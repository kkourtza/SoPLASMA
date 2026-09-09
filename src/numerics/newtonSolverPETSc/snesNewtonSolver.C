/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

\*---------------------------------------------------------------------------*/

#include "snesNewtonSolver.H"
#include "snesBridge.H"
#include "blockMatrixCOO.H"
#include "addToRunTimeSelectionTable.H"
#include <chrono>
#include <sstream>

#include "fvCFD.H"
#include "electromagneticsModel.H"
#include "singleRegionPoisson.H"
#include "plasmaSpecies.H"
#include "plasmaTransport.H"
#include "plasmaTransportModel.H"
#include "driftDiffusion.H"
#include "plasmaEnergy.H"
#include "localEnergyEnergyModel.H"

namespace Foam
{
    defineTypeNameAndDebug(snesNewtonSolver, 0);
    addToRunTimeSelectionTable(plasmaNewtonSolver, snesNewtonSolver, dictionary);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

// TEMPORARY (2026-09-09): L2 norm over this rank's own cells, for the
// residual-vs-terms diagnostic in residualCallback.
Foam::scalar blockNorm(const Foam::scalarField& f, const Foam::label nc)
{
    Foam::scalar s = 0;
    for (Foam::label c = 0; c < nc; ++c) { s += f[c]*f[c]; }
    return Foam::sqrt(Foam::returnReduce(s, Foam::sumOp<Foam::scalar>()));
}

// Global RMS -- a PER-CELL magnitude, so it does not drift with mesh size the
// way an L2 norm does. gAverage is a collective and is called unconditionally
// (rule 31: never behind a local guard).
Foam::scalar rmsOf(const Foam::scalarField& p)
{
    return Foam::sqrt(Foam::gAverage(Foam::sqr(p)));
}

// A scale must never be zero (a cold all-zero field would divide by it), and
// falling back to 1 makes the scaling the identity rather than a trap.
Foam::scalar safeScale(const Foam::scalar s)
{
    return (s > Foam::SMALL ? s : 1.0);
}

// Everything the residual callback needs, gathered ONCE in solveOuterStep()
// -- not reconstructed per call, since SNES calls the residual repeatedly
// (real Newton steps and internal matrix-free Jacobian-vector products).
struct ResidualContext
{
    Foam::electromagneticsModel* em;
    Foam::singleRegionPoisson* srp;     // same object as em, checked+cast once
    Foam::plasmaSpecies* species;
    Foam::plasmaTransport* transport;
    Foam::localEnergyEnergyModel* lmea; // nullptr under LFA
    Foam::label nCellsLocal;
    Foam::label nSpecies;
    bool hasEnergy;

    // PER-BLOCK SCALING, fixed for the whole Newton solve (recomputing it
    // per residual call would move the target and the problem would stop
    // being stationary). Block order matches the DOF layout:
    // [0]=Poisson, [1..nSpecies]=species, [last]=energy.
    //
    // Measured 2026-09-09, and this is why it exists: unscaled, the species
    // ddt blocks carry ||F|| ~ 2e24 while the Poisson block carries 1.2e-3.
    // A single combined L2 norm is then 100% species -- with rtol 1e-8 the
    // convergence target sits 19 orders of magnitude ABOVE the entire
    // Poisson residual, so SNES could not see Poisson at all. This is the
    // `typ u` per-component scaling Knoll & Keyes call for in section 2.3.1
    // (eq. 13) of J. Comput. Phys. 193 (2004) 357-397.
    const Foam::scalarList* sX; // state scale:    x_petsc = x_phys / sX
    const Foam::scalarList* sF; // residual scale: F_petsc = F_phys / sF
};

// Block layout, one instance per rank: [ePotential][n_0]..[n_{nSpecies-1}]
// [nEps_e if hasEnergy] -- flat concatenation, no VecNest, exactly the
// pattern the 2-field proof of concept validated in both serial and
// parallel (doc/newton-outer-solver-design.md).
void residualCallback
(
    int n,
    const double* x,
    double* F,
    void* userDataVoid
)
{
    using namespace Foam;

    // TEMPORARY diagnostic (2026-09-09): rule 15/27 check on my OWN residual.
    // ||F|| in isolation says nothing -- the species ddt term alone is
    // n/dt ~ 1e17/1e-12 = 1e29, so a ||F|| of 1e24 may be a healthy 1e-5
    // RELATIVE imbalance. Prints each block's ||F|| next to the norms of its
    // OWN terms, for the first few calls only.
    static label callCount = 0;
    ++callCount;
    const bool diagThisCall = (callCount <= 3);
    if (diagThisCall)
    {
        Pout<< "[diag] residualCallback call #" << callCount << endl;
    }

    // TEMPORARY per-phase TIMERS (2026-09-09, REMOVE AFTER). One Newton
    // timestep costs ~97 s against Picard's 44 steps/s, and JFNK needs one
    // residual per GMRES iteration. Each call does FOUR full matrix
    // discretisations plus a per-cell stiff chemistry ODE -- the chemistry is
    // the SUSPECT, not the conclusion, so this measures instead of assuming.
    static double tUnpack = 0, tDerived = 0, tTransport = 0,
                  tChem = 0, tEnergy = 0, tResid = 0;
    const auto clockNow = []{ return std::chrono::steady_clock::now(); };
    const auto secSince = [](const std::chrono::steady_clock::time_point t0)
    {
        return std::chrono::duration<double>
        (
            std::chrono::steady_clock::now() - t0
        ).count();
    };
    auto tPhase = clockNow();

    ResidualContext& ctx = *static_cast<ResidualContext*>(userDataVoid);
    const label nc = ctx.nCellsLocal;
    electromagneticsModel& em = *ctx.em;
    singleRegionPoisson& srp = *ctx.srp;
    plasmaSpecies& species = *ctx.species;
    plasmaTransport& transport = *ctx.transport;

    // ---- 1. Unpack the trial state into the REAL, live OpenFOAM fields.
    // These are the SAME fields the Picard branch reads/writes -- nothing
    // is reconstructed, so every existing accessor (mu(), D(), chemP(),
    // Psrc(), ...) reflects THIS trial state once step 2-5 below refresh
    // the models built on top of them.
    const scalarList& sX = *ctx.sX;
    const scalarList& sF = *ctx.sF;

    volScalarField& ePotential = em.ePotentialRef();
    {
        scalarField& f = ePotential.primitiveFieldRef();
        for (label c = 0; c < nc; ++c) { f[c] = x[c]*sX[0]; }
        ePotential.correctBoundaryConditions();
    }

    // INTERNAL VALUES ONLY HERE. The species and nEps boundary conditions are
    // corrected LATER, at step 4b, and the ordering is load-bearing:
    // electronDDWallFluxMixed and the ddWallFlux family read the LIVE phiE and
    // patch mobility, so correcting them here -- before updateDerivedFields()
    // rebuilds phiE from the new ePotential and before transportModel.correct()
    // refreshes mu -- gave patch VALUES from iterate k-1. That makes
    // F = F(u_k, u_{k-1}) rather than F(u), which is the same purity defect as
    // the meanE ordering and the LFA seed, and it corrupts the matrix-free
    // Jacobian in the same way. Found by review, 2026-09-09.
    for (label s = 0; s < ctx.nSpecies; ++s)
    {
        volScalarField& ns = species.numberDensity(s);
        scalarField& f = ns.primitiveFieldRef();
        const label off = (1 + s)*nc;
        for (label c = 0; c < nc; ++c) { f[c] = x[off + c]*sX[1 + s]; }
    }

    if (ctx.hasEnergy)
    {
        volScalarField& nEps = ctx.lmea->nEpsRef();
        scalarField& f = nEps.primitiveFieldRef();
        const label off = (1 + ctx.nSpecies)*nc;
        const scalar sxE = sX[1 + ctx.nSpecies];
        for (label c = 0; c < nc; ++c) { f[c] = x[off + c]*sxE; }
    }

    // ---- 1b. chargeDensity is DERIVED from the species densities, so it must
    // be rebuilt from the trial state before the Poisson residual reads it.
    // Without this the Poisson block sees a chargeDensity left over from
    // whenever the Picard path last updated it, which does not merely corrupt
    // the JACOBIAN's Poisson<->species coupling -- it corrupts F itself, so
    // Newton converges to the root of a Poisson-with-LAGGED-source system.
    // That is precisely the segregated coupling this solver exists to remove,
    // and it is why agreeing with Picard was NOT evidence of correctness.
    // Found 2026-09-09.
    species.updateChargeDensity();

    tUnpack += secSince(tPhase); tPhase = clockNow();

    // ---- 2. Refresh E, Emag, phiE, reducedE from the new ePotential --
    // exactly what singleRegionPoisson::solve() does internally before
    // building its own matrix.
    srp.updateDerivedFields();

    tDerived += secSince(tPhase); tPhase = clockNow();

    // ---- 3. Refresh meanE (and T_, and the energy coefficients) FIRST, so
    // that the transport and rate lookups below key on THIS trial state.
    //
    // ORDER IS LOAD-BEARING, and it was wrong until 2026-09-09. meanE_ is
    // recomputed only by updateDerived(), reachable only from
    // localEnergyEnergyModel::correct(). That call used to sit AFTER the
    // transport and chemistry refreshes, so mu/D and every reaction rate --
    // all keyed on the "meanE" field -- were evaluated against the PREVIOUS
    // residual call's mean energy. F was therefore a function of (u_k,
    // u_{k-1}), not of u, which silently corrupts every matrix-free
    // difference quotient: the Jacobian-vector product differences a
    // quantity that depends on evaluation HISTORY.
    //
    // eEqn()/updateSources() still runs AFTER the chemistry refresh (step 5),
    // because Psrc_/Lsrc_ depend on the chemistry sources under
    // `energySource chemistry`. Only correct() moved up here.
    if (ctx.hasEnergy)
    {
        ctx.lmea->correct();
    }

    // ---- 4. Refresh each species' transport coefficients (mu, D -- looked
    // up against the just-refreshed reducedE/meanE) from the new densities.
    for (label s = 0; s < ctx.nSpecies; ++s)
    {
        transport.transportModel(s).correct();
    }

    // ---- 4b. NOW correct the species/nEps boundary conditions, with phiE
    // (step 2), meanE (step 3) and mu/D (step 4) all rebuilt from THIS trial
    // state. See the note at step 1 for why this cannot happen earlier.
    for (label s = 0; s < ctx.nSpecies; ++s)
    {
        species.numberDensity(s).correctBoundaryConditions();
    }
    if (ctx.hasEnergy)
    {
        ctx.lmea->nEpsRef().correctBoundaryConditions();
    }

    tTransport += secSince(tPhase); tPhase = clockNow();

    // ---- 4. Refresh chemP_/chemL_ for the new species densities. Needs
    // real fvScalarMatrix objects to hand to refreshChemistrySources() (the
    // function it wraps writes species source terms into them the same way
    // the real Picard branch does) -- built here and simply discarded
    // afterward, never .solve()'d or .residual()'d. .nEqn().ptr() transfers
    // ownership out of the tmp<> (NOT `new fvScalarMatrix(nEqn())`, which
    // copy-constructs from a temporary and left a dangling internal
    // reference -- a real bug caught via gdb, 2026-09-09).
    {
        List<autoPtr<fvScalarMatrix>> eqns(ctx.nSpecies);
        for (label s = 0; s < ctx.nSpecies; ++s)
        {
            eqns[s].reset(transport.transportModel(s).nEqn().ptr());
        }
        const volScalarField& ne = species.numberDensity(species.electronSpeciesID());

        // NOTE, measured 2026-09-09 and recorded so it is not re-tried:
        // freezing chemistry here (holding chemP_/chemL_ across the Newton
        // solve) was tested as a DIAGNOSTIC for whether this path makes F
        // history-dependent -- plasmaTransport.C:4310-4327 does increment
        // chemOuterCount_ and overwrite chemSrcPrev_ on every call. It makes
        // NO measurable difference: the KSP residual history matched to 8
        // digits (35.95884133842 frozen vs 35.95884334902 live). So chemistry
        // history dependence is NOT the cause of the linear-solve failure at
        // a developed discharge state. The switch was removed again rather
        // than left as a code path nothing needs.
        transport.refreshChemistrySources(eqns, ne, em.Emag());
    }

    tChem += secSince(tPhase); tPhase = clockNow();

    // ---- 5. Refresh the energy model's Psrc_/Lsrc_, which depend on the
    // chemistry sources just refreshed above. eEqn()'s updateSources() call
    // does it; the returned matrix is discarded, only that side effect is
    // used, exactly as the chemistry step discards its matrices.
    //
    // correct() is NOT called here any more -- it moved to step 3, because
    // it owns the meanE refresh that steps 3-4 must key on. See there.
    if (ctx.hasEnergy)
    {
        localEnergyEnergyModel& lmea = *ctx.lmea;
        tmp<fvScalarMatrix> tDiscard = lmea.eEqn();
    }

    tEnergy += secSince(tPhase); tPhase = clockNow();

    // ---- 6. THE RESIDUAL. Explicit fvc:: evaluation throughout -- never
    // fvm::+.residual(), confirmed broken in parallel for that usage
    // independent of PETSc (see
    // fvmatrix-residual-broken-in-parallel-use-fvc-instead memory).

    // Poisson. Handles BOTH schemes -- the case's own em.PoissonScheme()
    // decides, so a case running `explicit` gets the right equation too,
    // not just the semiImplicit form this was designed against.
    {
        const dimensionedScalar& epsilon = em.epsilon();
        const volScalarField& chargeDensity = em.chargeDensity();

        if (em.PoissonScheme() == "semiImplicit")
        {
            const dimensionedScalar dt = ePotential.mesh().time().deltaT();
            const volScalarField effEps
            (
                epsilon + dt*transport.electricalConductivity()
            );
            const volScalarField rhsSource
            (
                -chargeDensity - dt*transport.diffusiveChargeSource()
            );

            tmp<volScalarField> tLap = fvc::laplacian(effEps, ePotential);
            const scalarField& lap = tLap().primitiveField();
            const scalarField& rhs = rhsSource.primitiveField();
            for (label c = 0; c < nc; ++c) { F[c] = (lap[c] - rhs[c])/sF[0]; }

            if (diagThisCall)
            {
                scalarField fBlock(nc);
                for (label c = 0; c < nc; ++c) { fBlock[c] = F[c]; }
                Pout<< "[diag]   Poisson: |lap|=" << blockNorm(lap, nc)
                    << " |rhs|=" << blockNorm(rhs, nc)
                    << " sF=" << sF[0]
                    << " |F_scaled|=" << blockNorm(fBlock, nc) << endl;
            }
        }
        else
        {
            tmp<volScalarField> tLap = fvc::laplacian(epsilon, ePotential);
            const scalarField& lap = tLap().primitiveField();
            const scalarField& rho = chargeDensity.primitiveField();
            for (label c = 0; c < nc; ++c) { F[c] = (lap[c] + rho[c])/sF[0]; }

            if (diagThisCall)
            {
                scalarField fBlock(nc);
                for (label c = 0; c < nc; ++c) { fBlock[c] = F[c]; }
                Pout<< "[diag]   Poisson(explicit): |lap|=" << blockNorm(lap, nc)
                    << " |rho|=" << blockNorm(rho, nc)
                    << " sX=" << sX[0] << " sF=" << sF[0]
                    << " |F_scaled|=" << blockNorm(fBlock, nc) << endl;
            }
        }
    }

    // Species. R_s = ddt(n_s) + div(phi_s, n_s) - lap(D_s, n_s)
    //              - (chemP_s - chemL_s * n_s)
    // phi_s reconstructed EXACTLY as driftDiffusion::convectivePhi() builds
    // it (Z * interpolate(mu) * phiE) -- same formula, same coefficient
    // fields (mu()/D(), just refreshed at step 3), not a re-derivation.
    for (label s = 0; s < ctx.nSpecies; ++s)
    {
        const volScalarField& ns = species.numberDensity(s);
        const plasmaTransportModel& model = transport.transportModel(s);
        const scalar Z = species.speciesChargeNumber(s);
        const word sName = species.speciesNames()[s];

        // phi_s is the CARRIER flux (mobility x field), not yet the total
        // particle flux -- fvc::div(phi_s) alone would silently drop the
        // n_s multiplication entirely (a real bug caught before this ever
        // ran on a real case). fvc::flux(phi, n, name) is the genuine
        // explicit equivalent of fvm::div(phi, n): it interpolates n to
        // faces with the SAME scheme fvm::div would use and forms the
        // actual total flux, which fvc::div then integrates.
        tmp<surfaceScalarField> tPhi
        (
            Z*fvc::interpolate(model.mu())*em.phiE()
        );
        tmp<surfaceScalarField> tTotalFlux = fvc::flux
        (
            tPhi(), ns, "div(phi_" + sName + ",n_" + sName + ")"
        );
        tmp<volScalarField> tDdt = fvc::ddt(ns);
        tmp<volScalarField> tDiv = fvc::div(tTotalFlux());
        tmp<volScalarField> tLap = fvc::laplacian(model.D(), ns);

        const scalarField& ddtF = tDdt().primitiveField();
        const scalarField& divF = tDiv().primitiveField();
        const scalarField& lapF = tLap().primitiveField();
        const scalarField& P = transport.chemP(s);
        const scalarField& Lr = transport.chemL(s);
        const scalarField& nsF = ns.primitiveField();

        const label off = (1 + s)*nc;
        const scalar sfS = sF[1 + s];
        for (label c = 0; c < nc; ++c)
        {
            F[off + c] =
                (ddtF[c] + divF[c] - lapF[c] - P[c] + Lr[c]*nsF[c])/sfS;
        }

        if (diagThisCall)
        {
            scalarField fBlock(nc), lossF(nc);
            for (label c = 0; c < nc; ++c)
            {
                fBlock[c] = F[off + c];
                lossF[c] = Lr[c]*nsF[c];
            }
            Pout<< "[diag]   n_" << sName
                << ": |n|=" << blockNorm(nsF, nc)
                << " |ddt|=" << blockNorm(ddtF, nc)
                << " |div|=" << blockNorm(divF, nc)
                << " |lap|=" << blockNorm(lapF, nc)
                << " |chemP|=" << blockNorm(P, nc)
                << " |chemL*n|=" << blockNorm(lossF, nc)
                << " sX=" << sX[1 + s] << " sF=" << sfS
                << " |F_scaled|=" << blockNorm(fBlock, nc) << endl;
        }
    }

    // Electron energy (LMEA). R = ddt(nEps) + div(phiEps,nEps)
    //                            - lap(DEpsEff,nEps) + Lsrc*nEps - Psrc
    // matching eEqn()'s own assembly (localEnergyEnergyModel.C) exactly,
    // reusing muEpsEff()/DEpsEff() (already energyFactor-scaled) and
    // Psrc()/Lsrc() (already chemistry-ODE-aware, from step 5).
    if (ctx.hasEnergy)
    {
        localEnergyEnergyModel& lmea = *ctx.lmea;
        const volScalarField& nEps = lmea.nEps();
        const scalar Ze = species.speciesChargeNumber(species.electronSpeciesID());

        // Explicit scheme name, matching eEqn()'s own convention: muEps has
        // no fvSchemes entry of its own, so the energy model's real
        // assembly borrows the electron mobility's ("interpolate(mu_e)")
        // rather than have every case declare one for a field it never
        // asked for -- see localEnergyEnergyModel.C's own comment on this.
        tmp<surfaceScalarField> tPhiEps
        (
            Ze*fvc::interpolate(lmea.muEpsEff(), "interpolate(mu_e)")*em.phiE()
        );
        // Same fvc::flux() correction as the species loop above -- phiEps
        // is the carrier flux, not the total one.
        tmp<surfaceScalarField> tTotalFluxEps = fvc::flux
        (
            tPhiEps(), nEps, "div(phi_e,n_e)"
        );
        tmp<volScalarField> tDdt = fvc::ddt(nEps);
        tmp<volScalarField> tDiv = fvc::div(tTotalFluxEps());
        // Explicit override, matching eEqn()'s own convention: DEpsEff is
        // registered as "DEps" (no matching fvSchemes entry of its own),
        // so the real assembly borrows the electron diffusivity's key.
        tmp<volScalarField> tLap = fvc::laplacian
        (
            lmea.DEpsEff(), nEps, "laplacian(D_e,n_e)"
        );

        const scalarField& ddtF = tDdt().primitiveField();
        const scalarField& divF = tDiv().primitiveField();
        const scalarField& lapF = tLap().primitiveField();
        const scalarField& Psrc = lmea.Psrc().primitiveField();
        const scalarField& Lsrc = lmea.Lsrc().primitiveField();
        const scalarField& nEpsF = nEps.primitiveField();

        const label off = (1 + ctx.nSpecies)*nc;
        const scalar sfE = sF[1 + ctx.nSpecies];
        for (label c = 0; c < nc; ++c)
        {
            F[off + c] =
                (ddtF[c] + divF[c] - lapF[c] + Lsrc[c]*nEpsF[c] - Psrc[c])
               /sfE;
        }

        if (diagThisCall)
        {
            scalarField fBlock(nc), lossF(nc);
            for (label c = 0; c < nc; ++c)
            {
                fBlock[c] = F[off + c];
                lossF[c] = Lsrc[c]*nEpsF[c];
            }
            Pout<< "[diag]   nEps_e"
                << ": |nEps|=" << blockNorm(nEpsF, nc)
                << " |ddt|=" << blockNorm(ddtF, nc)
                << " |div|=" << blockNorm(divF, nc)
                << " |lap|=" << blockNorm(lapF, nc)
                << " |Psrc|=" << blockNorm(Psrc, nc)
                << " |Lsrc*nEps|=" << blockNorm(lossF, nc)
                << " sX=" << sX[1 + ctx.nSpecies] << " sF=" << sfE
                << " |F_scaled|=" << blockNorm(fBlock, nc) << endl;
        }
    }

    tResid += secSince(tPhase);

    if (callCount % 2000 == 0 || callCount <= 3)
    {
        const double tot = tUnpack + tDerived + tTransport + tChem
                         + tEnergy + tResid;
        Pout<< "[timing] calls=" << callCount
            << "  total=" << tot << " s"
            << "  | unpack=" << tUnpack
            << " derived=" << tDerived
            << " transport=" << tTransport
            << " CHEM=" << tChem
            << " energy=" << tEnergy
            << " residual=" << tResid
            << "  (per call: " << tot/callCount << " s)" << endl;
    }
}

// Physics-based (Knoll & Keyes) preconditioner: approximate J^-1*b by ONE
// block Gauss-Seidel sweep -- Poisson, then species, then energy, each a
// native linear solve of the equation's OWN linearization (frozen
// coefficients at the current outer iterate, refreshed every residual
// call already), exactly mirroring the REAL Picard sweep's own order and
// operators. This is deliberately the SAME sequence
// soPlasmaFoam.C's Picard branch runs -- the old segregated solver
// repurposed as the preconditioner, not duplicated logic -- except each
// equation solves for a CORRECTION (homogeneous BCs, zeroed before every
// application) against a synthetic right-hand side from GMRES, not the
// real field against the real source.
struct PCContext
{
    Foam::electromagneticsModel* em;
    Foam::singleRegionPoisson* srp;
    Foam::plasmaSpecies* species;
    Foam::plasmaTransport* transport;
    Foam::localEnergyEnergyModel* lmea;
    Foam::label nCellsLocal;
    Foam::label nSpecies;
    bool hasEnergy;

    // Persistent scratch correction fields, constructed ONCE in
    // solveOuterStep() (mirroring dphi/dphi2 in the 2-field proof of
    // concept) -- each copies its REAL counterpart as a BC prototype (so
    // it inherits the right patch TYPES) then gets zeroed every
    // application: a correction must vanish wherever the real field is
    // Dirichlet-pinned.
    Foam::volScalarField* dePotential;
    Foam::PtrList<Foam::volScalarField>* dSpecies;
    Foam::volScalarField* dEnergy; // nullptr if !hasEnergy

    // Same per-block scaling the residual uses (see ResidualContext). GMRES
    // hands b in SCALED residual units and expects y in SCALED state units,
    // while the equations solved in between are PHYSICAL -- so b is taken
    // back to physical units on the way in, and the correction is scaled on
    // the way out.
    const Foam::scalarList* sX;
    const Foam::scalarList* sF;
};

void pcApplyCallback
(
    int n,
    const double* b,
    double* y,
    void* userDataVoid
)
{
    using namespace Foam;

    // TEMPORARY diagnostic (2026-09-09), paired with residualCallback's own.
    // Capped: this is called once per GMRES iteration, so an uncapped print
    // floods a long run.
    static label pcCallCount = 0;
    ++pcCallCount;
    if (pcCallCount <= 10)
    {
        Pout<< "[diag] pcApplyCallback call #" << pcCallCount << endl;
    }

    PCContext& ctx = *static_cast<PCContext*>(userDataVoid);
    const label nc = ctx.nCellsLocal;

    // THE VECTOR MAY BE SHORTER THAN THE FULL FIELD LAYOUT, and if it is we
    // must not touch it as if it were not.
    //
    // SNESVINEWTONRSLS is the REDUCED-SPACE active-set method: components
    // sitting at their bound are removed from the Newton system, so the KSP --
    // and therefore this preconditioner -- is handed a vector covering only
    // the INACTIVE components (virs.c). This routine indexes the full
    // [ePotential][n_0..n_{N-1}][nEps] layout, so on a reduced vector it wrote
    // past the end. Found 2026-09-09 by PETSc's own `-malloc_debug`, which
    // named it in one run: "error detected in PCApply_Shell() ... is corrupted
    // (probably write past end of array)". A backtrace of the eventual abort
    // was useless, because corruption is only noticed later, at free().
    //
    // Falling back to the IDENTITY is safe and legitimate -- a preconditioner
    // need only be an approximation, and M = I is the weakest valid one, so
    // the solve continues unpreconditioned for that application instead of
    // corrupting memory. It is NOT the real answer: the proper route is an
    // assembled Pmat with PCFIELDSPLIT, which SNESVI supports natively via
    // PCFieldSplitRestrictIS. See doc/newton-outer-solver-design.md.
    const label nExpected = (1 + ctx.nSpecies + (ctx.hasEnergy ? 1 : 0))*nc;
    if (label(n) != nExpected)
    {
        static bool warned = false;
        if (!warned)
        {
            warned = true;
            Pout<< "snesNewtonSolver: preconditioner received a REDUCED system"
                << " (" << n << " of " << nExpected << " unknowns -- the"
                << " bounded solver's active set). Falling back to the"
                << " identity for those applications; the physics-based"
                << " preconditioner cannot be applied to a restricted"
                << " vector. See doc/newton-outer-solver-design.md." << endl;
        }
        for (int i = 0; i < n; ++i) { y[i] = b[i]; }
        return;
    }
    // Note: unlike residualCallback, the PC correction fields never need
    // E/phiE/reducedE refreshed (they are pure linear corrections, not
    // trial states the real derived fields must track), so ctx.srp is
    // unused here.
    electromagneticsModel& em = *ctx.em;
    plasmaSpecies& species = *ctx.species;
    plasmaTransport& transport = *ctx.transport;

    // Unpack b into per-block RHS fields. The VALUES are plain raw doubles
    // (residualCallback's own convention -- packed via .primitiveField(),
    // dimension tag discarded), but their MAGNITUDES are those of the real
    // dimensioned residual quantities, so each block's field must carry the
    // matching real dimension for the fvm:: equations below (built from
    // REAL dimensioned coefficients: epsilon, mu, D, ...) to dimension-check.
    // Poisson residual's dimension is chargeDensity's (see residualCallback:
    // F = lap(effEps,ePotential) - rhsSource, both chargeDensity-dimensioned).
    volScalarField bPoisson
    (
        IOobject("bPoisson", em.ePotential().time().timeName(),
                 em.ePotential().mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
        em.ePotential().mesh(), dimensionedScalar(em.chargeDensity().dimensions(), Zero)
    );
    const scalarList& sX = *ctx.sX;
    const scalarList& sF = *ctx.sF;
    {
        scalarField& f = bPoisson.primitiveFieldRef();
        for (label c = 0; c < nc; ++c) { f[c] = b[c]*sF[0]; }
    }

    PtrList<volScalarField> bSpecies(ctx.nSpecies);
    for (label s = 0; s < ctx.nSpecies; ++s)
    {
        bSpecies.set
        (
            s,
            new volScalarField
            (
                IOobject("bSpecies_" + Foam::name(s), em.ePotential().time().timeName(),
                         em.ePotential().mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                em.ePotential().mesh(),
                dimensionedScalar
                (
                    species.numberDensity(s).dimensions()/dimTime, Zero
                )
            )
        );
        scalarField& f = bSpecies[s].primitiveFieldRef();
        const label off = (1 + s)*nc;
        const scalar sfS = sF[1 + s];
        for (label c = 0; c < nc; ++c) { f[c] = b[off + c]*sfS; }
    }

    autoPtr<volScalarField> bEnergy;
    if (ctx.hasEnergy)
    {
        bEnergy.reset
        (
            new volScalarField
            (
                IOobject("bEnergy", em.ePotential().time().timeName(),
                         em.ePotential().mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
                em.ePotential().mesh(),
                dimensionedScalar(ctx.lmea->nEps().dimensions()/dimTime, Zero)
            )
        );
        scalarField& f = bEnergy->primitiveFieldRef();
        const label off = (1 + ctx.nSpecies)*nc;
        const scalar sfE = sF[1 + ctx.nSpecies];
        for (label c = 0; c < nc; ++c) { f[c] = b[off + c]*sfE; }
    }

    // ---- Step 1: Poisson correction. Same operator/coefficient as the
    // residual (effEps for semiImplicit, epsilon for explicit) -- already
    // genuinely linear in ePotential, so no further linearisation needed.
    volScalarField& dePotential = *ctx.dePotential;
    dePotential.primitiveFieldRef() = Zero;
    dePotential.correctBoundaryConditions();
    {
        const dimensionedScalar& epsilon = em.epsilon();
        if (em.PoissonScheme() == "semiImplicit")
        {
            const dimensionedScalar dt = em.ePotential().mesh().time().deltaT();
            const volScalarField effEps
            (
                epsilon + dt*transport.electricalConductivity()
            );
            // Explicit scheme name: d_ePotential is a differently-named
            // field, so the auto-derived key
            // "laplacian((epsilon+...),d_ePotential)" does not match the
            // case's own "...,ePotential)" entry -- borrow it explicitly,
            // same pattern as the energy equation's muEps/DEps overrides.
            fvScalarMatrix eqn
            (
                fvm::laplacian
                (
                    effEps, dePotential,
                    "laplacian((epsilon+(deltaT*electricalConductivity)),ePotential)"
                )
             == bPoisson
            );
            eqn.solve(em.ePotential().mesh().solution().subDict("solvers").subDict("ePotential"));
        }
        else
        {
            fvScalarMatrix eqn
            (
                fvm::laplacian(epsilon, dePotential, "laplacian(epsilon,ePotential)")
             == bPoisson
            );
            eqn.solve(em.ePotential().mesh().solution().subDict("solvers").subDict("ePotential"));
        }
    }

    // ---- Step 2: species corrections, Gauss-Seidel (using the JUST-solved
    // dePotential's contribution to each species' convective term, exactly
    // as the real Picard sweep sees the freshly-solved ePotential before
    // building the species matrices).
    for (label s = 0; s < ctx.nSpecies; ++s)
    {
        volScalarField& dn = (*ctx.dSpecies)[s];
        dn.primitiveFieldRef() = Zero;
        dn.correctBoundaryConditions();

        const plasmaTransportModel& model = transport.transportModel(s);
        const scalar Z = species.speciesChargeNumber(s);
        const word sName = species.speciesNames()[s];

        // Off-diagonal Jacobian contribution from dePotential: the species
        // convective flux depends on phiE, which the just-solved correction
        // perturbs. -Z*interpolate(mu)*fvc::snGrad(dePotential)*magSf is the
        // linearised flux PERTURBATION (same construction as em's own
        // phiE(), applied to the CORRECTION instead of the real potential).
        tmp<surfaceScalarField> tDPhiE
        (
            -fvc::snGrad(dePotential, "snGrad(ePotential)")
            *em.ePotential().mesh().magSf()
        );
        tmp<surfaceScalarField> tDPhi
        (
            Z*fvc::interpolate(model.mu())*tDPhiE()
        );
        tmp<surfaceScalarField> tDFlux = fvc::flux
        (
            tDPhi(), species.numberDensity(s),
            "div(phi_" + sName + ",n_" + sName + ")"
        );

        const volScalarField rhs_s
        (
            IOobject("rhs_" + sName, em.ePotential().time().timeName(),
                     em.ePotential().mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
            bSpecies[s] - fvc::div(tDFlux())
        );

        tmp<surfaceScalarField> tPhiFrozen
        (
            Z*fvc::interpolate(model.mu())*em.phiE()
        );
        // fvm::Sp needs a volScalarField coefficient, not the raw
        // scalarField chemL(s) returns (List<scalarField>, no GeometricField
        // wrapper) -- wrap it. chemL is a rate coefficient (chemL*n must
        // match ddt(n)'s dimension), so it needs 1/time, not dimensionless,
        // to dimension-check against the ddt/div/laplacian terms below.
        volScalarField chemLField
        (
            IOobject("chemLField_" + sName, em.ePotential().time().timeName(),
                     em.ePotential().mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
            em.ePotential().mesh(), dimensionedScalar(dimless/dimTime, Zero)
        );
        chemLField.primitiveFieldRef() = transport.chemL(s);

        fvScalarMatrix eqn
        (
            fvm::ddt(dn)
          + fvm::div(tPhiFrozen(), dn, "div(phi_" + sName + ",n_" + sName + ")")
          - fvm::laplacian
            (
                model.D(), dn, "laplacian(D_" + sName + ",n_" + sName + ")"
            )
          + fvm::Sp(chemLField, dn)
         ==
            rhs_s
        );
        eqn.solve(em.ePotential().mesh().solution().subDict("solvers").subDict("n_" + sName));
    }

    // ---- Step 3: energy correction, using the LATEST species corrections'
    // contribution to Joule heating is NOT included (Psrc/Lsrc are treated
    // as frozen, matching how the real energy source already comes from a
    // separately-integrated stiff ODE, not a simple algebraic function of
    // n_e this preconditioner could usefully differentiate through).
    if (ctx.hasEnergy)
    {
        localEnergyEnergyModel& lmea = *ctx.lmea;
        volScalarField& dEnergy = *ctx.dEnergy;
        dEnergy.primitiveFieldRef() = Zero;
        dEnergy.correctBoundaryConditions();

        const scalar Ze = species.speciesChargeNumber(species.electronSpeciesID());
        tmp<surfaceScalarField> tPhiEpsFrozen
        (
            Ze*fvc::interpolate(lmea.muEpsEff(), "interpolate(mu_e)")*em.phiE()
        );
        fvScalarMatrix eqn
        (
            fvm::ddt(dEnergy)
          + fvm::div(tPhiEpsFrozen(), dEnergy, "div(phi_e,n_e)")
          - fvm::laplacian(lmea.DEpsEff(), dEnergy, "laplacian(D_e,n_e)")
          + fvm::Sp(lmea.Lsrc(), dEnergy)
         ==
            *bEnergy
        );
        // Same fallback as the real plasmaEnergy::solve() (plasmaEnergy.C):
        // an explicit "nEps_<species>" entry wins if present, otherwise
        // borrow "n_<species>"'s settings -- cases are not required to
        // declare a solver for an equation enabled with one keyword.
        {
            const word eName = species.speciesNames()[species.electronSpeciesID()];
            const dictionary& solvers =
                em.ePotential().mesh().solution().subDict("solvers");
            const word key =
                solvers.found("nEps_" + eName, keyType::REGEX)
              ? word("nEps_" + eName)
              : word("n_" + eName);
            eqn.solve(solvers.subDict(key));
        }
    }

    // Pack the corrections back into y, in SCALED state units (the
    // corrections were solved for in physical units).
    {
        const scalarField& f = dePotential.primitiveField();
        for (label c = 0; c < nc; ++c) { y[c] = f[c]/sX[0]; }
    }
    for (label s = 0; s < ctx.nSpecies; ++s)
    {
        const scalarField& f = (*ctx.dSpecies)[s].primitiveField();
        const label off = (1 + s)*nc;
        const scalar sxS = sX[1 + s];
        for (label c = 0; c < nc; ++c) { y[off + c] = f[c]/sxS; }
    }
    if (ctx.hasEnergy)
    {
        const scalarField& f = ctx.dEnergy->primitiveField();
        const label off = (1 + ctx.nSpecies)*nc;
        const scalar sxE = sX[1 + ctx.nSpecies];
        for (label c = 0; c < nc; ++c) { y[off + c] = f[c]/sxE; }
    }
}

// Homogeneous BC-type list for a correction field: fixedValue-0 wherever
// the real field's patch currently fixes a value (Dirichlet-like), plain
// zeroGradient elsewhere. Deliberately NEVER copies the real field's own
// boundary CONDITION OBJECTS -- custom physics BCs (e.g. ddWallFluxMixed)
// derive their species identity from the FIELD'S OWN registered name
// ("n_e" -> "e" via a prefix strip), which a "d_"-prefixed correction
// field breaks outright (a real crash found and fixed, 2026-09-09:
// "d_e not found in table" from plasmaSpecies::speciesID(), thrown by
// ddWallFluxMixedFvPatchScalarField::updateCoeffs() when
// correctBoundaryConditions() ran on the correction field). A correction
// has no wall-flux physics of its own -- only its Dirichlet/Neumann
// STRUCTURE needs to match, never the real BC's behaviour.
wordList homogeneousPatchTypes(const volScalarField& realField)
{
    wordList types(realField.boundaryField().size());
    forAll(realField.boundaryField(), patchi)
    {
        types[patchi] =
        (
            realField.boundaryField()[patchi].fixesValue()
          ? fixedValueFvPatchScalarField::typeName
          : zeroGradientFvPatchScalarField::typeName
        );
    }
    return types;
}

} // namespace

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

Foam::snesNewtonSolver::snesNewtonSolver
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    rtol_(dict.getOrDefault<scalar>("rtol", 1e-8)),
    maxIt_(dict.getOrDefault<label>("maxIt", 50)),
    mffdErr_(dict.getOrDefault<scalar>("mffdErr", 1e-5)),
    bounded_(dict.getOrDefault<bool>("bounded", true))
{
    // Lazy, ONCE-only: soPlasmaFoam's main() never calls initPetsc() itself
    // (this library is optionally loaded, so soPlasmaFoam must stay
    // PETSc-free when it is not), and by the time THIS constructor runs
    // (from inside plasmaSimulationControls reading, well after
    // setRootCase.H/createTime.H) MPI is already initialised by OpenFOAM
    // itself -- PetscInitialize attaches to it rather than calling MPI_Init
    // again, the normal way to embed PETSc in an already-MPI-parallel host.
    static bool petscInitialised = false;
    if (!petscInitialised)
    {
        initPetsc();
        petscInitialised = true;
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * //

void Foam::snesNewtonSolver::solveOuterStep
(
    electromagneticsModel& em,
    plasmaSpecies& species,
    plasmaTransport& transport,
    plasmaEnergy* energy
)
{
    if (!isA<singleRegionPoisson>(em))
    {
        FatalErrorInFunction
            << "outerSolver newton (type SNES) supports only"
            << " singleRegionPoisson so far -- got '"
            << em.type() << "'." << nl
            << "    multiRegionPoisson (dielectric regions) is not yet"
            << " supported; see doc/newton-outer-solver-design.md."
            << nl << exit(FatalError);
    }
    singleRegionPoisson& srp = refCast<singleRegionPoisson>(em);

    const label nSpecies = species.nSpecies();

    // COLD ALL-ZERO START: refuse it, loudly, rather than pretend to converge.
    //
    // Measured 2026-09-09 on grubert2009 from t=0, where every field starts
    // `internalField uniform 0`: every residual term is then exactly zero, so
    // SNES reports `CONVERGED_FNORM_ABS, 0 iterations` having done NOTHING --
    // a convergence that means the opposite of what it says. The step's own
    // clamp then jams the densities from 0 to the floor, a state change made
    // OUTSIDE the equations, and at the next step F is entirely the
    // time-derivative of that jump, which nothing in the equations can
    // balance (DIVERGED_LINEAR_SOLVE).
    //
    // This is rule 30's principle: the answer to a mechanism that fails on a
    // bad input is a GUARD that refuses the input, not a parallel mechanism.
    // The real fix is a formulation in which positivity is structural so no
    // clamp is needed -- see doc/newton-outer-solver-design.md.
    {
        bool allSpeciesZero = true;
        for (label s = 0; s < nSpecies; ++s)
        {
            if (rmsOf(species.numberDensity(s).primitiveField()) > SMALL)
            {
                allSpeciesZero = false;
            }
        }

        // Only a concern for the UNBOUNDED solve. With `bounded` on, x=0 is
        // INFEASIBLE, so SNESVI projects the initial guess into the feasible
        // region and the solution simply sits ON the constraint -- which is
        // exactly what an active-set method is for, and is the physically
        // right answer at t=0 (the density stays at its floor until
        // ionisation builds it).
        if (allSpeciesZero && !bounded_)
        {
            FatalErrorInFunction
                << "outerSolver newton cannot start from an all-zero state."
                << nl << nl
                << "    Every species number density is identically zero, so"
                << " the residual F is exactly zero and SNES would report"
                << " CONVERGED having solved nothing. The density clamp would"
                << " then move the state from outside the equations, and the"
                << " next step's residual would be unsolvable." << nl << nl
                << "    WHAT TO DO: either set `bounded true` in the newtonSolver"
                << " dict (the density floor then becomes a CONSTRAINT inside"
                << " the solve), or reach a physical state with `outerSolver"
                << " picard` first, then restart in newton mode from it --"
                << " `startFrom latestTime` in system/controlDict. The"
                << " densities must sit above the floor, not on it." << nl
                << "    See doc/newton-outer-solver-design.md (defect A)."
                << nl << exit(FatalError);
        }
    }

    for (label s = 0; s < nSpecies; ++s)
    {
        if (!isA<driftDiffusion>(transport.transportModel(s)))
        {
            FatalErrorInFunction
                << "outerSolver newton (type SNES) supports only"
                << " driftDiffusion transport so far -- species '"
                << species.speciesNames()[s] << "' uses '"
                << transport.transportModel(s).modelName() << "'." << nl
                << "    See doc/newton-outer-solver-design.md."
                << nl << exit(FatalError);
        }
    }

    localEnergyEnergyModel* lmea = nullptr;
    if (energy)
    {
        const localEnergyEnergyModel* clmea = energy->lmeaModel();
        if (!clmea)
        {
            FatalErrorInFunction
                << "outerSolver newton (type SNES) supports only LMEA"
                << " (localEnergyEnergyModel) so far -- this case's energy"
                << " model is not LMEA." << nl
                << "    See doc/newton-outer-solver-design.md."
                << nl << exit(FatalError);
        }
        // The underlying object is not actually const -- plasmaEnergy owns
        // it in a mutable PtrList<plasmaEnergyModel>; lmeaModel() only
        // returns a const view because plasmaEnergy exposes it through a
        // const method. Needed non-const here to call correct()/eEqn().
        lmea = const_cast<localEnergyEnergyModel*>(clmea);
    }

    const label nCellsLocal = mesh_.nCells();
    const label nFields = 1 + nSpecies + (lmea ? 1 : 0);
    const label nLocalTotal = nFields*nCellsLocal;

    // ---- BURN ANY ONE-SHOT MODEL INITIALISATION *BEFORE* THE SOLVE.
    //
    // localEnergyEnergyModel::correct() carries a one-shot LFA seed
    // (`if (seedFromLFA_) { seedFromLFA_ = false; nEps_ == meanE0*n_e; }`).
    // The Newton path's first correct() is INSIDE residualCallback, so that
    // seed fired on residual call #1 and not on calls #2+: F was literally a
    // DIFFERENT FUNCTION at the base point than at every perturbed point, so
    // the matrix-free product differenced F_2(u+hv) against F_1(u) and J*v was
    // meaningless. GMRES's own recurrence then could not agree with the
    // explicitly computed b - A*x, which is the "impossible" discrepancy
    // recorded in the design doc and wrongly attributed first to nonlinear
    // stiffness and then to conditioning.
    //
    // Measured 2026-09-09 at t=2e-8: the nEps block's scaled residual was
    // 172.46 on call #1 against 0.00086 on call #2 -- a factor of 2e5, and
    // ~97% of the reported initial norm. Calling correct() once here makes
    // the seed part of the STATE the solve starts from, which is what a
    // one-shot initialisation is supposed to be.
    //
    // NOT THE WHOLE STORY, and the rest is a PRODUCTION bug, not a Newton one:
    // that seed is gated on `time().timeIndex() == time().startTimeIndex()`
    // evaluated in the CONSTRUCTOR, before the time loop, so it is ALWAYS
    // true and every LMEA restart discards its stored energy state at the
    // first correct() -- on the Picard path too. Recorded in the design doc
    // for its own fix; deliberately not patched from here.
    if (lmea)
    {
        lmea->correct();
    }

    // ---- PER-BLOCK SCALES, computed ONCE for this whole Newton solve.
    // See ResidualContext's own comment for the measurement that forced
    // this. Each block's state scale is the RMS of its own field; each
    // block's residual scale is the RMS of that block's DOMINANT term,
    // measured 2026-09-09 to be ddt for every transported field and the
    // laplacian for Poisson.
    scalarList sX(nFields, 1.0);
    scalarList sF(nFields, 1.0);
    {
        const dimensionedScalar dt = mesh_.time().deltaT();

        sX[0] = safeScale(rmsOf(em.ePotential().primitiveField()));
        if (em.PoissonScheme() == "semiImplicit")
        {
            const volScalarField effEps
            (
                em.epsilon() + dt*transport.electricalConductivity()
            );
            sF[0] = safeScale
            (
                rmsOf(fvc::laplacian(effEps, em.ePotential())().primitiveField())
            );
        }
        else
        {
            sF[0] = safeScale
            (
                rmsOf
                (
                    fvc::laplacian(em.epsilon(), em.ePotential())()
                        .primitiveField()
                )
            );
        }

        for (label s = 0; s < nSpecies; ++s)
        {
            sX[1 + s] = safeScale(rmsOf(species.numberDensity(s).primitiveField()));
            sF[1 + s] = safeScale(sX[1 + s]/dt.value());
        }

        if (lmea)
        {
            sX[1 + nSpecies] = safeScale(rmsOf(lmea->nEps().primitiveField()));
            sF[1 + nSpecies] = safeScale(sX[1 + nSpecies]/dt.value());
        }
    }

    ResidualContext ctx
    {
        &em, &srp, &species, &transport, lmea,
        nCellsLocal, nSpecies, lmea != nullptr,
        &sX, &sF
    };

    // Persistent PC scratch fields -- constructed ONCE per outer step, each
    // copying its real counterpart as a BC prototype (correct patch TYPES)
    // then zeroed on every PC application (see pcApplyCallback).
    // Real DIMENSIONS (the correction combines with real-dimensioned
    // coefficients -- mu, D, epsilon, ... -- in pcApplyCallback's
    // off-diagonal terms), but homogeneous, physics-free BC TYPES only --
    // never the real field's own BC OBJECTS (see homogeneousPatchTypes()'s
    // comment for the crash that caused this).
    volScalarField dePotential
    (
        IOobject("d_ePotential", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_, dimensionedScalar(em.ePotential().dimensions(), Zero),
        homogeneousPatchTypes(em.ePotential())
    );
    PtrList<volScalarField> dSpecies(nSpecies);
    for (label s = 0; s < nSpecies; ++s)
    {
        const volScalarField& realN = species.numberDensity(s);
        dSpecies.set
        (
            s,
            new volScalarField
            (
                IOobject
                (
                    "d_" + species.speciesNames()[s], mesh_.time().timeName(),
                    mesh_, IOobject::NO_READ, IOobject::NO_WRITE
                ),
                mesh_, dimensionedScalar(realN.dimensions(), Zero),
                homogeneousPatchTypes(realN)
            )
        );
    }
    autoPtr<volScalarField> dEnergy;
    if (lmea)
    {
        const volScalarField& realNEps = lmea->nEps();
        dEnergy.reset
        (
            new volScalarField
            (
                IOobject("d_nEps_e", mesh_.time().timeName(), mesh_,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh_, dimensionedScalar(realNEps.dimensions(), Zero),
                homogeneousPatchTypes(realNEps)
            )
        );
    }

    PCContext pcCtx
    {
        &em, &srp, &species, &transport, lmea,
        nCellsLocal, nSpecies, lmea != nullptr,
        &dePotential, &dSpecies, dEnergy.get(),
        &sX, &sF
    };

    // Packed in SCALED state units -- everything PETSc sees is scaled.
    List<double> xBuf(nLocalTotal, Zero);
    {
        const scalarField& ePotF = em.ePotential().primitiveField();
        for (label c = 0; c < nCellsLocal; ++c) { xBuf[c] = ePotF[c]/sX[0]; }

        for (label s = 0; s < nSpecies; ++s)
        {
            const scalarField& nf = species.numberDensity(s).primitiveField();
            const label off = (1 + s)*nCellsLocal;
            const scalar sxS = sX[1 + s];
            for (label c = 0; c < nCellsLocal; ++c) { xBuf[off + c] = nf[c]/sxS; }
        }

        if (lmea)
        {
            const scalarField& nef = lmea->nEps().primitiveField();
            const label off = (1 + nSpecies)*nCellsLocal;
            const scalar sxE = sX[1 + nSpecies];
            for (label c = 0; c < nCellsLocal; ++c) { xBuf[off + c] = nef[c]/sxE; }
        }
    }

    // Physics-based (Knoll & Keyes) preconditioned matrix-free JFNK -- see
    // pcApplyCallback above.
    //
    // SUPERSEDED BY the -mat_mffd_err finding below (2026-09-09): the earlier
    // note here said bare/unpreconditioned JFNK "reliably hits
    // SNES_DIVERGED_LINEAR_SOLVE on the real system's first genuine Newton
    // iteration -- expected, matching the design doc's diagnosis". That
    // attribution was WRONG. The breakdown was the differencing step being
    // built from machine epsilon while F is only evaluable to ~1e-4, and it
    // occurred IDENTICALLY with the physics-based PC and with PCNONE -- which
    // is what proved the operator, not the preconditioner, was at fault. What
    // still holds: the PC is worth having (it converges the Krylov system in
    // 3-7 iterations). What does not: that bare JFNK's failure demonstrated
    // the system's stiffness.
    // setEnv("PETSC_OPTIONS",...) does NOT work here: PetscInitialize()
    // already ran once, in this object's CONSTRUCTOR (see there), which
    // consumed whatever PETSC_OPTIONS existed at PROGRAM START -- setting
    // the env var afterward from inside solveOuterStep() is never seen by
    // PETSc. setPetscOptions() inserts directly into PETSc's live options
    // database instead, which SNESSetFromOptions() (inside solveWithSNES)
    // picks up fresh on every call.
    const word petscOptions =
        "-snes_rtol " + Foam::name(rtol_)
      + " -snes_max_it " + Foam::name(maxIt_)
      + " -snes_monitor -snes_converged_reason"
      // -ksp_monitor deliberately NOT set: it prints per GMRES iteration and
      // floods a long run. -ksp_max_it caps the Krylov work per Newton step.
      + " -ksp_converged_reason -ksp_max_it 100"
      + " -mat_mffd_type wp"
      // Differencing step: PETSc's own default unless a case overrides it.
      // It is deliberately NOT set to a large value by default -- see
      // mffdErr_'s declaration for why that would have papered over a real
      // ordering bug rather than accommodating genuine noise in F.
      + (mffdErr_ > 0 ? " -mat_mffd_err " + Foam::name(mffdErr_) : word(""));
    setPetscOptions(petscOptions.c_str());

    // ---- LOWER BOUNDS, in SCALED units (everything PETSc sees is scaled).
    //
    // The species floor is the SAME number clampNumberDensities() applies
    // (plasmaSpecies::clampNumberDensity, active only where the floor is
    // positive), but imposed as a CONSTRAINT inside the solve instead of as a
    // mutation of the state afterwards. That difference is what makes F(u)=0
    // reachable: a post-hoc clamp moves the state to somewhere the residual
    // does not call a root, so Newton cannot converge to it.
    //
    // ePotential is genuinely unbounded (it is signed). nEps is floored at
    // zero: a negative energy density is meaningless, but its physical floor
    // is n_floor*meanE_min rather than a stated constant, so zero is the
    // defensible bound rather than an invented one.
    List<double> lowerBounds(nLocalTotal, snesNoBound);
    if (bounded_)
    {
        for (label c = 0; c < nCellsLocal; ++c)
        {
            lowerBounds[c] = snesNoBound;      // ePotential: free
        }

        for (label s = 0; s < nSpecies; ++s)
        {
            const scalar floorS = species.speciesMinNumberDensity(s);
            const label off = (1 + s)*nCellsLocal;
            const scalar sxS = sX[1 + s];
            const double lb =
                (floorS > 0 ? double(floorS/sxS) : 0.0);
            for (label c = 0; c < nCellsLocal; ++c)
            {
                lowerBounds[off + c] = lb;
            }
        }

        if (lmea)
        {
            const label off = (1 + nSpecies)*nCellsLocal;
            for (label c = 0; c < nCellsLocal; ++c)
            {
                lowerBounds[off + c] = 0.0;
            }
        }
    }

    // ---- ASSEMBLED PRECONDITIONING MATRIX (Pmat), in COO form.
    //
    // The diagonal blocks are the SAME fvm:: linearisations pcApplyCallback
    // solves, so the physics exists once, in OpenFOAM's language (rule 40).
    // What is NEW here is the CROSS-COUPLING, which the block Gauss-Seidel
    // shell could not represent at all and which PCFIELDSPLIT's Schur
    // complement needs.
    //
    // Assembled ONCE per outer step and held fixed for the solve -- a lagged
    // preconditioner. The matrix-free Amat still tracks the current iterate,
    // so this costs nothing in the converged answer (K&K 3.1/5.1).
    // PRIME THE MODELS FIRST. Everything the blocks below read -- mu, D,
    // chemP/chemL, Psrc/Lsrc, phiE, meanE -- is refreshed (and in the
    // chemistry's case SIZED) inside residualCallback. Assembling before any
    // residual evaluation read fields that did not exist yet and segfaulted
    // in the species block (localised with gdb + markers, 2026-09-09).
    //
    // One extra residual evaluation costs ~3 ms and buys a guarantee worth
    // more than that: Pmat and F then describe the SAME state.
    {
        List<double> primeF(nLocalTotal, Zero);
        residualCallback(int(nLocalTotal), xBuf.data(), primeF.data(), &ctx);
    }

    blockMatrixCOO coo(nCellsLocal, nFields);
    {
        // fvm:: matrices are VOLUME-INTEGRATED; the residual is per unit
        // volume. Convert with 1/V so Pmat approximates dF/dx and not
        // V*dF/dx, which would be wrong by a per-row factor on a graded mesh.
        scalarField rV(nCellsLocal);
        {
            const scalarField& Vc = mesh_.V().field();
            forAll(rV, c) { rV[c] = 1.0/Vc[c]; }
        }

        // --- Poisson diagonal block (field 0)
        {
            const dimensionedScalar& epsilon = em.epsilon();
            const scalar sc = sX[0]/sF[0];

            if (em.PoissonScheme() == "semiImplicit")
            {
                const dimensionedScalar dt = mesh_.time().deltaT();
                const volScalarField effEps
                (
                    epsilon + dt*transport.electricalConductivity()
                );
                fvScalarMatrix pEqn
                (
                    fvm::laplacian
                    (
                        effEps, dePotential,
                        "laplacian((epsilon+(deltaT*electricalConductivity)),ePotential)"
                    )
                );
                coo.addFvMatrix(0, pEqn, sc, rV);
            }
            else
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian
                    (
                        epsilon, dePotential, "laplacian(epsilon,ePotential)"
                    )
                );
                coo.addFvMatrix(0, pEqn, sc, rV);
            }
        }

            // --- THE COUPLING d(Poisson residual)/d(n_i), which is why any of
        // this exists. SIGN DERIVED FROM OUR OWN RESIDUAL, not copied: the
        // Poisson residual is F_0 = lap + chargeDensity + ... and
        // chargeDensity = sum_i n_i * speciesCharges_[i]
        // (plasmaSpecies::updateChargeDensity), so dF_0/dn_i = +speciesCharge(i).
        // It is DIAGONAL because chargeDensity is a pointwise sum, and only
        // CHARGED species contribute.
        for (const label id : species.chargedSpeciesIDs())
        {
            const scalar q = species.speciesCharge(id).value();
            const scalarField coeff(nCellsLocal, q);
            coo.addDiagonalBlock(0, 1 + id, coeff, sX[1 + id]/sF[0]);
        }

            // --- Species diagonal blocks, mirroring pcApplyCallback's operators
        for (label s = 0; s < nSpecies; ++s)
        {
            const plasmaTransportModel& model = transport.transportModel(s);
            const scalar Z = species.speciesChargeNumber(s);
            const word sName = species.speciesNames()[s];
            volScalarField& dn = dSpecies[s];

            tmp<surfaceScalarField> tPhi
            (
                Z*fvc::interpolate(model.mu())*em.phiE()
            );
            volScalarField chemLField
            (
                IOobject("chemLcoo_" + sName, mesh_.time().timeName(), mesh_,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh_, dimensionedScalar(dimless/dimTime, Zero)
            );
            // chemP_/chemL_ are SIZED by the chemistry timestep preamble,
            // which runs inside residualCallback -- i.e. AFTER this point on
            // the first outer step of a run. Reading them unsized overran a
            // List and segfaulted (found with gdb, 2026-09-09). Omitting the
            // loss term from the PRECONDITIONER until it exists is safe: a
            // preconditioner need only approximate, and it is refreshed on
            // every subsequent step.
            if (transport.chemL(s).size() == nCellsLocal)
            {
                chemLField.primitiveFieldRef() = transport.chemL(s);
            }

            fvScalarMatrix sEqn
            (
                fvm::ddt(dn)
              + fvm::div(tPhi(), dn, "div(phi_" + sName + ",n_" + sName + ")")
              - fvm::laplacian
                (
                    model.D(), dn,
                    "laplacian(D_" + sName + ",n_" + sName + ")"
                )
              + fvm::Sp(chemLField, dn)
            );
            coo.addFvMatrix(1 + s, sEqn, sX[1 + s]/sF[1 + s], rV);
        }

            // --- Energy diagonal block
        if (lmea)
        {
            const scalar Ze =
                species.speciesChargeNumber(species.electronSpeciesID());
            volScalarField& dE = *dEnergy;

            tmp<surfaceScalarField> tPhiEps
            (
                Ze*fvc::interpolate(lmea->muEpsEff(), "interpolate(mu_e)")
               *em.phiE()
            );
            fvScalarMatrix eEqnP
            (
                fvm::ddt(dE)
              + fvm::div(tPhiEps(), dE, "div(phi_e,n_e)")
              - fvm::laplacian(lmea->DEpsEff(), dE, "laplacian(D_e,n_e)")
              + fvm::Sp(lmea->Lsrc(), dE)
            );
            const label fE = 1 + nSpecies;
            coo.addFvMatrix(fE, eEqnP, sX[fE]/sF[fE], rV);
        }
    }

    SnesPmatCOO pmatCOO;
    pmatCOO.n = int(coo.nEntries());
    pmatCOO.rows = coo.rows();
    pmatCOO.cols = coo.cols();
    pmatCOO.vals = coo.vals();

    int its = 0;
    const int reason = solveWithSNES
    (
        int(nLocalTotal),
        xBuf.data(),
        &residualCallback,
        &ctx,
        &pcApplyCallback,
        &pcCtx,
        bounded_ ? lowerBounds.cdata() : nullptr,
        &pmatCOO,
        &its
    );

    // Copy the converged (or not) state back -- ALWAYS, even on failure,
    // so the caller's own diagnostics see what SNES actually left behind
    // rather than a silently-unwound trial state.
    {
        scalarField& ePotF = em.ePotentialRef().primitiveFieldRef();
        for (label c = 0; c < nCellsLocal; ++c) { ePotF[c] = xBuf[c]*sX[0]; }
        em.ePotentialRef().correctBoundaryConditions();

        for (label s = 0; s < nSpecies; ++s)
        {
            scalarField& nf = species.numberDensity(s).primitiveFieldRef();
            const label off = (1 + s)*nCellsLocal;
            const scalar sxS = sX[1 + s];
            for (label c = 0; c < nCellsLocal; ++c) { nf[c] = xBuf[off + c]*sxS; }
            species.numberDensity(s).correctBoundaryConditions();
        }

        if (lmea)
        {
            scalarField& nef = lmea->nEpsRef().primitiveFieldRef();
            const label off = (1 + nSpecies)*nCellsLocal;
            const scalar sxE = sX[1 + nSpecies];
            for (label c = 0; c < nCellsLocal; ++c) { nef[c] = xBuf[off + c]*sxE; }
            lmea->nEpsRef().correctBoundaryConditions();
        }
    }
    srp.updateDerivedFields();

    // Required bookkeeping the Picard sequence does every corrector that
    // is NOT itself "the solve" -- omitting these left densities at their
    // raw (possibly zero/uniform) initial values with no floor applied,
    // and meanE_/T_ stale at their PRE-solve value (found via all-zero
    // diagnostics after the first successful run, 2026-09-09; see
    // doc/newton-outer-solver-design.md):
    //   - species_.clampNumberDensities() (plasmaTransport.C:1309, inside
    //     solve(), which this path never calls) applies minNumberDensity.
    //   - energyModels_[i].correct() a SECOND time (plasmaEnergy.C:171,
    //     "Re-derive AFTER the solve") re-derives meanE_/T_ from the
    //     CONVERGED nEps -- the step-5 correct() during the residual
    //     evaluation used the LAST TRIAL state, not necessarily the final
    //     converged one.
    species.clampNumberDensities();
    if (lmea) { lmea->correct(); }

    Info<< "outerSolver newton (SNES): reason " << reason
        << " (positive = converged), iterations " << its << endl;

    if (reason <= 0)
    {
        FatalErrorInFunction
            << "SNES did not converge (reason " << reason << ")."
            << nl << exit(FatalError);
    }
}

// ************************************************************************* //
