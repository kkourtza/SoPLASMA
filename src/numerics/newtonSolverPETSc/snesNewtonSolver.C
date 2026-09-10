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
#include "ddWallFluxMixedFvPatchScalarField.H"
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
// For the d(Psrc)/d(phi) convection block: built directly rather than through
// fvm::div so it needs no per-case fvSchemes entry and cannot inherit a
// LIMITED scheme. See jouleJacobian_ in solveOuterStep().
#include "gaussConvectionScheme.H"
#include "upwind.H"

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

    // PER-CELL scales, one array per block. ALWAYS used by the residual and
    // state paths; filled with the uniform sX/sF value when perCellScaling is
    // off, which is what makes "switch off changes nothing" bit-exact rather
    // than merely intended. See perCellScaling_ in the header.
    const Foam::List<Foam::scalarField>* sXc;
    const Foam::List<Foam::scalarField>* sFc;
};

// Block layout, one instance per rank: [ePotential][n_0]..[n_{nSpecies-1}]
// [nEps_e if hasEnergy] -- flat concatenation, no VecNest, exactly the
// pattern the 2-field proof of concept validated in both serial and
// parallel (docs/design/newton-outer-solver-design.md).
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
        const scalarField& sxSc = (*ctx.sXc)[1 + s];
        for (label c = 0; c < nc; ++c) { f[c] = x[off + c]*sxSc[c]; }
    }

    if (ctx.hasEnergy)
    {
        volScalarField& nEps = ctx.lmea->nEpsRef();
        scalarField& f = nEps.primitiveFieldRef();
        const label off = (1 + ctx.nSpecies)*nc;
        const scalarField& sxEc = (*ctx.sXc)[1 + ctx.nSpecies];
        for (label c = 0; c < nc; ++c) { f[c] = x[off + c]*sxEc[c]; }
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
    // chemP_/chemL_ ARE NOT SIZED ON EVERY CHEMISTRY PATH, and reading them
    // unsized is a SEGV inside PETSc's residual evaluation with no usable
    // stack -- which is exactly how it presented (2026-09-10, on
    // positiveStreamer_fixedMesh). plasmaTransport's own comment states the
    // rule: computeChemistrySources() runs only for `solver ode`, while
    // `adaptive`, `adaptiveError` and `implicitRate` go through
    // mechanismSourceTerms(); `explicitSource` does NEITHER and adds its
    // source straight to the matrix, so there is nothing here to read.
    //
    // A guard that merely SKIPPED the source would be worse than the crash:
    // it would silently solve a streamer with no chemistry at all. So this
    // refuses, and names the fix.
    // Either record will do: chemP/chemL from `ode`/`adaptive`/
    // `adaptiveError`/`implicitRate`, or the net source that
    // `explicitSource` now records (with no loss coefficient, because it
    // has none). Both reproduce the term the Picard assembly applies.
    const bool havePL  = transport.chemistrySourcesAvailable();
    const bool haveNet = transport.chemNetSourceAvailable();
    if (!havePL && !haveNet)
    {
        FatalErrorInFunction
            << "outerSolver newton: the chemistry production/loss fields"
            << " chemP/chemL are not available for this case." << nl
            << "    This solver builds its own species residual and needs"
            << " them; `ode`, `adaptive`, `adaptiveError` and `implicitRate`"
            << " populate chemP/chemL, and `explicitSource` records a net"
            << " source instead. This case provided NEITHER." << nl
            << "    WHAT TO DO: use one of those chemistry solvers, or"
            << " outerSolver picard." << nl
            << exit(FatalError);
    }

    const scalarField zeroLoss(nc, Zero);

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
        // L is identically zero on the explicitSource path -- that path has
        // no implicit loss coefficient, so -P alone IS its whole source.
        // zeroLoss is SIZED (hoisted above the loop): a default-constructed
        // scalarField is EMPTY, and Lr[c] on it is the same out-of-bounds
        // read that made this callback segfault in the first place.
        const scalarField& P =
            havePL ? transport.chemP(s) : transport.chemNetSource(s);
        const scalarField& Lr =
            havePL ? transport.chemL(s) : zeroLoss;
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

    // Per-cell state scales, same arrays the residual callback uses. The PC
    // maps physical -> petsc, so it must apply the SAME scaling the
    // residual's inverse mapping does or the preconditioner would solve a
    // differently-scaled operator than the one being preconditioned.
    const Foam::List<Foam::scalarField>* sXc;
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
    // PCFieldSplitRestrictIS. See docs/design/newton-outer-solver-design.md.
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
                << " vector. See docs/design/newton-outer-solver-design.md." << endl;
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
        const scalarField& sxS = (*ctx.sXc)[1 + s];
        for (label c = 0; c < nc; ++c) { y[off + c] = f[c]/sxS[c]; }
    }
    if (ctx.hasEnergy)
    {
        const scalarField& f = ctx.dEnergy->primitiveField();
        const label off = (1 + ctx.nSpecies)*nc;
        const scalarField& sxE = (*ctx.sXc)[1 + ctx.nSpecies];
        for (label c = 0; c < nc; ++c) { y[off + c] = f[c]/sxE[c]; }
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
    bounded_(dict.getOrDefault<bool>("bounded", false)),
    petscOptions_(dict.getOrDefault<string>("petscOptions", string::null)),
    rebalanceScales_(dict.getOrDefault<bool>("rebalanceScales", true)),
    chemJacobian_(dict.getOrDefault<bool>("chemJacobian", true)),
    chemCrossJacobian_(dict.getOrDefault<bool>("chemCrossJacobian", false)),
    jouleJacobian_(dict.getOrDefault<bool>("jouleJacobian", false)),
    schurOnPhi_(dict.getOrDefault<bool>("schurOnPhi", false)),
    perCellScaling_(dict.getOrDefault<bool>("perCellScaling", false)),
    perCellScaleFloor_
    (
        dict.getOrDefault<scalar>("perCellScaleFloor", 1e-6)
    ),
    jouleJacobianRatioMax_
    (
        dict.getOrDefault<scalar>("jouleJacobianRatioMax", 10.0)
    ),
    adaptiveForcing_(dict.getOrDefault<bool>("adaptiveForcing", true)),
    extrapolateGuess_(dict.getOrDefault<bool>("extrapolateGuess", false)),
    extrapolatePotential_(dict.getOrDefault<bool>("extrapolatePotential", true))
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
            << " supported; see docs/design/newton-outer-solver-design.md."
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
    // clamp is needed -- see docs/design/newton-outer-solver-design.md.
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
                << "    See docs/design/newton-outer-solver-design.md (defect A)."
                << nl << exit(FatalError);
        }
    }

    for (label s = 0; s < nSpecies; ++s)
    {
        // The requirement is COEFFICIENTS, not a specific model type: this
        // solver assembles the species equation itself (fvc::flux + the
        // case's div/laplacian schemes), so any model that can hand over
        // mu() and D() works. `immobile` qualifies by returning ZERO for
        // both, which reduces the same assembly to ddt(n) == sources with
        // no branching here -- see immobile.H.
        //
        // Was `isA<driftDiffusion>` until 2026-09-10, which refused the
        // positiveStreamer benchmark outright because it PINS
        // `ionTransport immobile`.
        if (!transport.transportModel(s).providesTransportCoefficients())
        {
            FatalErrorInFunction
                << "outerSolver newton (type SNES) needs a transport model"
                << " that exposes mu() and D() -- species '"
                << species.speciesNames()[s] << "' uses '"
                << transport.transportModel(s).modelName()
                << "', which does not." << nl
                << "    Supported today: driftDiffusion, immobile." << nl
                << "    See docs/design/newton-outer-solver-design.md."
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
                << "    See docs/design/newton-outer-solver-design.md."
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

    // ---- PER-CELL SCALES. Declared here because the state PACK below already
    // needs them, and sX is final at this point. sF is NOT -- rebalanceScales_
    // rescales it further down -- so sFcell is refreshed after that pass.
    //
    // When OFF these hold the uniform sX/sF value in every cell, so the
    // residual and state arithmetic is byte-for-byte what the scalar path did.
    // That is the regression gate; see perCellScaling_ in the header.
    //
    // When ON, a cell's state scale is its own magnitude, floored at a fraction
    // of the block maximum: xhat ~ 1 everywhere, so a cell orders below the
    // peak is no longer invisible to a Krylov norm. The floor is not optional --
    // a cell AT the density floor would otherwise get scale ~0 and its scaled
    // residual would explode, the same pathology reversed.
    //
    // phi keeps a UNIFORM scale even when this is on: it is O(100 V) across the
    // whole gap and has no spread problem. Only the transported blocks, whose
    // spread is the reason this exists, go per-cell.
    List<scalarField> sXcell(nFields);
    List<scalarField> sFcell(nFields);
    for (label b = 0; b < nFields; ++b)
    {
        sXcell[b].setSize(nCellsLocal, sX[b]);
        sFcell[b].setSize(nCellsLocal, sF[b]);
    }
    if (perCellScaling_)
    {
        auto fillFrom = [&](const label b, const scalarField& v)
        {
            const scalar mx = max(gMax(mag(v)), VSMALL);
            const scalar flo = perCellScaleFloor_*mx;
            forAll(sXcell[b], c) { sXcell[b][c] = max(mag(v[c]), flo); }
        };
        for (label sp = 0; sp < nSpecies; ++sp)
        {
            fillFrom(1 + sp, species.numberDensity(sp).primitiveField());
        }
        if (lmea) { fillFrom(1 + nSpecies, lmea->nEps().primitiveField()); }
    }

    ResidualContext ctx
    {
        &em, &srp, &species, &transport, lmea,
        nCellsLocal, nSpecies, lmea != nullptr,
        &sX, &sF,
        &sXcell, &sFcell
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
        &sX, &sF,
        &sXcell
    };

    // Packed in SCALED state units -- everything PETSc sees is scaled.
    List<double> xBuf(nLocalTotal, Zero);
    {
        const scalarField& ePotF = em.ePotential().primitiveField();
        // EXTRAPOLATION RATIO. deltaT0() is the PREVIOUS step, so this is
        // (dt/dt_old) and it is 1 for a constant step. Guarded: on the first
        // step there is no meaningful history and r collapses to 0, which
        // recovers the old zeroth-order guess exactly.
        const scalar dtNow = mesh_.time().deltaTValue();
        const scalar dtOld = mesh_.time().deltaT0Value();
        const scalar r =
            (extrapolateGuess_ && dtOld > VSMALL && mesh_.time().timeIndex() > 2)
          ? dtNow/dtOld
          : 0.0;

        {
            // The potential has no ddt, so its extrapolation is a separate
            // question from the transported fields' -- see
            // extrapolatePotential_.
            const scalar rPhi = (extrapolatePotential_ ? r : 0.0);
            const scalarField& p0 = em.ePotential().oldTime().primitiveField();
            for (label c = 0; c < nCellsLocal; ++c)
            {
                xBuf[c] = (ePotF[c] + rPhi*(ePotF[c] - p0[c]))/sX[0];
            }
        }

        for (label s = 0; s < nSpecies; ++s)
        {
            const volScalarField& ns = species.numberDensity(s);
            const scalarField& nf = ns.primitiveField();
            const scalarField& n0 = ns.oldTime().primitiveField();
            const label off = (1 + s)*nCellsLocal;
            const scalarField& sxS = sXcell[1 + s];
            // FLOOR the extrapolation: an undershoot to a negative density
            // breaks the rate-table lookups the residual performs.
            const scalar fl = max(species.speciesMinNumberDensity(s), scalar(0));
            for (label c = 0; c < nCellsLocal; ++c)
            {
                const scalar e = nf[c] + r*(nf[c] - n0[c]);
                xBuf[off + c] = max(e, fl)/sxS[c];
            }
        }

        if (lmea)
        {
            const volScalarField& ne = lmea->nEps();
            const scalarField& nef = ne.primitiveField();
            const scalarField& ne0 = ne.oldTime().primitiveField();
            const label off = (1 + nSpecies)*nCellsLocal;
            const scalarField& sxE = sXcell[1 + nSpecies];
            for (label c = 0; c < nCellsLocal; ++c)
            {
                const scalar e = nef[c] + r*(nef[c] - ne0[c]);
                xBuf[off + c] = max(e, scalar(0))/sxE[c];
            }
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
      // FGMRES, NOT GMRES, and the restart raised to match -ksp_max_it so no
      // restart happens inside one solve.
      //
      // THE PRECONDITIONER VARIES BETWEEN KRYLOV ITERATIONS: PCFIELDSPLIT with
      // a Schur complement solves its blocks ITERATIVELY, so its action is not
      // a fixed linear operator. Standard GMRES assumes a constant right
      // preconditioner and loses the Arnoldi relation when that assumption
      // fails; FGMRES stores the preconditioned vectors precisely so a varying
      // PC is admissible. Knoll & Keyes (JCP 193 (2004) 357-397) section 3.5
      // makes exactly this recommendation for JFNK with a variable PC.
      //
      // MEASURED, 2026-09-10, and this is why it is not a stylistic choice:
      // with adaptive dt enabled on grubert2009, the Krylov solve failed with
      // DIVERGED_BREAKDOWN at ITERATION 30 -- GMRES's default restart length --
      // once dt had grown to ~6.5e-10 (650x the fixed step Phase D used). The
      // breakdown is at the restart boundary, which is the classic signature of
      // a varying PC breaking GMRES's assumption rather than of a bad matrix:
      // no NaN was present, and the SNES norms were falling cleanly (3.5e2 ->
      // 1.3e-9 in 5 iterations) on the steps before it.
      + " -ksp_type fgmres -ksp_gmres_restart 100"
      // Inexact Newton: see adaptiveForcing_ for the measurement.
      + (adaptiveForcing_ ? " -snes_ksp_ew" : word(""))
      + " -mat_mffd_type wp"
      // FIELDSPLIT defaults, baked in rather than left to the environment so
      // a long unattended run is reproducible. Schur is the only fieldsplit
      // type that uses the off-diagonal species->Poisson block, which is the
      // entire reason the Pmat is assembled; `full` matched an exact LU of
      // Pmat in the survey's measurements (7 KSP iterations against 36 for
      // the block Gauss-Seidel shell). A case may still override any of these.
      + " -pc_fieldsplit_schur_fact_type full"
      // THE INNER SOLVES, which were left at PETSc's defaults until 2026-09-10
      // and were costing an order of magnitude. Found with -ksp_view, which is
      // the only way to see them -- nothing in this solver's own output reports
      // the inner KSP configuration.
      //
      // The A00 (phi) solve is EXACT. Not an optimisation: the Schur complement
      // applies an inner phi solve every time it is applied, so an ITERATIVE
      // inner solve taking a varying number of iterations (measured: 55, 65,
      // 71, 71, 82, 85 on successive applications) makes S's action a
      // NON-FIXED operator, and the outer Krylov method then breaks down.
      // That is what crashed grubert_steady with a SIGFPE -- see
      // docs/design/newton-ignition-experiments.md section 25c. LU on the phi block is
      // affordable precisely because it is ONE field (2000 rows, nnz 9190,
      // factor fill 3.1 on the reference case).
      // SCALABLE, not exact -- rule 43. An earlier version of this default was
      // `preonly` + `lu`, justified as affordable at 2000 rows. That is true at
      // 2000 rows and useless at 1e6: sparse LU is ~O(n^1.5) in 2-D and worse
      // in 3-D, and this solver exists for large 2-D and eventually 3-D cases.
      //
      // The REQUIREMENT is that this solve be a FIXED LINEAR OPERATOR, because
      // the Schur complement applies it on every application and a varying one
      // makes S non-fixed, which breaks the outer Krylov method (that is what
      // raised the SIGFPE of section 25c). FIXEDNESS IS A PROPERTY OF THE
      // ITERATION, NOT OF EXACTNESS: a FIXED number of AMG cycles with a frozen
      // setup is both a fixed linear operator and O(n).
      //
      // `richardson` with max_it 2 and the convergence test SKIPPED, so it
      // never exits early and every application performs identical arithmetic.
      // NOTE that gmres/cg CANNOT be used here at all, however tempting: their
      // Krylov polynomial depends on the right-hand side, so they are not
      // linear in b. Richardson and Chebyshev are.
      //
      // A small case that wants the exact reference can still ask for it --
      // case petscOptions are applied LAST and win:
      //     petscOptions "-fieldsplit_phi_ksp_type preonly -fieldsplit_phi_pc_type lu";
      + " -fieldsplit_phi_ksp_type richardson"
        " -fieldsplit_phi_ksp_max_it 2"
        " -fieldsplit_phi_ksp_convergence_test skip"
        " -fieldsplit_phi_pc_type hypre"
      // The Schur (transport) solve: FGMRES for the same varying-PC reason as
      // the outer level, and three limits that PETSc's defaults get wrong for
      // this problem.
      //
      // rtol 1e-2, NOT the default 1e-5: the OUTER solve runs at rtol 0.3
      // under Eisenstat-Walker, so an inner solve driven to 1e-5 is four
      // orders tighter than anything the outer iteration can use. FGMRES
      // exists to tolerate a variable-quality preconditioner; this is exactly
      // the licence it grants.
      //
      // max_it 200, NOT the default 10000: a solve needing more than this is
      // better handed to retryStep than ground out. The default let 377 solves
      // reach 10000 iterations, wasting 3.8 million of them.
      //
      // restart 100 to match the outer -ksp_gmres_restart above; the default
      // 30 meant 333 restarts in a capped solve, which is how GMRES stagnates.
      //
      // MEASURED 2026-09-10 against the PETSc defaults, same case, same
      // window, normalised per unit SIMULATED time (not per step):
      //   inner iterations/solve  median 230 -> 53, mean 260 -> 59
      //   worst single solve      7226 -> 166 iterations
      //   solves hitting the cap  377 -> 2
      //   total inner Krylov work 9.3x LESS
      + " -fieldsplit_transport_ksp_type fgmres"
        " -fieldsplit_transport_ksp_rtol 1e-2"
        " -fieldsplit_transport_ksp_max_it 200"
        " -fieldsplit_transport_ksp_gmres_restart 100"
      // Differencing step: PETSc's own default unless a case overrides it.
      // It is deliberately NOT set to a large value by default -- see
      // mffdErr_'s declaration for why that would have papered over a real
      // ordering bug rather than accommodating genuine noise in F.
      + (mffdErr_ > 0 ? " -mat_mffd_err " + Foam::name(mffdErr_) : word(""));

    // THE CASE'S OWN OPTIONS GO LAST, so they WIN: PetscOptionsInsertString
    // applies settings in order and a later one replaces an earlier one. That
    // ordering is the whole point -- the defaults above are a starting point a
    // case can tune from, not a ceiling it has to rebuild the library to pass.
    const string allPetscOptions =
        petscOptions_.empty()
      ? string(petscOptions)
      : string(petscOptions) + " " + petscOptions_;

    if (!petscOptions_.empty())
    {
        static bool reported = false;
        if (!reported)
        {
            reported = true;
            Info<< "outerSolver newton: case PETSc options appended (they"
                << " override the defaults): " << petscOptions_ << endl;
        }
    }

    setPetscOptions(allPetscOptions.c_str());

    // TRIAL-ITERATE TOLERANCE, for the duration of this solve only.
    //
    // Newton evaluates F at iterates the line search and the matrix-free
    // Jacobian propose, and some of those are unphysical. A boundary condition
    // that raises FatalError at such a point kills the run over a state that
    // was about to be REJECTED. Measured 2026-09-10: grubert2009_ballast400
    // died with `wall-flux condition on patch anode has become singular` INSIDE
    // an unfinished SNES solve, 828 residual evaluations into the step.
    //
    // Restored below, so an ACCEPTED state still gets the hard, well-diagnosed
    // error -- there the failure is real physics and the diagnosis is worth
    // having.
    const bool tolSaved = ddWallFluxMixedFvPatchScalarField::tolerateSingular_;
    ddWallFluxMixedFvPatchScalarField::tolerateSingular_ = true;

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
            const scalarField& sxS = sXcell[1 + s];
            for (label c = 0; c < nCellsLocal; ++c)
            {
                // INSIDE the loop now: the floor is a PHYSICAL density, so it
                // must be divided by the same per-cell scale the state is, and
                // that scale differs per cell. Hoisting it out (as this did
                // when the scale was one constant per block) would apply one
                // cell's bound to every cell.
                lowerBounds[off + c] =
                    (floorS > 0 ? double(floorS/sxS[c]) : 0.0);
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

        // PURITY CHECK: is F a FUNCTION OF u AT ALL?
        //
        // A line search cannot increase ||F|| if the direction is a descent
        // direction and F is a genuine function of u. Measured 2026-09-11 on
        // grubert2009_ballast400, the gnorm OSCILLATES through the line search
        // (13.18 -> 13.10 -> 13.23 -> ... -> 16.85 -> ... -> 13.94) and
        // `Cubic step no good` appears 193 times in 52 steps. That is the
        // signature of F not being a pure function of u -- the residual
        // evaluation MUTATES state (updateDerived() writes a clamp back into
        // nEps_; the chemistry caches; the LFA seed), so re-evaluating at the
        // same u gives a DIFFERENT F and Newton is chasing a moving target.
        //
        // This evaluates F twice at the SAME x and reports the difference. A
        // non-zero result is a BUG, not a tolerance: it means the Jacobian
        // PETSc differences is the derivative of something that is not a
        // function.
        {
            List<double> pf2(nLocalTotal, Zero);
            residualCallback(int(nLocalTotal), xBuf.data(), pf2.data(), &ctx);

            // HYSTERESIS: evaluate somewhere ELSE, then come BACK to x. This
            // is what a line search does, and it is the form of impurity that
            // matters -- F can be deterministic for two consecutive calls at
            // the same x and STILL depend on where it was evaluated before,
            // if the evaluation mutates internal state (derived fields,
            // caches, clamps written back into the state).
            {
                List<double> xAway(xBuf);
                for (label i = 0; i < nLocalTotal; ++i)
                {
                    xAway[i] = xBuf[i]*1.01 + 1e-3;   // a real excursion
                }
                List<double> fAway(nLocalTotal, Zero);
                residualCallback(int(nLocalTotal), xAway.data(), fAway.data(), &ctx);

                List<double> pf3(nLocalTotal, Zero);
                residualCallback(int(nLocalTotal), xBuf.data(), pf3.data(), &ctx);

                scalar hsum = 0, hmax = 0; label hworst = -1;
                for (label b = 0; b < nFields; ++b)
                {
                    scalar bd = 0;
                    for (label c = 0; c < nCellsLocal; ++c)
                    {
                        const label i = b*nCellsLocal + c;
                        const scalar d = mag(pf3[i] - primeF[i]);
                        hsum += d*d;
                        bd = max(bd, d);
                    }
                    if (bd > hmax) { hmax = bd; hworst = b; }
                }
                const scalar hG = returnReduce(hsum, sumOp<scalar>());
                const scalar hM = returnReduce(hmax, maxOp<scalar>());

                static label hystCount = 0;
                if (hystCount < 3)
                {
                    ++hystCount;
                    Info<< "  [hysteresis] F(x), F(x'), F(x) again: "
                        << "||dF|| = " << Foam::sqrt(hG)
                        << ", worst block = " << hworst
                        << " (max |dF| = " << hM << ")" << nl
                        << "    NONZERO means an intervening evaluation"
                        << " CHANGED F at the same x -- the line search is"
                        << " chasing a moving target." << endl;
                }
            }

            scalar dsum = 0, fsum = 0, dmax = 0;
            label worstBlock = -1;
            for (label b = 0; b < nFields; ++b)
            {
                scalar bd = 0;
                for (label c = 0; c < nCellsLocal; ++c)
                {
                    const label i = b*nCellsLocal + c;
                    const scalar d = mag(pf2[i] - primeF[i]);
                    dsum += d*d;
                    fsum += sqr(primeF[i]);
                    bd = max(bd, d);
                }
                if (bd > dmax) { dmax = bd; worstBlock = b; }
            }
            const scalar dG = returnReduce(dsum, sumOp<scalar>());
            const scalar fG = returnReduce(fsum, sumOp<scalar>());
            const scalar mG = returnReduce(dmax, maxOp<scalar>());

            static label purityCount = 0;
            if (purityCount < 3)
            {
                ++purityCount;
                Info<< "  [purity] F evaluated TWICE at the same x: "
                    << "||dF|| = " << Foam::sqrt(dG)
                    << ", ||F|| = " << Foam::sqrt(fG)
                    << ", relative = "
                    << Foam::sqrt(dG)/(Foam::sqrt(fG) + VSMALL)
                    << ", worst block = " << worstBlock
                    << " (max |dF| = " << mG << ")" << nl
                    << "    ANY nonzero here means F is NOT a function of u --"
                    << " the matrix-free Jacobian is differencing noise."
                    << endl;
            }
        }

        // REBALANCE sF FROM THE MEASURED RESIDUAL, so every block starts the
        // solve at |F_scaled| ~ 1.
        //
        // `sF = sX/dt` for the transported blocks is an A PRIORI estimate of
        // how big the residual OUGHT to be, and it is not a good one. Measured
        // 2026-09-10 on grubert2009_ballast400 at the pre-ignition state, the
        // per-block scaled residuals at the first Newton call were:
        //
        //     Poisson  44.7    n_e  35.6    n_Ar2p  0.147
        //     n_Arp     0.0077 nEps_e 111.2
        //
        // a spread of 14,000x across blocks that are all supposed to be O(1).
        // A Krylov method minimises the NORM OF THE WHOLE VECTOR, so in that
        // state it is fitting nEps_e and the Poisson block and is nearly blind
        // to n_Arp: the ion block contributes ~1e-4 of the residual norm and
        // therefore almost nothing to the Krylov space, however wrong it is.
        //
        // Normalising by the ACTUAL initial residual costs nothing -- the
        // priming evaluation above already had to happen so that Pmat and F
        // describe the same state -- and makes the blocks commensurate, which
        // is the property the scaling was introduced to provide in the first
        // place (Knoll & Keyes 2.3.1: the `typ u` scaling exists so that no
        // component dominates the norm merely through its units).
        //
        // Done BEFORE the Pmat assembly below, which reads sX/sF to scale its
        // entries, so the matrix and the residual stay consistent.
        if (rebalanceScales_)
        {
            for (label b = 0; b < nFields; ++b)
            {
                const label off = b*nCellsLocal;

                scalar ss = 0;
                for (label c = 0; c < nCellsLocal; ++c)
                {
                    ss += sqr(scalar(primeF[off + c]));
                }

                // Collectives on every rank, same number of times (rule 31).
                const scalar ssG = returnReduce(ss, sumOp<scalar>());
                const label nG = returnReduce(nCellsLocal, sumOp<label>());

                const scalar rms = (nG > 0 ? Foam::sqrt(ssG/nG) : 0);

                // A block whose residual is already zero needs no rebalancing,
                // and dividing by it would be a zero-scale catastrophe.
                if (rms > SMALL)
                {
                    sF[b] *= rms;
                }
            }
        }
    }

    // sF was rescaled by rebalanceScales_ above, so its per-cell copy is stale.
    // The RESIDUAL scale stays per BLOCK even when perCellScaling is on: sF is
    // already normalised by the measured priming residual, and dividing the
    // residual per cell as well would rescale the quantity the convergence test
    // reads, changing what `rtol` means from cell to cell.
    for (label b = 0; b < nFields; ++b)
    {
        sFcell[b] = sF[b];
    }

    // ---- COLUMN SCALES AS volScalarFields, for the per-cell Pmat path.
    //
    // Needed only when perCellScaling is on. A volScalarField rather than a
    // bare array because a processor-interface entry's COLUMN is a cell on
    // another rank: correctBoundaryConditions() exchanges the patch values, so
    // addFvMatrixBlockPerCell can read the neighbour's scale directly. See
    // blockMatrixCOO.H.
    //
    // Only the SPECIES and ENERGY columns are built. phi keeps a uniform scale,
    // so every block whose COLUMN is phi (the drift-coupling blocks) still uses
    // the original scalar path unchanged.
    PtrList<volScalarField> colScale(perCellScaling_ ? nFields : 0);
    if (perCellScaling_)
    {
        for (label b = 1; b < nFields; ++b)
        {
            colScale.set
            (
                b,
                new volScalarField
                (
                    IOobject
                    (
                        "colScale" + Foam::name(b), mesh_.time().timeName(),
                        mesh_, IOobject::NO_READ, IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("one", dimless, 1.0),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
            colScale[b].primitiveFieldRef() = sXcell[b];
            colScale[b].correctBoundaryConditions();
        }
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
            if (perCellScaling_)
            {
                // addDiagonalBlock takes a scalar multiplier, so the per-cell
                // COLUMN scale is folded into the coefficient array instead.
                // This block is cell-local (no off-diagonals), so there is no
                // neighbour column to worry about.
                scalarField cs(coeff);
                forAll(cs, c) { cs[c] *= sXcell[1 + id][c]; }
                coo.addDiagonalBlock(0, 1 + id, cs, 1.0/sF[0]);
            }
            else
            {
                coo.addDiagonalBlock(0, 1 + id, coeff, sX[1 + id]/sF[0]);
            }
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

            // -d(chemP)/d(n), approximated as -P/n. See chemJacobian_.
            // Sized-guarded exactly like chemL above: chemP_ is sized by the
            // chemistry preamble inside residualCallback, so on the very first
            // outer step of a run it may not exist yet.
            volScalarField chemPoverN
            (
                IOobject("chemPoverN_" + sName, mesh_.time().timeName(), mesh_,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh_, dimensionedScalar(dimless/dimTime, Zero)
            );
            if (chemJacobian_ && transport.chemP(s).size() == nCellsLocal)
            {
                const scalarField& P = transport.chemP(s);
                const scalarField& nsf = species.numberDensity(s).primitiveField();
                const scalar nFloor = max(species.speciesMinNumberDensity(s), SMALL);

                scalarField& pn = chemPoverN.primitiveFieldRef();
                forAll(pn, c)
                {
                    // NEGATIVE: F carries -P, so dF/dn gets -dP/dn.
                    pn[c] = -P[c]/max(nsf[c], nFloor);
                }
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
              + fvm::Sp(chemPoverN, dn)
            );
            if (perCellScaling_)
            {
                // rowFactor carries 1/V and 1/sF; the column factor is the
                // species' own per-cell state scale. Reduces to the scalar
                // call below when the scales are uniform -- but NOT bitwise,
                // which is why this is a branch and not a replacement.
                scalarField rf(rV);
                forAll(rf, c) { rf[c] /= sF[1 + s]; }
                coo.addFvMatrixBlockPerCell
                (
                    1 + s, 1 + s, sEqn, rf, colScale[1 + s]
                );
            }
            else
            {
                coo.addFvMatrix(1 + s, sEqn, sX[1 + s]/sF[1 + s], rV);
            }

            // --- THE COUPLING d(species residual)/d(ePotential), WITHOUT
            // WHICH THE SCHUR COMPLEMENT IS INERT.
            //
            // The drift term is div(Z*mu_f*phiE*n) and
            // phiE = -snGrad(ePotential)*magSf
            // (singleRegionPoisson.C:39), so perturbing the potential by
            // d(phi) perturbs the flux by -Z*mu_f*n_f*snGrad(dphi)*magSf, and
            //
            //     d/dphi [ div(Z*mu*phiE*n) ] = -laplacian(Z*mu*n, dphi)
            //
            // because fvm::laplacian(G, psi) IS div(G_f*snGrad(psi)*magSf).
            //
            // WHY THIS IS THE BUG AND NOT AN OMISSION. PCFIELDSPLIT with a
            // Schur complement on splits {phi, transport} forms
            //
            //     S = A_tt - A_tp * A_pp^-1 * A_pt
            //
            // A_pt (d(Poisson)/d(n), the charge density) was assembled above.
            // A_tp is THIS block, and it was absent -- so A_tp = 0, hence
            // S = A_tt EXACTLY and the Schur complement contributed nothing at
            // all. The assembled Pmat, the two ISs and the Schur factorisation
            // were doing the work of a block-triangular preconditioner that
            // knows the densities move the field but not that the field moves
            // the densities.
            //
            // That is invisible while the coupling is weak and fatal once it
            // is not, which is exactly the observed behaviour (measured
            // 2026-09-10, docs/design/newton-ignition-experiments.md): Newton runs
            // pre-ignition and stalls with DIVERGED_ITS at ignition, and NO
            // sub-preconditioner helps -- hypre, bjacobi, selfp, 5x the Krylov
            // budget, Eisenstat-Walker, three line searches and both
            // differencing steps all failed, because the missing physics is
            // not in the matrix for any of them to precondition.
            //
            // Z == 0 for a neutral species: it does not drift, the block is
            // identically zero, and assembling it would only add explicit
            // zeros.
            if (mag(Z) > SMALL)
            {
                volScalarField driftCoeff
                (
                    IOobject
                    (
                        "driftCoeff_" + sName, mesh_.time().timeName(), mesh_,
                        IOobject::NO_READ, IOobject::NO_WRITE
                    ),
                    Z*model.mu()*species.numberDensity(s)
                );

                fvScalarMatrix cEqn
                (
                  - fvm::laplacian
                    (
                        driftCoeff, dePotential,
                        "laplacian(epsilon,ePotential)"
                    )
                );

                coo.addFvMatrixBlock
                (
                    1 + s, 0, cEqn, sX[0]/sF[1 + s], rV
                );
            }

            // --- d(species residual)/d(n_e), the CROSS-SPECIES chemistry
            // coupling. Diagonal in cells (a source term is pointwise), so it
            // uses addDiagonalBlock rather than an fvMatrix. See
            // chemCrossJacobian_ for the approximation and its limits.
            if
            (
                chemCrossJacobian_
             && s != species.electronSpeciesID()
             && transport.chemP(s).size() == nCellsLocal
            )
            {
                const label eID = species.electronSpeciesID();
                const scalarField& P = transport.chemP(s);
                const scalarField& ne =
                    species.numberDensity(eID).primitiveField();
                const scalar neFloor =
                    max(species.speciesMinNumberDensity(eID), SMALL);

                scalarField coeff(nCellsLocal);
                forAll(coeff, c)
                {
                    // F carries -P, so dF/dn_e gets -dP/dn_e.
                    coeff[c] = -P[c]/max(ne[c], neFloor);
                }

                coo.addDiagonalBlock
                (
                    1 + s, 1 + eID, coeff, sX[1 + eID]/sF[1 + s]
                );
            }
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
            if (perCellScaling_)
            {
                scalarField rf(rV);
                forAll(rf, c) { rf[c] /= sF[fE]; }
                coo.addFvMatrixBlockPerCell(fE, fE, eEqnP, rf, colScale[fE]);
            }
            else
            {
                coo.addFvMatrix(fE, eEqnP, sX[fE]/sF[fE], rV);
            }

            // --- d(energy residual)/d(ePotential), the same missing coupling
            // as the species block above and for the same reason: the energy
            // drift flux is div(Ze*muEps_f*phiE*nEps) and
            // phiE = -snGrad(ePotential)*magSf, so
            //
            //     d/dphi [ div(Ze*muEps*phiE*nEps) ]
            //         = -laplacian(Ze*muEps*nEps, dphi)
            //
            // The (phi, energy) block is NOT its mirror and is deliberately
            // absent: the Poisson residual depends on the species only through
            // chargeDensity = sum_i q_i*n_i, and nEps carries no charge, so
            // d(F_phi)/d(nEps) is identically zero. The coupling here is
            // genuinely one-way.
            //
            // d(Psrc)/d(phi) -- JOULE HEATING's response to the field -- is
            // added below under jouleJacobian_, no longer missing. See there
            // for the derivation and for why an analytic Joule form is
            // legitimate even under `energySource chemistry`.
            {
                volScalarField driftCoeffE
                (
                    IOobject
                    (
                        "driftCoeff_nEps", mesh_.time().timeName(), mesh_,
                        IOobject::NO_READ, IOobject::NO_WRITE
                    ),
                    Ze*lmea->muEpsEff()*lmea->nEps()
                );

                fvScalarMatrix cEqnE
                (
                  - fvm::laplacian
                    (
                        driftCoeffE, dePotential,
                        "laplacian(epsilon,ePotential)"
                    )
                );

                coo.addFvMatrixBlock(fE, 0, cEqnE, sX[0]/sF[fE], rV);
            }

            // d(Psrc)/d(phi): JOULE HEATING's response to the potential.
            //
            // DERIVATION. updateSources() builds
            //     Psrc = jouleHeating = mu_e * n_e * |E|^2 / (1 V) = C |grad phi|^2,
            //     C = mu_e n_e / (1 V) >= 0
            // and E = -grad(phi), so
            //     d(Psrc) = 2 C grad(phi) . grad(dphi) = -2 C E . grad(dphi).
            // The energy residual carries MINUS Psrc, hence
            //     d(F_eps)/d(phi) = +2 C E . grad(dphi) = a . grad(dphi),
            //     a = 2 C E.
            //
            // a.grad(dphi) is NOT laplacian-shaped -- that is why it needed its
            // own derivation rather than a copy of the drift term above. Using
            // a.grad(u) = div(a u) - u div(a) it becomes a convection operator
            // plus a zeroth-order term, and phiE = E.Sf is exactly the face flux
            // the convection operator wants:
            //     a.Sf = interpolate(2C) * phiE
            //
            // ADAPTIVE TO THE USER'S `energySource` CHOICE, WITH NO BRANCH.
            // Psrc is Joule heating only under `energySource model`; under
            // `energySource chemistry` the per-cell ODE supplies it from
            // chemPeps and no closed form for d(chemPeps)/d(phi) exists. Rather
            // than assume one path, the analytic derivative is scaled by the
            // model's OWN per-cell ratio
            //
            //     s = Psrc / jouleHeating
            //
            // and the field-dependence is taken to remain ~|E|^2, giving
            // d(Psrc)/d(phi) ~ s * d(jouleHeating)/d(phi). This is exact where
            // it can be and corrected where it cannot:
            //
            //   `energySource model`     : Psrc IS jouleHeating, so s == 1
            //                              identically and the derivative is EXACT.
            //   `energySource chemistry` : s carries the ODE's actual departure
            //                              from the bare Joule term, per cell and
            //                              per step, instead of ignoring it.
            //
            // One formula serves both, so nothing has to be kept in sync with
            // the energy model's source logic -- the ratio reads whatever that
            // logic produced, including the optional Coulomb-heating term.
            //
            // s is clipped: jouleHeating -> 0 wherever E -> 0, and there the
            // derivative is negligible anyway, so an unclipped ratio would be
            // pure noise amplification. Outside the clip the term is bounded
            // rather than dropped, which keeps the block from vanishing exactly
            // where a sheath makes it interesting.
            //
            // AND THIS IS THE PRECONDITIONER MATRIX, not the Jacobian. The
            // Jacobian action stays matrix-free (mffd), so Pmat only has to
            // APPROXIMATE dF/dx -- being in the right direction and the right
            // order of magnitude is all a preconditioner needs.
            //
            // SWITCHABLE AND MEASURED, following chemJacobian_/chemCrossJacobian_:
            // the second of those was my own idea and turned out to cost 29%, so
            // a plausible-sounding Jacobian block does not get to be a default on
            // theory alone. See docs/design/newton-ignition-experiments.md.
            if (jouleJacobian_)
            {
                const label eIDj = species.electronSpeciesID();

                volScalarField jouleCoeff
                (
                    IOobject
                    (
                        "jouleCoeff_dPsrc_dphi", mesh_.time().timeName(), mesh_,
                        IOobject::NO_READ, IOobject::NO_WRITE
                    ),
                    2.0*lmea->muEff()*species.numberDensity(eIDj)
                   /dimensionedScalar("oneVolt", dimensionSet(1,2,-3,0,0,-1,0), 1.0)
                );

                // s = Psrc/jouleHeating, per cell: EXACTLY 1 under
                // `energySource model`, the ODE's own departure from the bare
                // Joule term under `energySource chemistry`. See the note above.
                {
                    const scalarField& psrcF = lmea->Psrc().primitiveField();
                    const scalarField& jhF =
                        lmea->jouleHeating().primitiveField();
                    scalarField& jc = jouleCoeff.primitiveFieldRef();

                    // Reference scale for "is jouleHeating meaningful here",
                    // taken from the field itself so it needs no tuning.
                    const scalar jhRef =
                        SMALL*gMax(mag(jhF)) + VSMALL;

                    forAll(jc, c)
                    {
                        const scalar s =
                            (mag(jhF[c]) > jhRef ? psrcF[c]/jhF[c] : 1.0);

                        // Bounded, not dropped: keeps the block alive in a
                        // sheath while refusing a runaway ratio.
                        jc[c] *= min(max(s, 0.0), jouleJacobianRatioMax_);
                    }
                }
                jouleCoeff.correctBoundaryConditions();

                // a.Sf, with phiE = E.Sf already carrying the field's sign.
                //
                // linearInterpolate, NOT fvc::interpolate: this coefficient is
                // internal to the preconditioner and exists in no case
                // dictionary, so fvc::interpolate would demand an
                // `interpolate(jouleCoeff_dPsrc_dphi)` entry in every
                // fvSchemes in existence (these cases list interpolation
                // schemes per field with no `default`, so it is a hard
                // FatalIOError, not a fallback -- observed immediately).
                // linear is also the right choice on its merits: this is a
                // convective flux coefficient, where the harmonic scheme used
                // for DIFFUSIVITIES would be wrong.
                const surfaceScalarField aFlux
                (
                    linearInterpolate(jouleCoeff)*em.phiE()
                );

                // UPWIND, and built directly rather than via fvm::div, for two
                // independent reasons.
                //
                // 1. fvm::div would look up
                //    `div((interpolate(jouleCoeff...)*phiE),d_ePotential)` in
                //    divSchemes. These cases declare div schemes per field with
                //    no `default`, so that is a hard FatalIOError in every case
                //    in the tree (observed). An internal preconditioner operator
                //    must not require a new entry in every user's fvSchemes.
                //
                // 2. The declared drift schemes are all `Gauss ROUNDF` -- BOUNDED
                //    schemes, chosen because a density must stay positive. This
                //    operator acts on a POTENTIAL CORRECTION, which is signed and
                //    has no positivity requirement, and a limited scheme is
                //    state-dependent: its coefficients change with the field it
                //    is applied to. A preconditioner matrix must be a fixed
                //    linear operator -- assembling one that is not is precisely
                //    the defect that caused the SIGFPE this session (section
                //    25c), one level down.
                //
                // Upwind is also the right choice on its merits: it is
                // diagonally dominant, which is what makes a convection block
                // useful to a preconditioner in the first place.
                fvScalarMatrix jEqn
                (
                    fv::gaussConvectionScheme<scalar>
                    (
                        mesh_,
                        aFlux,
                        tmp<surfaceInterpolationScheme<scalar>>
                        (
                            new upwind<scalar>(mesh_, aFlux)
                        )
                    ).fvmDiv(aFlux, dePotential)
                  - fvm::Sp(fvc::div(aFlux), dePotential)
                );

                coo.addFvMatrixBlock(fE, 0, jEqn, sX[0]/sF[fE], rV);
            }
        }
    }

    // NON-FINITE PMAT SCAN, added 2026-09-10.
    //
    // A NaN or inf anywhere in the preconditioner matrix does NOT fail
    // cleanly. It reaches PETSc, and the first MatMult inside
    // PCApply_FieldSplit_Schur raises SIGFPE -- whose backtrace names only
    // PETSc internals (MatMult_SeqAIJ / PCApply / KSPFGMRESCycle) and gives no
    // hint which of the Jacobian blocks assembled above produced it. That is
    // exactly how grubert_steady died at t=1.954e-6, and finding it took a gdb
    // run over a 14 ns restart. One pass over the COO arrays costs nothing
    // measurable against a SNES solve and names the BLOCK, the CELL and the
    // VALUE instead.
    //
    // Deliberately FATAL, not a step rejection: a non-finite Jacobian entry is
    // a defect in the assembly above, not the physical stiffness that
    // retryStep exists to absorb. Masking it as a failed step would hide it.
    {
        const label nEnt = coo.nEntries();
        const label* rr = coo.rows();
        const label* cc = coo.cols();
        const scalar* vv = coo.vals();
        const label lStart = coo.globalRow(0, 0);

        // TWO failure modes, not one. A non-finite entry is the obvious one.
        // The other is an entry that is FINITE BUT ENORMOUS: MatMult computes
        // sum(a_ij*x_j), so it can raise FE_OVERFLOW -- the same SIGFPE, at
        // the same place in the same backtrace -- with every input finite.
        // Checking only isfinite() would report "clean" and leave the crash
        // unexplained. After the sX/sF scaling these entries should be O(1);
        // 1e200 is far past anything a scaled Jacobian can legitimately hold
        // and cannot be multiplied by anything above ~1e108 without
        // overflowing.
        const scalar hugeVal = 1e200;

        label nBad = 0;
        label firstBad = -1;
        label maxAt = -1;
        scalar maxAbs = 0;

        for (label e = 0; e < nEnt; ++e)
        {
            const scalar av = Foam::mag(vv[e]);

            if (std::isfinite(vv[e]) ? (av > hugeVal) : true)
            {
                if (firstBad < 0) { firstBad = e; }
                ++nBad;
            }

            // NaN fails every comparison, so it never becomes the max; that is
            // fine, the non-finite branch above already caught it.
            if (std::isfinite(vv[e]) && av > maxAbs) { maxAbs = av; maxAt = e; }
        }

        const label nBadG = returnReduce(nBad, sumOp<label>());
        const scalar maxAbsG = returnReduce(maxAbs, maxOp<scalar>());

        // Name the block an entry lives in. Row layout is
        // localStart + field*nCells + cell (blockMatrixCOO::globalRow).
        const auto blockName = [&](const label f) -> word
        {
            if (f == 0) { return word("ePotential"); }
            if (f >= 1 && f <= nSpecies)
            {
                return species.speciesNames()[f - 1];
            }
            if (lmea && f == 1 + nSpecies) { return word("nEps_e"); }
            return word("field" + Foam::name(f));
        };

        const auto describe = [&](const label e) -> string
        {
            if (e < 0 || nCellsLocal <= 0) { return string("<none>"); }
            const label rLoc = rr[e] - lStart;
            const label cLoc = cc[e] - lStart;
            return "d(" + blockName(rLoc/nCellsLocal) + ")/d("
                 + blockName(cLoc/nCellsLocal) + ") at local cell "
                 + Foam::name(rLoc % nCellsLocal);
        };

        if (nBadG > 0)
        {
            const bool nonFinite =
                (firstBad >= 0 && !std::isfinite(vv[firstBad]));

            FatalErrorInFunction
                << "The Newton preconditioner matrix contains " << nBadG
                << " unusable entr" << (nBadG == 1 ? "y" : "ies")
                << " (non-finite, or |value| > " << hugeVal << ")." << nl
                << "    First on this rank: " << describe(firstBad).c_str()
                << nl
                << "      value " << (firstBad >= 0 ? vv[firstBad] : 0)
                << "   -- " << (nonFinite ? "NON-FINITE" : "finite but huge")
                << nl
                << "      (row " << (firstBad >= 0 ? rr[firstBad] : -1)
                << ", col " << (firstBad >= 0 ? cc[firstBad] : -1)
                << ", of " << nEnt << " entries; localStart " << lStart
                << ", nCells " << nCellsLocal << ")" << nl
                << "    Largest finite |entry| anywhere: " << maxAbsG
                << "   at " << describe(maxAt).c_str() << nl
                << nl
                << "    This is an assembly defect, not physical stiffness."
                << " Handing it to PETSc" << nl
                << "    lets the first MatMult inside PCApply raise SIGFPE"
                << " -- by overflow if the" << nl
                << "    entries are merely huge -- with a backtrace naming"
                << " only PETSc internals," << nl
                << "    so it is refused here instead." << nl
                << exit(FatalError);
        }

        // The magnitude TREND is the diagnostic, not just the final value: a
        // fatal threshold alone cannot distinguish "my entries were fine, the
        // overflow happened inside PETSc" from "my entries were 1e150 and the
        // threshold was set too high". Printed every outer step so the
        // approach to a crash is visible in the log. One line against the
        // ~550 this solver already writes per step.
        Info<< "    Pmat: " << nEnt << " entries, max|a| " << maxAbsG
            << " at " << describe(maxAt).c_str() << endl;
    }

    // THE SEMI-IMPLICIT POISSON OPERATOR AS THE SCHUR PRECONDITIONER.
    //
    // Assembled only when schurOnPhi_ is set, because it preconditions S_f --
    // the Schur complement on the PHI block -- which exists only when the
    // transport block is the one eliminated (see SnesPmatCOO::schurOnPhi).
    //
    // div((eps + dt*sigma) grad .) IS S_f to leading order in dt. Derivation in
    // docs/design/schur-semiimplicit-poisson-preconditioner.md; in one line,
    //     A_ft inv(A_tt) A_tf ~ q * dt * (-div(Z mu n grad .))
    //                         = -dt div(sigma grad .)
    // because sigma = q Z mu n is exactly the electrical conductivity. So the
    // operator that is INCONSISTENT inside a Newton residual (experiment log
    // section 23, measured 0.3278 ceiling) is precisely the right object here.
    // Wrong place, not wrong idea.
    //
    // SAME scaling (sX[0]/sF[0]) and same 1/V conversion as the Poisson
    // diagonal block above: PETSc uses this matrix IN PLACE OF S_f, so it has
    // to live in the same scaled units as the rest of Pmat or it preconditions
    // a differently-normalised operator.
    blockMatrixCOO schurCoo(nCellsLocal, 1);
    if (schurOnPhi_)
    {
        scalarField rVs(nCellsLocal);
        {
            const scalarField& Vc = mesh_.V().field();
            forAll(rVs, c) { rVs[c] = 1.0/Vc[c]; }
        }

        const volScalarField effEpsSchur
        (
            em.epsilon()
          + mesh_.time().deltaT()*transport.electricalConductivity()
        );

        // The scheme name is the PLAIN one the drift blocks already use, not
        // "laplacian((epsilon+(deltaT*electricalConductivity)),ePotential)":
        // that second name is only declared by cases that actually run the
        // semi-implicit Poisson, and the Newton path forces the scheme to
        // explicit, so it is absent exactly where this is needed.
        fvScalarMatrix sEqn
        (
            fvm::laplacian
            (
                effEpsSchur, dePotential, "laplacian(epsilon,ePotential)"
            )
        );

        schurCoo.addFvMatrix(0, sEqn, sX[0]/sF[0], rVs);
    }

    SnesPmatCOO pmatCOO;
    pmatCOO.n = int(coo.nEntries());
    pmatCOO.rows = coo.rows();
    pmatCOO.cols = coo.cols();
    pmatCOO.vals = coo.vals();
    pmatCOO.schurOnPhi = schurOnPhi_;
    if (schurOnPhi_)
    {
        pmatCOO.nSchur = int(schurCoo.nEntries());
        pmatCOO.schurRows = schurCoo.rows();
        pmatCOO.schurCols = schurCoo.cols();
        pmatCOO.schurVals = schurCoo.vals();
    }

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
        int(nFields),
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
            const scalarField& sxS = sXcell[1 + s];
            for (label c = 0; c < nCellsLocal; ++c) { nf[c] = xBuf[off + c]*sxS[c]; }
            species.numberDensity(s).correctBoundaryConditions();
        }

        if (lmea)
        {
            scalarField& nef = lmea->nEpsRef().primitiveFieldRef();
            const label off = (1 + nSpecies)*nCellsLocal;
            const scalarField& sxE = sXcell[1 + nSpecies];
            for (label c = 0; c < nCellsLocal; ++c) { nef[c] = xBuf[off + c]*sxE[c]; }
            lmea->nEpsRef().correctBoundaryConditions();
        }
    }
    srp.updateDerivedFields();

    // Required bookkeeping the Picard sequence does every corrector that
    // is NOT itself "the solve" -- omitting these left densities at their
    // raw (possibly zero/uniform) initial values with no floor applied,
    // and meanE_/T_ stale at their PRE-solve value (found via all-zero
    // diagnostics after the first successful run, 2026-09-09; see
    // docs/design/newton-outer-solver-design.md):
    //   - species_.clampNumberDensities() (plasmaTransport.C:1309, inside
    //     solve(), which this path never calls) applies minNumberDensity.
    //   - energyModels_[i].correct() a SECOND time (plasmaEnergy.C:171,
    //     "Re-derive AFTER the solve") re-derives meanE_/T_ from the
    //     CONVERGED nEps -- the step-5 correct() during the residual
    //     evaluation used the LAST TRIAL state, not necessarily the final
    //     converged one.
    ddWallFluxMixedFvPatchScalarField::tolerateSingular_ = tolSaved;

    species.clampNumberDensities();
    if (lmea) { lmea->correct(); }

    // REFRESH THE CONVECTIVE FLUXES, or the species Courant limiter is INERT.
    //
    // plasmaTransport populates convectiveFlux_ inside its own solve(), from
    // fvMatrix::flux() -- "Placed AFTER the solve, because fvMatrix::flux() is
    // only defined once the matrix has been solved". The Newton path REPLACES
    // that solve and never calls it, so the field kept its zero-initialised
    // value and `limitSpeciesCo` silently protected nothing.
    //
    // MEASURED 2026-09-11 on grubert2009_ballast400: `Co_conv (e)` reported
    // EXACTLY 0 for entire runs -- a single distinct value across 1300+ steps
    // -- while `Co_conv (energy)`, which the LMEA model populates because
    // Newton does call lmea->correct(), read 6-27 on the same steps. That
    // asymmetry is what exposed it. Harmless on a case whose cap is 1500, and
    // a silent loss of protection on any case that sets a real one.
    //
    // Computed from STATE, not from a matrix: the drift face flux IS
    // phi_s*n_f with phi_s = Z*mu_f*phiE, which is exactly what the Courant
    // number needs (convRate ~ 0.5*sum|phi_s|/V, i.e. v/dx). No solved matrix
    // is required, so this does not resurrect the dependency that put the
    // original call inside solve().
    for (const label i : species.mobileSpeciesIDs())
    {
        const scalar Z = species.speciesChargeNumber(i);
        const plasmaTransportModel& m = transport.transportModel(i);

        transport.convectiveFlux(i) =
            Z*fvc::interpolate(m.mu())*em.phiE()
           *fvc::interpolate(species.numberDensity(i));
    }

    Info<< "outerSolver newton (SNES): reason " << reason
        << " (positive = converged), iterations " << its << endl;

    // REPORTED, NOT FATAL. A step Newton cannot solve is a step that was too
    // large, and the response to that is a shorter step -- exactly what Picard
    // already does under `outerCoupling/onNonConvergence retryStep`. Killing
    // the run instead denied Newton the one thing it needs to find its own
    // stability ceiling, and made the Picard->Newton handover fatal by
    // construction: Newton inherits whatever dt the easy Picard phase wound up
    // to and meets it head-on.
    //
    // soPlasmaFoam reads lastSolveConverged() and, when it is false, discards
    // the step and retries at a shorter dt. If the case's policy is `fatal`,
    // plasmaTimeControl raises that itself -- the policy lives in ONE place
    // rather than being duplicated here.
    lastConverged_ = (reason > 0);

    if (!lastConverged_)
    {
        WarningInFunction
            << "SNES did not converge (reason " << reason << ") at dt = "
            << mesh_.time().deltaTValue() << "." << nl
            << "    The step is being reported as non-converged so the outer"
            << " loop can DISCARD and retry it at a shorter dt." << nl
            << "    A persistent failure here means the step is above this"
            << " problem's Newton stability ceiling." << endl;
    }
}

// ************************************************************************* //
