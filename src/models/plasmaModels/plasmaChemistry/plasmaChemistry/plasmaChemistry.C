/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
\*---------------------------------------------------------------------------*/

#include "plasmaChemistry.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

Foam::label Foam::plasmaChemistry::stateIndex(const word& name) const
{
    forAll(species_, i)
    {
        if (species_[i] == name) return i;
    }
    return -1;
}


void Foam::plasmaChemistry::readMechanism
(
    const fileName& path,
    const dictionary& dict
)
{
    IFstream is(path);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open mechanism dictionary " << path << nl
            << exit(FatalError);
    }
    const dictionary mech(is);

    const word electron = mech.getOrDefault<word>("electronSpecies", "Electron");
    const word caseElectron = dict.getOrDefault<word>("electronName", "e");

    // Reference composition, for species the case does not transport. Their
    // densities are parameters, not state.
    if (mech.found("composition"))
    {
        const scalar nGas = dict.getOrDefault<scalar>("backgroundDensity", 0.0);
        for (const entry& e : mech.subDict("composition"))
        {
            background_.set(e.keyword(), readScalar(e.stream())*nGas);
        }
    }

    auto resolve = [&](const word& raw) -> word
    {
        return (raw == electron) ? caseElectron : raw;
    };

    // Cantera-side names of the carried species. The mechanism writes the
    // OpenFOAM-safe spelling everywhere (`N2p`, `O2m`) because a `+` or `-` is
    // not a legal OpenFOAM word; `speciesAlias` maps back to what Cantera
    // parses. Species absent from the table are spelt identically in both.
    canteraNames_ = species_;
    if (mech.found("speciesAlias"))
    {
        for (const entry& e : mech.subDict("speciesAlias"))
        {
            alias_.set(e.keyword(), word(string(e.stream()), false));
        }
        forAll(species_, i)
        {
            canteraNames_[i] = alias_.lookup(species_[i], species_[i]);
        }
    }

    // Fill a reaction's stoichiometry, folding untransported reactants into a
    // fixed density product so the right-hand side does not re-multiply them.
    auto fill = [&](const wordList& names, labelList& out, scalar& fixed)
    {
        DynamicList<label> idx;
        forAll(names, i)
        {
            const word nm = resolve(names[i]);
            const label s = stateIndex(nm);
            idx.append(s);
            if (s < 0)
            {
                // Not carried: use its fixed background density, or zero,
                // which switches the reaction off rather than pretending the
                // species is absent from the rate law.
                fixed *= background_.lookup(nm, 0.0);
            }
        }
        out = labelList(idx);
    };

    // ---- electron-impact processes ----------------------------------------
    if (mech.found("processes"))
    {
        const List<dictionary> procs(mech.lookup("processes"));
        forAll(procs, i)
        {
            const dictionary& p = procs[i];
            plasmaReactionSpec rx;
            rx.tabulated = nTabulated_++;
            tabulatedIds_.append(p.get<word>("id"));

            rx.fixedReactantDensity = 1.0;
            fill(p.get<wordList>("reactants"), rx.reactants,
                 rx.fixedReactantDensity);
            scalar dummy = 1.0;
            fill(p.get<wordList>("products"), rx.products, dummy);

            if (p.found("collider"))
            {
                const word c = resolve(p.get<word>("collider"));
                rx.collider = stateIndex(c);
                if (rx.collider < 0)
                {
                    rx.colliderFixedDensity = background_.lookup(c, 0.0);
                }
            }
            reactions_.append(rx);
        }
    }

    // ---- heavy reactions ---------------------------------------------------
    //
    // Skipped entirely for `reactions electronImpact`, which is what the
    // streamer validation was produced with. Filtering HERE rather than at the
    // source-term stage means every downstream consumer -- the ODE, the
    // Jacobian, the heat release, the diagnostics -- sees one consistent
    // reaction set, instead of each having to remember to exclude the same
    // subset.
    const bool withHeavy =
        dict.getOrDefault<bool>("includeHeavyReactions", true);

    if (mech.found("heavyReactions") && withHeavy)
    {
        const List<dictionary> heavy(mech.lookup("heavyReactions"));
        forAll(heavy, i)
        {
            const dictionary& h = heavy[i];
            plasmaReactionSpec rx;
            rx.A  = h.get<scalar>("A");
            rx.b  = h.get<scalar>("b");
            rx.Ta = h.get<scalar>("Ta");
            rx.deltaH = h.getOrDefault<scalar>("deltaH", 0.0);

            rx.fixedReactantDensity = 1.0;
            fill(h.get<wordList>("reactants"), rx.reactants,
                 rx.fixedReactantDensity);
            scalar dummy = 1.0;
            fill(h.get<wordList>("products"), rx.products, dummy);

            if (h.found("collider"))
            {
                const word c = resolve(h.get<word>("collider"));
                rx.collider = stateIndex(c);
                if (rx.collider < 0)
                {
                    rx.colliderFixedDensity = background_.lookup(c, 0.0);
                }
            }
            reactions_.append(rx);
        }
    }

    if (reactions_.empty())
    {
        FatalErrorInFunction
            << "Mechanism " << path << " yielded no reactions." << nl
            << exit(FatalError);
    }
}


Foam::scalar Foam::plasmaChemistry::checkBackends(const scalar Tgas) const
{
    if (!cantera_) return 0.0;

    // A reference state with every carried species present. The comparison is
    // only meaningful where a reaction actually proceeds, and a state built
    // from the case's own initial condition would leave most of the mechanism
    // at zero rate and therefore untested. Fractions of the background density
    // rather than equal values, so the test spans the range of magnitudes the
    // solver really sees.
    scalar nGas = 0;
    forAllConstIters(background_, it) nGas += it.val();
    if (nGas <= 0) nGas = 2.5e25;

    scalarField n(species_.size(), Zero);
    forAll(n, i) n[i] = nGas*Foam::pow(10.0, -3.0 - 3.0*scalar(i % 4));

    // Native heavy only: the electron-impact rates are tabulated and are not
    // part of what Cantera evaluates, so they must be out of both sides.
    scalarField kZero(kTab_.size(), Zero);
    kTab_ = kZero;

    const bool savedHeavy = ode_->heavy();
    const plasmaChemistryCantera* savedCt = nullptr;
    ode_->setCantera(nullptr);
    ode_->setHeavy(true);
    ode_->setTgas(Tgas);

    scalarField Pn, Ln;
    ode_->productionLoss(n, Pn, Ln);

    ode_->setHeavy(savedHeavy);
    ode_->setCantera(savedHeavy ? savedCt : cantera_.get());

    scalarField Pc, Lc;
    cantera_->productionLoss(n, Tgas, Pc, Lc);

    // Compared on the NET rate, not on P and L separately: the two backends
    // may legitimately split a reversible reaction differently between
    // production and loss while agreeing on dn/dt, which is the only thing
    // that reaches the solution.
    scalar worst = 0.0, scale = 0.0;
    scalarField dn(n.size(), Zero);
    forAll(n, i)
    {
        dn[i] = (Pn[i] - Ln[i]*n[i]) - (Pc[i] - Lc[i]*n[i]);
        scale = max(scale, mag(Pn[i] - Ln[i]*n[i]));
    }
    if (scale <= VSMALL) return 0.0;

    forAll(dn, i) worst = max(worst, mag(dn[i])/scale);
    return worst;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasmaChemistry::plasmaChemistry
(
    const fileName& mechanismDict,
    const wordList& species,
    const scalarField& charge,
    const dictionary& dict
)
:
    species_(species),
    charge_(charge),
    odeDict_(dict.subOrEmptyDict("odeCoeffs"))
{
    readMechanism(mechanismDict, dict);

    kTab_.setSize(max(nTabulated_, 1), Zero);
    ode_.reset(new plasmaChemistryODE(reactions_, kTab_, charge_,
                                      species_.size()));

    // Stiff by default. Plasma rates span twenty orders of magnitude, and an
    // explicit solver will appear to work while quietly demanding timesteps
    // far below the transport limit -- which shows up as a slow run, not as a
    // wrong answer, and is therefore easy to misattribute.
    if (!odeDict_.found("solver"))
    {
        // DEFAULT: rodas23. Why this one, out of the dozen ODESolver offers.
        //
        // Plasma chemistry is stiff by a wide margin -- rate coefficients here
        // span roughly twenty orders of magnitude, from three-body attachment
        // at 1e-43 m^6/s to electron-impact ionisation at 1e-14 m^3/s -- and
        // it is integrated ONCE PER CELL, millions of times per timestep. The
        // choice is therefore governed by three things at once:
        //
        //   STABILITY.  Anything explicit (Euler, RKCK45, RKDP45, RKF45) is
        //     limited by the fastest reaction, not by the transport. It will
        //     appear to work and then demand timesteps orders below the CFD
        //     limit, which presents as a slow run rather than a wrong answer
        //     and is therefore easy to misattribute. Ruled out on principle,
        //     not by benchmark.
        //
        //   COST PER STEP.  rodas23 is a Rosenbrock method: it needs ONE
        //     Jacobian and ONE LU factorisation per step and takes no Newton
        //     iterations at all. There is consequently no nonlinear
        //     convergence failure to handle -- valuable when the same code
        //     runs unattended over a million cells. seulex reaches higher
        //     order by extrapolation, but pays with several factorisations per
        //     step, and the extra order buys little here because the accuracy
        //     that matters is set by the splitting (first order, see
        //     docs/numerics-chemistry-coupling.md), not by the substep.
        //
        //   LINEAR INVARIANTS.  Charge conservation is a linear invariant of
        //     the reaction set. In exact arithmetic any linearly-implicit
        //     method preserves it, because q.f = 0 implies q.J = 0. In
        //     floating point they do not, and they differ: measured on the
        //     same case, the residual of the integrated source was 3.5e-06
        //     with seulex against 5.0e-08 with rodas23, a factor of 70. That
        //     residual feeds the Poisson equation through the net charge, so
        //     it is not a cosmetic difference. (It is projected out
        //     afterwards regardless -- but starting closer is worth having.)
        //
        // rodas23 is also L-stable and stiffly accurate, so it damps the fast
        // transients instead of ringing on them, which is the behaviour wanted
        // when a substep is much longer than the fastest chemical timescale.
        //
        // `seulex` remains a reasonable choice if substep accuracy ever
        // becomes the limiting error; `EulerSI` is there for cheap, low-accuracy
        // work. Both are selectable through odeCoeffs/solver.
        odeDict_.add("solver", word("rodas23"));
    }
    odeDict_.add("absTol", odeDict_.getOrDefault<scalar>("absTol", 1e-6), true);
    odeDict_.add("relTol", odeDict_.getOrDefault<scalar>("relTol", 1e-4), true);

    // SUBSTEP BUDGET PER CELL. OpenFOAM's default is 10000, and exceeding it
    // raises a FatalError -- which plasmaTransport catches (it enables
    // FatalError.throwExceptions() around the cell loop) and counts as an ODE
    // failure, so it is recoverable rather than fatal here.
    //
    // The problem is the COST of reaching that budget. MEASURED 2026-08-21: on
    // a diverged LMEA step, one rank sat in rodas23 for over 18 minutes at
    // full CPU while the other seven waited in an MPI collective. The run
    // looked alive and was making no progress, and because the outer loop
    // never reached its corrector cap, `onNonConvergence reduceDeltaT` was
    // never reached either -- the recovery path was unreachable while a single
    // cell could grind.
    //
    // 2000 keeps ample headroom for genuinely stiff cells (the healthy steps
    // on that case use tens of substeps) while turning a pathological cell
    // into a counted failure in a fraction of the time. Raise it per case
    // through odeCoeffs/maxSteps if a mechanism genuinely needs more.
    odeDict_.add
    (
        "maxSteps", odeDict_.getOrDefault<label>("maxSteps", 2000), true
    );

    solver_ = ODESolver::New(*ode_, odeDict_);

    // ---- heavy-chemistry backend -------------------------------------------
    backend_ = dict.getOrDefault<word>("chemistryBackend", "native");

    if (backend_ == "cantera")
    {
        // Cantera reads <mech>.heavy.yaml, written by mechc beside the .foam
        // dictionary from the same master mechanism. Deriving the name rather
        // than asking for it keeps the two from being paired wrongly: they are
        // two projections of one file and must never be mixed across versions.
        fileName yaml(mechanismDict.lessExt() + ".heavy.yaml");

        cantera_.reset
        (
            new plasmaChemistryCantera(yaml, species_, canteraNames_)
        );

        forAllConstIters(background_, it)
        {
            // ONLY species the case does not transport. background_ holds a
            // reference density for every species in the mechanism's
            // composition, including ones that are carried -- the native path
            // consults it solely when stateIndex() fails, so a carried species
            // never reads it. Pushing those to Cantera as well would add the
            // reference density on top of the transported one and count the
            // bulk gas twice, which is most of the mixture.
            if (stateIndex(it.key()) >= 0) continue;

            // background_ is keyed on the OpenFOAM spelling; translate.
            cantera_->setBackground(alias_.lookup(it.key(), it.key()), it.val());
        }

        // How the heavy Jacobian is obtained. `auto` takes Cantera's analytic
        // composition derivative where the mechanism's kinetics manager
        // implements it, which is the fast and exact path; `finiteDifference`
        // forces the fallback, at nSpecie+1 kinetics evaluations per Jacobian.
        // There is no option to borrow the native Jacobian: it is wrong in
        // exactly the cases Cantera was selected for.
        const word jac(dict.getOrDefault<word>("canteraJacobian", "auto"));
        if (jac == "finiteDifference")
        {
            cantera_->useFiniteDifferenceJacobian();
        }
        else if (jac != "auto")
        {
            FatalErrorInFunction
                << "Unknown canteraJacobian " << jac << nl
                << "    Valid: auto, finiteDifference" << nl
                << exit(FatalError);
        }

        // Heavy reactions are now Cantera's job; evaluating them natively too
        // would double every heavy rate.
        ode_->setHeavy(false);
        ode_->setCantera(cantera_.get());

        const wordList missing(cantera_->unmatched());

        // The two backends are cross-checked on the reactions they share. This
        // is no longer load-bearing for the Jacobian -- Cantera supplies its
        // own analytic composition derivative, so `ode` and `adaptive` are
        // exact under either backend -- but it stays because it is the only
        // thing that turns a silent factor-of-N error into a message at t=0.
        // A units block, a third body counted twice, a mismatched alias and a
        // dropped negative activation energy have all been caught by it.
        //
        // A LARGE value is not necessarily a defect: a mechanism using falloff,
        // PLOG or Chebyshev is one the native parser cannot express, and
        // Cantera's answer is then the correct one. It is reported either way,
        // because "the backends differ and here is by how much" is information
        // the user needs and cannot get anywhere else.
        backendMismatch_ = checkBackends(dict.getOrDefault<scalar>("Tgas", 300.0));

        Info<< "plasmaChemistry: heavy chemistry from Cantera, " << yaml.name()
            << nl
            << "    species Cantera does not carry: "
            << (missing.empty() ? wordList({word("none")}) : missing) << nl
            << "    native-vs-Cantera heavy rate agreement: "
            << backendMismatch_ << " (relative)" << nl
            << "    heavy Jacobian: "
            << (cantera_->analyticJacobian()
                  ? "Cantera analytic (netProductionRates_ddCi)"
                  : "finite-difference fallback -- nSpecie+1 evaluations per"
                    " Jacobian, only in cells the adaptive switch integrates")
            << endl;

        if (backendMismatch_ > 1e-6)
        {
            WarningInFunction
                << "The native and Cantera heavy rates differ by "
                << backendMismatch_ << " at the reference state." << nl
                << "    Expected when the mechanism uses reaction forms the"
                << " native parser cannot represent (falloff, PLOG, Chebyshev,"
                << " reversible): Cantera's rates are then the correct ones and"
                << " this number is the size of what the native path was"
                << " missing." << nl
                << "    Unexpected otherwise -- for a mechanism of plain"
                << " Arrhenius and three-body reactions the two evaluate the"
                << " same expressions and agree to roundoff (~1e-15). Anything"
                << " larger is a defect in one of them." << nl
                << "    Either way both `implicitRate` and `ode` remain valid:"
                << " the Jacobian is Cantera's own where Cantera supplies the"
                << " rates." << endl;
        }
    }
    else if (backend_ != "native")
    {
        FatalErrorInFunction
            << "Unknown chemistryBackend " << backend_ << nl
            << "    Valid: native, cantera" << nl
            << exit(FatalError);
    }

    Info<< "plasmaChemistry: " << reactions_.size() << " reactions ("
        << nTabulated_ << " electron-impact, "
        << reactions_.size() - nTabulated_ << " heavy), "
        << species_.size() << " species, ODE solver "
        << odeDict_.get<word>("solver") << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plasmaChemistry::setBackground(const word& name, const scalar n)
{
    background_.set(name, n);
}


void Foam::plasmaChemistry::enableEnergyEquation
(
    const label electronIndex,
    const scalar epsMin,
    const scalar epsMax,
    const scalar neFloor,
    std::function<scalar(scalar)> muE,
    std::function<scalar(scalar)> PlossN
)
{
    ode_->enableEnergy
    (
        electronIndex, epsMin, epsMax, neFloor,
        std::move(muE), std::move(PlossN)
    );

    // REBUILD. ODESolver allocates its workspace from odes_.nEqns() in its
    // constructor, so a solver built before this call has every internal array
    // one element short -- an out-of-bounds write inside the integrator, which
    // is silent in an optimised build. Rebuilding is cheap and happens once.
    solver_ = ODESolver::New(*ode_, odeDict_);

    Info<< "plasmaChemistry: electron ENERGY added to the integrated state"
        << " (nEqns " << ode_->nEqns() << ", energy index "
        << ode_->energyIndex() << "); ODE solver rebuilt for the new size."
        << endl;
}


void Foam::plasmaChemistry::setEnergyCell
(
    const scalar Emag,
    const scalar Ngas
) const
{
    ode_->setEnergyCell(Emag, Ngas);
}


bool Foam::plasmaChemistry::energyIntegrated() const
{
    return ode_->energyActive();
}


void Foam::plasmaChemistry::integrate
(
    scalarField& n,
    const scalarField& kTab,
    const scalar Tgas,
    const scalar dt,
    const scalarField* ext,
    const scalarField* extSlope
) const
{
    kTab_ = kTab;
    ode_->setTgas(Tgas);
    ode_->setExternal(ext);
    ode_->setExternalSlope(extSlope, dt);

    // The state the transport rate was MEASURED at, so a removing cross term
    // can be carried as a sink proportional to what is left rather than as a
    // constant drain that empties the cell and keeps going. Taken here rather
    // than passed in: `n` is the start-of-substep state by definition, and the
    // solver overwrites it in place. See plasmaChemistryODE::extRef_.
    if (ext)
    {
        ode_->setExternalRef(n);
    }

    scalar dtTry = dt;
    solver_->solve(0.0, dt, n, dtTry);

    // Densities cannot be negative. An ODE solver can return a small negative
    // value for a species being consumed to near-exhaustion, and letting it
    // through turns the next step's rate law into nonsense.
    ode_->setExternal(nullptr);
    ode_->setExternalSlope(nullptr, 0);
    ode_->clearExternalRef();

    forAll(n, i)
    {
        n[i] = max(n[i], scalar(0));
    }
}


void Foam::plasmaChemistry::derivatives
(
    const scalarField& n, const scalarField& kTab,
    const scalar Tgas, scalarField& dndt
) const
{
    kTab_ = kTab;
    ode_->setTgas(Tgas);
    dndt.setSize(species_.size());
    ode_->derivatives(0.0, n, dndt);
}


void Foam::plasmaChemistry::jacobian
(
    const scalarField& n, const scalarField& kTab,
    const scalar Tgas, scalarField& dfdx, scalarSquareMatrix& dfdy
) const
{
    kTab_ = kTab;
    ode_->setTgas(Tgas);
    dfdx.setSize(species_.size());
    ode_->jacobian(0.0, n, dfdx, dfdy);
}


void Foam::plasmaChemistry::productionLoss
(
    const scalarField& n, const scalarField& kTab,
    const scalar Tgas, scalarField& P, scalarField& L
) const
{
    kTab_ = kTab;
    ode_->setTgas(Tgas);
    ode_->setExternal(nullptr);
    ode_->productionLoss(n, P, L);

    if (cantera_)
    {
        cantera_->productionLoss(n, Tgas, Pc_, Lc_);
        P += Pc_;
        L += Lc_;
    }
}


Foam::scalar Foam::plasmaChemistry::heavyHeatRelease
(
    const scalarField& n,
    const scalar Tgas
) const
{
    ode_->setTgas(Tgas);
    return ode_->heavyHeatRelease(n);
}


Foam::scalar Foam::plasmaChemistry::chargeResidual
(
    const scalarField& n,
    const scalarField& kTab,
    const scalar Tgas
) const
{
    kTab_ = kTab;
    ode_->setTgas(Tgas);
    return ode_->chargeResidual(n);
}


// ************************************************************************* //
