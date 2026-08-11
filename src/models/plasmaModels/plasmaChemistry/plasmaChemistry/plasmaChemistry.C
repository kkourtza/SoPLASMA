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
    if (mech.found("heavyReactions"))
    {
        const List<dictionary> heavy(mech.lookup("heavyReactions"));
        forAll(heavy, i)
        {
            const dictionary& h = heavy[i];
            plasmaReactionSpec rx;
            rx.A  = h.get<scalar>("A");
            rx.b  = h.get<scalar>("b");
            rx.Ta = h.get<scalar>("Ta");

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

    solver_ = ODESolver::New(*ode_, odeDict_);

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


void Foam::plasmaChemistry::integrate
(
    scalarField& n,
    const scalarField& kTab,
    const scalar Tgas,
    const scalar dt,
    const scalarField* ext
) const
{
    kTab_ = kTab;
    ode_->setTgas(Tgas);
    ode_->setExternal(ext);

    scalar dtTry = dt;
    solver_->solve(0.0, dt, n, dtTry);

    // Densities cannot be negative. An ODE solver can return a small negative
    // value for a species being consumed to near-exhaustion, and letting it
    // through turns the next step's rate law into nonsense.
    ode_->setExternal(nullptr);

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
