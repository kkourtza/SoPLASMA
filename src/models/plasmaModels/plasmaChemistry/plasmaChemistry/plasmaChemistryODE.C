/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
\*---------------------------------------------------------------------------*/

#include "plasmaChemistryODE.H"
#include "scalarMatrices.H"
#include "plasmaChemistryCantera.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasmaChemistryODE::plasmaChemistryODE
(
    const List<plasmaReactionSpec>& reactions,
    const scalarField& kTab,
    const scalarField& charge,
    const label nSpecie
)
:
    reactions_(reactions),
    kTab_(kTab),
    Tgas_(300.0),
    charge_(charge),
    nSpecie_(nSpecie),
    rate_(reactions.size(), Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plasmaChemistryODE::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    // The transport contribution enters as a rate LINEAR in x, so the
    // trajectory the chemistry follows is the one the cell actually takes
    // rather than the one it would take if it were closed. `x` was unnamed
    // here while the term was constant.
    // Linear across the substep, centred so the mean stays *ext_. A constant
    // cross term caps the coupling at first order however well converged --
    // measured, see plasmaChemistryODE.H. crossTerm() applies that centring
    // and splits the result: an ADDING transport rate is a constant source, a
    // REMOVING one a linear sink proportional to what is left, so the species
    // cannot be driven through zero by transport at any Courant number.
    dydx = Zero;
    if (ext_)
    {
        forAll(dydx, i)
        {
            scalar rate, sink;
            crossTerm(i, x, rate, sink);
            dydx[i] = rate + sink*max(y[i], scalar(0));
        }
    }

    forAll(reactions_, r)
    {
        const plasmaReactionSpec& rx = reactions_[r];
        if (!heavy_ && rx.tabulated < 0) continue;

        // Rate coefficient. Electron-impact rates were interpolated from the
        // EEDF tables before the substep began and are frozen; heavy rates are
        // Arrhenius in the gas temperature.
        //
        // Ta is tested against zero, not against positive: a NEGATIVE
        // activation temperature is physical and common in this mechanism --
        // O2(b1) + O2 has Ta = -241 K -- and skipping the exponential for it
        // dropped a factor of 2.2 from that rate. Found by the Cantera
        // cross-check, which is the argument for having two backends at all.
        const scalar k = rateCoeff(rx, r);

        // Rate = k * PRODUCT(n_reactant), with repetition giving the order.
        // A negative density cannot occur physically but does occur
        // transiently inside an ODE step; clipping at zero here keeps the RHS
        // finite instead of letting an odd-order product change sign and drive
        // the solver further negative.
        scalar q = k*rx.fixedReactantDensity;
        forAll(rx.reactants, i)
        {
            const label s = rx.reactants[i];
            if (s >= 0)
            {
                q *= max(y[s], scalar(0));
            }
        }

        if (rx.collider >= 0)
        {
            q *= max(y[rx.collider], scalar(0));
        }
        else if (rx.colliderFixedDensity > 0)
        {
            q *= rx.colliderFixedDensity;
        }

        rate_[r] = q;

        forAll(rx.reactants, i)
        {
            const label s = rx.reactants[i];
            if (s >= 0) dydx[s] -= q;
        }
        forAll(rx.products, i)
        {
            const label s = rx.products[i];
            if (s >= 0) dydx[s] += q;
        }
    }

    // Heavy contribution from the external backend, as production minus the
    // loss it was split into. Evaluated at the CURRENT state y, so it tracks
    // the trajectory through the substep exactly as the native terms do --
    // freezing it at the entry state would make the heavy chemistry a constant
    // source and cap the step at first order.
    if (ct_)
    {
        ct_->productionLoss(y, Tgas_, Pc_, Lc_);
        forAll(dydx, s)
        {
            dydx[s] += Pc_[s] - Lc_[s]*max(y[s], scalar(0));
        }
    }
}


void Foam::plasmaChemistryODE::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    // EXPLICIT TIME DEPENDENCE, and it must be reported here.
    //
    // With a linear cross term the right-hand side does depend on x, and
    // rodas23 is a ROSENBROCK method: it uses df/dx directly. Leaving this
    // zero while derivatives() varies with x costs ORDER, not merely
    // accuracy -- which is the whole point of the change. d(T)/dx is the
    // slope itself; the chemistry part still has none, since k is frozen
    // over the substep and Tgas is constant.
    //
    // Everything below the slope term is the original reasoning:
    // d/dt has no explicit time dependence: k is frozen over the substep and
    // Tgas is constant, so df/dx is identically zero. Saying so is not an
    // approximation -- an implicit solver that assumed otherwise would be
    // integrating a term that does not exist.
    dfdx = Zero;
    if (ext_ && extSlope_ && extDt_ > 0)
    {
        dfdx = *extSlope_;
    }
    dfdy = Zero;

    // The REMOVING part of the cross term depends on the state, so it belongs
    // in dfdy as well: d/dy (sink y) = sink, a NEGATIVE diagonal. Supplying it
    // is what turns transport-driven depletion from an explicit source the
    // integrator must resolve with tiny substeps into a stiff decay a
    // Rosenbrock method takes in one -- the failure mode this repairs.
    //
    // Where the sink is active the derivative is sink*y rather than the raw
    // rate, so dfdx carries the slope scaled the same way; an inconsistent
    // dfdx costs ORDER here for the reason recorded above.
    if (ext_ && !extRef_.empty())
    {
        forAll(y, i)
        {
            scalar rate, sink;
            crossTerm(i, x, rate, sink);
            if (sink != 0)
            {
                dfdy(i, i) += sink;

                // The sink is constant across the substep by construction, so
                // this species' cross term has NO explicit time dependence and
                // the slope must not be reported for it. See crossTerm().
                dfdx[i] = 0;
            }
        }
    }

    // Under `chemistryBackend cantera` the heavy reactions are SKIPPED here
    // and their block is contributed by Cantera below, from its own analytic
    // composition derivative. Taking the heavy Jacobian from the native
    // Arrhenius data instead would be exact only while the two backends
    // describe the same reactions -- and inexact in precisely the cases
    // Cantera is selected for (falloff, PLOG, Chebyshev, reversible). A
    // Rosenbrock method is not Newton-iterated, so that costs ORDER, not just
    // convergence rate; the whole scheme would quietly drop to first.

    forAll(reactions_, r)
    {
        const plasmaReactionSpec& rx = reactions_[r];
        if (!heavy_ && rx.tabulated < 0) continue;

        const scalar k = rateCoeff(rx, r);

        scalar base = k*rx.fixedReactantDensity;
        if (rx.collider >= 0)
        {
            base *= max(y[rx.collider], scalar(0));
        }
        else if (rx.colliderFixedDensity > 0)
        {
            base *= rx.colliderFixedDensity;
        }

        // d(rate)/d(n_s) for each transported reactant, by differentiating the
        // product one factor at a time. Written analytically rather than by
        // finite differences: a numerical Jacobian on a system whose rates
        // span twenty orders of magnitude loses the small entries entirely,
        // and those are the stiff ones the implicit solver needs.
        forAll(rx.reactants, i)
        {
            const label si = rx.reactants[i];
            if (si < 0) continue;

            scalar d = base;
            forAll(rx.reactants, j)
            {
                if (j == i) continue;
                const label sj = rx.reactants[j];
                if (sj >= 0) d *= max(y[sj], scalar(0));
            }

            forAll(rx.reactants, m)
            {
                const label sm = rx.reactants[m];
                if (sm >= 0) dfdy(sm, si) -= d;
            }
            forAll(rx.products, m)
            {
                const label sm = rx.products[m];
                if (sm >= 0) dfdy(sm, si) += d;
            }
        }

        // ...and with respect to a transported collider, which multiplies the
        // rate without appearing on either side.
        if (rx.collider >= 0)
        {
            scalar d = k*rx.fixedReactantDensity;
            forAll(rx.reactants, j)
            {
                const label sj = rx.reactants[j];
                if (sj >= 0) d *= max(y[sj], scalar(0));
            }

            forAll(rx.reactants, m)
            {
                const label sm = rx.reactants[m];
                if (sm >= 0) dfdy(sm, rx.collider) -= d;
            }
            forAll(rx.products, m)
            {
                const label sm = rx.products[m];
                if (sm >= 0) dfdy(sm, rx.collider) += d;
            }
        }
    }

    // Heavy block from the external backend, in the same state ordering. Its
    // entries need no unit conversion: the rate and the state carry the same
    // factor of Avogadro's number and it cancels in the derivative.
    if (ct_)
    {
        ct_->jacobian(y, Tgas_, dfdy);
    }
}


Foam::scalar Foam::plasmaChemistryODE::chargeResidual(const scalarField& y) const
{
    // Measured on the REACTIONS only. Transport moves charge between cells, so
    // including it here would report a non-zero residual for a perfectly
    // balanced mechanism.
    const scalarField* saved = ext_;
    const_cast<plasmaChemistryODE*>(this)->ext_ = nullptr;

    scalarField dydx(nSpecie_, Zero);
    derivatives(0.0, y, dydx);

    const_cast<plasmaChemistryODE*>(this)->ext_ = saved;

    scalar net = 0.0, traffic = 0.0;
    forAll(dydx, s)
    {
        net     += charge_[s]*dydx[s];
        traffic += mag(charge_[s]*dydx[s]);
    }
    return (traffic > VSMALL) ? mag(net)/traffic : 0.0;
}


// ************************************************************************* //


void Foam::plasmaChemistryODE::reportNegativeRateGuard
(
    const wordList& tabulatedIds
) const
{
    if (negativeRateWarned_) return;

    // Collective: the offending cells may live on any rank.
    label count = negativeRateCount_;
    reduce(count, sumOp<label>());

    if (count == 0) return;

    negativeRateWarned_ = true;

    // Lowest offending reaction index anywhere. -1 means "none on this rank"
    // and must not win the minimum.
    label first = (negativeRateFirst_ < 0) ? labelMax : negativeRateFirst_;
    reduce(first, minOp<label>());

    word what("reaction " + Foam::name(first));
    if (first >= 0 && first < reactions_.size())
    {
        const label t = reactions_[first].tabulated;
        what += (t >= 0 && t < tabulatedIds.size())
              ? " (tabulated channel " + tabulatedIds[t] + ")"
              : " (heavy/Arrhenius)";
    }

    WarningInFunction
        << "a rate coefficient evaluated NEGATIVE and was clamped to zero."
        << nl
        << "    occurrences so far: " << count
        << ", first at " << what << nl
        << "    A negative k makes q = k*PROD(n) negative, which turns a"
        << " production term into a sink and puts an ANTI-DAMPING coefficient"
        << " on the implicit diagonal via fvm::Sp." << nl
        << "    Most likely cause: a rate table read with `outOfBounds"
        << " extrapolate`, linearly extrapolated past the top of a"
        << " non-monotonic k(E/N) until it crosses zero." << nl
        << "    Clamping keeps the run stable, but the MECHANISM is what needs"
        << " fixing: extend the table's range, or switch that table to"
        << " `outOfBounds clamp`." << nl
        << "    This warning is printed once per run." << endl;
}


void Foam::plasmaChemistryODE::productionLoss
(
    const scalarField& y,
    scalarField& P,
    scalarField& L
) const
{
    P.setSize(nSpecie_); P = Zero;
    L.setSize(nSpecie_); L = Zero;

    forAll(reactions_, r)
    {
        const plasmaReactionSpec& rx = reactions_[r];
        if (!heavy_ && rx.tabulated < 0) continue;

        const scalar k = rateCoeff(rx, r);

        scalar q = k*rx.fixedReactantDensity;
        forAll(rx.reactants, i)
        {
            const label s = rx.reactants[i];
            if (s >= 0) q *= max(y[s], scalar(0));
        }
        if (rx.collider >= 0)
        {
            q *= max(y[rx.collider], scalar(0));
        }
        else if (rx.colliderFixedDensity > 0)
        {
            q *= rx.colliderFixedDensity;
        }

        // Products are production. Reactants are loss, expressed as a
        // coefficient by dividing out the species' own density -- which is
        // exact here, because the species appears as a factor of q.
        forAll(rx.products, i)
        {
            const label s = rx.products[i];
            if (s >= 0) P[s] += q;
        }
        forAll(rx.reactants, i)
        {
            const label s = rx.reactants[i];
            if (s < 0) continue;
            // The loss COEFFICIENT is q/n, so a vanishing density makes it
            // unbounded. VSMALL (1e-300) as the floor produces L ~ 1e300,
            // which overflows the matrix and raises SIGFPE in the linear
            // solver -- found by a stress test, not in review.
            //
            // The floor is physical rather than numerical: below one particle
            // per cubic metre there is nothing to lose, and a sink for an
            // absent species is meaningless.
            const scalar ns = max(y[s], scalar(0));
            if (ns > 1.0)
            {
                L[s] += q/ns;
            }
            else
            {
                P[s] -= q;
            }
        }
    }
}


Foam::scalar Foam::plasmaChemistryODE::heavyHeatRelease
(
    const scalarField& y
) const
{
    scalar Q = 0.0;

    forAll(reactions_, r)
    {
        const plasmaReactionSpec& rx = reactions_[r];

        // Heavy only: tabulated channels are accounted by PgasN. See the
        // header for why mixing the two double-counts.
        if (rx.tabulated >= 0 || rx.deltaH == 0.0) continue;

        const scalar k = rateCoeff(rx, r);

        scalar q = k*rx.fixedReactantDensity;
        forAll(rx.reactants, i)
        {
            const label s = rx.reactants[i];
            if (s >= 0) q *= max(y[s], scalar(0));
        }
        if (rx.collider >= 0)
        {
            q *= max(y[rx.collider], scalar(0));
        }
        else if (rx.colliderFixedDensity > 0)
        {
            q *= rx.colliderFixedDensity;
        }

        // deltaH is products minus reactants, so exothermic is negative and
        // the heat delivered to the gas is -deltaH.
        Q += q*(-rx.deltaH);
    }

    return Q;
}
