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
    const scalar,
    const scalarField& y,
    scalarField& dydx
) const
{
    dydx = Zero;

    forAll(reactions_, r)
    {
        const plasmaReactionSpec& rx = reactions_[r];

        // Rate coefficient. Electron-impact rates were interpolated from the
        // EEDF tables before the substep began and are frozen; heavy rates are
        // Arrhenius in the gas temperature.
        scalar k =
            (rx.tabulated >= 0)
          ? kTab_[rx.tabulated]
          : rx.A*Foam::pow(Tgas_, rx.b)
              *(rx.Ta > 0 ? Foam::exp(-rx.Ta/Tgas_) : 1.0);

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
}


void Foam::plasmaChemistryODE::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    // d/dt has no explicit time dependence: k is frozen over the substep and
    // Tgas is constant, so df/dx is identically zero. Saying so is not an
    // approximation -- an implicit solver that assumed otherwise would be
    // integrating a term that does not exist.
    dfdx = Zero;
    dfdy = Zero;

    forAll(reactions_, r)
    {
        const plasmaReactionSpec& rx = reactions_[r];

        scalar k =
            (rx.tabulated >= 0)
          ? kTab_[rx.tabulated]
          : rx.A*Foam::pow(Tgas_, rx.b)
              *(rx.Ta > 0 ? Foam::exp(-rx.Ta/Tgas_) : 1.0);

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
}


Foam::scalar Foam::plasmaChemistryODE::chargeResidual(const scalarField& y) const
{
    scalarField dydx(nSpecie_, Zero);
    derivatives(0.0, y, dydx);

    scalar net = 0.0, traffic = 0.0;
    forAll(dydx, s)
    {
        net     += charge_[s]*dydx[s];
        traffic += mag(charge_[s]*dydx[s]);
    }
    return (traffic > VSMALL) ? mag(net)/traffic : 0.0;
}


// ************************************************************************* //
