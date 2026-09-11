/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::plasmaChemistryCantera

Description
    Heavy chemistry evaluated by Cantera instead of by the native reaction list.

    A sibling of the native path, not a replacement: it answers the same
    question -- production and loss per species at a state -- so the caller and
    the equation assembly do not change. `chemistryBackend cantera`.

    WHAT IT EVALUATES, and what it does not. Cantera is given a HEAVY-ONLY
    mechanism (`<mech>.heavy.yaml`, written by mechc). The electron-impact
    rates continue to come from the Boltzmann tables. That division is what
    makes the backend cheap and correct at once:

      * evaluating the electron-impact reactions here as well would count them
        twice, since the tabulated rates are added by the same assembly;
      * Cantera evaluates an electron-impact rate from the EEDF, which would
        have to be injected PER CELL -- the expense the tables exist to avoid.

    WHAT IT BUYS. Reaction forms the native parser does not implement:
    pressure-dependent falloff (Troe, SRI, Lindemann), PLOG, Chebyshev,
    reversible reactions with equilibrium constants from the thermo data, and
    surface chemistry. A mechanism using any of those cannot be expressed
    natively today, and this is the path for it.

    WHAT IT COSTS. A libcantera dependency in the OpenFOAM build, and a state
    set plus a kinetics evaluation per cell. Cantera works in kmol/m^3; the
    conversion to number density is applied here so the caller sees the same
    units either way.

SourceFiles
    plasmaChemistryCantera.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaChemistryCantera_H
#define plasmaChemistryCantera_H

#include "scalarField.H"
#include "wordList.H"
#include "fileName.H"
#include "labelList.H"
#include "scalarMatrices.H"

#include <memory>

namespace Cantera
{
    class Solution;
}

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class plasmaChemistryCantera Declaration
\*---------------------------------------------------------------------------*/

class plasmaChemistryCantera
{
    // Private Data

        std::shared_ptr<Cantera::Solution> sol_;

        //- State-vector index -> Cantera species index; -1 when Cantera does
        //  not carry that species (the electron, typically, which takes no
        //  part in heavy chemistry)
        labelList toCantera_;

        //- The inverse: Cantera species index -> state index, -1 for species
        //  Cantera carries that the case does not transport. Held explicitly
        //  because the Jacobian is scattered FROM Cantera's ordering, and
        //  searching toCantera_ per non-zero would make that quadratic.
        labelList fromCantera_;

        //- Species the caller carries, in state order
        wordList species_;

        //- NUMBER DENSITIES [m^-3] of species Cantera knows but the case does
        //  not transport, held fixed. The background is most of the gas, so
        //  omitting it would not merely perturb the rates -- a three-body
        //  reaction would lose its collider entirely.
        scalarField background_;

        //- Working arrays, sized once
        mutable scalarField X_;
        mutable scalarField creation_;
        mutable scalarField destruction_;

        //- Whether Cantera's analytic composition derivative is available for
        //  this mechanism. Probed once at construction rather than trusted:
        //  netProductionRates_ddCi() is declared on the base Kinetics class and
        //  throws NotImplementedError for managers that do not implement it, so
        //  the only honest test is to call it.
        bool analyticJacobian_ = true;

        //- Working state for the finite-difference fallback
        mutable scalarField nPert_, Pp_, Lp_, P0_, L0_;


    // Private Member Functions

        //- Push a state onto Cantera's phase. Factored out because the
        //  Jacobian needs the state set but not the rates evaluated, and
        //  reaching it through productionLoss() would pay for a kinetics
        //  evaluation it then discards.
        //
        //  Returns false when there is nothing there to react.
        bool setState(const scalarField& n, const scalar Tgas) const;


public:

    // Constructors

        //- Construct from the heavy-only mechanism and the carried species
        plasmaChemistryCantera
        (
            const fileName& heavyMechanism,
            const wordList& species,
            const wordList& canteraNames
        );

        ~plasmaChemistryCantera();


    // Member Functions

        //- Production and loss coefficient from the HEAVY reactions at this
        //  state, in the same units and the same form as the native path:
        //  dn_s/dt = P_s - L_s n_s, number density per second.
        //
        //  Added to whatever the native electron-impact evaluation produced,
        //  never replacing it.
        void productionLoss
        (
            const scalarField& n,
            const scalar Tgas,
            scalarField& P,
            scalarField& L
        ) const;

        //- Fixed number density of an untransported species [m^-3]. Silently
        //  ignored when Cantera does not carry the species, which is the
        //  correct behaviour: it cannot then appear in any heavy reaction.
        void setBackground(const word& canteraName, const scalar n);

        //- d(dn_i/dt)/d(n_j) of the HEAVY reactions, ADDED to dfdy.
        //
        //  Cantera returns d(omega_i)/d(C_j) in kmol/m^3 per kmol/m^3 per
        //  second. Both our rate and our state carry the same factor of
        //  Avogadro's number, so the two cancel exactly and the entries
        //  transfer with NO unit conversion -- unlike the rates themselves.
        //
        //  Analytic where Cantera implements it, finite-difference otherwise.
        //  Never omitted and never taken from the other backend: a Rosenbrock
        //  method is not Newton-iterated, so an inexact Jacobian costs ORDER,
        //  not just convergence rate.
        void jacobian
        (
            const scalarField& n,
            const scalar Tgas,
            scalarSquareMatrix& dfdy
        ) const;

        //- Force the finite-difference path even where the analytic
        //  derivative is available. `canteraJacobian finiteDifference`.
        //  Exists so the fallback is reachable and therefore testable: a
        //  branch that only runs on kinetics managers nobody here uses is a
        //  branch nobody has ever run.
        void useFiniteDifferenceJacobian() { analyticJacobian_ = false; }

        //- True when the Jacobian above is Cantera's analytic derivative,
        //  false when it falls back to finite differences
        bool analyticJacobian() const { return analyticJacobian_; }

        //- Species Cantera did not recognise, for reporting
        wordList unmatched() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
