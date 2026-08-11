/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::plasmaChemistryODE

Description
    The plasma chemistry of one cell, as an ODE system.

    Derives from Foam::ODESystem and nothing else, which is the whole point.
    OpenFOAM's chemistryModel / reactingFoam stack is the obvious base and the
    wrong one: it is built on psiReactionThermo, so it carries mass fractions
    and a single temperature where this needs number densities and a separate
    electron temperature; every ReactionRate type in it is k(T) Arrhenius,
    where electron-impact rates come from the EEDF as k(E/N); and its species
    are bound to a thermo package that would have to be fabricated for N2(A3).

    ODESystem asks for two things -- derivatives() and jacobian() -- and knows
    nothing about thermodynamics. Implementing it gives Rosenbrock12/23/34,
    SIBS, EulerSI and the adaptive Runge-Kutta family through
    ODESolver::New(dict), with none of that baggage.

    State vector: species NUMBER DENSITIES [m^-3], in the order the caller
    supplies. Not mass fractions -- an electron's mass fraction is 1e-9 of a
    plasma that is entirely controlled by it, and the resulting stiffness is
    numerical rather than physical.

    Two rate families, one assembly:

      electron-impact   k_j(E/N) from the Boltzmann tables, passed in already
                        interpolated. Frozen for the duration of the substep:
                        under operator splitting the field does not change
                        while the chemistry integrates, so re-interpolating
                        inside the ODE would be cost without meaning.

      heavy             k = A T^b exp(-Ta/T) in the GAS temperature, evaluated
                        here. No EEDF involved, which is why these need no
                        Boltzmann solve and can run without Cantera.

    Charge conservation is a property of the RIGHT-HAND SIDE, not of the
    result: if every reaction balances charge then sum(q_i dn_i/dt) is zero
    identically, for any state and any timestep. checkCharge() asserts that on
    the assembled derivative, which is a stronger statement than comparing
    densities before and after a step.

SourceFiles
    plasmaChemistryODE.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaChemistryODE_H
#define plasmaChemistryODE_H

#include "ODESystem.H"
#include "scalarField.H"
#include "labelList.H"
#include "wordList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class plasmaReactionSpec Declaration
\*---------------------------------------------------------------------------*/

//- One reaction, in the form the ODE evaluates.
//
//  Stoichiometry is held as index lists WITH REPETITION, exactly as the
//  mechanism dictionary writes it: `2 O` is (O O). Repetition rather than
//  coefficients keeps the rate law and the stoichiometry consistent by
//  construction -- the reaction order IS the length of the reactant list, so
//  the two cannot disagree.
struct plasmaReactionSpec
{
    //- Indices into the state vector; -1 marks a species that is not carried
    //  (background gas), whose density is a fixed parameter instead.
    labelList reactants;
    labelList products;

    //- Fixed density product for untransported reactants [m^-3]^n. Folded in
    //  once at setup rather than per RHS evaluation.
    scalar fixedReactantDensity = 1.0;

    //- Third body, if any. -1 when transported, index otherwise; the density
    //  multiplies the rate without the species appearing on either side.
    label collider = -1;
    scalar colliderFixedDensity = 0.0;

    //- k = A T^b exp(-Ta/T) for heavy reactions.
    scalar A = 0.0, b = 0.0, Ta = 0.0;

    //- Index into the pre-interpolated electron-impact rate array, or -1 for
    //  a heavy reaction evaluated from A, b, Ta.
    label tabulated = -1;

    //- Charge carried by each side, for the balance assertion
    scalar chargeImbalance = 0.0;
};


/*---------------------------------------------------------------------------*\
                       Class plasmaChemistryODE Declaration
\*---------------------------------------------------------------------------*/

class plasmaChemistryODE
:
    public ODESystem
{
    // Private Data

        const List<plasmaReactionSpec>& reactions_;

        //- Electron-impact rate coefficients at this cell, already
        //  interpolated at the local field and frozen for the substep
        const scalarField& kTab_;

        //- Gas temperature [K], for the Arrhenius rates
        scalar Tgas_;

        //- Charge number per state entry, for the conservation assertion
        const scalarField& charge_;

        label nSpecie_;

        //- Reaction rates, sized once and reused: derivatives() is called
        //  several times per step by every ODE solver worth using.
        mutable scalarField rate_;


public:

    // Constructors

        plasmaChemistryODE
        (
            const List<plasmaReactionSpec>& reactions,
            const scalarField& kTab,
            const scalarField& charge,
            const label nSpecie
        );


    // Member Functions

        //- Number of equations: one per transported species
        label nEqns() const override { return nSpecie_; }

        //- Gas temperature used by the Arrhenius rates
        void setTgas(const scalar T) { Tgas_ = T; }

        void derivatives
        (
            const scalar x,
            const scalarField& y,
            scalarField& dydx
        ) const override;

        void jacobian
        (
            const scalar x,
            const scalarField& y,
            scalarField& dfdx,
            scalarSquareMatrix& dfdy
        ) const override;

        //- Largest |sum(q_i dn_i/dt)| relative to the total charge traffic.
        //
        //  Zero to roundoff for a charge-conserving mechanism, at ANY state --
        //  it is a property of the reactions, not of the solution. Checking it
        //  on the derivative catches an unbalanced reaction immediately and
        //  unambiguously, where checking densities after a step confounds it
        //  with integration error.
        scalar chargeResidual(const scalarField& y) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
