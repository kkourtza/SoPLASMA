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

class plasmaChemistryCantera;

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

    //- Reaction enthalpy [eV], products minus reactants, from the mechanism's
    //  own thermo. EXOTHERMIC is NEGATIVE, so the heat released to the gas is
    //  -deltaH per event -- emitted with the thermodynamic sign so it cannot be
    //  mistaken for "heat out".
    //
    //  This is where fast gas heating comes from: quenching of the electronic
    //  states by O2 releases most of what the electrons stored, tens of
    //  nanoseconds after they stored it. Computed by mechc from species
    //  enthalpies, never apportioned by hand.
    scalar deltaH = 0.0;

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

        //- Whether the HEAVY (Arrhenius) reactions are evaluated here.
        //
        //  False under `chemistryBackend cantera`: the same reactions are then
        //  evaluated by Cantera and added by the caller, and evaluating them
        //  in both places would double every heavy rate. The reactions stay in
        //  the list either way, because the Jacobian is still taken from them.
        bool heavy_ = true;

        //- Heavy backend, when it is not this class. Non-null exactly when
        //  heavy_ is false: the heavy reactions have to enter the integrated
        //  trajectory from SOMEWHERE, and leaving them out of derivatives()
        //  would silently drop recombination and quenching from the ODE path
        //  while the implicit-rate path still had them.
        const plasmaChemistryCantera* ct_ = nullptr;

        //- Cantera production and loss, sized on first use
        mutable scalarField Pc_, Lc_;

        //- Constant external rate added to every derivative [m^-3 s^-1].
        //
        //  The TRANSPORT contribution over the step, held fixed while the
        //  chemistry integrates. Without it the ODE follows a trajectory in
        //  which the cell is closed -- no drift in, no diffusion out -- and
        //  the mean rate along that trajectory is only a first-order
        //  approximation to the mean along the true one. Including it, and
        //  iterating it to convergence in the outer loop, is what makes the
        //  coupling second order.
        //
        //  MEASURED CORRECTION (2026-08-18). Holding it CONSTANT is not
        //  enough: with the outer loop converged to 1e-12 and ext_ at its
        //  fixed point, `solver ode` still measured p = 0.988 on the 40x40
        //  study. The transport rate VARIES across the step, and integrating
        //  a nonlinear f against a constant stand-in perturbs the trajectory
        //  at O(dt) however exact that constant is. So the claim above --
        //  that iterating a constant ext_ to convergence gives second order
        //  -- is false, and the slope below is what repairs it.
        //
        //  Null when not supplied, which is the closed-cell case.
        const scalarField* ext_ = nullptr;

        //- Time derivative of that external rate [m^-3 s^-2], so the cross
        //  term is LINEAR across the substep rather than constant:
        //
        //      T(x) = ext_ + extSlope_*(x - dt/2)
        //
        //  Centred on dt/2 so the step MEAN is exactly ext_ -- the quantity
        //  that was already validated -- and only the variation is new. A
        //  slope lagged by one timestep is O(dt) accurate, and it enters
        //  multiplied by (x - dt/2) ~ dt, so T is O(dt^2) and the coupling
        //  reaches second order without needing the slope to be exact.
        const scalarField* extSlope_ = nullptr;

        //- Substep length, needed to centre the slope. Set by integrate().
        scalar extDt_ = 0;

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

        //- Constant external (transport) rate for this cell; null to disable
        void setExternal(const scalarField* e) { ext_ = e; }
        void setExternalSlope(const scalarField* s, const scalar dt)
        { extSlope_ = s; extDt_ = dt; }

        //- Evaluate the heavy reactions here (true) or leave them to the
        //  Cantera backend (false). Affects derivatives() and
        //  productionLoss(); jacobian() always includes them, see there.
        void setHeavy(const bool b) { heavy_ = b; }

        //- Heavy backend to consult in derivatives(); null for native
        void setCantera(const plasmaChemistryCantera* c) { ct_ = c; }
        bool heavy() const { return heavy_; }

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

        //- Heat released to the gas by the HEAVY reactions [eV m^-3 s^-1].
        //
        //  sum over heavy reactions of  q_r * (-deltaH_r)
        //
        //  This is fast gas heating. The electronic states store what the
        //  electrons gave them, and quenching by O2 hands it to the gas tens
        //  of nanoseconds later; the delay is physical and comes out of the
        //  kinetics rather than being modelled.
        //
        //  DELIBERATELY EXCLUDES the electron-impact reactions, and the
        //  exclusion is the whole reason this is a separate function rather
        //  than a loop over `reactions_`. Their contribution to gas heating is
        //  the tabulated PgasN, which is a sum over EVERY electron-impact
        //  channel including the `energy-only` ones that carry no species and
        //  never reach this reaction list at all. Adding a per-reaction term
        //  here on top of the table would count the resolved channels twice
        //  and still miss the lumped ones.
        scalar heavyHeatRelease(const scalarField& y) const;

        //- Split the instantaneous reaction rate into production and a LOSS
        //  COEFFICIENT: dn_s/dt = P_s - L_s * n_s.
        //
        //  The loss is returned as a coefficient rather than a rate so the
        //  caller can treat it implicitly (fvm::Sp). Two things follow from
        //  that, and both matter for a stiff system driven by an explicit
        //  source:
        //
        //    * POSITIVITY. An implicit sink cannot drive a density negative,
        //      whatever the timestep. An explicit one can, and does, as soon
        //      as a loss timescale falls below dt.
        //    * CONVERGENCE. The Picard iteration over the outer loop is a
        //      contraction on the loss term once it is implicit; left
        //      explicit, a fast sink makes it diverge.
        //
        //  Evaluated at the state given -- not averaged over a step -- because
        //  BDF2 requires the right-hand side at t^{n+1}. Supplying a mean
        //  instead is a Lie composition and caps the scheme at first order.
        void productionLoss
        (
            const scalarField& y,
            scalarField& P,
            scalarField& L
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
