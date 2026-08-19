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

        //- NEGATIVE RATE COEFFICIENT guard.
        //
        //  Rate tables may be read with `outOfBounds extrapolate`, and linear
        //  extrapolation past the top of a NON-MONOTONIC k(E/N) can cross
        //  zero. A negative k is not merely inaccurate: q = k*PROD(n) then
        //  goes negative, which turns a production term into a sink and, worse,
        //  puts a NEGATIVE coefficient on the implicit diagonal via
        //  `fvm::Sp(L, n)` -- anti-damping, the same destabilising shape as the
        //  transport cross-term sink that had to be abandoned.
        //
        //  The code already defended the UNBOUNDED direction (the VSMALL ->
        //  L ~ 1e300 -> SIGFPE story below) but not this one. Clamping to zero
        //  is the conservative repair: the reaction simply stops rather than
        //  running backwards. It is REPORTED once, naming the reaction, because
        //  a table extrapolating negative is a mechanism defect to fix at
        //  source, not something to silently absorb.
        //
        //  Below-range lookups are safe: the tables start at zero.
        mutable label negativeRateCount_ = 0;
        mutable label negativeRateFirst_ = -1;
        mutable bool  negativeRateWarned_ = false;

        //- The rate coefficient for one reaction, guarded against the above.
        //  Single point of truth: derivatives(), jacobian() and
        //  productionLoss() all evaluated k independently, so a clamp applied
        //  to only some of them would leave the three disagreeing about the
        //  same reaction.
        inline scalar rateCoeff
        (
            const plasmaReactionSpec& rx,
            const label r
        ) const;

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

        //- Start-of-substep densities, used to make a REMOVING cross term
        //  proportional to the state instead of constant [m^-3].
        //
        //  A constant negative rate removes the same amount per unit time
        //  whatever is left in the cell, so it drives the species through zero
        //  once |ext| dt exceeds the content -- which is exactly the species
        //  Courant number exceeding one. The ODE then integrates a trajectory
        //  with negative densities, the rate laws turn to nonsense, and
        //  rodas23 grinds into its step limit. MEASURED: both Co 2.5 streamer
        //  arms died at step 14 with "Integration steps greater than maximum
        //  10000" and negative entries in the state vector, and neither did
        //  before the cross term reached the adaptive stiff path.
        //
        //  So the REMOVING part is carried as a linear sink instead:
        //
        //      e < 0:   dn/dt += e n/n_0     (rather than  dn/dt += e)
        //
        //  which equals e exactly at n = n_0, the state the transport rate was
        //  measured at, so nothing changes where the old form was valid. Away
        //  from it the removal decays with what is left and the species
        //  approaches zero asymptotically instead of crossing it. The Jacobian
        //  gains a negative diagonal e/n_0, which a Rosenbrock method uses:
        //  the term becomes stiff-but-well-posed rather than explicit.
        //
        //  This is the same trade the wall-loss model already makes, for the
        //  same reason -- exponential decay rather than n -= nu n dt, so a step
        //  longer than 1/nu decays instead of going negative.
        //
        //  ADDING cross terms (e > 0) stay constant: inflow cannot make a
        //  density negative, and scaling it would suppress filling of an empty
        //  cell, which is how a streamer front advances.
        //
        //  Empty when not supplied, in which case the cross term is used raw.
        scalarField extRef_;

        //- Cap on the sink coefficient, in units of 1/dt, from
        //  `chemCrossSinkLimit`. ZERO DISABLES THE SINK ENTIRELY, and that is
        //  the default, because the sink was MEASURED to cost accuracy in
        //  every form tried.
        //
        //  Streamer at Co 1.5, peak n_e at 2 ns against the Co 0.25 reference
        //  (9.453e19), all with the cross term on:
        //
        //      no sink      +42%
        //      cap 1/dt     +83%
        //      cap 3/dt     +67%
        //      cap 10/dt    +95%
        //
        //  The hypothesis was a removal deficit -- a cap of nu removes only
        //  1 - exp(-nu dt) of what transport takes -- which predicts the error
        //  falling MONOTONICALLY as the cap rises. It does not. The prediction
        //  failed, so the mechanism is refuted and this is not a knob to tune;
        //  the sink is kept only so the experiment can be repeated.
        //
        //  Positivity is instead protected where it actually mattered -- the
        //  Co 2.5 crash -- by catching the integrator's failure and falling
        //  back to the linearised source. See plasmaTransport.C.
        //
        //  Caveat on all four numbers: variants differing only in this one
        //  constant scatter over 28 percentage points, while front position
        //  agrees to within a cell. Head peak density at 2 ns is an
        //  ill-conditioned metric and none of these differences should be read
        //  as precise.
        scalar extSinkLimit_ = 0.0;

        //- Cross term for species i at substep time x, split into the part
        //  that is added as a rate and the linear-sink coefficient [1/s].
        //
        //  rate is the ADDING part, sink the coefficient multiplying y[i].
        inline void crossTerm
        (
            const label i,
            const scalar x,
            scalar& rate,
            scalar& sink
        ) const
        {
            rate = 0;
            sink = 0;
            if (!ext_) return;

            // The step MEAN decides which branch this species takes, and the
            // sink coefficient is FIXED across the substep.
            //
            // Letting the sink vary with x instead cost the Jacobian its
            // bound: d f/dx then carries the slope scaled by y/n_0, and the
            // slope here is (ext - extPrev)/dt ~ 1e42 in streamer-head cells.
            // Multiplying that by an unbounded ratio collapses the Rosenbrock
            // step controller -- MEASURED: the integrator stalled at
            // dxDid = 8.9e-31 against a step of 1.7e-12 and hit its step
            // limit. The removal's own time variation is a second-order
            // detail next to its proportionality to what is left, which is
            // the part that matters, so only the ADDING branch keeps the
            // slope.
            const scalar e0 = (*ext_)[i];

            // GUARDED, not unconditional. The sink is only needed where the
            // constant form would actually take the species through zero,
            // i.e. where the removal over the step is a significant fraction
            // of what the cell holds. Everywhere else -- the overwhelming
            // majority of cells -- the old constant form is safe, cheaper,
            // and bit-for-bit what it always was.
            //
            // Applying it to every removing species instead was measured at
            // 564 s to step 10 of the Co 1.5 streamer, against 16-20 s/step
            // for the runs that did not carry a sink at all: it turns an
            // explicit source into a stiff decay with a Jacobian diagonal in
            // cells that never needed one.
            //
            // theta = 1/2: at most half the content removed over the step.
            if (extRef_.empty())
            {
                rate = e0;
                if (extSlope_ && extDt_ > 0)
                {
                    rate += (*extSlope_)[i]*(x - 0.5*extDt_);
                }
                return;
            }

            const bool safeExplicit =
                (extDt_ > 0) && (-e0*extDt_ <= 0.5*extRef_[i]);

            if (e0 >= 0 || safeExplicit || extSinkLimit_ <= 0)
            {
                rate = e0;
                if (extSlope_ && extDt_ > 0)
                {
                    rate += (*extSlope_)[i]*(x - 0.5*extDt_);
                }
                return;
            }

            // n_0 can be zero where transport is removing from a cell the
            // linear solve had already emptied -- an artifact, not a rate.
            // Capping the coefficient keeps that from becoming an infinitely
            // fast decay the integrator would have to resolve.
            //
            // The cap is applied by COMPARISON, never by dividing and then
            // clamping: with n_0 at VSMALL the quotient overflows to infinity
            // before any min() sees it, and SIGFPE is raised inside
            // derivatives(). Measured -- that is exactly how the first version
            // of this died, on step 2 of the Co 1.5 streamer.
            //
            // THE CAP IS 1/dt, and that value is the whole cost story.
            //
            // Transport wanting to remove more than the cell holds is not a
            // fast decay -- the most that can leave is what is there. Capping
            // at 10/dt instead treated it as one, and since a far-field cell
            // holds n_0 ~ 0 of most excited species, ANY negative ext there
            // exceeds the cap and the whole quiescent mesh acquired a decay
            // ten times faster than the step. nu dt = 10 is stiff by
            // construction and rodas23 pays for it in substeps: MEASURED at
            // 50-90 s/step against 9-12 s/step for the build without any of
            // this, and guarding which cells got the sink did NOT help,
            // because those cells all failed the guard.
            //
            // At 1/dt the term removes ~63% of the content over the step and
            // adds no stiffness whatsoever. The real removal is the linearised
            // transport solve's job in any case; ext is a cross-term estimate,
            // not the transport operator.
            const scalar n0 = extRef_[i];
            const scalar nuMax =
                (extDt_ > 0) ? extSinkLimit_/extDt_ : 1e30;

            sink = (-e0 >= nuMax*n0) ? -nuMax : e0/n0;
        }

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

        //- Warn ONCE if any rate coefficient came back negative. Collective:
        //  call it on every rank. See negativeRateCount_ for why a negative k
        //  matters rather than merely being inaccurate.
        void reportNegativeRateGuard(const wordList& reactionNames) const;

        //- Constant external (transport) rate for this cell; null to disable
        void setExternal(const scalarField* e) { ext_ = e; }
        void setExternalSlope(const scalarField* s, const scalar dt)
        { extSlope_ = s; extDt_ = dt; }

        //- Start-of-substep state that a REMOVING cross term is scaled
        //  against; clear it to carry the cross term raw. See extRef_.
        void setExternalRef(const scalarField& r) { extRef_ = r; }
        void clearExternalRef() { extRef_.clear(); }

        //- Cap on the sink coefficient, in units of 1/dt. Default 1.
        void setCrossSinkLimit(const scalar c) { extSinkLimit_ = c; }

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


// * * * * * * * * * * * * * * Inline Definitions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::plasmaChemistryODE::rateCoeff
(
    const plasmaReactionSpec& rx,
    const label r
) const
{
    const scalar k =
        (rx.tabulated >= 0)
      ? kTab_[rx.tabulated]
      : rx.A*Foam::pow(Tgas_, rx.b)
          *(rx.Ta != 0 ? Foam::exp(-rx.Ta/Tgas_) : 1.0);

    // One predictable branch in a hot loop; see negativeRateCount_.
    if (k < 0)
    {
        ++negativeRateCount_;
        if (negativeRateFirst_ < 0)
        {
            negativeRateFirst_ = r;
        }
        return 0;
    }

    return k;
}


#endif

// ************************************************************************* //
