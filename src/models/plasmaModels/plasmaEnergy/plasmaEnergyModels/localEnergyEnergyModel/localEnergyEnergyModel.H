/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::localEnergyEnergyModel

Description
    Local Mean Energy Approximation (LMEA): a TRANSPORTED electron energy
    equation, replacing the LFA closure in which the mean energy is a local
    function of E/N.

    Solved CONSERVATIVELY for the electron energy density

        n_eps = n_e * eps_bar        [eV/m^3]

    not for eps_bar or T_e directly. Solving the intensive variable is
    non-conservative and breaks the identity eps_bar = n_eps/n_e that the
    source linearisation depends on.

        d(n_eps)/dt + div(Gamma_eps) = S_eps
        Gamma_eps = (5/3) mu_e E n_eps - (5/3) D_e grad(n_eps)
        S_eps     = -e Gamma_e . E  -  P_loss

    The 5/3 factors are the energy-flux moments of the two-term expansion
    (Hagelaar & Kroesen 2000, Eqs. 10-13).

    PUBLISHES `meanE` [eV] as a registered field. That is the whole
    integration story: `TabulatedProperty1D` looks its coordinate up from the
    registry BY NAME, so mobilities, diffusivities and rate coefficients move
    from LFA to LMEA by changing `lookupVariable reducedE` to the published
    mean-energy field and pointing at the `*_vs_meanE` tables that the
    Boltzmann solver already emits. No new evaluator code is involved.

    PHASE A (this file): the loss term is implicit through fvm::Sp, which is
    unambiguous because it is a genuine sink; Joule heating is explicit. The
    full Hagelaar bracket -- which is sign-INDEFINITE through its E.grad(n_e)
    diffusion term -- is deliberately NOT adopted here. See
    docs/electron-energy-lmea-plan.md Sect. 3.0: a converged PIMPLE outer loop
    already treats implicitly the couplings whose explicit evaluation causes
    the oscillation his scheme exists to damp, so the machinery may not be
    needed at all. Measure before adopting.

SourceFiles
    localEnergyEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef localEnergyEnergyModel_H
#define localEnergyEnergyModel_H

#include "plasmaEnergyModel.H"
#include <functional>
#include "plasmaPropertyEvaluator.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class localEnergyEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class localEnergyEnergyModel
:
    public plasmaEnergyModel
{
private:

    //- Electron energy density n_eps = n_e eps_bar [eV/m^3]. THE transported
    //  variable: conservative, so the equation can be written in divergence
    //  form and the loss kept implicit.
    volScalarField nEps_;

    //- Mean energy [eV], published for the tabulated lookups. Derived, not
    //  solved: eps_bar = n_eps / n_e.
    volScalarField meanE_;

    //- Electron temperature [K], for the T() interface. (2/3) eps_bar / k_B.
    volScalarField T_;

    //- Electron mobility and diffusivity, evaluated at the LOCAL MEAN ENERGY
    //  rather than at E/N -- which is what makes this LMEA rather than LFA.
    //- PARTICLE mobility mu_e. Feeds the JOULE term (-e Gamma_e . E acts on
    //  the electrons) and its derivative in the source sensitivity.
    autoPtr<plasmaPropertyEvaluator> muE_;

    //- ENERGY mobility mu_eps and diffusivity D_eps. Feed the energy FLUX
    //  only.
    //
    //  These are SEPARATE fields on purpose. Serving both roles from one
    //  evaluator was tried and was wrong by exactly mu_eps/mu_e = 1.512: the
    //  Joule heating came out 51% high, which pushed the equilibrium mean
    //  energy up with it. Measured, then fixed:
    //      correct J/L 0.8505 x 1.512 = 1.286, the value the code reported.
    autoPtr<plasmaPropertyEvaluator> muEps_;
    autoPtr<plasmaPropertyEvaluator> DE_;

    //- Power lost to elastic and inelastic collisions, per gas density
    //  [eV m^3/s], tabulated against mean energy by the Boltzmann solver.
    autoPtr<plasmaPropertyEvaluator> Pelastic_;
    autoPtr<plasmaPropertyEvaluator> Pinelastic_;

    //- Working fields for the above, so they are not reallocated per step.
    volScalarField muEf_;      // particle mobility
    volScalarField muEpsf_;    // energy mobility
    volScalarField DEf_;       // energy diffusivity
    volScalarField PlossN_;

    //- SOURCE SENSITIVITY dS/deps, Hagelaar's linearisation. Implicit in how
    //  the source RESPONDS to the mean energy, not merely in the magnitude of
    //  the loss.
    //
    //  MEASURED in 0-D (validation/lmea/dt_stability.csv, air at 200 Td):
    //  linearising only the loss magnitude survives to dt ~ 5 tau_eps and then
    //  oscillates 2.3-8.8 eV at 18 tau; adding the sensitivity holds
    //  4.85-5.36 eV there. About 9x on the oscillation-free step against
    //  explicit -- NOT the ~100x Hagelaar reports, and the difference is
    //  stated rather than papered over.
    //
    //  Enabled by default because ROBUSTNESS is the point: our streamer runs
    //  at dt/tau ~ 0.7, comfortably inside even the explicit limit, but a case
    //  we have not run has no such guarantee.
    bool sourceSensitivity_ = true;

    //- True when the case supplies the real muEpsN/DEpsN tables, false when
    //  falling back to (5/3) x the particle coefficients. Reported, because
    //  the difference is up to 60% in diffusion on this air set.
    bool energyCoeffsTabulated_ = false;

    //- Seed from the LFA equilibrium on the FIRST correct(), never at
    //  construction: there the Poisson solve has not run and reducedE is
    //  zero, so the lookup returns the bottom of the table.
    mutable bool seedFromLFA_ = false;
    autoPtr<plasmaPropertyEvaluator> seedEval_;

    //- Convect and diffuse the energy. DIAGNOSTIC SWITCH, default true.
    //
    //  Off makes every cell a local 0-D energy balance, which MUST then
    //  reproduce meanEnergy_vs_reducedE -- that is the local-equilibrium limit
    //  by definition. So it separates a wrong SOURCE from a wrong TRANSPORT
    //  term in a single short run, which arguing about a hot cell cannot do:
    //  LMEA is legitimately allowed to exceed the local equilibrium wherever
    //  energy is convected in, so a high value is not by itself a defect.
    bool transportEnergy_ = true;
    mutable bool chemSourceReported_ = false;

    //- Top of the tabulated mean-energy range, and whether the
    //  beyond-the-table warning has been issued. Extrapolation is ALLOWED --
    //  clamping would freeze the physics at the table edge, which is less
    //  robust, not more -- but it is reported.
    scalar tableMaxMeanE_ = 0;
    bool extrapWarned_ = false;
    mutable bool negativeWarned_ = false;
    volScalarField dSdEps_;

    //- DENSITY FLOOR for eps_bar = n_eps/n_e, which is 0/0 in the far field.
    //
    //  Flagged as a risk in the plan before any of this was written, and it is
    //  the same shape as the near-empty-cell pathology that killed the
    //  transport cross-term sink: a quantity divided by a density that tends
    //  to zero in exactly the cells that do not matter.
    scalar nEfloor_;

    //- Mean energy is clamped to the range the tables cover, so a transient
    //  excursion cannot be extrapolated into a negative rate coefficient --
    //  the failure the chemistry rate guard exists to catch.
    scalar meanEmin_;
    scalar meanEmax_;

    //- Report the mean-energy RANGE every N steps. The single most useful
    //  diagnostic this model has: a monotone climb to meanEmax_ is a runaway,
    //  while an alternating high/low is the period-2 oscillation Hagelaar's
    //  Fig. 1 shows and the 0-D sweep reproduced. The two need opposite
    //  responses, and the range distinguishes them at a glance where a single
    //  number cannot.
    label reportInterval_;
    label lastReport_;
    scalar meanEprev_;
    label clampedSteps_;

    //- State of the first cell to need the upper clamp, so the cause can be
    //  told from the symptom. See updateDerived().
    mutable label clampCell_ = -1;
    mutable scalar clampNe_ = GREAT;
    mutable scalar clampNEps_ = 0;
    mutable scalar clampRaw_ = 0;

    //- Outer iterations within the CURRENT step that needed the clamp, and
    //  the step they belong to. A transient inside the outer loop that later
    //  correctors recover is not the same failure as a converged solution
    //  outside the tabulated range, and only the second is a defect.
    mutable label lastClampTimeIndex_ = -1;
    mutable label clampedItersThisStep_ = 0;

    localEnergyEnergyModel(const localEnergyEnergyModel&) = delete;
    void operator=(const localEnergyEnergyModel&) = delete;

    //- Refresh meanE_ and T_ from n_eps and n_e.
    void updateDerived();


public:

    //- Runtime type information
    TypeName("localEnergy");

    // Constructors

        localEnergyEnergyModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex,
            const volVectorField& E
        );

    // Member functions

        virtual void correct() override;

        //- The electron energy equation. Unlike every other energy model,
        //  this returns a REAL matrix: that is what LMEA is.
        virtual tmp<fvScalarMatrix> eEqn() const override;

        virtual bool isUniform() const override { return false; }
        virtual bool isConstant() const override { return false; }

        virtual tmp<volScalarField> T() const override;
        virtual const dimensionedScalar& Tvalue() const override;

        //- Largest energy-transport rate in the domain [1/s], from the
        //  ACTUAL energy flux rather than 5/3 x the species one.
        //
        //      rate = 0.5 SUM_faces |phi_eps| / V
        //
        //  Co_energy = rate * dt, the same construction the species Courant
        //  limiter uses. The 5/3 estimate it replaces was doubly wrong: it
        //  assumed the Maxwellian ratio (measured 1.24-1.57, not 1.667) AND
        //  it used mu evaluated at reducedE rather than at meanE.
        scalar maxEnergyRate() const;

        //- Largest ENERGY RELAXATION rate in the domain, 1/tau_eps =
        //  P_loss/n_eps [1/s]. A stiffness constraint on the step, distinct
        //  from the convective rate above.
        scalar maxEnergyRelaxationRate() const;

        //- OPTION 4 hand-off. True when this model's mobility and loss
        //  properties can be evaluated at a scalar eps, which the chemistry
        //  ODE needs in order to integrate the energy alongside the species.
        bool canExportEnergyCoefficients() const;

        //- The transported energy density field.
        const volScalarField& nEps() const { return nEps_; }

        //- mu*N as a function of eps [1/(V m s)]. RAW table units -- the
        //  caller divides by the gas density, which changes with T_gas.
        std::function<scalar(scalar)> muNofEps() const;

        //- Total power loss per unit gas density as a function of eps
        //  [eV m^3/s], elastic plus inelastic.
        std::function<scalar(scalar)> PlossNofEps() const;

        //- The bounds and floor the ODE must apply, so the integrator and the
        //  field model cannot disagree about what eps means.
        scalar meanEnergyMin() const { return meanEmin_; }
        scalar meanEnergyMax() const { return meanEmax_; }
        scalar electronDensityFloor() const { return nEfloor_; }

        //- Mean electron energy [eV].
        const volScalarField& meanE() const { return meanE_; }
};


} // End namespace Foam

#endif

// ************************************************************************* //
