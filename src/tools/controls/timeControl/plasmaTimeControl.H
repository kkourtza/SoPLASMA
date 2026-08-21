/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::plasmaTimeControl

Description
    Class for automatic time-step control in plasma simulations.

    This class serves as a specialized replacement for the standard OpenFOAM
    'setDeltaT.H' functionality. It calculates the maximum permissible time
    step based on plasma-specific stability criteria.

    The class reads the 'plasmaTimeControl' dictionary (found within
    system/plasmaSimulationControl) to determine which criteria to apply and
    their limiting coefficients. It ultimately sets runTime.deltaT().

    If 'adjustTimeStep' is false, this class performs no action.

SourceFiles
    plasmaTimeControl.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaTimeControl_H
#define plasmaTimeControl_H

#include "fvCFD.H"
#include "IOdictionary.H"

#include "SoPLASMAConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class plasmaTransport;    
/*---------------------------------------------------------------------------*\
                      Class plasmaTimeControl Declaration
\*---------------------------------------------------------------------------*/

class plasmaTimeControl
{
private:

    // Private data

        //- Reference to the Time object
        Time& runTime_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Dictionary to hold control parameters
        dictionary dict_;

        // Control switches and limits

            //- Switch to enable/disable time step adjustment
            bool adjustTimeStep_;

            //- Maximum physical time step allowed
            scalar maxDeltaT_;

        // Plasma specific criteria

            //- Switch to limit Dielectric Relaxation Ratio
            Switch limitDielectricRelaxationRatio_;

            //- Switch to print Dielectric Relaxation Ratio
            Switch printDielectricRelaxationRatio_;

            //- Coefficient for Dielectric Relaxation ratio 
            scalar maxDielectricRelaxationRatio_;

            //- Switch to enable/disable Species Courant limiting
            Switch limitSpeciesCo_;

            //- Switch to print Species Courant limiting
            Switch printSpeciesCo_;

            //- Maximum allowed convective Courant for the tracked species
            scalar maxSpeciesConvectiveCo_;

            //- Maximum allowed diffusive Courant for the tracked species
            scalar maxSpeciesDiffusiveCo_;

            //- Name of the species to track for Courant limiting (default: "e")
            word courantSpeciesName_;

            //- Switch to limit by Chemistry (Source term)
            Switch limitChemistryCo_;

            //- Switch to print Chemistry source restriction
            Switch printChemistryCo_;

            //- Maximum allowed relative change due to chemistry
            scalar maxChemistryCo_;

            //- Switch to limit voltage rise rate
            Switch limitVoltageRiseRate_;

            //- Switch to print voltage rise rate info
            Switch printVoltageRiseRate_;

            //- Maximum allowed voltage rise rate [V/s]
            scalar maxVoltageRiseRate_;

            //- Name of the patch where voltage is applied
            word voltagePatchName_;

            //- Previous timestep voltage on the applied patch [V]
            scalar prevPatchVoltage_;

        // Outer (PIMPLE) coupling convergence
        //
        // The Poisson and species equations are lagged relative to each other
        // within a timestep. Converging that loop is what makes the scheme
        // second order -- with a single corrector it is FIRST order whatever
        // `ddtSchemes` says, because no amount of ddt accuracy repairs a lag.
        // Order is therefore not something a case selects; it is what the case
        // ACHIEVES, and the run has to say which it got.

            //- Chase the coupling to `tolerance` (true), or accept a single
            //  lagged corrector (false)
            bool outerChaseConvergence_;

            //- Cap on outer correctors. A cap, not a target: the loop exits as
            //  soon as residualControl is satisfied.
            label outerMaxCorrectors_;

            //- Coupling residual target. MUST exceed the linear-solver
            //  tolerance, or the criterion can never be met and the loop
            //  silently runs to the cap on every step.
            scalar outerTolerance_;

            //- Has prevPatchVoltage_ been given a real value yet? Until it
            //  has, differencing against it would report the whole applied
            //  voltage as a one-step rise. See the first-step guard.
            bool voltageSeeded_;

            //- warn | reduceDeltaT | retryStep | fatal
            //
            //  `reduceDeltaT` ACCEPTS the non-converged step and shrinks the
            //  NEXT one, so an unconverged solution enters the history and
            //  stays there. Measured on the streamer: 11% of steps at Co 1.5
            //  and 16% at Co 2.5, against 0% for the Co 0.25 reference.
            //
            //  `retryStep` discards it instead -- the fields are reset to
            //  oldTime(), deltaT is halved and the step is re-run, so no
            //  unconverged state is ever kept. Measured correctors per step
            //  are bimodal (converge within ~8, or not at all), so a retry
            //  after a LOW cap is cheaper than accepting after a high one.
            word outerOnNonConvergence_;

            //- Retries already spent on the current step, and the cap on them.
            //  Halving 5 times is deltaT/32; past that the step is accepted
            //  with a warning rather than the run deadlocking.
            label outerRetries_ = 0;
            label outerMaxRetries_ = 5;

            //- Time index the retry counter belongs to, so it resets itself
            //  when the step changes without needing a call from the solver.
            label outerRetryTimeIndex_ = -1;

            //- Time value at the END of the last ACCEPTED step, captured by
            //  the solver before it advances the clock.
            //
            //  It must be captured, NOT inferred. The obvious
            //  `value() - deltaT()` is wrong here: ++runTime advances using
            //  the PREVIOUS step's deltaT, while adjustDeltaT() then installs
            //  a new one, so the two differ by the ramp factor (1.2x) -- and
            //  across repeated retries the error compounds, because value()
            //  has already been corrected once. MEASURED: that arithmetic sent
            //  the clock NEGATIVE on the third retry (t = -6e-14).
            scalar tStepStart_ = 0;

            //- REJECTION MEMORY: the smallest deltaT known to have failed to
            //  converge, and the ceiling derived from it.
            //
            //  Without it the controller rediscovers the failure scale over
            //  and over. MEASURED on the Co 1.5 streamer: steps failed at
            //  ~1.0e-11, succeeded at ~5e-12, regrew at 1.2x/step back through
            //  the failure scale in about four steps, and failed again -- a
            //  stable sawtooth between 4.3e-12 and 8.5e-12 that spent one
            //  discarded attempt every fourth step for no accuracy gain,
            //  because the SUSTAINABLE step lies between the two.
            //
            //  So a failure is remembered: growth is capped just below it
            //  (retryCeilingFactor_), and the ceiling is relaxed only slowly
            //  (retryCeilingRelax_ per successful step) so the solver can
            //  still follow the physics -- on a timescale far longer than the
            //  failure cycle. Zero disables the memory.
            scalar dtFailCeiling_ = 0;
            scalar retryCeilingFactor_ = 0.9;
            scalar retryCeilingRelax_ = 1.05;

            //- Set for the step FOLLOWING a rejection, which may not grow at
            //  all. Standard in Hairer-style step control, and on its own it
            //  removes the tightest oscillation: a step that has just been
            //  halved has no business immediately climbing again.
            bool noGrowthNextStep_ = false;

            //- ABSOLUTE FLOOR on deltaT, and whether it was derived rather
            //  than set by the user.
            //
            //  Without a floor the retry ladder descends forever. MEASURED on
            //  a deliberately broken case (`maxCorrectors 3`, where nothing
            //  can converge): deltaT fell 1.44e-12 -> 5.4e-14 -> 2.025e-15 in
            //  three steps, because a step can lose 32x (five halvings) while
            //  regrowth is capped at 1.2x per step. The run then never reaches
            //  endTime -- for a user the worst failure mode there is: no
            //  crash, no message, just a job that never finishes.
            //
            //  The default is derived from the requirement that the run be
            //  able to FINISH: below (endTime - startTime)/minDeltaTMaxSteps_
            //  it cannot, whatever the physics. That is a genuine impossibility
            //  rather than a tuned magic number, which is why it is safe to
            //  apply automatically.
            //  TWO TERMS, because either alone is wrong:
            //
            //    tSpan/minDeltaTMaxSteps_   -- "this run can never finish"
            //    deltaT0_/minDeltaTStepFactor_ -- "this is far below the
            //                                  physical scale you started at"
            //
            //  The first scales with run LENGTH, which decouples from the
            //  physical timescale: a 1 ms discharge would derive 1e-10 s while
            //  its streamer still needs ~1e-12 s, so the floor would sit ABOVE
            //  the step the Courant limiter legitimately demands and -- since
            //  the floor is applied last and overrides the physical limiters
            //  -- would silently clamp dt up to an unstable value. The second
            //  term anchors to the scale the user themselves specified, so the
            //  floor stays orders of magnitude below any real step size.
            //  Taking the MINIMUM keeps the guard while removing that risk.
            scalar minDeltaT_ = 0;
            bool   minDeltaTAuto_ = true;
            scalar minDeltaTMaxSteps_ = 1e7;
            scalar minDeltaTStepFactor_ = 1e4;

            //- The deltaT this run STARTED from, captured once. It must not be
            //  re-sampled on a later read(): deltaT falls as the controller
            //  works, so re-sampling would ratchet the floor down in step with
            //  the very collapse it exists to stop.
            scalar deltaT0_ = 0;

            //- Outer-coupling stability margin from the adaptive (Aitken)
        //  relaxation: min(omega) over the last step. 1 when adaptive
        //  relaxation is off, so neither threshold below can bind.
        scalar relaxMargin_ = 1.0;

        //- Largest factor this step, reported beside the minimum
        scalar relaxOmegaMax_ = 1.0;

        //- Below this, deltaT is not allowed to GROW
        scalar relaxMarginHold_ = 0.5;

        //- Below this, deltaT is actively SHRUNK by relaxMarginBackoff_
        scalar relaxMarginShrink_ = 0.2;

        //- Shrink factor applied below relaxMarginShrink_
        scalar relaxMarginBackoff_ = 0.7;

        //- Whether the margin governed deltaT this step, for reporting
        bool relaxMarginBound_ = false;

        //- Steps accepted WITHOUT the outer loop converging, because the
            //  dt lever was exhausted (at the floor, or out of retries).
            //  Graceful degradation: the run completes and the user is told
            //  exactly which steps are suspect, rather than the solver either
            //  lying or deadlocking.
            label degradedSteps_ = 0;
            label degradedRun_ = 0;
            label degradedConsecutive_ = 0;
            label maxDegradedConsecutive_ = 10;

            //- Times the floor actually clamped a requested reduction.
            label floorHits_ = 0;

            //- TREND detection, so a run ratcheting down is visible BEFORE it
            //  becomes terminal. Absolute floors only fire at the very end of
            //  a collapse; the ratio against the largest step the run has
            //  sustained catches it while there is still something to do
            //  about it.
            scalar dtHighWater_ = 0;
            scalar dtCollapseRatio_ = 100;

            //- ENERGY transport RATE [1/s] when an LMEA energy equation is
            //  transported, 0 when there is none. Co_energy = rate*dt.
            //
            //  The electron ENERGY convects and diffuses at 5/3 the particle
            //  rates (Hagelaar & Kroesen Eqs. 11-13), so its Courant and
            //  diffusion numbers are 5/3 of the species ones computed here.
            //  A case switching to LMEA at its existing
            //  maxSpeciesConvectiveCo is therefore silently 67% over that
            //  limit -- MEASURED: at Co 1.5 the energy ran at ~2.5 and the
            //  convective term drove n_eps negative on step 1, in cells whose
            //  SOURCE was strongly positive.
            //
            //  Set from the energy model rather than by the user, so enabling
            //  LMEA cannot leave the step silently unsafe.
            scalar energyRate_ = 0;
            bool   collapseWarned_ = false;

            //- ENERGY RELAXATION rate, 1/tau_eps = P_loss/n_eps [1/s].
            //
            //  A SEPARATE constraint from energyRate_ above, which is a
            //  CONVECTIVE (Courant) rate. The two are limited by different
            //  coefficients on purpose: a Courant number of order 1 is
            //  correct for transport and badly wrong for a stiff relaxation.
            //
            //  MEASURED on the 1.15M-cell streamer, 2026-08-21: at the cells
            //  that reached the 100 eV clamp, tau_eps = n_eps/P_loss = 1.65 ps
            //  while the step chosen by the SPECIES Courant number was
            //  5.16 ps -- dt/tau = 3.1. The mean energy then oscillated
            //  between the clamp and ~7 eV and the run died on SIGFPE. A
            //  reference run pinned at dt = 0.25 ps (dt/tau = 0.15) passed
            //  through the same window with ZERO clamps and a smooth monotone
            //  settle, which is what fixes the ratio below at 0.2.
            //
            //  Set from the energy model, like energyRate_, so enabling LMEA
            //  cannot leave the step silently unsafe.
            scalar energyRelaxRate_ = 0;
            bool   limitEnergyRelaxation_ = true;
            scalar maxEnergyRelaxationRatio_ = 0.2;
            bool   energyRelaxWarned_ = false;

            //- Set when the last step exhausted the cap without converging;
            //  read by adjustDeltaT() so the response is a time-step decision
            //  alongside the Courant and dielectric criteria, not a separate
            //  mechanism.
            bool outerHitCap_;

            //- Correctors actually used on the last step, peak-held between
            //  reports
            label outerItersUsed_;

            //- Steps since the last outer-coupling report
            label outerReportCounter_;

        // Private Member Functions

            //- Read controls from dictionary
            void read();

            //- Compute average voltage on the applied patch
            scalar patchVoltageAvg(const plasmaTransport& transport) const;

public:

    // Constructors

        //- Construct from Time and Mesh
        plasmaTimeControl(Time& runTime, const fvMesh& mesh);

    //- Destructor
    virtual ~plasmaTimeControl() {}

    // Member Functions

        //- Adjust deltaT using the specific plasmaTransport object
        //  Note: We pass it by const reference
        void adjustDeltaT(const plasmaTransport& transport);

        //- Set the first time step of the simulation
        void setInitialDeltaT(const plasmaTransport& transport);

        //- Translate the `outerCoupling` block into PIMPLE's own
        //  nOuterCorrectors and residualControl.
        //
        //  STATIC, and called BEFORE pimpleControl is constructed, because
        //  that is when PIMPLE reads its dictionary. Driving OpenFOAM's loop
        //  rather than replacing it keeps `finalIter()` exact -- and the
        //  chemistry trailing half-step and the gas-energy solve are both
        //  gated on it, so a hand-rolled loop that got it subtly wrong would
        //  advance the chemistry by the wrong multiple of dt.
        static void configureOuterCoupling(fvMesh& mesh);

        //- Record what the outer loop did on this step.
        //  `used` is the number of correctors executed. Convergence is
        //  inferred from the loop exiting BEFORE the cap, which is exactly
        //  OpenFOAM's own residualControl decision rather than a second,
        //  possibly disagreeing, criterion of ours.
        void noteOuterLoop(const label used);

        //- Should the step just solved be DISCARDED and re-run smaller?
        //
        //  True only under `onNonConvergence retryStep`, when the loop hit its
        //  cap and retries remain. The caller must then reset the solution to
        //  oldTime() and call prepareRetry() before solving again.
        bool stepRejected() const;

        //- Record the time the step is starting FROM. Call it in the solver
        //  loop BEFORE ++runTime, which is the only place the value is
        //  unambiguous. See tStepStart_ for why it cannot be inferred later.
        void noteStepStart() { tStepStart_ = runTime_.value(); }

        //- Declare that an LMEA energy equation is transported, so the step
        //  is limited on ITS Courant number too. See energyCourantFactor_.
        void setEnergyRate(const scalar r)
        {
            energyRate_ = r;
        }

        //- Declare the ENERGY RELAXATION rate 1/tau_eps = P_loss/n_eps, so
        //  the step is limited by how fast the energy equilibrates as well as
        //  by how fast it convects. See energyRelaxRate_ for the measurement
        //  that motivates it.
        void setEnergyRelaxRate(const scalar r)
        {
            energyRelaxRate_ = r;
        }

        //- Halve deltaT and move the clock back onto the shortened step.
        //
        //  The time VALUE has already advanced by the old deltaT, so it is
        //  corrected here rather than rewound: t_prev + deltaT_new. The time
        //  INDEX is deliberately left alone -- every per-step cache in the
        //  solver is keyed on it, and oldTime() rotates on it, so holding it
        //  fixed is what makes a retry land on the same start state.
        void prepareRetry();

        //- Is deltaT sitting on the floor, so that halving it again would
        //  achieve nothing? Retrying costs a full solve, so there is no point
        //  spending one on a step size the controller is not allowed to use.
        bool atMinDeltaT() const;

        //- Record that a step was accepted WITHOUT the outer loop converging.
        //  Warns per occurrence and, if degradation becomes systemic rather
        //  than isolated, stops the run with an actionable diagnostic.
        void noteDegradedStep();

        //- Did the last step exhaust the corrector cap without converging?
        inline bool outerHitCap() const
        {
            return outerHitCap_;
        }

        //- Number of non-converged steps accepted over the whole run.
        inline label degradedRun() const
        {
            return degradedRun_;
        }

    // Accessors (Getters)

        //- Is adjustable time stepping enabled?
        inline bool adjustTimeStep() const
        {
            return adjustTimeStep_;
        }

        //- Is dielectric relaxation limiting enabled?
        inline bool limitDielectricRelaxationRatio() const
        {
            return limitDielectricRelaxationRatio_;
        }

        //- Get the max allowed dielectric relaxation ratio
        inline scalar maxDielectricRelaxationRatio() const
        {
            return maxDielectricRelaxationRatio_;
        }

        //- Is species Courant number limiting enabled?
        inline bool limitSpeciesCo() const
        {
            return limitSpeciesCo_;
        }

        //- Get the max allowed convective species Courant number
        inline scalar maxSpeciesConvectiveCo() const
        {
            return maxSpeciesConvectiveCo_;
        }

        //- Get the max allowed diffusive species Courant number
        inline scalar maxSpeciesDiffusiveCo() const
        {
            return maxSpeciesDiffusiveCo_;
        }

        //- Get the species name used for Courant number limiting
        inline const word& courantSpeciesName() const
        {
            return courantSpeciesName_;
        }

        //- Is chemistry limiting enabled?
        inline bool limitChemistryCo() const
        {
            return limitChemistryCo_;
        }

        //- Get the max allowed chemistry limit
        inline scalar maxChemistryCo() const
        {
            return maxChemistryCo_;
        }

        //- Is voltage rise rate limiting enabled?
        inline bool limitVoltageRiseRate() const
        {
            return limitVoltageRiseRate_;
        }

        //- Get the maximum allowed voltage rise rate [V/dt]
        inline scalar maxVoltageRiseRate() const
        {
            return maxVoltageRiseRate_;
        }

        //- Get the name of the patch where voltage is applied
        inline const word& voltagePatchName() const
        {
            return voltagePatchName_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
