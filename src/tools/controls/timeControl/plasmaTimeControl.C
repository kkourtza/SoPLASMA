/*---------------------------------------------------------------------------*\
  File: plasmaTimeControl.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaTimeControl.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaTimeControl.H"
#include "plasmaTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaTimeControl::plasmaTimeControl(Time& runTime, const fvMesh& mesh)
:
    runTime_(runTime),
    mesh_(mesh),
    dict_(dictionary::null),
    adjustTimeStep_(false),
    maxDeltaT_(GREAT),
    limitDielectricRelaxationRatio_(false),
    printDielectricRelaxationRatio_(false),
    maxDielectricRelaxationRatio_(1.0),
    limitSpeciesCo_(false),
    printSpeciesCo_(false),
    maxSpeciesConvectiveCo_(1.0),
    maxSpeciesDiffusiveCo_(1.0),
    courantSpeciesName_("e"),
    limitChemistryCo_(false),
    printChemistryCo_(false),
    maxChemistryCo_(1.0),
    limitVoltageRiseRate_(false),
    printVoltageRiseRate_(false),
    maxVoltageRiseRate_(GREAT),
    voltagePatchName_(""),
    prevPatchVoltage_(0.0),
    voltageSeeded_(false),
    outerChaseConvergence_(true),
    outerMaxCorrectors_(20),
    outerTolerance_(1e-8),
    outerOnNonConvergence_("reduceDeltaT"),
    outerHitCap_(false),
    outerItersUsed_(0),
    outerReportCounter_(0)
{
    read();

    Info << "Plasma time control succesfully constructed." << nl << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void plasmaTimeControl::read()
{
    IOdictionary plasmaDict
        (
            IOobject
            (
                "plasmaSimulationControls",
                runTime_.system(),
                runTime_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );


    if (plasmaDict.found("plasmaTimeControl"))
    {
        dict_ = plasmaDict.subDict("plasmaTimeControl");

        adjustTimeStep_ =
            dict_.lookupOrDefault<Switch>("adjustTimeStep", false);

        maxDeltaT_ =
            dict_.lookupOrDefault<scalar>("maxDeltaT", GREAT);

        // ABSOLUTE FLOOR. Derived from two independent impossibilities, the
        // smaller winning -- see minDeltaT_ in the header for why one term
        // alone is unsafe. Deriving it from what CANNOT be true, rather than
        // guessing a working step size, is what makes it safe to apply
        // without the user having chosen it.
        minDeltaTMaxSteps_ =
            dict_.lookupOrDefault<scalar>("minDeltaTMaxSteps", 1e7);

        minDeltaTStepFactor_ =
            dict_.lookupOrDefault<scalar>("minDeltaTStepFactor", 1e4);

        maxDegradedConsecutive_ =
            dict_.lookupOrDefault<label>("maxConsecutiveDegradedSteps", 10);

        dtCollapseRatio_ =
            dict_.lookupOrDefault<scalar>("deltaTCollapseRatio", 100);

        minDeltaTAuto_ = !dict_.found("minDeltaT");

        const scalar tSpan =
            max(runTime_.endTime().value() - runTime_.startTime().value(),
                SMALL);

        // Captured ONCE -- see deltaT0_. On a restart this is the step the
        // restart began from, which may already be reduced; that only lowers
        // the floor, which is the safe direction.
        if (deltaT0_ <= 0)
        {
            deltaT0_ = runTime_.deltaTValue();
        }

        const scalar floorFromSpan =
            tSpan/max(minDeltaTMaxSteps_, scalar(1));

        const scalar floorFromStep =
            (deltaT0_ > 0)
          ? deltaT0_/max(minDeltaTStepFactor_, scalar(1))
          : GREAT;

        const scalar derivedFloor = min(floorFromSpan, floorFromStep);

        minDeltaT_ = minDeltaTAuto_
            ? derivedFloor
            : dict_.get<scalar>("minDeltaT");

        // AUTO-CORRECT the unambiguous configuration mistakes. Each one has a
        // single defensible repair, so repairing it beats refusing to run --
        // but never silently.
        if (minDeltaT_ <= 0)
        {
            WarningInFunction
                << "plasmaTimeControl/minDeltaT = " << minDeltaT_
                << " is not positive, so it cannot bound the retry ladder."
                << nl << "    Using the derived floor "
                << derivedFloor << " s instead." << endl;

            minDeltaT_ = derivedFloor;
            minDeltaTAuto_ = true;
        }

        if (minDeltaT_ > maxDeltaT_)
        {
            WarningInFunction
                << "plasmaTimeControl/minDeltaT (" << minDeltaT_
                << ") exceeds maxDeltaT (" << maxDeltaT_
                << "), which leaves no admissible step size." << nl
                << "    Lowering minDeltaT to maxDeltaT/1000 = "
                << maxDeltaT_/1000 << " s." << endl;

            minDeltaT_ = maxDeltaT_/1000;
        }

        if (adjustTimeStep_ && runTime_.deltaTValue() < minDeltaT_)
        {
            WarningInFunction
                << "the initial deltaT (" << runTime_.deltaTValue()
                << ") is already below minDeltaT (" << minDeltaT_ << ")."
                << nl << "    Raising the initial step to the floor."
                << endl;

            runTime_.setDeltaT(minDeltaT_);
        }

        if (adjustTimeStep_)
        {
            Info<< "plasmaTimeControl: deltaT floor " << minDeltaT_ << " s";

            if (minDeltaTAuto_)
            {
                Info<< " (derived, "
                    << (floorFromStep < floorFromSpan
                          ? "from the initial step"
                          : "from endTime")
                    << ")";
            }
            else
            {
                Info<< " (set by minDeltaT)";
            }

            Info<< nl
                << "    below it the run stops with a diagnostic rather than"
                << " descending forever" << endl;
        }

        // Dielectric relaxation
        limitDielectricRelaxationRatio_ =
         dict_.lookupOrDefault<Switch>("limitDielectricRelaxationRatio", false);

        printDielectricRelaxationRatio_ =
         dict_.lookupOrDefault<Switch>("printDielectricRelaxationRatio", false);

        maxDielectricRelaxationRatio_ =
            dict_.lookupOrDefault<scalar>("maxDielectricRelaxationRatio", 1.0);

        // Species Courant
        limitSpeciesCo_ = 
            dict_.lookupOrDefault<Switch>("limitSpeciesCo", false);

        printSpeciesCo_ = 
            dict_.lookupOrDefault<Switch>("printSpeciesCo", false);

        if (limitSpeciesCo_ || printSpeciesCo_)
        {
            maxSpeciesConvectiveCo_ =
                dict_.lookupOrDefault<scalar>("maxSpeciesConvectiveCo", 1.0);

            maxSpeciesDiffusiveCo_ =
                dict_.lookupOrDefault<scalar>("maxSpeciesDiffusiveCo", 1.0);

            courantSpeciesName_ =
                dict_.lookupOrDefault<word>("courantSpeciesName", "e");
        }

        // Read OUTSIDE the limitSpeciesCo_ block: the energy relaxation limit
        // is a stiffness constraint, not a Courant choice, and switching the
        // Courant limiter off must not silently disable it.
        limitEnergyRelaxation_ =
            dict_.lookupOrDefault<Switch>("limitEnergyRelaxation", true);

        maxEnergyRelaxationRatio_ =
            dict_.lookupOrDefault<scalar>("maxEnergyRelaxationRatio", 0.2);

        // Chemistry Courant
        limitChemistryCo_ = 
            dict_.lookupOrDefault<Switch>("limitChemistryCo", false);

        printChemistryCo_ = 
            dict_.lookupOrDefault<Switch>("printChemistryCo", false);

        maxChemistryCo_ = 
            dict_.lookupOrDefault<scalar>("maxChemistryCo", 1.0);

        // Outer-coupling settings. Kept in step with configureOuterCoupling()
        // so the reported cap is the cap PIMPLE was actually given.
        {
            const dictionary& oc = plasmaDict.subOrEmptyDict("outerCoupling");
            outerChaseConvergence_ =
                (oc.getOrDefault<word>("target", "converged") == "converged");
            outerMaxCorrectors_ =
                oc.getOrDefault<label>("maxCorrectors", 20);
            outerTolerance_ =
                oc.getOrDefault<scalar>("tolerance", 1e-8);
            outerOnNonConvergence_ =
                oc.getOrDefault<word>("onNonConvergence", "reduceDeltaT");
            outerMaxRetries_ =
                oc.getOrDefault<label>("maxRetries", 5);

            // Rejection memory. `retryCeilingFactor 0` disables it and
            // restores the plain 1.2x-per-step regrowth.
            retryCeilingFactor_ =
                oc.getOrDefault<scalar>("retryCeilingFactor", 0.9);
            retryCeilingRelax_ =
                oc.getOrDefault<scalar>("retryCeilingRelax", 1.05);

            // Stability margin from the ADAPTIVE outer-loop relaxation.
            // plasmaTransport reports min(omega) over the step; ~1 means the
            // Picard coupling was benign at this deltaT, small means it was
            // only stable because Aitken damped it hard.
            //
            // Without this feedback a damped run simply walks deltaT back up
            // to the new stability boundary and sits on it -- MEASURED: a
            // fixed factor of 0.9 passed 6.19 ps, which let the controller
            // grow deltaT to 8.03 ps, where the same limit cycle returned.
            // Nothing else in this controller knows that mode exists.
            //
            // Inert unless adaptiveRelaxation is on: relaxationMargin() then
            // returns 1 and neither threshold can bind.
            // THRESHOLDS: "the relaxation is running out of room", NOT
            // "the coupling would be unstable undamped".
            //
            // omega ~= 1/(1+g), so omega < 0.5 is exactly g > 1 -- the
            // undamped instability threshold. Keying the governor there was
            // MEASURED (2026-08-21) to throttle deltaT constantly on the LMEA
            // case, whose omega naturally sits at 0.37-0.71: it reached only
            // 5.16 ps where a fixed factor of 0.8 ran stably at 7.65 ps. But
            // g > 1 is the situation relaxation EXISTS to handle, so it is the
            // wrong trigger.
            //
            // deltaT only needs limiting when the relaxation can no longer
            // cope, i.e. when omega approaches its own floor (relaxOmegaMin,
            // default 0.05). Hence thresholds just above it.
            relaxMarginHold_ =
                oc.getOrDefault<scalar>("relaxMarginHold", 0.15);
            relaxMarginShrink_ =
                oc.getOrDefault<scalar>("relaxMarginShrink", 0.08);
            relaxMarginBackoff_ =
                oc.getOrDefault<scalar>("relaxMarginBackoff", 0.7);

            if (retryCeilingRelax_ < 1)
            {
                FatalErrorInFunction
                    << "outerCoupling/retryCeilingRelax = "
                    << retryCeilingRelax_ << " is below 1, so the ceiling"
                    << " would shrink on every step and deltaT could never"
                    << " recover." << exit(FatalError);
            }

            if
            (
                outerOnNonConvergence_ != "warn"
             && outerOnNonConvergence_ != "reduceDeltaT"
             && outerOnNonConvergence_ != "retryStep"
             && outerOnNonConvergence_ != "fatal"
            )
            {
                FatalErrorInFunction
                    << "outerCoupling/onNonConvergence = `"
                    << outerOnNonConvergence_ << "` is not recognised." << nl
                    << "    Use warn | reduceDeltaT | retryStep | fatal."
                    << exit(FatalError);
            }

            if (outerOnNonConvergence_ == "retryStep" && !adjustTimeStep_)
            {
                FatalErrorInFunction
                    << "outerCoupling/onNonConvergence = `retryStep` needs"
                    << " adjustTimeStep on: the retry works by shortening"
                    << " deltaT." << exit(FatalError);
            }
        }

        // Voltage rise rate
        limitVoltageRiseRate_ =
            dict_.lookupOrDefault<Switch>("limitVoltageRiseRate", false);

        printVoltageRiseRate_ =
            dict_.lookupOrDefault<Switch>("printVoltageRiseRate", false);

        if (limitVoltageRiseRate_ || printVoltageRiseRate_)
        {
            // VOLTS PER STEP, despite the name -- the limiter solves for the
            // deltaT that makes the applied voltage change by this much in one
            // step, and the report prints it as `dV/step [V]`. It is NOT a
            // rate in V/s. Kept under the old key so existing cases still
            // read, but documented here because the name misleads.
            //
            // 100 V/step is a sane ceiling for a kV-scale pulse: fine enough
            // that the rise is resolved, coarse enough never to bind while the
            // Courant and dielectric limits are doing their job. It is a
            // SAFETY NET, not the primary control.
            maxVoltageRiseRate_ =
                dict_.lookupOrDefault<scalar>("maxVoltageRisePerStep", 100.0);
            maxVoltageRiseRate_ =
                dict_.lookupOrDefault<scalar>
                (
                    "maxVoltageRiseRate", maxVoltageRiseRate_
                );

            voltagePatchName_ =
                dict_.lookupOrDefault<word>("voltagePatchName", "");

            if (voltagePatchName_.empty())
            {
                FatalIOErrorInFunction(dict_)
                    << "'voltagePatchName' must be specified when "
                    << "'limitVoltageRiseRate' or 'printVoltageRiseRate' "
                    << "is true."
                    << nl << exit(FatalIOError);
            }
        }
    }
}

scalar plasmaTimeControl::patchVoltageAvg
(
    const plasmaTransport& transport
) const
{
    const volScalarField& phi = transport.species().em().ePotential();
    const label patchi = mesh_.boundaryMesh().findPatchID(voltagePatchName_);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch '" << voltagePatchName_ << "' not found in mesh."
            << nl << exit(FatalError);
    }

    const scalarField& phiPatch = phi.boundaryField()[patchi];

    // Compute local contributions (zero on ranks with no faces on this patch)
    scalar localSum   = phiPatch.empty() ? 0.0 : sum(phiPatch);
    label  localCount = phiPatch.size();

    // All ranks participate in both reductions
    reduce(localSum,   sumOp<scalar>());
    reduce(localCount, sumOp<label>());

    return localCount > 0 ? localSum / scalar(localCount) : 0.0;
}

void plasmaTimeControl::adjustDeltaT(const plasmaTransport& transport)
{
    read();

    const scalar currentDeltaT = runTime_.deltaTValue();
    const scalar eps0 = constant::plasma::epsilon0.value();
    scalar newDeltaT = maxDeltaT_;

    scalar maxSigma = 0.0;
    scalar maxConvFluxRate  = 0.0;
    scalar maxDiffFluxRate  = 0.0;
    // scalar meanConvFluxRate  = 0.0;
    // scalar meanDiffFluxRate  = 0.0;
    scalar maxKeff = 0.0;
    scalar voltageRiseRate  = 0.0;

    //  Dielectric relaxation (tau = epsilon / sigma)
    if (limitDielectricRelaxationRatio_ || printDielectricRelaxationRatio_)
    {
        tmp<volScalarField> tSigma = transport.electricalConductivity();
        const volScalarField& sigma = tSigma();
        
        maxSigma = gMax(mag(sigma)().primitiveField());
        scalar dielectricLimit = 
            (maxDielectricRelaxationRatio_ * eps0) / (maxSigma + VSMALL);

        if (limitDielectricRelaxationRatio_)
        {
            newDeltaT = min(newDeltaT, dielectricLimit);
        }
        tSigma.clear();
    }
    
    //  Species Courant Limit (convective and diffusive)
    if (limitSpeciesCo_ || printSpeciesCo_)
    {
        label sIdx = transport.species().speciesID(courantSpeciesName_);
        const volScalarField& n = transport.species().numberDensity(sIdx);
        const scalarField& nField = n.primitiveField();
        const scalarField& V = mesh_.V().field();

        // Convective
        const surfaceScalarField& convFlux = transport.convectiveFlux(sIdx);
        const scalarField sumConv
        (
            fvc::surfaceSum(mag(convFlux))().primitiveField()
        );
        const scalarField convRate
        (
            0.5 * sumConv / (nField * V + VSMALL)
        );
        maxConvFluxRate = gMax(convRate);
        if (limitSpeciesCo_)
        {
            newDeltaT = min
            (
                newDeltaT,
                maxSpeciesConvectiveCo_ / (maxConvFluxRate + VSMALL)
            );

            // The ENERGY rides the same flux at 5/3 the speed, so it reaches
            // the same Courant number at 3/5 of the step. Applied to the same
            // cap rather than a separate one: a user who has chosen a
            // convective Courant number has chosen it for the physics, and
            // the energy equation should honour that choice rather than need
            // a second key that can disagree with it.
            if (energyRate_ > 0)
            {
                newDeltaT = min
                (
                    newDeltaT,
                    maxSpeciesConvectiveCo_ / (energyRate_ + VSMALL)
                );
            }
        }

        // ENERGY RELAXATION LIMIT -- a SEPARATE constraint, deliberately
        // outside the limitSpeciesCo_ block and with its own coefficient.
        //
        // The convective limit above shares maxSpeciesConvectiveCo_ because a
        // Courant number chosen for the physics should apply to whatever rides
        // the same flux. This one must NOT: it bounds dt against the stiff
        // source relaxation time tau_eps = n_eps/P_loss, where a coefficient
        // of order 1 is not merely inaccurate but unstable. MEASURED
        // 2026-08-21: dt/tau = 3.1 oscillated between the 100 eV clamp and
        // 7 eV and died on SIGFPE; dt/tau = 0.15 settled smoothly with zero
        // clamps. Hence a default ratio of 0.2, not 1.0.
        if (limitEnergyRelaxation_ && energyRelaxRate_ > 0)
        {
            const scalar dtRelax =
                maxEnergyRelaxationRatio_ / (energyRelaxRate_ + VSMALL);

            // Report only when this limiter is the one that BINDS, so the log
            // says which constraint is setting the step rather than leaving
            // the user to infer it from a number that moved.
            if (dtRelax < newDeltaT && !energyRelaxWarned_)
            {
                energyRelaxWarned_ = true;
                Info<< "plasmaTimeControl: the ENERGY RELAXATION limit is now"
                    << " setting deltaT (tau_eps = "
                    << 1.0/(energyRelaxRate_ + VSMALL) << " s, ratio "
                    << maxEnergyRelaxationRatio_ << " -> " << dtRelax
                    << " s)." << nl
                    << "    The energy equilibrates faster than the transport"
                    << " Courant condition resolves. This is expected once the"
                    << " mean energy is high, where P_loss rises steeply."
                    << endl;
            }

            newDeltaT = min(newDeltaT, dtRelax);
        }

        // Diffusive
        const surfaceScalarField& diffFlux = transport.diffusiveFlux(sIdx);
        const scalarField sumDiff
        (
            fvc::surfaceSum(mag(diffFlux))().primitiveField()
        );
        const scalarField diffRate
        (
            0.5 * sumDiff / (nField * V + VSMALL)
        );
        maxDiffFluxRate = gMax(diffRate);

        if (limitSpeciesCo_)
        {
            newDeltaT = min
            (
                newDeltaT,
                maxSpeciesDiffusiveCo_ / (maxDiffFluxRate + VSMALL)
            );
        }
    }

    // Chemistry Limit
    if (limitChemistryCo_ || printChemistryCo_)
    {
        // Co_chem bounds how far the chemistry moves the state in one step.
        // It is NOT a stability limit: `implicitRate` is unconditionally
        // stable and the ODE path substeps under its own error control. What
        // it constrains is the SPLITTING error -- the field, and the rate
        // coefficients read from it, are frozen while the chemistry
        // integrates over dt. So it is a coupling constraint.
        //
        // Two measures, because neither alone is sufficient:
        //
        //  * k_eff is the NET electron rate (S_ion - S_att)/n_e. It goes to
        //    zero where ionisation and attachment nearly cancel, even though
        //    both processes are fast -- blind to the regime it should
        //    constrain.
        //  * chemL_ is the per-species loss coefficient the chemistry
        //    actually assembles, with no such cancellation. It is live in all
        //    chemistry modes (computeChemistrySources() sizes it for pure
        //    `ode` too), but it is 0 before the first chemistry evaluation and
        //    with `reactions none`.
        //
        // Taking the larger uses whichever is informative and degrades to
        // k_eff alone whenever P/L has not been assembled yet.
        //
        // NOT mag(k_eff): a NEGATIVE net rate is attachment-dominated decay in
        // the far field, which is benign and stable. Limiting dt on it
        // throttled the step for no accuracy gain.
        // Co_chem bounds how far the chemistry moves the STATE in one step.
        // It is NOT a stability limit: `implicitRate` is unconditionally
        // stable and the ODE path substeps under its own error control. What
        // it constrains is the SPLITTING error -- the field, and the rate
        // coefficients read from it, are frozen while the chemistry integrates
        // over dt. So it is a coupling constraint.
        //
        // The source is the per-species fractional net rate of change (see
        // plasmaTransport::maxChemStateRate), which is the mechanism-derived
        // generalisation of the legacy air-fitted k_eff: for electrons the two
        // are equal, so the 0.9 calibration carries over.
        //
        // k_eff remains the fallback for the legacy Townsend path, which
        // assembles no P/L. Growth only, never mag(): a negative net rate is
        // attachment-dominated decay, which is benign and stable.
        const scalar stateRate = transport.maxChemStateRate();

        if (stateRate > 0)
        {
            maxKeff = stateRate;
        }
        else
        {
            const volScalarField& keff = transport.k_eff();
            maxKeff = max(gMax(keff.primitiveField()), scalar(0));
        }

        if (limitChemistryCo_)
        {
            newDeltaT = min
            (
                newDeltaT,
                maxChemistryCo_ / (maxKeff + VSMALL)
            );
        }
    }

    // Voltage rise rate
    if (limitVoltageRiseRate_ || printVoltageRiseRate_)
    {
        const scalar currentVoltage = patchVoltageAvg(transport);

        // FIRST STEP: there is no previous voltage to difference against, and
        // seeding prevPatchVoltage_ with 0 makes dV the WHOLE applied voltage
        // -- kilovolts in one step. The inferred rate is then enormous and the
        // limiter pins deltaT at ~5e-15 s, which is a start-up artefact, not
        // physics. So the first pass only records the voltage.
        if (!voltageSeeded_)
        {
            voltageSeeded_ = true;
            prevPatchVoltage_ = currentVoltage;
        }

        const scalar dV = currentVoltage - prevPatchVoltage_;
        voltageRiseRate = dV / (runTime_.deltaT0Value() + VSMALL);

        if (limitVoltageRiseRate_)
        {
            const scalar dtLimit =
                maxVoltageRiseRate_ / (mag(voltageRiseRate) + VSMALL);

            newDeltaT = min(newDeltaT, dtLimit);
        }

        prevPatchVoltage_ = currentVoltage;
    }

    // OUTER-COUPLING STABILITY MARGIN, from the adaptive (Aitken) relaxation.
    //
    // Applied here with the other per-step limits so the smallest still wins.
    // It MUST live in adjustDeltaT and not in setInitialDeltaT: the latter
    // runs once, so a governor placed there is dead for the whole run -- which
    // is exactly the bug this replaces, caught because the reported margin sat
    // at 1 while the banner showed Aitken active with omegaMax = 0.8.
    //
    // Deliberately keyed on omega rather than on a non-convergence verdict: a
    // diverging step can grind indefinitely inside the stiff chemistry
    // integrator and never reach the corrector cap, so that verdict may never
    // arrive. omega is available every corrector.
    relaxMargin_ = transport.relaxationMargin();
    relaxOmegaMax_ = transport.relaxationOmegaMax();
    relaxMarginBound_ = false;
    if (relaxMargin_ < relaxMarginShrink_)
    {
        newDeltaT = min(newDeltaT, currentDeltaT*relaxMarginBackoff_);
        relaxMarginBound_ = true;
    }
    else if (relaxMargin_ < relaxMarginHold_)
    {
        newDeltaT = min(newDeltaT, currentDeltaT);
        relaxMarginBound_ = true;
    }

    // The outer loop failing to converge is a TIME-STEP verdict, so it is
    // applied here with the others and the smallest limit still wins. A
    // separate mechanism that set deltaT on its own could override the
    // dielectric or Courant limits and trade an accuracy problem for a
    // stability one.
    if (outerHitCap_ && outerOnNonConvergence_ == "reduceDeltaT")
    {
        // Halve. The coupling residual is not a smooth function of deltaT the
        // way a Courant number is -- there is no limit to solve for -- so this
        // backs off and lets the next step report whether that was enough.
        newDeltaT = min(newDeltaT, 0.5*currentDeltaT);
    }

    // Apply and clamp
    if (adjustTimeStep_)
    {
        reduce(newDeltaT, minOp<scalar>());

        if (newDeltaT > currentDeltaT)
        {
            // (2) A step that was just halved may not grow again immediately.
            if (noGrowthNextStep_)
            {
                newDeltaT = currentDeltaT;
            }
            else
            {
                newDeltaT = min(newDeltaT, currentDeltaT * 1.2);
            }
        }
        noGrowthNextStep_ = false;

        // (1) Rejection memory: stay below the scale that is known to fail,
        // and let that ceiling drift up only slowly. See dtFailCeiling_.
        if (dtFailCeiling_ > 0)
        {
            newDeltaT = min(newDeltaT, dtFailCeiling_);
            dtFailCeiling_ *= retryCeilingRelax_;
        }

        // (3) ABSOLUTE FLOOR, applied LAST so nothing can push below it. Note
        // this deliberately overrides the physical limiters too: if the
        // Courant or dielectric criterion genuinely wants a step this small
        // the run is not viable anyway, and stopping with a diagnostic serves
        // the user better than grinding towards a result that never arrives.
        if (newDeltaT < minDeltaT_)
        {
            ++floorHits_;
            newDeltaT = minDeltaT_;
        }

        runTime_.setDeltaT(newDeltaT);

        // (4) TREND. Warn once, while there is still time to act, when the
        // step has collapsed far below the largest this run has sustained.
        dtHighWater_ = max(dtHighWater_, newDeltaT);

        if
        (
            !collapseWarned_
         && dtCollapseRatio_ > 1
         && dtHighWater_ > 0
         && newDeltaT*dtCollapseRatio_ < dtHighWater_
        )
        {
            collapseWarned_ = true;

            WarningInFunction
                << "deltaT has COLLAPSED by more than " << dtCollapseRatio_
                << "x, from " << dtHighWater_ << " s to " << newDeltaT
                << " s." << nl
                << "    The run is losing ground faster than the 1.2x/step"
                << " regrowth can recover it, and may not reach endTime." << nl
                << "    The floor is " << minDeltaT_ << " s ("
                << (minDeltaTAuto_ ? "derived" : "user-set") << ")." << nl
                << "    Non-converged steps accepted so far: " << degradedRun_
                << ", floor hits: " << floorHits_ << nl
                << "    Watch the next few steps: if deltaT keeps falling,"
                << " stop and see the guidance printed when the floor is"
                << " reached." << endl;
        }
    }

    // Report
    const scalar actualDeltaT = runTime_.deltaTValue();

    word   bindingName  = "maxDeltaT";

    auto fmtLine = [](
        const std::string& label,
        scalar value,
        bool   limited,
        scalar limitVal,
        bool   isBinding,
        const std::string& maxLabel = "max",
        int    lw = 26) -> std::string
    {
        std::string line = "  " + label + ":";
        while (static_cast<int>(line.size()) < lw) line += ' ';
        line += Foam::name(value).c_str();

        if (limited)
        {
            while (static_cast<int>(line.size()) < lw + 14) line += ' ';
            line += "[ " + maxLabel + ": "
                  + std::string(Foam::name(limitVal).c_str()) + " ]";
        }

        if (isBinding) line += "  <--";
        return line;
    };

    Info<< "  Time step control" << nl
        << "  " << std::string(52, '-').c_str() << nl
        << "  current deltaT:           " << actualDeltaT << nl;

    if (limitDielectricRelaxationRatio_ || printDielectricRelaxationRatio_)
    {
        const scalar val = (actualDeltaT * maxSigma) / eps0;
        const bool lim  = limitDielectricRelaxationRatio_ && adjustTimeStep_;
        const bool binding =
            lim && val >= maxDielectricRelaxationRatio_ * 0.99;

        Info<< "  " << fmtLine(
                "Diel. relax. ratio",
                val, lim, maxDielectricRelaxationRatio_,
                binding
            ).c_str() << nl;
    }

    if (limitSpeciesCo_ || printSpeciesCo_)
    {
        const scalar convVal = maxConvFluxRate * actualDeltaT;
        const scalar diffVal = maxDiffFluxRate * actualDeltaT;
        const bool lim = limitSpeciesCo_ && adjustTimeStep_;
        const bool convBinding = lim && convVal >= 
                                                maxSpeciesConvectiveCo_ * 0.99;
        const bool diffBinding = lim && diffVal >= 
                                                 maxSpeciesDiffusiveCo_ * 0.99;

        Info<< "  " << fmtLine(
                "Co_conv (" + courantSpeciesName_ + ")",
                convVal, lim, maxSpeciesConvectiveCo_,
                convBinding
            ).c_str() << nl
            << "  " << fmtLine(
                "Co_diff (" + courantSpeciesName_ + ")",
                diffVal, lim, maxSpeciesDiffusiveCo_,
                diffBinding
            ).c_str() << nl;
    }

    if (limitChemistryCo_ || printChemistryCo_)
    {
        const scalar val = maxKeff * actualDeltaT;
        const bool lim  = limitChemistryCo_ && adjustTimeStep_;
        const bool binding = lim && val >= maxChemistryCo_ * 0.99;

        Info<< "  " << fmtLine(
                "Co_chem",
                val, lim, maxChemistryCo_,
                binding
            ).c_str() << nl;

        // DIAGNOSTIC, not a limiter input. The fastest single loss channel is
        // worth seeing -- it is what `chemStiffnessLimit` routes cells on --
        // but it must not steer the step: P ~ L n means the channels can be
        // enormous while the state barely moves, and limiting on it was
        // measured to double the step count for nothing. See
        // plasmaTransport::maxChemStateRate().
        const scalar Lpk = transport.maxChemLossRate();
        if (Lpk > 0)
        {
            Info<< "    (fastest loss channel:  " << Lpk
                << " 1/s, tau = " << 1.0/Lpk
                << " s -- diagnostic, does not limit dt)" << nl;
        }
    }

    if (limitVoltageRiseRate_ || printVoltageRiseRate_)
    {
        const scalar dVPerStep = voltageRiseRate * actualDeltaT;
        const bool lim  = limitVoltageRiseRate_ && adjustTimeStep_;
        const bool binding =
            lim && mag(dVPerStep) >= maxVoltageRiseRate_ * 0.99;

        Info<< "  " << fmtLine(
                "dV/step [V]",
                dVPerStep, lim, maxVoltageRiseRate_,
                binding, "abs max"
            ).c_str() << nl;
    }

    // ADAPTIVE-RELAXATION STABILITY MARGIN.
    //
    // Without this line there is no way to tell whether the Aitken factor is
    // doing anything at all -- the start-up banner only says it is enabled.
    // A margin of 1 means NOTHING was damped at this deltaT: either the
    // coupling is genuinely benign, or the global factor is being averaged
    // over a domain where only a few cells oscillate and cannot see them.
    // Those want opposite fixes, so the number has to be visible.
    if (relaxMargin_ < 1.0 - SMALL || relaxMarginBound_)
    {
        Info<< "    omega [coupling margin]:  min " << relaxMargin_
            << "  max " << relaxOmegaMax_
            << (relaxMarginBound_ ? "   <--  deltaT GOVERNED" : "")
            << nl;
    }
    else
    {
        Info<< "    omega [coupling margin]:  min 1  (nothing damped)" << nl;
    }

    // OUTER-LOOP HALVING, reported beside the other constraints.
    //
    // Without this line the report is actively misleading: it shows every
    // Courant number comfortably below its cap and gives no reason why deltaT
    // is not larger. The reason is here -- the coupling did not converge, so
    // adjustDeltaT() halved the step, and that limit wins over the Courant
    // ones. Working this out from the outside took a long detour, from
    // "Co_conv 1.43 [max 2.5]" to the conclusion that the controller was
    // misbehaving, when it was doing exactly what it was told.
    //
    // The corrector count is printed too, because it says whether raising
    // maxCorrectors would help: a loop that used every one of them and was
    // still moving is slow, one that stalls is diverging, and the two want
    // opposite responses.
    if (outerHitCap_)
    {
        Info<< "  outer loop:               NOT converged in "
            << outerItersUsed_ << " correctors (cap "
            << outerMaxCorrectors_ << ")";

        if (outerOnNonConvergence_ == "reduceDeltaT" && adjustTimeStep_)
        {
            Info<< "  <--  deltaT HALVED" << nl
                << "                            step is not second order;"
                << " raise outerCoupling/maxCorrectors if this persists";
        }
        Info<< nl;
    }

    // TREND, per report interval. A run ratcheting down shows up here long
    // before it hits the floor -- a cluster of degraded steps can lose ground
    // far faster than 1.2x-per-step regrowth recovers it.
    if (degradedSteps_ > 0 || floorHits_ > 0)
    {
        Info<< "  non-converged accepted:   " << degradedSteps_
            << " this interval, " << degradedRun_ << " total";

        if (floorHits_ > 0)
        {
            Info<< nl << "  deltaT floor:             hit " << floorHits_
                << " times (minDeltaT " << minDeltaT_ << " s, "
                << (minDeltaTAuto_ ? "derived" : "user-set") << ")";
        }
        Info<< nl;
    }
    degradedSteps_ = 0;

    Info<< "  " << std::string(52, '-').c_str() << nl << endl;
}

void plasmaTimeControl::setInitialDeltaT(const plasmaTransport& transport)
{
    if (!adjustTimeStep_) return;

    scalar newDeltaT = maxDeltaT_;
    const scalar currentDeltaT = runTime_.deltaTValue();
    const scalar eps0 = constant::plasma::epsilon0.value();

    // Dielectric relaxation (tau = epsilon / sigma)
    if (limitDielectricRelaxationRatio_)
    {
        tmp<volScalarField> tSigma = transport.electricalConductivity();
        const volScalarField& sigma = tSigma();
        const scalar maxSigma = gMax(mag(sigma)().primitiveField());
        newDeltaT = min
        (
            newDeltaT,
            (maxDielectricRelaxationRatio_ * eps0) / (maxSigma + VSMALL)
        );
        tSigma.clear();
    }

    // Species Courant Limit (convective and diffusive)
    if (limitSpeciesCo_)
    {
        const label sIdx = transport.species().speciesID(courantSpeciesName_);
        const volScalarField& n = transport.species().numberDensity(sIdx);
        const scalarField& nField = n.primitiveField();
        const scalarField& V = mesh_.V().field();

        // Convective
        {
            const surfaceScalarField& convFlux =
                transport.convectiveFlux(sIdx);
            const scalarField sumConv
            (
                fvc::surfaceSum(mag(convFlux))().primitiveField()
            );
            const scalar maxRate =
                gMax(0.5 * sumConv / (nField * V + VSMALL));
            newDeltaT = min
            (
                newDeltaT,
                maxSpeciesConvectiveCo_ / (maxRate + VSMALL)
            );
        }

        // Diffusive
        {
            const surfaceScalarField& diffFlux =
                transport.diffusiveFlux(sIdx);
            const scalarField sumDiff
            (
                fvc::surfaceSum(mag(diffFlux))().primitiveField()
            );
            const scalar maxRate =
                gMax(0.5 * sumDiff / (nField * V + VSMALL));
            newDeltaT = min
            (
                newDeltaT,
                maxSpeciesDiffusiveCo_ / (maxRate + VSMALL)
            );
        }
    }

    // Chemistry Courant Limit
    if (limitChemistryCo_)
    {
        // Growth only, not mag() -- see adjustDeltaT(). chemL_ is necessarily
        // empty here (the chemistry has not run yet), so this is the k_eff
        // measure alone by construction.
        const volScalarField& keff = transport.k_eff();
        const scalar maxKeff = max(gMax(keff.primitiveField()), scalar(0));
        newDeltaT = min
        (
            newDeltaT,
            maxChemistryCo_ / (maxKeff + VSMALL)
        );
    }

    // Voltage rise rate
    if (limitVoltageRiseRate_ || printVoltageRiseRate_)
    {
        prevPatchVoltage_ = patchVoltageAvg(transport);
    }

    // Reduction and setting
    reduce(newDeltaT, minOp<scalar>());
    if (newDeltaT > currentDeltaT)
    {
        newDeltaT = min(newDeltaT, currentDeltaT * 1.2);
    }
    runTime_.setDeltaT(newDeltaT);   
}

void plasmaTimeControl::configureOuterCoupling(fvMesh& mesh)
{
    // Read from the same file as the rest of the plasma controls, so a user
    // does not have to know that this one happens to be implemented through
    // PIMPLE.
    IOdictionary controls
    (
        IOobject
        (
            "plasmaSimulationControls",
            mesh.time().system(),
            mesh.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    const dictionary& oc = controls.subOrEmptyDict("outerCoupling");

    const word target = oc.getOrDefault<word>("target", "converged");
    if (target != "converged" && target != "lagged")
    {
        FatalErrorInFunction
            << "Unknown outerCoupling/target '" << target << "'" << nl
            << "Valid: converged | lagged" << nl
            << exit(FatalError);
    }

    dictionary& pimpleDict = mesh.solution().subDict("PIMPLE");

    if (target == "lagged")
    {
        // A DELIBERATE first-order mode, and it says so. The cost of one
        // corrector is not just a looser tolerance: the Poisson and species
        // equations are lagged relative to each other and the scheme drops to
        // first order in time whatever ddtSchemes says. Measured on the
        // shipped benchmark: p = 0.87 at one corrector, p = 1.94 converged.
        pimpleDict.set("nOuterCorrectors", 1);
        Info<< "plasmaTimeControl: outerCoupling target `lagged`" << nl
            << "    ONE outer corrector. The Poisson-species coupling is not"
            << " converged, so the" << nl
            << "    scheme is FIRST ORDER in time regardless of ddtSchemes."
            << " Use for scoping runs." << endl;
        return;
    }

    const label maxCorr = oc.getOrDefault<label>("maxCorrectors", 20);
    const scalar tol = oc.getOrDefault<scalar>("tolerance", 1e-8);

    if (maxCorr < 2)
    {
        FatalErrorInFunction
            << "outerCoupling/maxCorrectors is " << maxCorr
            << ", so the loop cannot iterate." << nl
            << "    `target converged` needs room to converge; use"
            << " `target lagged` if one corrector is what you want." << nl
            << exit(FatalError);
    }

    // THE TOLERANCE MUST BE LOOSER THAN THE LINEAR SOLVER'S. Otherwise the
    // outer criterion asks for something the inner solve never delivers, the
    // loop runs to its cap on every step, and `nOuterCorrectors` silently
    // becomes a fixed iteration count again. That is exactly the state the
    // shipped streamer case was in: outer tolerance 1e-10 against a linear
    // tolerance of 1e-10, 4 of 4 correctors used on all 3000 steps.
    scalar worstInner = 0;
    const dictionary& solvers = mesh.solution().subDict("solvers");
    for (const entry& e : solvers)
    {
        if (!e.isDict()) continue;
        worstInner =
            max(worstInner, e.dict().getOrDefault<scalar>("tolerance", 0));
    }

    if (worstInner > 0 && tol <= 10*worstInner)
    {
        FatalErrorInFunction
            << "outerCoupling/tolerance (" << tol << ") is not comfortably"
            << " above the linear-solver" << nl
            << "    tolerance (" << worstInner << ")." << nl << nl
            << "    The outer loop measures how much the coupled solution"
            << " still moves between" << nl
            << "    correctors. It cannot resolve a change smaller than the"
            << " inner solves' own" << nl
            << "    convergence, so the criterion would never be met and the"
            << " loop would run to" << nl
            << "    maxCorrectors every step -- a fixed iteration count"
            << " wearing the name of a" << nl
            << "    convergence test." << nl << nl
            << "    Use at least 10x the linear-solver tolerance, i.e. >= "
            << 10*worstInner << "." << nl
            << exit(FatalError);
    }

    pimpleDict.set("nOuterCorrectors", maxCorr);

    // `outerCoupling` is the single user-facing knob, so it OWNS
    // residualControl and replaces whatever the case wrote. Say so out loud:
    // a setting that is silently ignored is its own trap.
    {
        const dictionary& prior = pimpleDict.subOrEmptyDict("residualControl");
        if (!prior.empty())
        {
            Info<< "plasmaTimeControl: PIMPLE/residualControl in the case"
                << " (" << flatOutput(prior.toc()) << ")" << nl
                << "    is SUPERSEDED by the outerCoupling block and can be"
                << " deleted from fvSolution." << endl;
        }
    }

    // WHICH FIELDS GATE THE LOOP -- and why the default is ePotential ALONE.
    //
    // This block used to build the pattern `"n_.*"` through a `word` that kept
    // its own quote characters: `word::stripInvalid()` is DEBUG-GATED
    // (etc/controlDict ships `word 0;`), so the quotes survived and
    // solutionControl compiled a pattern that could never match a field called
    // `n_e`. For the whole life of that code the species half of this
    // criterion was DEAD, and the outer loop exited on ePotential alone while
    // n_e was still at 1e-2.
    //
    // The patterns below are therefore built as explicit REGEX keyTypes. But
    // the DEFAULT deliberately remains `(ePotential)`, i.e. exactly the
    // behaviour that was actually in force, because that is what every
    // validated result in this project was produced under -- the shipped
    // benchmark's p = 1.936 order measurement included. Silently switching the
    // default to "all species and the energy" changes the accept/reject
    // verdict of every existing case, and it was MEASURED to make the LMEA
    // case degrade 3 of its first 4 steps (2026-08-21, `co-lmea-o4-rcfix`).
    //
    // It degrades for a reason worth recording before anyone raises the list:
    // OpenFOAM's normFactor is built from the deviation of the field about its
    // own average, so for a NEARLY UNIFORM field -- which n_e is before the
    // streamer develops -- it collapses and the normalised residual saturates
    // near 1 however well converged the coupling is. An absolute 1e-8 on n_e
    // is then unreachable in principle, not merely in practice.
    //
    // So this is a knob, not a new default. Widen it per case, knowing that.
    const wordList gateFields
    (
        oc.getOrDefault<wordList>("gateFields", wordList({"ePotential"}))
    );

    dictionary rc;

    // Resolve every requested pattern against the fields actually registered
    // NOW and report the mapping. A pattern matching nothing is exactly the
    // failure that hid here for so long, so it is an ERROR: a criterion that
    // cannot fire is decoration, not a check.
    const wordList regFields(mesh.sortedNames<volScalarField>());

    label nMatched = 0;

    Info<< "plasmaTimeControl: outerCoupling target `converged`" << nl
        << "    residual " << tol << ", up to " << maxCorr
        << " correctors (a CAP -- the loop exits when converged)" << nl
        << "    fields gating the loop:" << endl;

    for (const word& req : gateFields)
    {
        // Anything carrying regex metacharacters is a pattern; a plain name is
        // matched literally. Built here rather than passed through a quoted
        // `word`, which is what broke before.
        const bool isPattern =
            (req.find_first_of(".*[]()^$+?|\\") != std::string::npos);

        const keyType key
        (
            req,
            isPattern ? keyType::REGEX : keyType::LITERAL
        );

        DynamicList<word> matched;
        for (const word& f : regFields)
        {
            if (key.match(f))
            {
                matched.append(f);
            }
        }

        if (matched.empty())
        {
            FatalErrorInFunction
                << "outerCoupling/gateFields lists `" << req
                << "`, which matches NO solved field" << nl
                << "    on mesh `" << mesh.name() << "`." << nl << nl
                << "    A convergence criterion that cannot fire is not a"
                << " criterion: the outer loop" << nl
                << "    would exit on the remaining fields alone and report a"
                << " diverging solution" << nl
                << "    as converged. That is a real bug this project has"
                << " already been bitten by." << nl << nl
                << "    Registered volScalarFields: "
                << flatOutput(regFields) << nl
                << exit(FatalError);
        }

        dictionary fd;
        fd.set("tolerance", tol);
        fd.set("relTol", scalar(0));
        rc.set(key, fd);

        nMatched += matched.size();

        Info<< "    " << req << " -> " << flatOutput(matched) << endl;
    }

    pimpleDict.set("residualControl", rc);

    Info<< "    " << nMatched << " field(s) must meet " << tol
        << " for the step to be accepted" << endl;
}


bool plasmaTimeControl::atMinDeltaT() const
{
    // A retry costs a full solve, so it is only worth spending one if the
    // step can actually be made meaningfully smaller.
    return runTime_.deltaTValue() <= minDeltaT_*(1 + SMALL);
}


void plasmaTimeControl::noteDegradedStep()
{
    ++degradedSteps_;
    ++degradedRun_;
    ++degradedConsecutive_;

    WarningInFunction
        << "step accepted WITHOUT the outer loop converging, at t = "
        << runTime_.timeName() << ", deltaT = " << runTime_.deltaTValue()
        << " s." << nl
        << "    The time-step lever is exhausted ("
        << (atMinDeltaT() ? "deltaT is at the floor" : "retries spent")
        << "), so the best available iterate was kept." << nl
        << "    This step is NOT converged to residualControl -- treat"
        << " results from it as suspect. Consecutive: "
        << degradedConsecutive_ << " of " << maxDegradedConsecutive_
        << " allowed." << endl;

    // Isolated degradation is a fair trade: the run finishes and the user
    // knows which steps to distrust. SYSTEMATIC degradation is not -- it means
    // the solver has driven deltaT to the floor and is still failing, so
    // every subsequent result is meaningless and continuing only burns CPU.
    if (degradedConsecutive_ >= maxDegradedConsecutive_)
    {
        FatalErrorInFunction
            << nl
            << "==============================================" << nl
            << " The solver cannot converge this case." << nl
            << "==============================================" << nl << nl
            << degradedConsecutive_ << " consecutive time steps were accepted"
            << " without the outer (PIMPLE) loop reaching its residual"
            << " tolerance, with deltaT already at " << runTime_.deltaTValue()
            << " s" << nl
            << "(floor " << minDeltaT_ << " s, "
            << (minDeltaTAuto_ ? "derived" : "set by you")
            << "). Reducing the time step is the solver's last automatic"
            << " lever and it is used up, so it is stopping rather than"
            << " producing results that are not solutions." << nl << nl
            << "State at the failure:" << nl
            << "    time                    " << runTime_.timeName() << nl
            << "    deltaT                  " << runTime_.deltaTValue()
            << " s" << nl
            << "    largest deltaT sustained " << dtHighWater_ << " s" << nl
            << "    correctors on last step " << outerItersUsed_
            << " (cap " << outerMaxCorrectors_ << ")" << nl
            << "    retries spent this step " << outerRetries_
            << " of " << outerMaxRetries_ << nl
            << "    non-converged steps     " << degradedRun_
            << " over the run" << nl
            << "    floor hits              " << floorHits_ << nl << nl
            << "What to check, in the order most likely to help:" << nl << nl
            << "  1. MESH QUALITY. Run checkMesh. Non-orthogonality above ~70"
            << " or high skewness makes the Poisson and drift-diffusion"
            << " solves refuse to converge no matter how small the step is."
            << " This is the most common cause and the solver cannot fix it."
            << nl << nl
            << "  2. BOUNDARY CONDITIONS. An over-specified or inconsistent"
            << " set (for example a fixed potential facing a fixed charge"
            << " flux) has no steady solution for the outer loop to find."
            << nl << nl
            << "  3. NEGATIVE OR RUNAWAY DENSITIES. Check the chemistry report"
            << " above for negative-density clips and ODE failures. A species"
            << " going negative makes rate coefficients meaningless and the"
            << " coupling diverge." << nl << nl
            << "  4. TIME-STEP CRITERIA. Lower maxSpeciesCo (try halving it)"
            << " so the step is limited by physics before the coupling"
            << " fails, rather than after." << nl << nl
            << "  5. FIELD/CHARGE COUPLING. If this is a fast-ionisation"
            << " transient, check maxDielectricRelaxationRatio -- the"
            << " semi-implicit Poisson needs it near or below 1 when space"
            << " charge is growing quickly." << nl << nl
            << "  6. Only if the above are all sound: raise maxCorrectors, or"
            << " loosen outerCoupling/residual. Note the outer loop here is"
            << " bimodal -- it converges in a few correctors or not at all --"
            << " so raising the cap rarely helps." << nl << nl
            << "To let the run continue anyway and inspect the damage, set"
            << " plasmaTimeControl/maxConsecutiveDegradedSteps to a larger"
            << " value. The results from those steps are not converged."
            << nl
            << exit(FatalError);
    }
}


bool plasmaTimeControl::stepRejected() const
{
    if (outerOnNonConvergence_ != "retryStep" || !outerHitCap_) return false;

    // Halving a step that is already on the floor buys nothing but a wasted
    // solve; let the caller accept and degrade instead.
    if (atMinDeltaT()) return false;

    // The counter belongs to a time index, so a new step resets it without the
    // solver having to remember to.
    const label ti = runTime_.timeIndex();
    if (ti != outerRetryTimeIndex_) return true;

    return outerRetries_ < outerMaxRetries_;
}


void plasmaTimeControl::prepareRetry()
{
    const label ti = runTime_.timeIndex();
    if (ti != outerRetryTimeIndex_)
    {
        outerRetryTimeIndex_ = ti;
        outerRetries_ = 0;
    }
    ++outerRetries_;

    const scalar dtOld = runTime_.deltaTValue();

    // Never descend through the floor: stepRejected() already refuses to
    // retry AT it, and this keeps the last admissible rung on it rather than
    // below it.
    const scalar dtNew = max(0.5*dtOld, minDeltaT_);

    if (dtNew > 0.5*dtOld)
    {
        ++floorHits_;
    }

    // Remember the scale that failed, so the controller stops climbing back
    // into it every few steps, and forbid growth on the next step.
    if (retryCeilingFactor_ > 0)
    {
        const scalar ceiling = retryCeilingFactor_*dtOld;
        dtFailCeiling_ = (dtFailCeiling_ > 0)
                       ? min(dtFailCeiling_, ceiling)
                       : ceiling;
    }
    noGrowthNextStep_ = true;

    // Land the clock on the shortened step, measured from the CAPTURED start
    // of this step rather than from anything reconstructed out of the current
    // time and deltaT -- see tStepStart_. setTime keeps the INDEX, which is
    // what every per-step cache and oldTime() rotation is keyed on.
    runTime_.setDeltaT(dtNew);
    runTime_.setTime(tStepStart_ + dtNew, ti);

    Info<< "  outer loop did not converge -- DISCARDING this step and"
        << " re-running it" << nl
        << "    retry " << outerRetries_ << " of " << outerMaxRetries_
        << ", deltaT " << dtOld << " -> " << dtNew
        << ", t = " << runTime_.value() << endl;

    // The step is being redone, so the verdict from the attempt just thrown
    // away must not also shrink the NEXT step in adjustDeltaT().
    outerHitCap_ = false;
}


void plasmaTimeControl::noteOuterLoop(const label used)
{
    outerItersUsed_ = max(outerItersUsed_, used);

    // Converged == exited before the cap. That is OpenFOAM's own
    // residualControl verdict; re-deriving it here from residuals would risk a
    // second criterion that disagrees with the one steering the loop.
    const bool hitCap = outerChaseConvergence_ && (used >= outerMaxCorrectors_);
    outerHitCap_ = hitCap;

    // A converged step breaks the run of degraded ones. Only a CONSECUTIVE
    // run is evidence that the case itself cannot be solved -- isolated
    // failures around a fast transient are expected and recoverable.
    if (!hitCap)
    {
        degradedConsecutive_ = 0;
    }

    if (hitCap)
    {
        if (outerOnNonConvergence_ == "fatal")
        {
            FatalErrorInFunction
                << "The outer coupling did not converge in "
                << outerMaxCorrectors_ << " correctors." << nl
                << "    The Poisson-species coupling is unconverged, so this"
                << " step is first order at" << nl
                << "    best. Raise maxCorrectors, loosen tolerance, or reduce"
                << " deltaT." << nl
                << exit(FatalError);
        }

        // `reduceDeltaT` is handled in adjustDeltaT(), where every other
        // time-step criterion lives, so the smallest of them wins rather than
        // this one overriding the Courant or dielectric limits.
        WarningInFunction
            << "outer coupling hit the " << outerMaxCorrectors_
            << "-corrector cap without converging; this step is not"
            << " second order." << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
