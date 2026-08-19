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
        const volScalarField& keff = transport.k_eff();
        maxKeff = gMax(mag(keff)().primitiveField());
        
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

        runTime_.setDeltaT(newDeltaT);
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
        const volScalarField& keff = transport.k_eff();
        const scalar maxKeff = gMax(mag(keff)().primitiveField());
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

    // Applied to ePotential and to every transported species, which together
    // ARE the coupling: the lag is between the field and the densities.
    dictionary rc;
    for (const word& f : wordList({"ePotential", "\"n_.*\""}))
    {
        dictionary fd;
        fd.set("tolerance", tol);
        fd.set("relTol", scalar(0));
        rc.set(f, fd);
    }
    pimpleDict.set("residualControl", rc);

    Info<< "plasmaTimeControl: outerCoupling target `converged`" << nl
        << "    residual " << tol << ", up to " << maxCorr
        << " correctors (a CAP -- the loop exits when converged)" << endl;
}


bool plasmaTimeControl::stepRejected() const
{
    if (outerOnNonConvergence_ != "retryStep" || !outerHitCap_) return false;

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
    const scalar dtNew = 0.5*dtOld;

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
