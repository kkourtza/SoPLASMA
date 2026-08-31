/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026 Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Application
    soPlasmaFoam

Description
    Transient solver for coupled gas (plasma) and dielectric domains developed
    for plasma simulation purposes.

Usage
    \b soPlasmaFoam [OPTIONS]

    Example:
        soPlasmaFoam -case testCase

Author
    Rention Pasolari
    Contact: r.pasolari@gmail.com
\*---------------------------------------------------------------------------*/

#include <iomanip>

#include "IFstream.H"
#include "plasmaOuterRelaxation.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "regionProperties.H"
#include "fvOptions.H"
#include "loopControl.H"
#include "fvSolution.H"
#include "solutionControl.H"
#include "mappedPatchBase.H"
#include "pimpleControl.H"

#include "plasmaConstants.H"
#include "plasmaTimeControl.H"
#include "electromagneticsModel.H"
#include "multiRegionPoisson.H"
#include "dielectricRegion.H"
#include "plasmaSpecies.H"
#include "plasmaBoltzmann.H"
#include "plasmaTransport.H"
#include "plasmaEnergy.H"
#include "driftDiffusion.H"
#include "plasmaSimulationDiagnostics.H"
#include "plasmaDischargeCurrent.H"
#include "plasmaSimulationProfiler.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for coupled gas (plasma) and dielectric domains"
        " developed for plasma simulation purposes."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"


// //- TEMPORARY
//     // MUST_READ: this is the prescribed Gaussian from 0/S_iz
//     volScalarField ionizationSourceField
//     (
//         IOobject
//         (
//             "ionizationSourceField",
//             runTime.timeName(),
//             gasMesh(),
//             IOobject::MUST_READ,
//             IOobject::AUTO_WRITE
//         ),
//         gasMesh()
//     );



    //- Create the electromagnetics model
    autoPtr<electromagneticsModel> em =
        electromagneticsModel::New(gasMesh(), dielectricFvMeshes);

    //- Solve the EEDF and write the rate/transport tables, BEFORE anything
    //  reads them.
    //
    //  Ordering matters and is not obvious: the species transport models are
    //  built inside plasmaSpecies, and a `fromMechanism` mobility looks for its
    //  table in its constructor. plasmaReactionRates -- which owns the
    //  mechanism -- is not built until plasmaTransport, several lines later. So
    //  the sweep has to be kicked off here, from the dictionary, rather than
    //  from the class that happens to own the mechanism afterwards.
    //
    //  Idempotent: plasmaReactionRates calls the same thing and finds the
    //  tables already current.
    {
        IOdictionary transportDict
        (
            IOobject
            (
                "plasmaTransportProperties",
                runTime.constant(),
                gasMesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );

        if (transportDict.found("chemistry"))
        {
            const dictionary& chem = transportDict.subDict("chemistry");
            const fileName mechFile = chem.get<fileName>("mechanism");

            // THE MECHANISM HASH, read from the mechanism dictionary.
            //
            // This used to pass `word::null`, and plasmaBoltzmann::ensureTables
            // only reuses an existing table set when
            // `!expectedHash.empty() && head.find(expectedHash)`. An empty
            // hash therefore made the reuse branch UNREACHABLE: this call
            // re-solved the entire Boltzmann sweep on EVERY launch -- 261
            // points, minutes of wall clock, for every user of the solver --
            // and then plasmaReactionRates, which does pass the real hash,
            // found the freshly written tables and reused them.
            //
            // Read the same way plasmaReactionRates reads it
            // (plasmaReactionRates::readMechanism): the hash may begin with a
            // digit, which an OpenFOAM word may not, so mechc quotes it and it
            // must be read as a string.
            word mechHash = word::null;
            {
                IFstream mis(mechFile);
                if (mis.good())
                {
                    dictionary mdict(mis);
                    if (mdict.found("mechanismHash"))
                    {
                        mechHash =
                            word(string(mdict.get<string>("mechanismHash")));
                    }
                }
                // A missing hash is NOT fatal here: it simply restores the old
                // always-rebuild behaviour, which is correct if slow.
                // plasmaReactionRates fails loudly if the mechanism is
                // unreadable, so this does not mask a broken case.
            }

            plasmaBoltzmann::ensureTables
            (
                chem,
                chem.getOrDefault<fileName>
                (
                    "manifest", mechFile.lessExt() + ".mech.json"
                ),
                chem.getOrDefault<fileName>
                (
                    "tableDir", "constant/plasmaTables"
                ),
                mechHash
            );
        }
    }

    //- Create the plasmaSpecies model
    plasmaSpecies species(gasMesh(), em());

    // Gas energy. Constructed BEFORE plasmaTransport and registered on the
    // mesh, so the transport model can find it by lookup: plasmaEnergy owns
    // T_gas and the energy equation, plasmaTransport supplies the heat source.
    //
    // ONLY WHEN THE CASE ASKS FOR IT, on design grounds: constructing it
    // builds a per-species energy model for every species, and a feature that
    // is switched off should cost nothing, construction included.
    //
    // Not for the reason first written here. A segfault at the first timestep
    // was blamed on unconditional construction; the actual cause was an ABI
    // mismatch -- plasmaTransport.H gained members, and libplasmaTools, which
    // takes a const plasmaTransport&, was still compiled against the old
    // layout after a partial rebuild. Changing a public header's member layout
    // invalidates every library that links it, not only the one that defines
    // it.
    autoPtr<plasmaEnergy> energy;
    {
        // Ask plasmaEnergy whether it is needed, rather than re-deriving the
        // condition here. Gas heating and the per-species energy models are
        // INDEPENDENT physics: a case may want LMEA without solving the gas
        // temperature, or gas heating without LMEA. Testing only
        // `energy { solve }` here meant enabling LMEA also switched on the
        // gas-temperature equation, which then demanded a conduction scheme
        // and a T_gas linear solver for an equation the case never asked for.
        if (plasmaEnergy::required(species))
        {
            energy.reset(new plasmaEnergy(species, gasMesh(), em().E()));
        }
    }

    //- Translate the `outerCoupling` block into PIMPLE's nOuterCorrectors and
    //  residualControl.
    //
    //  BEFORE plasmaTransport, not merely before pimpleControl: the chemistry
    //  reads nOuterCorrectors at construction to check that `implicitRate` and
    //  `adaptive` have an outer loop to converge. Configuring afterwards would
    //  make that check fire on the case's raw value while the run went on to
    //  use the overridden one -- a guard reading a number nothing used.
    plasmaTimeControl::configureOuterCoupling(gasMesh());

    //- Create the plasmaTransport model
    plasmaTransport transport(gasMesh(), species);

    //- Create the PIMPLE loop control
    pimpleControl pimple(gasMesh());

    //- Create the timeControl manager and set initial time-step
    plasmaTimeControl timeControl(runTime, gasMesh());

    // BEFORE setInitialDeltaT: an LMEA energy equation convects at 5/3 the
    // particle rates, so the very first step must already respect its Courant
    // number. Declaring it afterwards would leave step 1 unguarded, which is
    // exactly where the energy was measured going negative.
    if (energy)
    {
        timeControl.setEnergyRate(energy->maxEnergyRate());
        timeControl.setEnergyRelaxRate(energy->maxEnergyRelaxationRate());
    }

    timeControl.setInitialDeltaT(transport);

    //- Create the plasmaSimulationDiagnostics manager
    plasmaSimulationDiagnostics diagnostics(runTime, transport);

    //- External circuit current by Sato's equation. Default OFF, so no
    //  existing case changes behaviour or cost. Its geometric weighting field
    //  is solved ONCE here, at construction.
    IOdictionary plasmaControlsDict
    (
        IOobject
        (
            "plasmaSimulationControls",
            runTime.system(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    plasmaDischargeCurrent dischargeCurrent
    (
        gasMesh(), plasmaControlsDict, em()
    );

    #include "reportSimulationSummary.H"

    runTime.writeNow();

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.run())
    {
        // BEFORE the clock moves: the only point at which the time this step
        // starts from is unambiguous. A retry needs it to place the shortened
        // step, and it cannot be reconstructed afterwards -- ++runTime uses
        // the previous deltaT while adjustDeltaT then installs a new one.
        // The energy rate follows the solution, so refresh it each step --
        // unlike the old constant 5/3 factor, which never changed.
        if (energy) timeControl.setEnergyRate(energy->maxEnergyRate());
        if (energy)
        {
            timeControl.setEnergyRelaxRate(energy->maxEnergyRelaxationRate());
        }

        timeControl.noteStepStart();

        // CHOOSE deltaT BEFORE ADVANCING THE CLOCK.
        //
        // It used to be chosen after ++runTime, which meant the clock advanced
        // by the PREVIOUS step's deltaT while fvm::ddt then discretised over
        // the newly chosen one. Under the 1.2x ramp those differ by 20%, so
        // the state written at `Time = t` had actually been advanced over a
        // different interval than the label implies. It does not accumulate --
        // the two running sums telescope, so the mismatch stays bounded at
        // about one step -- but the label and the discretisation should agree,
        // and a retry that has to place a shortened step cannot reason about a
        // clock that moved by a different amount than it is about to solve.
        //
        // adjustDeltaT reads only the transport state, which the advance does
        // not touch, so moving it earlier changes nothing it depends on.
        timeControl.adjustDeltaT(transport);

        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;
        gasMesh().update();

        // RETRY LOOP.
        //
        // Under `outerCoupling/onNonConvergence retryStep` a step whose outer
        // loop did not converge is DISCARDED rather than accepted: the fields
        // go back to oldTime(), deltaT is halved, and the step is solved
        // again. Under every other setting stepRejected() is false and this
        // loop runs exactly once, so existing cases are untouched.
        //
        // The time INDEX is held fixed across retries. That is the whole
        // trick: oldTime() rotates on the index, and so do the solver's
        // per-step caches, so a retry starts from precisely the state the
        // first attempt did -- with no snapshot buffers and no rewinding of
        // the BDF2 history.
        for (bool stepDone = false; !stepDone; /*retry*/)
        {
        // Correctors executed this step. The loop exiting before the cap is
        // OpenFOAM's own residualControl verdict that the Poisson-species
        // coupling has converged -- which is what makes the step second order.
        label nOuter = 0;

        // Semi-implicit Poisson branch
        if (em->PoissonScheme() == "semiImplicit")
        {
            while (pimple.loop())
            {
                ++nOuter;

                // Solve electromagnetics
                em->solve
                (
                    transport.electricalConductivity(), 
                    transport.diffusiveChargeSource()
                );

                // Solve transport equations
                transport.solve(pimple.finalIter());

                // Electron energy (LMEA), AFTER the species: the conservative
                // form recovers eps_bar = n_eps/n_e, so it needs n_e^{k+1}.
                // That is Hagelaar & Kroesen's ordering requirement, and it is
                // why this sits here rather than beside the Poisson solve.
                // A no-op for the LFA family, whose eEqn() returns nullptr.
                if (energy) energy->solveSpeciesEnergy();

                if (pimple.finalIter())
                {
                    // Update charge density
                    species.updateChargeDensity();

                    // Update surface charge
                    transport.updateSurfaceCharge();
                }
            }
        }
        else //Explicit Poisson branch
        {
            while (pimple.loop())
            {
                ++nOuter;

                plasmaSimulationProfiler::start("Electromagnetics");
                em->solve();
                plasmaSimulationProfiler::stop("Electromagnetics");

                // Solve transport equations
                plasmaSimulationProfiler::start("Plasma Transport");
                transport.solve(pimple.finalIter());
                plasmaSimulationProfiler::stop("Plasma Transport");

                // Electron energy (LMEA); see the semi-implicit branch above.
                if (energy) energy->solveSpeciesEnergy();

                // Update charge density
                plasmaSimulationProfiler::start("Update charge density");
                species.updateChargeDensity();
                plasmaSimulationProfiler::stop("Update charge density");

                // Update surface charge
                plasmaSimulationProfiler::start("Update surface charge");
                transport.updateSurfaceCharge();
                plasmaSimulationProfiler::stop("Update surface charge");


            }
        }

        // AFTER the loop, BEFORE the next adjustDeltaT: a step that exhausted
        // the cap is not second order, and the response is a smaller step.
        timeControl.noteOuterLoop(nOuter);

        if (timeControl.stepRejected())
        {
            // ================= THE STEP-DISCARD INVARIANT =================
            //
            // EVERY piece of solution state must be idempotent WITHIN a
            // timestep: this step is about to be thrown away and re-run, and
            // the answer must not depend on the attempt being discarded.
            //
            // Anything added to this solver in future must satisfy it one of
            // three ways:
            //   * recomputed from oldTime() each pass, and restored here on
            //     rejection -- REQUIRED for anything carrying fvm::ddt
            //     (species densities, T_gas, n_eps today);
            //   * baselined so a re-run recomputes instead of accumulating
            //     (e_vib does `ev = ev_start + ...`, never `ev += ...`;
            //     surfCharge resets from its own oldTime before accumulating);
            //   * keyed on timeIndex, for caches.
            // Per-pass counters must be ZEROED, or the discarded attempt is
            // counted alongside the kept one.
            //
            // THIS IS NOT THEORETICAL. n_eps was missed when retryStep was
            // written -- it lives in plasmaEnergy, not plasmaTransport -- so a
            // retried LMEA step re-solved with n_e at t^n and n_eps at the
            // DISCARDED t^{n+1}, and meanE = n_eps/n_e mixed two instants.
            // Measured: SIGFPE at step 6 after ONE discard. It survived
            // review because all four retryStep validation beds are LFA and
            // never exercised an electron-energy equation.
            //
            // A new energy model that forgets to override discardStep() fails
            // SILENTLY -- the base is a no-op. Test against a case that
            // actually retries (co1.5-forceretry forces it with
            // maxCorrectors 3).
            // ==============================================================
            transport.discardStep();
            if (energy) energy->discardStep();

            // The relaxation coordinator must forget the abandoned attempt's
            // iterate history too, or the retry starts from a residual that
            // describes a step that was thrown away.
            //
            // It also gets the chance to ESCALATE. Under `outerScheme
            // adaptive` the first retry of a step switches to Anderson and
            // re-runs at the SAME deltaT: try a stronger solver before
            // conceding the step is too long. A second failure returns
            // escalated == false and the retry ladder shortens the step as
            // usual, so the mechanism that always works is never bypassed.
            bool escalated = false;
            {
                plasmaOuterRelaxation* r =
                    plasmaOuterRelaxation::lookup(gasMesh());
                if (r)
                {
                    r->discardStep();
                    escalated = r->escalateForRetry();
                }
            }
            timeControl.prepareRetry(escalated);
            continue;                  // solve this step again
        }

        // Accepting a step the outer loop never converged: either the retries
        // are spent or deltaT is on its floor, so the solver has no lever
        // left. Keep the best iterate, but say so -- and stop the run if this
        // stops being an isolated event.
        if (timeControl.outerHitCap())
        {
            timeControl.noteDegradedStep();
        }

        stepDone = true;
        }   // end retry loop

        diagnostics.report();

        // AFTER the retry loop, so a discarded step never contributes: the
        // current is a diagnostic of the state that was actually accepted.
        dischargeCurrent.update(transport, species, em());

        // TEMPORAL ERROR MEASUREMENT (Phase 1: report only).
        //
        // AFTER the retry loop, for the same reason dischargeCurrent.update()
        // is: the estimate must describe the state that was ACCEPTED, not a
        // discarded attempt.
        //
        // The field list is assembled HERE because the solver is what knows
        // which models are active -- the same place that calls both
        // discardStep()s. Each model appends what IT transports, so a case is
        // covered whatever the user selected: species always, nEps_e under
        // LMEA, T_gas under gas heating.
        if (timeControl.reportTemporalError())
        {
            DynamicList<const volScalarField*> errFields;
            DynamicList<word> errNames;

            transport.appendTransportedFields(errFields, errNames);
            if (energy) energy->appendTransportedFields(errFields, errNames);

            timeControl.measureTemporalError(errFields, errNames);
        }


        runTime.write();
        runTime.printExecutionTime(Info);
    }

    plasmaSimulationProfiler::report();
    Info<< "End\n" << endl;

    return 0;
}
