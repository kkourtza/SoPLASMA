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
                word::null
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
        const dictionary& bg = species.backgroundDict();
        if (bg.subOrEmptyDict("energy").getOrDefault<bool>("solve", false))
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
    timeControl.setInitialDeltaT(transport);

    //- Create the plasmaSimulationDiagnostics manager
    plasmaSimulationDiagnostics diagnostics(runTime, transport);

    #include "reportSimulationSummary.H"

    runTime.writeNow();

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.run())
    {
        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;
        gasMesh().update();

        timeControl.adjustDeltaT(transport);

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

        diagnostics.report();

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    plasmaSimulationProfiler::report();
    Info<< "End\n" << endl;

    return 0;
}
