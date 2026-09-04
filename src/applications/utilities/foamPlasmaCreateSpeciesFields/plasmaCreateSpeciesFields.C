/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025 Rention Pasolari

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
    plasmaCreateSpeciesFields

Description
    Utility to generate template number-density fields (n_speciesName) in the
    0/ directory for all species listed in the plasmaSpecies dictionary.

    The tool:
      - detects and uses the gas region (if present in a multi-region case),
      - reads species from the plasmaSpecies dictionary,
      - loads the mesh of the region,
      - creates volScalarField files with default values and zeroGradient
        boundary conditions for all patches.

    These generated fields are **only templates**.
    The user must review and modify the internalField values, boundary
    conditions, and any other field settings according to the needs of the
    specific plasma simulation.

Usage
    \b plasmaDielectricFoam [OPTIONS]

    Example:
        plasmaCreateSpeciesFields -case plasmaCase

        or 

        mpirun -np 4 plasmaCreateFields -parallel

Author
    Rention Pasolari
    Contact: r.pasolari@gmail.com
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "regionProperties.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "plasmaSpecies.H"
#include "boundaryRoleLibrary.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "wallPolyPatch.H"
#include "OFstream.H"
#include <sstream>
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- What a patch IS, as far as species transport is concerned.
//
//  DERIVED, never asked for. The potential's own boundary conditions and the
//  region declarations already say what every patch is; restating it per
//  species in changeDictionary is how a case comes to have a wall flux on n_e
//  and zeroGradient on nEps_e -- measured on needleDBD 2026-09-04, from a
//  hand-written 0.orig/nEps_e template that still named the STREAMER case's
//  patches (axis, wedge_0, grounded_electrode).
enum patchRole
{
    prMechanical,     //!< empty / wedge / symmetry / cyclic / processor
    prOpen,           //!< an open boundary: nothing is absorbed
    prSurfaceCharging,//!< a dielectric surface: absorbs AND accumulates sigma
    prConductor       //!< metal: absorbs, but the charge leaves via the circuit
};


//- Classify every patch of a region from the CASE'S OWN DESCRIPTION.
//
//  THE DIRECTION MATTERS. This used to read 0/<region>/ePotential and map the
//  ELECTROSTATICS boundary-condition type back to a semantic class --
//  `thinDielectricPotential` meant a charging surface, `uniformFixedValue`
//  meant a conductor. That is deriving one module's answer from another
//  module's answer, and it is backwards under G2: both modules must derive from
//  the DESCRIPTION. It also forced an ordering (the potential had to be
//  generated first) that no longer exists.
//
//  Three things are still derived rather than declared, and a case may not
//  state them: MECHANICAL patches come from the mesh, region INTERFACES from
//  the topology, and everything else from configuration/boundaries via the
//  shared role library.
static void classifyPatches
(
    const fvMesh& mesh,
    const Time& runTime,
    HashTable<patchRole>& roles
)
{
    wordList farFields, dielectrics, gases;
    {
        IOdictionary rp
        (
            IOobject
            (
                "regionProperties",
                runTime.constant(),
                runTime,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        if (rp.found("regions"))
        {
            HashTable<wordList> regions;
            rp.readEntry("regions", regions);

            if (regions.found("farField"))   farFields   = regions["farField"];
            if (regions.found("dielectric")) dielectrics = regions["dielectric"];
            if (regions.found("gas"))        gases       = regions["gas"];
        }
    }

    const dictionary decl(boundaryRoleLibrary::caseDeclaration(runTime));

    forAll(mesh.boundary(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        const word& name = pp.name();

        // 1. MECHANICAL -- from the MESH. A mechanical constraint is not a
        //    physical description, and a flux condition on one is meaningless.
        if
        (
            isA<emptyPolyPatch>(pp)
         || isA<wedgePolyPatch>(pp)
         || isA<symmetryPolyPatch>(pp)
         || isA<symmetryPlanePolyPatch>(pp)
         || isA<cyclicPolyPatch>(pp)
         || isA<processorPolyPatch>(pp)
        )
        {
            roles.insert(name, prMechanical);
            continue;
        }

        // 2. A REGION INTERFACE -- from the TOPOLOGY.
        bool isInterface = false;

        for (const word& ff : farFields)
        {
            if (name == mesh.name() + "_to_" + ff)
            {
                // Fictitious AIR, not a surface: nothing is absorbed and
                // nothing may accumulate. Charge deposited here would invent a
                // dielectric surface in the middle of the gas and shield the
                // very field the region was added to resolve.
                roles.insert(name, prOpen);
                isInterface = true;
            }
        }
        if (isInterface) continue;

        for (const word& d : dielectrics)
        {
            if (name == mesh.name() + "_to_" + d)
            {
                roles.insert(name, prSurfaceCharging);
                isInterface = true;
            }
        }
        for (const word& g : gases)
        {
            if (name == mesh.name() + "_to_" + g)
            {
                roles.insert(name, prSurfaceCharging);
                isInterface = true;
            }
        }
        if (isInterface) continue;

        // 3. DECLARED by the case, and translated by the role library.
        if (!decl.found(name))
        {
            FatalErrorInFunction
                << "Patch `" << name << "` of region `" << mesh.name()
                << "` is not described." << nl << nl
                << "    Add a block to configuration/boundaries saying what it"
                   " IS. A patch is NOT" << nl
                << "    declared only when it is a mechanical constraint"
                   " (empty, wedge, symmetry," << nl
                << "    processor) or a region interface -- both are derived."
                   " This one is neither." << nl << nl
                << "    Available kinds: " << boundaryRoleLibrary::kinds()
                << nl
                << "    `plasmaSetupBoundaries -listKinds` describes them."
                << nl
                << exit(FatalError);
        }

        const word kind(decl.subDict(name).get<word>("kind"));

        const word surf
        (
            boundaryRoleLibrary::surfaceClass
            (
                kind,
                "configuration/boundaries, patch `" + name + "`"
            )
        );

        roles.insert
        (
            name,
            surf == "chargingSurface" ? prSurfaceCharging
          : surf == "conductor"       ? prConductor
          : surf == "mechanical"      ? prMechanical
          :                             prOpen
        );
    }
}


//- The boundary-condition entry for one patch of one field.
static void writePatchEntry
(
    Ostream& os,
    const fvMesh& mesh,
    const label patchi,
    const patchRole role,
    const word& kind,          //!< "electron" | "ion" | "energy" | "neutral"
    const word& fluxFamily,    //!< "Mixed" | "Implicit"
    const word& material,      //!< the surface's material, or word::null
    const dictionary& emission //!< emission mechanisms, possibly empty
)
{
    const polyPatch& pp = mesh.boundaryMesh()[patchi];

    os  << "    " << pp.name() << nl << "    {" << nl;

    if (role == prMechanical)
    {
        // The mechanical constraint IS the boundary condition here.
        word t("zeroGradient");
        if (isA<emptyPolyPatch>(pp))              t = "empty";
        else if (isA<wedgePolyPatch>(pp))         t = "wedge";
        else if (isA<symmetryPolyPatch>(pp))      t = "symmetry";
        else if (isA<symmetryPlanePolyPatch>(pp)) t = "symmetryPlane";
        else if (isA<cyclicPolyPatch>(pp))        t = "cyclic";
        else if (isA<processorPolyPatch>(pp))     t = "processor";

        os << "        type            " << t << ";" << nl;
    }
    else if (role == prOpen || kind == "neutral")
    {
        // Neutrals are not absorbed by a wall flux condition: the drift-
        // diffusion wall flux is built on the charged-particle thermal and
        // drift fluxes. An open boundary absorbs nothing either.
        os << "        type            zeroGradient;" << nl;
    }
    else
    {
        // A SOLID SURFACE. Every charged species gets a drift-diffusion wall
        // flux, and so does the electron energy density -- thermal-only is
        // SINGULAR at a sharp electrode, and setting the drift flux on one of
        // the pair leaves the other singular, so a partial fix looks like a
        // failure.
        const word type
        (
            kind == "electron" ? "electronDDWallFlux" + fluxFamily
          : kind == "energy"   ? "energyDDWallFlux"   + fluxFamily
          :                      "ionDDWallFlux"      + fluxFamily
        );

        os  << "        type            " << type << ";" << nl
            << "        value           uniform 0;" << nl;

        // `T` is a FIELD NAME, not a temperature: the field the wall flux
        // evaluates the thermal speed from. Electrons and their energy are at
        // the ELECTRON temperature; ions are at the gas temperature.
        os  << "        T               "
            << (kind == "ion" ? "T_gas" : "T_e") << ";" << nl;

        os  << "        includeDriftFlux true;" << nl;

        // SIGMA IS FOR CHARGED SPECIES ONLY.
        //
        // nEps_e is an ENERGY density, not a charge density. Counting it would
        // add q_e times an energy flux to a surface charge. It cannot happen
        // by accident today -- the accumulation loops the SPECIES list and
        // nEps_e is not a species, it belongs to plasmaEnergy -- so the entry
        // would be inert rather than wrong; it is omitted so that nobody reads
        // it and concludes otherwise.
        if (kind != "energy")
        {
            os  << "        enableSurfaceCharging "
                << (role == prSurfaceCharging ? "true" : "false") << ";" << nl;
        }

        // THE PATCH'S MATERIAL, passed through so an emission model can
        // resolve a work function, a secondary yield or a band structure from
        // the shipped libraries instead of being handed a number. This is the
        // route that lets a patch say what it is MADE OF.
        //
        // Written on the ELECTRON condition only: emission is electron
        // emission, and an ion or energy condition has nothing to do with it.
        if (kind == "electron" && !material.empty())
        {
            os  << "        material        " << material << ";" << nl;
        }

        // EMISSION MECHANISMS, verbatim from the boundary role (or a per-patch
        // override). Empty for every shipped kind -- see the note in
        // etc/boundaryRoles for why switching one on by default would be the
        // wrong default.
        if (kind == "electron" && !emission.empty())
        {
            os  << "        emission" << nl << "        {" << nl;

            OStringStream eos;
            emission.write(eos, false);

            // Re-indent the block so the generated file stays readable.
            const string body(eos.str());
            std::string line;
            std::istringstream is(body);
            while (std::getline(is, line))
            {
                if (line.find_first_not_of(" \t") == std::string::npos) continue;
                os << "        " << line.c_str() << nl;
            }

            os  << "        }" << nl;
        }
    }

    os << "    }" << nl << nl;
}


//- Write one scalar field file with DERIVED boundary conditions.
//
//  Written TEXTUALLY rather than by constructing a volScalarField with real
//  boundary conditions: the wall-flux conditions look up `plasmaTransport` in
//  the object registry, which does not exist in a generation utility.
static void writeDerivedField
(
    const fvMesh& mesh,
    const Time& runTime,
    const word& fieldName,
    const dimensionSet& dims,
    const scalar internalValue,
    const word& kind,
    const word& fluxFamily,
    const HashTable<patchRole>& roles,
    const dictionary& decl
)
{
    IOobject io
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    mkDir(io.path());

    OFstream os(io.objectPath());

    io.writeHeader(os, volScalarField::typeName);

    os  << "dimensions      " << dims << ";" << nl << nl
        << "internalField   uniform " << internalValue << ";" << nl << nl
        << "boundaryField" << nl << "{" << nl;

    forAll(mesh.boundary(), patchi)
    {
        const word& name = mesh.boundaryMesh()[patchi].name();

        // The patch's material and emission block, from its DECLARATION and
        // its KIND's role. A per-patch `emission` in configuration/boundaries
        // overrides the kind's, so a single experiment does not require
        // editing the shipped library.
        word material(word::null);
        dictionary emission;

        if (decl.found(name))
        {
            const dictionary& pd = decl.subDict(name);

            material = pd.getOrDefault<word>("material", word::null);

            if (pd.found("kind"))
            {
                const dictionary& role = boundaryRoleLibrary::lookup
                (
                    pd.get<word>("kind"),
                    "configuration/boundaries, patch `" + name + "`"
                );

                if (role.found("emission")) emission = role.subDict("emission");
            }

            if (pd.found("emission")) emission = pd.subDict("emission");
        }

        writePatchEntry
        (
            os, mesh, patchi,
            roles.found(name) ? roles[name] : prOpen,
            kind, fluxFamily, material, emission
        );
    }

    os << "}" << nl << nl;

    io.writeEndDivider(os);
}

} // End namespace Foam



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Create species fields in 0 folder"
    );

    argList::noBanner();
    argList::noJobInfo();
    argList::noFunctionObjects();

    argList::addOption
    (
        "dict",
        "file",
        "Alternative plasmaSpecies dictionary"
    );

    argList args(argc, argv);

    // Create minimal Time
    Time runTime
    (
        Time::controlDictName,
        args.rootPath(),
        args.caseName()
    );

    word gasRegionName;

    fileName globalConstantDir = runTime.globalPath()/"constant";
    fileName regionPropsFile = runTime.constant()/"regionProperties";

    // Find the gas region, or default region in single region cases
    if (isFile(regionPropsFile))
    {
        if (Pstream::master())
        {
            Info << "Found regionProperties: " << regionPropsFile << nl;
        }

        regionProperties rp(runTime);

        if (rp.found("gas"))
        {
            const wordList& names = rp["gas"];

            if (names.size() == 0)
            {
                FatalErrorInFunction
                    << "A 'gas' region type is defined in regionProperties, "
                    << "but no region names are listed under it." << nl
                    << exit(FatalError);
            }

            if (names.size() > 1)
            {
                FatalErrorInFunction
                    << "Multiple gas regions detected in "
                    << "constant/regionProperties:" << nl
                    << "  gas: " << names << nl
                    << "This utility supports ONLY ONE gas region." << nl
                    << exit(FatalError);
            }

            // Exactly one gas region
            gasRegionName = names[0];

            if (Pstream::master())
            {
                Info << "Detected gas region: " << gasRegionName << nl;
            }
        }
    }

    // Report region used
    if (Pstream::master())
    {
        if (gasRegionName.empty())
            Info << "No gas region found → using constant/ as path." << nl;
        else
            Info << "Using region: " << gasRegionName << nl;
    }

    fileName plasmaSpeciesPath;

    // Find the plasmaSpecies file in the case
    if (!gasRegionName.empty())
    {
        // Multi-region path
        fileName regionDir = globalConstantDir/gasRegionName;

        if (isFile(regionDir/"plasmaSpeciesProperties"))
        {
            plasmaSpeciesPath = regionDir/"plasmaSpeciesProperties";
        }
        else
        {
            FatalErrorInFunction
                << "Region '" << gasRegionName << "' exists but has no "
                << "plasmaSpeciesProperties file:" << nl
                << "  " << regionDir/"plasmaSpeciesProperties" << nl
                << exit(FatalError);
        }
    }
    else
    {
        // Single-region path
        if (isFile(runTime.constant()/"plasmaSpeciesProperties"))
        {
            plasmaSpeciesPath = runTime.constant()/"plasmaSpeciesProperties";
        }
        else
        {
            FatalErrorInFunction
                << "No 'constant/plasmaSpeciesProperties' file found." << nl 
                << exit(FatalError);
        }
    }

    // Read plasmaSpecies dict
    if (Pstream::master())
    {
        Info << "Reading plasmaSpecies from: " << plasmaSpeciesPath << nl;
    }

    IOdictionary speciesDict
    (
        IOobject
        (
            "plasmaSpeciesProperties",
            plasmaSpeciesPath.path(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (!speciesDict.found("activeSpecies"))
    {
        FatalErrorInFunction
            << "The plasmaSpecies dictionary has no 'activeSpecies' entry!"
            << exit(FatalError);
    }

    // `activeSpecies fromMechanism;` is resolved through plasmaSpecies, which
    // is also what the solver uses. This utility must create EXACTLY the fields
    // the solver will later look for, so the two cannot be allowed to derive
    // the list independently.
    wordList species;

    // CHARGE DECIDES THE BOUNDARY CONDITION, so it is read here rather than
    // guessed from the name: an electron gets a condition with secondary
    // emission, every other charged species gets one without, and a neutral
    // gets none. `n_Om` and `n_O2m` are ions despite being negative.
    HashTable<scalar> charges;
    {
        ITstream& is = speciesDict.lookup("activeSpecies");
        token firstToken(is);
        is.rewind();

        if (firstToken.isWord() && firstToken.wordToken() == "fromMechanism")
        {
            species = plasmaSpecies::speciesFromMechanism
        (
            speciesDict, &charges
        );
        }
        else
        {
            is >> species;
        }
    }

    if (Pstream::master())
    {
        Info << "Found " << species.size() << " active species." << nl;
    }

    for (const word& s : species)
    {
        Info << "  - " << s << nl;
    }

    // Determine region path for mesh creation
    word regionToLoad = gasRegionName;

    // Create the mesh
    fvMesh mesh
    (
        IOobject
        (
            regionToLoad,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Field Creation Loop
    dimensionSet dimDensity(0, -3, 0, 0, 0, 0, 0);

    // WHAT EVERY PATCH IS, derived once from the potential's own boundary
    // conditions and the region declarations. See classifyPatches().
    HashTable<patchRole> roles;
    classifyPatches(mesh, runTime, roles);

    const dictionary decl(boundaryRoleLibrary::caseDeclaration(runTime));

    // The wall-flux family. `Mixed` imposes the flux through the mixed
    // condition's valueFraction; `Implicit` writes it into the matrix. Mixed is
    // the default because it is the one the validated multi-region cases use.
    word fluxFamily("Mixed");
    {
        IOdictionary controls
        (
            IOobject
            (
                "plasmaSimulationControls",
                runTime.system(),
                runTime,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        if (controls.found("wallFluxFamily"))
        {
            fluxFamily = controls.get<word>("wallFluxFamily");

            if (fluxFamily != "Mixed" && fluxFamily != "Implicit")
            {
                FatalIOErrorInFunction(controls)
                    << "Unknown `wallFluxFamily` `" << fluxFamily
                    << "`. Valid: Mixed | Implicit." << nl
                    << exit(FatalIOError);
            }
        }
    }

    if (Pstream::master())
    {
        Info<< nl << "Boundary conditions DERIVED per patch (wall-flux family "
            << fluxFamily << "):" << nl;

        forAll(mesh.boundary(), patchi)
        {
            const word& name = mesh.boundaryMesh()[patchi].name();
            const patchRole r = roles.found(name) ? roles[name] : prOpen;

            Info<< "    " << name << "  ->  "
                << (
                       r == prMechanical      ? "mechanical (kept as-is)"
                     : r == prOpen            ? "open: zeroGradient, no charging"
                     : r == prSurfaceCharging ? "dielectric surface: wall flux + surface charging"
                     :                          "conductor: wall flux, NO local sigma"
                   ) << nl;
        }
        Info<< endl;
    }

    for (const word& s : species)
    {
        const word fieldName = "n_" + s;

        const scalar q = charges.found(s) ? charges[s] : 0.0;

        // The electron is the negatively charged species the mechanism calls
        // `e` (or `E`). Everything else charged is an ion, and gets a
        // condition WITHOUT secondary emission -- an ion striking a surface
        // may release an electron, but an ion arriving is not itself an
        // emission event.
        const word kind
        (
            mag(q) < SMALL   ? "neutral"
          : (s == "e" || s == "E" || s == "el") ? "electron"
          :                  "ion"
        );

        if (Pstream::master())
        {
            Info << "Creating field: " << fieldName
                 << "  (charge " << q << ", " << kind << ")" << endl;
        }

        writeDerivedField
        (
            mesh, runTime, fieldName, dimDensity, 0.0,
            kind, fluxFamily, roles, decl
        );
    }

    // THE ELECTRON ENERGY DENSITY IS GENERATED TOO, under LMEA.
    //
    // It used to be a HAND-WRITTEN 0.orig/nEps_e per case, and on needleDBD
    // that template was a copy of the STREAMER case's -- it named `axis`,
    // `wedge_0`, `wedge_1` and `grounded_electrode`, none of which exist in
    // that mesh, and had no entry for the barrier surface at all. The result
    // was a wall flux on n_e and zeroGradient on nEps_e, which is exactly the
    // pairing that leaves the energy equation singular at a sharp electrode.
    // Measured 2026-09-04.
    {
        const word energyModel
        (
            speciesDict.getOrDefault<word>("electronEnergyModel", "none")
        );

        if (energyModel == "LMEA")
        {
            const dimensionSet dimEnergyDensity(0, -3, 0, 0, 0, 0, 0);

            if (Pstream::master())
            {
                Info << "Creating field: nEps_e  (electron energy density,"
                     << " LMEA)" << endl;
            }

            writeDerivedField
            (
                mesh, runTime, "nEps_e", dimEnergyDensity, 0.0,
                "energy", fluxFamily, roles, decl
            );
        }
    }

    if (Pstream::master())
    {
        Info << nl << "Fields created successfully." << nl;
    }

    return 0;
}
