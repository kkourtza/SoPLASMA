/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
\*---------------------------------------------------------------------------*/

#include "plasmaBoltzmann.H"

#include "IFstream.H"
#include "OSspecific.H"
#include "Switch.H"
#include "error.H"

// The Boltzmann solver. Only this header, so the OpenFOAM build needs no Eigen
// include path -- the mechanism is constructed inside the library.
#include "MechTables.H"

#include <fstream>
#include <sstream>
#include <string>

// * * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * //

namespace
{

//- Sweep settings from the `boltzmann` sub-dictionary.
//
//  Every default here is the same one genMechTables uses, so a case that says
//  nothing gets exactly the table set the CLI would have produced. The two must
//  agree, which is why the sweep itself lives in the library and not in either.
Boltzmann::MechTableOptions readOptions
(
    const Foam::dictionary& chem,
    const Foam::fileName& manifest,
    const Foam::fileName& tableDir
)
{
    const Foam::dictionary b = chem.subOrEmptyDict("boltzmann");

    Boltzmann::MechTableOptions o;
    o.manifest = manifest;
    o.outDir   = tableDir;

    o.EN_min      = b.getOrDefault<Foam::scalar>("ENmin", 0.1);
    o.EN_max      = b.getOrDefault<Foam::scalar>("ENmax", 2000.0);
    o.nPoints     = b.getOrDefault<Foam::label>("nPoints", 200);
    o.T_gas       = b.getOrDefault<Foam::scalar>("Tgas", 300.0);
    o.T_exc       = b.getOrDefault<Foam::scalar>("Texc", -1.0);
    o.gridPoints  = b.getOrDefault<Foam::label>("gridPoints", 200);
    o.eedfPoints  = b.getOrDefault<Foam::label>("eedfPoints", 200);
    o.thermalFloor = b.getOrDefault<Foam::scalar>("thermalFloor", 0.01);
    o.writeEEDF   = b.getOrDefault<Foam::Switch>("writeEEDF", true);
    o.eedfNative  = b.getOrDefault<Foam::Switch>("writeEEDFNative", false);

    // Growth model. Defaults to temporal, and deliberately so: with no growth
    // model the mean energy runs away wherever ionisation is strong -- 91 eV at
    // 2000 Td in dry air, which looks like a solver failure and is not.
    o.growth = b.getOrDefault<Foam::word>("growthModel", "temporal");

    // Pressure in atm, for density-scaled (three-body) processes. 0 leaves them
    // out and keeps the tables density-independent; a positive value includes
    // them, and the resulting tables are valid only near that density because a
    // three-body process breaks E/N similarity.
    //
    // The default is 1 atm rather than 0: e + O2 + M -> O2- + M is the dominant
    // electron loss channel in atmospheric air below ~50 Td, and excluding it
    // by default would quietly bias the electron density high in exactly the
    // cases this solver is aimed at.
    o.pressure_atm = b.getOrDefault<Foam::scalar>("pressureAtm", 1.0);

    return o;
}


//- First line of a table file, where the mechanism hash is recorded.
Foam::string firstLineOf(const Foam::fileName& path)
{
    std::ifstream f(path);
    std::string line;
    if (f && std::getline(f, line)) return Foam::string(line);
    return Foam::string::null;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::plasmaBoltzmann::status Foam::plasmaBoltzmann::ensureTables
(
    const dictionary& chem,
    const fileName& manifest,
    const fileName& tableDir,
    const word& expectedHash
)
{
    if (!chem.getOrDefault<Switch>("generateTables", true))
    {
        Info<< "plasmaBoltzmann: generateTables off; using the tables in "
            << tableDir << " as found" << endl;
        return disabled;
    }

    if (!isFile(manifest))
    {
        FatalErrorInFunction
            << "Cannot find the mechanism manifest " << manifest << nl
            << "    It is produced by mechc alongside the .foam dictionary."
            << " Either point `manifest` at it, or set `generateTables no`"
            << " to use a pre-built table set." << nl
            << exit(FatalError);
    }

    // Reuse only if an existing table carries the manifest's hash. Checking one
    // representative table is enough because the whole set is written together
    // by one sweep -- a half-written set is not a state that occurs.
    const fileName probe = tableDir/"muN_vs_reducedE";
    if (isFile(probe))
    {
        const string head = firstLineOf(probe);
        if (!expectedHash.empty() && head.find(expectedHash) != std::string::npos)
        {
            Info<< "plasmaBoltzmann: tables in " << tableDir
                << " match mechanism [" << expectedHash.c_str()
                << "]; reusing" << endl;
            return reused;
        }

        Info<< "plasmaBoltzmann: tables in " << tableDir
            << " do not carry mechanism [" << expectedHash.c_str()
            << "]; rebuilding" << endl;
    }

    rebuild(chem, manifest, tableDir);
    return generated;
}


void Foam::plasmaBoltzmann::rebuild
(
    const dictionary& chem,
    const fileName& manifest,
    const fileName& tableDir
)
{
    mkDir(tableDir);

    const Boltzmann::MechTableOptions o = readOptions(chem, manifest, tableDir);

    Info<< "plasmaBoltzmann: solving the EEDF over " << o.nPoints
        << " points, E/N " << o.EN_min << " .. " << o.EN_max << " Td"
        << ", T_gas " << o.T_gas << " K"
        << ", growth " << o.growth.c_str() << endl;

    // The sweep writes its own progress. Sent to a stringstream and echoed at
    // the end so that it cannot interleave with OpenFOAM's output under MPI,
    // where every rank runs this.
    std::ostringstream log;
    int rc = 1;
    try
    {
        rc = Boltzmann::buildMechTablesFromManifest(o, log);
    }
    catch (const std::exception& e)
    {
        FatalErrorInFunction
            << "The Boltzmann sweep failed: " << e.what() << nl
            << log.str().c_str() << nl << exit(FatalError);
    }

    if (rc != 0)
    {
        // Non-zero means some E/N points did not converge. Fatal rather than a
        // warning: an unconverged row is a silently wrong rate coefficient at
        // that field, and the run would interpolate straight through it.
        FatalErrorInFunction
            << "The Boltzmann sweep reported unconverged points." << nl
            << log.str().c_str() << nl
            << "    Refine `gridPoints`, narrow the E/N range, or fix the"
            << " mechanism." << nl << exit(FatalError);
    }

    Info<< "plasmaBoltzmann: wrote tables to " << tableDir << endl;
}


// ************************************************************************* //
