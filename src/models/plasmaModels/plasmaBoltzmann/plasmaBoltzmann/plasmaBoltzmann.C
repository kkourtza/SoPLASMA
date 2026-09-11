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
#include "OFstream.H"
#include "OSspecific.H"
#include "Switch.H"
#include "error.H"

// The Boltzmann solver. Only this header, so the OpenFOAM build needs no Eigen
// include path -- the mechanism is constructed inside the library.
#include "MechTables.H"

#include <fstream>
#include <sstream>
#include "pressureUnits.H"
#include "IOdictionary.H"
#include "objectRegistry.H"
#include "Time.H"
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
    const Foam::fileName& tableDir,
    const Foam::scalar pGasPa
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

    // DERIVED from the gas the case actually runs at -- never stated here.
    //
    // Density-scaled (three-body) processes break E/N similarity: their
    // contribution to the EEDF scales with N, so the tables are valid only near
    // the density they were solved at. There is therefore exactly ONE correct
    // density for the sweep, the case's own, and nothing to choose.
    //
    // The library still wants atm because that is the convention of the
    // cross-section literature. That conversion is internal; a user only ever
    // types Pa.
    o.pressure_atm = Foam::constant::plasma::atmFromPa(pGasPa);

    // REJECTED, not silently migrated.
    //
    // `pressureAtm` was a second, independent spelling of the gas pressure, in
    // a different unit, in a different file from `backgroundGas/pressure`. It
    // is now derived. Reading it would be worse than ignoring it and ignoring
    // it would be worse than stopping, because -- as sweepStamp below says --
    // changing the sweep pressure changes EVERY table while leaving the
    // mechanism hash untouched. A case whose stated value disagreed with its
    // own gas would silently get a table set built for a different gas.
    if (b.found("pressureAtm"))
    {
        FatalIOErrorInFunction(b)
            << "`pressureAtm` in the `boltzmann` block is no longer read." << Foam::nl
            << Foam::nl
            << "    The sweep pressure is DERIVED from `backgroundGas/pressure`"
            << " in plasmaSpeciesProperties, which for this case is" << Foam::nl
            << "        "
            << Foam::constant::plasma::pressureStr(pGasPa).c_str() << Foam::nl
            << Foam::nl
            << "    Delete the entry. Pressure is stated ONCE, in Pa, in"
            << " `backgroundGas/pressure`." << Foam::nl
            << "    Three-body processes are included whenever the mechanism"
            << " has them; there is nothing to switch." << Foam::nl
            << exit(Foam::FatalIOError);
    }

    return o;
}


//- A stamp of the sweep settings, written beside the tables.
//
//  The mechanism hash alone is NOT enough to decide whether a table set is
//  current: changing ENmax, Tgas, pressureAtm or the growth model changes every
//  table while leaving the mechanism -- and therefore its hash -- untouched. A
//  case that raised ENmax from 500 to 2000 Td would have silently kept the old
//  500 Td tables and extrapolated above them.
Foam::string sweepStamp(const Boltzmann::MechTableOptions& o)
{
    std::ostringstream ss;
    ss << "EN=" << o.EN_min << ":" << o.EN_max << ":" << o.nPoints
       << " Tgas=" << o.T_gas << " Texc=" << o.T_exc
       << " grid=" << o.gridPoints
       << " pPa=" << o.pressure_atm*Foam::constant::plasma::PaPerAtm
       << " growth=" << o.growth << " floor=" << o.thermalFloor;
    return Foam::string(ss.str());
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

Foam::scalar Foam::plasmaBoltzmann::gasPressurePa
(
    const Foam::objectRegistry& obr
)
{
    IOdictionary sp
    (
        IOobject
        (
            "plasmaSpeciesProperties",
            obr.time().constant(),
            obr,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    // 101325 Pa if the case says nothing, matching plasmaSpecies, which closes
    // the number density from the same key with the same default. The two must
    // agree or the tables are solved for a different gas than the transport.
    return sp.subOrEmptyDict("backgroundGas")
             .getOrDefault<scalar>("pressure", constant::plasma::PaPerAtm);
}


Foam::plasmaBoltzmann::status Foam::plasmaBoltzmann::ensureTables
(
    const dictionary& chem,
    const fileName& manifest,
    const fileName& tableDir,
    const word& expectedHash,
    const scalar pGasPa
)
{
    if (!chem.getOrDefault<Switch>("generateTables", true))
    {
        // `generateTables no` MUST STILL CHECK THE STAMP.
        //
        // This branch used to return here, before the stamp was ever read, so
        // the one path that cannot self-correct was also the one with no check:
        // take a pre-built table set from another case, run at a different
        // pressure or a different ENmax, and every rate coefficient and
        // mobility comes from the wrong sweep -- converging, plausible, wrong.
        // "as found" was doing a great deal of quiet work.
        //
        // `no` is an explicit instruction NOT to rebuild, so the only honest
        // options are to stop or to lie. It stops.
        const fileName stampFile = tableDir/"sweep.stamp";
        const string want = sweepStamp(readOptions(chem, manifest, tableDir, pGasPa));

        if (!isFile(stampFile))
        {
            FatalErrorInFunction
                << "`generateTables no`, but " << tableDir
                << " carries no sweep.stamp." << nl
                << nl
                << "    Without it there is no way to tell what conditions"
                << " those tables were solved at, and a table set" << nl
                << "    built for a different gas or a different E/N range is"
                << " indistinguishable from a correct one." << nl
                << "    A missing stamp is therefore treated as a mismatch, not"
                << " as permission." << nl
                << nl
                << "    This case needs:  " << want.c_str() << nl
                << nl
                << "    Set `generateTables yes` to solve the sweep here, or"
                << " point `tableDir` at a set built for this case." << nl
                << exit(FatalError);
        }

        const string have = firstLineOf(stampFile);

        if (have != want)
        {
            // Name the field that moved. The stamp has eight, and "they differ"
            // would leave the reader diffing two long strings by eye.
            string diffs;
            {
                std::istringstream hs(have), ws(want);
                std::string ht, wt;
                while (hs >> ht && ws >> wt)
                {
                    if (ht != wt)
                    {
                        diffs += "        " + ht + "   ->   " + wt + "\n";
                    }
                }
            }

            FatalErrorInFunction
                << "`generateTables no`, but the tables in " << tableDir
                << " were not built for this case." << nl
                << nl
                << "    tables were built at:  " << have.c_str() << nl
                << "    this case needs:       " << want.c_str() << nl
                << nl
                << "    differing:" << nl
                << diffs.c_str()
                << nl
                << "    Every rate coefficient and mobility would come from the"
                << " wrong sweep. Set `generateTables yes` to" << nl
                << "    re-solve it here, or point `tableDir` at a set built for"
                << " this case." << nl
                << exit(FatalError);
        }

        Info<< "plasmaBoltzmann: generateTables off; the tables in " << tableDir
            << " match this case's sweep settings." << endl;
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
    const string stamp = sweepStamp(readOptions(chem, manifest, tableDir, pGasPa));
    const fileName stampFile = tableDir/"sweep.stamp";

    const fileName probe = tableDir/"muN_vs_reducedE";
    if (isFile(probe))
    {
        const string head = firstLineOf(probe);
        const bool stampOk =
            isFile(stampFile) && firstLineOf(stampFile) == stamp;

        if (!stampOk && isFile(probe))
        {
            Info<< "plasmaBoltzmann: sweep settings differ from the tables in "
                << tableDir << "; rebuilding" << nl
                << "    was: " << firstLineOf(stampFile).c_str() << nl
                << "    now: " << stamp.c_str() << endl;
        }

        if (stampOk && !expectedHash.empty()
         && head.find(expectedHash) != std::string::npos)
        {
            Info<< "plasmaBoltzmann: tables in " << tableDir
                << " match mechanism [" << expectedHash.c_str()
                << "]; reusing" << endl;
            return reused;
        }

        if (stampOk)
        {
            Info<< "plasmaBoltzmann: tables in " << tableDir
                << " do not carry mechanism [" << expectedHash.c_str()
                << "]; rebuilding" << endl;
        }
    }

    rebuild(chem, manifest, tableDir, pGasPa);
    return generated;
}


void Foam::plasmaBoltzmann::rebuild
(
    const dictionary& chem,
    const fileName& manifest,
    const fileName& tableDir,
    const scalar pGasPa,
    const HashTable<scalar>& composition,
    const scalar Tgas
)
{
    mkDir(tableDir);

    Boltzmann::MechTableOptions o = readOptions(chem, manifest, tableDir, pGasPa);

    if (Tgas > 0)
    {
        o.T_gas = Tgas;
        // T_exc follows T_gas unless the case pinned it. Leaving a stale
        // excitation temperature behind while the gas temperature moves would
        // silently change the superelastic balance.
        if (!chem.subOrEmptyDict("boltzmann").found("Texc"))
        {
            o.T_exc = -1;
        }
    }

    forAllConstIters(composition, it)
    {
        o.composition[it.key()] = it.val();
    }

    if (!composition.empty())
    {
        Info<< "plasmaBoltzmann: composition";
        forAllConstIters(composition, it)
        {
            Info<< " " << it.key() << "=" << it.val();
        }
        Info<< endl;
    }

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

    // Written LAST, so an interrupted sweep leaves no stamp and the next run
    // rebuilds rather than trusting a half-written set.
    {
        OFstream os(tableDir/"sweep.stamp");
        os << sweepStamp(o).c_str() << endl;
    }

    Info<< "plasmaBoltzmann: wrote tables to " << tableDir << endl;
}


// ************************************************************************* //
