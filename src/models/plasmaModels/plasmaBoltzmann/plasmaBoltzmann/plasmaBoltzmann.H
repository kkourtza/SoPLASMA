/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::plasmaBoltzmann

Description
    Runs the two-term electron Boltzmann solver in process and writes the rate
    and transport tables the run will interpolate.

    This removes the offline step. A case used to require someone to run
    `genMechTables` by hand with the right arguments and copy the result in,
    which meant the tables and the case could disagree about gas temperature,
    pressure or growth model with nothing to detect it. Here the sweep reads the
    same dictionary the solver does.

    Running it at start-up is affordable rather than merely tolerable: the sweep
    is OpenMP-parallel and takes ~0.45 s for 200 E/N points over 49 processes,
    against minutes for the mesh alone. It is also what makes an EEDF refresh
    during a run realistic (see `eedfRefresh`).

    Tables are always written to FILES, never held only in memory, for two
    reasons. They are the record of what the run actually used, and they are
    what you plot when a rate looks wrong. The format is unchanged, so anything
    that read them before still does.

    Staleness is decided by the mechanism hash, not by timestamps. A table set
    whose header carries the manifest's hash is reused; anything else is
    rebuilt. Timestamps would rebuild after a harmless touch and, worse, would
    NOT rebuild if a manifest were restored from an older copy.

SourceFiles
    plasmaBoltzmann.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaBoltzmann_H
#define plasmaBoltzmann_H

#include "objectRegistry.H"
#include "dictionary.H"
#include "fileName.H"
#include "word.H"
#include "HashTable.H"
#include "scalar.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class plasmaBoltzmann Declaration
\*---------------------------------------------------------------------------*/

class plasmaBoltzmann
{
public:

    //- Outcome of an ensureTables() call, for reporting
    enum status { generated, reused, disabled };

    //- Ensure `tableDir` holds a table set built from `manifest`.
    //
    //  Returns `reused` when an existing set already carries the manifest's
    //  mechanism hash, `generated` when the sweep ran, `disabled` when the
    //  dictionary asks for pre-built tables (`generateTables no`).
    //
    //  `dict` is the case's `chemistry` sub-dictionary. Sweep settings are read
    //  from its `boltzmann` sub-dictionary; the gas temperature and pressure
    //  come from `Tgas` and `pressure` there, defaulting to the values a case
    //  would otherwise state in `backgroundGas`.
    //- The gas pressure [Pa], from `backgroundGas/pressure` in
    //  constant/<region>/plasmaSpeciesProperties.
    //
    //  ONE OWNER, ONE UNIT. The sweep needs the gas density because three-body
    //  processes break E/N similarity -- their contribution to the EEDF scales
    //  with N, so the tables are valid only near the density they were solved
    //  at. That density is not an independent setting: it is the one the case
    //  runs at.
    //
    //  It used to be stated separately as `boltzmann/pressureAtm`, in ATM, so a
    //  case wrote the same physical quantity twice in two units and they could
    //  disagree. Every shipped case set it to the hand-computed 0.98692 against
    //  a pressure of 1e5 Pa; nothing used the `0` sentinel it also carried.
    static scalar gasPressurePa(const objectRegistry& obr);

    static status ensureTables
    (
        const dictionary& dict,
        const fileName& manifest,
        const fileName& tableDir,
        const word& expectedHash,
        const scalar pGasPa
    );

    //- Force a rebuild regardless of the existing table set. Used by the EEDF
    //  refresh path, where the mechanism has not changed but the conditions it
    //  was solved at have.
    //
    //  `composition` overrides the mechanism's reference mole fractions, and
    //  `Tgas` (when positive) the dictionary's gas temperature. Both default to
    //  "unchanged", which is what a start-up build wants. Passing the CURRENT
    //  values is the entire point of a refresh: re-solving at the composition
    //  the tables were already built at reproduces them exactly, at full cost.
    static void rebuild
    (
        const dictionary& dict,
        const fileName& manifest,
        const fileName& tableDir,
        const scalar pGasPa,
        const HashTable<scalar>& composition = HashTable<scalar>(),
        const scalar Tgas = -1
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
