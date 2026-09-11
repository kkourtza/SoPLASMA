/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Description
    ONE UNIT FOR PRESSURE INPUT, PASCAL, AND ONE WAY TO REPORT IT.

    Every user-facing pressure in a SoPLASMA case is stated in Pa. There is no
    second unit to convert, no second key, and nothing to keep consistent.

    That was not always so. Until 2026-09-02 one physical quantity was read from
    THREE keys in TWO files in TWO units:

        backgroundGas/pressure         [Pa]   plasmaSpecies      -> gas density
        backgroundGas/energy/pressure  [Pa]   plasmaEnergy,      -> each divided
                                              plasmaTransport       by 101325
        boltzmann/pressureAtm          [atm]  plasmaBoltzmann    -> the density
                                              the three-body tables are built at

    A case therefore wrote 1e5 twice and 0.98692 once, for one gas. Every
    shipped case that set `pressureAtm` set it to the same hand-computed
    0.98692, and none used the `0` sentinel it also carried -- so the key bought
    nothing and cost a number that could silently disagree with the other two.

    REPORTING IS THE MIRROR OF INPUT. Pa is authoritative and comes first, with
    the conversions in parentheses, so a value can be read against whichever
    unit the reader's literature uses without anyone having to convert by hand:

        1e+05 Pa (0.98692 atm, 750.06 Torr, 1 bar)

    Use `pressureStr()` for anything printed. Use `toPa()`/`atmFromPa()` where a
    conversion is genuinely needed internally -- a Boltzmann sweep's three-body
    channels are tabulated per atmosphere because that is the convention of the
    cross-section literature, not because a user should ever type it.

\*---------------------------------------------------------------------------*/

#ifndef pressureUnits_H
#define pressureUnits_H

#include "scalar.H"
#include "string.H"
#include <sstream>

namespace Foam
{
namespace constant
{
namespace plasma
{

//- Standard atmosphere [Pa]. The ONE definition; do not re-spell 101325.
constexpr scalar PaPerAtm  = 101325.0;

//- Torr [Pa] = atm/760
constexpr scalar PaPerTorr = PaPerAtm/760.0;

//- Bar [Pa]
constexpr scalar PaPerBar  = 1.0e5;

inline scalar atmFromPa(const scalar pPa)  { return pPa/PaPerAtm; }
inline scalar torrFromPa(const scalar pPa) { return pPa/PaPerTorr; }
inline scalar barFromPa(const scalar pPa)  { return pPa/PaPerBar; }

//- Pressure for HUMAN OUTPUT: Pa first, conversions in parentheses.
//
//  Every Info<< of a pressure goes through this, so the format cannot drift
//  between subsystems and a reader never has to convert to compare against a
//  paper quoted in Torr.
inline string pressureStr(const scalar pPa)
{
    std::ostringstream os;
    os.precision(6);
    os << pPa << " Pa ("
       << atmFromPa(pPa)  << " atm, "
       << torrFromPa(pPa) << " Torr, "
       << barFromPa(pPa)  << " bar)";
    return string(os.str());
}

} // End namespace plasma
} // End namespace constant
} // End namespace Foam

#endif

// ************************************************************************* //
