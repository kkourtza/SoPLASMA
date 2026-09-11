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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::plasmaWallBC

Description
    Abstract wall boundary condition interface. Serves as a base for all
    SoPLASMA wall boundary conditions that need to support surface charging.
    All wall BCs derive from this so they can be found via dynamic_cast
    using this type.
    
\*---------------------------------------------------------------------------*/

#ifndef plasmaWallBC_H
#define plasmaWallBC_H

#include "HashSet.H"
#include "Ostream.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class plasmaWallBC Declaration
\*---------------------------------------------------------------------------*/

class plasmaWallBC
{
protected:

    //- THE KEYS THE CASE ACTUALLY STATED, as opposed to those the condition
    //  defaulted. Shared by BOTH wall-flux families, so the concept has one
    //  name and one definition (G3).
    //
    //  WHY THIS EXISTS. `soPlasmaFoam` calls `runTime.writeNow()` at start-up,
    //  which overwrites `0/` with the fully-RESOLVED boundary dictionaries. If
    //  `write()` emits a defaulted value, that default is baked into the
    //  case's own `0/`, and a dictionary value beats a code default -- so A
    //  CHANGED DEFAULT CAN NEVER AGAIN REACH ANY CASE THAT HAS BEEN RUN ONCE.
    //  Such a case keeps the old physics while looking identical to one that
    //  got the new.
    //
    //  MEASURED 2026-09-04: exactly this hid the eq. (6.14) -> (6.15) energy
    //  weight change. A needleDBD tree carried `fluxEnergyFactor
    //  1.33333333333` in its own `0/gas/nEps_e`, written by an earlier run, so
    //  three successive runs all used 4/3 and a 5/3-vs-4/3 comparison came out
    //  at a ratio of exactly 1.0000. The generator never emits the key and no
    //  tracked tutorial carries it, so the contamination came purely from
    //  having been run.
    //
    //  This is G1: "a setting stated in two places is a defect even when both
    //  copies agree, because they will not agree later." The default lived in
    //  the code AND in `0/`, and they disagreed. So an optional entry is
    //  written back ONLY if the case supplied it. Restart is unaffected: an
    //  absent key defaults exactly as it does for a fresh case.
    //
    //  THE TRADE-OFF, stated rather than buried: a written time directory is
    //  no longer a complete record of what was used. That is what fixes the
    //  bug; the record moves to the start-up log, which cannot go stale.
    wordHashSet suppliedKeys_;

public:

    //- Destructor
    virtual ~plasmaWallBC() = default;

    // Wall properties
    virtual bool enableSurfaceCharging() const
    {
        return false;
    }

    //- Did the case state this key?
    inline bool wasSupplied(const word& k) const
    {
        return suppliedKeys_.found(k);
    }

    //- Write an OPTIONAL entry only if the case stated it. See suppliedKeys_
    //  for why writing a defaulted value is the defect.
    template<class T>
    void writeEntryIfSupplied(Ostream& os, const word& k, const T& v) const
    {
        if (wasSupplied(k))
        {
            os.writeEntry(k, v);
        }
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
