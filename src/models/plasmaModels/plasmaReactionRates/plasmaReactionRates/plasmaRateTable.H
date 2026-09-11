/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::plasmaRateTable

Description:
    One-dimensional lookup for plasma rate coefficients, with a user-selectable
    response beyond the tabulated range.

    Exists because OpenFOAM's interpolationTable offers only clamp, warn or
    error outside its range -- there is no extrapolation mode -- and clamping a
    rate coefficient is not the neutral choice it appears to be. It replaces the
    curve with a flat plateau and puts a slope discontinuity at the table edge,
    which for an avalanche is a systematic bias rather than a conservative one.

    The mode governs the UPPER end only. Below the table the value is always
    held at the first entry, whatever the mode: the first row is the zero-field
    value, so there is nothing below it to extrapolate towards, and a downward
    power law would diverge as x -> 0 rather than vanish. Extrapolating at the
    low end would be actively wrong, not merely uncertain.

    Modes (dictionary entry `outOfBounds`), applying above the table:

      clamp        hold the end value. Predictable, and defensible when the
                   tabulated range genuinely covers the physics.

      extrapolate  continue the trend as a power law, i.e. linear in
                   log(y) against log(x), from the last two tabulated points.
                   Rate coefficients are smooth and monotone in that space over
                   any reasonable range, so this degrades gracefully instead of
                   lying flat. This is the default.

    Neither is "safe" above the table: the two-term approximation itself is
    unreliable much above ~1000 Td. The choice is between two kinds of
    uncertainty, and extrapolation is the better-behaved one -- which is why the
    class also counts how many lookups left the range, so the cost of the choice
    is visible rather than assumed.

    Interpolation INSIDE the range is linear in x, matching what
    interpolationTable does, so switching to this class does not silently change
    results where the table already covered the field.

SourceFiles:
    plasmaRateTable.C
\*---------------------------------------------------------------------------*/

#ifndef plasmaRateTable_H
#define plasmaRateTable_H

#include "scalarField.H"
#include "fileName.H"
#include "dictionary.H"
#include "word.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class plasmaRateTable Declaration
\*---------------------------------------------------------------------------*/

class plasmaRateTable
{
public:

    //- Response outside the tabulated range
    enum boundsHandling { bhClamp, bhExtrapolate };

    //- Read a mode name; FatalError on anything else
    static boundsHandling boundsFromWord(const word& w);


private:

    // Private Data

        //- Tabulated abscissae, strictly increasing
        scalarField x_;

        //- Tabulated values
        scalarField y_;

        boundsHandling bounds_ = bhExtrapolate;

        //- Source file, for diagnostics
        fileName path_;

        //- Precomputed log-log slope of the last interval, used above the
        //  table. Zero when the end values are not both positive, in which case
        //  the top is clamped regardless of mode -- a power law through a zero
        //  is not defined, and inventing one would be worse than holding.
        scalar topSlope_ = 0;
        bool   topIsPowerLaw_ = false;

        //- Lookups that fell outside, counted so the cost of the bounds choice
        //  can be reported rather than assumed
        mutable label nBelow_ = 0;
        mutable label nAbove_ = 0;


public:

    // Constructors

        plasmaRateTable() = default;

        //- Read from an OpenFOAM-format table file: ( (x y) (x y) ... )
        plasmaRateTable(const fileName& path, boundsHandling bounds);


    // Member Functions

        //- Value at x
        scalar value(const scalar x) const;

        //- Fill out with value(x) for every entry of in
        void value(const scalarField& in, scalarField& out) const;

        //- Tabulated range
        scalar xMin() const { return x_.size() ? x_.first() : 0; }
        scalar xMax() const { return x_.size() ? x_.last()  : 0; }

        label size() const { return x_.size(); }
        const fileName& path() const { return path_; }
        boundsHandling bounds() const { return bounds_; }

        //- Out-of-range lookup counts since the last reset
        label nBelow() const { return nBelow_; }
        label nAbove() const { return nAbove_; }
        void resetCounts() const { nBelow_ = 0; nAbove_ = 0; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
