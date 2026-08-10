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
\*---------------------------------------------------------------------------*/

#include "plasmaRateTable.H"
#include "error.H"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

Foam::plasmaRateTable::boundsHandling
Foam::plasmaRateTable::boundsFromWord(const word& w)
{
    if (w == "clamp")       return bhClamp;
    if (w == "extrapolate") return bhExtrapolate;

    FatalErrorInFunction
        << "Unknown outOfBounds '" << w << "'" << nl
        << "Valid: clamp | extrapolate" << nl
        << exit(FatalError);

    return bhClamp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasmaRateTable::plasmaRateTable
(
    const fileName& path,
    boundsHandling bounds
)
:
    bounds_(bounds),
    path_(path)
{
    // Parsed with strtod rather than OpenFOAM's stream operators: the file
    // contains a "//" header and a bare "(" line, and the stream operator
    // raises a FatalIOError on a token it cannot read rather than simply
    // failing. A tolerant parse is required -- reading a table must not be able
    // to abort a run on a comment line.
    std::ifstream f(path);
    if (!f)
    {
        FatalErrorInFunction
            << "Cannot open rate table " << path << nl << exit(FatalError);
    }

    std::vector<scalar> xs, ys;
    std::string line;
    while (std::getline(f, line))
    {
        const auto o = line.find('(');
        if (o == std::string::npos) continue;

        const char* p = line.c_str() + o + 1;
        char* end = nullptr;
        const double a = std::strtod(p, &end);
        if (end == p) continue;                  // header or bare "("
        const char* q = end;
        const double b = std::strtod(q, &end);
        if (end == q) continue;                  // only one number: not a row

        xs.push_back(a);
        ys.push_back(b);
    }

    if (xs.size() < 2)
    {
        FatalErrorInFunction
            << "Rate table " << path << " has " << xs.size()
            << " usable rows; at least 2 are required" << nl << exit(FatalError);
    }

    x_.setSize(xs.size());
    y_.setSize(ys.size());
    forAll(x_, i) { x_[i] = xs[i]; y_[i] = ys[i]; }

    for (label i = 1; i < x_.size(); ++i)
    {
        if (x_[i] <= x_[i-1])
        {
            FatalErrorInFunction
                << "Rate table " << path << " is not strictly increasing at row "
                << i << " (" << x_[i-1] << " -> " << x_[i] << ")" << nl
                << exit(FatalError);
        }
    }

    // Power-law slope of the final interval, for extrapolation above the table.
    // Requires both ends positive: a power law through zero is undefined, and
    // inventing one would be worse than holding the value.
    const label n = x_.size();
    if (y_[n-1] > 0 && y_[n-2] > 0 && x_[n-1] > 0 && x_[n-2] > 0)
    {
        topSlope_ =
            (::log(y_[n-1]) - ::log(y_[n-2]))
          / (::log(x_[n-1]) - ::log(x_[n-2]));
        topIsPowerLaw_ = true;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::plasmaRateTable::value(const scalar x) const
{
    const label n = x_.size();

    // Below the table: hold. Rate coefficients fall to zero at low field and
    // the first row is already the zero-field value, so there is nothing
    // meaningful to extrapolate towards -- and a power law downwards would
    // diverge rather than vanish.
    if (x <= x_[0])
    {
        if (x < x_[0]) ++nBelow_;
        return y_[0];
    }

    if (x >= x_[n-1])
    {
        if (x > x_[n-1]) ++nAbove_;

        if (bounds_ == bhClamp || !topIsPowerLaw_ || x <= 0)
        {
            return y_[n-1];
        }

        // y = y_last * (x/x_last)^slope  -- linear in log-log, continuing the
        // trend of the final tabulated interval.
        return y_[n-1]*::exp(topSlope_*(::log(x) - ::log(x_[n-1])));
    }

    // Inside: linear in x, matching interpolationTable so that switching to
    // this class changes nothing where the table already covered the field.
    label lo = 0, hi = n - 1;
    while (hi - lo > 1)
    {
        const label mid = (lo + hi)/2;
        if (x_[mid] > x) hi = mid; else lo = mid;
    }
    const scalar t = (x - x_[lo])/(x_[hi] - x_[lo]);
    return y_[lo] + t*(y_[hi] - y_[lo]);
}


void Foam::plasmaRateTable::value
(
    const scalarField& in,
    scalarField& out
) const
{
    out.setSize(in.size());
    forAll(in, i)
    {
        out[i] = value(in[i]);
    }
}


// ************************************************************************* //
