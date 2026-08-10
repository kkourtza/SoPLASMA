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

#include "plasmaReactionRates.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "dictionary.H"
#include <fstream>
#include <cstdlib>

// * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //

void Foam::plasmaReactionRates::readMechanism(const fileName& dictPath)
{
    IFstream is(dictPath);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open mechanism dictionary " << dictPath << nl
            << "It is written by the BoltzmannSolver compiler:" << nl
            << "    mechc.py <mechanism>.yaml -o <dir>/" << nl
            << exit(FatalError);
    }

    dictionary dict(is);

    mechanismName_ = dict.get<word>("mechanismName");
    // Read as a STRING, not a word: the hash is hex and may begin with a
    // digit, which an OpenFOAM word may not. mechc quotes it for the same
    // reason -- a quoted token is a string, and asking for a word here would
    // reject the very file that was written to be readable.
    mechanismHash_ = string(dict.get<string>("mechanismHash"));

    // Required, not optional. Defaulting it would let a mechanism compiled by
    // an older mechc load against code that assumes the declaration exists,
    // and the failure mode is an electron that silently matches nothing.
    electronSpecies_ = dict.get<word>("electronSpecies");

    const label nProcesses = dict.get<label>("nProcesses");

    // `processes` lists only the entries that are fluid reactions; the manifest
    // also contains energy-only channels (real EEDF energy sinks with no
    // species) and density-scaled curves. So this list is SHORTER than
    // nProcesses, and the indices in it are the canonical manifest indices --
    // deliberately not contiguous. Nothing here may assume otherwise.
    const List<dictionary> entries(dict.lookup("processes"));

    reactions_.setSize(entries.size());
    forAll(entries, i)
    {
        const dictionary& e = entries[i];
        reaction& r = reactions_[i];
        r.index     = e.get<label>("index");
        r.id        = e.get<word>("id");
        r.target    = e.get<word>("target");
        r.reactants = e.get<wordList>("reactants");
        r.products  = e.get<wordList>("products");
        r.collider  = e.getOrDefault<word>("collider", word::null);
        r.rateScale = e.getOrDefault<scalar>("rateScale", 1.0);
        r.threshold = e.getOrDefault<scalar>("threshold", 0.0);

        if (r.index < 0 || r.index >= nProcesses)
        {
            FatalErrorInFunction
                << "Reaction '" << r.id << "' has manifest index " << r.index
                << " outside [0," << nProcesses << ")" << nl
                << exit(FatalError);
        }
    }

    // Charges, so that discarding a charged product can be made fatal. Absent
    // in older dictionaries, in which case the check degrades to a warning.
    if (dict.found("speciesCharge"))
    {
        const dictionary& cd = dict.subDict("speciesCharge");
        for (const entry& e : cd)
        {
            charge_.insert(e.keyword(), readLabel(e.stream()));
        }
    }

    Info<< "plasmaReactionRates: mechanism " << mechanismName_
        << " [" << mechanismHash_.c_str() << "], "
        << reactions_.size() << " electron-impact reactions of "
        << nProcesses << " processes" << endl;
}


Foam::label Foam::plasmaReactionRates::chargeOf(const word& name) const
{
    const auto it = charge_.cfind(name);
    return it.good() ? *it : 0;
}


void Foam::plasmaReactionRates::checkHash(const fileName& tablePath) const
{
    // genMechTables writes the mechanism name and hash into the first line of
    // every table as a comment. Reading it is what makes a stale pair
    // detectable: the tables and the dictionary must come from the same master
    // file, or the indices in one do not mean what the other says they mean.
    IFstream is(tablePath);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open rate table " << tablePath << nl
            << "Generate the tables with:" << nl
            << "    genMechTables <mech>.mech.json --EN <lo>:<hi>:<n> -o <dir>/"
            << nl << exit(FatalError);
    }

    string firstLine;
    is.getLine(firstLine);

    if (firstLine.find(mechanismHash_) == std::string::npos)
    {
        FatalErrorInFunction
            << "Mechanism hash mismatch." << nl
            << "    dictionary : " << mechanismHash_ << nl
            << "    table      : " << tablePath << nl
            << "    header     : " << firstLine << nl << nl
            << "The tables were generated from a different mechanism than the"
            << " dictionary describes." << nl
            << "Process indices would not agree. Regenerate both from the same"
            << " master file." << nl
            << exit(FatalError);
    }
}


void Foam::plasmaReactionRates::buildEvaluators(const fileName& tableDir)
{
    rate_.setSize(reactions_.size());
    k_.setSize(reactions_.size());
    alpha_.setSize(reactions_.size());
    alphaN_.setSize(reactions_.size());

    const dimensionSet rateDims (0, 3, -1, 0, 0, 0, 0);   // m^3/s
    const dimensionSet alphaDims(0, 2,  0, 0, 0, 0, 0);   // m^2, reduced

    forAll(reactions_, i)
    {
        const reaction& r = reactions_[i];
        const fileName tablePath = tableDir/("k_" + r.id + "_vs_" + tableKey_);

        checkHash(tablePath);
        rate_.set(i, new plasmaRateTable(tablePath, bounds_));

        k_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "k_" + r.id, mesh_.time().timeName(), mesh_,
                    IOobject::NO_READ, IOobject::NO_WRITE
                ),
                mesh_, dimensionedScalar(rateDims, Zero)
            )
        );

        // Townsend coefficient where the reaction has one. Presence is decided
        // by whether genMechTables wrote the file -- ionising and attaching
        // reactions only -- rather than by re-deriving the reaction kind here,
        // so the two tools cannot disagree about which is which.
        for (const word& pre : {word("alpha_"), word("eta_")})
        {
            const fileName aPath = tableDir/(pre + r.id + "_vs_" + tableKey_);
            if (!isFile(aPath)) continue;

            checkHash(aPath);
            alpha_.set(i, new plasmaRateTable(aPath, bounds_));
            alphaN_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        pre + r.id, mesh_.time().timeName(), mesh_,
                        IOobject::NO_READ, IOobject::NO_WRITE
                    ),
                    mesh_, dimensionedScalar(alphaDims, Zero)
                )
            );
            break;
        }
    }
}


void Foam::plasmaReactionRates::readRefreshControl(const dictionary& dict)
{
    // Default is everyStep: correct by construction, and the only sensible
    // default when the cost of a refresh is not yet known for the user's case.
    const word mode = dict.getOrDefault<word>("eedfRefresh", "never");

    if (mode == "never")
    {
        refreshMode_ = rmNever;
    }
    else if (mode == "interval")
    {
        refreshMode_ = rmInterval;
        refreshInterval_ = dict.get<label>("refreshInterval");
        if (refreshInterval_ < 1)
        {
            FatalErrorInFunction
                << "refreshInterval must be >= 1, got " << refreshInterval_
                << exit(FatalError);
        }
    }
    else if (mode == "time")
    {
        refreshMode_ = rmTime;
        refreshDeltaT_ = dict.get<scalar>("refreshDeltaT");
        if (refreshDeltaT_ <= 0)
        {
            FatalErrorInFunction
                << "refreshDeltaT must be > 0, got " << refreshDeltaT_
                << exit(FatalError);
        }
    }
    else if (mode == "relativeChange")
    {
        refreshMode_ = rmRelativeChange;
        refreshTolerance_ =
            dict.getOrDefault<scalar>("refreshTolerance", 0.05);
        if (refreshTolerance_ <= 0)
        {
            FatalErrorInFunction
                << "refreshTolerance must be > 0, got " << refreshTolerance_
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown eedfRefresh mode '" << mode << "'" << nl
            << "Valid: never | interval | time | relativeChange" << nl
            << exit(FatalError);
    }

    Info<< "plasmaReactionRates: eedfRefresh " << mode;
    if (refreshMode_ == rmInterval)      Info<< " (every " << refreshInterval_ << " steps)";
    if (refreshMode_ == rmTime)          Info<< " (every " << refreshDeltaT_ << " s)";
    if (refreshMode_ == rmRelativeChange) Info<< " (tolerance " << refreshTolerance_ << ")";
    Info<< endl;
}


bool Foam::plasmaReactionRates::eedfRefreshDue()
{
    // `never` is the default and the normal mode: tables built once, offline,
    // and reused. Nothing below applies.
    if (refreshMode_ == rmNever)
    {
        return false;
    }

    const label idx = mesh_.time().timeIndex();

    // A mesh change forces one: a cached field sized to a previous mesh is not
    // merely stale, it is wrong.
    if (lastRefreshIndex_ < 0 || lastLookup_.size() != mesh_.nCells())
    {
        return true;
    }

    switch (refreshMode_)
    {
        case rmNever:
            return false;

        case rmInterval:
            return (idx - lastRefreshIndex_) >= refreshInterval_;

        case rmTime:
            return (mesh_.time().value() - lastRefreshTime_) >= refreshDeltaT_;

        case rmRelativeChange:
        {
            // Refresh when the keyed field has moved anywhere by more than the
            // tolerance. The comparison is relative to the LOCAL value, so a
            // quiescent region does not force refreshes on account of a single
            // active cell elsewhere being large -- and reduced across
            // processors, because the decision must be the same on all of them
            // or they fall out of step.
            const volScalarField& lut =
                mesh_.lookupObject<volScalarField>(lookupVariable_);

            scalar worst = 0;
            forAll(lut, c)
            {
                const scalar ref = max(mag(lastLookup_[c]), SMALL);
                worst = max(worst, mag(lut[c] - lastLookup_[c])/ref);
            }
            reduce(worst, maxOp<scalar>());
            return worst >= refreshTolerance_;
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasmaReactionRates::plasmaReactionRates
(
    const fvMesh& mesh,
    const fileName& mechanismDict,
    const fileName& tableDir,
    const word& lookupVariable,
    const word& tableKey,
    const dictionary& refreshDict
)
:
    mesh_(mesh),
    lookupVariable_(lookupVariable),
    tableKey_(tableKey)
{
    readMechanism(mechanismDict);
    bounds_ = plasmaRateTable::boundsFromWord
    (
        refreshDict.getOrDefault<word>("outOfBounds", "extrapolate")
    );
    buildEvaluators(tableDir);
    readRefreshControl(refreshDict);
    // The first correct() happens here, when the lookup field is still
    // zero-initialised. Suppress the range report for it -- reporting on a
    // field that has not been computed yet is worse than not reporting: it
    // once printed a reassuring "within range" while every lookup was in fact
    // clamping to zero. Re-arm afterwards so the first REAL step reports.
    rangeReported_ = true;
    correct();
    rangeReported_ = false;
    eedfRefreshed();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plasmaReactionRates::reportRange() const
{
    if (rangeReported_ || rate_.empty()) return;
    if (!mesh_.foundObject<volScalarField>(lookupVariable_)) return;

    const scalarField& f =
        mesh_.lookupObject<volScalarField>(lookupVariable_).primitiveField();
    if (f.empty()) return;

    const scalar lo = gMin(f), hi = gMax(f);
    const scalar tLo = rate_[0].xMin(), tHi = rate_[0].xMax();

    label nAbove = 0, nBelow = 0;
    forAll(f, c)
    {
        if (f[c] > tHi) ++nAbove;
        if (f[c] < tLo) ++nBelow;
    }
    reduce(nAbove, sumOp<label>());
    reduce(nBelow, sumOp<label>());
    label nTot = f.size();
    reduce(nTot, sumOp<label>());

    Info<< "plasmaReactionRates: " << lookupVariable_ << " spans ["
        << lo << ", " << hi << "], tables cover [" << tLo << ", " << tHi << "]"
        << endl;

    // Report the FRACTION out of range, not the extremes. A single singular
    // cell reporting a huge maximum is alarming and useless; "30 cells of 1.1
    // million" is what tells you whether to care.
    if (nAbove || nBelow)
    {
        Info<< "    out of range: " << nAbove << " above, " << nBelow
            << " below, of " << nTot << " cells ("
            << 100.0*scalar(nAbove + nBelow)/max(nTot, label(1)) << " %)" << nl
            << "    handling: "
            << (bounds_ == plasmaRateTable::bhExtrapolate
                    ? "extrapolate (power law from the last tabulated interval)"
                    : "clamp (hold the end value)")
            << endl;

        if (nBelow && lo > 0 && tLo > 0 && lo < 1e-3*tLo)
        {
            WarningInFunction
                << "the field minimum is orders of magnitude below the table."
                << " Suspect a units mismatch: reducedE is in V m^2, while E/N"
                << " is conventionally quoted in Td (1 Td = 1e-21 V m^2)."
                << endl;
        }
    }
    rangeReported_ = true;
}


void Foam::plasmaReactionRates::correct()
{
    // Unthrottled by design: k_j must track the field every step, or the
    // coefficients stay frozen at t=0 while the discharge moves. This is an
    // interpolation, not a Boltzmann solve -- it is the cheap half.
    if (!mesh_.foundObject<volScalarField>(lookupVariable_)) return;

    const scalarField& lut =
        mesh_.lookupObject<volScalarField>(lookupVariable_).primitiveField();

    forAll(reactions_, i)
    {
        rate_[i].value(lut, k_[i].primitiveFieldRef());
        k_[i].correctBoundaryConditions();
        if (alpha_.set(i))
        {
            alpha_[i].value(lut, alphaN_[i].primitiveFieldRef());
            alphaN_[i].correctBoundaryConditions();
        }
    }

    reportRange();
}


void Foam::plasmaReactionRates::eedfRefreshed()
{
    lastRefreshIndex_ = mesh_.time().timeIndex();
    lastRefreshTime_  = mesh_.time().value();
    ++nRefresh_;

    // Snapshot the keyed field so relativeChange has a baseline next step.
    if (mesh_.foundObject<volScalarField>(lookupVariable_))
    {
        lastLookup_ =
            mesh_.lookupObject<volScalarField>(lookupVariable_).primitiveField();
    }
}


Foam::label Foam::plasmaReactionRates::indexOf(const word& id) const
{
    forAll(reactions_, i)
    {
        if (reactions_[i].id == id)
        {
            return i;
        }
    }
    return -1;
}


const Foam::volScalarField&
Foam::plasmaReactionRates::k(const word& id) const
{
    const label i = indexOf(id);
    if (i < 0)
    {
        // Deliberately fatal. A mistyped id would otherwise contribute nothing
        // and leave a mechanism that looks complete but is quietly missing a
        // reaction.
        FatalErrorInFunction
            << "No reaction with id '" << id << "' in mechanism "
            << mechanismName_ << nl << "Available:" << nl;
        forAll(reactions_, j)
        {
            FatalError << "    " << reactions_[j].id << nl;
        }
        FatalError << exit(FatalError);
    }
    return k_[i];
}


// ************************************************************************* //
