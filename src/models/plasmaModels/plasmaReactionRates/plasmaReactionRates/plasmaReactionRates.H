/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

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
    Foam::plasmaReactionRates

Description:
    Electron-impact reaction rate coefficients, tabulated against the reduced
    electric field and evaluated per cell.

    Replaces the hard-coded Townsend fits in plasmaTransport::solve() with
    coefficients generated from a mechanism, so the chemistry is defined by a
    data file rather than by the source.

    The tables are written by the BoltzmannSolver tool `genMechTables` and the
    mechanism dictionary by its compiler `mechc`; both derive from ONE master
    mechanism file, so nothing here can drift from what the EEDF was solved
    with. See <BoltzmannSolver>/docs/soplasma_workflow.md.

    One plasmaPropertyEvaluator is held per reaction rather than a bespoke
    interpolator: `tabulated1D` already wraps OpenFOAM's interpolationTable and
    already accepts `reducedE` as a lookup variable, which is exactly what these
    tables are keyed on.

    UNITS -- the trap worth stating in the header, because getting it wrong is
    invisible. Rate coefficients here are UNWEIGHTED: k_j in m^3/s per unit
    TARGET density, the BOLSIG+/Cantera convention. The reaction rate is

        R_j = k_j * n_e * n_target                                   [m^-3 s^-1]

    so the caller supplies the target density. Do NOT multiply by a mole
    fraction: the tables already exclude it. (An earlier version of the solver
    returned x_j*k_j, which agrees with this to within a factor that looks
    entirely plausible -- 0.79 for N2 -- and is therefore easy to miss.)

    Staleness: the mechanism dictionary carries a hash of the master file, and
    every table carries the same hash in its header comment. They are compared
    on construction, because a table generated from a different mechanism has
    process indices that no longer mean what the dictionary says they mean.

SourceFiles:
    plasmaReactionRates.C
\*---------------------------------------------------------------------------*/

#ifndef plasmaReactionRates_H
#define plasmaReactionRates_H

#include "fvMesh.H"
#include "volFields.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "wordList.H"
#include "HashTable.H"
#include "plasmaRateTable.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class plasmaReactionRates Declaration
\*---------------------------------------------------------------------------*/

class plasmaReactionRates
{
public:

    //- One electron-impact reaction, as declared by the mechanism dictionary
    struct reaction
    {
        //- Canonical index in the mechanism manifest
        label index = -1;

        //- Mechanism process id, e.g. "EI_N2_ION_N2". The binding contract is
        //  this string; the index is derived and only used for fast access.
        word id;

        //- Target species name. This is NOT the whole left-hand side: it is
        //  the heavy collision partner, and it exists separately from
        //  `reactants` because the rate coefficient is per unit TARGET
        //  density, which is a question about the rate rather than about
        //  stoichiometry.
        word target;

        //- Full left-hand side, electron included. Explicit precisely because
        //  it used to be implicit: consumers were expected to remember that
        //  every electron-impact process consumes one electron, and the
        //  source-term assembly did not, so each excitation manufactured an
        //  electron from nothing.
        wordList reactants;

        //- Full right-hand side. The electron appears under the mechanism's
        //  declared name (see electronSpecies()), not as a bare "e".
        wordList products;

        //- Third body, for three-body processes; empty otherwise. A catalyst:
        //  it is on both sides and cancels, so it is deliberately NOT in
        //  reactants/products (which must balance on their own). It enters only
        //  as a density factor, which is what makes the rate third order.
        word collider;

        //- Unit conversion applied to the tabulated coefficient, so that
        //  k*rateScale is in m^6/s for a three-body process. 1 otherwise.
        //  Emitted as a number by mechc so this side never has to know what a
        //  `density-cm3` curve is.
        scalar rateScale = 1.0;

        //- True when the rate is third order and needs the collider density
        bool threeBody() const { return !collider.empty(); }

        //- Energy the electron loses [eV]; 0 where the mechanism gives none
        //  (attachment removes the electron rather than degrading it)
        scalar threshold = 0.0;
    };


    //- The electron under the name the mechanism uses in `reactants` and
    //  `products`. A case names its own electron species freely, so the two
    //  must be matched through this declaration rather than by assuming a
    //  spelling; code that tested for a literal "e" worked by coincidence.
    const word& electronSpecies() const { return electronSpecies_; }

    //- Charge of a mechanism species in units of e, or 0 if unknown.
    //  OpenFOAM cannot know the charge of a species it does not transport, so
    //  this comes from the mechanism dictionary. It is what makes discarding a
    //  charged product detectable.
    label chargeOf(const word& name) const;

    //- Charges of every species the mechanism touches
    const HashTable<label>& speciesCharge() const { return charge_; }


private:

    // Private Data

        const fvMesh& mesh_;

        //- Kept so a refresh can re-run the sweep with the same settings the
        //  start-up build used, changing only the conditions.
        fileName manifest_;
        fileName tableDir_;
        dictionary chemDict_;

        //- Mechanism name and the hash of the master file it was compiled from
        word mechanismName_;
        string mechanismHash_;

        //- Reactions, in canonical manifest order
        List<reaction> reactions_;

        //- Charge of each mechanism species, from the dictionary
        HashTable<label> charge_;

        //- Mechanism's name for the electron, from `electronSpecies` in the
        //  mechanism dictionary
        word electronSpecies_ = "Electron";

        //- One rate-coefficient table per reaction, same order.
        //  plasmaRateTable rather than the tabulated1D evaluator because
        //  OpenFOAM's interpolationTable cannot extrapolate, and clamping a
        //  rate coefficient above the tabulated range is a systematic bias
        //  rather than the neutral choice it looks like.
        PtrList<plasmaRateTable> rate_;

        //- Cached rate coefficient fields [m^3/s], same order
        PtrList<volScalarField> k_;

        //- Per-reaction Townsend coefficients [m^2], for `sourceModel townsend`.
        //  Only ionising and attaching reactions have one; the entry is unset
        //  for the rest, which is what hasTownsend() reports.
        //
        //    alpha_j/N = nu_j/N / v_drift = x_target * k_j / v_drift
        //
        //  Note these ARE mole-fraction weighted (they are nu/N) while k_ is
        //  not. That is precisely what makes the two formulations equal:
        //  (alpha_j/N) N v_drift n_e == k_j n_target n_e, verified to ratio
        //  1.000000 over 100-2000 Td.
        PtrList<volScalarField> alphaN_;
        PtrList<plasmaRateTable> alpha_;

        //- Field the tables are keyed on. "reducedE" today; switching to the
        //  mean electron energy is a change of this word plus the table set,
        //  because genMechTables writes both keys from the same sweep.
        word lookupVariable_;

        //- Suffix of the table FILES, which is the key they are tabulated
        //  against ("EN" or "meanE") -- distinct from lookupVariable_, which is
        //  the name of the OpenFOAM FIELD holding that key ("reducedE"). They
        //  are different namespaces and conflating them looks for
        //  k_<id>_vs_reducedE, which genMechTables never writes.
        word tableKey_;


    // EEDF refresh control
    //
    // TWO DIFFERENT THINGS, and conflating them is expensive in one direction
    // and wrong in the other:
    //
    //   correct()       re-interpolates k_j from the CURRENT lookup field.
    //                   Cheap, and it must happen every step -- otherwise k_j
    //                   stays frozen at its t=0 value while the discharge
    //                   moves. NOT throttled.
    //
    //   refreshEEDF()   re-solves the EEDF and regenerates the tables for the
    //                   current gas state and composition. Expensive: a
    //                   Boltzmann solve, not an interpolation. THIS is what the
    //                   control below governs.
    //
    // Default is `never`: the tables are built once, offline, from the master
    // mechanism, and reused for the whole run. That is the intended normal mode
    // of operation and it is exact as long as the composition and gas
    // temperature stay near the ones the tables were built at.
    //
    // The other modes exist because that assumption eventually fails -- the
    // EEDF is a functional of the whole mixture, so once excited states and
    // ions build up appreciably the tables drift out of validity:
    //
    //   never            build once, never re-solve (default)
    //   interval         re-solve every N timesteps
    //   time             re-solve every deltaT seconds of simulated time
    //   relativeChange   re-solve when the lookup field has moved anywhere on
    //                    the mesh by more than a given fraction
    //
    // In every mode other than `never` the first evaluation still happens at
    // construction, and a mesh change forces one, because a cached field sized
    // to a previous mesh is not stale but wrong.

        enum refreshMode { rmNever, rmInterval, rmTime, rmRelativeChange };

        refreshMode refreshMode_ = rmNever;

        //- Steps between refreshes (rmInterval)
        label refreshInterval_ = 1;

        //- Simulated time between refreshes (rmTime) [s]
        scalar refreshDeltaT_ = 0.0;

        //- Fractional change in the lookup field forcing a refresh
        //  (rmRelativeChange)
        scalar refreshTolerance_ = 0.05;

        //- Bookkeeping for the criteria above
        label  lastRefreshIndex_ = -1;
        scalar lastRefreshTime_ = -GREAT;
        scalarField lastLookup_;

        //- Refreshes actually performed; reported so the saving (or its
        //  absence) is visible rather than assumed
        label nRefresh_ = 0;

        //- Response outside the tabulated range, from `outOfBounds` in the
        //  chemistry dictionary: clamp | extrapolate (default extrapolate)
        plasmaRateTable::boundsHandling bounds_ =
            plasmaRateTable::bhExtrapolate;

        //- The one-time range report is deferred to the first call AFTER
        //  construction. Reporting from the constructor measures the
        //  zero-initialised field and reassuringly says "within range" about a
        //  field that has not been computed yet -- which is exactly what it
        //  did, while every lookup was silently clamping to zero.
        //- mutable: reportRange() is const (it only observes) but has to
        //  remember that it has already fired.
        mutable bool rangeReported_ = false;
        void reportRange() const;

        //- Read the refresh sub-dictionary
        void readRefreshControl(const dictionary& dict);

        //- Whether the criterion says a refresh is due now
        bool refreshDue();


    // Private Member Functions

        //- Read <mech>.foam: reaction list, ids, hash
        void readMechanism(const fileName& dictPath);

        //- Build one evaluator per reaction from its k_<id>_vs_<key> table
        void buildEvaluators(const fileName& tableDir);

        //- Compare the hash in a table header against the mechanism's, and
        //  FatalError if they differ -- a stale pair is not recoverable at
        //  runtime and must not be silently used.
        void checkHash(const fileName& tablePath) const;

        //- No copying
        plasmaReactionRates(const plasmaReactionRates&) = delete;
        void operator=(const plasmaReactionRates&) = delete;


public:

    // Constructors

        //- Construct from the mechanism dictionary and the table directory
        plasmaReactionRates
        (
            const fvMesh& mesh,
            const fileName& mechanismDict,
            const fileName& tableDir,
            const word& lookupVariable = "reducedE",
            const word& tableKey = "EN",
            const dictionary& refreshDict = dictionary::null
        );


    //- Destructor
    ~plasmaReactionRates() = default;


    // Member Functions

        //- Re-interpolate every rate coefficient from the current lookup
        //  field. Cheap; call once per timestep after that field is updated.
        //  Never throttled -- see the note on refresh control above.
        void correct();

        //- Whether the EEDF should be re-solved now, per the configured
        //  criterion. Always false in the default `never` mode.
        bool eedfRefreshDue();

        //- Re-solve the EEDF at the given composition and gas temperature,
        //  regenerate the tables and reload the interpolators. Returns false
        //  when refresh is `never`.
        //
        //  The caller supplies the composition because it owns the species
        //  densities; this class owns only the mechanism and its tables. It
        //  should be weighted by ELECTRON DENSITY rather than by volume: the
        //  rate coefficients matter where the electrons are, and a
        //  volume-average over a mostly undisturbed domain would report the
        //  initial gas no matter what the discharge had done to it.
        //
        //  Both halves are required and are done here: regenerating the tables
        //  without reloading them leaves the run using the old coefficients
        //  from memory -- a refresh at full cost that changes nothing.
        bool refreshEEDF(const HashTable<scalar>& composition, const scalar Tgas);

        //- Record that an EEDF refresh has just been performed, resetting the
        //  criterion. Call after acting on eedfRefreshDue().
        void eedfRefreshed();

        //- EEDF refreshes performed so far (0 in `never` mode)
        label nRefresh() const { return nRefresh_; }

        //- Number of electron-impact reactions
        label size() const { return reactions_.size(); }

        //- Reaction metadata, in canonical manifest order
        const List<reaction>& reactions() const { return reactions_; }

        //- Rate coefficient field for reaction i [m^3/s], UNWEIGHTED
        const volScalarField& k(const label i) const { return k_[i]; }

        //- Whether reaction i has a Townsend coefficient (ionising/attaching)
        bool hasTownsend(const label i) const { return alpha_.set(i); }

        //- Townsend coefficient for reaction i [m^2], reduced and
        //  mole-fraction weighted. Only valid where hasTownsend(i).
        const volScalarField& townsend(const label i) const { return alphaN_[i]; }

        //- Rate coefficient by mechanism id; FatalError if unknown, because a
        //  mistyped id would otherwise silently contribute nothing
        const volScalarField& k(const word& id) const;

        //- Index of a mechanism id, or -1
        label indexOf(const word& id) const;

        //- Mechanism identity, for reporting and for staleness checks
        const word& mechanismName() const { return mechanismName_; }
        const string& mechanismHash() const { return mechanismHash_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
