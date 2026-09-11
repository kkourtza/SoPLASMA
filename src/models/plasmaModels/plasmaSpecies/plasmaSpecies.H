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
    Foam::plasmaSpecies
    
Description
    Core class for plasma species data management.

    This class provides the minimal data structures required to represent
    plasma species in OpenFOAM-based plasma simulations. It defines which
    species exist, their basic properties (e.g. charge number), and
    constructs and owns the corresponding volume fields such as number
    density.

    plasmaSpecies does not implement any physical models (transport,
    reactions, electron energy, etc.). Instead, it acts as a central storage
    and access point for all species-related quantities. Other physics
    modules operate on these fields but do not create their own copies.

    Only a single plasmaSpecies object should exist on a given mesh.

SourceFiles
    plasmaSpeciesI.H
    plasmaSpecies.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaSpecies_H
#define plasmaSpecies_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "PtrList.H"
#include "labelList.H"
#include "wordList.H"
#include "HashTable.H"

#include "SoPLASMAConstants.H"
#include "electromagneticsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class plasmaSpecies Declaration
\*---------------------------------------------------------------------------*/

class plasmaSpecies
:
    public IOdictionary
{
private:

    //- Disallow copy
    plasmaSpecies(const plasmaSpecies&) = delete;
    void operator=(const plasmaSpecies&) = delete;

protected:

    // Protected Data

        //- Mesh reference
        const fvMesh& mesh_;

        //- Reference to the electromagnetics model
        electromagneticsModel& em_;

        //- Active spcecies

            //- Number of active species
            label nSpecies_;

            //- Active species names
            wordList speciesNames_;

            //- Active species index-name mapping
            HashTable<label> speciesIDs_;

            //- Mass per active species [kg]
            PtrList<dimensionedScalar> speciesMasses_;

            //- Charge per active species [C]
            PtrList<dimensionedScalar> speciesCharges_;

            //- Charge number Z_i per active species (-1, +1, etc)
            List<scalar> speciesChargeNumbers_;

            //- Number density per active species [1/m3]
            PtrList<volScalarField> numberDensities_;

            //- Minimum number density clamp per active species [1/m3]
            List<scalar> speciesMinNumberDensities_;

            //- Property dictionaries per active species
            HashTable<dictionary> speciesDicts_;

            //- Default active species properties
            dictionary defaultSpeciesDict_;

        //- Background gas

            //- Background gas name
            word backgroundName_;

            //- Background gas number density N [1/m3]
            dimensionedScalar backgroundDensity_;

            //- True when activeSpecies was derived from the mechanism
            bool fromMechanism_ = false;

            //- Charge number and mass [kg] per species, as the mechanism
            //  declares them. Used only where the case does not state its own,
            //  so an explicit value always wins.
            HashTable<scalar> mechCharge_;
            HashTable<scalar> mechMass_;

            //- How derived ions are transported: immobile | driftDiffusion
            word ionTransport_ = "immobile";
            word ionFluxScheme_ = "standard";

            //- Default transport for derived NEUTRAL species. Same vocabulary
            //  as ionTransport deliberately: `immobile` and `diffusion` mean
            //  the same thing whichever class of species they are applied to,
            //  and a second set of words for the same concepts would be one
            //  more thing to look up.
            word neutralTransport_ = "immobile";

            //- D0 [m^2/s] per neutral species at Tref, pref, from the
            //  mechanism's `speciesDiffusivity` block
            HashTable<scalar> mechDiff_;

            //- Temperature exponent per species, fitted by mechc. Per species
            //  because Omega_D(T*) differs between them; a single shared value
            //  would be right for at most one.
            HashTable<scalar> mechDiffExp_;

            //- Reference state for the above
            scalar diffTref_ = 300.0;
            scalar diffPref_ = 1.0e5;

            //- Total neutral number density field [1/m3] 
            //- Sum of nBackground_ + all active neutrals
            volScalarField totalNeutralDensity_;

            //- Background gas dictionary
            dictionary backgroundDict_;

            //- The electron energy closure, `LFA` or `LMEA`. RESOLVED ONCE,
            //  here, from the global `electronEnergyModel` key (with the
            //  legacy per-species spelling accepted for one release).
            //
            //  Everything downstream -- whether the n_eps equation is solved,
            //  which field the electron transport coefficients are keyed on,
            //  which field the reaction rates are keyed on -- ASKS this and
            //  never re-derives the condition. Three places used to derive it
            //  independently, which is how a run came to be half-LMEA.
            //
            //  ORTHOGONAL to gas heating (`backgroundGas/energy/solve`):
            //  all four combinations are supported and neither implies the
            //  other. An earlier version keyed LMEA off `energy { solve }`,
            //  so enabling LMEA silently switched on the gas temperature.
            word electronEnergyModel_;

        //- Groups

            //- Index for the electron species
            label electronSpeciesID_;

            //- Indices for all ions (positive + negative)
            labelList ionSpeciesIDs_;

            //- Indices for positive ions
            labelList positiveIonSpeciesIDs_;

            //- Indices for negative ions
            labelList negativeIonSpeciesIDs_;

            //- Indices for active neutrals
            labelList neutralSpeciesIDs_;

            //- Indices for all charged species (ions + electrons)
            labelList chargedSpeciesIDs_;

            //- Indices for mobile species
            labelList mobileSpeciesIDs_;

            //- Indices for immobile species
            labelList immobileSpeciesIDs_;

            // REMOVED 2026-09-01: the five per-species TEMPERATURE groupings
            // (constantTemperature / dynamicTemperature /
            // followBackgroundTemp / solveEnergy / fieldTemperature).
            //
            // They were populated from a per-species `energyModel` key whose
            // accepted values -- isothermal / backgroundGas / localField /
            // solveEnergy -- were a SECOND, incompatible vocabulary for the
            // same key name that plasmaEnergy and plasmaTransport read with
            // {gasTemperature, isothermal, localEnergy, localField}. Only two
            // (`localField` was REMOVED 2026-09-01: it interpolated a
            // hand-supplied T(E/N) table and published T_<sp>, a field
            // nothing reads -- plasmaEnergy::T(label) has no callers --
            // and the sweep already writes meanEnergy_vs_reducedE. Under
            // LMEA T_e follows from the solved meanE; under LFA no
            // electron temperature is wanted at all.
            //
            // `isothermal` was REMOVED the same day and for the same
            // reason -- a constant T_<sp> that nothing reads -- and it is
            // NOT the mechanism for vibrational non-equilibrium: that is
            // the dynamic eVib reservoir in plasmaEnergy, drained by
            // Millikan-White/Popov V-T relaxation. A hand-set constant
            // could not represent an evolving T_vib in any case, and no
            // T_vib field or consumer exists. Selectable models are now
            // `gasTemperature` and `localEnergy`.)
            // values were common, and `gasTemperature` -- the documented
            // default -- matched no branch here, so a species using the
            // default fell into no group at all.
            //
            // All five accessors were called ZERO times outside their own
            // declarations, so the grouping decided nothing: it existed only
            // to confuse the person writing the case. Gas temperature is a
            // GLOBAL property, resolved once from backgroundGas/energy, and
            // the electron closure is `electronEnergyModel`. Do not
            // reintroduce a per-species temperature key without a consumer.

            //- Species that participate in reactions
            labelList reactingSpeciesIDs_;

            //- Species that are inert (transport only, no chemical reactions)
            labelList nonReactingSpeciesIDs_;


    // Private Member Functions

        //- Fill speciesNames_, mechCharge_ and mechMass_ from the compiled
        //  mechanism dictionary, for `activeSpecies fromMechanism;`
        void readMechanismSpecies();

        //- Resolve electronEnergyModel_ into `LFA` or `LMEA`.
        //
        //  MUST run AFTER electronSpeciesID_ is assigned (plasmaSpecies.C,
        //  in the species loop), because the legacy per-species spelling can
        //  only be read off the electron's own sub-dictionary. It is a second
        //  pass over the dictionary for exactly that reason.
        void resolveElectronEnergyModel();

        //- Point the ELECTRON's drift-diffusion coefficients at the field the
        //  closure implies, and reject a contradiction.
        //
        //  This closes the last silent half-LMEA hole. The chemistry key and
        //  the energy model both follow electronEnergyModel_, but the
        //  electron's own driftDiffusionCoeffs -- which drive the n_e flux --
        //  were hand-maintained, so LMEA could run its transport on `reducedE`
        //  with no error and plausible output.
        //
        //  ELECTRON ONLY. Ion coefficients are functions of E/N or |E| and are
        //  never touched; the ion blocks deliberately omit `lookupVariable` so
        //  that MechanismProperty's `reducedE` default applies.
        void deriveElectronTransportKey();


public:

    // Static Member Functions

        //- Species list implied by `activeSpecies fromMechanism;`.
        //
        //  Static and dictionary-only because plasmaCreateSpeciesFields needs
        //  the same list before a mesh exists, and it must create exactly the
        //  fields the solver will look for. Optionally returns the charge
        //  numbers and masses [kg] the mechanism declares.
        static wordList speciesFromMechanism
        (
            const dictionary& speciesDict,
            HashTable<scalar>* chargeOut = nullptr,
            HashTable<scalar>* massOut = nullptr,
            word* ionTrOut = nullptr,
            word* fluxOut = nullptr,
            word* neutralTrOut = nullptr,
            HashTable<scalar>* diffOut = nullptr,
            scalar* diffTrefOut = nullptr,
            scalar* diffPrefOut = nullptr,
            HashTable<scalar>* diffExpOut = nullptr
        );


        //- D0 [m^2/s] per neutral species at the mechanism's reference state.
        //  Empty for mechanisms compiled before mechc emitted it.
        inline const HashTable<scalar>& mechanismDiffusivity() const
        {
            return mechDiff_;
        }

        //- Reference temperature for mechanismDiffusivity()
        inline scalar diffusivityTref() const { return diffTref_; }

        //- Temperature exponent per species for mechanismDiffusivity()
        inline scalar diffusivityExponent(const word& s) const
        {
            return mechDiffExp_.lookup(s, 1.68);
        }

    //- Runtime type information
    TypeName("plasmaSpecies");

    // Constructors

        //- Construct from mesh
        explicit plasmaSpecies
        (
            const fvMesh& mesh,
            electromagneticsModel& em
        );

    //- Destructor
    virtual ~plasmaSpecies() = default;

    // Member functions
        
        //- Update only the volume charge density (instantaneous)
        void updateChargeDensity();

        //- INVARIANT: the chargeDensity field Poisson uses must equal
        //  sum_i q_i n_i for the CURRENT densities.
        //
        //  This is not a duplicate of the step audit. The audit catches a
        //  skipped CALL; this catches a wrong RESULT -- a missing charged
        //  species, a sign error, a clamp applied to a density AFTER the
        //  charge update, an aliased field. Either can freeze Poisson's
        //  source, and on 2026-09-06 a skipped call did, for 364670 steps.
        //
        //  Call it ONCE PER ACCEPTED STEP, after the last corrector: within
        //  the corrector loop the source is legitimately one iteration behind
        //  the densities, so a strict check there would fire on correct runs.
        //  A violation is FATAL and names the worst cell. CLAUDE.md rule 27.
        void verifyChargeDensity(const label everyN = 1) const;

        //- Clamp all number density fields to their minimum values
        void clampNumberDensities();

        //- Clamp a single species number density field
        void clampNumberDensity(label i);

    // Access

        //- Metadata access
        inline const fvMesh& mesh() const;
        inline label nSpecies() const;
        inline const wordList& speciesNames() const;
        inline const word& speciesName(label i) const;

        //- Electromagnetics model access
        inline electromagneticsModel& em() const;

        //- Get species index from name
        inline label speciesID(const word& name) const;

        //- Physical properties access
        inline const List<scalar>& speciesChargeNumbers() const;
        inline scalar speciesChargeNumber(label i) const;
        inline const PtrList<dimensionedScalar>& speciesCharges() const;
        inline const dimensionedScalar& speciesCharge(label i) const;
        inline const PtrList<dimensionedScalar>& speciesMasses() const;
        inline const dimensionedScalar& speciesMass(label i) const;
        inline const PtrList<volScalarField>& numberDensities() const;
        inline PtrList<volScalarField>& numberDensities();
        inline const volScalarField& numberDensity(label i) const;
        inline const List<scalar>& speciesMinNumberDensities() const;
        inline scalar speciesMinNumberDensity(label i) const;
        inline volScalarField& numberDensity(label i);

        //- Dictionary access
        inline const dictionary& speciesDict(label i) const;
        inline const dictionary& speciesDict(const word& name) const; 
        inline const HashTable<dictionary>& speciesDicts() const;
        inline const dictionary& defaultSpeciesDict() const;

        //- Background gas access
        inline const word& backgroundName() const;
        inline const dimensionedScalar& backgroundDensity() const;
        inline const volScalarField& totalNeutralDensity() const;
        inline const dictionary& backgroundDict() const;

        //- Group access
        inline label electronSpeciesID() const;

        //- The electron energy closure: `LFA` or `LMEA`.
        inline const word& electronEnergyModel() const;

        //- Is the local MEAN ENERGY approximation active? The one question
        //  every downstream consumer should ask; do not re-derive it from
        //  dictionary keys.
        inline bool isLMEA() const;

        //- The field the electron coefficients and reaction rates are keyed
        //  on: `meanE` under LMEA, `reducedE` under LFA. Determined entirely
        //  by the closure, so it is derived and never read from a case.
        inline const word& electronLookupKey() const;
        inline const labelList& ionSpeciesIDs() const;
        inline const labelList& positiveIonSpeciesIDs() const;
        inline const labelList& negativeIonSpeciesIDs() const;
        inline const labelList& neutralSpeciesIDs() const;
        inline const labelList& chargedSpeciesIDs() const;

        inline const labelList& mobileSpeciesIDs() const;
        inline const labelList& immobileSpeciesIDs() const;
        inline const labelList& reactingSpeciesIDs() const;
        inline const labelList& nonReactingSpeciesIDs() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaSpeciesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
