/*---------------------------------------------------------------------------*\
  File: ionInducedSEE.H
  Part of: SoPLASMA

Class
    Foam::emissionModels::ionInducedSEE

Description
    Secondary electron emission by ION IMPACT:

        Gamma_e = SUM_i gamma_i max(0, Gamma_i . n)

    over the POSITIVE ions with a wall flux. Only inward flux counts, and only
    positive ions: a negative ion arriving is not a neutralisation event of
    this kind.

    THIS IS THE MECHANISM THAT DECIDES SELF-SUSTAINMENT. In a barrier discharge
    gamma is what turns a single avalanche into a sustained discharge, so a
    value fifty times too large manufactures sustainment the case never asked
    for -- and the run looks entirely healthy. Everything below is arranged so
    that the value used is stated, sourced and reported, never inferred.

    THREE WAYS TO GET gamma, chosen by `yield`:

    yield constant   (default)
        One number for every ion, with optional per-species overrides. This is
        what the wall-flux conditions have always done through `defaultSEEC`,
        and it remains the default because it is the honest choice when the
        surface state is unknown -- which it usually is.

    yield hagstrum
        PHYSICS-RIGID, per ion and per surface. Potential (Auger) emission,
        which is the channel that matters below ~1 keV, i.e. in every
        collisional discharge this solver targets. For a METAL of work function
        phi:

            gamma_i ~ 0.016 (0.8 E_ion - 2 phi)

        and the process is FORBIDDEN, gamma = 0, when 0.8 E_ion < 2 phi -- Auger
        neutralisation cannot lift a second electron over the barrier. Raizer,
        "Gas Discharge Physics" (1991) sec. 4, after Hagstrum, Phys. Rev. 96
        (1954) 336.

        For an INSULATOR phi does not apply. Emission requires crossing the band
        gap AND the electron affinity, E_ion > Eg + 2 chi, and the yield is
        scaled by how far above that threshold the ion is. This is why MgO is
        exceptional and why the model must be told which kind of surface it is
        on -- `surface metal` or `surface insulator`.

        E_ion comes from $SoPLASMA_ETC/materials/ionisationEnergies, indexed by
        the ion's PARENT NEUTRAL: `N2p` -> `N2`. phi from workFunctions, Eg and
        chi from dielectrics. A species whose parent is not in the library is
        FATAL, not defaulted -- a fabricated ionisation energy is a fabricated
        discharge.

        NOTE WHAT THIS DOES NOT COVER: kinetic emission, which switches on above
        roughly 0.5-1 keV impact energy. Atmospheric-pressure ions are
        collisional and arrive with a few eV to tens of eV, so it is negligible
        there -- but it is NOT negligible in a low-pressure sheath, and this
        model does not represent it.

    yield table
        gamma from a tabulated file against reduced field E/N. THE CAPABILITY
        SHIPS WITH NO DATA, deliberately. Published gamma_eff(E/N) -- Phelps &
        Petrovic, PSST 8 (1999) R21 being the standard -- is FITTED PER
        GAS/CATHODE PAIR and lumps ion impact, fast-neutral impact,
        photoemission and electron backscatter into one effective coefficient.
        An argon-on-copper fit is not an air-on-acrylic coefficient, and
        shipping one as a default would be worse than a constant because it
        would look authoritative. Supply your own table for your own pair, and
        say where it came from.

SourceFiles
    ionInducedSEE.C

\*---------------------------------------------------------------------------*/

#ifndef ionInducedSEE_H
#define ionInducedSEE_H

#include "emissionModel.H"
#include "HashTable.H"

namespace Foam
{
namespace emissionModels
{

/*---------------------------------------------------------------------------*\
                          Class ionInducedSEE Declaration
\*---------------------------------------------------------------------------*/

class ionInducedSEE
:
    public emissionModel
{
    //- constant | hagstrum | table
    const word yield_;

    //- `yield constant`: the coefficient for every ion
    const scalar gamma_;

    //- `yield constant`: per-species overrides
    const dictionary speciesGamma_;

    //- `yield hagstrum`: metal | insulator
    const word surface_;

    //- Resolved gamma per species NAME, filled once and reported
    mutable HashTable<scalar> resolved_;

    //- Whether the resolution has run
    mutable bool resolvedOnce_;


    // Private Member Functions

        //- gamma for one ion species, by whichever route `yield` selects
        scalar gammaFor(const word& ionName) const;

        //- The parent neutral of an ion: `N2p` -> `N2`, `O2m` -> `O2`
        static word parentNeutral(const word& ionName);

        //- Fill resolved_ for every positive ion with a wall flux
        void resolve() const;


public:

    TypeName("ionInducedSEE");

    ionInducedSEE(const fvPatch& p, const dictionary& dict);

    virtual ~ionInducedSEE() = default;

    virtual tmp<scalarField> emittedFlux() const;

    virtual void report(Ostream& os) const;
};


} // End namespace emissionModels
} // End namespace Foam

#endif

// ************************************************************************* //
