/*---------------------------------------------------------------------------*\
  File: ionInducedSEE.C
  Part of: SoPLASMA
\*---------------------------------------------------------------------------*/

#include "ionInducedSEE.H"
#include "addToRunTimeSelectionTable.H"
#include "materialLibrary.H"
#include "plasmaTransport.H"
#include "surfaceFields.H"
#include "volFields.H"

namespace Foam
{
namespace emissionModels
{
    defineTypeNameAndDebug(ionInducedSEE, 0);
    addToRunTimeSelectionTable(emissionModel, ionInducedSEE, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emissionModels::ionInducedSEE::ionInducedSEE
(
    const fvPatch& p,
    const dictionary& dict
)
:
    emissionModel(p, dict),
    yield_(dict.getOrDefault<word>("yield", "constant")),
    // 0.001, the contaminated-oxide / barrier figure. The SAME value and the
    // same reasoning as the wall-flux `defaultSEEC` default: the surfaces this
    // solver targets are air-exposed electrodes and dielectric barriers, both
    // contaminated, and 0.05 is a clean-metal-in-vacuum figure at the wrong end
    // of the range.
    gamma_(dict.getOrDefault<scalar>("gamma", 0.001)),
    speciesGamma_(dict.subOrEmptyDict("speciesGamma")),
    surface_(dict.getOrDefault<word>("surface", "metal")),
    resolved_(),
    resolvedOnce_(false)
{
    if (yield_ != "constant" && yield_ != "hagstrum" && yield_ != "table")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown `yield` `" << yield_ << "` on patch " << p.name()
            << "." << nl
            << "    Valid: constant (default) | hagstrum | table." << nl
            << exit(FatalIOError);
    }

    if (surface_ != "metal" && surface_ != "insulator")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown `surface` `" << surface_ << "`." << nl
            << "    Valid: metal (default) | insulator. The Auger criterion is"
               " DIFFERENT for the two:" << nl
            << "    a metal uses the work function, an insulator the band gap"
               " and electron affinity." << nl
            << exit(FatalIOError);
    }

    if (yield_ == "table")
    {
        FatalIOErrorInFunction(dict)
            << "`yield table` is not implemented yet." << nl << nl
            << "    The CAPABILITY is intended; NO DATA ships with it, and that"
               " is deliberate." << nl
            << "    Published gamma_eff(E/N) -- Phelps & Petrovic, PSST 8"
               " (1999) R21 -- is fitted PER" << nl
            << "    GAS/CATHODE PAIR and lumps ion impact, fast neutrals,"
               " photoemission and electron" << nl
            << "    backscatter into one effective number. An argon-on-copper"
               " fit is not an" << nl
            << "    air-on-acrylic coefficient." << nl << nl
            << "    Use `yield constant` with a value you can defend, or"
               " `yield hagstrum` for the" << nl
            << "    physics-rigid potential-emission estimate." << nl
            << exit(FatalIOError);
    }

    if (gamma_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "`gamma` must not be negative, got " << gamma_ << "." << nl
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::emissionModels::ionInducedSEE::parentNeutral
(
    const word& ionName
)
{
    // THE MECHANISM'S NAMING CONVENTION: a positive ion is the neutral with a
    // trailing `p`, a negative ion with a trailing `m`. N2p -> N2, O2m -> O2.
    // Nothing else is stripped: `Om` is atomic O-, whose parent is O.
    if (ionName.size() > 1)
    {
        const char last = ionName[ionName.size() - 1];

        if (last == 'p' || last == 'm')
        {
            return word(ionName.substr(0, ionName.size() - 1));
        }
    }

    return ionName;
}


Foam::scalar Foam::emissionModels::ionInducedSEE::gammaFor
(
    const word& ionName
) const
{
    if (yield_ == "constant")
    {
        // A per-species override WINS, exactly as speciesSEEC did.
        return speciesGamma_.getOrDefault<scalar>(ionName, gamma_);
    }

    // --- yield hagstrum ----------------------------------------------------
    const word parent(parentNeutral(ionName));

    const string ctx
    (
        "ionInducedSEE `yield hagstrum` on patch `" + patch_.name()
      + "`, ion `" + ionName + "`"
    );

    const scalar Eion = materialLibrary::get
    (
        parent, "E", ctx, "ionisationEnergies"
    );

    if (surface_ == "metal")
    {
        const scalar phi = workFunction();

        // Raizer 1991 sec. 4, after Hagstrum. FORBIDDEN below the threshold:
        // Auger neutralisation cannot lift a second electron over the barrier,
        // and a negative gamma is not a small gamma -- it is no emission.
        return max(0.0, 0.016*(0.8*Eion - 2.0*phi));
    }

    // --- insulator ---------------------------------------------------------
    // phi is NOT the right quantity. Emission must cross the band gap AND the
    // electron affinity: E_ion > Eg + 2 chi.
    if (!dict_.found("material"))
    {
        FatalIOErrorInFunction(dict_)
            << "`surface insulator` needs a `material`, to get its band gap and"
            << " electron affinity." << nl
            << "    Patch `" << patch_.name() << "`." << nl
            << exit(FatalIOError);
    }

    const word mat(dict_.get<word>("material"));

    const scalar Eg  = materialLibrary::get(mat, "bandGap", ctx);
    const scalar chi = materialLibrary::get(mat, "electronAffinity", ctx);

    const scalar threshold = Eg + 2.0*chi;

    // The same LINEAR form as the metal case with the barrier replaced by the
    // insulator one. It is a first-order estimate, not a measured yield, and
    // the report says so.
    return max(0.0, 0.016*(0.8*Eion - threshold));
}


void Foam::emissionModels::ionInducedSEE::resolve() const
{
    if (resolvedOnce_) return;
    resolvedOnce_ = true;

    const objectRegistry& db = patch_.boundaryMesh().mesh().thisDb();

    if (!db.foundObject<plasmaTransport>("plasmaTransport")) return;

    const plasmaSpecies& sp =
        db.lookupObject<plasmaTransport>("plasmaTransport").species();

    for (const label i : sp.positiveIonSpeciesIDs())
    {
        const word& name = sp.speciesName(i);
        resolved_.set(name, gammaFor(name));
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::emissionModels::ionInducedSEE::emittedFlux() const
{
    auto tflux = tmp<scalarField>::New(patch_.size(), Zero);
    scalarField& flux = tflux.ref();

    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    if (!mesh.thisDb().foundObject<plasmaTransport>("plasmaTransport"))
    {
        return tflux;
    }

    const plasmaSpecies& sp =
        mesh.thisDb().lookupObject<plasmaTransport>("plasmaTransport").species();

    resolve();

    const scalarField& magSf = patch_.magSf();

    // EXACTLY the sum the wall-flux condition has always formed:
    //   SUM over POSITIVE ions of gamma_i * max(0, Gamma_i / |Sf|)
    // Only positive ions, and only INWARD flux -- an ion leaving the surface
    // is not a neutralisation event.
    for (const label i : sp.positiveIonSpeciesIDs())
    {
        const word fluxName("particleFlux_" + sp.speciesName(i));

        if (!mesh.foundObject<surfaceScalarField>(fluxName)) continue;

        const fvsPatchScalarField& phiI =
            patch_.lookupPatchField<surfaceScalarField, scalar>(fluxName);

        const scalar g = resolved_.lookup(sp.speciesName(i), 0.0);

        flux += g*max(scalar(0), phiI/magSf);
    }

    return tflux;
}


void Foam::emissionModels::ionInducedSEE::report(Ostream& os) const
{
    os  << "        ionInducedSEE      yield = " << yield_;

    if (yield_ == "constant")
    {
        os  << ", gamma = " << gamma_;
        if (!speciesGamma_.empty())
        {
            os << " (overridden for " << speciesGamma_.toc() << ")";
        }
        os << nl;
    }
    else
    {
        os  << ", surface = " << surface_ << nl;
    }

    resolve();

    if (resolved_.empty())
    {
        os  << "            NO POSITIVE ION HAS A WALL FLUX on this patch, so"
               " this mechanism emits" << nl
            << "            NOTHING. Ions must be `transportModel"
               " driftDiffusion` to reach a wall." << nl;
        return;
    }

    os  << "            resolved per ion:";
    forAllConstIters(resolved_, it)
    {
        os  << "  " << it.key() << " = " << it.val();
    }
    os  << nl;

    if (yield_ == "hagstrum")
    {
        os  << "            Hagstrum/Raizer POTENTIAL emission only -- no"
               " kinetic channel, which is" << nl
            << "            negligible below ~1 keV but NOT in a low-pressure"
               " sheath." << nl;

        bool anyZero = false;
        forAllConstIters(resolved_, it)
        {
            if (it.val() <= 0) anyZero = true;
        }
        if (anyZero)
        {
            os  << "            A zero above is PHYSICAL, not a failure: Auger"
                   " neutralisation cannot lift a" << nl
                << "            second electron over the barrier for that ion"
                   " on this surface." << nl;
        }
    }
}


// ************************************************************************* //
