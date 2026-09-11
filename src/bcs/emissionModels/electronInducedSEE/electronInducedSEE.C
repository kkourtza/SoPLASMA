/*---------------------------------------------------------------------------*\
  File: electronInducedSEE.C
  Part of: SoPLASMA
\*---------------------------------------------------------------------------*/

#include "electronInducedSEE.H"
#include "addToRunTimeSelectionTable.H"
#include "materialLibrary.H"
#include "surfaceFields.H"
#include "volFields.H"

namespace Foam
{
namespace emissionModels
{
    defineTypeNameAndDebug(electronInducedSEE, 0);
    addToRunTimeSelectionTable(emissionModel, electronInducedSEE, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emissionModels::electronInducedSEE::electronInducedSEE
(
    const fvPatch& p,
    const dictionary& dict
)
:
    emissionModel(p, dict),
    deltaMax_(required("deltaMax", "deltaMax")),
    energyDeltaMax_(required("energyDeltaMax", "energyDeltaMax")),
    // Vaughan's threshold. Overridable because it is the one number in the
    // curve that is a convention rather than a measurement.
    E0_(dict.getOrDefault<scalar>("thresholdEnergy", 12.5)),
    caveatReported_(false)
{
    if (deltaMax_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "`deltaMax` must be positive, got " << deltaMax_ << "." << nl
            << exit(FatalIOError);
    }

    if (energyDeltaMax_ <= E0_)
    {
        FatalIOErrorInFunction(dict)
            << "`energyDeltaMax` (" << energyDeltaMax_
            << " eV) must exceed the threshold (" << E0_ << " eV)." << nl
            << "    Vaughan's curve is parametrised by v = (E-E0)/(Emax-E0),"
               " which is meaningless" << nl
            << "    otherwise." << nl << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::emissionModels::electronInducedSEE::required
(
    const word& key,
    const word& libKey
) const
{
    if (dict_.found(key)) return dict_.get<scalar>(key);

    if (dict_.found("material"))
    {
        const word mat(dict_.get<word>("material"));

        const string ctx
        (
            "electronInducedSEE on patch `" + patch_.name() + "`"
        );

        // Look in BOTH libraries: a surface may be a conductor or a
        // dielectric, and the yield parameters live with whichever describes
        // it. Tried in that order; the error below lists both.
        for (const word& lib : {word("workFunctions"), word("dielectrics")})
        {
            if (materialLibrary::available(lib))
            {
                const dictionary& e =
                    materialLibrary::lookup(mat, ctx, lib);

                if (e.found(libKey)) return e.get<scalar>(libKey);
            }
        }
    }

    FatalIOErrorInFunction(dict_)
        << "electronInducedSEE on patch `" << patch_.name() << "` needs `"
        << key << "`." << nl << nl
        << "    Give it explicitly, or name a `material` whose library entry"
           " carries it." << nl << nl
        << "    IT IS NOT DEFAULTED, and for an INSULATOR it is deliberately"
           " not shipped: the" << nl
        << "    published spread for one insulator (MgO spans roughly 4 to 25)"
           " is far larger than" << nl
        << "    the difference between materials, so a point value would be"
           " authoritative-looking" << nl
        << "    and wrong. Measure it, or take a value from a source you can"
           " cite for YOUR surface." << nl
        << exit(FatalIOError);

    return 0;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::emissionModels::electronInducedSEE::vaughanYield
(
    const scalar E,
    const scalar deltaMax,
    const scalar Emax,
    const scalar E0
)
{
    // Below the threshold nothing is emitted. This is the branch an
    // atmospheric-pressure discharge lives in: electrons arrive with a few eV
    // against a 12.5 eV threshold.
    if (E <= E0) return 0.0;

    const scalar v = (E - E0)/(Emax - E0);

    if (v <= 0) return 0.0;

    // Vaughan, IEEE Trans. Electron Devices 36 (1989) 1963; 40 (1993) 830.
    if (v > 3.6)
    {
        return deltaMax*1.125*Foam::pow(v, -0.35);
    }

    const scalar k = (v < 1.0) ? 0.56 : 0.25;

    return deltaMax*Foam::pow(v*Foam::exp(1.0 - v), k);
}


Foam::tmp<Foam::scalarField>
Foam::emissionModels::electronInducedSEE::emittedFlux() const
{
    auto tflux = tmp<scalarField>::New(patch_.size(), Zero);
    scalarField& flux = tflux.ref();

    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    // The incident electron flux, from the same per-species flux field the
    // ion mechanism uses.
    if (!mesh.foundObject<surfaceScalarField>("particleFlux_e")) return tflux;

    // The incident ENERGY. `meanE` exists under LMEA; without it there is no
    // energy to evaluate a strongly non-linear yield at, and assuming one
    // would invent the answer.
    if (!mesh.foundObject<volScalarField>("meanE"))
    {
        FatalErrorInFunction
            << "electronInducedSEE on patch `" << patch_.name()
            << "` needs the electron mean energy `meanE`, which does not"
               " exist." << nl << nl
            << "    It exists under `electronEnergyModel LMEA`. Under LFA there"
               " is no electron energy" << nl
            << "    field to evaluate the yield at, and Vaughan's curve is"
               " strongly non-linear in it," << nl
            << "    so no reasonable single value can be assumed." << nl
            << exit(FatalError);
    }

    const fvsPatchScalarField& phiE =
        patch_.lookupPatchField<surfaceScalarField, scalar>("particleFlux_e");

    const scalarField& meanE =
        mesh.lookupObject<volScalarField>("meanE")
            .boundaryField()[patch_.index()];

    const scalarField& magSf = patch_.magSf();

    forAll(flux, i)
    {
        const scalar incident = max(scalar(0), phiE[i]/magSf[i]);

        if (incident <= 0) continue;

        flux[i] =
            vaughanYield(meanE[i], deltaMax_, energyDeltaMax_, E0_)*incident;
    }

    return tflux;
}


void Foam::emissionModels::electronInducedSEE::report(Ostream& os) const
{
    os  << "        electronInducedSEE delta_max = " << deltaMax_
        << " at " << energyDeltaMax_ << " eV, threshold " << E0_ << " eV"
        << nl
        << "            Vaughan's universal curve (1989, corr. 1993)." << nl
        << "            EVALUATED AT THE MEAN ENERGY, which is an"
           " APPROXIMATION: delta is strongly" << nl
        << "            non-linear, so the correct treatment integrates"
           " delta(E) over the electron" << nl
        << "            energy DISTRIBUTION at the wall. The two agree only"
           " where delta is locally" << nl
        << "            linear." << nl;

    if (deltaMax_ > 1.0)
    {
        os  << "            delta_max > 1: this surface can emit MORE"
               " electrons than it collects," << nl
            << "            which is how a dielectric charges POSITIVE. Expect"
               " that, and check it." << nl;
    }
}


// ************************************************************************* //
