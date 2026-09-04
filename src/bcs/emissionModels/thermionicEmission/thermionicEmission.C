/*---------------------------------------------------------------------------*\
  File: thermionicEmission.C
  Part of: SoPLASMA
\*---------------------------------------------------------------------------*/

#include "thermionicEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "plasmaConstants.H"

namespace Foam
{
namespace emissionModels
{
    defineTypeNameAndDebug(thermionicEmission, 0);
    addToRunTimeSelectionTable(emissionModel, thermionicEmission, dictionary);
}
}

// A0 = 4 pi m_e e k_B^2 / h^3, the universal Richardson constant.
// = 1.201735e6 A m^-2 K^-2   (120.1735 A cm^-2 K^-2)
// Exact in the fundamental constants; CRC Handbook 95th ed., sec. 12.
static const Foam::scalar RD_A0 = 1.201735e6;

// sqrt(e^3/(4 pi eps0)) expressed for phi in eV and F in V/m:
//   dphi[eV] = 3.794686e-5 sqrt(F[V/m])
static const Foam::scalar SCHOTTKY_C = 3.794686e-5;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emissionModels::thermionicEmission::thermionicEmission
(
    const fvPatch& p,
    const dictionary& dict
)
:
    emissionModel(p, dict),
    lambdaR_(dict.getOrDefault<scalar>("richardsonCorrection", 1.0)),
    schottky_(dict.getOrDefault<Switch>("schottky", true)),
    phi_(workFunction())
{
    if (lambdaR_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "`richardsonCorrection` must be positive, not " << lambdaR_
            << " (patch " << p.name() << ")." << nl << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::emissionModels::thermionicEmission::emittedFlux() const
{
    const scalarField T(patchTemperature());

    // Only looked up when actually needed, so a case with `schottky false`
    // does not require a field it never uses.
    tmp<scalarField> tE
    (
        schottky_ ? patchField() : tmp<scalarField>::New(patch_.size(), Zero)
    );
    const scalarField& E = tE();

    auto tflux = tmp<scalarField>::New(patch_.size(), Zero);
    scalarField& flux = tflux.ref();

    // k_B in eV/K, so phi and k_B T are in the same unit and the exponent is
    // dimensionless without a further conversion.
    const scalar kB_eV =
        constant::plasma::kappaBoltzmann.value()
      / constant::plasma::eCharge.value();

    forAll(flux, i)
    {
        if (T[i] <= SMALL) continue;

        scalar phiEff = phi_;

        if (schottky_)
        {
            // The barrier cannot be lowered below zero: at that point there is
            // no barrier and the emission is not thermionic at all. Clamped,
            // and the clamp is physical rather than numerical.
            phiEff = max(0.0, phi_ - SCHOTTKY_C*Foam::sqrt(max(E[i], 0.0)));
        }

        const scalar J = lambdaR_*RD_A0*sqr(T[i])
                       * Foam::exp(-phiEff/(kB_eV*T[i]));

        flux[i] = J/constant::plasma::eCharge.value();
    }

    return tflux;
}


void Foam::emissionModels::thermionicEmission::report(Ostream& os) const
{
    os  << "        thermionicEmission phi = " << phi_ << " eV"
        << ", lambda_R = " << lambdaR_
        << ", schottky = " << schottky_.c_str() << nl
        << "            temperature is the BOUNDARY VALUE of T_gas on this"
           " patch -- set that condition" << nl
        << "            to run a heated electrode." << nl;

    if (lambdaR_ == 1.0)
    {
        os  << "            lambda_R = 1 is the free-electron value; it is ~0.5"
               " for tungsten, so this" << nl
            << "            overestimates by ~2x -- small against the"
               " exponential, but not zero." << nl;
    }
}


// ************************************************************************* //
