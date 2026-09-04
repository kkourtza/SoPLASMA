/*---------------------------------------------------------------------------*\
  File: fieldEmission.C
  Part of: SoPLASMA
\*---------------------------------------------------------------------------*/

#include "fieldEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "plasmaConstants.H"

// * * * * * * * * * * * * * * * * Static Data  * * * * * * * * * * * * * * //

namespace Foam
{
namespace emissionModels
{
    defineTypeNameAndDebug(fieldEmission, 0);
    addToRunTimeSelectionTable(emissionModel, fieldEmission, dictionary);
}
}

// EXACT combinations of fundamental constants, not fits. Written out with their
// definitions so they can be checked rather than trusted.
//   a = e^3/(8 pi h)                = 1.541434e-6  A eV V^-2
//   b = (8 pi/3) sqrt(2 m_e)/(e h)  = 6.830890e9   eV^-3/2 V m^-1
// Forbes & Deane, Proc. R. Soc. A 463 (2007) 2907, Table 1.
static const Foam::scalar FN_a = 1.541434e-6;
static const Foam::scalar FN_b = 6.830890e9;

// e^3/(4 pi eps0) in eV^2 m/V, for the Schottky-Nordheim scaled barrier field
//   f = (e^3/4 pi eps0) F / phi^2
// = 1.439964 eV nm  ->  1.439964e-9 eV^2 m/V
static const Foam::scalar SN_c = 1.439964e-9;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emissionModels::fieldEmission::fieldEmission
(
    const fvPatch& p,
    const dictionary& dict
)
:
    emissionModel(p, dict),
    beta_(dict.getOrDefault<scalar>("fieldEnhancement", 1.0)),
    barrier_(dict.getOrDefault<word>("barrier", "schottkyNordheim")),
    phi_(workFunction())
{
    if (beta_ < 1.0)
    {
        FatalIOErrorInFunction(dict)
            << "`fieldEnhancement` = " << beta_ << " on patch " << p.name()
            << " is less than 1." << nl
            << "    Enhancement is a geometric amplification of the local"
               " field by surface topography;" << nl
            << "    it cannot reduce it. Use 1 for an ideally smooth surface."
            << nl << exit(FatalIOError);
    }

    if (barrier_ != "elementary" && barrier_ != "schottkyNordheim")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown `barrier` `" << barrier_ << "`." << nl
            << "    Valid: elementary | schottkyNordheim (default)." << nl
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::emissionModels::fieldEmission::emittedFlux() const
{
    const scalarField E(patchField());

    auto tflux = tmp<scalarField>::New(patch_.size(), Zero);
    scalarField& flux = tflux.ref();

    const scalar phi32 = Foam::pow(phi_, 1.5);

    forAll(flux, i)
    {
        const scalar F = beta_*E[i];

        // Below any sensible field the exponential underflows anyway; the
        // guard keeps the division and the log finite.
        if (F <= SMALL) continue;

        scalar v = 1.0;

        if (barrier_ == "schottkyNordheim")
        {
            const scalar f = SN_c*F/sqr(phi_);

            // v ~ 1 - f + (f/6) ln f is the standard approximation and is
            // valid for f < 1. At f >= 1 the barrier is suppressed entirely --
            // there is no tunnelling problem left, the electrons are simply
            // above it -- and the approximation goes negative, which would
            // turn the exponent's sign and produce a nonsensical current. So
            // it is clamped at 0, meaning "barrier gone".
            v = (f < 1.0 && f > SMALL)
              ? max(0.0, 1.0 - f + (f/6.0)*Foam::log(f))
              : (f >= 1.0 ? 0.0 : 1.0);
        }

        const scalar expo = -v*FN_b*phi32/F;

        // J [A/m2]
        const scalar J = (FN_a/phi_)*sqr(F)*Foam::exp(expo);

        // Electron flux = current density / elementary charge.
        flux[i] = J/constant::plasma::eCharge.value();
    }

    return tflux;
}


void Foam::emissionModels::fieldEmission::report(Ostream& os) const
{
    os  << "        fieldEmission      phi = " << phi_ << " eV"
        << ", beta = " << beta_
        << ", barrier = " << barrier_ << nl;

    if (beta_ == 1.0)
    {
        os  << "            beta = 1 (ideally smooth). Emission is EXPONENTIAL"
               " in the field, so a real" << nl
            << "            surface's roughness usually dominates it -- fit"
               " beta if this matters." << nl;
    }
}


// ************************************************************************* //
