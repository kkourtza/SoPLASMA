/*---------------------------------------------------------------------------*\
  File: immobile.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::immobile.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "fvm.H"

#include "immobile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(immobile, 0);
addToRunTimeSelectionTable
(
    plasmaTransportModel,
    immobile, 
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immobile::immobile
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    plasmaSpecies& species,
    const label specieIndex
)
:
    plasmaTransportModel
    (
        modelName, 
        dict, 
        mesh, 
        species, 
        specieIndex
    )
{}

namespace
{
    //- Stand-in for "no transport" in a coefficient that will be HARMONICALLY
    //  interpolated, where exactly zero is a division by zero. See
    //  immobile::mu() for the measurement and the magnitude argument.
    static const Foam::scalar immobileTinyCoeff = 1.0e-30;
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void immobile::correct()
{
    // No tranport coefficients to update
}

tmp<fvScalarMatrix> immobile::nEqn() const
{
    return fvm::ddt(species_.numberDensity(specieIndex_));
}

void immobile::updateFluxes
(
    const fvScalarMatrix& nEqnMatrix,
    surfaceScalarField&,
    surfaceScalarField&,
    surfaceScalarField&
) const
{
    FatalErrorInFunction
        << "updateFluxes() called for immobile species '"
        << species_.speciesName(specieIndex_) << "'." << nl
        << "Immobile species have no particle flux." << nl
        << "This is a logic error in plasmaTransport." << nl
        << abort(FatalError);
}

tmp<volScalarField> immobile::mu() const
{
    // NEGLIGIBLE, NOT EXACTLY ZERO -- and the difference is load-bearing.
    //
    // Exactly zero SIGFPEs. The streamer beds interpolate mobility with
    // `harmonic` (2ab/(a+b)) and diffusivity with `Gauss harmonic corrected`,
    // so an identically-zero coefficient makes the face value 0/0. Caught
    // under gdb 2026-09-10:
    //     #0 Foam::divide(Field&, const double&, const UList&)
    //     #2 Foam::fvc::interpolate<double>(volScalarField, word)
    //     #4 residualCallback
    // The case's own fvSchemes predicted it: its "interpolate\(mu_.*\)"
    // catch-all is annotated "Reached only when an ion is given a
    // transportModel other than `immobile`" -- which is exactly what handing
    // an immobile species to the Newton assembly does.
    //
    // 1e-30 is harmonic-safe (harmonic(e,e) = e, no division by zero, and
    // 1e-30 is a normal double, nowhere near the 2.2e-308 subnormal floor)
    // and physically inert: against an ion mobility of ~1e-4 m^2/(V.s) it is
    // 26 orders down, giving a drift flux ~1e-4 m^-2 s^-1 where the chemistry
    // sources are ~1e26 -- 30 orders below the terms it sits beside.
    //
    // PICARD IS UNAFFECTED EITHER WAY: immobile::nEqn() is fvm::ddt(n) alone,
    // so mu()/D() are reached only by a caller that assembles the species
    // equation itself, i.e. the Newton outer solver.
    return volScalarField::New
    (
        "mu_" + species_.speciesName(specieIndex_),
        mesh(),
        dimensionedScalar(dimensionSet(-1, 0, 2, 0, 0, 1, 0), immobileTinyCoeff)
    );
}


tmp<volScalarField> immobile::D() const
{
    // Negligible, not zero -- see mu() for the harmonic-interpolation reason.
    return volScalarField::New
    (
        "D_" + species_.speciesName(specieIndex_),
        mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0, 0, 0), immobileTinyCoeff)
    );
}


tmp<volScalarField> immobile::electricalConductivity() const
{
    FatalErrorInFunction
        << "electricalConductivity() called for immobile species '"
        << species_.speciesName(specieIndex_) << "'." << nl
        << "Immobile species have zero conductivity and should be "
        << "skipped by the caller (use mobileSpeciesIDs/chargedSpeciesIDs)." << nl
        << abort(FatalError);
    return nullptr; 
}

tmp<volScalarField> immobile::diffusiveChargeSource() const
{
    FatalErrorInFunction
        << "diffusiveChargeSource() called for immobile species '"
        << species_.speciesName(specieIndex_) << "'." << nl
        << "Immobile species have zero diffusive charge source and should be "
        << "skipped by the caller (use mobileSpeciesIDs/chargedSpeciesIDs)." << nl
        << abort(FatalError);
    return nullptr; 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
