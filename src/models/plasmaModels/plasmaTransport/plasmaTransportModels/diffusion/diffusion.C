/*---------------------------------------------------------------------------*\
  File: diffusion.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::diffusion.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "diffusion.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(diffusion, 0);
addToRunTimeSelectionTable
(
    plasmaTransportModel,
    diffusion,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

diffusion::diffusion
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
    ),
    diffusivityModel_(nullptr)
{
    const word& sName = species_.speciesName(specieIndex_);

    // UNCHARGED ONLY, because this model has no drift term and a charged
    // species drifts. Selecting it for one would silently delete the dominant
    // transport mechanism -- the field would still accelerate the species in
    // the Poisson solve while nothing moved it.
    if (mag(species_.speciesChargeNumber(specieIndex_)) > SMALL)
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << sName << "' is CHARGED (charge number "
            << species_.speciesChargeNumber(specieIndex_)
            << ") and cannot use `transportModel diffusion`." << nl << nl
            << "    This model carries no drift term, and a charged species"
            << " drifts." << nl << nl
            << "    Use `transportModel driftDiffusion`, which offers both"
            << " flux treatments:" << nl
            << "      fluxScheme standard           drift and diffusion as"
            << " separate terms, so the" << nl
            << "                                    convective scheme is yours"
            << " in fvSchemes (ROUNDF, TVD, ...)" << nl
            << "      fluxScheme ScharfetterGummel  the two combined in one"
            << " exponentially fitted flux," << nl
            << "                                    stable where drift"
            << " dominates" << nl << nl
            << "    Or `immobile` if the species does not move on the"
            << " timescale of interest." << nl
            << exit(FatalIOError);
    }

    if (!dict_.found("diffusivity"))
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << sName
            << "': `transportModel diffusion` needs a 'diffusivity'"
            << " sub-dictionary in '" << dict_.name() << "'." << nl
            << exit(FatalIOError);
    }

    const dictionary& diffusivityDict = dict_.subDict("diffusivity");

    diffusivityModel_ = plasmaDiffusivityModel::New
    (
        diffusivityDict.get<word>("type"),
        diffusivityDict,
        mesh_,
        species_,
        specieIndex_
    );
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void diffusion::correct()
{
    diffusivityModel_->correct();
}


tmp<fvScalarMatrix> diffusion::nEqn() const
{
    volScalarField& n = species_.numberDensity(specieIndex_);
    const volScalarField& D = diffusivityModel_->D();

    tmp<fvScalarMatrix> tEqn = fvm::ddt(n);
    tEqn.ref() -= fvm::laplacian(D, n);

    return tEqn;
}


void diffusion::updateFluxes
(
    const fvScalarMatrix& nEqnMatrix,
    surfaceScalarField& convectiveFlux,
    surfaceScalarField& diffusiveFlux,
    surfaceScalarField& particleFlux
) const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);
    const volScalarField& D = diffusivityModel_->D();

    // Sign convention follows driftDiffusion's `standard` branch, where the
    // Laplacian is subtracted from the equation and its flux negated.
    const fvScalarMatrix diffMat(fvm::laplacian(D, n));

    diffusiveFlux = -diffMat.flux();
    convectiveFlux = 0.0*diffusiveFlux;   // no drift: identically zero
    particleFlux = diffusiveFlux;
}


tmp<volScalarField> diffusion::electricalConductivity() const
{
    // Zero rather than fatal: plasmaTransport sums this over the MOBILE
    // species, which now includes uncharged diffusing ones. An uncharged
    // species contributes no conductivity, so zero is both the safe answer for
    // the caller and the physically correct one.
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "sigma_" + species_.speciesName(specieIndex_),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimensionSet(-1, -3, 3, 0, 0, 2, 0), Zero)
    );
}


tmp<volScalarField> diffusion::diffusiveChargeSource() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "diffChargeSrc_" + species_.speciesName(specieIndex_),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimensionSet(0, -3, 0, 0, 0, 1, 0), Zero)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
