/*---------------------------------------------------------------------------*\
  File: energyDDWallFluxMixedFvPatchScalarField.C
  Part of: SoPLASMA
  Copyright (C) 2026
  License: GNU General Public License v3 or later
\*---------------------------------------------------------------------------*/

#include "energyDDWallFluxMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

namespace Foam
{

defineTypeNameAndDebug(energyDDWallFluxMixedFvPatchScalarField, 0);

addToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    energyDDWallFluxMixedFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField>
energyDDWallFluxMixedFvPatchScalarField::calcAbsorptionVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal
) const
{
    // The electron condition, scaled. The energy leaves on the SAME electrons,
    // so reusing its velocity keeps the two conditions consistent by
    // construction: any change to the thermal velocity, the drift-flux option
    // or the surface treatment is inherited rather than duplicated here.
    tmp<scalarField> tVel =
        electronDDWallFluxMixedFvPatchScalarField::calcAbsorptionVelocity
        (
            m, T, uDriftNormal
        );

    tVel.ref() *= fluxEnergyFactor_;

    return tVel;
}


tmp<scalarField>
energyDDWallFluxMixedFvPatchScalarField::calcEffectiveWallVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal
) const
{
    tmp<scalarField> tVel =
        electronDDWallFluxMixedFvPatchScalarField::calcEffectiveWallVelocity
        (
            m, T, uDriftNormal
        );

    tVel.ref() *= fluxEnergyFactor_;

    return tVel;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

energyDDWallFluxMixedFvPatchScalarField::
energyDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    electronDDWallFluxMixedFvPatchScalarField(p, iF),
    fluxEnergyFactor_(4.0/3.0),
    secondaryElectronEnergy_(0.0)
{}


energyDDWallFluxMixedFvPatchScalarField::
energyDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    electronDDWallFluxMixedFvPatchScalarField(p, iF, dict),
    // 4/3, not 5/2 and not 1 -- see the header. 5/2 is the INFLOW value for a
    // Maxwellian reservoir and would over-drain at every wall; 1 is the
    // population mean and would under-drain, because faster electrons reach
    // the wall more often.
    fluxEnergyFactor_(dict.lookupOrDefault<scalar>("fluxEnergyFactor", 4.0/3.0)),
    // USER-DEFINED, with a default rather than a derivation, because there is
    // no rigid calculation available: the birth energy of a secondary electron
    // depends on the surface material and work function and on the incident
    // ion energy, none of which this solver models. 2 eV is a conventional
    // value for metals and dielectrics in the 1-5 eV range reported in the
    // literature. Anything claiming to compute it from first principles here
    // would be inventing the number, not deriving it.
    secondaryElectronEnergy_
    (
        dict.lookupOrDefault<scalar>("secondaryElectronEnergy", 2.0)
    )
{}


energyDDWallFluxMixedFvPatchScalarField::
energyDDWallFluxMixedFvPatchScalarField
(
    const energyDDWallFluxMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    electronDDWallFluxMixedFvPatchScalarField(ptf, p, iF, mapper),
    fluxEnergyFactor_(ptf.fluxEnergyFactor_),
    secondaryElectronEnergy_(ptf.secondaryElectronEnergy_)
{}


energyDDWallFluxMixedFvPatchScalarField::
energyDDWallFluxMixedFvPatchScalarField
(
    const energyDDWallFluxMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    electronDDWallFluxMixedFvPatchScalarField(ptf, iF),
    fluxEnergyFactor_(ptf.fluxEnergyFactor_),
    secondaryElectronEnergy_(ptf.secondaryElectronEnergy_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void energyDDWallFluxMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // The secondary electrons that enter here are the SAME ones the electron
    // condition already counted; each carries secondaryElectronEnergy_, so the
    // energy inflow is simply that inflow times that energy.
    //
    // WAS A DELTA, AND THE DELTA WAS WRONG. The previous version captured
    // refGrad BEFORE calling the parent and treated it as a baseline:
    //
    //     grad0 = refGrad;  parent();  refGrad = grad0 + s*(refGrad - grad0);
    //
    // That assumes the parent ACCUMULATES onto refGrad. It does not -- the
    // base `ddWallFluxMixed::updateCoeffs` RESETS `refGrad() = 0.0` before
    // adding the emission term. So grad0 was never a baseline; it was the
    // value left over from the PREVIOUS corrector, which the parent had
    // already discarded. The recursion that produced,
    //
    //     x_{n+1} = s*SEE_{n+1} + (1 - s)*x_n,
    //
    // does not settle on s*SEE at all: at the shipped s = 2 it alternates
    // sign every corrector and never converges. Nothing downstream can damp a
    // boundary condition that flips sign on each pass.
    //
    // After the parent returns, refGrad IS the secondary-electron particle
    // inflow (the base zeroed it first), so the energy version is one
    // multiplication. The ion loop, the per-species coefficients and the
    // diffusivity lookup are still inherited, not duplicated.
    // The parent has now built the wall flux from the ENERGY's own mu_eps and
    // D_eps (see patchMobility/patchDiffusivity below), and its refGrad is the
    // secondary-electron inflow expressed against D_eps. Those are the same
    // electrons the electron condition emitted; each carries
    // secondaryElectronEnergy_, so the energy inflow is that inflow times that
    // energy -- one multiplication, no diffusivity ratio anywhere.
    electronDDWallFluxMixedFvPatchScalarField::updateCoeffs();

    this->refGrad() *= secondaryElectronEnergy_;
}


tmp<scalarField> energyDDWallFluxMixedFvPatchScalarField::patchMobility
(
    const driftDiffusion& ddModel
) const
{
    return energyCoefficient("muEps", ddModel, true);
}


tmp<scalarField> energyDDWallFluxMixedFvPatchScalarField::patchDiffusivity
(
    const driftDiffusion& ddModel
) const
{
    return energyCoefficient("DEps", ddModel, false);
}


//- The energy's own coefficient, published by localEnergyEnergyModel from the
//  same expression its laplacian/div use. NOT derived from the electron
//  coefficient by a factor: the 5/3 is the Maxwellian limit, while these are
//  tabulated from the EEDF (the model applies a factor of 1 in that case).
tmp<scalarField> energyDDWallFluxMixedFvPatchScalarField::energyCoefficient
(
    const word& fieldName,
    const driftDiffusion& ddModel,
    const bool isMobility
) const
{
    if (!this->db().foundObject<volScalarField>(fieldName))
    {
        FatalErrorInFunction
            << "energyDDWallFluxMixed needs `" << fieldName << "`, which"
               " localEnergyEnergyModel registers." << nl
            << "    It was not found: this condition is only meaningful with"
               " an LMEA electron-energy" << nl
            << "    equation enabled." << nl
            << exit(FatalError);
    }

    return tmp<scalarField>::New
    (
        this->db().lookupObject<volScalarField>(fieldName)
            .boundaryField()[this->patch().index()]
    );
}


void energyDDWallFluxMixedFvPatchScalarField::write(Ostream& os) const
{
    electronDDWallFluxMixedFvPatchScalarField::write(os);

    os.writeEntry("fluxEnergyFactor", fluxEnergyFactor_);
    os.writeEntry("secondaryElectronEnergy", secondaryElectronEnergy_);
}


} // End namespace Foam

// ************************************************************************* //
