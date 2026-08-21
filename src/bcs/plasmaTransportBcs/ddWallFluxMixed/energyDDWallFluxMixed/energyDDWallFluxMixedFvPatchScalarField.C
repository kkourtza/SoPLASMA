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

    // The parent adds the SECONDARY-ELECTRON PARTICLE inflow to refGrad. The
    // same electrons carry secondaryElectronEnergy_ each into the energy
    // field, so the contribution here is that many times larger.
    //
    // Captured as a DELTA rather than recomputed, so the ion loop, the
    // per-species emission coefficients and the diffusivity lookup are
    // inherited from the electron condition rather than duplicated -- two
    // copies of that loop is exactly how the electron and energy conditions
    // would drift apart under a later change.
    const scalarField grad0(this->refGrad());

    electronDDWallFluxMixedFvPatchScalarField::updateCoeffs();

    this->refGrad() =
        grad0 + secondaryElectronEnergy_*(this->refGrad() - grad0);
}


void energyDDWallFluxMixedFvPatchScalarField::write(Ostream& os) const
{
    electronDDWallFluxMixedFvPatchScalarField::write(os);

    os.writeEntry("fluxEnergyFactor", fluxEnergyFactor_);
    os.writeEntry("secondaryElectronEnergy", secondaryElectronEnergy_);
}


} // End namespace Foam

// ************************************************************************* //
