/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.
    Copyright (C) 2026
    GNU General Public License v3 or later.

Class
    Foam::energyDDWallFluxMixedFvPatchScalarField

Description
    Wall boundary condition for the ELECTRON ENERGY DENSITY n_eps under LMEA.

    Follows the electron condition exactly -- thermal absorption velocity,
    optional drift flux, secondary emission, surface charging -- because the
    energy leaves the domain on the same electrons the electron condition
    already accounts for. Only the WEIGHT differs:

        Q_e . n = 2 k_B T_e Gamma_e^{out} - dE_se Gamma_se^{in}

    (Esgee VizGlow theory manual Sect. 2.4.)

    THE FACTOR IS 4/3, AND THAT IS THE ONLY SUBTLE PART. For a Maxwellian the
    FLUX-weighted mean energy is 2 k_B T_e, while the POPULATION mean energy is
    (3/2) k_B T_e = eps. So

        2 k_B T_e = (4/3) eps    =>    Gamma_eps = (4/3) eps Gamma_e
                                                = (4/3) n_eps u_abs

    which is the electron condition applied to n_eps with the velocities scaled
    by 4/3. It is NOT 5/2: that is the value for electrons ENTERING from a
    Maxwellian reservoir, and using it here would over-drain the energy at
    every wall. It is not 1 either -- taking the population mean would
    UNDER-drain, because faster electrons reach the wall more often.

    The secondary-electron inflow enters with the emission energy rather than
    with the local mean energy, since those electrons are born at the surface.

SourceFiles
    energyDDWallFluxMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef energyDDWallFluxMixedFvPatchScalarField_H
#define energyDDWallFluxMixedFvPatchScalarField_H

#include "electronDDWallFluxMixedFvPatchScalarField.H"

namespace Foam
{

class energyDDWallFluxMixedFvPatchScalarField
:
    public electronDDWallFluxMixedFvPatchScalarField
{
protected:

    //- Flux-weighted energy factor, 2k_B T_e / eps = 4/3 for a Maxwellian.
    //  Exposed because a strongly non-Maxwellian EEDF has a different ratio,
    //  and the two-term solver knows it even though this condition does not.
    scalar fluxEnergyFactor_;

    //- Energy of a secondary electron at birth [eV]. USER-DEFINED with a
    //  2 eV default: it depends on the surface material, its work function
    //  and the incident ion energy, none of which this solver models, so
    //  there is no rigid calculation to make automatically.
    scalar secondaryElectronEnergy_;

    //- The ENERGY's transport coefficients, not the electron's.

    virtual tmp<scalarField> patchMobility(const driftDiffusion&) const;

    virtual tmp<scalarField> patchDiffusivity(const driftDiffusion&) const;


    //- Fetch a coefficient published by localEnergyEnergyModel.

    tmp<scalarField> energyCoefficient

    (

        const word& fieldName,

        const driftDiffusion& ddModel,

        const bool isMobility

    ) const;


    virtual tmp<scalarField> calcAbsorptionVelocity
    (
        const dimensionedScalar& m,
        const scalarField& T,
        const scalarField& uDriftNormal
    ) const override;

    virtual tmp<scalarField> calcEffectiveWallVelocity
    (
        const dimensionedScalar& m,
        const scalarField& T,
        const scalarField& uDriftNormal
    ) const override;

public:

    TypeName("energyDDWallFluxMixed");

    energyDDWallFluxMixedFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    energyDDWallFluxMixedFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    energyDDWallFluxMixedFvPatchScalarField
    (
        const energyDDWallFluxMixedFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    energyDDWallFluxMixedFvPatchScalarField
    (
        const energyDDWallFluxMixedFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchField<scalar>> clone() const
    {
        return tmp<fvPatchField<scalar>>
        (
            new energyDDWallFluxMixedFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchField<scalar>> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchField<scalar>>
        (
            new energyDDWallFluxMixedFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
