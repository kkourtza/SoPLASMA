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

        Q . n = (1 - r) n_e w_w eps_w - Gc eps_s        [eq. (6.13) split]

    (Hagelaar HDR, Literature/, ch. 6 pp. 38-39.) The loss speed w_w and the
    reflection probability r are the ELECTRON condition's, unchanged -- the
    energy leaves the domain on the same electrons. Only the WEIGHT differs.

    THE WEIGHT IS 5/3, AND THAT IS THE ONLY SUBTLE PART. Eq. (6.15) averages
    the energy flux over a SHIFTED Maxwellian -- the same distribution
    assumption that produces the particle base (6.6) -- and gives the mean
    energy lost per electron as

        eps_w = ( 5/2 - sqrt(pi)*Gamma_w/(n_e vT) ) T_e
              = ( 5/2 - (Gamma_w/n_e)/A ) T_e            with A = vT/sqrt(pi)

    so with eps = (3/2) T_e for a Maxwellian population,

        eps_w/eps = 5/3 - (2/3)*(Gamma_w/n_e)/A .

    SUPERSEDED 2026-09-04: this condition previously used the CONSTANT 4/3,
    i.e. eps_w = 2 k_B T_e, which is eq. (6.14) -- the CENTRED half-Maxwellian
    that pairs with the particle base (6.3). Commit 7e6a4b3 moved the particle
    base to (6.6) and left the energy weight on (6.14), so the pair was
    mismatched; the consistent pairs are (6.3, 6.14) and (6.6, 6.15). What
    still holds from the old note: the weight is NOT 1 (the population mean
    would under-drain, because faster electrons reach the wall more often).
    What does not hold is the old file's argument that "it is NOT 5/2, that is
    the value for electrons entering from a Maxwellian reservoir" -- that
    conflated Hagelaar's (5/2)T_e with the (5/2)kT monatomic enthalpy: same
    number, different origin.

    TWO INDEPENDENT CHECKS, both from the text and both reproduced exactly on
    2026-09-04:
      - at Gamma_w = 0, (6.15) is a factor 5/4 above (6.14): 5/3 = (5/4)(4/3);
      - substituting the reflection-only Gamma_w = r n w_w and letting r -> 1
        gives eps_w -> 2 T_e, "as for a centred Maxwellian", which the text
        asserts and which the formula above satisfies identically.
    A third, structural: the HDR's own energy equation carries W = (5/3) n e w
    (p. 38), so a BC weight of 5/3 matches the bulk energy-flux definition
    while 4/3 does not.

    REFLECTED ELECTRONS CARRY THEIR ENERGY BACK. Eq. (6.13) as printed puts a
    single eps_s on the whole of Gamma_w; splitting Gamma_w by (6.2) into the
    reflected part r n w_w -- which returns at the INCIDENT mean energy eps_w,
    reflection being elastic -- and the created part Gc, which is born at the
    surface with eps_s, gives the form at the top. That split is the faithful
    reading of (6.2), and it makes the energy closure exactly parallel to the
    particle one.

    Both effects enter through the SAME gRatio = (Gamma_w/n_e)/A that the
    electron condition already computes, so there is one definition of
    Gamma_w in the pair and not two.

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

    //- Flux-weighted energy weight eps_w/eps at ZERO wall creation and zero
    //  reflection: 5/3 for a Maxwellian, from eq. (6.15). The Gamma_w
    //  correction of (6.15) is applied on top of this, per face.
    //
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


    //- ALWAYS the electron. This condition exists only for the


    //  electron-energy field, so the species is not something the


    //  user should have to state -- and the base class's inference


    //  (strip a leading "n_") gets it WRONG for `nEps_e`, failing


    //  with "nEps_e not found in table" and a list of species.


    //  An explicit `species` entry still wins, for a case that


    //  names its electron something unusual.


    virtual word resolveSpeciesName() const;



    //- n_e on this patch, NOT n_eps.
    //
    //  MEASURED DEFECT, 2026-09-04: this class inherited the electron
    //  routine's Gamma_w/n with `*this` as the density, so on the energy field
    //  it divided the emitted ELECTRON flux by n_eps = n_e*eps[eV]. Since
    //  nEps_e carries dimensions [0 -3 0 0 0 0 0] the result is still a
    //  velocity, so DIMENSIONAL ANALYSIS COULD NOT CATCH IT -- the eq. (6.8)
    //  creation correction was simply under-applied by a factor of eps in eV,
    //  i.e. by about 1-10 in a discharge. Dated because it is a live fix, not
    //  a historical note.
    virtual tmp<scalarField> normalisingDensity() const override;

    //- eps_w/eps per face, eq. (6.15), clamped at zero.
    //
    //  Hagelaar does not clamp (6.15); the clamp is added here for the same
    //  reason (6.6) carries one -- under strong creation the |w| << vT
    //  expansion behind (6.15) no longer holds, and a NEGATIVE mean energy
    //  per lost electron is not a physical state to impose.
    tmp<scalarField> energyWeight(const scalarField& gRatio) const;

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
