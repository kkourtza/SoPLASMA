/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026 Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::electronDDWallFluxMixedFvPatchScalarField

Description
    Concrete implementation of ddWallFluxMixedFvPatchScalarField for 
    electron species.

    This boundary condition calculates the electron flux to a solid surface 
    based on thermal motion. Given a switch (includeDriftFlux_) the user can
    include the drift flux for the calculation of the total wall flux (if the
    drift is towards the wall).

    The boundary condition accounts for the Secondary Electron Emission, based
    on the ion wall flux.

Usage
    Example boundary condition entry:
    \verbatim
    <patchName>
    {
        type            electronDDWallFluxMixed;

        // Electron temperature [K]
        T               Te; 

        // Optional: Enable surface charge accumulation
        enableSurfaceCharging false;

        // Optional: Include the drift flux
        includeDriftFlux      false;

        // Optional: Enable secondary electron emission
        enableSEE             true;
        defaultSEEC           0.05; // SEE coefficient
        speciesSEEC
        {
            ion1              0.01;
            ion2              0.04;
        }
    }
    \endverbatim
    
SourceFiles
    electronDDWallFluxMixedFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef electronDDWallFluxMixedFvPatchScalarField_H
#define electronDDWallFluxMixedFvPatchScalarField_H

#include "ddWallFluxMixedFvPatchScalarField.H"
#include "emissionModel.H"
#include "PtrList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class electronDDWallFluxMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class electronDDWallFluxMixedFvPatchScalarField
:
    public ddWallFluxMixedFvPatchScalarField
{
protected:

    // Protected data

       //- Switch to enable surface charging
        Switch enableSurfaceCharging_;

        //- Switch to include drift velocity in wall flux calculation
        Switch includeDriftFlux_;

        //- Switch to enable secondary electron emission
        Switch enableSEE_;

        //- Diagnostic: print, per patch, the ion flux that drives secondary
        //  emission and the emission it produces. Emission is expected at the
        //  CATHODE (positive ions are accelerated into it) and NOT at the
        //  anode (they are repelled), so a non-zero anode total is a finding,
        //  not a detail. Default off.
        Switch seeReport_;

        //- Latch so the "SEE cannot fire" warning is printed once, not once
        //  per corrector per patch.
        mutable bool seeInertReported_;

        //- 1 eV when the case names no temperature field. Electron walls are
        //  not at the gas temperature, and an LFA case running the
        //  `gasTemperature` energy model registers no T_e to follow.
        virtual dimensionedScalar defaultTValue() const;

        //- Default SEE coefficient
        scalar defaultSEEC_;

        //- Species-specific SEE coefficients
        dictionary speciesSEEC_;

        //- ADDITIONAL EMISSION MECHANISMS, as a sub-dictionary keyed by model:
        //
        //      emission
        //      {
        //          ionInducedSEE { yield hagstrum; surface metal; material copper; }
        //          fieldEmission { material copper; fieldEnhancement 50; }
        //          thermionicEmission { material tungsten; }
        //      }
        //
        //  Each contributes an electron flux, summed and added to refGradient
        //  as Gamma/D -- the same way the built-in `enableSEE` path does.
        //
        //  STORED AS A DICTIONARY and built on first use, not constructed
        //  here: the models are abstract, so a PtrList of them would make this
        //  condition's four copy constructors do deep clones of a runtime-
        //  selected hierarchy for no benefit. Lazy construction also lets
        //  ionInducedSEE reach `plasmaTransport` in the registry, which does
        //  not exist while boundary conditions are being built.
        dictionary emissionDict_;

        //- The models, built on first updateCoeffs()
        mutable PtrList<emissionModel> emission_;

        //- Whether the emission report has been printed
        mutable bool emissionReported_;

        //- List to keep the SEE coefficients for each ion species
        List<scalar> seec_;

        //- Switch to read the SEE entries only once
        Switch mapped_;

    // Protected Member Functions

        //- Calculate the electron-specific wall absorption velocity
        virtual tmp<scalarField> calcAbsorptionVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal
        ) const override;

        //- Calculate the electron-specific wall effective velocity
        virtual tmp<scalarField> calcEffectiveWallVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal
        ) const override;

public:

    //- Runtime type information
    TypeName("electronDDWallFluxMixed");

    // Constructors

        //- Construct from patch and internal field
        electronDDWallFluxMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        electronDDWallFluxMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        electronDDWallFluxMixedFvPatchScalarField
        (
            const electronDDWallFluxMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        electronDDWallFluxMixedFvPatchScalarField
        (
            const electronDDWallFluxMixedFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        electronDDWallFluxMixedFvPatchScalarField
        (
            const electronDDWallFluxMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new electronDDWallFluxMixedFvPatchScalarField(*this)
            );
        }

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new electronDDWallFluxMixedFvPatchScalarField(*this, iF)
            );
        }

        // Member Functions

            virtual void updateCoeffs();

            //- Return the surface charging switch value
            inline bool enableSurfaceCharging() const
            {
                return enableSurfaceCharging_;
            }

            //- Write dictionary entries
            virtual void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
