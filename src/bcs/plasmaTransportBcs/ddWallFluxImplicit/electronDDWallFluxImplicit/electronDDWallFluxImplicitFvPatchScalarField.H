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
    Foam::electronDDWallFluxImplicitFvPatchScalarField

Description
    Concrete implementation of ddWallFluxImplicitFvPatchScalarField for 
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
        type            electronDDWallFluxImplicit;

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
    electronDDWallFluxImplicitFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef electronDDWallFluxImplicitFvPatchScalarField_H
#define electronDDWallFluxImplicitFvPatchScalarField_H

#include "ddWallFluxImplicitFvPatchScalarField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class electronDDWallFluxImplicitFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class electronDDWallFluxImplicitFvPatchScalarField
:
    public ddWallFluxImplicitFvPatchScalarField
{
private:
    // Private Member Functions

        //- Build seec_ list from plasmaSpecies registry (called once)
        void buildSEECList() const;

        //- Compute total SEE flux field on this patch [1/m^2/s]
        tmp<scalarField> calcSEEFlux() const;

protected:

    // Protected data

       //- Switch to enable surface charging
        Switch enableSurfaceCharging_;

        //- Switch to include drift velocity in wall flux calculation
        Switch includeDriftFlux_;

        //- Switch to enable secondary electron emission
        Switch enableSEE_;

        //- Default SEE coefficient
        scalar defaultSEEC_;

        //- Species-specific SEE coefficients
        dictionary speciesSEEC_;

        //- List to keep the SEE coefficients for each ion species
        mutable List<scalar> seec_;

        //- Switch to read the SEE entries only once
        mutable Switch mapped_;

        //- Electron reflection probability r of Hagelaar eq. (6.2).
        //  SAME KEY as the Mixed family (G3): one concept, one name.
        scalar electronReflection_;

    // Protected Member Functions

        //- The three factors by which REFLECTION scales this condition.
        //
        //  Hagelaar HDR (Literature/) ch. 6. Eq. (6.2) builds Gamma_w from the
        //  total loss speed and (6.6) feeds it back, so the closure solved
        //  exactly and UNCLAMPED is
        //
        //      Gamma = [(1-r)/(1+r)]*n*(A + Dd)  -  [2/(1+r)]*Gc
        //
        //  with A the thermal base, Dd the drift addition and Gc the creation
        //  flux. Every factor here is a CONSTANT, so unlike the Mixed family's
        //  form this stays fully implicit in n -- no lag is introduced.
        //
        //  NOTE THE FACTOR 2 ON EMISSION AT r = 0. It is not a fudge: it is
        //  eq. (6.1) closed with (6.6),
        //      Gamma = n*w_w - Gamma_w = n*(A - Gamma_w/n) - Gamma_w
        //            = n*A - 2*Gamma_w .
        //  MEASURED DIVERGENCE, 2026-09-04: commit 70f6674 put that second
        //  factor into the Mixed family (inside the loss speed) and not here,
        //  so this family under-counted every emitted electron by exactly 2
        //  from 70f6674 until today.
        scalar reflectionFluxFactor() const
        {
            return (1.0 - electronReflection_)/(1.0 + electronReflection_);
        }

        scalar reflectionCreationFactor() const
        {
            return 2.0/(1.0 + electronReflection_);
        }

        //- Return the thermal wall loss speed, eq. (6.6), scaled by
        //  reflectionFluxFactor(). Drift is scaled by the same factor in
        //  valueInternalCoeffs(), and creation by
        //  reflectionCreationFactor() in gradientBoundaryCoeffs().
        //
        //  LIMITATION, STATED RATHER THAN HIDDEN: this family does NOT carry
        //  eq. (6.6)'s max(..., 0) clamp. The clamp needs Gamma_w/n with n at
        //  the wall, which would make an intentionally-linear condition
        //  lagged and nonlinear. Consequence: under emission strong enough
        //  that Gamma_w/n exceeds A, this family lets the net wall flux become
        //  an electron SOURCE, where the Mixed family clamps it to the
        //  drift-only value. Hagelaar's reason for the clamp -- published
        //  gamma values "have generally been deduced neglecting thermal
        //  electron loss to the cathode" -- applies to both. Use the Mixed
        //  family (the default) at a cathode with significant SEE.
        virtual tmp<scalarField> calcWallThermalVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T
        ) const override;

public:

    //- Runtime type information
    TypeName("electronDDWallFluxImplicit");

    // Constructors

        //- Construct from patch and internal field
        electronDDWallFluxImplicitFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        electronDDWallFluxImplicitFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        electronDDWallFluxImplicitFvPatchScalarField
        (
            const electronDDWallFluxImplicitFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        electronDDWallFluxImplicitFvPatchScalarField
        (
            const electronDDWallFluxImplicitFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        electronDDWallFluxImplicitFvPatchScalarField
        (
            const electronDDWallFluxImplicitFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new electronDDWallFluxImplicitFvPatchScalarField(*this)
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
                new electronDDWallFluxImplicitFvPatchScalarField(*this, iF)
            );
        }

        // Member Functions

            //- Coefficient for fvm::div diagonal
            virtual tmp<Field<scalar>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Coefficient for fvm::laplacian diagonal
            virtual tmp<Field<scalar>> gradientInternalCoeffs() const;

            //- RHS coefficient for fvm::laplacian
            virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const;

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
