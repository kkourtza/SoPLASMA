/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::thinDielectricPotentialFvPatchScalarField

Description
    *** WORK IN PROGRESS -- NOT BUILT. ***
    There is no .C and this file is deliberately NOT in src/bcs/Make/files, so
    nothing compiles or links it and the type is not registered. It is the
    INTERFACE DESIGN, kept in the code where an implementer will look for it;
    the physics is in
    BoltzmannSolver docs/models/poisson_equation/ (see also the memory note
    thin-dielectric-bc). Do not use it expecting a working boundary condition.

    Electric potential on an UNMESHED dielectric surface: an insulating sheet
    thin enough that meshing it is not worth it, optionally backed by a
    conductor.

    It fills a real gap. `coupledElectricPotential` is the only other boundary
    condition that consumes a surface charge, and it REQUIRES a neighbour
    region. On a patch with nothing meshed behind it there was therefore no way
    to feed accumulated charge back into the potential: charge accumulated in
    the field and did NOTHING -- no field, no shielding, no self-limiting. The
    obvious alternative, `zeroGradient`, is E.n = 0, a perfect insulator that
    ignores its own surface charge.

    THE CONDITION. Let n be the outward normal from the gas, and let a layer of
    thickness d and relative permittivity epsilonR sit outside the patch, backed
    by a conductor held at V_b. Thin means the field inside the layer is uniform,
    so the drop from the surface potential V_s to V_b over d gives a normal
    displacement C (V_s - V_b) with

        C = epsilon0 * epsilonR / d          surface capacitance per area

    Gauss's law at the gas/layer surface, with free surface charge sigma and n
    directed from the gas into the layer, reads D_layer.n - D_gas.n = sigma. With
    D_gas.n = -epsilon_gas dV/dn this gives a ROBIN condition:

        epsilon_gas dV/dn = sigma - C (V_s - V_b)

    THE TWO CASES ARE SEPARATE `kind`s at the case level, because one boundary
    condition with an optional thickness admits a plausible but MEANINGLESS
    configuration -- a free-standing sheet whose stated thickness changes
    nothing:

      thinDielectricSurface       FREE-STANDING, nothing behind it. There is no
                                  conductor to be a capacitor to, so C = 0:
                                      epsilon_gas dV/dn = sigma
                                  `thickness` and `backingPotential` are
                                  REJECTED here, which is the point.

      thinDielectricOnElectrode   BACKED -- a barrier bonded to an electrode,
                                  the usual DBD arrangement. `thickness` and
                                  `backingPotential` are both REQUIRED: they are
                                  meaningless apart, and requiring them together
                                  is what stops a case asking for a capacitance
                                  to nothing.

    LIMITS, which are the tests:

        C -> 0     epsilon_gas dV/dn = sigma      pure Neumann, sigma as source
        C -> inf   V_s -> V_b                     the layer vanishes and the
                                                  surface IS the backing
                                                  electrode (Dirichlet)
        sigma = 0, C -> inf                       V_s = V_b, a grounded
                                                  electrode

    So one implementation spans the whole range, and the thickness term is what
    gives the model a memory of d: without it, results cannot depend on layer
    thickness, while a real barrier's charging rate and saturation do.

    MAPPED ONTO mixedFvPatchField. With
    V_s = f*refValue + (1-f)*(V_c + refGrad/Delta) and Delta = deltaCoeffs,
    discretising dV/dn as (V_s - V_c)*Delta gives

        f        = C / (C + epsilon_gas*Delta)
        refValue = V_b + sigma/C
        refGrad  = 0

    and sigma/C diverges as C -> 0. The PRODUCT is finite --
    f*refValue -> sigma/(epsilon_gas*Delta) -- so the free-standing branch is
    implemented as pure refGradient instead, never forming refValue. See
    updateCoeffs.

Usage
    \table
        Property          | Description                        | Required
        epsilonR          | layer relative permittivity        | no (default 1)
        thickness         | layer thickness [m]                | only if backed
        backingPotential  | conductor potential behind it [V]  | only if backed
        surfCharge        | name of the surface-charge field   | no
    \endtable

    Free-standing:
    \verbatim
    <patchName>
    {
        type            thinDielectricPotential;
        epsilonR        4.6;
        surfCharge      surfCharge;
        value           uniform 0;
    }
    \endverbatim

    Backed by a grounded electrode:
    \verbatim
    <patchName>
    {
        type                thinDielectricPotential;
        epsilonR            4.6;
        thickness           100e-6;
        backingPotential    0;
        surfCharge          surfCharge;
        value               uniform 0;
    }
    \endverbatim

SourceFiles
    thinDielectricPotentialFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thinDielectricPotentialFvPatchScalarField_H
#define thinDielectricPotentialFvPatchScalarField_H

#include "mixedFvPatchFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class thinDielectricPotentialFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class thinDielectricPotentialFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Relative permittivity of the layer
        scalar epsilonR_;

        //- Layer thickness [m]. Negative means FREE-STANDING, i.e. no backing
        //  conductor and therefore no capacitance: the `thinDielectricSurface`
        //  case. Not defaulted to a number, because a guessed thickness would
        //  silently invent a capacitance.
        scalar thickness_;

        //- Potential of the conductor behind the layer [V]. Meaningful only
        //  when `thickness_` is positive.
        scalar backingPotential_;

        //- Name of the surface-charge field, or "none".
        word surfChargeName_;


    // Private Member Functions

        //- Surface capacitance per unit area, epsilon0*epsilonR/d.
        //  Zero for the free-standing case.
        scalar layerCapacitance() const;

        //- The surface charge on this patch, zero where none is named.
        tmp<scalarField> surfaceCharge() const;


public:

    //- Runtime type information
    TypeName("thinDielectricPotential");


    // Constructors

        thinDielectricPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        thinDielectricPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        thinDielectricPotentialFvPatchScalarField
        (
            const thinDielectricPotentialFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        thinDielectricPotentialFvPatchScalarField
        (
            const thinDielectricPotentialFvPatchScalarField&
        );

        thinDielectricPotentialFvPatchScalarField
        (
            const thinDielectricPotentialFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new thinDielectricPotentialFvPatchScalarField(*this)
            );
        }

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new thinDielectricPotentialFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- True when this condition consumes a surface charge.
        //
        //  Surface charging is DERIVED, not switched: charge accumulates on
        //  exactly those patches whose potential condition reads it. A
        //  dielectric surface therefore charges because something reads the
        //  charge, and a metal electrode does not because nothing does.
        bool consumesSurfaceCharge() const
        {
            return surfChargeName_ != "none";
        }

        virtual void updateCoeffs();

        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
