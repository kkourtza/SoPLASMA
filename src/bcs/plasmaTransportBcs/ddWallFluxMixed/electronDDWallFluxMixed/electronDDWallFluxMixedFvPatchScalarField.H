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

        //- ELECTRON REFLECTION COEFFICIENT r [-], default 0.
        //
        //  A wall does not absorb every electron that reaches it. With a
        //  probability r of elastic/inelastic backscatter, the NET thermal
        //  flux is reduced:
        //
        //      Gamma_net = (1 - r)/(1 + r) * (1/4) n v_th
        //
        //  which is the standard kinetic result (see e.g. Lieberman &
        //  Lichtenberg, "Principles of Plasma Discharges", ch. 6): of the
        //  incident flux only (1-r) is absorbed, while the near-wall density
        //  supports both the incident and reflected populations, hence (1+r).
        //
        //  r = 0 gives a factor of EXACTLY 1, so the default reproduces every
        //  previous result bit-for-bit. That is the regression this was
        //  verified against.
        //
        //  SCOPE, STATED BECAUSE IT IS PARTIAL. The factor is applied to the
        //  THERMAL flux only -- the term ¼ n v_th, in both the effective and
        //  the absorption velocity, so the mixed condition stays consistent.
        //  It is NOT applied to:
        //
        //    * the DRIFT term. Whether directed drift into a wall should carry
        //      the same reflection factor is a modelling choice, not an
        //      algebraic consequence, and treatments differ.
        //    * the EMISSION terms (ion-induced SEE and the emission models),
        //      which in a consistent treatment carry 2/(1+r). NOT APPLIED, by
        //      decision of 2026-09-04, and the derivation is recorded here so
        //      the omission is understood rather than inherited.
        //
        //      DERIVATION. Take the near-wall distribution as two
        //      half-Maxwellians at one temperature, amplitudes A_in (toward
        //      the wall) and A_out (away):
        //          Gamma_in = (1/4) A_in v_th,   n = (1/2)(A_in + A_out)
        //      With reflection r and an emitted flux Gamma_emit =
        //      (1/4) A_emit v_th,
        //          A_out = r A_in + A_emit
        //      so the net flux is (1/4) v_th [A_in (1-r) - A_emit]. Eliminating
        //      A_in = (2n - A_emit)/(1+r) gives
        //          Gamma = (1-r)/(1+r) (1/2) n v_th - 2/(1+r) Gamma_emit
        //
        //      BOTH FACTORS COME FROM THE SAME ELIMINATION, so they are a
        //      MATCHED PAIR: the emitted population contributes to the
        //      outgoing flux AND to the near-wall density, and that double
        //      role is where the 2 comes from. Applying one without the other,
        //      as here, is therefore INTERNALLY INCONSISTENT AT r > 0 -- named
        //      plainly so nobody has to rediscover it.
        //
        //      Note also that the derivation yields (1/2) n v_th, whereas this
        //      family uses (1/4) n v_th -- the full-Maxwellian convention,
        //      which ignores the depletion of the outgoing half. So even the
        //      thermal factor is a hybrid of two conventions.
        //
        //      AT r = 0 ALL OF THIS IS MOOT: both factors reduce to their
        //      identity values and the condition is exactly what it was.
        //
        //  So r > 0 is a PARTIAL model today and the report says so.
        scalar electronReflection_;

        //- THE SURFACE'S MATERIAL, declared once on the patch.
        //
        //  Inherited by every emission mechanism on this patch, so a surface
        //  says what it is MADE OF once and each mechanism resolves whatever
        //  it needs -- work function, secondary yield, band structure -- from
        //  the shipped libraries. A mechanism that names its own `material`
        //  overrides it.
        //
        //  This indirection is necessary rather than tidy: an emission model
        //  is handed only its OWN sub-dictionary, so without merging the
        //  material in it cannot see a `material` declared beside `emission`
        //  and fails asking for a work function the case did supply.
        word material_;

        //- THE EMITTED FLUX [1/m2/s], computed BEFORE the base class builds
        //  the loss speed, because eq. (6.8) needs it there.
        //
        //  Hagelaar eq. (6.1) is  Gamma = n w_w - Gamma_w, and eq. (6.6)/(6.8)
        //  put Gamma_w INSIDE the loss speed too:
        //      w_w = max( thermal - Gamma_w/n , 0 ) + max( drift_to_wall , 0 )
        //  This condition imposes n*(uDrift_n + uEff), i.e. n*w_w, so the
        //  separate -Gamma_w of eq. (6.1) must ALSO be applied -- which is what
        //  the refGradient term does. Emission therefore appears in BOTH
        //  places, and that is not double counting: it is exactly where the
        //  factor of 2 on emission comes from. Unclamped and with r = 0,
        //      Gamma = n*thermal + n*drift - 2*Gamma_w
        //  which is eq. (6.6) closed with (6.1), as derived.
        mutable scalarField emitted_;

        //- Whether emitted_ holds a value for the current evaluation
        mutable bool emittedValid_;

        //- The models, built on first updateCoeffs()
        mutable PtrList<emissionModel> emission_;

        //- Whether the emission report has been printed
        mutable bool emissionReported_;

        //- List to keep the SEE coefficients for each ion species
        List<scalar> seec_;

        //- Switch to read the SEE entries only once
        Switch mapped_;

    // Protected Member Functions

        //- The density that normalises Gamma_w/n in eqs (6.6), (6.8), (6.15).
        //
        //  For n_e that is THIS field. The energy override must supply n_e
        //  instead, because Gamma_w/n_eps is not the quantity in any of those
        //  equations -- see the measured defect recorded there.
        virtual tmp<scalarField> normalisingDensity() const;

        //- TOTAL wall loss speed W, eqs (6.1), (6.2), (6.6), (6.8) of
        //  Hagelaar's HDR (Literature/, ch. 6, pp. 36-37), solved exactly.
        //
        //  Eq. (6.2) builds Gamma_w from the TOTAL loss speed w_w, and (6.8)
        //  feeds Gamma_w back into w_w, so REFLECTION IS IMPLICIT, not a
        //  prefactor:
        //
        //      w_w = max(A - Gamma_w/n, 0) + Dd,   Gamma_w/n = r*w_w + Gc/n
        //
        //  with A the thermal base (6.6), Dd the drift addition of (6.8) and
        //  Gc the creation flux. The two branches of the max partition
        //  cleanly -- the unclamped condition reduces to exactly the negation
        //  of the clamped one -- so the implicit relation has a CLOSED FORM:
        //
        //      w_w = max( (A + Dd - Gc/n)/(1 + r), Dd )
        //      W   = (1 - r)*w_w                            [eq. (6.1)+(6.2)]
        //
        //  VERIFIED against the one reflection case Hagelaar closes
        //  explicitly: with Gc = 0 and Dd = 0 this gives
        //      Gamma.n = (1/sqrt(pi))*(1-r)/(1+r)*n*vT
        //  which is eq. (6.7) exactly, including the 1/sqrt(pi).
        //
        //  Combining r with (6.8)'s drift term is NOT quoted from the text --
        //  (6.7) is the only reflection case closed there. It follows from
        //  reading (6.2) literally, which uses the total w_w and therefore
        //  reflects the drift-driven incident flux too. Stated rather than
        //  presented as a citation. Derived 2026-09-04.
        //
        //  At r = 0 this reproduces the pre-2026-09-04 form, because
        //      max(A + Dd - Gc/n, Dd) == Dd + max(A - Gc/n, 0).
        //  ALGEBRAICALLY identical, NOT bitwise: the two forms associate the
        //  additions differently. Measured 2026-09-04 by testWallFlux over a
        //  372-point (Dd, Gc/n) sweep: max relative difference 3.7e-16, max
        //  absolute 2.3e-10 m/s against A = 4.7e5 m/s.
        //
        //  gRatio returns (Gamma_w/n)/A, the dimensionless creation-plus-
        //  reflection strength that eq. (6.15) needs for the energy weight.
        //  Zero where A is not positive.
        tmp<scalarField> wallLossSpeed
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal,
            scalarField& gRatio
        ) const;

        virtual tmp<scalarField> calcAbsorptionVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal
        ) const override;

        //- Calculate the electron-specific wall effective velocity
        //- The declared electron reflection coefficient (see the base).
        virtual scalar reflectionCoefficient() const
        {
            return electronReflection_;
        }

        virtual tmp<scalarField> calcEffectiveWallVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal
        ) const override;

public:

        //- THE CLOSURE, as a pure function of four scalars.
        //
        //  Public and static so a unit bed can exercise the SHIPPED algebra
        //  against analytic ground truth instead of a re-implementation of
        //  it. src/applications/utilities/testWallFlux does exactly that.
        //
        //      A       thermal base, eq. (6.6) = vT/sqrt(pi)      [m/s]
        //      Dd      drift addition of eq. (6.8), >= 0          [m/s]
        //      GcOverN creation flux over density, Gc/n           [m/s]
        //      r       reflection probability of eq. (6.2), [0,1)
        //      W       TOTAL loss speed (1-r)*w_w                 [m/s]
        //      gRatio  (Gamma_w/n)/A, for eq. (6.15)              [-]
        //
        //  See wallLossSpeed() for the derivation and the eq. (6.7) check.
        static void hagelaarClosure
        (
            const scalar A,
            const scalar Dd,
            const scalar GcOverN,
            const scalar r,
            scalar& W,
            scalar& gRatio
        )
        {
            const scalar wW = max((A + Dd - GcOverN)/(1.0 + r), Dd);

            W = (1.0 - r)*wW;

            // gRatio TAKES w_w FROM EQ. (6.6), NOT FROM wW ABOVE.
            //
            // Eq. (6.15), which hagelaarEnergyWeight() implements, ends
            // "where w_w from equation (6.6)" -- the DRIFT-FREE loss speed.
            // wW above is eq. (6.8), which adds the drift OUTSIDE the max.
            // Feeding eq. (6.8) here makes gRatio grow without bound in the
            // field, so max(0, 5/3 - (2/3)gRatio) reaches ZERO and the wall
            // stops carrying ANY energy away while still absorbing electrons.
            // MEASURED 2026-09-07 on grubert2009: at r = 0.36 that cliff is at
            // Dd/A = 8.44, and crossing it took Te,max from 11.6 eV to the
            // mean-energy clamp. Isolated by splitting r between the two
            // conditions: r on nEps_e alone reproduced it, r on n_e alone was
            // indistinguishable from r = 0.
            //
            // With eq. (6.6) and no creation this is r/(1+r) -- a CONSTANT,
            // and exactly the expression Hagelaar substitutes himself on p.39
            // to check the limit: Gamma_w = r n w_w = (vT/sqrt(pi)) r/(1+r)
            // into (6.15) with r -> 1 yields 2 Te. Verified: it does, to all
            // digits. gRatio is then bounded in [0, 1/2] for r in [0, 1], so
            // eps_w stays in [2 Te, 5/2 Te] and the max() cannot fire on
            // reflection alone -- which is right, because the text puts that
            // clamp there for strong wall CREATION (secondary emission at the
            // cathode), not for reflection.
            //
            // IDENTICAL to the old expression at Dd = 0 for every r, and at
            // r = 0 with no emission for every Dd. It therefore cannot move
            // any result this framework has already validated: every prior
            // case ran r = 0 with no anode SEE, and at a cathode the electron
            // drift is outward so Dd = max(0, uDriftNormal) = 0 anyway.
            const scalar wW66 = max((A - GcOverN)/(1.0 + r), scalar(0));

            gRatio = (A > VSMALL) ? (r*wW66 + GcOverN)/A : 0.0;
        }

        //- eps_w/eps of eq. (6.15), clamped at zero. Static for the same
        //  reason as hagelaarClosure(); f0 is the zero-creation weight, 5/3
        //  for a Maxwellian.
        static scalar hagelaarEnergyWeight(const scalar f0, const scalar gRatio)
        {
            return max(scalar(0), f0 - (2.0/3.0)*gRatio);
        }


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
