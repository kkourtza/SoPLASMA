/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::plasmaDischargeCurrent

Description
    External circuit current by Sato's equation (Morrow & Sato, J. Phys. D 32
    (1999) L20).

        I(t) = SUM_s q_s INT_V Gamma_s . e_hat dV  +  C_g dV_a/dt

        e_hat = E_L / V_a          [1/m]   geometric weighting field
        C_g   = INT_V eps |e_hat|^2 dV     gap capacitance

    `E_L` is the LAPLACIAN (space-charge-free) field. The total field is
    E = E_L + E', and the whole point of the derivation is that E' CANCELS out
    of the final expression -- so using the solved field in place of E_L would
    be wrong, not merely approximate.

    Because e_hat is purely geometric it is computed ONCE at construction and
    cached. It changes only if the mesh or the electrode set changes.

    This is the Shockley-Ramo theorem generalised from a single drifting charge
    to a continuum of species; e_hat is its "weighting field". The equivalent
    electrode surface integral,

        INT_electrode J_tot . dS = - INT_V e_hat . J_tot dV

    follows from div(J_tot) = 0 plus Green's reciprocity, and is implemented
    separately as a VALIDATION route only. It is the numerically worse of the
    two: it samples one layer of cells in the steepest-gradient region of the
    mesh, and it differences conduction against displacement current where the
    two nearly cancel. Sato's volume form integrates a smooth weight and
    collapses the displacement term to an analytic constant. See
    docs/discharge-current-sato-plan.md Sect. 5.1.

SourceFiles
    plasmaDischargeCurrent.C
    weightingField.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaDischargeCurrent_H
#define plasmaDischargeCurrent_H

#include "fvCFD.H"
#include "fvMesh.H"
#include "wordList.H"
#include "OFstream.H"
#include "plasmaSpecies.H"
#include "plasmaTransport.H"
#include "plasmaConstants.H"

namespace Foam
{

class plasmaSpecies;
class plasmaTransport;
class electromagneticsModel;

/*---------------------------------------------------------------------------*\
                  Class plasmaDischargeCurrent Declaration
\*---------------------------------------------------------------------------*/

class plasmaDischargeCurrent
{
    // Private Data

        const fvMesh& mesh_;

        //- Off unless a `dischargeCurrent` sub-dictionary asks for it, so no
        //  existing case changes behaviour or cost.
        bool enabled_ = false;

        //- Patch driven by the applied voltage, and the patches held at zero.
        word drivenPatch_;
        wordList groundedPatches_;

        //- Geometric weighting field e_hat = E_L/V_a [1/m] on the GAS mesh,
        //  computed once. The conduction term needs it only here, because the
        //  charged species live only in the gas.
        autoPtr<volVectorField> eHat_;

        //- The unit-voltage potential, per region. Kept because C_g integrates
        //  over the dielectrics too, and because they are worth writing out
        //  when a weighting field looks wrong.
        autoPtr<volScalarField> psiHatGas_;
        PtrList<volScalarField> psiHatDiel_;

        //- Gap capacitance INT eps |e_hat|^2 dV [F], a constant of the mesh.
        scalar Cg_ = 0;

        //- Per-species conduction columns in the CSV.
        bool perSpecies_ = false;

        //- Output cadence. INDEPENDENT of the field write interval and may be
        //  far shorter: the current is one scalar per step, so writing it every
        //  step costs nothing next to a field write, while a nanosecond pulse
        //  is destroyed by sampling at the field cadence. 1 = every step.
        label writeInterval_ = 1;
        label printInterval_ = 0;

        //- Previous applied voltage and time, for dV_a/dt.
        scalar prevVoltage_ = 0;
        scalar prevTime_ = 0;
        bool   seeded_ = false;

        autoPtr<OFstream> os_;

        //- Per-species conduction current [A], for the optional columns.
        List<scalar> perSpeciesCurrent_;

        //- CROSS-CHECK: the same current by the electrode surface integral,
        //  which is mathematically equivalent and numerically unrelated.
        //
        //      INT_electrode J_tot . dS = - INT_V e_hat . J_tot dV
        //
        //  (div J_tot = 0 plus Green's reciprocity; outward normal). The sign
        //  is FIXED BY THAT IDENTITY, not chosen to make the numbers agree.
        //
        //  Off by default: it needs dE/dt at the boundary, which is the noisy
        //  quantity Sato's form exists to avoid, and it costs a stored field.
        bool crossCheck_ = false;

        //- eps E . Sf on the driven patch at the previous step, for dE/dt.
        scalarField prevEFlux_;
        bool eFluxSeeded_ = false;

        //- The surface-integral current [A], conduction and displacement.
        void surfaceIntegralCurrent
        (
            const plasmaTransport& transport,
            const plasmaSpecies& species,
            const electromagneticsModel& em,
            scalar& Icond,
            scalar& Idisp
        );

        //- 2*pi/theta on a wedge mesh, 1 otherwise. Per MESH, because
        //  nothing guarantees a dielectric region was built with the same
        //  wedge angle as the gas.
        //
        //  mesh.V() on a wedge holds only the slice volumes, so every
        //  EXTENSIVE quantity -- charge, capacitance, current -- must be
        //  scaled to the full revolution. A 5 degree wedge would otherwise
        //  understate the current by 72x.
        static scalar wedgeRevolutionFactor(const fvMesh& m);

        //- Cached value of the above.
        scalar revolutionFactor_ = 1.0;

        //- Solve div(eps grad psi) = 0 with psi = 1 on the driven patch and 0
        //  on the grounded ones, and set eHat_ = -grad(psi), Cg_ from it.
        void computeWeightingField(const electromagneticsModel& em);

        //- Area-weighted mean potential on the driven patch [V].
        scalar appliedVoltage(const electromagneticsModel& em) const;


public:

    // Constructors

        plasmaDischargeCurrent
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const electromagneticsModel& em
        );


    // Member Functions

        inline bool enabled() const { return enabled_; }

        //- Gap capacitance [F]. Zero when disabled.
        inline scalar Cg() const { return Cg_; }

        //- The geometric weighting field [1/m].
        const volVectorField& eHat() const { return eHat_(); }

        //- Write psi_hat and e_hat per region, for inspecting a weighting
        //  field that looks wrong. Cheap and written on demand only.
        void writeFields() const;

        //- Accumulate and write this step's current.
        void update
        (
            const plasmaTransport& transport,
            const plasmaSpecies& species,
            const electromagneticsModel& em
        );
};


} // End namespace Foam

#endif

// ************************************************************************* //
