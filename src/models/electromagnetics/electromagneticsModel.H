/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::electromagneticsModel
    
Description
    Abstract base class for electromagnetic and electrostatic field solvers.

    Owns all fields common to plasma EM models: electric potential, electric
    field, charge density, and surface charge. These are persistent, registered
    in the mesh registry, and read from disk where present.

    For pure electrostatic runs, chargeDensity and surfCharge are zero or
    read from disk and never modified. For plasma runs, the calling solver
    updates them via mutable references before each solve.

    The class has no knowledge of plasma species or transport models.
    Conductivity and diffusive sources are passed as raw fields by the solver.

SourceFiles
    electromagneticsModelI.H
    electromagneticsModel.C

\*---------------------------------------------------------------------------*/

#ifndef electromagneticsModel_H
#define electromagneticsModel_H

#include "floatingElectrode.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvc.H"
#include "fvm.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "UPtrList.H"
#include "runTimeSelectionTables.H"

#include "SoPLASMAConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class electromagneticsModel Declaration
\*---------------------------------------------------------------------------*/

class electromagneticsModel
:
    public IOdictionary
{
private:

    //- Disallow copy
    electromagneticsModel(const electromagneticsModel&) = delete;
    void operator=(const electromagneticsModel&) = delete;

protected:

    //- Reference to the primary gas/plasma mesh
    const fvMesh& mesh_;

    //- Reference to the list of dielectric meshes
    const UPtrList<fvMesh>& dielectricMeshes_;

    //- Electric potential field [V]
    volScalarField ePotential_;

        //- The floating-electrode constraint, when the case has one.
        //
        //  Constructed LAZILY on the first solve, not in this constructor: on
        //  a multi-region model the dielectric regions do not exist until the
        //  derived constructor has run, and the classification scans them.
        mutable autoPtr<floatingElectrode> floating_;

        //- Whether the case has been searched for a floating conductor
        mutable bool floatingChecked_ = false;

    //- Electric field [V/m]
    volVectorField E_;

    //- Magnitude of the electric field [V/m]
    volScalarField Emag_;

    //- Flux of the electric field [V * m2]
    surfaceScalarField phiE_;

    //- Reduced electric field [V * m2]
    volScalarField reducedE_;

    //- Voumetric charge density in the gas region [C/m3]
    volScalarField chargeDensity_;

    //- Surface charge density (owned by gas) [C/m2]
    volScalarField surfCharge_;

    //- Dielectric permittivity in the gas region [F/m]
    dimensionedScalar epsilon_;

    //- Relative permittivity of the medium (gas region)
    scalar epsilonR_;

public:

    //- The GLOBAL Poisson numerics.
    //
    //  These are uniform over every region the potential is solved on -- the
    //  gas and every dielectric -- because the model holds exactly one of each
    //  and reads them once. They are therefore not per-region settings and do
    //  not belong in a per-region file.
    //
    //  They live in `system/plasmaSimulationControls` under `poisson`, beside
    //  `plasmaTimeControl` and `outerCoupling`, because they are SOLUTION
    //  CONTROL. That is OpenFOAM's own split: `constant/` holds physical
    //  properties, `system/` holds how they are solved. Permittivity is
    //  physical and lives per region under `constant/<region>/`; a
    //  discretisation scheme is not.
    struct poissonNumerics
    {
        //- `explicit` treats the charge as a pure source; `semiImplicit` also
        //  brings the conductivity in implicitly, which stabilises the
        //  space-charge / Poisson coupling. Note this is a SCHEME, not a
        //  linear solver -- the linear solver for ePotential is in fvSolution,
        //  and the old case-variable name `poissonSolver` invited exactly that
        //  confusion.
        word scheme = "semiImplicit";

        //- How E is formed from the potential: `grad` | `reconstruct`
        word EScheme = "reconstruct";

        //- Non-orthogonal correctors for the Poisson solve
        label nNonOrthogonalCorrectors = 0;

        //- Segregated-coupling Picard loop, used only when the interface is
        //  explicit. Irrelevant under monolithic coupling.
        label maxNonCoupledIterations = 100;
        scalar nonCoupledTolerance = 1e-6;
    };

    //- Read the global Poisson numerics.
    //
    //  `legacyCoeffs` is the old `<type>Coeffs` sub-dictionary of
    //  constant/electromagneticsProperties. It is consulted ONLY as a fallback,
    //  so an unmigrated case keeps running, and using it prints where the
    //  setting has moved to.
    static poissonNumerics readPoissonNumerics
    (
        const fvMesh& mesh,
        const dictionary& legacyCoeffs
    );

    //- One region's relative permittivity, from
    //  `constant/<region>/electricalProperties`.
    //
    //  Per-region, per OpenFOAM's own convention -- the same place and shape as
    //  `constant/<region>/thermophysicalProperties` in chtMultiRegionFoam. This
    //  replaces a GLOBAL file that carried a sub-dictionary per region and was
    //  then copied byte-identically into every region directory, so every copy
    //  held every region's data.
    //
    //  `legacyCoeffs` is consulted only as a fallback: for the gas its own
    //  `dielectricConstant`, for a dielectric the sub-dictionary named after
    //  the region.
    //
    //  The KIND is looked up here rather than passed in, because the caller
    //  iterating the dielectric mesh list has no way to know which of those
    //  meshes are `farField`. Defaulting differs by kind:
    //      gas        -> 1.0 (a gas is a vacuum to a few parts in 10^4)
    //      farField   -> 1.0 (that IS the definition of the kind)
    //      dielectric -> FATAL, since a default of 1 would be vacuum and would
    //                    make a real barrier electrically invisible
    static scalar readEpsilonR
    (
        const fvMesh& regionMesh,
        const dictionary& legacyCoeffs
    );

    //- The KIND declared for this mesh region in constant/regionProperties.
    //
    //  Returns `gas` when there is no regionProperties at all, which is the
    //  right answer for a single-region case.
    //
    //  Reads the `regions` entry into a HashTable<wordList> directly rather
    //  than using Foam::regionProperties, so this header needs no meshTools
    //  dependency.
    static word regionKind(const fvMesh& regionMesh);

    //- The region kinds this solver family implements.
    //
    //      gas         the plasma: species, chemistry, electron energy, Poisson
    //      dielectric  a solid barrier: Poisson only, permittivity REQUIRED
    //      farField    a FICTITIOUS air/vacuum region: Poisson only,
    //                  epsilonR defaults to 1.0. Its purpose is to extend the
    //                  ELECTROSTATIC domain without extending the plasma one --
    //                  Poisson is long-range, so it needs a large domain, while
    //                  species transport and photoionization do not.
    static const wordList& knownRegionKinds();

    //- Reject any region KIND that no solver here implements.
    //
    //  `kinds` is the top-level key list of constant/regionProperties, i.e.
    //  `rp.sortedToc()`. Takes a wordList rather than the regionProperties
    //  object so this header needs no meshTools dependency.
    //
    //  Lives here because THIS class is what degrades silently: it derives
    //  singleRegionPoisson vs multiRegionPoisson from whether any dielectric
    //  region exists, so a kind that matches nothing produces a case with no
    //  dielectrics rather than an error. Both solvers call it.
    static void validateRegionKinds(const wordList& kinds);

    //- Runtime type information
    TypeName("electromagneticsModel");

    //- Declare the Run-Time Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        electromagneticsModel,
        dictionary,
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes
        ),
        (mesh, dielectricMeshes)
    );


    // Constructors

        //- Construct from meshes
        electromagneticsModel
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes
        );

    // Selectors
    
        //- Selector
        static autoPtr<electromagneticsModel> New
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes
        );

    //- Destructor
    virtual ~electromagneticsModel() = default;

    // Solve functions
    
        //- Solve Poisson equation (for explicit Poisson)
        virtual void solve() = 0;

        //- Solve semi-implicit Poisson equation (for semiImplicit Poisson)
        virtual void solve
        (
            const volScalarField& electricalConductivity,
            const volScalarField& diffusiveChargeSource
        ) = 0;

    // Field access
    
        const volScalarField& ePotential() const { return ePotential_; }

        //- Non-const access to the potential.
        //
        //  Deliberately named apart from ePotential() so every mutation of the
        //  solved field is greppable. TWO legitimate users:
        //
        //   1. the floating-electrode constraint, which adds dV_f*psi AFTER
        //      the solve; that is an exact superposition, not a fix-up,
        //      because psi solves the homogeneous problem;
        //   2. plasmaExternalCircuit, which sets a driven electrode's
        //      BOUNDARY value from V_source - R*I before the next step. It
        //      touches the boundary field only, never the interior, so it
        //      changes the problem being solved rather than the solution --
        //      which is what an electrode fed through a ballast IS.
        //
        //  Added 2026-09-06; the comment said "one legitimate user" and would
        //  otherwise have read as a prohibition on the second.
        volScalarField& ePotentialRef() { return ePotential_; }

        //- Apply the floating-electrode constraint, if the case has one.
        //
        //  MUST be called after every Poisson solve and BEFORE the derived
        //  fields are formed, because it changes the potential. Pass the
        //  effective gas permittivity under `scheme semiImplicit` and nullptr
        //  under `explicit`, so psi is built from the same operator the solve
        //  used -- superposition is exact only then.
        //
        //  Called from BOTH Poisson models. A correction applied in only one
        //  of two sibling models is this project's most repeated defect, so
        //  the call sites are enumerated here: singleRegionPoisson::solve()
        //  and solve(sigma,...), multiRegionPoisson::solve() and
        //  solve(sigma,...).
        void correctFloatingElectrode(const volScalarField* effEpsGas);

        //- The floating electrode, or nullptr
        const floatingElectrode* floatingElectrodePtr() const
        {
            return floating_ ? &floating_() : nullptr;
        }

        //- The floating electrode, or nullptr. Non-const, for the charge
        //  ledger: plasmaTransport feeds I_plasma in from the species wall
        //  fluxes, which is the only place they exist.
        floatingElectrode* floatingElectrodePtr()
        {
            return floating_ ? &floating_() : nullptr;
        }

        //- Classify the potential's IMPOSED-VALUE patches, across EVERY
        //  region, into the three things a conductor can be.
        //
        //  ONE OWNER for this rule. Sato's discharge current needs the driven
        //  and grounded lists; the floating electrode needs the floating one
        //  AND the other two (as the patches its unit-potential problem holds
        //  at zero). Two copies of a classification rule is two rules.
        //
        //      driven    a Dirichlet condition that is TIME-VARYING, or a
        //                NON-ZERO constant
        //      grounded  a NON-time-varying Dirichlet equal to ZERO
        //      floating  a conductor whose potential is UNKNOWN
        //      (skipped) everything else -- zeroGradient, empty, processor,
        //                the region interfaces, thinDielectricPotential
        //
        //  The time-varying test is load-bearing: at t = 0 a ramp reads
        //  EXACTLY ZERO, so a value-only rule would call the driven electrode
        //  of every ramped case a ground.
        void classifyElectrodePatches
        (
            wordList& driven,
            wordList& grounded,
            wordList& floating
        ) const;
        const volVectorField& E() const { return E_; }
        const volScalarField& Emag() const { return Emag_; }
        const surfaceScalarField& phiE() const { return phiE_; }
        const volScalarField& reducedE() const { return reducedE_; }
        const volScalarField& chargeDensity() const { return chargeDensity_; }
        volScalarField& chargeDensity() { return chargeDensity_; }
        const volScalarField& surfCharge() const { return surfCharge_; }
        volScalarField& surfCharge() { return surfCharge_; }
        const dimensionedScalar& epsilon() const { return epsilon_; }
        scalar epsilonR() const { return epsilonR_; }

        //- Gas number density used to reduce E into E/N.
        //
        //  Pushed in by plasmaSpecies, which owns it. It must NOT be read
        //  independently from electromagneticsProperties: the chemistry and
        //  the transport coefficients are looked up at E/N, so a second
        //  definition of N means the rates are evaluated at a different gas
        //  state than the one the case declares -- silently, and by whatever
        //  factor the two happen to differ. That is what happened here: the
        //  Poisson model defaulted to 1 atm while the species model derived
        //  1 bar from the ideal gas law, and every rate lookup in the run was
        //  off by 1.3%.
        virtual void setBackgroundDensity(const dimensionedScalar& N) = 0;
        virtual void setBackgroundDensity(const volScalarField& N) = 0;

        virtual const word& PoissonScheme() const = 0;

        //- Change the Poisson discretisation AT RUN TIME.
        //
        //  Exists for exactly one caller: soPlasmaFoam forcing `explicit` when
        //  the Newton outer solver takes over. The semi-implicit form PREDICTS
        //  rho^{n+1} and folds it into the operator, which is what a SEGREGATED
        //  solver needs and what a Newton residual -- which already carries
        //  rho^{n+1} as an unknown -- must not have, because the charge
        //  relaxation is then counted twice and the coupled system has no
        //  consistent root (docs/design/newton-ignition-experiments.md section 23).
        //
        //  SAFE TO CHANGE MID-RUN, and the reason is specific to this field:
        //  ePotential has NO time derivative. It is recomputed from the charge
        //  density every step, so there is no history that a change of
        //  discretisation could leave inconsistent -- phi simply jumps by the
        //  O(dt*sigma/eps) difference between the two operators and is
        //  thereafter consistent. The same change to a TRANSPORTED field, which
        //  carries fvm::ddt and an old-time history, would NOT be safe.
        virtual void setPoissonScheme(const word& scheme) = 0;

        //- Refresh E, phiE, Emag, reducedE from the current ePotential.
        //
        //  Hoisted to the base 2026-09-11. Both concrete models already had
        //  it; only the CAST to singleRegionPoisson forced the Newton outer
        //  solver to know which model it had. It calls this and nothing else
        //  on the concrete type, so making it virtual is what lets Newton
        //  support multiRegionPoisson at all.
        virtual void updateDerivedFields() = 0;
        virtual label nNonOrthCorr() const = 0;
        virtual const word& EScheme() const = 0;

        //- Number of EXTRA Poisson-only regions beyond the gas mesh:
        //  dielectrics, and the fictitious far-field regions that solve
        //  Poisson alone to get the correct far-field distribution.
        //  ZERO for singleRegionPoisson.
        //
        //  These exist for the Newton outer solver, which packs them into the
        //  RAGGED TAIL of its DOF layout (phi lives in every region; the
        //  species do not, so the tail extends the Poisson block ONLY).
        //
        //  The inter-region coupling needs no monolithic assembly.
        //  coupledElectricPotential::updateCoeffs() sets valueFraction,
        //  refValue (the mapped neighbour potential) and refGrad (the
        //  surface-charge jump) REGARDLESS of useImplicit -- useImplicit only
        //  additionally fills source() for manipulateMatrix -- so an explicit
        //  fvc:: residual already sees the neighbouring region. Measured
        //  2026-09-11 on needleDBD: a +1 V bump of the dielectric potential
        //  moves the explicit gas residual by rms 1.549 on the 160 interface
        //  cells against 0.0035 in the interior. Newton's own outer iteration
        //  then converges the region coupling, which is exactly what the
        //  monolithic assembly does inside the LINEAR solve -- but at the
        //  nonlinear level, so the surface-charge and permittivity
        //  nonlinearity come with it.
        virtual label nExtraPoissonRegions() const { return 0; }

        //- Potential field of extra Poisson region i. Non-const: Newton
        //  writes the trial state into it.
        virtual volScalarField& extraPoissonPotential(const label i);

        //- Permittivity of extra Poisson region i.
        virtual const dimensionedScalar& extraPoissonEpsilon
        (
            const label i
        ) const;


    // Mesh and dictionary access

        //- Access to the main mesh
        inline const fvMesh& mesh() const;

        //- Access to the list of dielectric meshes
        inline const UPtrList<fvMesh>& dielectricMeshes() const;

        //- Access to the dictionary
        inline const dictionary& dict() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "electromagneticsModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
