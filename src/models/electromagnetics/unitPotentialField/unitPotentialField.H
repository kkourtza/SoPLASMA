/*---------------------------------------------------------------------------*\
  File: unitPotentialField.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

Class
    Foam::unitPotentialField

Description
    A UNIT-POTENTIAL Laplace solve, monolithic across the gas and every
    dielectric region.

        div(eps grad(psi)) = 0,   psi = 1 on one patch, psi = 0 on the others

    psi is dimensionless and depends ONLY on the geometry and the
    permittivities, which is what makes it reusable and, when the operator is
    constant, cacheable.

    TWO features need exactly this field, with different patches held at 1:

      * Sato's discharge current  -- psi_hat with the MEASURING electrode at 1.
        Its capacitance is the ENERGY integral, INT eps |grad psi|^2 dV.
      * The floating electrode    -- psi_hat with the FLOATING electrode at 1.
        Its capacitance is the SURFACE integral, INT eps grad(psi).n dA, i.e.
        the self-capacitance C_self that converts a known charge into a
        potential.

    Extracted from plasmaDischargeCurrent::computeWeightingField 2026-09-03 so
    the second feature could not fork a second copy. Two facts make sharing
    non-optional rather than tidy: the interface-BC RETARGETING below, and the
    Krylov solver defaults, were each found by a measured wrong answer, and a
    copy would have inherited neither.

SourceFiles
    unitPotentialField.C

\*---------------------------------------------------------------------------*/

#ifndef unitPotentialField_H
#define unitPotentialField_H

#include "fvCFD.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class electromagneticsModel;

/*---------------------------------------------------------------------------*\
                     Class unitPotentialField Declaration
\*---------------------------------------------------------------------------*/

class unitPotentialField
{
    // Private data

        //- The GAS mesh, which owns the assembly
        const fvMesh& mesh_;

        //- Field name, used for the registry and for the fvSolution lookup
        const word name_;

        //- The patch held at psi = 1
        const word unitPatch_;

        //- The patches held at psi = 0
        const wordList zeroPatches_;

        //- psi on the gas mesh
        autoPtr<volScalarField> gas_;

        //- psi on each dielectric mesh
        PtrList<volScalarField> dielectrics_;

        //- Whether the last solve was monolithic across regions
        bool coupled_;

        //- Number of dielectric regions seen by the last solve
        label nDiel_;

        //- Whether a solve has happened at all
        bool solved_;

        //- Caller-supplied advice appended to the coupled-path fatal, so the
        //  error tells the user how to act on THEIR feature rather than on a
        //  helper class they never configured.
        const string hint_;


    // Private Member Functions

        //- Clone the potential's field and REBUILD the boundary conditions for
        //  a unit-potential problem. See the .C -- the interface retargeting
        //  here is load-bearing and was worth a 3.76x error.
        autoPtr<volScalarField> clone
        (
            const volScalarField& ePot,
            const word& regionName
        ) const;

        //- Linear-solver settings for this one-shot cold solve
        dictionary solverDict() const;

        //- The shared solve. `effEpsGas` null means use the model's constant
        //  gas permittivity.
        void solveImpl
        (
            const electromagneticsModel& em,
            const volScalarField* effEpsGas
        );

        //- No copying: it owns registered fields
        unitPotentialField(const unitPotentialField&) = delete;
        void operator=(const unitPotentialField&) = delete;


public:

    // Constructors

        //- Construct from the gas mesh and the patch assignment
        unitPotentialField
        (
            const fvMesh& gasMesh,
            const word& name,
            const word& unitPatch,
            const wordList& zeroPatches,
            const string& hint = string::null
        );


    //- Destructor
    ~unitPotentialField() = default;


    // Member Functions

        //- Solve the unit-potential problem with the model's permittivities.
        //  Safe to call repeatedly: it re-solves in place.
        //
        //  Correct for an EXPLICIT Poisson scheme, where the operator is eps.
        void solve(const electromagneticsModel& em);

        //- Solve with an EFFECTIVE gas permittivity.
        //
        //  REQUIRED under `scheme semiImplicit`, where the real operator is
        //  eps + dt*sigma_electrical rather than eps. Superposition is exact
        //  only if psi solves THE SAME OPERATOR as the field it corrects, so
        //  passing the true eps there would silently break it -- and since
        //  both dt and sigma change every step, psi and its capacitance must
        //  be rebuilt every step too.
        //
        //  Dielectrics keep their CONSTANT eps: an insulator has no
        //  conductivity, so eps_eff = eps there. This mirrors
        //  multiRegionPoisson exactly, which builds effEps for the gas and
        //  uses reg.epsilon() for each dielectric.
        void solve
        (
            const electromagneticsModel& em,
            const volScalarField& effEpsGas
        );

        //- INT eps |grad psi|^2 dV summed over EVERY region, with each
        //  region's own wedge revolution factor applied [F].
        //  This is Sato's gap capacitance.
        scalar energyCapacitance(const electromagneticsModel& em) const;

        //- INT eps grad(psi).n dA over `patchName` on the GAS mesh [F].
        //  This is the self-capacitance of that patch: the charge it carries
        //  per volt of its own potential, with every other electrode at zero.
        scalar surfaceCapacitance
        (
            const electromagneticsModel& em,
            const word& patchName
        ) const;

        //- Wedge/axisymmetric scaling to a full 2*pi revolution for a mesh
        static scalar wedgeRevolutionFactor(const fvMesh& m);


    // Access

        inline const word& name() const { return name_; }
        inline const word& unitPatch() const { return unitPatch_; }
        inline const wordList& zeroPatches() const { return zeroPatches_; }
        inline bool solved() const { return solved_; }
        inline bool coupled() const { return coupled_; }
        inline label nDielectrics() const { return nDiel_; }

        //- psi on the gas mesh. Fatal before the first solve.
        const volScalarField& gas() const;

        //- psi on the dielectric meshes
        inline const PtrList<volScalarField>& dielectrics() const
        {
            return dielectrics_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
