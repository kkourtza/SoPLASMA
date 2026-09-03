/*---------------------------------------------------------------------------*\
  File: floatingElectrode.H
  Part of: SoPLASMA

Class
    Foam::floatingElectrode

Description
    The CONSTRAINT for a conductor whose charge is known and whose potential is
    not. The boundary condition floatingElectrodePotential only CARRIES the
    equipotential value; this class determines it.

    Poisson is LINEAR in V, so the constraint closes in CLOSED FORM -- no
    shooting, no relaxation parameter, no convergence criterion. Superposition
    for a fixed space charge and fixed other-electrode potentials gives

        V = V_rho + V_f * psi

    where V_rho solves the problem with the floating patch held at zero and psi
    solves the homogeneous problem with it at one. Then

        INT D.n dA = Q_rho + V_f * C_self,     C_self = INT eps grad(psi).n dA

    and setting that equal to the known Q(t),

        V_f = (Q(t) - Q_rho) / C_self.

    IMPLEMENTED AS A CORRECTION, not as two solves. Rather than re-solving with
    the patch at zero to get V_rho, the ordinary solve runs with the patch at
    its CURRENT value V_f_old, and the residual charge error is removed
    afterwards:

        dV_f = (Q_target - Q_measured) / C_self
        V_f  += dV_f
        V    += dV_f * psi          in EVERY region

    This is algebraically identical -- L(psi) = 0, so adding dV_f*psi keeps the
    field a solution while shifting the patch value by exactly dV_f -- and it
    leaves the main Poisson solve completely untouched. It is EXACT, not
    iterative: one correction lands on Q_target by construction, because the
    relation Q(V_f) is exactly linear.

    Cost: ONE extra Poisson solve whenever psi must be rebuilt, plus two
    surface integrals per step.

    WHEN psi AND C_self MUST BE REBUILT
    They depend only on the OPERATOR and the mesh, so they are cached -- but
    only when that is provably safe. Under `scheme semiImplicit`, which is the
    DEFAULT, the operator is eps + dt*sigma_electrical and BOTH factors change
    every step, so they are rebuilt every step. Under `scheme explicit` on a
    static mesh the operator is eps alone and the cache holds for the whole run.
    A mesh change (AMR, motion) invalidates it in either case.

    Theory, the three analytic tests and the expected physical behaviour:
    docs/models/poisson_equation/floating-electrode.md

SourceFiles
    floatingElectrode.C

\*---------------------------------------------------------------------------*/

#ifndef floatingElectrode_H
#define floatingElectrode_H

#include "fvCFD.H"
#include "unitPotentialField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class electromagneticsModel;

/*---------------------------------------------------------------------------*\
                       Class floatingElectrode Declaration
\*---------------------------------------------------------------------------*/

class floatingElectrode
{
    // Private data

        //- The gas mesh
        const fvMesh& mesh_;

        //- The floating patch, on the gas mesh
        word patchName_;

        //- Patch index on the gas mesh
        label patchi_;

        //- Initial charge Q0 [C], from the boundary condition
        scalar Q0_;

        //- Total charge Q(t) [C]
        scalar Q_;

        //- The floating potential [V], the result
        scalar Vf_;

        //- Self-capacitance [F]
        scalar Cself_;

        //- psi: unit on the floating patch, zero on every known electrode
        autoPtr<unitPotentialField> psi_;

        //- Whether psi/C_self are currently valid
        bool cacheValid_;

        //- The mesh state psi was built against, to detect a mesh change
        label builtAtMeshEvent_;

        //- Whether the start-up report has been printed
        bool reported_;


    // Private Member Functions

        //- Build psi and C_self against the CURRENT operator
        void build(const electromagneticsModel& em, const volScalarField* effEpsGas);

        //- Q on the floating patch from the present field [C]
        scalar measureCharge(const electromagneticsModel& em) const;

        floatingElectrode(const floatingElectrode&) = delete;
        void operator=(const floatingElectrode&) = delete;


public:

    // Constructors

        //- Construct for the single floating patch found on the model
        floatingElectrode
        (
            const fvMesh& gasMesh,
            const electromagneticsModel& em
        );


    //- Destructor
    ~floatingElectrode() = default;


    // Member Functions

        //- Whether the model has any floating conductor at all
        static bool present(const electromagneticsModel& em);

        //- Apply the constraint. Call AFTER the Poisson solve, passing the
        //  effective gas permittivity under semiImplicit (null under
        //  explicit), so psi is built from the same operator.
        void correct
        (
            electromagneticsModel& em,
            const volScalarField* effEpsGas
        );

        //- Add the plasma charge delivered over dt [C]
        void addCharge(const scalar dQ) { Q_ += dQ; }


    // Access

        inline const word& patchName() const { return patchName_; }
        inline scalar Vf() const { return Vf_; }
        inline scalar Q() const { return Q_; }
        inline scalar Q0() const { return Q0_; }
        inline scalar Cself() const { return Cself_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
