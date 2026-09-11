/*---------------------------------------------------------------------------*\
Application
    testSnesJFNK

Description
    PROOF OF CONCEPT for the JFNK outer-solver design
    (docs/design/newton-outer-solver-design.md), 2026-09-09. Not yet the real
    3-field system -- validates the CORE mechanism first: can a residual
    evaluated on an OpenFOAM field be driven to zero by PETSc's
    matrix-free SNES, correctly in parallel, with a physics-based
    (Knoll & Keyes) preconditioner built from the SAME equation via a
    native linear solve?

    The outer residual (FormFunction, below) is computed via EXPLICIT
    fvc:: evaluation, NOT via fvm::+fvMatrix::residual(). The latter was
    tried first (matching the design doc's original premise of reusing
    every equation's existing fvScalarMatrix "for free" via .residual())
    and was found, via an independent non-OpenFOAM finite-difference
    check, to silently and massively misreport in parallel: a field state
    whose true residual was ~6.5e9 was reported as ~1e-12 by BOTH PETSc's
    FormFunction and a standalone direct .residual() call on the identical
    matrix -- ruled out as a PETSc, preconditioning, or MPI-communicator
    issue, since none of those changed the (bit-identical) wrong outcome,
    while a completely native .solve() on that same matrix correctly
    detected the imbalance. fvc:: explicit evaluation does not have this
    defect (confirmed bit-identical between serial and 4-rank parallel).
    Practical implication for the real 3-field system: its outer residual
    must be written explicitly (fvc::), not harvested via .residual();
    the preconditioner, which uses fvm::+.solve() rather than .residual(),
    is unaffected and can still reuse existing equation assembly as
    planned.

    Manufactured solution: phi_exact(x) = sin(pi x / L) on the reused
    Grubert mesh (L = 0.01 m, the gap length), with homogeneous Dirichlet
    BCs at cathode/anode (both zero, matching sin(0) = sin(pi) = 0).

    PDE:      -D lap(phi) + phi^2 = f(x)
    Chosen so phi_exact solves it exactly:
              f(x) = D (pi/L)^2 sin(pi x/L) + sin(pi x/L)^2

    Residual reported to SNES:
              R(phi) = -D lap(phi) + phi^2 - f

    R(phi_exact) = 0 by construction; starting from phi = 0 everywhere
    (a genuinely bad initial guess, not "close to the answer") is a real
    test of whether SNES actually converges, not just holds still.

    This file includes fvCFD.H ONLY -- it never sees a PETSc header. All
    PETSc calls live in snesBridge.C, reached through the plain-types
    interface in snesBridge.H. See that header for why the split is
    required (a real ambiguous-overload clash between PETSc's and
    OpenFOAM's math headers when both are visible in one translation unit).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "snesBridge.H"

using namespace Foam;

namespace {

struct ResidualContext
{
    fvMesh* mesh;
    volScalarField* phi;
    const volScalarField* fSource;
    const dimensionedScalar* D;
};

void residualCallback
(
    int n,
    const double* x,
    double* F,
    void* userDataVoid
)
{
    ResidualContext& ctx = *static_cast<ResidualContext*>(userDataVoid);
    volScalarField& phi = *ctx.phi;
    const volScalarField& fSource = *ctx.fSource;
    const dimensionedScalar& D = *ctx.D;

    // Unpack the trial state into the OpenFOAM field, THEN exchange
    // processor-boundary (halo) values before assembling anything -- per
    // the parallel-correctness checklist: a residual assembled without a
    // fresh halo exchange is wrong at processor boundaries, not just
    // stale, and SNES calls this function repeatedly (real Newton steps
    // AND internal finite-differenced Jacobian-vector products), so this
    // must happen on every call.
    {
        scalarField& phiField = phi.primitiveFieldRef();
        for (int celli = 0; celli < n; ++celli)
        {
            phiField[celli] = x[celli];
        }
        phi.correctBoundaryConditions();
    }

    // R computed via EXPLICIT fvc:: evaluation, not fvm::+.residual().
    // fvMatrix::residual(), confirmed via an independent (non-OpenFOAM,
    // raw finite-difference-in-Python) calculation, can silently and
    // massively misreport when called this way on a decomposed/parallel
    // field: a state whose TRUE residual was ~6.5e9 was reported as
    // ~1e-12 by BOTH PETSc's FormFunction and a standalone direct call to
    // .residual() on the identical matrix -- yet .solve() on that SAME
    // matrix correctly detected the imbalance and needed 78 PCG iterations
    // to fix it. That rules out PETSc, preconditioning, and MPI-communicator
    // sharing as the cause; the defect is specific to .residual() itself
    // in this usage. fvc:: explicit evaluation was the ORIGINAL (v1)
    // approach here and was already confirmed bit-identical between serial
    // and 4-rank parallel, so it sidesteps the defect entirely. The
    // preconditioner below still uses fvm::+.solve() (not .residual()),
    // which is unaffected.
    tmp<volScalarField> tLap = fvc::laplacian(phi);
    const scalarField& lapField = tLap().primitiveField();
    const scalarField& phiField = phi.primitiveField();
    const scalarField& fSourceField = fSource.primitiveField();

    for (int celli = 0; celli < n; ++celli)
    {
        F[celli] =
            -D.value()*lapField[celli]
          + phiField[celli]*phiField[celli]
          - fSourceField[celli];
    }
}

struct PCContext
{
    fvMesh* mesh;
    volScalarField* phi;      // current outer iterate, read-only here
    volScalarField* dphi;     // scratch correction field, reused every call
    const dimensionedScalar* D;
};

// Physics-based (Knoll & Keyes) preconditioner: approximate J^-1*b by ONE
// native (non-matrix-free) linear solve of the true Jacobian of
// -D*lap(phi) + phi^2 -- i.e. -D*lap(dphi) + 2*phi*dphi == b -- using
// OpenFOAM's own parallel-correct linear solver (the same machinery every
// other OpenFOAM solve() call uses, confirmed correct in parallel via the
// native-Picard isolation test). Unpreconditioned matrix-free JFNK on this
// problem converged to a WRONG root in parallel while still satisfying the
// requested relative residual drop -- see docs/design/newton-outer-solver-design.md.
// dphi carries the SAME (homogeneous) BC types as phi: a correction must
// vanish wherever phi itself is fixed by a Dirichlet BC.
void pcApplyCallback
(
    int n,
    const double* b,
    double* y,
    void* userDataVoid
)
{
    PCContext& ctx = *static_cast<PCContext*>(userDataVoid);
    volScalarField& phi = *ctx.phi;
    volScalarField& dphi = *ctx.dphi;
    const dimensionedScalar& D = *ctx.D;

    volScalarField bField
    (
        IOobject("bField", phi.time().timeName(), phi.mesh(),
                 IOobject::NO_READ, IOobject::NO_WRITE),
        phi.mesh(),
        dimensionedScalar(dimless, Zero)
    );
    {
        scalarField& bf = bField.primitiveFieldRef();
        for (int celli = 0; celli < n; ++celli)
        {
            bf[celli] = b[celli];
        }
    }

    dphi.primitiveFieldRef() = Zero;
    dphi.correctBoundaryConditions();

    fvScalarMatrix eqn
    (
        -D*fvm::laplacian(dphi) + fvm::Sp(2.0*phi, dphi)
     ==
        bField
    );
    eqn.solve();

    const scalarField& dphiField = dphi.primitiveField();
    for (int celli = 0; celli < n; ++celli)
    {
        y[celli] = dphiField[celli];
    }
}

} // namespace


int main(int argc, char *argv[])
{
    // PETSc first, so it strips its own command-line flags (-snes_mf,
    // -snes_monitor, -snes_rtol, ...) before OpenFOAM's argList parses
    // whatever remains -- the standard pattern for combining PETSc's CLI
    // handling with another framework's.
    initPetsc(&argc, &argv);

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField phi
    (
        IOobject
        (
            "phi", runTime.timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        ),
        mesh
    );

    const scalar L = 0.01;      // gap length, matches the reused mesh
    // Length^2 dimensions: with phi/fSource dimensionless, lap(phi) carries
    // 1/length^2, so D needs length^2 to cancel back to dimensionless in
    // the preconditioner's fvm::-assembled linear solve below -- see
    // pcApplyCallback's comment. (The outer residual itself is computed
    // via raw fvc:: + scalarField arithmetic, bypassing dimension checking
    // entirely -- see residualCallback's comment.)
    const dimensionedScalar D("D", dimensionSet(0, 2, 0, 0, 0, 0, 0), 1.0);
    const scalar pi = Foam::constant::mathematical::pi;

    // Raw scalarField arithmetic here, not volScalarField: mesh.C() carries
    // LENGTH dimensions, and OpenFOAM's Foam::sin() correctly refuses a
    // dimensioned argument. This manufactured-solution setup is a plain
    // numerical test artifact with no physical dimensions to track, so
    // side-stepping GeometricField's dimension checking is the right tool,
    // not a workaround.
    const scalarField xCoordRaw(mesh.C().component(0)().primitiveField());
    scalarField phiExactRaw(xCoordRaw.size());
    forAll(phiExactRaw, celli)
    {
        phiExactRaw[celli] = Foam::sin((pi/L)*xCoordRaw[celli]);
    }
    const scalarField fSourceRaw
    (
        D.value()*sqr(pi/L)*phiExactRaw + phiExactRaw*phiExactRaw
    );

    volScalarField phiExact
    (
        IOobject("phiExact", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimless, Zero)
    );
    phiExact.primitiveFieldRef() = phiExactRaw;

    volScalarField fSource
    (
        IOobject("fSource", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimless, Zero)
    );
    fSource.primitiveFieldRef() = fSourceRaw;

    const label nCellsLocal = mesh.nCells();

    ResidualContext ctx{&mesh, &phi, &fSource, &D};

    // dphi copies phi as a BC prototype (same patch types/values -- the
    // Dirichlet value is 0 on both, so this gives the correct HOMOGENEOUS
    // BCs for a correction field with no extra bookkeeping), then gets
    // zeroed: a correction must vanish wherever phi itself is pinned.
    volScalarField dphi
    (
        IOobject("dphi", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        phi
    );
    dphi.primitiveFieldRef() = Zero;
    dphi.correctBoundaryConditions();

    PCContext pcCtx{&mesh, &phi, &dphi, &D};

    List<double> xBuf(nCellsLocal, Zero);
    {
        const scalarField& phiField = phi.primitiveField();
        forAll(phiField, celli)
        {
            xBuf[celli] = phiField[celli];
        }
    }

    int its = 0;
    const int reason = solveWithSNES
    (
        int(nCellsLocal),
        xBuf.data(),
        &residualCallback,
        &ctx,
        &pcApplyCallback,
        &pcCtx,
        &its
    );

    {
        scalarField& phiField = phi.primitiveFieldRef();
        forAll(phiField, celli)
        {
            phiField[celli] = xBuf[celli];
        }
        phi.correctBoundaryConditions();
    }

    const scalarField diff(phi.primitiveField() - phiExact.primitiveField());
    scalar sumSq = sum(diff*diff);
    scalar nTotal = scalar(mesh.nCells());
    reduce(sumSq, sumOp<scalar>());
    reduce(nTotal, sumOp<scalar>());
    const scalar l2Error = Foam::sqrt(sumSq/nTotal);

    Info<< "testSnesJFNK: SNESConvergedReason = " << reason
        << " (positive = converged), iterations = " << its << nl
        << "testSnesJFNK: L2 error vs manufactured solution = "
        << l2Error << nl
        << "testSnesJFNK: " << (reason > 0 ? "PASS" : "FAIL") << endl;

    phi.write();

    finalizePetsc();

    return reason > 0 ? 0 : 1;
}

// ************************************************************************* //
