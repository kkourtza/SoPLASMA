/*---------------------------------------------------------------------------*\
Application
    testSnesJFNK2Field

Description
    PROOF OF CONCEPT, step 2: validates the BLOCK/multi-field DOF layout
    and a multi-equation, Gauss-Seidel-style physics-based preconditioner,
    BEFORE wiring in the real Poisson+species+energy physics
    (doc/newton-outer-solver-design.md). Deliberately kept separate from
    the real system's complexity (Scharfetter-Gummel schemes, chemistry
    sources, LMEA transport lookups) so a bug in the BLOCK machinery itself
    cannot be confused with a bug in the real equations' physics -- the
    same separation of concerns that made the single-equation .residual()
    parallel defect (see testSnesJFNK.C's header) tractable to isolate.

    Two GENUINELY coupled nonlinear equations on the same reused mesh:
        -D*lap(phi)  + phi*phi2 = f1(x)
        -D*lap(phi2) + phi^2    = f2(x)
    with manufactured solution phi_exact = sin(pi x/L),
    phi2_exact = 2*sin(pi x/L) (both zero at the Dirichlet boundaries,
    phi2 scaled differently from phi so a DOF-layout mixup between the two
    fields' entries shows up as a real, non-degenerate error, not a
    coincidental match).

    Residual reported to SNES (per field, concatenated as
    [phi block][phi2 block] per rank -- the block DOF layout under test):
        R1(phi,phi2) = -D*lap(phi)  + phi*phi2 - f1
        R2(phi,phi2) = -D*lap(phi2) + phi^2    - f2
    computed via EXPLICIT fvc:: evaluation, per the fix in testSnesJFNK.C
    (fvm::+.residual() is NOT used for the outer residual anywhere in this
    project -- see doc/newton-outer-solver-design.md and the
    fvmatrix-residual-broken-in-parallel-use-fvc-instead memory).

    Preconditioner: ONE block Gauss-Seidel sweep, mirroring the real
    system's actual outer sequence (each equation solved once, in a fixed
    order, using the LATEST available values of the others -- exactly how
    soPlasmaFoam.C sequences Poisson -> species -> energy today). Given a
    right-hand side (b1,b2):
        1. solve  -D*lap(dphi)  + phi2*dphi        == b1   for dphi
        2. solve  -D*lap(dphi2)                    == b2 - 2*phi*dphi
           for dphi2, using the JUST-COMPUTED dphi (Gauss-Seidel, not
           Jacobi: the second solve sees the first's result within the
           SAME preconditioner application, matching how the real
           sequence lets species see the just-updated potential).
    Both linearizations use CURRENT (frozen) phi/phi2 as coefficients --
    Picard-style, not the exact Newton Jacobian -- exactly the "reuse the
    existing corrector sweep as a PCSHELL" plan.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "snesBridge.H"

using namespace Foam;

namespace {

struct ResidualContext
{
    volScalarField* phi;
    volScalarField* phi2;
    const volScalarField* fSource1;
    const volScalarField* fSource2;
    const dimensionedScalar* D;
    label nCellsLocal;
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
    volScalarField& phi2 = *ctx.phi2;
    const volScalarField& fSource1 = *ctx.fSource1;
    const volScalarField& fSource2 = *ctx.fSource2;
    const dimensionedScalar& D = *ctx.D;
    const label nc = ctx.nCellsLocal;

    // Block layout under test: this rank's n = 2*nc entries are
    // [phi block (nc)][phi2 block (nc)], contiguous, matching how
    // VecCreateMPI lays out nLocal per rank -- the flat-concatenation
    // approach the design doc's DOF-layout item resolves to in practice
    // (no VecNest needed: FormFunction/PCApply own the field semantics).
    {
        scalarField& phiField = phi.primitiveFieldRef();
        scalarField& phi2Field = phi2.primitiveFieldRef();
        for (label celli = 0; celli < nc; ++celli)
        {
            phiField[celli] = x[celli];
            phi2Field[celli] = x[nc + celli];
        }
        phi.correctBoundaryConditions();
        phi2.correctBoundaryConditions();
    }

    tmp<volScalarField> tLap1 = fvc::laplacian(phi);
    tmp<volScalarField> tLap2 = fvc::laplacian(phi2);
    const scalarField& lap1 = tLap1().primitiveField();
    const scalarField& lap2 = tLap2().primitiveField();
    const scalarField& phiField = phi.primitiveField();
    const scalarField& phi2Field = phi2.primitiveField();
    const scalarField& f1Field = fSource1.primitiveField();
    const scalarField& f2Field = fSource2.primitiveField();

    for (label celli = 0; celli < nc; ++celli)
    {
        F[celli] =
            -D.value()*lap1[celli]
          + phiField[celli]*phi2Field[celli]
          - f1Field[celli];

        F[nc + celli] =
            -D.value()*lap2[celli]
          + phiField[celli]*phiField[celli]
          - f2Field[celli];
    }
}

struct PCContext
{
    volScalarField* phi;
    volScalarField* phi2;
    volScalarField* dphi;
    volScalarField* dphi2;
    const dimensionedScalar* D;
    label nCellsLocal;
};

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
    volScalarField& phi2 = *ctx.phi2;
    volScalarField& dphi = *ctx.dphi;
    volScalarField& dphi2 = *ctx.dphi2;
    const dimensionedScalar& D = *ctx.D;
    const label nc = ctx.nCellsLocal;

    volScalarField b1Field
    (
        IOobject("b1Field", phi.time().timeName(), phi.mesh(),
                 IOobject::NO_READ, IOobject::NO_WRITE),
        phi.mesh(), dimensionedScalar(dimless, Zero)
    );
    volScalarField b2Field
    (
        IOobject("b2Field", phi.time().timeName(), phi.mesh(),
                 IOobject::NO_READ, IOobject::NO_WRITE),
        phi.mesh(), dimensionedScalar(dimless, Zero)
    );
    {
        scalarField& b1 = b1Field.primitiveFieldRef();
        scalarField& b2 = b2Field.primitiveFieldRef();
        for (label celli = 0; celli < nc; ++celli)
        {
            b1[celli] = b[celli];
            b2[celli] = b[nc + celli];
        }
    }

    // Step 1: dphi from the phi-equation's own linearization, freezing
    // phi2 (dropping the phi*dphi2 cross term -- the lower-triangular
    // block-Gauss-Seidel approximation).
    dphi.primitiveFieldRef() = Zero;
    dphi.correctBoundaryConditions();
    {
        fvScalarMatrix eqn
        (
            -D*fvm::laplacian(dphi) + fvm::Sp(phi2, dphi)
         ==
            b1Field
        );
        eqn.solve();
    }

    // Step 2: dphi2 from the phi2-equation's own linearization, now using
    // the JUST-COMPUTED dphi to fold in the coupling term 2*phi*dphi --
    // Gauss-Seidel, not Jacobi: this is what makes it "one sweep" in the
    // same sense as the real system's sequential Poisson->species->energy
    // pass, not two independent decoupled solves.
    dphi2.primitiveFieldRef() = Zero;
    dphi2.correctBoundaryConditions();
    {
        volScalarField rhs2
        (
            IOobject("rhs2", phi.time().timeName(), phi.mesh(),
                     IOobject::NO_READ, IOobject::NO_WRITE),
            b2Field - 2.0*phi*dphi
        );
        fvScalarMatrix eqn
        (
            -D*fvm::laplacian(dphi2)
         ==
            rhs2
        );
        eqn.solve();
    }

    const scalarField& dphiField = dphi.primitiveField();
    const scalarField& dphi2Field = dphi2.primitiveField();
    for (label celli = 0; celli < nc; ++celli)
    {
        y[celli] = dphiField[celli];
        y[nc + celli] = dphi2Field[celli];
    }
}

} // namespace


int main(int argc, char *argv[])
{
    initPetsc(&argc, &argv);

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    );
    volScalarField phi2
    (
        IOobject("phi2", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    );

    const scalar L = 0.01;
    const dimensionedScalar D("D", dimensionSet(0, 2, 0, 0, 0, 0, 0), 1.0);
    const scalar pi = Foam::constant::mathematical::pi;

    const scalarField xCoordRaw(mesh.C().component(0)().primitiveField());
    scalarField phiExactRaw(xCoordRaw.size());
    scalarField phi2ExactRaw(xCoordRaw.size());
    forAll(phiExactRaw, celli)
    {
        phiExactRaw[celli] = Foam::sin((pi/L)*xCoordRaw[celli]);
        phi2ExactRaw[celli] = 2.0*Foam::sin((pi/L)*xCoordRaw[celli]);
    }
    const scalarField f1Raw
    (
        D.value()*sqr(pi/L)*phiExactRaw + phiExactRaw*phi2ExactRaw
    );
    const scalarField f2Raw
    (
        D.value()*sqr(pi/L)*phi2ExactRaw + phiExactRaw*phiExactRaw
    );

    volScalarField phiExact
    (
        IOobject("phiExact", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimless, Zero)
    );
    phiExact.primitiveFieldRef() = phiExactRaw;

    volScalarField phi2Exact
    (
        IOobject("phi2Exact", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimless, Zero)
    );
    phi2Exact.primitiveFieldRef() = phi2ExactRaw;

    volScalarField fSource1
    (
        IOobject("fSource1", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimless, Zero)
    );
    fSource1.primitiveFieldRef() = f1Raw;

    volScalarField fSource2
    (
        IOobject("fSource2", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimless, Zero)
    );
    fSource2.primitiveFieldRef() = f2Raw;

    const label nCellsLocal = mesh.nCells();

    ResidualContext ctx{&phi, &phi2, &fSource1, &fSource2, &D, nCellsLocal};

    volScalarField dphi(IOobject("dphi", runTime.timeName(), mesh,
                                  IOobject::NO_READ, IOobject::NO_WRITE), phi);
    dphi.primitiveFieldRef() = Zero;
    dphi.correctBoundaryConditions();

    volScalarField dphi2(IOobject("dphi2", runTime.timeName(), mesh,
                                   IOobject::NO_READ, IOobject::NO_WRITE), phi2);
    dphi2.primitiveFieldRef() = Zero;
    dphi2.correctBoundaryConditions();

    PCContext pcCtx{&phi, &phi2, &dphi, &dphi2, &D, nCellsLocal};

    const label nLocalTotal = 2*nCellsLocal;
    List<double> xBuf(nLocalTotal, Zero);
    {
        const scalarField& phiField = phi.primitiveField();
        const scalarField& phi2Field = phi2.primitiveField();
        for (label celli = 0; celli < nCellsLocal; ++celli)
        {
            xBuf[celli] = phiField[celli];
            xBuf[nCellsLocal + celli] = phi2Field[celli];
        }
    }

    int its = 0;
    const int reason = solveWithSNES
    (
        int(nLocalTotal),
        xBuf.data(),
        &residualCallback,
        &ctx,
        &pcApplyCallback,
        &pcCtx,
        &its
    );

    {
        scalarField& phiField = phi.primitiveFieldRef();
        scalarField& phi2Field = phi2.primitiveFieldRef();
        for (label celli = 0; celli < nCellsLocal; ++celli)
        {
            phiField[celli] = xBuf[celli];
            phi2Field[celli] = xBuf[nCellsLocal + celli];
        }
        phi.correctBoundaryConditions();
        phi2.correctBoundaryConditions();
    }

    const scalarField diff1(phi.primitiveField() - phiExact.primitiveField());
    const scalarField diff2(phi2.primitiveField() - phi2Exact.primitiveField());
    scalar sumSq1 = sum(diff1*diff1);
    scalar sumSq2 = sum(diff2*diff2);
    scalar nTotal = scalar(mesh.nCells());
    reduce(sumSq1, sumOp<scalar>());
    reduce(sumSq2, sumOp<scalar>());
    reduce(nTotal, sumOp<scalar>());
    const scalar l2Error1 = Foam::sqrt(sumSq1/nTotal);
    const scalar l2Error2 = Foam::sqrt(sumSq2/nTotal);

    Info<< "testSnesJFNK2Field: SNESConvergedReason = " << reason
        << " (positive = converged), iterations = " << its << nl
        << "testSnesJFNK2Field: L2 error phi  = " << l2Error1 << nl
        << "testSnesJFNK2Field: L2 error phi2 = " << l2Error2 << nl
        << "testSnesJFNK2Field: "
        << (reason > 0 ? "PASS" : "FAIL") << endl;

    phi.write();
    phi2.write();

    finalizePetsc();

    return reason > 0 ? 0 : 1;
}

// ************************************************************************* //
