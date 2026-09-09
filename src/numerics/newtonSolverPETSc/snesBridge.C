// PETSc-only side of the bridge -- deliberately does NOT include fvCFD.H
// or any other OpenFOAM header. See snesBridge.H for why.
#include "snesBridge.H"
#include <petscsnes.h>

namespace {

struct CallbackContext
{
    ResidualCallback callback;
    void* userData;
};

struct PCContext
{
    PCApplyCallback callback;
    void* userData;
};

PetscErrorCode FormFunction(SNES /*snes*/, Vec x, Vec F, void* ctxVoid)
{
    CallbackContext& ctx = *static_cast<CallbackContext*>(ctxVoid);

    const PetscScalar* xArr;
    VecGetArrayRead(x, &xArr);

    PetscInt n;
    VecGetLocalSize(x, &n);

    PetscScalar* fArr;
    VecGetArray(F, &fArr);

    ctx.callback
    (
        int(n),
        static_cast<const double*>(xArr),
        static_cast<double*>(fArr),
        ctx.userData
    );

    VecRestoreArrayRead(x, &xArr);
    VecRestoreArray(F, &fArr);

    return 0;
}

PetscErrorCode PCApplyShell(PC pc, Vec b, Vec y)
{
    PCContext* ctx;
    PCShellGetContext(pc, reinterpret_cast<void**>(&ctx));

    const PetscScalar* bArr;
    VecGetArrayRead(b, &bArr);

    PetscInt n;
    VecGetLocalSize(b, &n);

    PetscScalar* yArr;
    VecGetArray(y, &yArr);

    ctx->callback
    (
        int(n),
        static_cast<const double*>(bArr),
        static_cast<double*>(yArr),
        ctx->userData
    );

    VecRestoreArrayRead(b, &bArr);
    VecRestoreArray(y, &yArr);

    return 0;
}

} // namespace

// A PRIVATE duplicate of MPI_COMM_WORLD for PETSc's own exclusive use.
// residualCallback (the FormFunction body) calls back into OpenFOAM code
// that does its OWN non-blocking Pstream communication (fvm::laplacian's
// processor-patch coupling, fvMatrix::residual()'s interface exchange) on
// MPI_COMM_WORLD, nested INSIDE a PETSc-driven call sequence that is ALSO
// doing MPI traffic (Vec norms, KSP reductions) on the same communicator.
// Two libraries' independent non-blocking sends/receives sharing one
// literal communicator can interleave and cross-match messages by tag --
// empirically, exactly this nesting (never present in the working
// native-solve path, which has no PETSc MPI activity to interleave with)
// was the only common factor across every JFNK variant that converged to
// a wrong-but-self-consistent root in parallel while giving the correct
// answer in serial. Giving PETSc a dup'd communicator makes its internal
// traffic use a distinct context, immune to collision with OpenFOAM's own.
MPI_Comm petscComm = MPI_COMM_NULL;

void initPetsc(int* argc, char*** argv)
{
    PetscInitialize(argc, argv, nullptr, nullptr);
    MPI_Comm_dup(MPI_COMM_WORLD, &petscComm);
}

void initPetsc()
{
    int argc = 0;
    char** argv = nullptr;
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    MPI_Comm_dup(MPI_COMM_WORLD, &petscComm);
}

void finalizePetsc()
{
    MPI_Comm_free(&petscComm);
    PetscFinalize();
}

void setPetscOptions(const char* options)
{
    PetscOptionsInsertString(nullptr, options);
}

int solveWithSNES
(
    int nLocal,
    double* xInOut,
    ResidualCallback callback,
    void* userData,
    PCApplyCallback pcCallback,
    void* pcUserData,
    int* itsOut
)
{
    Vec x, F;
    VecCreateMPI(petscComm, nLocal, PETSC_DETERMINE, &x);
    VecDuplicate(x, &F);

    {
        PetscScalar* xArr;
        VecGetArray(x, &xArr);
        for (int i = 0; i < nLocal; ++i)
        {
            xArr[i] = xInOut[i];
        }
        VecRestoreArray(x, &xArr);
    }

    CallbackContext ctx{callback, userData};
    PCContext pcCtx{pcCallback, pcUserData};

    SNES snes;
    SNESCreate(petscComm, &snes);
    SNESSetFunction(snes, F, FormFunction, &ctx);
    SNESSetType(snes, SNESNEWTONLS);
    // Matrix-free: no assembled Jacobian. PETSc finite-differences
    // FormFunction internally for Jacobian-vector products -- the JFNK
    // core idea. Unpreconditioned matrix-free JFNK on this problem was
    // empirically found to converge to a WRONG root in parallel (while
    // still satisfying the requested relative residual drop) -- see the
    // design doc. The physics-based (Knoll & Keyes) preconditioner below
    // is the fix: one native linear solve of the problem's own linearized
    // operator per GMRES preconditioner application.
    //
    // The second argument (mf_operator) tells SNES to ALSO derive its own
    // (dense, finite-differenced) preconditioning matrix -- wasteful and
    // beside the point when we supply a real PC via PCSHELL below, so it's
    // only enabled when no PC callback is given (falling back to whatever
    // PETSc's own default construction provides).
    SNESSetUseMatrixFree(snes, PETSC_TRUE, pcCallback ? PETSC_FALSE : PETSC_TRUE);

    if (pcCallback)
    {
        KSP ksp;
        SNESGetKSP(snes, &ksp);
        PC pc;
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCSHELL);
        PCShellSetContext(pc, &pcCtx);
        PCShellSetApply(pc, PCApplyShell);
        PCShellSetName(pc, "physicsBasedPicard");
    }

    SNESSetFromOptions(snes);

    SNESSolve(snes, nullptr, x);

    SNESConvergedReason reason;
    SNESGetConvergedReason(snes, &reason);
    PetscInt its;
    SNESGetIterationNumber(snes, &its);
    *itsOut = int(its);

    {
        const PetscScalar* xArr;
        VecGetArrayRead(x, &xArr);
        for (int i = 0; i < nLocal; ++i)
        {
            xInOut[i] = xArr[i];
        }
        VecRestoreArrayRead(x, &xArr);
    }

    VecDestroy(&x);
    VecDestroy(&F);
    SNESDestroy(&snes);

    return int(reason);
}
