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
    // COMPLETELY matrix-free: ONE MATMFFD operator used as BOTH Amat and Pmat,
    // with MatMFFDComputeJacobian registered as the Jacobian function. This is
    // the idiom PETSc's own documentation prescribes for this exact case --
    // "when using a completely matrix-free solver, that is the B matrix is
    // also the same matrix operator".
    //
    // Do NOT go back to SNESSetUseMatrixFree() here. Both of its spellings were
    // measured wrong for this solver on 2026-09-09:
    //
    // * mf_operator = TRUE cost a factor of ~1500. One Newton timestep took
    //   97 s, of which 92 s was 30,020 residual evaluations -- against the ~20
    //   that 3 Newton iterations x 3-7 GMRES iterations justify. That count is
    //   ~10,000 per Newton iteration, exactly the number of unknowns (5 fields
    //   x 2000 cells): the signature of finite-differencing a whole Jacobian
    //   COLUMN BY COLUMN. SNESSetUpMatrices() only takes its
    //   matrix-free-for-both branch when `snes->mf && !snes->mf_operator`;
    //   with mf_operator = TRUE it instead does DMCreateMatrix() for a real
    //   Pmat and fills it with SNESComputeJacobianDefault. PETSc's doc says it
    //   plainly -- mf_operator means "the user provided Pmat will continue to
    //   be used" -- and we never provide one, while the PCSHELL below ignores
    //   Pmat entirely. Pure waste. (Note `snes->mf = mf_operator ? TRUE : mf`,
    //   so passing mf = TRUE as well cannot rescue it.)
    //
    // * (mf_operator = FALSE, mf = TRUE) removed the waste -- 30,020 residual
    //   calls fell to ~20 -- but then GMRES hit DIVERGED_BREAKDOWN at 30
    //   iterations, IDENTICALLY with our PCSHELL and with PCNONE, which proves
    //   the operator itself was degenerate rather than the preconditioner. The
    //   cause is that nothing was refreshing the MFFD differencing BASE vector,
    //   so J*v was being differenced about a stale/zero state. Assembling the
    //   MATMFFD matrix is what pulls the current x from the SNES object, and
    //   MatMFFDComputeJacobian exists precisely to do that -- all it does is
    //   call MatAssemblyBegin/End on the operator.
    {
        Mat Jmf;
        MatCreateSNESMF(snes, &Jmf);
        MatSetFromOptions(Jmf); // so -mat_mffd_type is honoured
        SNESSetJacobian(snes, Jmf, Jmf, MatMFFDComputeJacobian, nullptr);
        MatDestroy(&Jmf);
    }

    {
        KSP ksp;
        SNESGetKSP(snes, &ksp);
        PC pc;
        KSPGetPC(ksp, &pc);

        if (pcCallback)
        {
            PCSetType(pc, PCSHELL);
            PCShellSetContext(pc, &pcCtx);
            PCShellSetApply(pc, PCApplyShell);
            PCShellSetName(pc, "physicsBasedPicard");
        }
        else
        {
            // Amat and Pmat are now BOTH the matrix-free operator, which no
            // factorisation-based default PC can touch. Unpreconditioned is
            // the only coherent choice here, and it must be explicit.
            PCSetType(pc, PCNONE);
        }
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
