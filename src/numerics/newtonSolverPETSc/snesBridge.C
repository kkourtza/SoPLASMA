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

const double snesNoBound = -1.0e308;

int solveWithSNES
(
    int nLocal,
    double* xInOut,
    ResidualCallback callback,
    void* userData,
    PCApplyCallback pcCallback,
    void* pcUserData,
    const double* lowerBounds,
    const SnesPmatCOO* pmat,
    int nFieldBlocks,
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

    // TYPE FIRST, THEN THE FUNCTION -- PETSc's canonical order. Setting the
    // type afterwards swaps the solver implementation underneath state that
    // has already been registered, and SNESVISetVariableBounds_VI() assigns
    // snes->vec_func itself (vi.c), so the ordering is load-bearing here in a
    // way it is not for plain NEWTONLS.
    SNESSetType(snes, lowerBounds ? SNESVINEWTONRSLS : SNESNEWTONLS);

    SNESSetFunction(snes, F, FormFunction, &ctx);

    // BOUNDED (variational-inequality) Newton when bounds are supplied.
    // SNESVINEWTONRSLS is the reduced-space active-set method: components at
    // their bound are removed from the Newton system for that iteration,
    // rather than being clamped afterwards.
    Vec xl = nullptr;
    if (lowerBounds)
    {
        // NOT PETSC_INFINITY. That is PETSC_MAX_REAL/4 ~ 4.5e307 -- finite,
        // but close enough to the top of the range that SNESVI's own
        // arithmetic on it (xu - xl is already 9e307) overflows to inf, and
        // OpenFOAM runs with FOAM_SIGFPE enabled, so an overflow TRAPS. That
        // was measured as an immediate SIGFPE at the first bounded solve,
        // 2026-09-09.
        //
        // A far smaller magnitude is both safe and sufficient because
        // everything PETSc sees here is SCALED to O(1) (see sX/sF): 1e30 is
        // effectively unbounded for an O(1) variable. The only cost is that
        // PETSc's `ntruebounds` bookkeeping compares against PETSC_INFINITY
        // exactly, so it will count these components as bounded -- that
        // figure is a heuristic, not a correctness input.
        const PetscScalar bigBound = 1.0e30;

        Vec xu;
        VecDuplicate(x, &xl);
        VecDuplicate(x, &xu);
        VecSet(xu, bigBound);
        {
            PetscScalar* lArr;
            VecGetArray(xl, &lArr);
            for (int i = 0; i < nLocal; ++i)
            {
                lArr[i] = (lowerBounds[i] <= snesNoBound)
                        ? -bigBound
                        : lowerBounds[i];
            }
            VecRestoreArray(xl, &lArr);
        }
        SNESVISetVariableBounds(snes, xl, xu);
        VecDestroy(&xu);
    }
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
    Mat Pmat = nullptr;
    {
        Mat Jmf;
        MatCreateSNESMF(snes, &Jmf);
        MatSetFromOptions(Jmf); // so -mat_mffd_type is honoured

        if (pmat && pmat->n > 0)
        {
            // A REAL assembled Pmat alongside the matrix-free Amat. PETSc's
            // manual is explicit that fieldsplit takes its blocks from Pmat,
            // not Amat, and a MATSHELL cannot supply them.
            MatCreate(petscComm, &Pmat);
            MatSetSizes(Pmat, nLocal, nLocal, PETSC_DETERMINE, PETSC_DETERMINE);
            MatSetType(Pmat, MATAIJ);
            MatSetFromOptions(Pmat);

            // Preallocation: count entries per row from the triplets. Without
            // this, assembly is quadratic in the worst case.
            {
                PetscInt* nnz = nullptr;
                PetscMalloc1(nLocal, &nnz);
                for (int i = 0; i < nLocal; ++i) { nnz[i] = 0; }
                for (int k = 0; k < pmat->n; ++k)
                {
                    const int r = pmat->rows[k];
                    if (r >= 0 && r < nLocal) { ++nnz[r]; }
                }
                MatSeqAIJSetPreallocation(Pmat, 0, nnz);
                MatMPIAIJSetPreallocation(Pmat, 0, nnz, 0, nullptr);
                PetscFree(nnz);
            }

            // ADD_VALUES so duplicate coordinates sum -- that is what lets the
            // diagonal blocks and the cross-coupling be appended separately.
            MatSetOption(Pmat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
            for (int k = 0; k < pmat->n; ++k)
            {
                const PetscInt r = pmat->rows[k];
                const PetscInt c = pmat->cols[k];
                const PetscScalar v = pmat->vals[k];
                MatSetValues(Pmat, 1, &r, 1, &c, &v, ADD_VALUES);
            }
            MatAssemblyBegin(Pmat, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(Pmat, MAT_FINAL_ASSEMBLY);

            // MatMFFDComputeJacobian assembles only the MFFD operator it is
            // given, leaving Pmat untouched -- exactly the lagged-Pmat
            // behaviour we want for one outer step.
            SNESSetJacobian(snes, Jmf, Pmat, MatMFFDComputeJacobian, nullptr);
        }
        else
        {
            SNESSetJacobian(snes, Jmf, Jmf, MatMFFDComputeJacobian, nullptr);
        }

        MatDestroy(&Jmf);
    }

    {
        KSP ksp;
        SNESGetKSP(snes, &ksp);
        PC pc;
        KSPGetPC(ksp, &pc);

        // BOUNDED AND FIELDSPLIT ARE MUTUALLY EXCLUSIVE TODAY, and the
        // reason is structural rather than a bug of ours.
        //
        // SNESVINEWTONRSLS calls SNESVIResetPCandKSP() on every active-set
        // change, which does KSPReset() + KSPResetFromOptions(); PCReset_FieldSplit
        // then DESTROYS the user-supplied index sets (ISDestroy(&ilink->is)).
        // The very next call, PCFieldSplitRestrictIS(), therefore operates on a
        // fieldsplit with no splits left, and PCSetUp fails with PETSc's
        // "Unhandled case, must have at least two fields, not 0".
        //
        // Note the `FromOptions` in that reset: splits defined by OPTIONS are
        // regenerated, but option-defined splits require an INTERLACED block
        // size, and our DOF layout is field-major. The real fix is a
        // DM-based field decomposition (-pc_fieldsplit_dm_splits, the pattern
        // in PETSc's own src/snes/tutorials/ex28.c), which lets PETSc rebuild
        // the splits after any reset. Until that exists, choose explicitly and
        // say so, rather than letting the run die inside PCSetUp.
        if (Pmat && nFieldBlocks > 1 && lowerBounds)
        {
            static bool warned = false;
            if (!warned)
            {
                warned = true;
                PetscPrintf
                (
                    petscComm,
                    "snesNewtonSolver: BOUNDED solve requested, so the "
                    "assembled-Pmat PCFIELDSPLIT preconditioner is disabled "
                    "and the physics-based shell is used instead.\n"
                    "    They are mutually exclusive today: the bounded "
                    "solver resets the KSP on every active-set change, which "
                    "destroys fieldsplit's index sets. A DM-based field "
                    "decomposition is the fix.\n"
                    "    Set `bounded false` in the newtonSolver dict to use "
                    "fieldsplit, at the cost of the density floor being a "
                    "post-solve clamp again.\n"
                );
            }
        }

        if (Pmat && nFieldBlocks > 1 && !lowerBounds)
        {
            // PCFIELDSPLIT on the ASSEMBLED Pmat.
            //
            // Two splits, because Schur requires EXACTLY two
            // (fieldsplit.c: "To use Schur complement preconditioner you must
            // have exactly 2 fields"): the elliptic potential against ALL the
            // transport fields. The Schur complement is the only fieldsplit
            // type that uses the off-diagonal coupling block, which is the
            // whole reason the Pmat was assembled.
            //
            // PCFieldSplitSetIS, NOT PCFieldSplitSetFields: our DOF layout is
            // FIELD-MAJOR (all of field 0, then all of field 1, ...), not
            // interlaced, so the block-size-based helper does not describe it.
            const PetscInt nc = nLocal/nFieldBlocks;

            IS isPhi, isTransport;
            ISCreateStride(petscComm, nc, 0, 1, &isPhi);
            ISCreateStride(petscComm, nLocal - nc, nc, 1, &isTransport);

            PCSetType(pc, PCFIELDSPLIT);
            PCFieldSplitSetIS(pc, "phi", isPhi);
            PCFieldSplitSetIS(pc, "transport", isTransport);
            PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);

            ISDestroy(&isPhi);
            ISDestroy(&isTransport);
        }
        else if (pcCallback)
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

    // SNES FIRST, THEN the vectors. SNES holds references to the solution and
    // function vectors, and SNESVI additionally keeps active-set work vectors
    // and index sets derived from the bounds, resetting its KSP whenever the
    // active set changes (virs.c). Freeing the vectors first left it operating
    // on released memory, which showed up as
    // "double free or corruption (out)" inside SNESDestroy -> KSPReset_GMRES.
    // Measured 2026-09-09, and only on the bounded path -- with `bounded
    // false` the same run is clean.
    SNESDestroy(&snes);
    if (Pmat) MatDestroy(&Pmat);
    VecDestroy(&x);
    VecDestroy(&F);
    if (xl) VecDestroy(&xl);

    return int(reason);
}
