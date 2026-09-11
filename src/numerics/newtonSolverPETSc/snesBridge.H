// Plain-types boundary between OpenFOAM (fvCFD.H) and PETSc (petscsnes.h).
// Neither header may be visible in the same translation unit as the other:
// PETSc's petscmath.h and OpenFOAM's own math overloads collide (ambiguous
// pow(), traced into OpenFOAM's error.H) once both are in scope anywhere in
// one TU. This header uses only int/double/void* so it can be included by
// both sides without pulling either heavy header across the boundary.
#ifndef snesBridge_H
#define snesBridge_H

// PETSc init/finalize also live here -- petscsys.h (pulled in by
// PetscInitialize's real declaration) is part of the same header chain
// that causes the pow() ambiguity, so even these two calls must stay on
// the PETSc-only side of the boundary.
void initPetsc(int* argc, char*** argv);

// No-args overload, for a caller (e.g. a runtime-loaded library invoked
// from deep inside soPlasmaFoam's own object construction) that has no
// access to the original argc/argv. Actual PETSc configuration in this
// project goes through setPetscOptions() below (NOT the PETSC_OPTIONS
// environment variable -- see its comment for why that doesn't work here),
// specifically to avoid argList rejecting unrecognised command-line flags
// -- so nothing is lost by not stripping recognised flags out of a real
// argv here.
void initPetsc();

void finalizePetsc();

// Inserts options into PETSc's global options database directly, callable
// any number of times (idempotent per call). Needed because
// PetscInitialize() only reads the PETSC_OPTIONS environment variable
// ONCE, at the moment it runs -- which for this project is inside
// snesNewtonSolver's CONSTRUCTOR, well before the first solveOuterStep()
// call that used to set that env var. Every rtol/maxIt (and any other)
// option set that way was silently never seen by PETSc -- a real bug
// found 2026-09-09 while investigating an apparently-unbounded KSP
// iteration count (PETSc's own default max_it, 10000, was in effect the
// entire time instead of anything this project configured).
void setPetscOptions(const char* options);

typedef void (*ResidualCallback)
(
    int n,
    const double* x,
    double* F,
    void* userData
);

// Physics-based (Knoll & Keyes) preconditioner action: given a right-hand
// side b, compute an APPROXIMATE J^-1*b into y -- in practice one native
// (non-matrix-free) linear solve of the problem's own linearized operator,
// reusing the same physics assembly the residual itself uses. Unpreconditioned
// matrix-free JFNK on this equation was empirically confirmed unstable in
// parallel (converges to a wrong root while still satisfying the requested
// relative residual drop -- see the design doc); this is the planned fix,
// not a bypass of it.
typedef void (*PCApplyCallback)
(
    int n,
    const double* b,
    double* y,
    void* userData
);

// Solves via PETSc's matrix-free SNES. xInOut is the initial guess on
// entry, the converged solution on exit (nLocal entries, this rank's
// local cells only -- global layout follows PETSc's own VecCreateMPI
// convention from the per-rank local sizes, which matches how OpenFOAM's
// mesh decomposition already partitions cells by rank).
// pcCallback/pcUserData may be nullptr for no preconditioning (PCNONE).
// Returns the SNESConvergedReason (positive = converged) and writes the
// iteration count to *itsOut.
// Sentinel meaning "this component has NO lower bound". Chosen rather than
// exposing PETSC_NINFINITY, because the OpenFOAM side of this boundary must
// not see a PETSc header (see the note at the top of this file).
extern const double snesNoBound;

// lowerBounds may be nullptr, in which case an ORDINARY unconstrained Newton
// solve is done (SNESNEWTONLS). When it is supplied, the solve becomes a
// BOUNDED one (SNESVINEWTONRSLS) enforcing lowerBounds[i] <= x[i], with no
// upper bound. Components set to snesNoBound are left free.
//
// WHY BOUNDS AT ALL: species number densities must stay positive. That was
// previously enforced by a clamp applied AFTER the solve, i.e. OUTSIDE the
// equations, which makes F(u)=0 unreachable -- no Newton method can converge
// to a state that its own residual does not call a root. Enforcing the same
// floor as a CONSTRAINT inside the solve is the principled version, and it is
// what makes a cold start possible at all (see docs/design/newton-outer-solver-design.md,
// defect A).
// ASSEMBLED PRECONDITIONING MATRIX, in coordinate (COO) form.
//
// When cooN > 0 a real MATAIJ Pmat is built from these triplets and handed to
// SNESSetJacobian as the SECOND matrix, with the matrix-free operator staying
// as Amat. That combination is the ONLY one that lets PCFIELDSPLIT work here:
// PETSc extracts field sub-blocks from Pmat, and it CANNOT extract them from a
// matrix-free (MATSHELL) operator -- it silently degrades every split to
// `PC type: none` instead of failing, which is why this exists.
//
// Duplicate (row, col) pairs are SUMMED, which is what lets separate physics
// blocks and the cross-coupling be appended independently.
//
// The matrix is assembled ONCE per outer step and held fixed for the whole
// SNES solve -- a deliberately lagged preconditioner (Knoll & Keyes 3.1/5.1).
// The matrix-free Amat still tracks the current iterate exactly, so this costs
// nothing in the converged answer.
struct SnesPmatCOO
{
    int n = 0;                      // number of triplets
    const int* rows = nullptr;
    const int* cols = nullptr;
    const double* vals = nullptr;

    //- RAGGED TAIL LENGTH: rows beyond nFieldBlocks*nCells that belong to the
    //  PHI split only. These are the extra Poisson-only regions (dielectric,
    //  far-field), where the potential is solved but no species exists.
    //
    //  It has to be carried here because the split index sets are built from
    //  nLocal/nFieldBlocks, and that division is WRONG the moment the layout
    //  stops being a clean nFields x nCells rectangle. ZERO restores the exact
    //  previous arithmetic and the exact previous ISCreateStride calls.
    int nTail = 0;

    //- WHICH BLOCK THE SCHUR COMPLEMENT ACTS ON. PETSc's fieldsplit eliminates
    //  split 0 and forms the Schur complement on split 1, and the choice is
    //  made purely by REGISTRATION ORDER -- there is no runtime option for it.
    //
    //  false (default, historic): register phi first, so
    //      S_t = A_tt - A_tf inv(A_ff) A_ft   acts on the TRANSPORT block.
    //  true: register transport first, so
    //      S_f = A_ff - A_ft inv(A_tt) A_tf   acts on the PHI block.
    //
    //  The second is worth having because S_f is a SCALAR ELLIPTIC operator on
    //  one field, which can be assembled explicitly and handed to AMG, whereas
    //  S_t is (nSpecies+1) x nCells and can only be APPLIED -- which is why it
    //  currently runs on ILU(0) of A11 and grinds. See
    //  docs/design/schur-semiimplicit-poisson-preconditioner.md.
    bool schurOnPhi = false;

    //- OPTIONAL explicitly assembled approximation of the phi Schur
    //  complement, supplied to PETSc as PC_FIELDSPLIT_SCHUR_PRE_USER.
    //
    //  The point of the exercise: the SEMI-IMPLICIT POISSON operator
    //  div((eps + dt*sigma) grad .) IS S_f to leading order in dt -- an
    //  identity, derived in the design note above, because
    //  A_ft inv(A_tt) A_tf ~ q * dt * (-div(Z mu n grad .)) = -dt div(sigma grad .)
    //  and sigma = q Z mu n is exactly the conductivity. So the operator that
    //  is INCONSISTENT inside a Newton residual (experiment log section 23) is
    //  precisely the right object in the preconditioner.
    //
    //  Indices are LOCAL TO THE PHI SPLIT (0 .. nCells-1), not global rows.
    //  nSchur == 0 means "not supplied", and PETSc's default (a11) is used.
    //  Only meaningful when schurOnPhi is true.
    int nSchur = 0;
    const int* schurRows = nullptr;
    const int* schurCols = nullptr;
    const double* schurVals = nullptr;
};

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
);

#endif
