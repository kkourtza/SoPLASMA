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
int solveWithSNES
(
    int nLocal,
    double* xInOut,
    ResidualCallback callback,
    void* userData,
    PCApplyCallback pcCallback,
    void* pcUserData,
    int* itsOut
);

#endif
