/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::snesNewtonSolver

Description
    Concrete `plasmaNewtonSolver` implementation using PETSc's matrix-free
    SNES -- see doc/newton-outer-solver-design.md for the full design and
    the proofs of concept this is built on.

    NEVER includes any PETSc header directly -- petscmath.h and OpenFOAM's
    own math operator overloads are ambiguous when both are visible in one
    translation unit (a real build obstacle found and fixed in the proof of
    concept). All PETSc calls live in `snesBridge.C`, reached through the
    plain-types interface in `snesBridge.H`, exactly as the proof of concept
    established.

    Physics-based (Knoll & Keyes) preconditioned matrix-free JFNK: one
    block Gauss-Seidel sweep (Poisson, then species, then energy -- the
    SAME order and operators the real Picard sequence uses) as the
    preconditioner, each equation solving a linearised CORRECTION against
    a synthetic GMRES right-hand side. Bare/unpreconditioned JFNK was
    tried first (matching how the single-equation proof of concept was
    validated) and found, 2026-09-09, to reliably hit
    SNES_DIVERGED_LINEAR_SOLVE on the real system's first genuine Newton
    iteration -- expected, matching the design doc's diagnosis that this
    coupled system is inherently stiff (fast electron transport against
    Poisson) independent of the later negative-differential-resistance
    region.

    Field block, one instance per rank, flat-concatenated in this fixed
    order: [ePotential][n_0]...[n_{nSpecies-1}][nEps_e if LMEA]. The outer
    residual for EVERY block is computed via EXPLICIT fvc:: evaluation,
    never fvm::+.residual() -- confirmed broken in parallel for that
    specific usage, independent of PETSc (see
    [[fvmatrix-residual-broken-in-parallel-use-fvc-instead]]).

SourceFiles
    snesNewtonSolver.C

\*---------------------------------------------------------------------------*/

#ifndef snesNewtonSolver_H
#define snesNewtonSolver_H

#include "plasmaNewtonSolver.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class snesNewtonSolver Declaration
\*---------------------------------------------------------------------------*/

class snesNewtonSolver
:
    public plasmaNewtonSolver
{
    const fvMesh& mesh_;

    //- `-snes_rtol`/`-snes_max_it`-equivalent, read from the dict rather
    //  than expecting PETSc command-line flags (which argList would reject
    //  as unrecognised -- a real obstacle found in the proof of concept).
    scalar rtol_;
    label maxIt_;

    //- Relative error to which F(u) can be EVALUATED, setting the matrix-free
    //  differencing step (PETSc `-mat_mffd_err`). NEGATIVE (the default) means
    //  do not set it at all and leave PETSc's own choice in place.
    //
    //  Knoll & Keyes section 2.3.1 gives the rule this exists for: when F can
    //  only be evaluated to eps_rel, the differencing parameter must be built
    //  from eps_rel rather than eps_mach. It is retained as an escape hatch
    //  for a genuinely noisy F -- for instance a case whose stiff chemistry
    //  ODE actually fires, which grubert2009's does not (its own report says
    //  "0 of 2000 cells integrated").
    //
    //  DEFAULT 1e-5, and it is a WORKAROUND, not a physical constant. PETSc's
    //  own default (PETSC_SQRT_MACHINE_EPSILON ~ 1.49e-8) does NOT work here.
    //  Measured 2026-09-09 over 20 steps from the 2e-9 Picard restart:
    //  1e-8, 1e-7 and 1e-6 all fail at step 1 with DIVERGED_BREAKDOWN;
    //  1e-5 converges on 20/20 steps in at most 3 Newton iterations.
    //
    //  What that means: F is still not a smooth function of u at the 1e-6
    //  level. Two identified, UNFIXED causes -- localEnergyEnergyModel's
    //  updateDerived() writes a CLAMP back into nEps_, so the state F is
    //  evaluated at is not the state PETSc perturbed; and the species/nEps
    //  boundary conditions are corrected before the derived fields and
    //  mobilities they depend on. Removing those should let this drop toward
    //  PETSc's default; until then a case that diverges here should try
    //  RAISING mffdErr before anything else.
    //
    //  SUPERSEDED, recorded because both earlier versions of this comment
    //  were wrong and one was committed:
    //  (1) It first defaulted to 1e-4, attributing the breakdown to noise
    //      from the adaptive chemistry ODE. Wrong -- that ODE never fires on
    //      this case ("0 of 2000 cells integrated", 3 reactions).
    //  (2) It then claimed the meanE ordering fix made 1e-8/1e-6/1e-4 all
    //      converge identically. That measurement was INVALID: mffdErr_ still
    //      defaulted to 1e-4 at the time, so setPetscOptions() inserted
    //      -mat_mffd_err 1e-4 AFTER the swept environment variable and
    //      overrode it. All three runs were the same configuration, which is
    //      why they agreed to the last digit. The ordering fix is still
    //      correct and still worth having -- it moved the usable threshold
    //      from ~1e-4 to ~1e-5 -- but it did not make F smooth.
    scalar mffdErr_;

    //- Enforce the species density floor as a CONSTRAINT inside the Newton
    //  solve (PETSc SNESVINEWTONRSLS) instead of as a post-hoc clamp.
    //  DEFAULT true. The clamp is applied outside the equations, which makes
    //  F(u)=0 unreachable -- see doc/newton-outer-solver-design.md, defect A,
    //  and the cold-start guard in solveOuterStep.
    bool bounded_;

    //- Extra PETSc options, appended AFTER this class's own defaults so they
    //  WIN (PetscOptionsInsertString lets a later setting replace an earlier
    //  one). Empty by default, so a case that says nothing gets exactly the
    //  behaviour it got before this existed.
    //
    //  Exists because preconditioner tuning is EXPERIMENTAL WORK, and without
    //  this every trial -- a different sub-PC, a different Schur
    //  approximation, a raised Krylov cap -- needed a library rebuild. The
    //  binding constraint on this solver today is the PC's behaviour as dt
    //  grows (measured 2026-09-10: KSP iterations 14 -> 15 -> 17 -> 20 -> 32 ->
    //  100(cap) as dt rose toward 3.5e-10), so making that cheap to iterate on
    //  is the difference between tuning it and not.
    //
    //  Example, in the case's newtonSolver dict:
    //      petscOptions "-fieldsplit_0_pc_type hypre -fieldsplit_1_pc_type ilu";
    string petscOptions_;

    //- Normalise each field block's residual scale by the MEASURED initial
    //  residual, so every block starts the solve at |F_scaled| ~ 1.
    //
    //  `sF = sX/dt` is an a priori guess at the residual size and is a poor
    //  one. Measured 2026-09-10 on grubert2009_ballast400 the per-block scaled
    //  residuals spanned 0.0077 (n_Arp) to 111 (nEps_e) -- 14,000x across
    //  blocks that are all supposed to be O(1). A Krylov method minimises the
    //  norm of the WHOLE vector, so it fits the loud blocks and is nearly
    //  blind to the quiet ones however wrong they are.
    //
    //  Default TRUE. Set false to reproduce the pre-2026-09-10 behaviour, or
    //  to A/B the effect.
    bool rebalanceScales_;

    //- Whether the LAST SNES solve converged. Reported through
    //  lastSolveConverged() instead of being turned into a FatalError, so a
    //  step Newton cannot solve is DISCARDED AND RETRIED at a shorter dt via
    //  the same `onNonConvergence retryStep` path Picard uses.
    //
    //  Measured 2026-09-10 on grubert2009_ballast400, and this is why it
    //  matters: on handover Newton inherits whatever dt Picard's easy phase
    //  wound up to (1.29e-9 there) and dies on contact, while the SAME state
    //  solves cleanly at dt <= 1e-11 (154 and 256 steps, zero failures). The
    //  step was merely too big, and a solver with no way to report that has no
    //  way to find its own ceiling.
    bool lastConverged_ = true;

public:

    TypeName("SNES");

    snesNewtonSolver(const fvMesh& mesh, const dictionary& dict);

    virtual ~snesNewtonSolver() = default;

    virtual void solveOuterStep
    (
        electromagneticsModel& em,
        plasmaSpecies& species,
        plasmaTransport& transport,
        plasmaEnergy* energy
    ) override;

    //- Did the last solveOuterStep() converge? A false here makes soPlasmaFoam
    //  DISCARD the step and retry it at a shorter dt through plasmaTimeControl,
    //  rather than the run dying on a FatalError.
    virtual bool lastSolveConverged() const override
    {
        return lastConverged_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
