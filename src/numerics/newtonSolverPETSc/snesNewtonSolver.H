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
    SNES -- see docs/design/newton-outer-solver-design.md for the full design and
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
    //  F(u)=0 unreachable -- see docs/design/newton-outer-solver-design.md, defect A,
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

    //- Include the CHEMISTRY PRODUCTION derivative d(chemP)/d(n) in the
    //  preconditioner's species diagonal.
    //
    //  The species residual is F = ddt + div - lap - P + L*n, so
    //  dF/dn = L - dP/dn. The Pmat carried fvm::Sp(chemL) -- the L term -- but
    //  NOT -dP/dn. For electron-impact ionisation P_e = k_iz*n_e*n_Ar, so
    //  dP_e/dn_e = k_iz*n_Ar IS THE IONISATION FREQUENCY: the avalanche growth
    //  rate, and the defining physics of ignition. Omitting it leaves the
    //  preconditioner believing the species equation is far more diagonally
    //  dominant -- far more stable -- than it actually is, exactly when it is
    //  not.
    //
    //  dP_s/dn_s is approximated as P_s/n_s, which is EXACT whenever the
    //  production is first order in the species itself (electron-impact
    //  ionisation is), and is a reasonable approximation otherwise. A
    //  preconditioner need only approximate.
    //
    //  This term is NEGATIVE on the diagonal and deliberately so: it makes the
    //  block less diagonally dominant, because an avalanche genuinely is
    //  unstable. A preconditioner that hides that is describing a different
    //  problem.
    //
    //  Default TRUE; set false to A/B it.
    bool chemJacobian_;

    //- Include the CROSS-SPECIES chemistry coupling d(chemP_s)/d(n_e) in the
    //  preconditioner, as an off-diagonal (species, electron) block.
    //
    //  The same missing derivative as chemJacobian_, one step out: electron
    //  impact ionisation produces an ION for every electron it produces, so
    //  d(P_ion)/d(n_e) is as large as d(P_e)/d(n_e) and is what tells the
    //  preconditioner that the ion population is slaved to the electron
    //  avalanche. Without it the ion blocks look like independent
    //  advection-diffusion problems with a source that happens to be there.
    //
    //  APPROXIMATED as d(P_s)/d(n_e) = P_s/n_e, i.e. ALL of a species'
    //  production is treated as first order in the electron density. That is
    //  right for electron-impact channels and wrong for the rest (ion-neutral
    //  chemistry, three-body recombination), so it is an approximation and not
    //  a Jacobian -- which is exactly what a preconditioner is allowed to be.
    //  It is separately switchable from chemJacobian_ so the two can be told
    //  apart in a measurement.
    //
    //  Default FALSE until measured on this case.
    bool chemCrossJacobian_;

    //- Include d(Psrc)/d(phi) -- JOULE HEATING's response to the potential --
    //  in the energy row of the preconditioner.
    //
    //  The last missing block in the (energy, phi) coupling. That row already
    //  carries the energy DRIFT term's field derivative, but Psrc ~ J.E is a
    //  direct and typically stronger dependence on the potential. It is not
    //  laplacian-shaped -- it linearises to a convection operator,
    //  a.grad(dphi) with a = 2*mu_e*n_e*E/(1 V) -- which is why it was
    //  deferred rather than copied from the drift term.
    //
    //  ADAPTIVE to `energySource`, with no branch: the analytic Joule
    //  derivative is scaled by the energy model's own per-cell ratio
    //  Psrc/jouleHeating, which is identically 1 under `energySource model`
    //  (making the derivative EXACT) and carries the ODE's actual departure
    //  from the bare Joule term under `energySource chemistry`.
    //
    //  DEFAULT FALSE -- MEASURED, and it did NOT win. Equal start, equal
    //  endTime, both arms reaching it, on grubert (`energySource chemistry`,
    //  so the ratio-corrected path rather than the exact one):
    //
    //      mean dt          3.461e-11  vs  3.307e-11   (+4.7% for ON)
    //      SNES iters total       1710  vs        1517  (+12.7% WORSE for ON)
    //      failures/retries          8  vs           9
    //
    //  The timestep gain does not pay for the extra Newton iterations. Kept
    //  because the term is real and because this case exercises only the
    //  APPROXIMATE path: under `energySource model` the derivative is EXACT
    //  (the ratio is identically 1), and it may well pay there. That is
    //  untested -- do not assume it either way.
    //
    //  Same discipline as chemCrossJacobian_, and the same outcome: a
    //  plausible, correctly-derived Jacobian block that measurement declined.
    bool jouleJacobian_;

    //- Form the Schur complement on the PHI block instead of the transport
    //  block, and precondition it with the SEMI-IMPLICIT POISSON operator.
    //
    //  PETSc's fieldsplit eliminates split 0 and forms S on split 1, chosen
    //  purely by registration order. The historic order (phi first) gives
    //      S_t = A_tt - A_tf inv(A_ff) A_ft
    //  which is (nSpecies+1) x nCells, can only be APPLIED, and therefore runs
    //  on PETSc's `a11` default -- ILU(0) of A_tt. MEASURED on grubert_long:
    //  167 inner solves per outer step, median 230 Krylov iterations each,
    //  2.6% grinding to the iteration cap.
    //
    //  Reversing the order gives
    //      S_f = A_ff - A_ft inv(A_tt) A_tf
    //  a SCALAR ELLIPTIC operator on one field -- and
    //  div((eps + dt*sigma) grad .) IS that operator to leading order in dt,
    //  so it can be ASSEMBLED and handed to AMG. See
    //  docs/design/schur-semiimplicit-poisson-preconditioner.md for the derivation.
    //
    //  Default FALSE. Two reasons, not one:
    //   * it inverts a working, measured configuration, and Tier 0 (the inner
    //     KSP limits, now baked into the default options) already bought 9.3x
    //     on its own -- this has to beat THAT, not the old defaults;
    //   * inv(A_tt) becomes the inner solve, and an ITERATIVE inner solve
    //     taking a varying number of iterations is exactly what broke down and
    //     raised the SIGFPE of section 25c. Anything enabling this must also
    //     pin the transport solve to a FIXED iteration count, e.g.
    //       -fieldsplit_transport_ksp_type richardson
    //       -fieldsplit_transport_ksp_max_it 5
    //       -fieldsplit_transport_ksp_convergence_test skip
    //     so every application performs identical arithmetic.
    bool schurOnPhi_;

    //- PER-CELL state/residual scaling instead of one constant per FIELD.
    //
    //  sX/sF are one scalar per block, so a block is normalised by ONE number
    //  however wide its own spread is. A Krylov method minimises the norm of
    //  the WHOLE vector, so it fits the cells at the top of the range and is
    //  nearly blind to those at the bottom -- the same argument that justified
    //  rebalanceScales_ ACROSS blocks, applied one level down, WITHIN one.
    //
    //  DEFAULT FALSE, and the expected payoff is smaller than when this was
    //  first proposed. It was queued on a measured spread of 1e11-5.7e18
    //  (8 orders) inside the Arp block. Re-measured 2026-09-10 in the regime
    //  where Newton actually RUNS -- a developed discharge -- the spread is
    //  1.7e4 for Ar2p, 4.1e3 for Arp, 8.9e2 for electrons: FOUR orders, not
    //  eight. The 8-order figure belongs to the COLD state, where every
    //  species sits on its floor and the start-up guard keeps Picard engaged,
    //  so Newton never sees it.
    //
    //  Kept switchable rather than dropped because four orders is still a real
    //  spread; enabled by nobody until an A/B says it pays.
    //
    //  BIT-IDENTITY WHEN OFF is the regression gate, and the implementation is
    //  shaped around it: the residual and state sites use per-cell arrays
    //  UNCONDITIONALLY, filled with the uniform value when this is false, so
    //  the arithmetic is unchanged; the Pmat sites BRANCH, because folding sF
    //  into the row factor and sX into a column factor reorders the
    //  multiplication and floating-point multiplication is not associative.
    bool perCellScaling_;

    //- Floor for the per-cell state scale, as a FRACTION of the block's own
    //  maximum. Without it a cell at the density floor gets a scale of ~0 and
    //  its scaled residual explodes -- the mirror of the problem this exists
    //  to fix. 1e-6 keeps six orders of genuine per-cell range.
    scalar perCellScaleFloor_;

    //- Ceiling on the Psrc/jouleHeating ratio above, so a cell where
    //  jouleHeating is small but nonzero cannot produce a runaway coefficient.
    //  Bounds the term rather than dropping it, which keeps the block alive in
    //  a sheath -- exactly where it matters most.
    scalar jouleJacobianRatioMax_;

    //- INEXACT NEWTON with an Eisenstat-Walker adaptive forcing term
    //  (PETSc `-snes_ksp_ew`): choose each linear solve's tolerance from how
    //  well the previous Newton step actually reduced ||F||, instead of
    //  solving every one to a fixed rtol.
    //
    //  Early Newton iterates are far from the root, so solving their linear
    //  systems tightly buys nothing -- the step is going to be superseded. This
    //  is standard JFNK practice and is what Knoll & Keyes mean by an inexact
    //  Newton method.
    //
    //  MEASURED 2026-09-10 on the grubert2009_ballast400 pre-ignition state,
    //  first 90 steps, everything else equal:
    //
    //      without   mean dt 1.917e-10   KSP failure rate 2.6%
    //      with      mean dt 3.186e-10   KSP failure rate 0.1%
    //                        +80%                         13x fewer
    //
    //  IT FAILED WHEN FIRST TRIED, and why it failed matters: that test ran
    //  against a Pmat that was missing both the Schur coupling block and the
    //  ionisation derivative. An adaptive forcing term reads the quality of the
    //  Newton steps it is given; with a Jacobian missing the dominant physics
    //  those steps carry no signal, so adapting to them is worthless. Fixing
    //  the matrix is what made this option pay.
    //
    //  Default TRUE. It is a real behaviour change for every existing case, so
    //  it is switchable and the measurement above is the justification.
    bool adaptiveForcing_;

    //- Start the Newton solve from a LINEAR EXTRAPOLATION of the history
    //  instead of from the previous step's solution.
    //
    //      u0 = u^n + (dt/dt_old)*(u^n - u^{n-1})
    //
    //  The initial guess was u^n -- a ZEROTH-order guess. Newton is only
    //  LOCALLY convergent, and as dt grows the root moves further from u^n
    //  until u^n leaves its basin of attraction; Newton then converges to a
    //  local minimum of ||F|| rather than to the root.
    //
    //  That is exactly what the failing steps look like (measured 2026-09-11 on
    //  grubert2009_ballast400): ||F|| falls to 0.29 of its initial value and
    //  FLATLINES, while successful steps reach 1e-8 in 4-10 iterations. And
    //  everything that would explain it differently has been eliminated by
    //  measurement -- F is provably exact (purity and hysteresis checks both
    //  return ||dF|| = 0), the differencing step is irrelevant (mffdErr 1e-5 /
    //  1e-7 / PETSc default all give the same dt to within 1%), the KSP
    //  converges, and relaxing errRtol 100x does not move dt.
    //
    //  Densities are floored at their own minNumberDensity: an extrapolated
    //  guess may undershoot to a negative value, which is not merely infeasible
    //  but breaks the rate-table lookups the residual performs.
    //
    //  Default FALSE until measured.
    bool extrapolateGuess_;

    //- Extrapolate the POTENTIAL too, or only the transported fields.
    //
    //  These are NOT the same case and the distinction is physical. Densities
    //  and nEps carry a ddt term: they genuinely evolve, so
    //  u^n + r*(u^n - u^{n-1}) is a proper second-order predictor and is
    //  standard practice. ePotential has NO time derivative -- it is set
    //  instantaneously by Poisson from the charge density -- so extrapolating
    //  it in time is not motivated the same way, and the extrapolated phi may
    //  be INCONSISTENT with the extrapolated rho. The consistent guess for phi
    //  would be the solution of Poisson AT the extrapolated densities.
    //
    //  Kept separately switchable so the question is answered by measurement
    //  rather than by argument. Default TRUE (i.e. extrapolate it) only
    //  because phi(t) is smooth in the cases at hand; set false to extrapolate
    //  the transported fields alone.
    bool extrapolatePotential_;

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
