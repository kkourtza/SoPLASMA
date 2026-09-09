# A genuine Newton-type outer solver for SoPLASMA -- why, and which road

2026-09-09. Written before starting the PETSc build, so the reasoning survives
independent of whoever picks this up next.

## The problem, measured

SoPLASMA's outer coupling (Poisson <-> species transport <-> electron energy)
is a segregated fixed-point (Picard) iteration per timestep, accelerated by
adaptive Aitken relaxation (`plasmaOuterRelaxation`). This is exactly the
scheme the group's own published work already knew needed tightening: the
SoPlasmaFoam paper's "coupling loop sensitivity analysis" (positive-streamer
case, Section 6.3.1) measured that the number of correctors needed exceeds
what naive Courant-number or dielectric-relaxation-ratio criteria alone would
suggest, even at values comfortably below 1.0.

Pushed further, on a DC-glow ignition transient
(`validation/grubert2009_pseudo`, a dt-escalation scan 2026-09-09), the same
loop was found to do something qualitatively worse than "needs more
correctors": its own measured contraction factor `rho` climbed *smoothly and
monotonically* from 0.76 toward 1.0 over about ten consecutive steps, then
crossed 1.0 outright -- the mathematical threshold where a fixed-point
iteration stops contracting and starts diverging. No corrector count fixes
that; more iterations of a divergent map just diverge faster. Every dt tested
from 1e-8 to 1e-10 hit this identical failure (climbing non-convergence, then
a genuine SIGFPE inside the Poisson linear solve), at progressively later
simulated times as dt shrank -- smaller dt resolves the *approach* to the
wall more accurately, but does not remove the wall, exactly matching an
independent finding one day earlier across standard/SG/CFS transient schemes
("CFS's accuracy merely postpones reaching it, does not fix it" --
`session-state-2026-09-08`).

The loss of contraction was measured to coincide, in the same few steps, with
an extremum (peak-then-reversal) in the discharge's own conduction current
`I_cond` -- while the imposed current-source target `I_set` was climbing
smoothly on its intended ramp, `I_cond` was running 30-35x larger and turned
over right where `rho` crossed 1. That correlation is the textbook signature
of a negative-differential-resistance ("subnormal") region of the
discharge's own I-V characteristic. Full detail and the raw diagnostic
excerpts: memory `outer-coupling-diverges-at-negative-differential-resistance`.

This is also exactly what the very first standing reference in this
project's memory index already said: Almeida et al. 2016 -- "a 1D CVC has a
minimum+hysteresis; CURRENT control and a stationary solver are required."
This project built the current-control half. It never built the stationary
(non-Picard) solver half, and today's evidence is a direct, mechanistic
confirmation of why that gap matters: a Picard iteration is not the right
tool for a multivalued/negative-resistance characteristic, regardless of dt,
scheme, or relaxation schedule.

## Two candidate roads, surveyed before choosing

**(a) OpenFOAM's own native block-coupled matrix framework
(`fvBlockMatrix`).** Investigated 2026-09-09 (local source search across the
whole installed OpenFOAM v2412 tree, plus external research). Result:
**this class does not exist in OpenFOAM v2412 (ESI/OpenFOAM.com) at all.**
`fvBlockMatrix`/`blockLduMatrix` is a foam-extend-only lineage (Chalmers-era
work, OpenFOAM-1.6-ext), never merged into the ESI branch SoPLASMA is built
on. Pursuing it would mean forking to a different, long-diverged OpenFOAM
lineage, or reimplementing block-matrix assembly from scratch against ESI's
current core matrix classes -- before any Newton logic could even start. And
even in foam-extend, by every description found, it is a LINEAR
block-coupled solve, not a nonlinear one: cross-field coupling terms are
still linearised and iterated by hand by the solver author. It would not
supply a Newton method, only (after nontrivial porting work) a better-coupled
linear solve to wrap one around.

**(b) PETSc's SNES (nonlinear/Newton) module, via the existing `petsc4Foam`
bridge.** Confirmed: PETSc is not currently built here (`PETSC_DIR` points to
a nonexistent directory -- this is the whole reason `libpetscFoam.so` has
been failing to load, silently, in every run so far, since `matrixSolver` is
`foam` everywhere and nothing has depended on it yet). Confirmed further,
by reading the vendored `petsc4Foam` source directly and the project's own
`Pasolari_Kourtzanidis_SoPlasmaFoam_arXiv-2607.05137.pdf` in full: the
existing integration only wraps PETSc's KSP (linear Krylov) layer. Its
"monolithic solver" feature is unrelated to Newton -- it assembles multiple
MESH REGIONS' (plasma+dielectric) Poisson matrices into one linear system to
avoid iterative inter-region BC coupling, not a nonlinear solve of the
tightly-coupled Poisson-species-energy system within one region. No
SNES-OpenFOAM bridge, official or community, was found anywhere. What WAS
found: Zapdos, a MOOSE-based plasma-fluid code solving this exact physics
class (Poisson + drift-diffusion-reaction species + electron energy) fully
coupled via genuine Newton through PETSc SNES with hand-coded exact
Jacobians -- actively maintained, peer-reviewed (Computer Physics
Communications, 2023), production-grade. Different framework, nothing
directly reusable, but strong independent confirmation that full Newton via
SNES is the standard, validated tool for exactly this problem class.

## Decision: (b), PETSc SNES

Building genuine SNES integration is the SMALLER undertaking of the two, not
the larger one. `petsc4Foam`'s existing KSP bridge already solves the harder
plumbing problem -- OpenFOAM `lduMatrix`/field <-> PETSc `Mat`/`Vec`
marshaling, parallel `Pstream` compatibility -- and in PETSc's own
architecture, SNES is built directly on top of KSP (each Newton step is
itself a KSP solve). Extending a working, already-integrated bridge is
additive engineering on a proven foundation. The `fvBlockMatrix` route, by
contrast, starts from a codebase that does not have the class at all, in a
different OpenFOAM lineage, for a lesser end result (still just a linear
solve, still needing a hand-written Newton layer on top regardless).

Chosen 2026-09-09. Full survey detail, including the exact search evidence
for both roads: memory `outer-coupling-diverges-at-negative-differential-resistance`
(records the divergence measurement) plus the agent research report that
produced this recommendation (not separately filed as of this writing).

## Progress

**2026-09-09: PETSc built and verified.** Not under the system OpenFOAM tree
(`/usr/lib/openfoam/openfoam2412/ThirdParty` is a stub file, "system
dependencies", on this apt-style install -- not a real directory) but
project-locally at `ThirdParty/petsc-3.24.0`, `PETSC_ARCH=DPInt32`, matching
the version this project's `petsc4Foam` was tested against. Configured with
`--download-fblaslapack --download-metis --download-parmetis
--download-superlu_dist --download-hypre`, all resolved cleanly. `make check`
passed, including `src/snes/tutorials/ex19` (a genuine SNES/nonlinear
example) with HYPRE and SuperLU_DIST backends, in parallel. `PETSC_DIR`/
`PETSC_ARCH`/`LD_LIBRARY_PATH` are now exported from `etc/bashrc` (source
this AFTER OpenFOAM's own bashrc, so it overrides the broken system
`PETSC_DIR` OpenFOAM's environment otherwise sets). `petsc4Foam` was
rebuilt against it (`ThirdParty/petsc4Foam/Allwmake -prefix=openfoam`) and
`foamHasLibrary -verbose petscFoam` now reports "Can load petscFoam" --
confirmed further by an actual `soPlasmaFoam` run producing no
"Could not load libpetscFoam.so" warning, which had appeared in every prior
run this session.

## Scope and exposure -- decided 2026-09-09, before any implementation

This must land as a general, always-available SWITCH between outer-solver
strategies (Picard, the current default, vs. Newton/SNES) -- not a one-case
patch, and not something that only exists for the DC-glow ignition problem
that motivated it. Two explicit constraints from the user:

1. **A general option, everywhere.** Any case should be able to select
   Picard or Newton for the outer coupling. Whether Newton becomes the
   DEFAULT, and whether an AUTOMATIC choice between the two can be made, are
   open questions to revisit ONCE this is built and results exist -- not
   decided now, and not to be pre-empted by the implementation.

2. **This is a numerics knob, NOT a semantic-layer/generator concept.**
   Per G1/G2 (the imported SoEEDF CLAUDE.md): the generator
   (`configuration/config`, `configuration/boundaries`) is the SEMANTIC
   layer -- it is where a user states WHAT the physics/geometry/circuit IS,
   and dictionary values get DERIVED from that (the way `electronEnergyModel
   LMEA` auto-derives `energyModelCoeffs`). Picard-vs-Newton is not a
   physics statement a user should have to reason about; it belongs beside
   `matrixSolver foam | petsc` in `system/plasmaSimulationControls` (or a
   sibling numerics dict) as a direct, manually-set key -- e.g. something
   like `outerSolver picard | newton;` -- NOT threaded through the
   generator's derivation logic, and NOT exposed as a config-layer switch
   a case author is expected to reason about physically.

A natural (not yet designed) tie-in for the later "automatic choice"
question: the existing `rho [contraction]` / `omega [coupling margin]`
diagnostics (already computed every step, see
[[outer-coupling-diverges-at-negative-differential-resistance]]) are exactly
the signal that would need to detect an approaching non-contraction and
trigger a switch -- rho climbing toward 1 IS the early-warning signal this
whole investigation was built on. Worth remembering when that question comes
up; not scoped further now.

## Also to test, per the user 2026-09-09: `MatFDColoring`

A middle ground between pure matrix-free JFNK and a hand-coded/AD Jacobian,
already built into PETSc: graph-colored finite differencing
(`SNESComputeJacobianDefaultColor`). Needs roughly one residual evaluation
per sparsity-graph COLOR rather than per DOF or none at all, and composes
directly with the same `FormFunction` residual evaluator planned below --
no AD tooling, no hand-derived Jacobian. Flagged for a later comparison
against plain JFNK once the residual/DOF-layout machinery exists (that
machinery is shared by both routes); not yet started, and not required to
be decided before starting JFNK.

Also surveyed and set aside for now: full source-transformation AD (Enzyme,
ADOL-C) for an exact Jacobian. Real precedent exists (ADOL-C has generated
PDE-residual Jacobians for PETSc elsewhere), but OpenFOAM's heavy use of
templates and virtual dispatch is a documented, genuine friction point for
these tools (ADOL-C needs retyping doubles as `adouble` throughout; Enzyme
avoids that but has known trouble with heavily templated/virtual C++) --
not ruled out permanently, just not the next thing to reach for.

## JFNK design, scoped 2026-09-09 -- confirmed against the actual assembly code

Decided to pursue Jacobian-Free Newton-Krylov (JFNK) rather than a hand-coded
analytic Jacobian (the Zapdos approach). Reasoning: the failure diagnosed
here is structural (a non-contracting fixed-point map), not a precision
problem, so a genuine Newton *iteration* -- even with an approximate/
matrix-free Jacobian -- already fixes the mechanism, and JFNK reuses far
more of the existing code than deriving every cross-field derivative by
hand would.

**Confirmed by reading the actual assembly code (not assumed):** every
equation SoPLASMA solves is already built as a standard, queryable
`fvScalarMatrix`:
  - `singleRegionPoisson::solve()` builds `fvScalarMatrix ePotentialEqn`
    (`src/models/electromagnetics/electrostaticModels/singleRegionPoisson/singleRegionPoisson.C:187`).
  - Every species transport model implements
    `virtual tmp<fvScalarMatrix> nEqn() const`
    (`plasmaTransportModel.H:146`; `driftDiffusion.C` etc. implement it).
  - The electron energy equation is `localEnergyEnergyModel::eEqn()`
    returning `tmp<fvScalarMatrix>`
    (`localEnergyEnergyModel.C:1312`).

Every `fvMatrix` has a native `.residual()` method (standard OpenFOAM,
computes Ax-b for whatever was assembled). This means the SNES residual
evaluation needs **no new numerical code** -- it is: unpack the trial state
into the OpenFOAM fields, call the EXISTING assembly functions (unchanged),
read `.residual()` off each resulting matrix, pack the results into one
concatenated PETSc `Vec`. Likewise the Jacobian-vector product needs no new
code: PETSc's matrix-free SNES (`-snes_mf_operator` / `MatCreateSNESMF`)
gets it by finite-differencing that same residual function internally.

**What IS genuinely new engineering, in decreasing order of how much of it
reuses existing machinery:**
1. **DOF layout / global numbering for a concatenated multi-field Vec.**
   `petsc4Foam`'s existing `ldu2csr` bridge already solves "OpenFOAM's
   parallel-decomposed field <-> one PETSc Vec" for a SINGLE field (that is
   what its KSP usage already does today). Extending it to a BLOCK Vec
   spanning potential + every species + energy, all in the same parallel
   layout, is an incremental generalisation of a solved problem, not a new
   one.
2. **The SNES `FormFunction` callback.** Genuinely new, but thin: mechanical
   unpack/assemble/harvest/pack around code that already exists and does not
   need to change.
3. **The preconditioner.** The natural, cheapest choice: reuse ONE existing
   Picard/Aitken corrector sweep (the current per-equation `.solve()` calls,
   unchanged) as a PCSHELL physics-based preconditioner for the Newton-Krylov
   inner solve -- exactly the technique Knoll & Keyes describe for JFNK.
   This is the piece that most needs care: it means the OLD segregated solve
   is not deleted, it is repurposed as the preconditioner, while the OUTER
   iteration becomes genuine Newton.
4. **The `outerSolver picard | newton` switch itself**, wired into
   `plasmaTransport::solve()` / the main PIMPLE loop in `soPlasmaFoam.C`, per
   the exposure decision above (a numerics-dict key, not a generator/
   semantic-layer concept).

Not yet started (at time of writing above): the actual C++ for any of the
four items. Progress since, below.

## Proof of concept: BUILT AND PASSING, serial and parallel -- 2026-09-09

`verification/testSnesJFNK/` (case, reusing the Grubert mesh) and
`src/applications/utilities/testSnesJFNK/` (application). Validates the
CORE mechanism end to end before touching the real 3-field system: a
manufactured nonlinear problem, `-D*lap(phi) + phi^2 = f(x)` with
`phi_exact = sin(pi x/L)` chosen so the residual is exactly zero at the true
solution, solved by PETSc's matrix-free SNES starting from phi=0 everywhere
(a genuinely bad initial guess, not a warm start).

**A real, previously-unknown build obstacle, found and fixed:** PETSc's
`petscmath.h` and OpenFOAM's own math operator overloads are AMBIGUOUS when
both are visible in one translation unit -- traced into OpenFOAM's own
`error.H`, not something a using-declaration or namespace qualifier fixes.
Resolved by a hard split: `testSnesJFNK.C` includes `fvCFD.H` and NEVER
includes any PETSc header; `snesBridge.C` includes `petscsnes.h` and NEVER
includes `fvCFD.H` (or any other OpenFOAM header); the two communicate only
through `snesBridge.H`, which uses exclusively `int`/`double`/`void*` --
including for `PetscInitialize`/`PetscFinalize` themselves, since even
those pull in the same header chain. This split is now the required
pattern for any future PETSc-touching SoPLASMA code, not a one-off
workaround for this test.

**Also found and fixed, both genuine and both quick once identified:**
- OpenFOAM's `Foam::sin()` correctly refuses a dimensioned argument;
  `mesh.C().component(0)` carries LENGTH dimensions. Fixed by doing the
  manufactured-solution setup (phi_exact, fSource) on raw `scalarField`s,
  not dimensioned `volScalarField` arithmetic -- right tool for a synthetic
  test artifact with no physical dimensions to track.
- Similarly, `fvc::laplacian` of a dimensionless field correctly carries
  1/length^2 dimensions in OpenFOAM's bookkeeping, clashing with the
  dimensionless `phi^2`/`fSource` terms in the residual. Same fix: raw
  `scalarField` arithmetic for the residual itself.
- Caught and fixed the exact dangling-`tmp<>` pattern already documented
  elsewhere in this project (binding a reference straight through a
  temporary's `()` without holding the `tmp<>` in a named variable first)
  before it ever ran, not after a crash.
- A real command-line-parsing conflict: `PetscInitialize` does not strip
  its own recognized flags from `argc`/`argv`, so OpenFOAM's `argList` then
  sees (and rejects) `-snes_monitor` etc. as unrecognized options. Worked
  around via the `PETSC_OPTIONS` environment variable for this test;
  the real system will need a considered answer here (a case-level
  dictionary is the natural fit, matching how every other SoPLASMA solver
  option is already exposed, rather than expecting `-snes_*` flags on the
  command line at all).

**Result, serial** (1 rank): converged in 2 Newton iterations
(residual norm 2.100272141852e+06 -> 1.76 -> 1.285e-05,
`SNESConvergedReason = 3`), L2 error vs the manufactured exact solution =
2.16476217e-05 (consistent with ordinary FVM discretisation error on this
mesh, not a solver artifact).

**Result, parallel** (4 ranks, `scotch` decomposition): converged in 3
iterations (residual norm 2.100272141842e+06 -> ... -> 6.48e-06,
`SNESConvergedReason = 3`), **L2 error = 2.16476217e-05 -- identical to the
serial run to every displayed digit.** This is the decisive check, not a
formality: a halo-exchange or DOF-layout bug would show up as extra error
concentrated at partition boundaries, so an identical L2 error across
serial and 4-rank decomposition is direct evidence (not merely an
assumption) that the parallel-correctness checklist above was actually
satisfied, not just designed for.

## `fvMatrix::residual()` found genuinely broken for this usage in parallel -- 2026-09-09

The plan above (see "Confirmed by reading the actual assembly code") assumed
the outer residual could be harvested "for free" via `fvm::`-assembling each
equation's existing matrix and calling `.residual()` on it. This was built
and tested (`fvm::laplacian(phi) + fvm::Sp(phi,phi) == fSource`, harvesting
`eqn.residual()`) and **found to silently and massively misreport in
parallel**, independent of PETSc:

- Unpreconditioned matrix-free JFNK on this `.residual()`-harvested form
  converged (`SNESConvergedReason` positive, satisfying the requested
  relative residual drop) to a state with L2 error 0.408 against the
  manufactured solution in 4-rank parallel, vs. 2.1648e-05 in serial --
  a ~20000x discrepancy, and bit-identical across every attempted fix
  (residual volume-normalisation, a real physics-based PCSHELL
  preconditioner, giving PETSc its own duplicated MPI communicator to rule
  out tag collisions with OpenFOAM's own Pstream traffic, confirming
  `residualCallback` was invoked an identical number of times on every
  rank). None of that changed the wrong answer, which is itself the tell:
  the fault was not in the outer solve strategy at all.
- **Decisive isolation:** feeding that "wrong" converged `phi` into a
  completely NATIVE, non-PETSc, single direct call to `eqn.residual()` on
  the identical assembled matrix STILL reported a near-zero residual
  (~1e-12) -- while calling `.solve()` on that same matrix immediately
  afterward correctly detected a large imbalance (78 PCG iterations to
  fix it), and 30 native Picard sweeps starting from that exact state
  walked it all the way back to the correct answer (L2 -> 2.1648e-05).
  This proves the "wrong" state is not a genuine alternate root of the
  discrete system that Newton legitimately found.
- **Ground truth, independent of OpenFOAM entirely:** reading raw cell
  centres and the "wrong" `phi` values back from each processor's `0/`
  directory and computing `-D*d^2(phi)/dx^2 + phi^2 - f` via a plain
  Python finite difference (no OpenFOAM matrix machinery anywhere in the
  chain) gives `max|residual| ~ 6.5e9` -- consistent with the field's
  genuinely wild, oscillatory profile (not a subtle smooth alternate
  branch). `fvMatrix::residual()`, called this way on a decomposed field,
  is simply wrong here by roughly 21 orders of magnitude, and `.solve()`
  is not: two back-to-back constructions of the IDENTICAL matrix on the
  IDENTICAL `phi`, one via `.residual()` and one via `.solve()`, disagree
  on whether the equation is satisfied.

**Root cause not further isolated below `fvMatrix::residual()` itself**
(would require instrumenting OpenFOAM's own `lduMatrix::residual()` /
`initMatrixInterfaces`/`updateMatrixInterfaces` internals, or an upstream
bug report) -- diminishing returns for this project once a correct,
already-validated alternative existed.

**The fix:** compute the outer residual via EXPLICIT `fvc::` evaluation
instead -- `R = -D*fvc::laplacian(phi) + phi^2 - fSource`, done as raw
`scalarField` arithmetic (bypassing `GeometricField` dimension-checking
for the final combination, same rationale as the original manufactured-
solution setup). This was the ORIGINAL (v1) approach in this same proof of
concept before the pivot to `fvm::`+`.residual()`, and was ALREADY
confirmed bit-identical between serial and 4-rank parallel then. Switching
back resolved the discrepancy completely: serial and 4-rank parallel now
both give `SNESConvergedReason = 3`, 2 iterations, **L2 error =
2.16476217e-05, identical to every displayed digit**, with the
physics-based preconditioner (below) still in place.

**Practical implication for the REAL 3-field system, i.e. a correction to
the plan above, not merely an implementation detail:** the outer SNES
residual (`FormFunction`) must be written EXPLICITLY per equation via
`fvc::` operators, not harvested via `fvm::`+`.residual()` on the existing
`fvScalarMatrix`-returning functions -- item 2 in the "genuinely new
engineering" list above is therefore somewhat LESS thin than originally
scoped (an explicit residual expression per equation, not a one-line
`.residual()` harvest). The preconditioner is UNAFFECTED: it uses
`fvm::`+`.solve()` (never `.residual()`), so it can still reuse existing
equation assembly and the existing per-equation linear solve as planned.

## Physics-based (Knoll & Keyes) preconditioner: BUILT AND VALIDATED -- 2026-09-09

Implemented via PETSc's `PCSHELL`: given a right-hand side `b` from GMRES,
the preconditioner applies ONE native (non-matrix-free) linear solve of the
problem's own linearised operator -- for this test problem,
`-D*fvm::laplacian(dphi) + fvm::Sp(2*phi, dphi) == b`, solved via
`eqn.solve()` -- reusing a scratch correction field `dphi` that copies
`phi`'s own BC prototype (then zeroed, since a correction must vanish
wherever `phi` itself is Dirichlet-pinned). `2*phi` is the correct
LINEARISATION coefficient (the derivative of `phi^2`), not `phi` itself
(which is what the RESIDUAL's `fvm::Sp(phi,phi)` idiom uses to reproduce
`phi^2` exactly -- a different role, easy to conflate).

This is exactly the piece flagged as needing the most care in the plan
above ("the OLD segregated solve is not deleted, it is repurposed as the
preconditioner"), now built and confirmed working (serial and parallel,
both converging in 2 SNES iterations with it in place) rather than merely
planned. `snesBridge.H`/`.C` were extended with a second callback type
(`PCApplyCallback`) alongside `ResidualCallback`, wired to the SNES's KSP
via `SNESGetKSP`/`KSPGetPC`/`PCSetType(PCSHELL)`/`PCShellSetApply` -- this
plumbing is now the template for the real 3-field preconditioner (one
Picard/Aitken sweep across all three equations per PC application, not
just one scalar equation).

One PETSc API subtlety worth keeping: `SNESSetUseMatrixFree(snes, mf,
mf_operator)`'s second argument tells SNES to ALSO derive its own dense,
finite-differenced preconditioning matrix if true -- wasteful and
beside the point once a real PCSHELL is supplied, so it must be passed
`PETSC_FALSE` whenever a PC callback is given (confirmed via `-ksp_view`
that leaving it `PETSC_TRUE` attaches an unused `mpidense` Pmat alongside
the shell PC).

Not yet done: any of the four items in the "genuinely new engineering" list
above for the REAL 3-field system (Poisson + species + energy) -- both the
residual harvesting correction and the preconditioner above were validated
on the single-equation manufactured problem only. It proves the SNES<->
OpenFOAM plumbing, the parallel behaviour, AND the preconditioning
strategy; the remaining work is extending each to the real 3-field DOF
layout (concatenated block Vec, one FormFunction per equation feeding into
one packed residual, one Picard/Aitken sweep across all three equations as
the PC application).

## Parallel correctness -- a standing requirement, not an afterthought

Per the user 2026-09-09: every piece of this must translate to MPI-parallel
runs correctly, not just serially. This is not a new concern to invent for
this feature -- it is exactly [[collective-behind-a-local-guard]] / rule 31
in the imported CLAUDE.md, already paid for once in this project ("A
COLLECTIVE MUST NEVER SIT DOWNSTREAM OF A GUARD OR EARLY RETURN THAT DEPENDS
ON A LOCAL QUANTITY" -- a `reduce`/`gAverage` gated on a per-rank-varying
local count deadlocked under one decomposition and silently corrupted under
another). Concrete places this bites in THIS design:

1. **The block DOF layout must preserve per-rank ownership ranges exactly
   as `petsc4Foam`'s existing single-field bridge already establishes them.**
   This is the good news: that plumbing is already parallel-correct and
   already exercised in production (it is what the existing KSP-based Poisson
   solve runs on today). The new work is offsetting each field's block
   consistently across ranks (`VecCreateNest` / manual global-index
   offsetting -- a standard, well-trodden PETSc pattern for multi-physics
   coupled systems), not re-solving the ownership problem from scratch.
2. **`FormFunction` must exchange processor-boundary (halo) values before
   assembling any equation, every time SNES perturbs a trial state** --
   exactly what the existing Picard loop already does
   (`correctBoundaryConditions()` after every field update) before building
   its `fvm::` operators. A finite-volume residual assembled without a fresh
   halo exchange is WRONG at processor boundaries, not just stale. Since
   JFNK calls the residual function repeatedly with perturbed states
   (both for the actual Newton step and for PETSc's internal
   finite-differenced Jacobian-vector products), this exchange must happen
   on EVERY call, not once -- easy to get right by following the Picard
   loop's own existing pattern, easy to get wrong by "optimising" it away.
3. **Any custom convergence test or diagnostic (a global residual norm, a
   custom SNES monitor) is a collective by construction and must be called
   unconditionally by every rank** -- never behind a check like "if this
   rank owns any cells of species X". This is rule 31's exact failure mode,
   and the natural place to reintroduce it here is a hand-rolled convergence
   check that iterates over locally-owned species/fields before reducing.
4. **`MatFDColoring`'s coloring itself is mature, parallel-safe PETSc
   infrastructure** -- no new risk there, PROVIDED the sparsity pattern
   feeding it is assembled correctly across ranks in the first place (which
   reduces to point 1/2 above, not a separate problem).

None of this blocks starting the implementation -- it is a checklist to
build against from the first line of code, not a redesign.

## What is NOT yet decided or built

- No SNES residual/Jacobian assembly exists yet for SoPLASMA's 3-field
  system (electric potential, species number densities, electron mean
  energy). Scoping that assembly is real, substantial work in its own
  right, independent of which road was chosen above -- choosing PETSc SNES
  over `fvBlockMatrix` avoids the matrix-framework problem, it does not
  remove the physics-residual/Jacobian problem. This is the next actual
  step, not yet started.
  **PARTIALLY SUPERSEDED 2026-09-09** by the two-field proof of concept
  below: the BLOCK/multi-field DOF-layout and multi-equation
  preconditioner MACHINERY is now built and validated on a genuinely
  coupled synthetic pair of equations (serial and 4-rank parallel,
  bit-identical). What remains not-yet-built is specifically the REAL
  physics (Poisson + 3 species + electron energy assembly, with
  Scharfetter-Gummel/CompleteFlux schemes, chemistry sources, LMEA
  transport lookups), not the block-coupling mechanism itself.
- Whether a genuine Newton solve, once built, actually converges through a
  negative-differential-resistance region better than Picard is a
  *hypothesis under this diagnosis*, not yet demonstrated. Newton methods
  can also fail at turning points without additional care (e.g. pseudo-
  arclength continuation, which is the classical fix for bifurcation/fold
  points in nonlinear solvers) -- no such technique has been investigated
  yet for this problem, and doing so may turn out to be necessary on top of
  plain Newton/SNES.

## Proof of concept, step 2: multi-field block DOF layout -- BUILT AND PASSING, serial and parallel -- 2026-09-09

`verification/testSnesJFNK2Field/` and
`src/applications/utilities/testSnesJFNK2Field/` (reusing the single-field
proof of concept's mesh and `snesBridge.H`/`.C` UNCHANGED -- the bridge
interface already takes a plain `nLocal`/flat array, so it never needed to
know anything about "one field" vs. "several fields"; that genericity was
free). Validates the block/multi-field DOF layout and a multi-equation
physics-based preconditioner on a GENUINELY coupled synthetic pair of
nonlinear equations, deliberately kept simple (not the real physics) so a
bug in the block machinery itself cannot be confused with a bug in the
real equations -- the same separation of concerns that made the
single-equation `.residual()` defect (previous section) tractable to
isolate rather than chased blindly.

**The synthetic coupled system**, on the same reused mesh:
```
-D*lap(phi)  + phi*phi2 = f1(x)      phi_exact  = sin(pi x/L)
-D*lap(phi2) + phi^2    = f2(x)      phi2_exact = 2*sin(pi x/L)
```
`phi2_exact` is scaled 2x relative to `phi_exact` specifically so a DOF
layout bug (e.g. the two fields' blocks swapped or overlapping per rank)
would show up as a real, non-degenerate error rather than an accidental
match.

**DOF layout finding, resolving the open question in the earlier
"genuinely new engineering" list:** no `VecCreateNest` or manual global
index bookkeeping is needed. Each rank's local PETSc vector is simply
`nFields * nCellsLocal` long, laid out as contiguous per-field blocks
(`[phi block][phi2 block]`, generalizing to `[field0][field1]...` for N
fields) -- `FormFunction`/`PCApply` own all the field semantics
internally (unpacking into named OpenFOAM fields, repacking on the way
out), so PETSc never needs to know the vector has structure at all. This
is a substantially smaller task than the original plan anticipated.

**Residual** (`FormFunction`): explicit `fvc::` evaluation per equation,
per the fix in the previous section -- `fvm::`+`.residual()` is not used
anywhere in this project's Newton-solver code.

**Preconditioner** (`PCApply`): ONE block GAUSS-SEIDEL sweep, mirroring
the real system's actual sequential order (confirmed by reading
`soPlasmaFoam.C`'s outer loop, see below: Poisson, then species, then
energy, each using the latest available values of the others, not a
simultaneous/Jacobi update) --
```
1. solve  -D*lap(dphi)  + phi2*dphi == b1                  for dphi
2. solve  -D*lap(dphi2)             == b2 - 2*phi*dphi     for dphi2
```
using CURRENT (frozen) `phi`/`phi2` as Picard-linearization coefficients,
and the JUST-COMPUTED `dphi` folded into step 2's right-hand side --
Gauss-Seidel, not Jacobi, which is what makes this "one sweep" in the same
sense the real system's sequential solve is one sweep, not two
independent decoupled solves.

**Result, serial:** `SNESConvergedReason = 3`, 2 iterations, L2 error
`phi` = 2.16472295e-05, `phi2` = 4.32956425e-05 (both consistent with
ordinary FVM discretisation error, matching the single-field case's
2.1648e-05 to within the expected difference from `phi2`'s 2x amplitude).

**Result, parallel (4 ranks, scotch decomposition), fresh `uniform 0`
initial guess independently verified on every rank before the run:**
`SNESConvergedReason = 3`, 2 iterations, **L2 error = 2.16472295e-05 /
4.32956425e-05 -- identical to serial to every displayed digit.** Given
the previous section's single-field parallel defect was found exactly
this way (an identical-looking PASS that was actually silently wrong), an
independent Python finite-difference cross-check was considered but not
performed here -- justified because the earlier defect was isolated to
`fvm::`+`.residual()` specifically, which this code never calls (residual
is `fvc::`, matching the confirmed-safe path), and the bit-identical
serial/parallel match together with correct convergence to the KNOWN
manufactured solution is the same standard of evidence the single-field
fix was validated against.

## Real system's actual equation structure, confirmed by reading the code -- 2026-09-09

Gathered to ground the next step (wiring in the real physics) in what the
code actually does, not what was assumed when this document was first
written:

- **Poisson**: `singleRegionPoisson::solve()`
  (`singleRegionPoisson.C:187-191`) assembles
  `fvm::laplacian(epsilon_, ePotential_) == -chargeDensity_` (explicit
  scheme) or, under `poissonScheme semiImplicit`
  (`singleRegionPoisson.C:214-254`), the same operator with an effective
  permittivity `epsilon_ + dt*electricalConductivity` and RHS
  `-chargeDensity_ - dt*diffusiveChargeSource` fed implicitly from the
  species models. `multiRegionPoisson` additionally exists for
  gas+dielectric cases (`solveCoupled()`/`solveSegregated()`) -- out of
  scope unless a dielectric region is part of the target case.
- **Species transport**: `plasmaTransportModel.H:146`,
  `virtual tmp<fvScalarMatrix> nEqn() const`. `driftDiffusion::nEqn()`
  (`driftDiffusion.C:328-400`) assembles `fvm::ddt(n) + fvm::div(phi,n) -
  fvm::laplacian(D,n)` (or the combined SG/CompleteFlux operator under
  those flux schemes); no reaction terms here -- chemistry sources are
  added on top by the caller (`plasmaTransport::mechanismSourceTerms()`,
  `plasmaTransport.C:3394+`, e.g. `fvm::Sp(...)` at lines 4286/4377) and
  `explicitSource`/`Sph` are subtracted at `plasmaTransport.C:1234-1260`.
  All species' `nEqn()`s are built together into
  `List<autoPtr<fvScalarMatrix>> eqns(nSpecies)`
  (`plasmaTransport.C:1095-1153`) then solved in a flat loop
  (`plasmaTransport.C:1266-1286`) -- not nested per-species Picard.
  `validation/grubert2009` transports 3 species (`e`, `Ar2p`, `Arp`) on a
  2000-cell mesh.
- **Electron energy**: `localEnergyEnergyModel::eEqn()`
  (`localEnergyEnergyModel.C:1312+`) assembles `fvm::ddt(nEps_) +
  <convection/diffusion, same SG/CompleteFlux structure> + fvm::Sp(Lsrc,
  nEps_) + <Newton-linearised chemistry source term> ==
  <Joule (+ Coulomb) heating> + <chemistry source>`. Its convective flux
  and Joule heating depend directly on the CURRENT `ePotential`-derived
  field and electron density; species mobility/diffusivity (in LMEA mode)
  are tabulated against the energy equation's own output -- the coupling
  is genuinely three-way (Poisson <-> species <-> energy), not a one-way
  chain.
- **The actual outer sequence**, `soPlasmaFoam.C`'s `while (pimple.loop())`
  (explicit-Poisson branch, `soPlasmaFoam.C:463-522`): `em->solve()` ->
  `transport.solve(...)` (all species) -> `energy->solveSpeciesEnergy()`
  -> `species.updateChargeDensity()` -> `transport.updateSurfaceCharge()`.
  `plasmaOuterRelaxation` (Aitken) is invoked FROM INSIDE the species and
  energy solves (`plasmaTransport.C:1311-1314`, `plasmaEnergy.C:156-167`),
  auto-triggering a JOINT relaxation across both once every enrolled field
  has contributed in that corrector. This whole sequence (with relaxation
  folded in) is exactly the "one sweep" the two-field proof of concept's
  Gauss-Seidel preconditioner is designed to generalize to.
- **Numerics exposure**: `system/plasmaSimulationControls`, an
  `IOdictionary` re-read every step (`MUST_READ_IF_MODIFIED`), read by
  `plasmaTimeControl::read()` (`plasmaTimeControl.C:75-109`) and containing
  sibling sub-dicts `plasmaTimeControl {...}`, `outerCoupling {target;
  adaptiveRelaxation; tolerance; maxCorrectors; onNonConvergence;}`,
  `dischargeCurrent {...}`, `poisson {scheme; EScheme;
  nNonOrthogonalCorrectors;}`. A new `outerSolver picard | newton;` key
  fits naturally as a sibling of `outerCoupling`, consistent with how
  `poisson.scheme` already selects between numerics variants as a plain
  dict key -- confirming the exposure decision made earlier in this
  document without having seen the actual dict structure.

**Next actual step, not yet started:** wire the REAL 5-field system
(`ePotential` + `e`/`Ar2p`/`Arp` + `nEps`, `grubert2009`'s ~2000-cell mesh,
~10,000 total DOF in serial) into this now-validated block/FormFunction/PC
machinery. The remaining genuine unknowns are physics-assembly questions,
not block-coupling questions: writing each equation's `fvc::` residual
explicitly (mechanically translating each existing `fvm::`-based `nEqn()`/
`eEqn()`/Poisson assembly into an explicit expression -- SG/CompleteFlux
schemes will need their own explicit `fvc::` equivalents, not yet written
for either), and whether the preconditioner's one-sweep-per-application
cost (a full existing Poisson+species+energy solve, per GMRES iteration,
possibly per Newton iteration multiple times) is affordable at the real
mesh/timestep scale -- not yet measured.

## Phase B: the real 5-field residual, built, crashed twice, fixed, FIRST SUCCESSFUL RUN -- 2026-09-09

Wired directly into `soPlasmaFoam` as a genuine runtime-selectable
`plasmaNewtonSolver` interface (`src/numerics/newtonSolver/`, PETSc-free,
always built) with a concrete `snesNewtonSolver` ("type SNES") implementation
in a SEPARATE, optionally-loaded library
(`src/numerics/newtonSolverPETSc/`, `libplasmaNewtonSolverPETSc.so`) --
mirroring how `petsc4Foam` is already loaded (case declares it in
`controlDict`'s `libs`, absent by default, `soPlasmaFoam` itself stays
PETSc-free). Selected via `outerCoupling.outerSolver picard | newton` +
`outerCoupling.newtonSolver { type SNES; }`, replacing the whole segregated
Picard sequence in BOTH the semi-implicit and explicit Poisson branches of
`soPlasmaFoam.C`'s `while (pimple.loop())`.

**Seven small, additive exposures were needed in production classes to make
the real 5-field residual buildable at all**, each following the exact same
pattern as the two in Phase A (promote a value already computed somewhere
into a public accessor/wrapper, zero logic change):
- `plasmaTransportModel::mu()`/`D()` (virtual, FatalErrors by default --
  the Newton path supports only `driftDiffusion` so far, checked and
  FatalError'd clearly at construction if any species uses something else).
- `plasmaTransport::transportModel(i)` (const AND mutable overloads --
  mutable needed to call `.correct()` and refresh mu/D from a new trial
  density).
- `plasmaTransport::refreshChemistrySources()`, a thin wrapper for the
  private `mechanismSourceTerms()`.
- `plasmaTransport::refreshChemistryTimestepState()` -- see the crash below;
  this one is NOT thin, it's the genuine per-timestep chemistry-state
  preamble extracted verbatim out of `solve()`.
- `localEnergyEnergyModel::muEpsEff()`/`DEpsEff()` (already
  energyFactor-scaled, exactly what `eEqn()` itself uses for `phiEps`/`DEps`).
- `singleRegionPoisson::updateDerivedFields()` promoted from private to
  public (pure access-level change) -- refreshes E/Emag/phiE/reducedE from
  a newly-written trial `ePotential`.

**Two real crashes, both root-caused with evidence, not guessed:**

1. **SEGV inside chemistry, found via `gdb -batch -ex run -ex "bt full"`**
   (no debug build of PETSc needed -- gdb's symbol-offset attribution on
   the crashing frame, `<...computeChemistrySources...>-49139`, was enough
   even though the frame LABEL gdb printed, `mechanismSourceTerms`, was
   misleading). Root cause: `plasmaTransport::solve()` has its own "step 0"
   preamble (chemN0_/chemExt_/chemExtPrev_/chemExtSlope_/chemExtEps_/...
   sizing, guarded on `chemTimeIndex_ != mesh_.time().timeIndex()`) that
   `mechanismSourceTerms()`'s Option-4 stiff-chemistry path depends on --
   and the Newton path never calls `solve()` at all, so those arrays were
   never sized (empty `List<scalarField>`s), and indexing into them
   segfaulted. The "recomputing the source each outer iteration is
   idempotent" comment on `mechanismSourceTerms()` turned out to describe
   idempotency ACROSS OUTER ITERATIONS WITHIN `solve()`'s own already-set-up
   context, not standalone callability bypassing `solve()`'s preamble
   entirely -- a real, non-obvious distinction this session's design missed
   on the first pass. Fixed by extracting that preamble verbatim into
   `plasmaTransport::refreshChemistryTimestepState()` (private, called from
   both `solve()` and the new `refreshChemistrySources()` wrapper, which
   calls it first every time -- safe/idempotent by the SAME guard).
2. **A real, separate bug caught by the SAME crash**, in this session's own
   new code: `eqns[s].reset(new fvScalarMatrix(transport.transportModel(s).nEqn()))`
   copy-constructs a NEW matrix from a `tmp<fvScalarMatrix>`'s dereferenced
   content -- the real code (`plasmaTransport.C:1152`) uses
   `eqns[i].reset(transportModels_[i].nEqn().ptr())` instead, transferring
   ownership out of the `tmp<>` rather than copy-constructing from it. Not
   the actual cause of the crash above (fixing it alone did not stop the
   SEGV), but a genuine defect fixed alongside it.
3. **Missing `fvc::flux()` in the explicit residual itself**, caught by a
   `FOAM FATAL IO ERROR` (not a crash) after the SEGV was fixed:
   `fvc::div(phi)` on a bare carrier flux (`Z*mu*phiE`, matching the
   REAL code's `phi`/`phiEps` variable, which is likewise not yet the total
   flux) silently omits the `n`/`nEps` multiplication that `fvm::div(phi,n)`
   performs internally -- there is no such thing as a genuine explicit
   equivalent of `fvm::div(phi, n)` via a bare `fvc::interpolate`+`fvc::div`
   pair. `fvc::flux(phi, n, schemeName)` is the actual OpenFOAM primitive
   for this (interpolates `n` to faces with the SAME scheme `fvm::div`
   would use, forming the real total flux, which `fvc::div` then
   integrates) -- used for both the species and energy convection terms
   now. Also needed explicit scheme-name overrides in three places
   (`interpolate(mu_e)`, `div(phi_e,n_e)`, `laplacian(D_e,n_e)`) for the
   energy equation specifically, exactly matching `eEqn()`'s own borrowed-
   scheme convention (`muEps`/`DEps` have no `fvSchemes` entry of their
   own) -- these are silent-if-missed correctness bugs (a case with a
   matching wildcard entry would not have caught them), not silent-if-
   missed crashes, so worth flagging as a class of risk for the remaining
   physics-translation work.

**Result: first fully successful run of the real 5-field Newton solver on
`grubert2009_pseudo`**, no crash, no FatalError, running cleanly through
80+ PIMPLE outer iterations (t = 1e-12 to ~8e-11 s) before being left
running longer in the background to see whether it reaches real (non-
trivial) Newton iterations as the discharge develops, and whether it
remains stable through the region
([[outer-coupling-diverges-at-negative-differential-resistance]], ~step
23,000) where Picard was shown to diverge. `SNESConvergedReason = 2`
(`FNORM_ABS`), 0 Newton iterations, on every step so far -- plausible at
this very early, still-quasi-static point in the discharge (S_iz~0,
meanE constant) rather than a bug, but not yet confirmed either way at
later times. Not yet done: the physics-based preconditioner for the real
system (Phase C, still bare/unpreconditioned JFNK here, matching the
proof-of-concept's own staged validation order), and re-verifying the
`plasmaTransport` refactor left the existing Picard path bit-identical
(deferred until this background run finishes, to avoid modifying
`plasmaSimulationControls` -- `MUST_READ_IF_MODIFIED` -- out from under a
live run reading it).

## Phase C: preconditioner built, and the TWO defects that made Newton meaningless until fixed -- 2026-09-09

The `SNESConvergedReason = 2, 0 iterations` reported above is now
understood, and it was NOT "plausible at an early quasi-static point". It
was a FAKE convergence, and finding out why exposed two separate real
defects. Both are fixed or worked around; the residual itself turned out
to be correct.

### The six build/wiring bugs, fixed first

Phase C's Gauss-Seidel preconditioner (Poisson -> species -> energy, each
solving a linearised CORRECTION against a synthetic GMRES right-hand
side) needed six fixes before it ran at all. Recorded because four of
them are traps any future block preconditioner in this codebase will hit:

1. **Correction-field DIMENSIONS.** `bPoisson`/`bSpecies`/`bEnergy` were
   built `dimless` while the equations they feed use REAL dimensioned
   coefficients (`epsilon`, `mu`, `D`). Fixed by giving each block's RHS
   the dimension its own residual actually carries -- Poisson's is
   `chargeDensity`'s, a transported block's is `field/dimTime`. `chemL`
   needed `dimless/dimTime` (a rate coefficient) so `fvm::Sp(chemL, dn)`
   matches `fvm::ddt(dn)`.
2. **Correction fields must NOT inherit the real field's BC OBJECTS.**
   Copying `species.numberDensity(s)` wholesale gave `d_e` the real
   `ddWallFluxMixed` boundary condition, whose `updateCoeffs()` resolves
   its species from the FIELD'S OWN REGISTERED NAME (`n_e` -> `e` by a
   prefix strip). On `d_e` that threw `d_e not found in table. Valid
   entries: 3(Arp e Ar2p)`. A correction has no wall-flux physics of its
   own: only its Dirichlet/Neumann STRUCTURE must match. Fixed with
   `homogeneousPatchTypes()` -- `fixedValue`-0 where the real patch
   `fixesValue()`, `zeroGradient` elsewhere, and never the real BC class.
3. **Scheme-name auto-derivation breaks on a renamed field.**
   `fvc::snGrad(dePotential)` derives `snGrad(d_ePotential)`, absent from
   `fvSchemes`. Every operator on a correction field needs the REAL
   field's key passed explicitly (`"snGrad(ePotential)"`,
   `"laplacian(epsilon,ePotential)"`, `"laplacian(D_<sp>,n_<sp>)"`).
4. **The energy equation's solver dict.** Hardcoding `subDict("nEps_e")`
   aborted: cases carry a regex `"n_.*"` that `nEps_e` does not match.
   The real code (`plasmaEnergy.C`) falls back to `n_<species>`'s
   settings, and the PC now replicates that same fallback rather than
   inventing a requirement on every case.
5. **`PETSC_OPTIONS` never reached PETSc AT ALL.** `PetscInitialize()`
   runs in `snesNewtonSolver`'s CONSTRUCTOR, which consumes the
   environment variable once, at program start -- long before
   `solveOuterStep()` set it. So every `rtol`/`maxIt` this project
   configured was silently ignored and PETSc ran on its own defaults
   (notably `ksp_max_it` 10000) for the whole of Phase B. Fixed by
   `setPetscOptions()` in `snesBridge`, which calls
   `PetscOptionsInsertString` into the live options database before each
   solve, where `SNESSetFromOptions()` picks it up.
6. A `tmp<>` ownership bug and a leftover nonsensical ternary, both
   self-caught.

### DEFECT A: the density clamp sits OUTSIDE the equations, so F=0 is unreachable

`validation/grubert2009_pseudo` starts `internalField uniform 0` for
EVERY field -- `n_e`, `n_Arp`, `n_Ar2p`, `nEps_e`, `ePotential`. The
density floor is a per-step CLAMP, not an initial value
([[density-floor-is-a-source]]). Measured consequence, with a per-block
residual-vs-terms probe:

    step 1 (t=1e-12): every term exactly 0 -> ||F|| = 0.000000000000e+00
                      -> SNES reports CONVERGED_FNORM_ABS, 0 iterations

That is the fake convergence: SNES converged on the ZERO state and did
nothing. `clampNumberDensities()` then jammed densities 0 -> 1e11, a
state change made outside the equations, and at the next step:

    n_e   : |ddt|=2.2361e24  |div|=1.72e21  |lap|=1.14e22  |F|=2.2354e24
    n_Arp : |ddt|=1.2496e19  ...                            |F|=2.2361e24
    nEps_e: |ddt|=8.671e22   ...                            |F|=8.667e22

`||F||` equals `||ddt||` to 0.03% in every transported block, with
`div`/`lap`/`chem` 100-1000x smaller. **The residual was ENTIRELY the
time-derivative of the clamp's own 0 -> 1e11 jump**, which nothing in the
equations can balance. That is why 18,400 residual evaluations bought one
Newton iteration and a factor-of-2 reduction.

A clamp applied after the solve means the discrete system being solved is
not `F(u) = 0`, so no Newton method can converge to it. Two responses:

* **Chosen now (interim):** do not cold-start Newton at all. Run Picard
  to an established state and restart Newton from it. Verified on the
  2e-9 restart: `n_e` 1.331e11-1.5e11, `n_Arp` 1.334e11, `nEps_e`
  1.97e11, and **0 cells pinned at the floor** for any species -- the
  clamp is inactive, so the residual is genuinely satisfiable.
* **Deferred (the real fix):** solve for `log(n)`, making `n = exp(psi) >
  0` structural so no clamp is needed. Standard for drift-diffusion, and
  Knoll & Keyes cite Gummel's semiconductor method in exactly this role
  (section 3.4). Touches the species residual, the preconditioner AND the
  real Picard path, so it is a project, not a patch. User asked for it to
  be kept on the list 2026-09-09 -- see [[deferred-action-items]].

### DEFECT B: unscaled blocks made SNES structurally blind to Poisson

With every block in its own physical units, one combined L2 norm is
meaningless. Measured at the same state:

    Poisson: |lap|=1.247e-3  |rhs|=5.85e-8   -> block ||F|| ~ 1.2e-3
    species: ||F|| ~ 2e24
    combined SNES function norm = 3.163016002050e+24

The combined norm is 100% species. With `rtol 1e-8` the convergence
target is 3.16e16 -- **nineteen orders of magnitude ABOVE the entire
Poisson residual**, which sits 27 decades below the norm's noise floor.
SNES could not see the Poisson equation at all, so "converged" said
nothing about it. This is precisely the per-component `typ u` scaling
Knoll & Keyes call for in section 2.3.1 (eq. 13) -- read in the paper
itself (`Literature/Knoll_1-s2.0-S0021999103004340-main.pdf`), not
inferred.

**Fix: per-block scaling, computed ONCE per `solveOuterStep`.** Fixed for
the whole Newton solve -- recomputing per residual call would move the
target and the problem would stop being stationary.

    x_petsc = x_phys / sX      F_petsc = F_phys / sF

    sX[block] = rms(that block's own field)
    sF[transported] = sX/dt            (ddt, the measured dominant term)
    sF[Poisson]     = rms(lap(effEps, ePotential))   (ditto)

Applied at four sites: packing `x`, unpacking `x` in the residual, writing
`F`, and in the PC -- where GMRES hands `b` in scaled RESIDUAL units and
expects `y` in scaled STATE units, so `b_phys = b*sF`, solve, `y =
d_phys/sX`. Each scale is floored (`safeScale`) so a cold all-zero field
cannot divide by zero.

### Result: Newton actually converges, and the residual is VALIDATED

From the 2e-9 Picard restart, with both defects addressed:

    0 SNES Function norm 5.287294435400e+02
    1 SNES Function norm 4.394442751829e+01     (12x)
    2 SNES Function norm 1.834287727213e-05     (2.4e6 x)

Superlinear, i.e. genuine Newton behaviour, and `pcApplyCallback` is now
invoked -- 4-5 GMRES iterations per Newton step, the flat low
GMRES/Newton ratio Knoll & Keyes report as this preconditioner class's
success signature (section 5.4, tables 4-5: 4-26 across grid
refinements).

**The residual is validated against ground truth**, which was the actual
question. At the Picard-converged state the scaled species residuals are
`n_e` 1.2e-3, `n_Arp` 2.6e-4, `n_Ar2p` 7.5e-3 -- small, as they must be
if the `fvc::` residual reproduces the equations Picard solved. Confirmed
independently by a temporary probe placed inside
`singleRegionPoisson::solve()` itself, immediately after
`ePotentialEqn.solve()`, evaluating the SAME `fvc::` expression on
`fvm::`'s just-produced solution:

    rms|fvc::lap|=1.3077e-09  rms|rhs|=1.3077e-09  rms|fvcRes|=3.70e-13

So `fvc::laplacian` does reproduce `fvm::laplacian`'s action, boundaries
included. (Probe REMOVED again -- it was in production code.)

**Two wrong explanations recorded, because both were stated before being
checked** (rule 22): the Poisson residual's apparent 21,300x mismatch was
first blamed on `effEps` (excluded by measurement: `rms(effEps)`
8.854190e-12 vs `epsilon` 8.854188e-12, identical to 6 digits, and
`sigma` only 2.9e-6 S/m), then on a circuit boundary update (excluded:
the cathode patch value in the Newton path is `-0.11415401995`, identical
to the restart file). The measurement that actually explains it:
`max|lap|` = 5.576e-4 at cells 0 and 2, `x = 6.77e-7` m -- the first cell
column against a Dirichlet boundary -- and those few cells carry the
whole norm. The cathode FACE is at -0.11415401995 while the first cell
CENTRE is at -0.114088468576, a 6.56e-5 V drop giving a near-wall field
of 96.8 V/m against 11.8 V/m in the interior; `effEps *
(96.8-11.8)/1.353e-6 = 5.6e-4`, reproducing the measured value exactly.
That 6.56e-5 V is the same order as the circuit's own `dV/step` of
-5.752e-5, so the leading explanation is that the written cathode BC is
ONE circuit increment ahead of the written interior field, making the
reloaded pair inconsistent in the first cell column
(cf. [[dt-ordering-clock-lead]]). NOT PROVEN -- stated as the leading
explanation only. It is benign for the solver either way: Newton removes
it in the first iteration, and write precision is excluded
(`writePrecision 12`, resolution ~5e-13, induced laplacian error ~2.5e-12
against a measured 5.6e-4).

### OPEN: the residual is too expensive -- ~30-60 s per timestep

Measured 2026-09-09: ONE Newton timestep takes 30-60 s (1 timestep in a
60 s window; the 90 s window completed 2 Newton iterations). Picard on the
same case runs 44 steps/s. Each residual evaluation does FOUR full matrix
discretisations (three species `nEqn()` plus the energy `eEqn()`, all
built and discarded purely for their side effects) AND a per-cell stiff
chemistry ODE integration, and JFNK needs one per GMRES iteration. Not
yet profiled -- the chemistry ODE is the suspect, not the conclusion.

Knoll & Keyes section 5.1 ("Jacobian lagging") is the relevant prior art:
expensive pieces can be lagged in the preconditioner while the true
operator stays in the residual. Note the asymmetry here -- lagging
chemistry in the RESIDUAL changes which system is being solved (it
becomes chemistry-lagged/semi-implicit, not fully coupled), whereas
lagging it only in the PRECONDITIONER costs nothing in the converged
answer. That distinction should decide the fix.

### Also still open

* **`species.updateChargeDensity()` is never called in the residual.**
  So `chargeDensity` does not follow the trial species densities, and the
  Poisson block's coupling to the species is MISSING from the Jacobian.
  Found by inspection while tracing the above, not yet fixed.
* The energy block's residual scale uses `ddt`, but at the restart state
  its dominant term is the LAPLACIAN (`|lap|` 1.2e26 vs `|ddt|` 8.4e24),
  leaving `|F_scaled|` ~ 527 for that block against ~1e-3 for species.
  A scale-choice imperfection, not a bug.
* Chemistry accuracy under Newton mode is UNVERIFIED: the
  `adaptiveError` `nOuterCorrectors > 1` guard is bypassed for
  `outerSolver newton` only (Picard's check is untouched and still
  fires). The user asked for simple tests of chemical accuracy here.
* The `plasmaTransport` refactor has still not been re-verified as
  bit-identical for the existing Picard path.
* Phase D -- the run past Picard's ~step 23,000 failure point -- has NOT
  been attempted. The per-step cost above must come down first.
