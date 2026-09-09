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

## The per-step cost RESOLVED: 97.88 s -> 0.195 s/step, and Newton now matches Picard's trajectory -- 2026-09-09

The "residual is too expensive" conclusion recorded above was WRONG about its
cause, and profiling said so immediately. SUPERSEDED BY the two findings here.

### The residual was never the problem: 3 ms, called 1500x too often

Per-phase timers inside residualCallback, one Newton timestep:

    calls=30020  total=92.09 s  | unpack=0.17 derived=2.27 transport=4.62
                                  CHEM=41.11 energy=23.59 residual=20.34
                                  (per call: 0.00307 s)

So a residual evaluation costs 3.07 ms -- cheap. The 92 s came from THIRTY
THOUSAND AND TWENTY of them in ONE timestep, against the ~20 that the actual
iteration counts justify (SNES took 3 Newton iterations; its KSP solves
converged in 3, 4 and 7 GMRES iterations). Chemistry is only ~45% of a call,
so lagging it -- the fix the section above proposed -- would have bought under
2x while the real defect was a factor of 1500.

30020/3 ~ 10,000 per Newton iteration is exactly the number of unknowns
(5 fields x 2000 cells): the signature of finite-differencing a whole Jacobian
COLUMN BY COLUMN.

### DEFECT 1: `mf_operator` made PETSc assemble a Jacobian it never used

`snesBridge.C` had

    SNESSetUseMatrixFree(snes, PETSC_TRUE, pcCallback ? PETSC_FALSE : PETSC_TRUE);

with a comment asserting the second argument controlled whether PETSc derives
its own preconditioning matrix. The argument SEMANTICS were backwards. From
PETSc's own source (`src/snes/interface/snes.c`):

    snes->mf = mf_operator ? PETSC_TRUE : mf;   // mf_operator forces mf on
    snes->mf_operator = mf_operator;

and its doc string: `mf_operator` = matrix-free "only for the Amat ... this
means the user provided Pmat WILL CONTINUE TO BE USED"; `mf` = matrix-free for
both, "both the Amat and Pmat set in SNESSetJacobian() will be ignored".
`SNESSetUpMatrices()` takes its matrix-free-for-both branch ONLY when
`snes->mf && !snes->mf_operator`; with `mf_operator = TRUE` it takes the other
branch, which calls `DMCreateMatrix()` for a real Pmat and fills it with
`SNESComputeJacobianDefault` -- one residual call per column. We never provide
a Pmat and the PCSHELL ignores it, so the entire assembly was waste. Note both
spellings of the old ternary set `mf_operator = TRUE`, so bare JFNK was doing
this too.

Replaced with the idiom PETSc's documentation prescribes for a completely
matrix-free solver where Pmat IS the operator:

    Mat Jmf;
    MatCreateSNESMF(snes, &Jmf);
    MatSetFromOptions(Jmf);
    SNESSetJacobian(snes, Jmf, Jmf, MatMFFDComputeJacobian, nullptr);
    MatDestroy(&Jmf);

### DEFECT 2: the differencing step assumed F was evaluable to machine epsilon

Removing defect 1 dropped 30,020 residual calls to ~20 -- and then GMRES hit
`DIVERGED_BREAKDOWN` at 30 iterations. `-ksp_monitor_true_residual` separated
the causes in one run:

    iter   preconditioned resid    true resid      ||r||/||b||
       0   2.332e+01               5.287e+02       1.00
       1   5.074e-01               1.055e+03       2.00
       3   2.235e-02               2.427e+09       4.59e+06
      29   4.305e-03               2.431e+09       4.60e+06

The PRECONDITIONED residual falls cleanly (23 -> 0.004), so the physics-based
PC works; the TRUE residual explodes 4.6 million-fold and flatlines, so the
Jacobian-vector product does not represent the true operator. Confirmed not to
be the preconditioner: the breakdown is IDENTICAL under PCNONE (verified via
`-pc_type none`, with zero pcApplyCallback invocations).

Cause: PETSc's default differencing parameter is built from sqrt(eps_mach),
which assumes F can be evaluated to machine precision. F here integrates a
PER-CELL ADAPTIVE STIFF CHEMISTRY ODE, whose internal step sequence changes
discontinuously with its input, so F carries noise many orders above eps_mach
and differencing it yields garbage. Knoll & Keyes state the rule in section
2.3.1: when F is only evaluable to eps_rel, the differencing parameter must be
built from eps_rel, not eps_mach.

Measured sweep of `-mat_mffd_err`, one Newton step from the 2e-9 restart:

    ~1.5e-8 (default)   DIVERGED_BREAKDOWN, true resid 2.4e9
    1e-6                DIVERGED_BREAKDOWN
    1e-4                CONVERGED: 528 -> 63.6 -> 8.76e-3 -> 3.33e-9 (3 its)
    1e-2                CONVERGED: 528 -> 73.9 -> 3.79e-4 -> 2.53e-10 (3 its)

So F is good to roughly 1e-4. Exposed as `newtonSolver/mffdErr`, DEFAULT 1e-4,
so no case has to know this (G1). The reasoning is recorded on the member's own
declaration in `snesNewtonSolver.H`.

### Result

    per timestep      residual calls   residual time   wall clock
    before                    30,020         92.09 s      97.88 s
    after                         16          0.05 s       0.195 s

20 consecutive steps from the 2e-9 restart: 320 residual calls (16/step),
2 Newton iterations per step, every step `reason 3`
(CONVERGED_FNORM_RELATIVE). Newton 7.07 s against Picard 1.92 s for the same
20 steps including start-up (Newton's start-up carries PetscInitialize), i.e.
~3.7x -- an ordinary price for a fully coupled solve, against 4000x before.

### THE PHYSICS GATE: Newton reproduces Picard's trajectory

Paired 20-step runs from the SAME 2e-9 restart, identical in every parameter
except `outerSolver` (rule 15: the control is stated). Relative difference of
the written fields at t = 2.02e-9:

    field           max|rel diff|   rms rel diff
    ePotential         1.075e-05      6.267e-06
    n_e                4.755e-06      1.134e-06
    n_Arp              7.459e-07      5.164e-08
    n_Ar2p             5.499e-07      3.852e-08
    nEps_e             8.339e-08      2.046e-08
    chargeDensity      8.197e-05      1.974e-05
    meanE / T_e        4.869e-06      1.161e-06

Agreement at the level of the two solvers' own tolerances (Picard outer 1e-6,
SNES rtol 1e-8). `chargeDensity` is the loosest, as expected: it is a small
difference of near-equal large numbers in a quasineutral plasma. This is the
"does the real-physics translation actually work" gate from the phase plan, and
it PASSES.

### Still open

* Phase D -- the run through and past Picard's ~step 23,000 failure point --
  is now UNBLOCKED and not yet attempted.
* `species.updateChargeDensity()` is still not called in the residual, so the
  Poisson block's coupling to the species is missing from the Jacobian. Newton
  converges regardless (2 iterations/step) because the PC supplies that
  coupling, but the Jacobian is not the true one. Fix before Phase D.
* The energy block's residual scale still uses `ddt` where the laplacian
  dominates; `|F_scaled|` ~ 527 for energy against ~1e-3 for species.
* Chemistry accuracy under Newton remains UNVERIFIED (the `adaptiveError`
  `nOuterCorrectors > 1` guard is bypassed for newton mode only).
* Temporary per-phase timers and per-block residual diagnostics are still in
  `residualCallback`; remove once Phase D is under way.

## Two CORRECTNESS defects in the residual, and a retracted measurement -- 2026-09-09

Found by a four-agent code study run against the working tree. Both defects are
fixed; one of my own measurements from the section above is RETRACTED.

### 0a. The residual never refreshed chargeDensity

`residualCallback` read `em.chargeDensity()` for the Poisson block but never
called `species.updateChargeDensity()` -- whose only call sites were
`soPlasmaFoam.C` (the Picard branches) and `plasmaTransport.C:550`. Verified by
grep, not inferred.

This is worse than a Jacobian defect. It corrupts F ITSELF, so Newton converged
to the root of a Poisson-with-LAGGED-source system -- exactly the segregated
coupling this solver exists to remove at NDR, with the enclosing
`pimple.loop()` supplying the outer lag. **It also means the "physics gate
passed" result above was not evidence of correctness**: agreeing with Picard is
precisely what a Newton solver that is secretly solving Picard's system does.
Fixed: `updateChargeDensity()` now runs right after the species are unpacked.

### 0b. meanE was refreshed AFTER the lookups that key on it

`meanE_` is recomputed only in `localEnergyEnergyModel::updateDerived()`,
reachable only from `correct()`. That call sat AFTER
`transportModel(s).correct()` and `refreshChemistrySources()` -- both of which
look up coefficients and reaction rates keyed on the `"meanE"` field. So F was
evaluated with the PREVIOUS residual call's mean energy: F = F(u_k, u_{k-1}),
not a function of u at all, which corrupts every matrix-free difference
quotient. Fixed by moving `correct()` ahead of the transport and chemistry
refreshes, while leaving `eEqn()`/`updateSources()` after them (Psrc_/Lsrc_
depend on the chemistry sources under `energySource chemistry`).

### RETRACTED: the claim that the ordering fix made F smooth

The section above recorded a sweep showing 1e-8, 1e-6 and 1e-4 all converging
"to the IDENTICAL final norm 2.247181640213e-10", read as insensitivity to the
differencing step and therefore as proof the operator was now correct.

THAT MEASUREMENT WAS INVALID. `mffdErr_` still defaulted to 1e-4 at the time, so
`setPetscOptions()` inserted `-mat_mffd_err 1e-4` AFTER the swept
`PETSC_OPTIONS` value and overrode it. All three runs were the same
configuration -- which is exactly why they agreed to the last digit. Rule 22:
the diagnostic was the first suspect and was not checked.

The valid measurement, with the default unset so the environment variable takes
effect, 20 steps from the 2e-9 restart:

    1e-8, 1e-7, 1e-6   fail at step 1, DIVERGED_BREAKDOWN
    1e-5               20/20 steps converged, <= 3 Newton iterations
    (PETSc default PETSC_SQRT_MACHINE_EPSILON ~ 1.49e-8 also fails)

So the ordering fix moved the usable threshold from ~1e-4 to ~1e-5 -- real, but
NOT a qualitative cure. **F is still not smooth at the 1e-6 level.** Two
identified, still-unfixed causes:

* `updateDerived()` writes a CLAMP back into `nEps_`, so the state F is
  evaluated at is not the state PETSc perturbed. A clamp inside the residual is
  the same class of defect as defect A above, one level down.
* species/`nEps` `correctBoundaryConditions()` run before `updateDerivedFields`
  and the mobility refresh, so `ddWallFluxMixed` (which reads live `phiE` and
  patch mobility) is one call stale.

`newtonSolver/mffdErr` therefore ships as a documented WORKAROUND, default
1e-5, with the reasoning on its declaration. A case that diverges should raise
it first.

### What the 20-step window can and cannot show

After both fixes, 20/20 steps converge (19 at 2 Newton iterations, 1 at 3;
8.13 s against Picard's 1.92 s). Three-way field comparison at t = 2.02e-9:

    field           fixed-vs-picard   fixed-vs-unfixed-newton
    ePotential            1.08e-05                  1.73e-06
    n_e                   4.76e-06                  3.05e-09
    nEps_e                7.85e-08                  4.89e-09
    chargeDensity         8.20e-05                  5.36e-08

The correctness fix moved the answer by ~1e-6 at most. **This window therefore
does NOT exercise what was fixed** (rule 23): over 1 ps in a near-quasi-static
region, a lagged Poisson source and a coupled one barely differ. The lag matters
at the negative-differential-resistance point, which is Phase D. Do not cite
these numbers as validation of the coupling.

### Remaining, in priority order

1. Phase D: run through Picard's ~step 23,000 divergence. This is the only test
   that exercises the coupling the fixes restored.
2. Remove the two remaining non-smoothness sources above, then check whether
   `mffdErr` can drop toward PETSc's default. That is the objective measure of
   whether F is finally a clean function of u.
3. VERIFY, not yet done: `plasmaTransport.C:4306-4327` mutates `chemSrcPrev_`
   and `chemOuterCount_` DURING residual evaluation. If anything consumes
   `chemPicardChange_` in Newton mode, the residual is history-dependent and
   that is a Tier-0 defect, not a performance note.
4. Correctness-free speedups, all identified and none applied: drop the three
   discarded `nEqn()` builds (every `eqns` access in `mechanismSourceTerms` is
   already guarded, so an empty list is legal today); drop the `eEqn()` build by
   exposing `updateSources()`; cache the preconditioner's coefficient fields per
   Newton step (~6 s of the 8 s run is PC+KSP, now the largest target); and get
   the MPI collectives and `Info<<` out of the matvec
   (`plasmaReactionRates.C:478-490`, `plasmaTransport.C:4075`,
   `driftDiffusion.C:387`).
5. `snesNewtonSolver` should FatalError, not silently solve a different system,
   on configurations the residual ignores: surface charging, photoionisation,
   and the legacy Townsend `!rates_` source. Verified absent for grubert2009
   ONLY -- re-check before offering Newton for needleDBD or the streamer cases.
6. Chemistry accuracy under Newton is still UNVERIFIED.
7. Remove the temporary per-phase timers and per-block residual diagnostics.

## Cold-start guard added, and PHASE D: Newton FAILS at a developed discharge -- 2026-09-09

### The guard (done)

`solveOuterStep` now refuses an all-zero start with an actionable message
instead of reporting CONVERGED after solving nothing. Verified BOTH ways
(rule 23 / [[silent-diagnostic-trap]]): it fires on `grubert2009` from t=0, and
stays silent on the 2e-9 restart, which still runs 20/20 steps converged.

### Phase D control arm: Picard does NOT crash -- it stops CONTRACTING

Picard, t=0 -> 4e-8 (40,000 steps, dt 1e-12 fixed), 953 s: ran to completion,
no FatalError, no reported convergence failure. But `rho [contraction]`
crossed 1.0 at step ~18,878 and climbed to ~1321 ("residual GREW"), while
`omega` stayed 0.8-1.0 and PIMPLE kept running its fixed 4 correctors.

So the documented "divergence" is not a crash in this configuration: the outer
loop simply stops contracting and the clock keeps advancing with the coupling
unconverged. Per [[deferred-action-items]] `rho` is independently known to cry
wolf, so this is NOT called a divergence on `rho` alone -- but it does sharpen
Phase D's question to: at a state past step 18,878, does Newton actually
CONVERGE the coupled system where Picard's outer loop no longer contracts?

### Phase D Newton arm: FAILS, and it is OUR solver, not the physics

Newton restarted at t=2e-8 (step 20,000): the very first step fails,
`DIVERGED_LINEAR_SOLVE`. Per-block scaled residuals there:

    Poisson   44.7   (= sqrt(2000), i.e. 1.0 per cell)
    n_e        0.089
    n_Ar2p     0.0025
    n_Arp      0.00011
    nEps_e   172.5    <-- dominates; |ddt| 2.98e25 vs |lap| 1.08e25

Ruled out BY MEASUREMENT, not by inspection:

* **The preconditioner.** Under `-pc_type none` (0 pcApplyCallback calls) the
  failure is identical.
* **The differencing step.** mffdErr 1e-4, 1e-3, 1e-2 all fail identically --
  no effect whatsoever, unlike at t=2e-9.
* **Chemistry history dependence.** Freezing chemP_/chemL_ across the solve
  changed the KSP residual history by less than one part in 1e8.
* **The mean-energy clamp.** Bounds are [0.039, 2644] eV; meanE is 0.81-1.62 eV,
  nowhere near binding.
* **GMRES restart length.** `fgmres` gives DIVERGED_DTOL at 30, `bcgs` and
  `gmres_restart 200` give DIVERGED_ITS at 100. So DIVERGED_BREAKDOWN was an
  artifact of the 30-vector restart; the real situation is a linear system that
  will not converge.

What remains, and it is the smoking gun: **under PCNONE the preconditioned and
true residual norms MUST be identical, and they are not** --

    iter 0   preconditioned 178.16   true 178.16     (equal, as required)
    iter 1   preconditioned  35.96   true 348.57     (impossible if A is linear)

GMRES's Arnoldi recurrence and the explicit b - A*x disagree, which cannot
happen for a consistent LINEAR operator. With the preconditioner, the
chemistry, the clamp and the Krylov method all excluded, the indicated cause is
that the finite-difference Jacobian is not acting linearly over the
perturbations GMRES uses -- i.e. genuine nonlinear stiffness at a developed
discharge, the regime Knoll & Keyes address in section 3.6 (nonlinear
preconditioning / ASPIN) and section 2.4 (globalization). NOT yet proven by a
direct linearity test; stated as the indicated cause.

Note this is consistent with everything else: it worked at t=2e-9 (weak,
quasi-static, chemistry inactive -- chemP was exactly 0) and fails at t=2e-8
(sheath forming, chemP ~ 3.4e9).

### What COMSOL does, read from the User's Guide (their Plasma Module)

Read 2026-09-09 from `Literature/COMSOL_PlasmaModuleUsersGuide.pdf` because
the log(n) question is really "should we copy COMSOL here". Findings, quoted:

* **Their DEFAULT is the log formulation**: "Finite element, log formulation
  (linear shape function) (the default) to solve the equations in logarithmic
  form". They solve for ln(n_e) AND ln(n_eps).
* **Their rationale**: "the electron number density can span 10 orders of
  magnitude over a very small distance... The best way of handling this from a
  numerical point of view is to solve for the log of the electron number and
  energy density. This also prevents a divide by zero."
* **Their admitted cost, in their own words**: "This makes it more numerically
  stable but INCREASES THE NONLINEARITY of the equation system, and as such the
  model might take slightly longer to solve."
* **It does not remove the zero problem, it relocates it**: "the solver can run
  into difficulties when the species mass fractions approach zero", handled by
  a "Source stabilization" term `R_k,tot = R_k + exp(-iota*ln(n_k))` with a
  USER-TUNED parameter (default 1, "if the plasma is high pressure
  (atmospheric) then it can help to lower this number to 0.25-0.5"). A
  per-regime tuning knob is a direct G1 violation in our terms.
* Zero becomes inexpressible: "specifying an electron density of zero is not
  allowed."
* **Their log form is FINITE ELEMENT.** "Finite volume (constant shape
  function)" is offered as a SEPARATE, non-log option -- they do not ship
  log + finite volume. For an FV/OpenFOAM code that is both a warning and an
  opening.
* Their heavy-species log form divides through by rho*w_k (the
  NON-conservative form): flux `V_k = D grad(W_k) + D grad(lnM) + D^T grad(lnT)
  - Z mu E`.

### Why this changes the log(n) recommendation

Two independent sources now point the same way. COMSOL's own guide says the log
formulation INCREASES nonlinearity; our own measurement says nonlinearity is
exactly what is blocking us at t=2e-8. **So log(n) would likely make our
CURRENT blocker worse, not better.** It addresses positivity, which for us only
bites at cold start (now guarded), and not the developed-state failure.

Candidate alternatives, to be planned properly rather than adopted by default:

1. **Pseudo-transient continuation / dt as continuation parameter** (K&K 2.4.2).
   Directly targets the measured blocker, cheap to try, no reformulation.
2. **Bounded Newton via PETSc's variational-inequality solvers**
   (SNESVINEWTONRSLS): enforce n >= n_floor as a CONSTRAINT inside the Newton
   solve instead of by a post-hoc clamp or a change of variable. Keeps the
   conservative FV form, the Scharfetter-Gummel scheme, every existing BC and
   every diagnostic; removes the clamp-outside-the-equations defect; already
   available in the PETSc we link. Smallest footprint, and genuinely distinct
   from COMSOL.
3. **Lean on exponential fitting we already own**: Scharfetter-Gummel is
   already exponentially fitted and positivity-friendly, which is the same
   physics insight log(n) encodes -- but inside the FLUX, keeping the unknown
   conservative. Strengthening that is differentiating rather than imitative.
4. **Nonlinear preconditioning (ASPIN-family, K&K 3.6)** for the unbalanced
   sheath-vs-bulk nonlinearity that is our actual failure. Research-grade and
   publishable (rule 33).

## The Phase D blocker is CONDITIONING, not nonlinearity: dt is the lever, and SG makes it WORSE -- 2026-09-09

SUPERSEDES the "nonlinear stiffness" attribution in the section above. That was
too pessimistic and is retracted.

### dt IS the lever (option 1 confirmed)

Newton restarted at t=2e-8, ONLY deltaT varied, 20 steps each, `standard` flux
scheme:

    dt        result
    1e-12     FAIL at step 1 (DIVERGED_LINEAR_SOLVE)   <- Picard's dt
    9e-13     FAIL at step 1
    7.5e-13   FAIL at step 1
    5e-13     20/20 converged, 3 Newton iterations
    1e-13     20/20 converged, 3 Newton iterations
    1e-14     20/20 converged, 3 Newton iterations

Threshold between 5e-13 and 7.5e-13, so the penalty against Picard's 1e-12 is
only about 1.5-2x -- NOT the 10x first guessed.

WHY, and it is a straightforward conditioning argument borne out by the block
norms at that state: for the species blocks |lap| exceeds |ddt| by ~1000x
(n_e: |ddt| 1.5e19, |div| 4.4e21, |lap| 1.4e22). The ddt term contributes
1/dt = 1e12 to the diagonal while diffusion contributes D/dx^2 ~ 5e14, so at
dt = 1e-12 the Jacobian is NOT diagonally dominant. Halving dt doubles the
diagonal and the linear system becomes solvable. This is a DIFFUSION-stiffness
limit, not a failure of Newton's method.

Note the honest comparison: Picard RAN at dt=1e-12 while Newton REFUSES there.
Picard was not solving the coupled system at that dt either (rho ~ 1321, not
contracting) -- it advanced the clock regardless. Newton failing is Newton
reporting the truth. But operationally it costs ~2x the steps at this state.

### Scharfetter-Gummel makes it WORSE (option 3 refuted)

Tested because SG is exponentially fitted and positivity-friendly, i.e. it
encodes log(n)'s insight inside the FLUX while keeping the unknown
conservative. Measured, only the scheme varied:

    dt = 1e-12   SG        FAIL   (as standard does)
    dt = 5e-13   standard  20/20 converged
    dt = 5e-13   SG        FAIL at step 1

So SG NARROWS the usable dt envelope rather than widening it. On reflection
this is expected: SG's Bernoulli-function flux introduces exponential
dependence on the local Peclet number, which makes F MORE nonlinear and so
harder to finite-difference. SG remains the right choice for positivity and
monotonicity -- it is the wrong lever for JFNK conditioning. Verified SG was
actually active (all 3 species and the LMEA energy report
`fluxScheme ScharfetterGummel`, "Discretizing transport with SG scheme...").

### Consequence for the log(n) question

This strengthens the case AGAINST log(n) further. COMSOL's guide says the log
form "increases the nonlinearity of the equation system"; we have now MEASURED
that a different exponentially-fitted reformulation (SG) does exactly that and
costs us dt headroom. Two independent lines of evidence now say that adding
exponential structure to the residual hurts the thing that currently limits us.

### Where this leaves the options

1. **dt-aware Newton (do this).** The framework already limits dt by `Co_diff`
   and already has `retryStep` for rejected steps. The natural fix is to let a
   failed SNES solve REJECT the step and retry at smaller dt, exactly as the
   temporal-error controller already does -- rather than FatalError. That turns
   a hard failure into an automatic, documented dt reduction of ~2x at the
   stiffest moments and needs no new mechanism (rule 30).
2. **Bounded Newton (SNESVINEWTONRSLS)** -- still worth doing, but for the
   COLD START and to remove the clamp-outside-the-equations defect. Measured
   NOT to be the Phase D blocker: zero out-of-range table lookups in either
   run, and the mean-energy clamp is not binding (bounds [0.039, 2644] eV
   against meanE 0.81-1.62).
3. SG for JFNK conditioning: REFUTED above, do not pursue for this purpose.
4. ASPIN / nonlinear preconditioning: the motivation is weaker than thought,
   since the blocker is linear conditioning rather than unbalanced
   nonlinearity. Keep as a research track, not the next step.

## THE REAL BLOCKER FOUND: a one-shot model initialisation firing INSIDE the residual -- 2026-09-09

SUPERSEDES both earlier attributions. "Nonlinear stiffness" was wrong; "dt is
the lever / conditioning" was also wrong (and its dt sweep was non-monotonic at
a second state, which no conditioning argument can produce).

### The defect

`localEnergyEnergyModel::correct()` carries a ONE-SHOT LFA seed:

    if (seedFromLFA_) { seedFromLFA_ = false; nEps_ == meanE0*n_e; }

The Newton path's first `correct()` sat INSIDE residualCallback (put there by
the meanE ordering fix), so the seed fired on residual call #1 and on no other.
F was therefore a DIFFERENT FUNCTION at the base point than at every perturbed
point: the matrix-free product differenced F_{call>=2}(u+hv) against
F_{call1}(u), so J*v was meaningless.

**That is the "impossible" measurement recorded above** -- under PCNONE the
preconditioned and true residual norms diverged from iteration 1, which cannot
happen for a consistent linear operator. It was never stiffness and never
conditioning; the operator was not a fixed linear map because the function
underneath it changed after the first evaluation.

Measured at t=2e-8, nEps block scaled residual: **172.46 on call #1 against
0.00086 on call #2** -- a factor of 2e5, and ~97% of the reported initial norm
(178.16 = sqrt(44.72^2 + 172.46^2)). It also explains why t=2e-9 worked: there
the stored nEps is still ~the LFA equilibrium, so the seed is nearly a no-op,
while by t=2e-8 the state has moved 3.4x away.

### Fixed, and verified three ways

`solveOuterStep` now calls `lmea->correct()` ONCE before the solve, making the
seed part of the state the solve starts from -- which is what a one-shot
initialisation is for. After the fix:

    seed message now precedes residual call #1 (log line 176 vs 177)
    nEps |F_scaled| calls #1/#2/#3: 26.9182 / 26.9183 / 26.9182   (was 172.46 / 0.00086)
    initial SNES norm: 52.20                                       (was 178.16)

F is now consistent across calls. The t=2e-8 / dt=1e-12 case STILL fails with
DIVERGED_LINEAR_SOLVE, so this was one blocker and not the only one -- but it
was the one invalidating every diagnosis built on the residual.

### PRODUCTION BUG, not a Newton one -- needs its own fix (rule 37)

The seed is gated on

    freshStart = mesh.time().timeIndex() == mesh.time().startTimeIndex()

evaluated in the CONSTRUCTOR, before the time loop -- so it is ALWAYS true,
restart or not. The comment beside it states an intent the test cannot
implement. Consequence: **every LMEA restart with `initialMeanEnergy` silently
discards its stored energy state at the first correct(), on the PICARD path
too.** Confirmed in the setup log. Deliberately NOT patched from the Newton
code; it needs its own change and its own verification.

### Also fixed this round

* **Bounded Newton** (`newtonSolver/bounded`, default true): the species floor
  is now a CONSTRAINT inside the solve (`SNESVINEWTONRSLS` +
  `SNESVISetVariableBounds`) instead of a post-hoc clamp. Verified no
  regression on the working 2e-9 restart (20/20 converged, same iteration
  counts). The cold-start guard stands down when `bounded` is on.
* **PETSC_INFINITY overflows under FOAM_SIGFPE.** It is PETSC_MAX_REAL/4 ~
  4.5e307, so xu - xl ~ 9e307 and any further arithmetic overflows to inf,
  which OpenFOAM TRAPS. Use a finite 1e30 -- effectively unbounded for the
  O(1) scaled variables. (PETSc's own docs recommend PETSC_INFINITY; that
  advice does not survive FPE trapping.)
* **PCSHELL wrote past the end of its output vector under the bounded solver.**
  SNESVINEWTONRSLS is REDUCED-SPACE: it removes components at their bound and
  hands the KSP a SHORTER vector, while pcApplyCallback indexed the full field
  layout. Now detected via the callback's own `n` and falling back to the
  identity for restricted applications, with a warning. That is safe, not
  right: the proper route is an assembled Pmat with PCFIELDSPLIT, which
  SNESVI supports natively via PCFieldSplitRestrictIS.
* **BC ordering**: species/nEps `correctBoundaryConditions()` moved to AFTER
  phiE/meanE/mu are rebuilt. Structurally correct (the wall-flux BCs read live
  phiE and patch mobility), though it changed the base-point residual not at
  all here.

### Reviewed and CLEARED (do not re-investigate)

* **The ddt-at-restart anomaly is NOT a bug.** `ddtSchemes backward` is
  three-level, so at u = u^n, ddt = 0.5(n^oo - n^o)/dt, i.e. minus half the
  PREVIOUS step's backward difference -- not zero. Verified numerically to 7
  significant figures against the restart files. Old-time handling is correct
  (`storeOldTimes()` runs before `primitiveFieldRef()` writes).
* **chemSrcPrev_/chemOuterCount_/chemPicardChange_** are mutated per residual
  call but NOTHING in F reads them; `finalOuterIteration()`, the only
  behavioural consumer, is dead code. Consistent with the measured freeze test
  (no difference to 8 digits).
* Scaling verified consistent at all four sites; signs verified block by block
  against the residual; the chargeDensity and meanE fixes are complete.

### Latent, worth fixing (rule 31)

* `plasmaReactionRates::reportRange()` has `if (f.empty()) return;` -- a LOCAL
  guard -- with gMin/gMax and three `reduce()` downstream. Rule 31's dangerous
  spelling. Reached from the residual. Bites only on a zero-cell rank.
* `plasmaTransport.C:4073` `gMax(chemL_[sp])` behind `if (!chemL_[sp].empty())`
  -- same pattern, plus one global reduction per species per matvec.
* ~10 global reductions per matvec overall. Benign in serial, a real cost in
  parallel.
* Peak-hold diagnostics (clampRaw_, advisoryLpeak_, chemStiffnessPeak_,
  chemPicardPeak_) are mutated once per matvec and are therefore MEANINGLESS in
  Newton mode. The earlier Phase D exclusion "the clamp is not binding" came
  from one of these and is not sound evidence.

### Next step, now well specified (from a verified PETSc survey)

Assemble a Pmat and switch to PCFIELDSPLIT Schur. Measured on a coupled toy
(MFFD Amat), KSP iterations over 3 Newton steps: PCNONE 77, additive 64,
multiplicative (~= our current PCSHELL) 36, **Schur with fact_type full 7** --
equal to an exact LU of Pmat.

* PCFIELDSPLIT CANNOT work off a matrix-free Pmat, and it FAILS SILENTLY:
  `-ksp_view` shows every split as `PC type: none` because PETSc cannot extract
  sub-blocks from a MATSHELL. Verified by running. `-fieldsplit_*_pc_type
  hypre` appears to work but densely probes the block -- ~10,000 residual
  evaluations per setup, i.e. the 1500x disaster reintroduced.
* Supported route: keep `MatCreateSNESMF` as Amat, pass a real MATAIJ Pmat as
  SNESSetJacobian's second argument (PETSc manual: blocks come from Pmat, not
  Amat). Never set `-pc_use_amat` or the fieldsplit `*_use_amat` variants.
* Our layout is field-major, so use `PCFieldSplitSetIS` with stride index sets,
  NOT `PCFieldSplitSetFields`. Schur requires EXACTLY two splits, so group
  [phi] against [n_e, n_Arp, n_Ar2p, nEps_e].
* The missing block is cheap and diagonal:
  **A_01 = d(Poisson residual)/dn_i = -Z_i e**.
* SNESVINEWTONRSLS works matrix-free (verified) and has first-class fieldsplit
  support. SNESVINEWTONSSLS and SNESNEWTONTR both HARD-FAIL matrix-free
  (they need MatMultTranspose, which MATMFFD lacks). NGMRES/QN/NCG/Anderson
  ignore the Jacobian and the physics PC entirely -- skip.
* Precedent to copy: `src/snes/tutorials/ex28.c` (elliptic coefficient coupled
  to a second field, run with `-snes_mf_operator -pc_type fieldsplit`) and
  `ex19.c` for AMG on the elliptic block.

## Pmat + PCFIELDSPLIT: the design, and what petsc4Foam does and does not give us -- 2026-09-09

Rule 40 (work inside OpenFOAM's structures; confine third-party libraries
behind a thin MECHANICAL boundary) decides the shape of this, and the first
question it forces is whether the framework already does the job.

### What petsc4Foam provides -- checked, not assumed

Upstream (https://gitlab.com/petsc/petsc4foam) README lists:

    "ldu2csr: convert OpenFOAM lduMatrix to PETSc mataij (csr)"
    "Selection of solvers and preconditioners available in PETSc and in its
     external libraries", GPU support, dictionary/rc configuration, caching.

It is a LINEAR-SOLVER REPLACEMENT AT THE SEGREGATED FIELD LEVEL -- selected by
naming `petsc` as the solver for one equation in `fvSolution`. It provides NO
SNES/nonlinear layer, NO block or coupled multi-field assembly, and NO
PCFIELDSPLIT. So the conversion primitive is reusable; the nonlinear layer and
the block assembly are genuinely ours to build, and building them is not
duplicating the framework (rule 30 satisfied).

NOTE, for future upstream merges: our vendored copy at `ThirdParty/petsc4Foam`
is LOCALLY PATCHED, not pristine -- see commit 215d032 ("Corrected petsc4foam
to run for implicit and no coupling in the same case"), and it carries
`useCoupledAssembly_` / `lduPrimitiveMeshAssembly` handling that the upstream
README does not mention.

### What we reuse from it: the CONVENTIONS, not the function

`petscSolver::buildMat` (petscSolver.C:525) is the ldu2csr implementation, and
it is exported. It cannot be called directly for our purpose because it sizes
the matrix to ONE field -- `MatSetSizes(Amat, nrows_, nrows_, ...)` with
`nrows_ = lduAddr.size()`. There is no field-block offset.

But its COO layout IS the pattern to follow, and it ALREADY uses a row/column
offset `off` for the multi-region coupled-assembly case, which is structurally
the same thing we need with fields in place of regions:

    diagonal   (celli + off,   celli + off)
    upper      (low[f] + off,  upp[f] + off)
    lower      (upp[f] + off,  low[f] + off)
    + one entry per processor-interface face

CRITICAL CONVENTION, and easy to get wrong: it copies `matrix_.diag()`
DIRECTLY. That is correct only because `fvMatrix::solveSegregated` has already
folded the boundary contribution into the diagonal with `addBoundaryDiag()`
before handing the matrix to the lduMatrix solver. Converting an `fvScalarMatrix`
ourselves, we must call `addBoundaryDiag(d, 0)` on a COPY of the diagonal, or
every boundary cell's diagonal is wrong.

### The design

1. **Physics stays in OpenFOAM.** Each diagonal block is assembled with the
   SAME `fvm::` expressions the preconditioner and the real equations already
   use -- `fvm::ddt + fvm::div - fvm::laplacian + fvm::Sp` -- so the physics
   exists once, in one language.
2. **Extraction is mechanical**: read `diag()` (+ `addBoundaryDiag`), `upper()`,
   `lower()` (falling back to `upper()` when `!hasLower()`), and the lduAddressing,
   into plain `rows/cols/vals` arrays with the field-block offset applied.
3. **The boundary stays plain-typed**: those arrays cross into `snesBridge`,
   which builds the MATAIJ. No translation unit sees both `petscsnes.h` and
   `fvCFD.H` -- the existing rule that made this library buildable at all.
4. **THE COUPLING BLOCK**, which is the entire point and is cheap:
   `A_01 = d(Poisson residual)/d n_i`. DERIVE THE SIGN FROM OUR OWN RESIDUAL,
   do not copy it from elsewhere: residualCallback computes
   `F_0 = lap(effEps, phi) - rhsSource` with
   `rhsSource = -chargeDensity - dt*diffusiveChargeSource`, i.e.
   `F_0 = lap + chargeDensity + dt*(...)`, and `chargeDensity = sum_i Z_i e n_i`.
   So **`dF_0/dn_i = +Z_i e`**, a DIAGONAL block. (An external survey quoted
   `-Z_i e`; that is a different sign convention and would be silently wrong here.)
5. **SCALING IS NOT OPTIONAL.** Everything PETSc sees is scaled
   (`x_p = x/sX`, `F_p = F/sF`), so the assembled Jacobian must be too:
   every entry in block (row r, col c) is multiplied by **`sX[c]/sF[r]`**.
   Getting this wrong yields a Pmat that is a valid matrix of the WRONG system,
   which fails in the least obvious way possible.
6. **Then** `PCFIELDSPLIT`: two splits via `PCFieldSplitSetIS` (our layout is
   field-major, so NOT `PCFieldSplitSetFields`), `[phi]` against
   `[n_e, n_Arp, n_Ar2p, nEps_e]`, `-pc_fieldsplit_type schur
   -pc_fieldsplit_schur_fact_type full`. Retire the PCSHELL -- which also
   removes its incompatibility with the bounded solver's reduced vectors.

### Verify the Pmat BEFORE wiring any preconditioner onto it

The check that makes this safe, and it is cheap: compare `Pmat*v` against the
matrix-free `Amat*v` for a few random v. They should agree to the differencing
error. If they do not, the assembly, the boundary handling, the coupling sign
or the scaling is wrong, and finding that out through preconditioner behaviour
instead would be miserable. Do this first.

### Traps already paid for, do not rediscover

* A matrix-free Pmat makes PCFIELDSPLIT degrade SILENTLY to `PC type: none`.
* `-fieldsplit_*_pc_type hypre` on a shell densely probes it -- ~10,000 residual
  evaluations per setup, i.e. the 1500x defect reintroduced.
* Never set `-pc_use_amat` or the fieldsplit `*_use_amat` variants.
* Schur requires EXACTLY two splits.
* `PETSC_INFINITY` overflows under `FOAM_SIGFPE`; use a finite 1e30.
