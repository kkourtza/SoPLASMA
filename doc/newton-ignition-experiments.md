# Newton at ignition — experiment log

**Live log. Every idea tried is recorded here, including the ones that failed
and the ones that turned out to be measuring nothing.** Started 2026-09-10.

## SUMMARY — read this first

| # | idea | verdict | the number that decides it |
|---|---|---|---|
| 5 | **Newton retries a failed step** | **WORKS — the one that mattered** | dt 1e-11 (hard cap) → **1.5–4e-10 self-regulated**, 15–40× |
| 3 | FGMRES instead of GMRES | works, kept | removed `DIVERGED_BREAKDOWN` at the restart boundary |
| 8 | rebalance residual scales | correct, neutral | blocks 0.0077–111 → all 44.72; KSP unchanged |
| 2 | preconditioner (hypre/bjacobi/selfp/500 its) | **no effect** | all `DIVERGED_ITS`, PC verified by `-ksp_view` |
| 4 | 8 PETSc options (EW, 3 line searches, mffd_err) | **all failed** | none survives the dt baseline dies at |
| 9 | **bounded Newton** (`SNESVINEWTONRSLS`) | **FAILS BADLY** | dt collapses to **5e-15**, 10⁵× worse than unbounded |
| 1 | dt cap sweep | diagnostic | ceiling between 1e-11 (works) and 1e-10 (dies) |

**The one-line conclusion so far:** the binding constraint was never the
preconditioner or any PETSc knob — it was that **Newton had no way to back off
from a step that was too large.** Giving it the retry path Picard already had
raised its sustained dt by 15–40× and made the Picard→Newton handover survivable
at all.

**Two traps this log exists to stop anyone repeating:** read the ACHIEVED dt,
never the step count (bounded Newton's 360 steps advanced 4e-11); and confirm
the knob moved before believing a comparison (four separate inert-knob
incidents, now rule 42).

## The case, fixed for every experiment below

`grubert2009_ballast400_picard` — Grubert DC glow, ballasted electrode
(`seriesResistor`, R = 1e8), driven by a SHARP ramp to −400 V over 4 µs
(100 V/µs, ~10× the standard bed's 10 V/µs, so breakdown at −121 V arrives at
t ≈ 1.21 µs instead of 12.1 µs). `adjustTimeStep true`.

It **ignites at t ≈ 1.97e-6**: `n_e` goes from the 1e11 seed to 4.7e18, `Arp`
to 5.7e18, `Emag` to 3.7e6 V/m.

**Picard's failure is a dt COLLAPSE, not a crash** — measured on the continuous
run:

    t = 1.97144e-06   dt = 4.6e-13
    t = 1.97156e-06   dt = 8.7e-13
    t = 1.97167e-06   dt = 7.2e-13
    t = 1.97177e-06   dt = 7.2e-15   <-- collapse
    t = 1.97177e-06   dt ~ 1e-15, stuck, 31-42 correctors/step

That collapse is the thing Newton has to beat. Everything below is aimed at it.

**Restart points.** Snapshots at 2e-7 spacing exist from 2.0e-7 to 1.8e-6;
`ballast400_fine` is re-running Picard from 1.80014e-6 with 1e-8 spacing to get
a restart point immediately before the collapse, because 1.8e-6 → ignition at
dt 1e-11 is ~17,000 steps.

## Established, with the verification that makes it trustworthy

Rule 42 applies throughout: a sweep is only reported once the swept parameter
is shown to have MOVED.

### 1. Newton has a hard dt ceiling near 1e-11 on this state — CONFIRMED

Restart from t=1.80014e-6, Newton on, only `maxDeltaT` varied. dt achieved was
confirmed different per arm, and SNES confirmed engaged:

| maxDeltaT | dt actually used | steps | fatal | SNES mentions |
|---|---|---|---|---|
| 1e-10 | 1e-10 | 1 | **1** | 6 |
| 1e-11 | 1e-11 | 154 | 0 | 769 |
| 1e-12 | 1e-12 | 256 | 0 | 1274 |

**`maxDeltaT` must go in the `plasmaTimeControl` dict** inside
`system/plasmaSimulationControls`. In `system/controlDict` it is read by
OpenFOAM's `Time` and then silently overridden by `plasmaTimeControl`'s own
lookup (default `GREAT`).

### 2. The preconditioner is NOT the binding constraint — CONFIRMED NEGATIVE

All at the inherited dt = 1.29e-9, i.e. >100× above the ceiling in (1), which
is why none of it could have helped. PC verified with `-ksp_view`,
`-options_left` showed nothing unused.

| arm | PC actually built | result |
|---|---|---|
| baseline | ilu / ilu | DIVERGED_ITS 100 |
| `-ksp_max_it 500` | ilu / ilu | DIVERGED_ITS 500 |
| `-fieldsplit_phi_pc_type hypre` | **hypre** / ilu | DIVERGED_ITS 100 |
| hypre + bjacobi + 500 its | **hypre** / **bjacobi** | DIVERGED_ITS 500 |
| `-pc_fieldsplit_schur_precondition selfp` | ilu / ilu | DIVERGED_ITS 500 |

AMG on the elliptic block plus 5× the Krylov budget changed nothing.

### 3. FGMRES fixed a genuine breakdown — CONFIRMED, and kept

With adaptive dt the Krylov solve failed `DIVERGED_BREAKDOWN` at **iteration
30**, GMRES's default restart length. `PCFIELDSPLIT`+Schur solves its blocks
iteratively, so the preconditioner VARIES between Krylov iterations and plain
GMRES loses the Arnoldi relation (Knoll & Keyes §3.5). `-ksp_type fgmres
-ksp_gmres_restart 100` removed the breakdown; the failure mode became the
honest `DIVERGED_ITS`. Now a baked-in default.

## Measurements that turned out to be MEASURING NOTHING

Recorded because the failure mode recurred four times in one day and is now
rule 42.

| what was "tested" | why it was inert |
|---|---|
| `-fieldsplit_0_pc_type hypre` | `PCFieldSplitSetIS` names the splits `phi` and `transport`; there is no split `0`. Two "PC tests" compared a config against itself. |
| `maxDeltaT` in `system/controlDict` | `plasmaTimeControl` reads its own dict; three arms ran at the same inherited dt = 1.29e-9 and returned byte-identical output. |
| `sed s/outerSolver.*/newton/` | The ballast case has NO `outerSolver` key (it descends from `ballast_clean`), so the sed was a no-op and two "Newton" arms ran **Picard** for 24,000 steps. |
| Phase D "40,000-step Newton success" | `anySpeciesOnFloor()` was re-tested every step; Newton was off for 21,174 of them. |

## Two solver bugs found and fixed along the way

**The handover guard asked `gMin`.** `anySpeciesOnFloor()` tested
`gMin(n) <= floor` — "is ANY cell at the floor?" — which is true FOREVER in a
real discharge because the quiescent far field always rests on the floor. On
this very case `Arp` has `minNumberDensity 1e11` and its min sits at exactly
1e11 for the whole run while its max reaches 5.7e18, so Newton could never be
handed to on the one case it was wanted for. Now `gMax(n) <= floor` — "is this
species at its floor EVERYWHERE", the actual cold-start degeneracy.

**The handover did not latch.** It flipped back to Picard mid-run whenever a
species touched its floor, silently (one-shot message). Now latched.

## Ideas queued, with what each is expected to address

The ceiling in (1) is the target. Each is tested at **dt = 1e-10, where
baseline Newton fails** — a binary, cheap, attributable outcome.

- **Eisenstat–Walker adaptive forcing** (`-snes_ksp_ew`). Standard JFNK
  practice: early Newton iterations do not need a tight linear solve. Directly
  targets `DIVERGED_ITS`, because it LOOSENS the tolerance the KSP must reach.
- **Line search** (`-snes_linesearch_type` `bt`/`l2`/`cp`/`basic`). A poor line
  search can stall an otherwise healthy Newton.
- **Differencing step** (`-mat_mffd_err`). The ignited state has densities
  spanning 1e11–1e18; the matrix-free directional derivative may be badly
  scaled.
- **Bounded Newton** (`SNESVINEWTONRSLS`). The density clamp is applied OUTSIDE
  the equations, so `F(u)=0` is not reachable while it is active; bounds put it
  inside.
- **Per-block scaling adequacy.** A single scalar scale per field cannot
  condition a block whose own values span seven decades, which is exactly what
  an ignited discharge has.
- **Newton retry on non-convergence.** Newton `FatalError`s where Picard has
  `onNonConvergence retryStep`. With a retry path Newton would find its own
  ceiling instead of dying at it — likely a precondition for the automatic
  Picard↔Newton switch.

## Results

### 4. EIGHT PETSc-option ideas, ALL FAILED — the ceiling is STRUCTURAL

Every arm: restart from the saved t=1.80014e-6 state, `outerSolver newton`,
`maxDeltaT 1e-10` (the level at which baseline Newton dies after one step),
exactly one option changed. dt achieved and SNES engagement confirmed per arm.

| idea | option | verdict |
|---|---|---|
| baseline (control) | — | failed |
| Eisenstat–Walker | `-snes_ksp_ew` | failed |
| line search: basic | `-snes_linesearch_type basic` | failed |
| line search: L2 | `-snes_linesearch_type l2` | failed |
| line search: critical point | `-snes_linesearch_type cp` | failed |
| differencing step, coarse | `-mat_mffd_err 1e-3` | failed |
| differencing step, fine | `-mat_mffd_err 1e-7` | failed |
| EW + L2 combined | both | failed |

Taken with (2), **no PETSc knob rescues dt = 1e-10.** Not the preconditioner,
not the Krylov budget, not the forcing term, not the line search, not the
differencing step. The ceiling is a property of the FORMULATION, not of the
solver configuration — which is what makes the remaining ideas (scaling, the
out-of-equation clamp, Pmat quality) the ones worth pursuing.

Eisenstat–Walker failing is the most informative of these: it LOOSENS the
linear tolerance, so if the problem were merely "the KSP cannot reach a tight
rtol in the budget", EW would have papered over it. It did not, so the linear
system is not merely expensive — the Newton step it defines is wrong or the
operator is genuinely intractable at that dt.

### 5. RETRY ON NON-CONVERGENCE — THE ONE THAT WORKED

Implemented (see below) and tested with **no dt cap at all**, so Newton
inherits the dt = 1.29e-9 that was previously an immediate FatalError:

    retry 1 of 5, deltaT 1.29279923982e-09 -> 6.46399619909e-10
    retry 2 of 5, deltaT 6.46399619909e-10 -> 3.23199809954e-10
    ...
    62 steps, 0 fatals, 15 retries, 556 SNES calls
    sustained dt ~ 1.2e-10 .. 2.3e-10

**Self-regulation beats any fixed cap.** The hard-capped ceiling measured in
(1) was 1e-11; letting Newton attempt large steps and pay for the occasional
rejection lands it at **1.5e-10 – 2.3e-10, i.e. 15–23× larger**, with no
failures. A fixed cap has to be conservative enough for the worst step in the
run; a retry ladder only has to be right on average.

This also removes the structural reason the Picard→Newton handover was fatal:
Newton inherits whatever dt the easy Picard phase wound up to, and now backs
off from it instead of dying on contact.

**Implementation, three pieces:**
* `plasmaNewtonSolver::lastSolveConverged()` — new virtual, default `true`.
  `snesNewtonSolver` reports the SNES reason through it and WARNS instead of
  raising `FatalError`.
* `plasmaTimeControl::noteOuterFailure()` — the Newton equivalent of hitting
  the corrector cap, so `onNonConvergence` means ONE thing whichever outer
  solver runs, and the policy lives in one place. Deliberately not expressed as
  `noteOuterLoop(maxCorrectors)`, which would fake a corrector count and
  corrupt the reported `outerItersUsed_`.
* A separate `outerSolveFailed_` flag rather than reusing `outerHitCap_`:
  `noteOuterLoop()` assigns that one unconditionally and runs AFTER the Newton
  call, so a shared flag would be erased by the very next statement. Cleared at
  each step start and on each retry.

### 11. THROUGH-IGNITION A/B, RUNNING OVERNIGHT (2026-09-10 ~02:05)

Two processes, identical case, identical restart point, **only the Schur
coupling block of (10) differs**:

| case | library | started from |
|---|---|---|
| `validation/newton_retry_ignition` | OLD (no coupling block) — **control** | t=1.80014e-6, still running |
| `validation/newton_ignition_fixed` | FIXED (coupling block) | t=1.86023e-6 |

The control was launched at 01:44:16, before the 01:59:11 rebuild, so it kept
the old `.so` mapped — which makes it a genuine control rather than a wasted
run. **Verified at the INODE level** rather than assumed, because a rebuilt
shared library is exactly the kind of thing that silently does not take effect:

    on-disk (fixed) inode  622918
    pid 4332  newton_retry_ignition  inode 616626   <= OLD lib
    pid 13964 newton_ignition_fixed  inode 622918   <= FIXED lib

The fixed run restarts from the control's own latest snapshot, so the 280 steps
of progress are kept — rule 41 paying for itself the same day it was written.

**What to compare in the morning.** The control's log already covers
t=1.86023e-6 onward, so extract its behaviour from that instant and compare
like-for-like against the fixed run:
* did either reach ignition (~1.97e-6), and did dt survive it?
* dt sustained through the breakdown — **the number that decides everything**.
  Picard collapses to ~1e-15 there;
* retries per step, and DIVERGED_ITS count per solve;
* wall clock per ns of simulated time.

### 6. Earlier through-ignition run — superseded by (11)

`validation/newton_retry_ignition`: Newton + retry, no dt cap, from
t=1.80014e-6 through the breakdown at ~1.97e-6 to 2.05e-6. This is the
experiment the whole log is aimed at: **does Newton hold a usable dt where
Picard collapses to 1e-15?**

### 7. THE BLOCKS ARE NOT COMMENSURATE — a 14,000× spread, MEASURED

The residual callback's own diagnostic, at the first Newton call of a failing
step (`idea_baseline`, dt = 1e-10):

    Poisson   |F_scaled| =  44.7
    n_e       |F_scaled| =  35.6
    n_Ar2p    |F_scaled| =   0.147
    n_Arp     |F_scaled| =   0.0077
    nEps_e    |F_scaled| = 111.2

Every one of these is supposed to be O(1) — that is the entire purpose of the
`typ u` scaling (Knoll & Keyes 2.3.1: no component should dominate the norm
merely through its units). They span **0.0077 to 111**.

**Why that is fatal rather than untidy:** a Krylov method minimises the norm of
the WHOLE vector. In this state the norm is essentially `nEps_e` plus Poisson;
`n_Arp` contributes ~1e-4 of it and therefore almost nothing to the Krylov
space, no matter how wrong the ion block is. GMRES is fitting two blocks and is
nearly blind to a third.

**The cause is in the code, not the physics.** `sF = sX/dt` for the transported
blocks is an A PRIORI estimate of how large the residual OUGHT to be. The
measurement above is what it actually is.

**Fix implemented: normalise `sF` by the MEASURED initial residual** (new
`rebalanceScales`, default true). It costs nothing — the priming residual
evaluation already had to happen so that Pmat and F describe the same state —
and it is done BEFORE the Pmat assembly, which reads `sX`/`sF`, so matrix and
residual stay consistent. A/B results below.

Note this also explains why every PETSc knob failed: no choice of
preconditioner, forcing term or line search repairs a residual in which one
block is invisible. That is a property of the vector being handed to the Krylov
method, not of how the Krylov method is run.

### 8. RESIDUAL REBALANCING — works exactly as designed, does NOT fix convergence

Verified it moved the knob (the first A/B printed identical numbers for both
arms, because the diagnostic fires on the PRIMING call which runs BEFORE the
rebalancing — a rule-42 near-miss; the later calls are the ones that show it):

| block | `rebalanceScales true` | `false` |
|---|---|---|
| Poisson | 44.7214 | 44.7190 |
| n_e | **44.7214** | 35.6396 |
| n_Ar2p | **44.7215** | 0.1465 |
| n_Arp | **44.7214** | 0.00766 |
| nEps_e | **44.7213** | 111.168 |

Every block now contributes equally. **And convergence did not improve**: KSP
counts 24/49/100 against 25/49/100, one retry each. So block imbalance was A
defect, not THE defect. Kept because it is correct, cheap, and removes a
confound from every later measurement.

### 9. BOUNDED NEWTON (`SNESVINEWTONRSLS`) — CLEAR FAILURE, dt collapses to 5e-15

The hypothesis was that the density clamp, applied OUTSIDE the equations, makes
`F(u)=0` unreachable, and that imposing the floor as a CONSTRAINT inside the
solve would fix it. It does the opposite.

    360 steps, 0 fatals, 2383 SNES calls   <-- looks like a success
    final t = 1.80018e-06                  <-- advanced 4e-11 in 360 steps
    dt      = 5.0e-15                      <-- COLLAPSED
    135 retries out of 360 steps

**The step count was the trap.** 360 steps of nothing. Bounded Newton grinds at
dt ~5e-15, which is *the same collapse Picard suffers at ignition* and 10⁵×
worse than unbounded Newton on the identical state (2–4e-10). Read the achieved
dt, never the step count.

Why it likely fails: `SNESVINEWTONRSLS` is a reduced-space active-set method.
At this state a large fraction of cells sit exactly on the floor, so the active
set is huge and changes every iteration, and the reduced space it solves in is
both small and unstable between iterations.

**Verdict: do not pursue bounded Newton for this problem.** The clamp is still a
real defect, but the fix is to remove the need for it (log(n)), not to
re-express it as a constraint.

### 10. ROOT CAUSE FOUND — the Schur complement was INERT

**The Pmat was missing d(species)/d(ePotential), so PCFIELDSPLIT's Schur
complement contributed nothing at all.**

Every block the Pmat actually contained:

| block | content | present? |
|---|---|---|
| (φ, φ) | `laplacian(ε, dφ)` | yes |
| (φ, nᵢ) | `q` — charge density, species → Poisson | yes |
| (nₛ, nₛ) | species diagonal | yes |
| (ε, ε) | energy diagonal | yes |
| **(nₛ, φ)** | **d(drift)/d(φ) — Poisson → species** | **NO** |
| **(ε, φ)** | **d(energy drift)/d(φ)** | **NO** |

PCFIELDSPLIT with a Schur complement on splits {φ, transport} forms

    S = A_tt − A_tφ · A_φφ⁻¹ · A_φt

`A_φt` (charge density) was assembled. `A_tφ` was absent, i.e. **zero** — so
**S = A_tt exactly**. The assembled Pmat, the two index sets and the Schur
factorisation were all doing the work of a block-TRIANGULAR preconditioner that
knows the densities move the field but not that the field moves the densities.
The entire reason the Pmat exists was defeated.

**This explains every observation in this log:**
* harmless while the coupling is weak (pre-ignition, a linear capacitive ramp),
  fatal once space charge makes the field and the field makes the ionisation;
* no sub-preconditioner could ever help — hypre, bjacobi, `selfp`, 5× the
  Krylov budget, Eisenstat–Walker, three line searches and both differencing
  steps all failed, because **the missing physics is not in the matrix for any
  of them to precondition**;
* residual rebalancing was correct and yet neutral, for the same reason.

**The derivative.** `phiE = -snGrad(ePotential)*magSf`
(`singleRegionPoisson.C:39`) and the drift term is `div(Z·μ_f·phiE·n)`, so

    d/dφ [ div(Z·μ·phiE·n) ] = -laplacian(Z·μ·n, dφ)

since `fvm::laplacian(Γ,ψ)` IS `div(Γ_f·snGrad(ψ)·magSf)`. The energy row gets
the same term with `Ze·μ_eps·nEps`.

**The (φ, energy) block is NOT the mirror of this and is correctly absent:** the
Poisson residual depends on species only through `chargeDensity = Σ qᵢnᵢ`, and
`nEps` carries no charge, so `d(F_φ)/d(nEps)` is identically zero. That coupling
is genuinely one-way.

**Implemented** as `blockMatrixCOO::addFvMatrixBlock(rowField, colField, ...)`,
generalising `addFvMatrix` (which is now the `rowField == colField` case) and
carrying the processor-interface handling across unchanged.

**Still missing from the energy row, deliberately and recorded rather than
skipped:** `d(Psrc)/d(φ)`, the response of JOULE HEATING to the field. `Psrc ~
J·E` is a direct and probably stronger dependence on the potential than the
drift term retained here, but it is not laplacian-shaped and needs its own
derivation.

## Still untested (next, in priority order)

1. **Per-CELL scaling.** `sX` is one scalar per FIELD. An ignited discharge has
   `Arp` spanning 1e11–5.7e18 WITHIN one block (7.7 decades), so cells at the
   floor scale to ~1e-8 of the block's own scale and their Jacobian columns are
   differenced at the noise level. Item (7) fixes the imbalance BETWEEN blocks;
   this is the same disease WITHIN one.
2. **Bounded Newton** (`SNESVINEWTONRSLS`, `bounded true`). The density clamp
   is applied OUTSIDE the equations, so while it is active `F(u)=0` is not
   reachable at all. At ignition, cells sit on the floor next to cells at 1e18.
   (Test running.)
3. **Pmat quality at the ignited state.** `-snes_test_jacobian` measured 6.75e-3
   pre-ignition; re-measure it here. A Pmat that was an adequate approximation
   in the quiescent state need not be one across a streamer head. (Test
   running.)
4. **log(n) formulation.** The principled fix for 1 and 2 at once, and already
   on the deferred list. Positivity becomes intrinsic, the clamp goes, and the
   seven-decade dynamic range becomes an O(1) range in the unknown. Everything
   measured today argues for it: the failures are all dynamic-range failures.
