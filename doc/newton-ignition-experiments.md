# Newton at ignition — experiment log

**Live log. Every idea tried is recorded here, including the ones that failed
and the ones that turned out to be measuring nothing.** Started 2026-09-10.

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

### 6. Through-ignition run — IN PROGRESS

`validation/newton_retry_ignition`: Newton + retry, no dt cap, from
t=1.80014e-6 through the breakdown at ~1.97e-6 to 2.05e-6. This is the
experiment the whole log is aimed at: **does Newton hold a usable dt where
Picard collapses to 1e-15?**

## Still untested (next, in priority order)

1. **Per-block scaling.** `sX`/`sF` are ONE SCALAR PER FIELD. An ignited
   discharge has `n_e` spanning 1e11–1e18 WITHIN the electron block, so a
   single scale cannot condition it, and the Newton step is computed in badly
   scaled variables. Strongest remaining candidate, and consistent with every
   PETSc knob failing.
2. **Bounded Newton** (`SNESVINEWTONRSLS`, `bounded true`). The density clamp
   is applied OUTSIDE the equations, so while it is active `F(u)=0` is not
   reachable at all. At ignition, cells sit on the floor next to cells at 1e18.
3. **Pmat quality at the ignited state.** `-snes_test_jacobian` measured 6.75e-3
   pre-ignition; re-measure it here. A Pmat that was an adequate approximation
   in the quiescent state need not be one across a streamer head.
4. **log(n) formulation.** The principled fix for both 1 and 2 at once, and
   already on the deferred list. Positivity becomes intrinsic, the clamp goes,
   and the seven-decade dynamic range becomes an O(1) range in the unknown.
