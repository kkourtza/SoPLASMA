# Newton at ignition — experiment log

**Live log. Every idea tried is recorded here, including the ones that failed
and the ones that turned out to be measuring nothing.** Started 2026-09-10.

## SUMMARY — read this first

| # | idea | verdict | the number that decides it |
|---|---|---|---|
| 5 | **Newton retries a failed step** | **WORKS** | dt 1e-11 (hard cap) → **1.5–4e-10 self-regulated**, 15–40× |
| 15 | **Eisenstat–Walker (inexact Newton)** | **BIGGEST OPTION WIN — now default** | mean dt **+80%**, KSP failures **13× fewer** |
| 12 | **ionisation derivative `dP/dn` in Pmat** | **WORKS — now default** | **+41%** advance/step; it IS the avalanche growth rate |
| 10 | **Schur coupling block `dF_s/dφ`** | **real bug — Schur was INERT** | `A_tφ = 0` ⇒ `S = A_tt` exactly; +7% |
| 17 | *(consequence)* | **regime change** | solver-driven rejections **38% → 8.6%** of steps |
| 3 | FGMRES instead of GMRES | works, kept | removed `DIVERGED_BREAKDOWN` at the restart boundary |
| 8 | rebalance residual scales | correct, neutral | blocks 0.0077–111 → all 44.72; KSP unchanged |
| 16 | `chemCrossJacobian` (my own idea) | **HURTS** | **−29%**; approximation too crude |
| 2 | preconditioner (hypre/bjacobi/selfp, retested) | **no effect** | −9%, −1%; verified by `-ksp_view` |
| 9 | **bounded Newton** (`SNESVINEWTONRSLS`) | **FAILS BADLY** | dt collapses to **5e-15**, 10⁵× worse |
| 14 | `ngmres`, Schur `upper`, quadratic backtracking | **all fail** | singular wall flux; 45%; 20% KSP failure |

**The one-line conclusion:** the preconditioner options were never the problem
and could never have been, *while the matrix was missing physics*. Two Jacobian
blocks were absent — the field→species drift coupling (which made the Schur
complement mathematically inert) and the ionisation derivative (the avalanche
growth rate itself). Once the Pmat was right, the option that ADAPTS to Jacobian
quality — Eisenstat–Walker, which had *failed* on the broken matrix — nearly
doubled the timestep and cut linear-solve failures 13×. Newton is now limited by
the temporal-accuracy controller rather than by its own convergence.

**Three traps this log exists to stop anyone repeating:** read the ACHIEVED dt,
never the step count (bounded Newton's 360 steps advanced 4e-11); confirm the
knob moved before believing a comparison (four inert-knob incidents, now rule
42); and **re-test options after fixing the model they were tested against** —
EW, `hypre`, `maxit`, and the line searches were all dismissed on evidence
gathered with a broken Pmat, and one of them turned out to be the biggest win
available.

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

**AND PICARD DOES NOT MERELY GRIND — IT DIES.** The reference run
`ballast400_fine` (pure Picard, from t=1.80014e-6) aborted at ignition with the
solver's own give-up message:

     The solver cannot converge this case.
    10 consecutive time steps were accepted without the outer (PIMPLE) loop
    reaching its residual tolerance, with deltaT already at 1.7986e-14 s
    (floor 1.7986e-14 s, derived).

        time                     1.971771e-06
        deltaT                   1.7986e-14 s
        largest deltaT sustained 1.29279923982e-09 s
        correctors on last step  150 (cap 150)
        n_e max                  4.61e+18

**dt collapsed 72,000× — 1.29e-9 down to the 1.8e-14 floor — and the run
stopped.** "Reducing the time step is the solver's last automatic lever and it
is used up."

THIS IS THE BENCHMARK. Newton has to get past t = 1.971771e-06 at a usable dt.
Everything below is aimed at it.

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

### 12. THE IONISATION DERIVATIVE — the biggest Pmat win so far

Second missing Jacobian term, found by asking the same question as (10) about a
different operator. `F = ddt + div - lap - P + L*n`, so `dF/dn = L - dP/dn`. The
Pmat carried `fvm::Sp(chemL)` — the `L` — and **not `-dP/dn`**. For
electron-impact ionisation `P_e = k_iz·n_e·n_Ar`, so

    dP_e/dn_e = k_iz·n_Ar = THE IONISATION FREQUENCY

the avalanche growth rate, i.e. the defining physics of ignition. Without it the
preconditioner believed the species equation was far more diagonally dominant —
far more stable — than it is, exactly where it is not.

Approximated `dP_s/dn_s ≈ P_s/n_s`, exact whenever production is first order in
the species itself. The term is NEGATIVE on the diagonal deliberately: an
avalanche genuinely is unstable, and a preconditioner that hides that is
describing a different problem.

| `chemJacobian` | advance/step | mean dt (90) |
|---|---|---|
| off | 9.50e-11 | 1.580e-10 |
| **on** | **1.344e-10 (+41%)** | **1.855e-10 (+17%)** |

For scale, the drift-coupling block of (10) was worth ~7%.

### 13. WHY IT STILL FAILS — the linear solve, every time

Diagnosed rather than guessed, on the best current configuration:

    retry causes      10 of 10 : SNES reason -3 (DIVERGED_LINEAR_SOLVE)
    SNES reasons      3 (37x), 4 (14x), -3 (10x)
    KSP on failure    100 iterations = the cap, every time

**Every single failure is the Krylov solve running out of iterations.** Not the
line search, not the nonlinear iteration, not a bad state. That makes the linear
solve the sole remaining target, and it retroactively justifies retesting the
options that were dismissed in (2) and (4) — those all ran against the BROKEN
Pmat, so they were never a fair test of "help the KSP".

Consistent with that reading, the two options that help most are exactly the two
that make the KSP's job easier rather than better-preconditioned:

| on top of `chemJacobian on` | mean dt (90) |
|---|---|
| — | 1.855e-10 |
| `-mat_mffd_type ds` | 1.936e-10 |
| `-ksp_rtol 1e-3` | **1.960e-10** |

### 14. Failures from the overnight batch, recorded

* `-snes_type ngmres` — **fails hard**, drives the state to a singular
  wall-flux condition ("the drift into the wall now exceeds the thermal speed")
  within 2 steps.
* `-pc_fieldsplit_schur_fact_type upper` — 45% KSP failure rate against 2.9%
  for `full`. Stopped early.
* `-snes_linesearch_type bt -snes_linesearch_order 2` — 20% KSP failure rate.
  Stopped early.

### 15. EISENSTAT–WALKER — +80% on dt, and a 13× drop in linear-solve failures

The single largest gain of the night, and it is the option that **failed** in
the first sweep (4). It failed there because it ran against the BROKEN Pmat: an
adaptive forcing term chooses the linear tolerance from how well the previous
Newton step actually reduced `||F||`, so it needs a Jacobian whose steps mean
something. With the Schur coupling (10) and the ionisation derivative (12) in
place, it finally has one.

All arms: same restart state, `chemJacobian true` unless stated, first 90 steps.

| arm | mean dt | vs baseline | KSP failure |
|---|---|---|---|
| `chemJacobian off` (baseline) | 1.773e-10 | — | 2.2% |
| `chemJacobian on` | 1.917e-10 | +8% | 2.7% |
| `-ksp_rtol 1e-3` | 2.015e-10 | +14% | 2.1% |
| `-mat_mffd_type ds -ksp_rtol 1e-3` | 2.310e-10 | +30% | 1.6% |
| **`-snes_ksp_ew`** | **3.190e-10** | **+80%** | **0.2%** |
| `-fieldsplit_transport_ksp_type gmres` (5 its) | 1.856e-10 | +5% | 4.0% |
| `-fieldsplit_phi_pc_type hypre` (retest) | 1.618e-10 | **−9%** | 2.6% |
| `-pc_fieldsplit_schur_precondition selfp` (retest) | 1.750e-10 | −1% | 3.2% |

Note that hypre and `selfp`, retested against the CORRECTED Pmat, are still no
help — so their earlier dismissal stands on better evidence than it did.

**This is the pattern of the whole night:** the preconditioner options were
never the problem, and were never going to be, while the matrix was missing
physics. Once the Pmat was right, the one option that adapts to the Jacobian's
actual quality paid immediately.

### 16. `chemCrossJacobian` — MY OWN IDEA, AND IT HURTS

Recorded prominently because it was a plausible follow-up to (12) and is wrong.
Approximating ALL of a species' production as first order in `n_e` gives

    mean dt 1.358e-10 against 1.904e-10 with it off   -- MINUS 29%

The approximation is right for electron-impact channels and wrong for
ion-neutral and recombination ones, and evidently the wrong part dominates.
Default stays FALSE. Do not turn it on without deriving the real per-reaction
derivative.

### 17. THE BOTTLENECK MOVED — Newton is now ACCURACY-limited, not solver-limited

The most important result of the night, and it is a change of regime rather
than another percentage.

`plasmaTimeControl` states its own health test in the source:

> *"This is the ACCURACY layer, and it is meant to be the primary setter of
> deltaT. The Courant and stiffness limits below remain as CAPS: **if one of
> them still names itself in `deltaT set by`, the controller is not in charge
> and that is the thing to investigate.**"*

Measured by that criterion:

| config | `temporal` (accuracy) | **`rejection`** (solver failed) | `Co_conv` |
|---|---|---|---|
| `chemJacobian off` | 76 | **49 — 38%** | 3 |
| `chemJacobian on` | 82 | **38 — 30%** | 4 |
| **`+ -snes_ksp_ew`** | 59 | **7 — 8.6%** | 14 |

**Solver-driven step rejections fell from 38% of steps to 8.6%.** The thing
setting the timestep is now the temporal-error controller, which is precisely
what the design intends. Further solver work has diminishing returns on this
case; the next lever is the accuracy tolerance itself, which is a physics
trade-off and the user's call, not a numerics fix.

And the FAILURE MODE moved with it. Re-diagnosed on the EW run:

    without EW   every retry: SNES reason -3 (DIVERGED_LINEAR_SOLVE), KSP at its cap
    with EW      -6 LINE_SEARCH (5x), -5 MAX_IT (4x), -3 LINEAR_SOLVE (2x)
                 and 9 solves hitting the maxIt = 50 cap

So the remaining failures are in the NONLINEAR iteration, not the linear one.
That makes `maxIt` and the line-search type relevant again — both were tested
early and dismissed, but against the broken Pmat AND without EW, so neither was
ever a fair test of this regime. Under test now.

### 18. IMPORTANT QUALIFICATION — with EW on, `chemJacobian` adds nothing

Measured at EQUAL SAMPLE SIZE, which is what exposed it. Both arms verified to
carry the setting they claim and `-snes_ksp_ew` confirmed reaching PETSc:

| arm | mean dt @16 | mean dt @30 |
|---|---|---|
| `chemJacobian on`, **no EW** | 1.998e-10 | 1.912e-10 |
| EW + `chemJacobian on` | 3.653e-10 | 3.325e-10 |
| EW + `chemJacobian OFF` | **3.653e-10** | **3.316e-10** |

**Identical.** The +41% attributed to the ionisation derivative in (12) was
measured WITHOUT EW; once EW is on, EW subsumes it. That makes sense
mechanically — an adaptive forcing term chooses the linear tolerance from how
well the previous step reduced `||F||`, so it compensates for a worse Jacobian
by simply not over-solving against it.

**`chemJacobian` stays default TRUE anyway**, for reasons that are not
performance: it is a genuinely missing Jacobian term, it is nearly free, it is
worth +41% whenever EW is off or unavailable, and it is never harmful. But the
honest attribution is that **EW is doing the heavy lifting**, and the earlier
commit message overstates the ionisation derivative's standalone importance in
the final configuration.

This is the unequal-sample trap for the fourth time in one session: `ew_nochem`
first appeared to BEAT `ew_v2` (3.416e-10 against 2.955e-10) purely because its
average covered 34 steps and the other's covered 93. Always compare at equal N.

### 19. Options stacked on top of EW — only the differencing type adds anything

Equal-sample (N=16), all with EW:

| on top of EW | mean dt |
|---|---|
| — | 3.653e-10 |
| **`-mat_mffd_type ds`** | **3.926e-10 (+7%)** |
| `-snes_linesearch_type l2` | 3.657e-10 (0%) |
| `maxIt 50 → 200` | 3.074e-10 (−16%) |
| `-snes_linesearch_type cp` | 2.782e-10 (−24%) |

Raising `maxIt` HURTS, which is counter-intuitive and worth stating: letting a
struggling Newton solve grind for 200 iterations delays the retry that would
have found a workable step sooner. The retry ladder is a better response to a
hard step than more iterations on it.

### 20. EW IS STATE-DEPENDENT — +80% on the easy state, ~neutral on the hard one

**RESOLVED by a clean isolation, and it corrects TWO of my own readings.**
Identical library, identical restart at t=1.90026e-6, one switch:

| equal N=30 | mean dt | advanced | retries | MAX_IT | LINE_SEARCH |
|---|---|---|---|---|---|
| **EW ON** (default) | **2.536e-10** | **5.39e-09** | 9 | **5** | 4 |
| EW OFF | 2.421e-10 | 5.11e-09 | 6 | 1 | 4 |

So at the deep-avalanche state EW is **+4.7% on dt** — marginally BETTER, not
worse — while causing **5× more `DIVERGED_MAX_IT`** and 50% more retries.

**Both of my earlier readings were wrong, in opposite directions:**
* the "+80%, now default" headline (15) was measured at the PRE-IGNITION state
  and does not generalise — at the hard state the gain is ~5%;
* the "EW may be hurting" flag was based on (a) a CONFOUNDED three-way
  comparison differing in three changes at once, and (b) a 4–5 step sample of
  the clean test, which was noise. At N=30 it reverses.

**Verdict: keep `adaptiveForcing true`.** It is a large win over most of a run
(the long pre-ignition phase is where most steps are), neutral-to-slightly-
positive at ignition, and never actually harmful. But the caveat is real and
belongs next to the default: **EW trades linear-solve effort for nonlinear
iterations, and at a hard state that trade is roughly break-even** — which is
why `MAX_IT` failures multiply even as dt holds up.

That is the same shape as the `maxIt` result: an option tuned on the easy phase
is not automatically right for the hard one. **Measure both states before making
anything a default** — I did not, and got a headline number that does not
generalise.

### 21. THE VERDICT RUN — Newton did NOT pass the benchmark, and WHY matters

`newton_ignition_fixed`, restarted at t=1.86023e-6, ran **14,301 steps with
1,534 retries** and stopped at **t = 1.949675e-06** — short of Picard's death
point at 1.971771e-06.

**But it failed for a completely different reason.** Not convergence:

    wall-flux condition on patch anode for field n_e has become singular.
    The drift into the wall now exceeds the thermal speed by more than the
    near-wall diffusive velocity D/delta, so D/delta + uEff <= 0.
    Present settings on this patch: reflection r = 0.36, fluxScheme `standard`.

That is a **boundary-condition model limit**, and the solver's own error text
names the remedy: `fluxScheme ScharfetterGummel`, whose denominator
`D/delta*Bern(Pe) + uAbs` has both terms non-negative for any `r`, so it cannot
invert — and which is "the better scheme for the drift-dominated cell this
failure happens in" anyway.

**Newton is NOT immune to the stiffening, only more graceful about it:**

| | Picard | Newton |
|---|---|---|
| died at | t = 1.971771e-06 | t = 1.949675e-06 |
| cause | **dt floor, 150/150 correctors** | **wall-flux BC singular** |
| dt at death | 1.7986e-14 | **9.74e-14 (5.4× larger)** |
| dt collapse | **72,000×** | **~2,000×** |
| `n_e` max | 4.61e+18 | **1.50e+19** |

So Newton's dt collapsed too — a 2,000× collapse is not a triumph — but 36×
less severely, and it carried the discharge to 3× the electron density before
stopping on a different failure entirely.

**HONEST READING: this is not "Newton beats Picard".** Newton stopped earlier in
TIME. What it does establish is that the *solver* is no longer the thing that
stops the run — a boundary-condition scheme is — and that is a different and
more tractable blocker.

### 22. `fluxScheme ScharfetterGummel` — UNDER TEST

And a rule-42 lesson on the way in. The first attempt edited `fluxScheme` into
the **patch dictionaries** of the restart time directory. That is the wrong
place and the edit was silently stripped at startup: the BC reads the scheme
from the TRANSPORT MODEL (`ddModel.fluxScheme()`), not from the patch entry. The
run looked like it was testing SG and was testing nothing.

The real knob is per-species, and this case already parameterises it:

    constant/plasmaSpeciesProperties:  driftDiffusionCoeffs { fluxScheme $driftDiffusionFluxScheme; }
    configuration/config:              driftDiffusionFluxScheme  standard;   ->  ScharfetterGummel;

Confirmed in force by the run's own species table before trusting it. Note this
changes the transport DISCRETISATION everywhere, not just at the wall, so it is
a physics-affecting change and its accuracy needs checking against the existing
2 ns results — not just "does it survive".

### 23. THE SEMI-IMPLICIT POISSON IS INCONSISTENT INSIDE A NEWTON RESIDUAL

**The user's insight, and it is the largest effect measured in this whole log.**

The semi-implicit Poisson exists to relax the dielectric-relaxation constraint.
Under Newton it was doing the OPPOSITE -- it was the thing enforcing it.

**Why.** The semi-implicit derivation PREDICTS the new charge density,
rho^{n+1} ~ rho^n - dt*div(sigma*E), and substitutes that into Gauss's law. That
substitution is what produces the operator

    div( (eps + dt*sigma) grad(phi) ) = -rho^n - dt*(diffusive)

and its whole purpose is to avoid needing rho^{n+1}, which a SEGREGATED solver
does not have. **A Newton residual already has rho^{n+1}** -- it is an unknown of
the coupled system, rebuilt from the trial densities at step 1b. So the
`dt*sigma` term accounts for the charge relaxation a SECOND time. The coupled
system then has no consistent root, which is exactly a residual floor, and the
size of the spurious term is `dt*sigma/eps` -- the dielectric relaxation ratio
itself.

**That predicted the ceiling quantitatively, and the data matched.** Across 3167
steps in two runs, NO step whose dielectric relaxation ratio exceeded 0.3278
ever converged. In one run the separation was perfect: 0 of 161 successes at or
above the lowest failing ratio.

**The fix is to use the plain Poisson residual**, since Newton makes the
semi-implicit device redundant:

| | semiImplicit | **explicit** |
|---|---|---|
| SNES failures | 4.92% | **0% (0 of 19)** |
| retries | 10 | **0** |
| dt | ~3.4e-12 | **7.40e-11 (22x)** |
| max dielectric relaxation ratio | 0.60 | **19.52** |

**Newton now steps over the dielectric relaxation time by ~19x** -- which is the
claim JFNK was adopted for and had never once demonstrated.

**Why every other hypothesis failed to find this**, and it is worth recording as
a lesson: the defect was not in the solver, the preconditioner, the Jacobian, the
initial guess, the tolerances or the chemistry -- all six were eliminated by
measurement. It was in the EQUATION the residual encoded. A Newton solver
converges to the root of whatever residual it is given, and this residual had no
root. Six increasingly sophisticated solver diagnostics could not see that,
because they all presuppose the equation is right.

**CONFIRMED over a longer run, and THE ACCURACY IS VERIFIED.** 152 steps, 151
converged, 0 failures, 0 retries, dielectric relaxation ratio reaching 308. It
also passed t = 1.949675e-06, the wall-flux singularity that killed the
semiImplicit run, and did so genuinely -- 0 clamp firings, 0 warnings.

The accuracy check that matters: explicit Poisson with dt FREE against explicit
Poisson with dt capped at 2e-12 (417 steps, 0 failures), compared at a common
time:

| | dt <= 2e-12 (reference) | dt free (~22x) | difference |
|---|---|---|---|
| `n_e` max | 1.581e+17 | 1.572e+17 | **0.57%** |
| `n_e` min | 6.749e+15 | 6.661e+15 | **1.3%** |
| `Emag` max | 9.707e+04 | 9.657e+04 | **0.52%** |

**The 22x larger step gives the same answer to ~1%.** The speed is not being
bought out of the solution.

**And that resolves the trajectory question.** The two EXPLICIT runs agree with
each other, so the large difference from Picard is NOT a large-dt error -- it is
that the semi-implicit form modifies the Poisson equation by `O(dt*sigma/eps)`,
which was **20-60%** in the Picard runs. Explicit Poisson is the true Gauss law;
the semi-implicit one is the approximation. That strengthens the case, but it
should still be checked against the PUBLISHED benchmark before being claimed as
better physics rather than merely different.

### 24. Eliminated by measurement, in order (all with the knob verified)

Recorded so none of these is re-run. Every one was a plausible candidate.

| hypothesis | how it died |
|---|---|
| accuracy / `errRtol` limit | 100x relaxation changed dt by 0.98-1.04x; `rejection memory` took over (2 -> 8 -> 43) |
| adaptive chemistry (`adaptiveError`) | `implicitRate` gave 6.9% failures against 4.5% |
| F not a function of u | purity check: `||dF|| = 0` exactly |
| F has hysteresis | F(x), F(x'), F(x): `||dF|| = 0` exactly |
| matrix-free differencing step | `mffdErr` 1e-5 / 1e-7 / PETSc default: dt ratios 1.00 / 1.00 / 0.99 |
| initial guess outside the basin | linear extrapolation, equal N: ratios 1.00 / 1.01 / 1.01 |
| "just needs more iterations" | most failing solves need 1,000-53,000 more; a minority need ~33 |

### 25. `grubert_steady` died of a SIGFPE, and it is IN THE PRECONDITIONER MATRIX

The run launched to answer "is steady state now reachable" (the user's step 1)
**did not get there.** It aborted at `t = 1.954433e-6` — close to, but not the
same as, Picard's death at `1.971771e-6`. This is a different failure from
Picard's: dt was **6.9e-11 and set by the PI temporal controller**, so the
explicit-Poisson gain held right up to the crash. Not a dt collapse.

**Located, by gdb over a 14 ns restart from the `1.940229e-06` snapshot** (the
snapshot existed because of rule 41 — this is the rule paying for itself):

```
#0  MatMult_SeqAIJ            <- the arithmetic fault
#1  MatMult
#2  PCApply_FieldSplit_Schur
#3  PCApply
#4  KSPFGMRESCycle
#5  KSPSolve_FGMRES
...
#10 solveWithSNES
#11 snesNewtonSolver::solveOuterStep
```

**This eliminates the whole class of hypotheses I would otherwise have spent
the day on.** It is NOT in the residual evaluation, so it is not a trial
iterate driving a density negative, not `Te = nEps/n_e` at a vanishing `n_e`,
and not a rate-table lookup out of bounds — every one of which was a live
candidate before the backtrace, and `bounded false` made the first of them
look likely. `MatMult` only multiplies and adds; it cannot divide. So the
assembled **preconditioner matrix already carried a non-finite or overflowing
entry** before PETSc ever touched it.

Two candidate sources remain, and they need different fixes:

1. one of the Jacobian blocks assembled in `solveOuterStep` produces inf/NaN;
2. the value is manufactured INSIDE the Schur application — `S = A_tt −
   A_tφ·A_φφ⁻¹·A_φt` — most plausibly the inner `A_φφ⁻¹` solve reaching inf
   before the `MatMult` by `A_φt`. Note only `-pc_fieldsplit_schur_fact_type
   full` is set; every inner KSP/PC is at PETSc defaults.

**A permanent scan now discriminates them** (`snesNewtonSolver.C`, immediately
before the `SnesPmatCOO` handoff): one pass over the COO arrays, reporting the
BLOCK (`d(rowField)/d(colField)`), the cell and the value, and refusing
fatally. Deliberately fatal rather than a step rejection — a non-finite
Jacobian entry is an assembly defect, not the physical stiffness `retryStep`
exists to absorb, and rejecting the step would have hidden it. Had this scan
existed, the gdb run would not have been necessary.

**MEASURED, and it eliminates candidate 1.** The fault reproduced (step 241,
t = 1.9537e-6) with the scan reporting **zero non-finite entries** — the scan
string verified present in the library that ran, per the build trap below.

Then the magnitude, because `MatMult` can raise `FE_OVERFLOW` with every input
finite and an isfinite-only test would report "clean" either way:

```
Pmat: 88710 entries, max|a| 5.1e5 ... 2.7e7 ... 1.3e7
                     ALWAYS at d(e)/d(ePotential), local cell 3
```

`max|a| ~ 1e7`, stable, over hundreds of steps. **So the matrix is not the
problem in either mode.** Overflow at `a ~ 1e7` needs a vector element above
`~1e301`, so the enormous quantity is the VECTOR, manufactured inside the
Krylov/Schur solve — candidate 2. An arm with
`-fieldsplit_phi_ksp_converged_reason` (split names are `phi`/`transport`, NOT
`0`/`1` — rule 42) is measuring whether an inner solve diverges before the
fault; through step 11 all 37,396 `phi` and 1,796 `transport` inner solves
report `CONVERGED_RTOL`.

Incidental, and worth keeping: the largest Jacobian entry in this problem is
always the electron/potential drift coupling at the same cell. That is the
block whose absence made the Schur complement inert (section 10).

**A one-line `max|a|` report per outer step is now permanent.** A fatal
threshold alone cannot distinguish "my entries were fine" from "the threshold
was set too high" — only the trend can, and it is one line against the ~550
this solver already writes per step.

**Note on the pre-existing `diag/V` SIGFPE diagnostic** in
`localEnergyEnergyModel.C`: it is for a DIFFERENT fault. Its own criterion is
"a diag/V near zero is the SIGFPE"; here it read 3.08e11 and 3.55e11. Not this.

### 25b. Explicit vs semi-implicit Poisson on GRUBERT: 110-400x the timestep

Six arms, all restarted from the same `1.940229e-06` snapshot. The Poisson
scheme was confirmed **from each run's own log line**, not from the config
(`"Model: singleRegionPoisson Poisson scheme: ..."`), because these cases set
it through a `$poissonScheme` variable and rule 42 applies.

| Poisson | case | steps | reached t | final dt |
|---|---|---|---|---|
| **explicit** | `grubert_steady` | 261 | 1.954433e-06 | **6.90e-11** |
| **explicit** | `poisson_std_long` | 214 | 1.953210e-06 | **6.14e-11** |
| semiImplicit | `newton_sg` | 5639 | 1.947073e-06 | 5.64e-13 |
| semiImplicit | `newton_cfs` | 5037 | 1.946779e-06 | 4.76e-13 |
| semiImplicit | `newton_tolBC` | 5396 | 1.947275e-06 | 3.25e-13 |
| semiImplicit | `extrap_off` | 11433 | 1.949345e-06 | 1.70e-13 |

**An independent reconfirmation of section 23 on a different case** — and much
larger here (110-400x) than the 22x measured on the streamer. Note the step
counts: the explicit arms got FURTHER in 214-261 steps than the semi-implicit
ones did in 5,000-11,400.

**It also explains why only the explicit arms hit the SIGFPE.** The fault is at
a STATE (t ~ 1.9535-1.9545e-6); the semi-implicit arms simply never advanced
far enough to reach it. `poisson_std_long` was killed at 1.953210e-06, a few
hundred steps short. So the FPE is not evidence against explicit Poisson — it
is a defect the explicit arms are the first to be fast enough to expose.

All six were KILLED, not crashed (zero `FOAM FATAL`, zero PETSc FPE, all
stopping within 4 s of each other at 11:20:5x).

### 26. Two of my own claims were wrong, both now corrected in the tree

Recorded because both were committed, and both would have misled the next
reader rather than merely being private mistakes.

| claim | what is actually true |
|---|---|
| the legacy Townsend fit is calibrated for a REDUCED field, so it underflows at real discharge fields (comment in `plasmaTransport.C`, twice) | **Raw V/m is exactly what its constants expect**, for air at 1 atm: α = 19 cm⁻¹ at the 30 kV/cm breakdown field, 1.06e3 cm⁻¹ at 100 kV/cm, `mu = 2.398·E^-0.26` → 0.036 m²/V/s at 1e7 V/m, η = 3.4 cm⁻¹. All correct. The fit is gas- and pressure-SPECIFIC, not miscalibrated. The 2026-09-09 measurement that motivated gating it off (S_iz ~250 orders low) stands, but because it was taken in a 100 Pa argon glow — outside the fit's regime — not because of a field normalisation error. `d4db586` |
| the Newton photoionization refusal (`8c7e572`) was verified | It was keyed on `photoionization_.valid()`, and that autoPtr is **always** valid — a case with no photoionization holds the null object `noPhotoionization`, TypeName `"none"`. The predicate was a constant true, so the guard **refused every case in existence** under Newton. Caught by it refusing `grubert_steady`, a 100 Pa argon glow with no photoionization key at all. `e07d0bb` |

**The lesson, and it is the same one as rule 42:** `8c7e572` verified only that
the guard FIRES (by stripping a chemistry dict until it did). It never checked
that the guard stays SILENT when it should — which is the half that was broken.
A guard has two directions and both are part of the test. Both are now verified
for the photoionization one: `none` → silent, Newton runs; `nTermHelmholtz` →
fires at step 0.

**A build trap found the same way:** `wmake` in `src/numerics` builds
`libplasmaNumerics` and silently does NOT build `libplasmaNewtonSolverPETSc`,
which has its own `Make/files` one directory down. It exits 0 having done
nothing to the Newton library. Caught only by checking the library mtime and
grepping the binary for the new string. **For this library, verify the artefact,
not the exit code.**

## Still untested / next, in priority order (as of 2026-09-11 ~03:40)

1. **`maxIt` is state-dependent — resolve it.** Raising `maxIt` 50→200 HURT at
   the easy pre-ignition state (−16%), but the deep-avalanche run
   (`newton_ignition_ew`, t≈1.90e-6) is failing with `DIVERGED_MAX_IT` at
   `maxIt 50`, which is the opposite signal. The two are not contradictory —
   more iterations waste time on a step that should be retried, EXCEPT when the
   step is genuinely solvable and just needs them. Worth an `maxIt` sweep AT the
   hard state rather than the easy one.
2. **Per-CELL scaling.** `sX` is still one scalar per FIELD while `Arp` spans
   1e11–5.7e18 within its block. Deliberately NOT attempted overnight: 42 usage
   sites, and a wrong refactor would silently corrupt every result above. Do it
   with a switch and an A/B, in daylight.
3. **`d(Psrc)/dφ`, Joule heating's response to the field**, still missing from
   the energy row. Direct and probably strong, but not laplacian-shaped, so it
   needs its own derivation rather than a copy of the drift term.
4. **Pmat staleness.** It is assembled once per outer step and held for the
   whole SNES solve. At ignition the state moves fast within a step; a mid-solve
   refresh is untested.
5. **log(n).** Everything measured points the same way — the failures are
   dynamic-range failures — and it fixes (2) and the clamp at once.

## What is running (overnight, 2026-09-11)

| case | configuration | purpose |
|---|---|---|
| `newton_ignition_ew` | coupling + ionisation deriv + EW | best config, through ignition |
| `newton_ignition_fixed` | coupling only | control |
| `newton_retry_ignition` | oldest library | control |
| `ballast400_fine` | PICARD | the reference trajectory + fine snapshots |
| `hyp_ew_best` | EW + `-mat_mffd_type ds` | best option stack, more samples |

Ignition is at t ≈ 1.97e-6. Picard's own behaviour there is already measured and
is the thing to beat: **dt collapses to ~1e-15 with 31–42 correctors/step.**

**Compare at EQUAL SAMPLE SIZE and from a COMMON start time.** Four separate
misreadings in this session came from comparing runs at unequal step counts or
different start points, in both directions.
