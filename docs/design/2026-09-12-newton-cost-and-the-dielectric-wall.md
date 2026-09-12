# Newton's cost, and the wall at the dielectric relaxation time

**Status: MEASURED. Two fixes ready to ship. THE dt/tau FRAMING IN THIS NOTE'S TITLE AND
SECTION 3 WAS WRONG AND HAS BEEN CORRECTED IN PLACE (2026-09-12); the retraction is in
`PROGRESS.md` section 5.**
Dated 2026-09-12. Supersedes nothing; it is the measurement record behind
`newton-outer-solver-design.md`, which remains the design record.

Companion: `2026-09-12-preconditioning-coupled-poisson-drift-diffusion-literature.md`
carries what the literature already says. This note carries what WE measured.
Keep them apart — that one is citations, this one is our own numbers.

---

## 1. The headline

Two independent things were established, and they point in opposite directions.

**(a) Newton is ~2x more expensive than it needs to be, and the cause is a PETSc
default nobody overrode. Two options fix it, no code change.**

**(b) The larger timestep Newton exists to enable was measured to be worth
NOTHING on our bed — and the literature says that is exactly what should happen
under our conditions.** If that holds, a better outer solver buys ROBUSTNESS,
not speed, and the benchmark's motivation has to be restated.

---

## 2. Newton's convergence rate IS its Eisenstat-Walker forcing term

**Tier 2 (instrumented convergence study).** Bed: `grubert2009` 2000 cells, LMEA,
dt=1e-9, `kspMaxIt 1000`, `adjustTimeStep false`.

`||F_{k+1}||/||F_k||` **equals the reduction the LINEAR solve actually
delivered**, to 1.5% median over 3.5 decades. Measured over 88 full-step Newton
iterations (lambda=1, so no line-search confound):

| requested eta | achieved eta | KSP its | SNES ratio |
|---|---|---|---|
| 0.00013 | 0.00011 | 120 | 0.00011 |
| 0.00219 | 0.00149 | 116 | 0.00150 |
| 0.03281 | 0.03093 | 17 | 0.03105 |

**So no quadratic term is being "lost" — it is unreachable by construction.** An
inexact Newton converges linearly at rate `eta` by definition. The observed
~0.85 is simply `eta`.

**Why eta sits there:** PETSc's EW `rtol_max` default is **0.9**, and
`snesNewtonSolver.C:1664` sets only `-snes_ksp_ew`, overriding no parameters.
eta is PINNED at that ceiling for **2,454 of 3,573 (68.7%)** Newton iterations,
because EW v2 recomputes it from the LINE-SEARCH-DAMPED norm ratio: one damped
step resets eta to the ceiling and the ladder must be re-descended.

For damped steps the same law holds in first-order form, `ratio = 1 -
lambda(1-eta)`, median error 7.0e-4 across FIVE decades of lambda — which also
proves **the Newton direction is a correct first-order descent direction.**

### The fix

    petscOptions "-snes_ksp_ew_version 3 -snes_linesearch_minlambda 1e-3";

* **EW version 3** (Chacon) caps rtol at `rtol_0` = 0.3 instead of `rtol_max` =
  0.9 and adds an anti-oversolving guard, so a damped step cannot reset the ladder.
* **minlambda 1e-3**: PETSc's default is 1e-12 with `bt`'s Armijo constant 1e-4,
  so at lambda=1e-5 the required decrease is 1e-9 relative — **below F's own
  measured noise floor of 2.5e-9**. About 10% of line-search decisions were being
  made on noise. The derived floor is ~2.5e-5; 1e-3 carries a 40x margin.

Verified INDEPENDENTLY (two separate counts, within 2%), both arms at an
identical 110 steps:

| | SNES its/step | KSP/step | cost/step |
|---|---|---|---|
| shipped (EW v2) | 33.5 | 251 | **376** |
| v3 + minlambda | 12.4 | 129 | **169 (-55%)** |

**NOT a blanket win.** v3 ALONE oversolves once the line search collapses (619
vs 383 cost/step); it needs `minlambda` with it. And it costs 3 converged steps
in the healthy window — under `adjustTimeStep true` those become `retryStep` +
a dt cut, which is a **PREDICTION, not a measurement** (the bed was fixed-dt).
**Test that before shipping it as a default.**

---

## 3. The wall, and what it is NOT

Bed: coarse streamer, 81,640 cells, warm-started at t=1e-09 from a DEVELOPED
streamer, all segregated-era limiters OFF, common endpoint t=2e-09.

| dt | Picard | Newton |
|---|---|---|
| 1e-11 | reached 2e-09, 100/100 converged in 2 correctors | — |
| 5e-11 | reached 2e-09 | — |
| 1e-10 | reached 2e-09, 10/10 converged in **2 correctors** | converges, \|\|F\|\| 700 -> 1.2e-4 in 32 its |
| 2e-10 | **SIGFPE step 1** | **no linear solve converges** |
| 5e-10 | **SIGFPE step 1** | **no linear solve converges** |

**THE WALL IS BETWEEN dt = 1e-10 AND 2e-10 FOR BOTH SOLVERS, AND WE DO NOT KNOW
WHAT IT IS.** It is NOT the dielectric relaxation time, NOT accuracy, and not any
single Courant number — see below. Do not fill this gap with a story; two
successive explanations have already failed.

### It is NOT the dielectric relaxation time

An earlier version of this note claimed the wall sat at dt = tau "because tau is
where dt*sigma overtakes eps0". **That was wrong**, and wrong for a mundane
reason: tau was computed from an ASSUMED electron mobility (0.04 m^2/Vs) instead
of the solver's own `maxSigma`, which it prints every step. The error was ~20x.

Read from the instrument instead (`Diel. relax. ratio` = `deltaT*maxSigma/eps0`,
`plasmaTimeControl.C:1818`), on the arm that SURVIVES:

    dt = 1e-10, a later step:
      Diel. relax. ratio : 4.79     <- ~5x PAST the dielectric limit
      Co_conv (e)        : 2.13     <- past the drift CFL
      Co_chem            : 5.10     <- past the chemistry limit
      temporal err       : 0.528  [target 1]
      "dt would be: 1.237e-10"      <- the accuracy controller would allow 24% MORE

**Picard runs with all three Picard-era limiters exceeded simultaneously —
dielectric x4.8, drift CFL x2.1, chemistry x5.1 — converging in TWO correctors,
with temporal error at half its target.** That is itself the benchmark's central
experiment ("remove the Picard-era limiters and see how far each solver can
actually step") and its answer: **those limiters are conservative by roughly
2-5x, and the segregated solver is fine without them.**

sigma is also not constant: the ratio reads 17.3 at the first step of that arm and
falls to 4.8 as the streamer evolves, so any single dt/tau number is a snapshot,
not a property of the bed.

### Three preconditioners, same wall

| preconditioner | result at dt = 2e-10 |
|---|---|
| fieldsplit + default `a11` Schur | `DIVERGED_ITS` @ 1000 |
| physics-based PCSHELL (`assembledPmat false`) | `DIVERGED_ITS` @ 1000, residual frozen |
| `schurOnPhi true` | 66 min, 566 evals, ZERO completed solves, residual unmoved |

**Supplying the Schur operator did not rescue it either.** But see the caveat
below before reading that as evidence about the operator.

### Two caveats on the Schur story, both from the literature

* `div((eps + dt*sigma) grad .)` is the Schur complement **only under
  `A_tt^-1 ~ dt I`**, which Chacon & Knoll 2003 p.581 condition on
  `dt <~ dt_A`, the ADVECTIVE CFL limit. Calling it "the exact Schur complement"
  or "the CORRECT S_f" — as earlier versions of this note did — overstates it.
* **The prescribed mitigation was never applied to the transport split.** Our own
  `schur-semiimplicit-poisson-preconditioner.md` lines 142-156 pre-registered this
  failure — "an inner solve taking a varying number of iterations makes the Schur
  operator NOT a fixed linear operator" — and prescribed `richardson` / `max_it 5`
  / `convergence_test skip`. The shipped code applies that to the PHI split and
  leaves transport on `fgmres`/`rtol 1e-2`/`max_it 200`. Under `schurOnPhi` it is
  the TRANSPORT block being inverted inside the Schur operator, so the arm may
  have failed for a pre-registered implementation reason rather than a physical
  one. **Re-test with the mitigation before concluding anything about S_f.**

## 4. The premise in doubt

At a COMMON physical window (t in [4.4e-8, 5.5e-8]), same bed, dt=1e-9 vs
dt=2e-10: **cost per ns of simulated time identical to ~1%**, while the big step
converged **8.3%** of its steps and the small one **100%**.

**The larger timestep bought nothing** — the extra Newton work exactly cancelled
the longer step.

**CAVEAT, and it is serious:** those two arms converged 8.3% vs 100% of steps, so
that comparison is across two different algorithms, not one algorithm at two dt.
It needs redoing at equal convergence before it is quotable.

**The literature predicts it** (Teunissen PSST 29 (2020) 015010 §3.6): beating
tau is beneficial only when `tau < tau_CFL`. The solver reports both itself, and
on this bed at dt=1e-11: `Diel. relax. ratio` 0.0927 (so tau = 1.08e-10) against
`Co_conv (e)` 0.1726 (so tau_CFL = 5.8e-11). **tau is ~2x LARGER than tau_CFL —
the CFL binds first, and beating tau cannot help.**

**THE DECISIVE EXPERIMENT, NOT YET RUN:** repeat the dt sweep on a bed where
`tau < tau_CFL` (coarser cells raise tau_CFL; `tau_CFL ~ dx`). A 16,900-cell bed
is staged at `$HOME/streamer-cfl` for exactly this. Until that is done, the
benchmark's central claim — that a larger dt is worth having — is UNTESTED under
the one condition where the literature says it should pay.

---

## 5. Refuted by one-variable test — do not re-propose

Each knob was PROVEN to have moved before its null result was believed.

* **A non-smooth clamp in the residual.** Was the leading hypothesis. The
  mean-energy floor was moved 0.0388 -> **1e-9 eV** (the banner confirms the
  7.6-order move) and the result was IDENTICAL IN EVERY DIGIT over 41 steps.
* **Flux-scheme non-smoothness** (`Gauss ROUNDF` -> `Gauss linear`).
* **The chemistry per-cell routing branch** (`adaptiveError` -> `implicitRate`):
  identical in every digit.
* **Tightening the matrix-free differencing.** `mffdErr` 1e-8 makes it WORSE (58
  `DIVERGED_LINE_SEARCH`): matvec relative error becomes `eps_F/err` = 25%.
* **`adaptiveForcing false`** relocates rather than fixes: PETSc's default KSP
  rtol 1e-5 is BELOW F's 2.5e-9 matvec noise floor. Usable band: eta ~2.5e-4 .. 0.9.
* **A "preconditioner cliff"/bimodal KSP distribution.** An artefact of
  `-ksp_max_it 100` TRUNCATING the distribution. True distribution: median 4,
  p90 14, max 115.

---

## 6. Traps this work paid for

* **`-snes_ksp_ew_monitor` DOES NOT EXIST in PETSc 3.24.** An arm built on it was
  byte-identical to its control (222,906 lines both) and measured nothing. The
  forcing term is visible only via `PETSC_OPTIONS="-info :snes"` in the
  ENVIRONMENT — PetscInfo is consumed at `PetscInitialize`, so it cannot go in
  the case dict.
* **A malformed option SEGVs inside PETSc.** `-ksp_monitor_true_residual
  ::ascii_info_detail` crashed every arm; removing it, the same arm ran 8 solves
  and 152 residual evaluations. **Verify with `-options_left` AND confirm the
  monitor prints.** Three unverified PETSc options cost real time in one day.
* **`pdftotext` silently corrupts equations in scanned papers.** Read
  equation-bearing pages as rendered IMAGES.
* **`writeControl runTime` with `writeInterval 10` means TEN SECONDS of simulated
  time.** A run ending at 1e-09 s wrote ZERO snapshots while every grep-the-key
  check passed.
* **Absolute dt ceilings do NOT transfer between beds.** The drift Courant number
  scales with `1/dx`, so the 449k bed (2.35x finer) fails at a 2.35x smaller dt
  than the 81k bed. Only the RATIO transfers.
* **This bed is 1.4-7.6 GB resident per solver on a 30 GB machine.** Ten
  concurrent 449k arms froze WSL twice.
* **`validation/diag_*` cold-started is an INVALID bed for solver diagnosis:**
  n_e pinned at `minNumberDensity` 1e11 all run (clamp outside the equations, so
  F=0 unreachable — `newton-outer-solver-design.md` DEFECT A), Co_conv(e) 20-37,
  Co_conv(energy) 30-55. Its 42% failure rate is a property of the bed.

---

## 7. What to do next, in order

1. **Ship the two EW options** — after testing them under `adjustTimeStep true`,
   which is the one untested condition.
2. **Run the `tau < tau_CFL` bed** (`$HOME/streamer-cfl`, 16,900 cells). This
   decides whether a larger dt is worth anything at all, and therefore whether
   any further preconditioner work is motivated by SPEED or only by ROBUSTNESS.
3. **Redo the dt cost comparison at equal convergence** — the 8.3%-vs-100%
   confound above.
4. Only then, preconditioning above tau. `schurOnPhi` is the one candidate with a
   principled reason to differ, and its arm has not yet completed a solve.
