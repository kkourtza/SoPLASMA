# PROGRESS — the live state. READ THIS FIRST.

**Last updated: 2026-09-12.**

This is the ONE orientation file (rule R0). It holds what is *happening*.
`docs/CAPABILITIES.md` holds what *exists* — the inventory, the tooling, the full
refuted list — and this file points at it rather than restating it. Keep both true
(D5); they are only worth re-reading if they are current.

`/orient` reads this for you. `/save-state` updates it.

---

## 1. WHERE WE ARE RIGHT NOW

**One thread is open: does the Newton/JFNK outer solver pay against the segregated
Picard sweep?** Overnight 2026-09-12 the "why is Newton slow" half was ANSWERED, and
it is not a bug in SoPlasma — it is a PETSc default nobody overrode.

**Newton's convergence rate IS its Eisenstat-Walker forcing term eta**, measured:
`||F_{k+1}||/||F_k||` equals the reduction the LINEAR solve delivered, to 1.5% median
over 3.5 decades. So no quadratic term is being lost; an inexact Newton converges
linearly at rate eta by definition. eta sits at ~0.85-0.9 because PETSc's EW
`rtol_max` default is **0.9** and `snesNewtonSolver.C:1664` sets only `-snes_ksp_ew`,
overriding none of the parameters. eta is PINNED at that ceiling for **68.7%** of
Newton iterations, because EW v2 recomputes it from the LINE-SEARCH-DAMPED norm
ratio, so one damped step resets the ladder.

**TWO OPTIONS, NO CODE CHANGE, -55% COST** (verified twice, independently, both arms
at an identical 110 steps): `petscOptions "-snes_ksp_ew_version 3
-snes_linesearch_minlambda 1e-3"` takes 33.5 -> 12.4 SNES its/step and 376 -> 169
cost/step. See [[newton-rate-is-the-EW-forcing-term]] for the caveats — it is NOT a
blanket win and its behaviour under `adjustTimeStep true` is a PREDICTION, untested.

**AND THE BENCHMARK'S CENTRAL ECONOMIC CLAIM NOW HAS A DIRECT MEASUREMENT.** At a
common physical window, dt=1e-9 against dt=2e-10 on the same bed: **cost per ns of
simulated time is identical to ~1%**, while the big step converges 8.3% of its steps
and the small one 100%. **The larger timestep buys nothing** — the extra Newton work
exactly cancels the longer step.

**NEXT CONCRETE ACTION: finish the warm-started ladder on the streamer bed.** Picard's
half is DONE (below). Newton's half needs arms that actually produce a result — the
`kspMaxIt 200` arms gave `DIVERGED_ITS` at 200 with 0 usable steps. A one-variable
pair (shipped vs EW-v3+minlambda) is running at dt=2e-10 with `kspMaxIt 1000`.

### The warm-started ladder, Picard half — MEASURED 2026-09-12

Coarse bed (81,640 cells, `$HOME/streamer-warm`), warm start t=1e-09 from a developed
streamer (peak n_e 1.19e19, tau = eps0/(e mu_e n_e) = 1.16e-10 s), limiters OFF,
common endpoint t=2e-09:

| dt | Picard |
|---|---|
| 1e-11 | reached 2e-09, 100/100 converged in 2 correctors |
| 5e-11 | reached 2e-09 |
| 1e-10 | reached 2e-09, 10/10 converged in **2 correctors** |
| 2e-10 | **SIGFPE step 1** |
| 5e-10 | **SIGFPE step 1** |

**The wall is between dt=1e-10 and 2e-10 for BOTH solvers, and we do NOT know what it
is.** Read from the solver's own instrument at a later step of the SURVIVING arm
(dt=1e-10): `Diel. relax. ratio` **4.79**, `Co_conv (e)` **2.13**, `Co_chem` **5.10**,
temporal error **0.528** against a target of 1, and the accuracy controller asking for
a **24% LARGER** step. So Picard runs with all three Picard-era limiters exceeded at
once, converging in two correctors — **those limiters are conservative by ~2-5x**,
which is itself the benchmark's central experiment answered. It is NOT the dielectric
time, NOT accuracy, and not any single Courant number. Do not fill the gap with a
story; two explanations have already failed.

At 449k Picard fails at 5e-11 instead — that bed is 2.35x finer. **Absolute dt ceilings
do NOT transfer between beds.** The final benchmark number must come from 449k.

## 2. IN FLIGHT

**Nothing is running.** `ps` shows no `soPlasmaFoam`/`mpirun` (checked 2026-09-11).

**But there is ~1 h of measurement on disk that nobody has read**, written after the
last physics commit (`aed5dbe`, 17:50) at 18:47. It is in no commit, no `COMPARE.md`
verdict and no memory. Two arms were cut off MID-RUN:

| arm | state when it stopped |
|---|---|
| `validation/lad_newton_dt1e-10_k1000` | died at step 1, `MPI_ABORT` errorcode 59 |
| `validation/lad_newton_dt2e-11_k1000` | t=1.2e-10, 7 steps, 1 `DIVERGED_ITS`→`DIVERGED_LINEAR_SOLVE` |
| `validation/lad_newton_dt5e-11_k1000` | t=1e-10, 3 steps |
| `validation/pcfix_maxit_dt1e-9` | complete, unrecorded — see §4 Task 1 |
| `validation/pcfix_maxit_dt1e-8` | complete, unrecorded |

All on the 449,413-cell bed, limiters off, fixed dt. **Read these before re-running
anything** — a completed measurement that gets silently redone is the cost R0 exists
to prevent.

---

## 3. BLOCKED

| item | blocked on | what unblocks it |
|---|---|---|
| Task 1, Newton-vs-Picard verdict | every prior verdict invalidated by the ksp cap; and the cold-start ladder cannot reach the stiff regime | the warm-started ladder (§4) |
| Task 3, the AP / semi-implicit-Poisson paper | Task 1's dt data — its key measurement IS the AP claim's central experiment | Task 1 |
| The 1.15M-cell wall under Newton | cause unknown; the one proposed cause is refuted (§5) | a preconditioner answer, not a case change |
| PCSHELL becoming the default | its 17 failures are all `-6 DIVERGED_LINE_SEARCH`, undiagnosed | diagnose the line search |
| Dielectric cases under Newton | the per-patch surface-charge GUARD is unwritten — a dielectric case does not refuse, it **silently drops the surface charge** | write the guard (correctness hole, not an enhancement) |

**SoPlasma is now PUSHED — resolved 2026-09-12.** All 208 commits are on `origin` at
`37b6b6e`; the tree is 0 ahead / 0 behind. The Newton/JFNK solver is off this machine.

**SoEEDF is now pushed too — resolved 2026-09-12.** `fix/mechc-audit-unitarity-bound` had
never been pushed at all; it now tracks `origin/fix/mechc-audit-unitarity-bound`, and
`master` went up with it. No rewrite was needed there — its largest blob is a 46 MB
Literature PDF, well under the limit. **Both trees are now 0 ahead / 0 behind with clean
status.**

---

## 4. TASK LEDGER

Newest first. `[x]` carries the date it was completed.

### Task 1 — Newton-vs-Picard benchmark on `positiveStreamer_fixedMesh` *(user's order, agreed 2026-09-10)*
Status: **IN FLIGHT, unanswered.** Needs a `COMPARE.md`. Already parametrised; no new code.

- [x] 2026-09-11 — `kspMaxIt` exposed as a `newtonSolver` dict key, default 100 unchanged so nothing moves silently (`f17a489`).
- [x] 2026-09-11 — cap measured on `grubert2009_pseudo` (2000 cells, dt=1e-10, 400 steps, 398 SNES solves, identical otherwise): **100 → 42/398 failures (10.6%), 38 of them `DIVERGED_LINEAR_SOLVE`; 1000 → 4/398 (1.0%), ZERO linear-solve failures, +2.5% wall clock.** The 39 solves needing >100 needed 104–115 and **all converged** — 0.92% of 4,235 solves, each 4–15 iterations short.
- [x] 2026-09-11 — the "preconditioner cliff"/bimodal-KSP diagnosis **RETRACTED**: the distribution looked bimodal only because the cap truncated it. True distribution median 4, p90 14, max 115 (`aed5dbe`).
- [x] 2026-09-11 — PCSHELL (`assembledPmat false`) measured **better**: KSP median 3 / p90 6 / max 22 vs fieldsplit 4 / 14 / 115 (5x tighter tail), 8% fewer residual evaluations (113,208 vs 123,234).
- [ ] **Read the five unrecorded arms in §2 and write their verdict into a `COMPARE.md`.**
- [ ] **THE KEY EXPERIMENT — the warm-started ladder.** Measure the largest dt each solver survives at fixed accuracy with the Picard-era limiters off.
- [ ] Re-test needleDBD multi-region with the cap raised (its run predates the key).

**Standing instructions on this task:** leave `changeDictionary` authoritative (it
produced the validated 2 ns LFA/LMEA results and both arms share it, so the comparison
stays fair) and do **not** resolve the BC fork here. Measure **per ns of simulated time
at equal accuracy, never per step** — Newton costs ~4–5x per step, so a per-step
comparison is rigged for Picard (C3).

### Task 2 — surface charge / dielectrics under Newton, in `tutorials/.../needleDBD`
Status: **unblocked as of 2026-09-11, not finished.** A correctness hole, not an enhancement.

- [x] 2026-09-11 — multi-region blocker REMOVED (`10b7b2c`): phi lives in a ragged COO tail, **no monolithic assembly needed**. Measured why: the explicit `fvc::` residual already sees the neighbouring region through `coupledElectricPotential`'s Robin condition (interface-cell d(res) rms 1.549 vs 0.00354 one cell in), so Newton's own outer iteration converges the region coupling.
- [x] 2026-09-11 — five gates passed, three defects caught that were invisible without running: empty Pmat tail rows (zero diagonal → FPE in the first `PCApply`); a **fifth** `fvSchemes` catch-all, `"snGrad(n_.*)"`; and parallel refused loudly rather than silently building a preconditioner from disconnected subdomains.
- [ ] **Write the per-patch surface-charge guard.** Until it exists a dielectric case under `outerSolver newton` silently drops the surface charge.
- [ ] Re-run `validation/needle_mrgate` with `kspMaxIt` raised — its `DIVERGED_LINEAR_SOLVE` at 17:35 predates the key landing at 17:49.
- [ ] Handover problem, still open: `anySpeciesOnFloor()` gates Picard→Newton on every species being off its floor, and in needleDBD all 13 start AT their floor — so Newton is unreachable until ignition. The gate seeds all 13 above their floors, which is **MACHINERY, NOT PHYSICS** and must never be quoted for DBD results.

### Task 3 — the AP proof / semi-implicit-Poisson paper
Status: **blocked on Task 1.** Publishable claim: the standard semi-implicit AP
reformulation is INCONSISTENT inside a fully-coupled Newton residual, and the coupled
solve attains AP on its own. Villa paper in `SoEEDF/Literature`. (E4: research is a
valid track, not a lower-priority one.)

### Environment / infrastructure *(this session, 2026-09-11)*
- [x] 2026-09-11 — `validation/` run output gitignored on the `smoke/` principle; **`git status` 1310 lines / 182 s → 1 line / 0.6 s**. Enabled `core.untrackedCache` + `feature.manyFiles`.
- [x] 2026-09-11 — `verification/` TRACKED (was untracked in its entirety — one `git clean` from losing the only analytic-ground-truth suite).
- [x] 2026-09-11 — `ROUNDW.{C,H}` committed (untracked source, in `Make/files`, 2nd-best of five ROUND variants).
- [x] 2026-09-11 — needleDBD `Allrun-serial` portability fix committed.
- [x] 2026-09-11 — CLAUDE.md rewritten and moved here; SoEEDF's is now a stub importing it.
      **475 → 369 lines while gaining 12 rules and the command card**, because the procedures
      moved into skills. All 31 rule ids and all 44 old numbers still resolve (checked by grep).
- [x] 2026-09-12 — `.claude/` built: 9 skills, 10 agents, 3 hooks, and `.claude/mined/` (the
      adversarially-verified extraction of 129 memories and ~90 documents across nine
      dimensions). The verification pass rejected **19 of 212** extracted commands — including
      one copied verbatim FROM a memory, so that memory is wrong too (`Starting time loop` is
      never printed, and the `sed` range then matches the whole log).
- [x] 2026-09-12 — memory store SHARED: the `soplasma-scratch` memory directory is now a symlink
      to the `Projects-SoEEDF` one, so 129 memories resolve from both roots and the base cannot
      fork. `SoEEDF/.claude/{skills,agents,hooks,mined}` symlink here too.
- [x] 2026-09-12 — **the push was unblocked.** It had been failing because **7 blobs exceeded
      GitHub's 100 MB hard limit** (up to 317.5 MB), all `validation/*/logs/log.soPlasmaFoam`,
      baked into 151 of the 208 unpushed commits. The `.gitignore` work stopped them
      accumulating but git pushes HISTORY, so the whole push was rejected. Purged all run
      output from the unpushed range with `git filter-repo` in a throwaway clone:
      **2094 MB → 39.7 MB, 0 blobs over 50 MB.** The rewritten chain still descends from
      origin's tip, so it went up as a **normal fast-forward, no `--force`**, and the working
      tree was never touched — the solver running at the time was undisturbed.
      **The purge set is derived, not hand-written**: purge = output-shaped AND not present in
      the current tree, so tree identity holds BY CONSTRUCTION. A first attempt that
      hand-wrote the path list silently dropped 414 tutorial `plasmaTables` files (breaking
      `positiveStreamer_LMEA_fast`, the 2 s debugging bed) and 14 `testSnesJFNK` mesh files —
      caught only because the tree hash was compared before and after. Verified identical:
      `ce3db1a6…` before and after.
      152 commits got new SHAs; the 102 hash citations across 42 docs and memories were
      rewritten from filter-repo's commit map in the same change (D1/D2), leaving 0 stale.
- [x] 2026-09-12 — SoEEDF pushed: the feature branch (53 commits, never pushed) and `master`
      (3 ahead). Four uncommitted items committed first, and **one deliberately NOT committed**:
      `Cross section_IST_Lisbon_He.txt`. `data/lxcat/cross-sections/` is gitignored because
      **LXCat does not authorise redistribution** (`data/README.md`), and their policy names
      "commercial interests, in particular" — which is exactly what SoPhy is headed for. The
      two sets already tracked are grandfathered exceptions; do not add more. The file was
      moved into `cross-sections/` so it is correctly ignored instead of showing as noise.
- [x] 2026-09-12 — **first real run of `/regression-gate`, and the 2D bed contradiction settled.**
      Tier 1 (analytic unit beds) PASS, all five exact strings. Tier 2 (electrostatics) PASS,
      `6 ok, 0 failed`, 16.6 s / 91 MB. Tier 3: BOTH baselines were STALE, not regressions —
      each verified four ways (mtimes, an `nCorr` field the baseline lacked, the error having
      IMPROVED, and two runs bit-identical) then refreshed as separate stated acts (`bcd7288`,
      `f110bac`). **No regressions anywhere.**
      The 2D bed's headline — "SG AND CFS DO NOT CONVERGE ON A NON-ORTHOGONAL MESH" — is
      **REFUTED by its own data**: SG order 0.98, CFS 1.01 against the control's 0.97. The
      original 2026-09-08 finding was correct AND HAS BEEN FIXED (`ScharfetterGummel.H:140`
      now uses `nonOrthDeltaCoeffs` + explicit `snGrad`, justified by the Bernoulli
      factorisation rather than patched); the README documented a defect that no longer
      existed. Its "DO NOT use SG or CFS on skewed meshes" line was **actively misleading
      guidance** and is corrected. Still UNVERIFIED for genuinely unstructured or graded
      meshes — one uniform 11.3° shear is not a claim about a tetrahedral mesh.
      Also found: the 2D `Allrun` hardcoded its scheme loop and **could not reproduce the
      `standard` control row its own README quoted** (A1). `SCHEMES` now defaults to all three;
      the control reproduces to every digit.
- [x] 2026-09-12 — **the WSL freezes are memory, measured not guessed.** Newton on the 449k bed
      is **5.2–7.1 GB RSS per arm**; Picard on the same bed is **1.56 GB** — Newton costs ~4x the
      memory, which is a genuine C3 cost result and not just a wall-clock question. Four arms =
      20.3 GB of a 30 GB WSL cap (no `.wslconfig`, so the default ~50% of RAM). The eight arms
      running during both freezes ≈ 40 GB, i.e. over the cap, and the VM swap-thrashes rather
      than cleanly OOM-killing — which is why the terminal froze and would not reopen.
      **Practical limit: at most 4 large-bed arms, no more than 3 of them Newton.**
      The regression suite was NOT the cause: 91 MB peak, 16.6 s idle.
- [~] **The regression gate is HALF done.** `/regression-gate` now carries the procedure and the
      classification B5 requires (REGRESSION vs INTENDED IMPROVEMENT vs STALE BASELINE), but it is
      a skill I execute — **there is still no script and no CI hook**, so nothing compares
      `results.txt` to `results.baseline.txt` unattended. Finishing it means a runnable comparator
      with a relative tolerance. Until then the gate only fires when someone invokes it.

### Also open, unordered *(CAPABILITIES §6)*
- [ ] `grubert_1d_I`, relaunched cold 2026-09-11 — the **first genuine test of current control**, because `particleFlux_` was fossilised under Newton so the regulator had never actually worked. Its earlier state came from ~3900 steps of UNREGULATED ramping and is not trustworthy as physics.
- [ ] The negative-`L` chemistry guard (`max(L,0)`) — cheap; fold in when next touching `plasmaChemistryODE`.
- [ ] Validation-suite regeneration, and the live BC fork (below).
- [ ] `plasmaChemistry0D.C:953` still hardcodes `100.0` as `meanEnergyMax` (`:1004` clamps to it) after the model-side default was derived from the table range.
- [ ] The second memory-consolidation pass.
- [ ] The Grubert physics thread (paused, not resolved) — the model reaches the Townsend→glow transition and fails at **sheath formation**. Excluded by measurement: the circuit, the density floor, the energy tables, the meanE clamp, resolution/grading. Remaining suspects: the wall-flux closure at a forming sheath, and the `n_e`/`nEps_e` coupling (implicated three times).
- [ ] **Outer-coupling conditioning — the largest measured performance lever**, promoted from optimisation to blocker for any steady solver. On `grubert2009_ballast_clean`: 22 correctors/step (tail 141, cap 150), Aitken omega 0.148, and the dt governor keyed on that same omega clamped dt to **1/83** of what accuracy allowed — compounded, ~1000x more work. Temporal error names `nEps_e` as the worst field, so suspect the ENERGY coupling first.

**THE BC FORK IS LIVE AND DELIBERATELY FROZEN.** In the four `positiveStreamer` beds
`etc/changeDictionary` runs AFTER the generators and overrides them, and the two paths
disagree on PHYSICS: charged species at the electrodes are `electron/ionDDWallFluxMixed`
(ABSORBING) under layer 1 vs `zeroGradient` (NON-absorbing) under changeDictionary;
`n_e` at `far` is `zeroGradient` vs `inletOutlet inletValue 1e13`. Frozen on 2026-09-10
because changeDictionary produced the validated 2 ns results. Do not "fix" it inside
Task 1.

---

## 5. FAILED APPROACHES — DO NOT RETRY

### RETRACTED 2026-09-12: "the wall sits at the dielectric relaxation time"

Claimed repeatedly on 2026-09-12, in `PROGRESS.md`, both design notes, a memory, and
in the brief sent to the literature review as a possible NOVEL result. **It is wrong.**

`tau` was computed from an ASSUMED electron mobility (0.04 m^2/Vs) rather than the
solver's own `maxSigma`, which `plasmaTimeControl.C:1818` prints every step as
`Diel. relax. ratio = deltaT*maxSigma/eps0`. The error was **~20x**. Quoted dt/tau
values of 0.86 / 1.7 / 4.3 were really **17.3 / 34.6 / 86**.

**Picard converges in TWO correctors at a dielectric ratio of 4.8-17** — i.e. it steps
5-17x past the dielectric relaxation time with EXPLICIT Poisson coupling — so that
limit is simply not the binding constraint for an outer loop that ITERATES the
coupling to convergence. A second explanation (quasi-neutrality suppressing the net
charge where sigma is largest) is ALSO unsupported: peak |chargeDensity| is 44% of
q*n_e, not a quasi-neutral channel.

**What survives:** the Eisenstat-Walker finding and the 3.7x pseudo-transient speedup,
which rest on SNES/KSP counts and not on tau. And the literature verdicts, which stand
on their own texts.

**The lesson, and it is the expensive one:** the solver was PRINTING the correct
quantity every step while I quoted a derived one I never checked against it. That is
A1's "the control is NAMED... report every magnitude against the physical scale it
should be judged by", failed on a number I then built four turns of argument on.


### Newton's linear convergence — four one-variable refutations (2026-09-12)

All against `validation/diag_A0_control` at an identical endpoint. Each knob was
PROVEN to move before its null result was believed.

- **A non-smooth clamp in the residual.** The mean-energy floor was moved
  0.0388 -> **1e-9 eV** (the banner confirms the 7.6-order move) and the result was
  IDENTICAL IN EVERY DIGIT over 41 steps. This was the leading hypothesis and it is
  dead. The `nEps_` write-back at `localEnergyEnergyModel.C:719-723` IS a real
  non-smoothness and still fires — it is a symptom of a bad step, not the cause.
- **Flux-scheme non-smoothness** — `Gauss ROUNDF` -> `Gauss linear`: refuted; the
  collapse step moves 44 -> 42.
- **The chemistry per-cell routing branch** — `adaptiveError` -> `implicitRate`:
  identical in every digit.
- **Tightening the matrix-free differencing** — `mffdErr` 1e-8 makes it WORSE (58
  `DIVERGED_LINE_SEARCH`): the matvec relative error becomes eps_F/err = 25%. The
  usable band is eta ~2.5e-4 .. 0.9.
- **`adaptiveForcing false`** does not fix it, it RELOCATES it: PETSc's default KSP
  rtol 1e-5 is BELOW F's 2.5e-9 matvec noise floor, so the solve cannot converge and
  you get `DIVERGED_LINEAR_SOLVE` instead.

**`-snes_ksp_ew_monitor` DOES NOT EXIST in PETSc 3.24.** An arm built on it was
byte-identical to its control (222,906 lines both) and measured nothing. The forcing
term is visible only via `PETSC_OPTIONS="-info :snes"` in the ENVIRONMENT — PetscInfo
is consumed at `PetscInitialize`, so it cannot go in the case dict. **Third unverified
PETSc option to cost real time in one day; verify with `-options_left` and confirm the
monitor actually prints before building anything on it.**

**`validation/diag_*` cold-started is an INVALID bed for solver diagnosis.** n_e is
pinned at `minNumberDensity` 1e11 for the whole run — the clamp sits outside the
equations so F=0 is unreachable (`newton-outer-solver-design.md` DEFECT A) — and at
dt=1e-9 it runs at Co_conv(e) 20-37, Co_conv(energy) 30-55 with `adjustTimeStep
false`. Its 42% failure rate is a property of the bed, not of Newton.


Each is measured. Re-proposing one costs the measurement again. Full detail and dates
in `docs/CAPABILITIES.md` §4/4b/4c.

**Preconditioning and scaling** — *tried X, didn't work because Y:*
- **ILU(1) on the transport split** — helped at 2k (−33% steps), *hurt* at 20k (+3.5%). A 2,000-cell artefact.
- **`selfp` Schur approximation** — same shape: 2k −22%, 20k KSP median 11→12.
- **ILU(1) + `selfp` together** — they INTERFERE; the combination is exactly baseline.
- **Per-cell residual/state scaling** — 2k −46% steps, 20k +1.4% mean dt, KSP median 13 vs 13.
- **`sourceAwareScaling`** — REMOVED 2026-09-11. `rms(chemP)/(rms(n)/dt)` is 0.833 for all seven source-dominated species: the source NEVER exceeds ddt, so `max(ddt, src)` is a no-op. A species produced from nothing has n ~ P·dt by construction, so n/dt *is* the source scale.
- **More `-fieldsplit_phi_ksp_max_it` for the 1.15M wall** — 8 cycles vs 2 at 449k: −5% SNES, −4% Krylov, **+14% wall clock**. A net loss.
- **`extrapolateGuess` for block imbalance** — 38 vs 40 converged solves over 11 steps. No effect here.
- **Ion mobility as the cause of Newton's conditioning trouble** — once the chemJacobian sign bug was fixed, mobile and immobile arms are identical digit for digit. The earlier apparent advantage WAS the bug.

**Grubert's DC-glow operating point** — every time-marched route is closed:
- **Voltage + ballast** — at the CVC minimum dR/dI = 0, so the load line is TANGENT and selects nothing (Almeida & Benilov 2017, plus ~10 of our own arms). All died at 150/150 correctors with dt collapsing ~25,000x. Includes "ballast-limited ignition". *Switched to: current control.*
- **Current control, time-marched** — sound but insufficient; all six `grubert2009_iset*` arms stalled at t~9.5e-7, one after 80,983 steps with 23,575 discards.
- **A gentler ignition ramp** — helps, does not solve: ignited quasi-statically at −137.7 V and cut rejections 30,547→1,727, then still overshot to 718x setpoint.
- **`ddtSchemes steadyState` as a case setting** — DOES NOT EXIST. The solver refuses it: ddt IS the diagonal, and removing it leaves rows with no diagonal.
- **`relaxationFactors` by hand** — BREAKS the working recipe, 9472/9472 converged → 0/10. SoPLASMA already runs adaptive Aitken outer relaxation; fixed factors fight it. *Switched to: the pseudo-transient recipe (BDF2 + `adjustTimeStep false` + fixed dt + currentSource, NO relaxationFactors), which works today with no code — but has only ever been run to 9.5 ns of a μs–ms problem, so it is STABLE BUT UNPROVEN, not closed.*

**Numerics and diagnostics:**
- **`fvMatrix::residual()` standalone in parallel** — misreports by ~21 orders of magnitude, independent of PETSc. *Switched to: explicit `fvc::`.*
- **`relativeChange` outer criterion** — REMOVED 2026-09-06, now fatal. It divides by a nearly-uniform field's deviation about its own mean; and its raw `break` skipped the loop contract, so `updateChargeDensity()` ran on 151 of 364,670 steps (0.04%).
- **GMRES (not FGMRES) on the transport split** — `DIVERGED_BREAKDOWN` from a varying operator. *FGMRES tolerates it, and is mandatory.*
- **`ROUND*01` schemes for densities** — they clamp the field to **[0,1]**, catastrophic for 1e16 m⁻³ — and "bounded" is exactly what one reaches for.
- **Tightening the Poisson tolerance for the lateral asymmetry** — improves phi/Ex, leaves `Ey` (the driving component) unchanged, and GAMG at 1e-16 costs 2000 iterations (286x) for nothing. *Switched to: `NY = 1` for 1-D, normal refinement otherwise.*
- **A mesh-symmetrisation utility or symmetry check** — rejected with the user: refinement already cures it (`Ey ~ NY^-1.7`), and an unstructured triangular mesh has no mirror symmetry at all, so a check would fire constantly and mean nothing.
- **`chemJacobian`'s P/n for ALL species** — exact only for the electron; for a species produced by electron impact on something else the true derivative is ~0 while P/n is enormous, and it inverted the diagonal's SIGN. Now clamped and electron-only by default.
- **The negative-`L` extrapolation mechanism** — RETIRED 2026-09-11: `plasmaRateTable.C` fits a POWER LAW to the last interval, which cannot cross zero.
- **The ngspice bridge** — decided against.

**PERMANENTLY OFF THE LIST BY USER INSTRUCTION: solve for log(n) instead of n.**
*"remove log(ne) from the deferred list or demote it to last and only if i say so"*
(2026-09-10). Do not propose it, schedule it, or treat it as a prerequisite. The interim
workaround — restart Newton from an established Picard state, where the clamp is
inactive — works.

**Already answered, do not re-run:** the 1.15M-cell mesh under Newton "to see if it
works". Answered both ways — Newton converges at 40k, 211k, 449k and does not at 1.15M,
proven to be MESH SIZE by running it warm (10 h CPU/rank, still step 1) and cold (3 h
CPU/rank), both at 99.9% CPU, i.e. computing, not deadlocked.

---

## 6. DISCOVERED TASKS

Found during other work; not yet scheduled.

- [~] **The regression gate is half-built as of 2026-09-12** — `/regression-gate` carries the
      procedure and the three-way classification, but no script exists, so this entry stays open.
      Originally found as: **the regression gate does not exist.** `verification/fluxScheme1D/` and `fluxScheme2Dnonortho/` each ship a `results.baseline.txt` and **nothing in the tree ever compares them** (`grep -rn baseline verification/` finds no caller). `fluxScheme2Dnonortho` already differs from its baseline — and it is a **STALE BASELINE, not a regression** (`nCorr` 2→15/16, errors ~2x lower, written three minutes later). Nothing on disk can tell those apart. *(Found 2026-09-11. This is what B5 requires and `/regression-gate` is meant to be.)*
- [ ] **`DANGLING: simulationType` in needleDBD** — a real live defect. The dangling-key detector prints exactly one line on that case today. A `$key` nothing references is inert, and a dangling `appliedVoltage` once invalidated a whole "low-field" arm that ran at the original 62.1 Td and exited 0. *(Found 2026-09-11.)*
- [ ] **`./Allwmake` is NOT a full build and says "Build complete." anyway.** `src/models/plasmaModels/Allwmake` builds 4 of 8 model directories; `src/applications/Allwmake` builds 4 of ~15. `./build-all.sh` is the only correct build. Either fix Allwmake or make it refuse. *(Found 2026-09-11.)*
- [ ] **The repo's `build.log` is a stale 2026-08-10 artefact** with 14 `OK` lines against a `BUILD_DIRS` list that now has 31 entries. Anything grepping it reads a month-old build as current. Delete it or regenerate it. *(It is gitignored, so it only misleads locally.)*
- [ ] **`report.py` silently drops 72 of 504 rows** in the fluxScheme1D study — `std:ROUNDW` and `CompleteFlux` are absent from its line-29 scheme list. Add them.
- [ ] **`check-no-running-solvers.sh` guards only 7 binaries** and misses `singleRegionElectrostaticFoam`, `multiRegionElectrostaticFoam`, `plasmaSetupBoundaries`, `testEmission`, `testWallFlux`, `testDischargeCurrent`. A build over any of those is unguarded.
- [ ] **`testAitken` has no pass/fail machinery** — it always exits 0 and prints no verdict. It is an exploratory bed, but it sits among the test utilities where its exit code looks like a pass.
- [ ] Push both trees. 205 unpushed commits here; SoEEDF's branch has no upstream.

---

## 7. RUNNING NOTES (when stuck)

Dated attempt log, so the next session inherits the search rather than repeating it.
Clear entries into §4 or §5 once resolved.

*(Empty — the ksp-cap thread resolved into §4. Open a dated block here the moment an
investigation takes more than two attempts, and record each refuted hypothesis as you
refute it, not afterwards: on 2026-09-11 four hypotheses were refuted one after another
and only the instruments settled it — A6.)*

---

## 8. POINTERS

| for | read |
|---|---|
| what EXISTS — inventory, tooling, full refuted list | `docs/CAPABILITIES.md` |
| the rules, and the command card | `CLAUDE.md` (already in context) |
| why a rule is trusted — the measured post-mortems | `../Projects/SoEEDF/docs/rules-postmortems.md` |
| every boundary role and circuit type | `etc/boundaryRoles` |
| every dictionary key | `docs/reference/` |
| a model's physics and options | `docs/models/**` |
| design notes and status — **check the `Status:` line, several are DESIGN only** | `docs/design/*.md` |
| the mined knowledge base, by dimension | `.claude/mined/*.md` |
| ~129 accumulated memories | the memory `MEMORY.md` index; newest for this thread are `newton-vs-picard-benchmark-state` and `ksp-cap-was-the-bug`, which are NOT `session-state-*` files |

**Memory ordering is not filename ordering.** For the Grubert thread,
`session-state-2026-09-07c` supersedes `-09-07b` supersedes `-09-07` supersedes `-09-06`,
and `-09-08` continues `-09-07c`. The newest state of the *current* (Newton) thread is
not a session-state file at all.
