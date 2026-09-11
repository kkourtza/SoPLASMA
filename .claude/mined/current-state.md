# dim6

## Summary
This dimension is the live state of work as of 2026-09-11 evening: what is done, what was cut off mid-flight, what is blocked, and what a PROGRESS.md must never let anyone re-attempt. The project is in the middle of ONE thread — does the Newton/JFNK outer solver pay against the segregated Picard sweep — and that thread was re-opened hours before the session ended: a hardcoded `-ksp_max_it 100` was found (commits f17a489, aed5dbe, 2026-09-11 17:49/17:50) to have produced EVERY negative Newton verdict on record, so the 449k head-to-head ladder, the pseudo-transient dt ceiling and the needleDBD multi-region result are all INVALIDATED and awaiting re-test. Five re-test arms exist on disk (`validation/pcfix_maxit_dt1e-9`, `_dt1e-8`, `lad_*_k1000`) whose results are NOT in any commit, doc or memory — two of them were still writing at 18:47 and were cut off mid-run. The second live fact is mechanical: `git status` in the live tree takes 182 s and prints 1310 lines, of which 1109 are experiment output, so every "small testable commits" rule is dead on arrival until `validation/` run output is ignored the way `smoke/` already is. The Grubert DC-glow thread is settled-as-closed for every time-marched route except the pseudo-transient recipe, which has only ever been run to 9.5 ns of a microsecond-to-millisecond problem. No PROGRESS.md, TODO.md or equivalent exists in either tree; `docs/CAPABILITIES.md` (§4/4b/4c/5/6) is the de-facto inventory and PROGRESS.md must point at it rather than restate it.

## Facts

### NOTHING IS RUNNING and the newest work is uncommitted. `ps` shows no soPlasmaFoam/mpirun process. Last commit is aed5dbe at 2026-09-11 17:50; the newest run output was written at 18:47 (`validation/lad_newton_dt2e-11_k1000.log`, `lad_newton_dt5e-11_k1000.log`). So ~1 h of measurement exists on disk with no commit, no COMPARE.md verdict and no memory entry — and two of those arms were cut off MID-RUN (7 and 3 steps).
**Evidence:** `ps -eo pid,etime,cmd | grep soPlasmaFoam` → empty; `git log -1 --format='%h %ad'` → aed5dbe 2026-09-11 17:50; `ls -lt /home/kkourtza/soplasma-scratch/validation/*.log | head -2` → 18:47

**Rule:** PROGRESS.md must open with a dated 'AS OF' block listing: nothing running, last commit hash+time, newest run-output mtime. On session start, run the three commands in exact_commands #1-#3 and reconcile — output newer than the last commit means a result nobody has read.

**Cost:** The ksp-cap re-test — the measurement that decides whether every negative Newton verdict stands — gets silently redone or, worse, the stale verdict gets quoted.


### THE CAP WAS THE BUG, and it is the single most important live fact. A hardcoded `-ksp_max_it 100` (`snesNewtonSolver.C:1530`) killed near-converged linear solves. On `grubert2009_pseudo` (2000 cells, dt=1e-10, 400 steps, 398 SNES solves, identical otherwise): kspMaxIt 100 → 42/398 failures (10.6%), 38 of them DIVERGED_LINEAR_SOLVE; kspMaxIt 1000 → 4/398 (1.0%), ZERO linear-solve failures, +2.5% wall clock. The 39 solves that needed >100 needed 104-115 and ALL converged — 0.92% of 4,235 solves, each 4-15 iterations short.
**Evidence:** I reproduced the counts from the raw logs: `grep -c 'Nonlinear solve did not converge' validation/pcfix_base_dt1e-10.log` → 42 (converged 356); `pcfix_maxit_dt1e-10.log` → 4 (394); `pcfix_shell_dt1e-10.log` → 17 (381). Memory `ksp-cap-was-the-bug.md`; commits f17a489, aed5dbe; docs/CAPABILITIES.md:401-472

**Rule:** Treat every Newton/JFNK verdict recorded before 2026-09-11 17:49 as INVALID until re-measured with `kspMaxIt 1000`. Never quote a failure COUNT without the `due to <REASON>` breakdown.

**Cost:** A whole benchmark's negative conclusion (Newton does not pay) rests on an expiring budget misread as a solver failure; the measured correction is 10.5x fewer failures for +2.5% wall clock.


### THE PCSHELL IS THE BETTER PRECONDITIONER and its remaining defect is named but undiagnosed: `assembledPmat false` gives KSP median 3 / p90 6 / max 22 against the fieldsplit's 4 / 14 / 115 (5x tighter tail) and 8% fewer residual evaluations (113,208 vs 123,234), but its 17 failures are ALL `-6 DIVERGED_LINE_SEARCH`. Fix the line search and it should beat everything.
**Evidence:** `grep -o 'due to [A-Z_]*' validation/pcfix_shell_dt1e-10.log | sort | uniq -c` → 5722 CONVERGED_RTOL, 284 CONVERGED_FNORM_RELATIVE, 97 CONVERGED_SNORM_RELATIVE, 17 DIVERGED_LINE_SEARCH, zero DIVERGED_LINEAR_SOLVE. Commit aed5dbe; memory ksp-cap-was-the-bug.md

**Rule:** Put 'diagnose the PCSHELL line-search failures (`assembledPmat false`)' on the task list as the cheapest remaining Newton win. It was previously unreachable in production because `snesBridge` prefers PCFIELDSPLIT whenever a Pmat is supplied.

**Cost:** The stronger preconditioner stays off by default because its headline failure RATE looked worse (4.3% vs 1.0%) while its linear algebra was better.


### UNRECORDED RESULT ON DISK #1 — the pseudo-transient dt ceiling re-test. With kspMaxIt 1000 at dt=1e-9 (400 steps to t=4e-07, 2000 cells): 168/398 SNES failures (42.2%), and the failure mode has MOVED — 125 DIVERGED_MAX_IT + 41 DIVERGED_LINE_SEARCH but only 2 DIVERGED_LINEAR_SOLVE. At dt=1e-8: 56/67 failures (83.6%), 70 steps to t=7e-07. So raising the cap did NOT simply buy 10x on dt; it converted a linear-solve failure into a genuine NONLINEAR one (SNES `maxIt 50`).
**Evidence:** `grep -c 'Nonlinear solve did not converge' /home/kkourtza/soplasma-scratch/validation/pcfix_maxit_dt1e-9.log` → 168 (converged 230); reason histogram as quoted; case config `validation/pcfix_maxit_dt1e-9/configuration/config:96 deltaT 1e-9`, `:109 outerSolver newton`; mesh note `nCells:2000`

**Rule:** Record this in PROGRESS.md as a MEASURED-BUT-UNWRITTEN result and re-read the two logs before proposing any pseudo-transient dt. The memory's claim that 'dt=1e-9 is now worth testing' HAS been tested; the answer is 42% at 1e-9, not 1%.

**Cost:** Someone re-runs a 400-step arm that already exists, or plans a ~1 h steady run on a 10x-dt assumption the data does not support.


### UNRECORDED RESULT ON DISK #2 — the 449k streamer ladder re-run with kspMaxIt 1000 was IN FLIGHT when the session ended. `lad_newton_dt1e-10_k1000` died at step 1 with MPI_ABORT errorcode 59; `lad_newton_dt2e-11_k1000` reached t=1.2e-10 (7 steps, 1 DIVERGED_ITS→DIVERGED_LINEAR_SOLVE); `lad_newton_dt5e-11_k1000` reached t=1e-10 (3 steps). All on the 449,413-cell bed with limiters off and fixed dt.
**Evidence:** `validation/lad_newton_dt*_k1000.log` mtimes 18:47/18:47/17:59; `ls -d validation/lad_newton_dt2e-11_k1000/[0-9]*` → 2e-11 4e-11 6e-11 8e-11 1e-10 1.2e-10; `constant/polyMesh/owner` note `nCells:449413`; config line 90 `deltaT 2e-11`, 71/74 `limitSpeciesCo false` `limitChemistryCo false`, 67 `adjustTimeStep false`

**Rule:** Do not read these three arms as a completed ladder — they are truncated. Either finish them or delete them, and write the COMPARE.md first (none of the three has one).

**Cost:** A truncated ladder gets quoted as a verdict, repeating exactly the error that CAPABILITIES 4c already had to retract once.


### THE BENCHMARK'S KEY EXPERIMENT HAS STILL NEVER BEEN RUN. Every ladder so far is COLD-started and therefore cannot answer the question by construction: cold start reaches t=2e-10 with peak n_e still at the 1.3e13 seed against a streamer head's ~1e20, so tau=eps0/(e*mu_e*n_e) ~1e-4 s and dt/tau ~ 2e-7 — seven orders BELOW the stiffness Newton exists to step over. Picard confirms it: 2-4 correctors against a cap of 20, Aitken INACTIVE. The required experiment is the WARM-STARTED ladder from a developed streamer.
**Evidence:** memory newton-vs-picard-benchmark-state.md; docs/CAPABILITIES.md:361-400; CAPABILITIES 3b 'Newton needs a successful PICARD step before it hands over… handover LATCHES'

**Rule:** Next concrete action for task 1: warm-start both solvers from the established t=1e-9 streamer state on `positiveStreamer_fixedMesh`/`scale_r4`, kspMaxIt 1000, limiters off, and measure wall clock PER NS OF SIMULATED TIME AT EQUAL ACCURACY — never per step (Newton costs ~4-5x/step, so per-step is rigged for Picard).

**Cost:** Two full ladders (seven arms + three k1000 arms) have already been run that answer nothing; a third would too.


### THE COST BAR NEWTON MUST CLEAR, and the caveat on it. Measured at 449k, dt 1e-11, streamer bed: Picard 14.6 s/step vs Newton 142 s/step (~10x). Warm-started Picard at 1.15M survives dt=5e-11, so Newton must hold ~5e-10 to break even — BUT Picard's COLD ceiling at 449k is 2e-11 (SIGFPE step 1 at both 5e-11 and 1e-10), a different bed and start, so the bar must be restated against a like-for-like control. Newton is NOT less robust in dt: at 5e-11 Picard dies outright and Newton completes 4 steps.
**Evidence:** docs/CAPABILITIES.md:314-360; memory newton-vs-picard-benchmark-state.md

**Rule:** State the break-even bar in PROGRESS.md with its control attached ('~10x, measured on 449k warm vs 1.15M warm — NOT like-for-like'). Any new bar must name the bed, the start (cold/warm) and the dt.

**Cost:** A break-even number quoted against the wrong control decides a solver architecture.


### THE 1.15M WALL IS MESH SIZE, and its one proposed cause is REFUTED. Newton converges at 40k, 211k, 449k and does not at 1.15M — proven both warm (10 h CPU/rank, still step 1) and cold (3 h CPU/rank, stuck after handover), both at 99.9% CPU. The 'starved Schur complement' hypothesis (more `-fieldsplit_phi_ksp_max_it`) was tested at 449k, 10 steps at dt=1e-11: 2 cycles = 77 SNES/387 Krylov/121 s per step; 8 cycles = 73/373/138 s — -5% SNES, -4% Krylov, +14% wall clock, a net LOSS. Cause still unknown.
**Evidence:** docs/CAPABILITIES.md:314-338 and the REFUTED row; commit abcd8d2 'Refute the phi-split hypothesis for the 1.15M wall' (2026-09-11 15:15)

**Rule:** Never re-run 1.15M under Newton 'to see if it works' — answered twice, both ways. The open question is the PRECONDITIONER, not the case. Keep the inner phi count FIXED (`convergence_test skip`): a convergence-tested inner solve makes the Schur operator non-linear in b, which is what raised the original SIGFPE.

**Cost:** A 10 h/rank run that has already been done twice.


### `git status` COSTS 182 SECONDS AND PRINTS 1310 LINES in the live SoPlasma tree; `git status --porcelain -uno` costs 0.58 s and prints 193; scoped to source (`src docs tools etc build-all.sh`) it costs 0.00 s and prints 0 — the SOURCE TREE IS CLEAN and every dirty path but 8 is experiment output.
**Evidence:** `/usr/bin/time -f '%e' git status --porcelain` → 181.97 s, 1310 lines; `-uno` → 0.58 s, 193 lines; scoped → 0.00 s, 0 lines (measured 2026-09-11 in /home/kkourtza/soplasma-scratch)

**Rule:** Put `git status --porcelain -uno` in CLAUDE.md as THE status command for this tree, and scope any untracked check to a path. A bare `git status` exceeds the default 120 s Bash tool timeout and will look like a hung tool.

**Cost:** Every session pays 3 minutes or a timeout for a status check, and the 'commit in small testable units' rule is unusable because the signal is 0.6% of the output.


### `validation/` IS NOT GITIGNORED — only `smoke/` run output is. 3,228 files are TRACKED under validation/: 1,650 `constant/plasmaTables`, 325 under `0/`, 223 `constant/polyMesh`, 154 inside other time directories, 154 `postProcessing`/`logs`/`log.*`, and just 26 COMPARE.md. On top of that 1,109 untracked paths and 129 tracked-but-deleted files.
**Evidence:** `cat /home/kkourtza/soplasma-scratch/.gitignore` (smoke/* block only); `git ls-files validation/ | wc -l` → 3228; breakdown by `awk -F/ '$3=="0"'` etc.; `git status --porcelain | awk '{print $1}' | sort | uniq -c` → 1117 ??, 129 D, 64 M

**Rule:** THE FIX, mirroring the existing smoke/ block exactly: add `validation/*/0/`, `validation/*/[0-9]*.[0-9]*/`, `validation/*/[0-9]*e-[0-9]*/`, `validation/*/logs/`, `validation/*/postProcessing/`, `validation/*/processor*/`, `validation/*/log.*`, `validation/*.log`, `validation/*/constant/polyMesh/`, `validation/*/system/fvSolution`; then `git rm -r --cached` the ~850 already-tracked run-output paths and commit the 129 deletions in the same pass. KEEP tracked: `COMPARE.md`, `configuration/config`, `configuration/boundaries`, `0.orig/`. Decide `constant/plasmaTables` deliberately (1,650 files; regenerable by sweep, but shipping them keeps cases fast — the same ship-vs-sweep inconsistency deferred-action-items already flags for the tutorials).

**Cost:** 1,310-line status, 182 s per check, and the .gitignore's own recorded history: things invisible to `git status` rot (all five smoke beds were broken and unnoticed by 2026-09-02).


### THE `[0-9]*` GLOB TRAP IS ALREADY PAID FOR — twice. The .gitignore says in its own comment that time-directory patterns deliberately do NOT use `[0-9]*`, which matches `0.orig` and would ignore the initial-field TEMPLATES: 'That glob has already destroyed 0.orig twice here; the failure surfaces much later as "cannot find file .../0/ePotential"'.
**Evidence:** /home/kkourtza/soplasma-scratch/.gitignore, smoke-bed block comment

**Rule:** Any new ignore rule for time directories must use the paired forms `[0-9]*.[0-9]*/` and `[0-9]*e-[0-9]*/` plus an explicit `0/`, never bare `[0-9]*`.

**Cost:** Destroyed 0.orig templates, surfacing much later as a misleading missing-field error.


### 197 COMMITS ARE UNPUSHED in the live SoPlasma tree and 52 in SoEEDF. The entire Newton/JFNK outer solver exists only on this machine plus these local commits. Untracked-and-unbacked-up on top of that: `verification/` (4 beds including `fluxScheme1D`, a 504-run/14-scheme sweep with REPORT.md and figures, and `testSnesJFNK`/`testSnesJFNK2Field`), 4 `ThirdParty/libROUNDSchemes/src/ROUNDW.*` files (the ROUNDW scheme added from upstream 2026-09-08), and `ThirdParty/petsc-3.24.0/`.
**Evidence:** `git rev-list --left-right --count origin/feature/chemistry_integration_and_BoltzmannSolver...HEAD` → 0 197; SoEEDF `origin/master...HEAD` → 0 52; `git status --porcelain | grep -v validation/` → 8 paths; `ls verification/` → fluxScheme1D fluxScheme2Dnonortho testSnesJFNK testSnesJFNK2Field

**Rule:** PROGRESS.md needs a standing 'BACKUP STATE' line (unpushed commit count + untracked new directories), refreshed with the fast `-uno` status. `ThirdParty/petsc-3.24.0/` should be ignored, never committed; `verification/` should be tracked (case definitions + REPORT.md) with its run output ignored.

**Cost:** Measured 2026-09-09: two days with no commit left ~6,900 lines of new solver source untracked with no backup of any kind, and the user had to catch it (memory prompt-to-commit-proactively).


### THE REAL UNCOMMITTED CONTENT IS INVISIBLE IN THE NOISE: of the 64 modified tracked files, 26 are `configuration/` case definitions and 14 are COMPARE.md — the pre-registered questions and verdicts that are the project's intellectual record. The 129 deletions are pruned snapshots in `ballastCo_1000` (64), `ballastCo_10000` (64) and `grubert2009_steady` (1), i.e. tracked run output deleted by the rule-41 ten-snapshot pruning.
**Evidence:** `git status --porcelain -uno` → 14 lines matching COMPARE.md (grubert2009_R1e9, R5e9, ballast_std, iset_All150/All1500/En150/En1500, n11_clamp100, n11_clampDerived, steady, …); deletion grouping by `cut -d/ -f2 | sort | uniq -c`

**Rule:** Commit the 14 COMPARE.md + 26 configuration/ changes as content (they are the record of what each arm was asked); resolve the 129 deletions by untracking, not by `git checkout` (restoring stale snapshots would resurrect deleted run output).

**Cost:** The verdicts of ~14 experiments live only in the working tree while 1,109 untracked field files hide them.


### TASK 1 (user's agreed order, 2026-09-10) — Newton-vs-Picard benchmark on `positiveStreamer_fixedMesh`, with a COMPARE.md. STATUS: in flight, unanswered. BLOCKED ON: (a) every prior verdict invalidated by the ksp cap; (b) the cold-start ladder cannot reach the regime. NEXT CONCRETE ACTION: warm-started ladder, kspMaxIt 1000, limiters off, speed measured per ns at equal accuracy. Standing instruction: leave `changeDictionary` authoritative (it produced the validated 2 ns LFA/LMEA results; both arms share it so the comparison stays fair) and do NOT resolve the BC fork here.
**Evidence:** docs/CAPABILITIES.md:618-631; memory deferred-action-items.md 'NEWTON vs PICARD BENCHMARK'

**Rule:** Seed PROGRESS.md's task list in this exact order (1 benchmark, 2 surface charge under Newton, 3 AP/semi-implicit-Poisson paper) and carry the 'blocked-on' + 'next action' fields verbatim; the order was agreed with the user and rule 25 says an 'ok/continue' means THAT thread.

**Cost:** Thread-switching. Memory `never-switch-threads` exists because it happened.


### TASK 2 — surface charge / dielectrics under Newton in `tutorials/plasma/soPlasmaFoam/needleDBD`. The multi-region blocker was REMOVED on 2026-09-11 (commit 10b7b2c, 17:33: phi in a ragged COO tail, no monolithic assembly), and `validation/needle_mrgate` then ran: handover 'Picard warm-up COMPLETE at t = 2e-12 … PERMANENT', then SNES reason -3 (DIVERGED_LINEAR_SOLVE). That run is at 17:35, i.e. BEFORE the `kspMaxIt` key landed at 17:49, and the memory records that its second linear solve hit 100 having already dropped ||F|| 4.2x on the first. It has NOT been re-run with the cap raised. Separately, the per-patch surface-charge GUARD is still unwritten: a dielectric case under `outerSolver newton` does not refuse — it silently drops the surface charge.
**Evidence:** `grep -in 'multiRegion|handover|reason -3' validation/needle_mrgate/log.mrgate2` lines 29/229/963; commit 10b7b2c 17:33 vs f17a489 17:49; memory surface-charge-under-newton-test-in-needledbd.md; docs/CAPABILITIES.md §6 item 2 'dielectric cases under Newton are currently unguarded — a correctness hole, not an enhancement'

**Rule:** Next concrete action: re-run `needle_mrgate` with `kspMaxIt 1000` before concluding anything about multi-region Newton; and land the cheap half (a boundary walk that refuses charging surfaces under Newton) independently of the residual work.

**Cost:** A silently wrong DBD result — surface charge dropped with no error — and a re-litigated 'Newton can't do multi-region' conclusion that was a 100-iteration budget.


### TASK 3 — the AP proof / semi-implicit-Poisson paper. STATUS: blocked on task 1's dt data ('its KEY measurement … is also the AP claim's central experiment, so it must precede the paper').
**Evidence:** docs/CAPABILITIES.md:632-634; memory ap-proof-newton-removes-dielectric-constraint.md

**Rule:** Do not start the paper before the warm-started ladder produces the largest-stable-dt numbers; the paper's central claim IS that measurement.

**Cost:** A publishable claim written ahead of its evidence (rule 33 / rule 11).


### ALSO OPEN, UNORDERED (CAPABILITIES §6): `grubert_1d_I` relaunched COLD 2026-09-11 — the FIRST genuine test of current control under Newton, because `particleFlux_` was fossilised so the regulator had never actually worked under Newton (its earlier state came from ~3900 steps of UNREGULATED ramping and is not trustworthy as physics); the negative-`L` chemistry guard (`max(L,0)`, cheap, fold in when next touching `plasmaChemistryODE`); the validation-suite REGENERATION with its live BC fork; the second memory-consolidation pass; and `plasmaChemistry0D.C:953` still hardcodes 100.0 as meanEnergyMax (`:1004` clamps to it) after the model-side default was derived from the table range.
**Evidence:** docs/CAPABILITIES.md:636-648; `grep -n '100\.0' src/applications/utilities/plasmaChemistry0D/plasmaChemistry0D.C` → 953, 1004; memory session-state-2026-09-07 'Still open, deliberately'

**Rule:** Carry all five into PROGRESS.md with their one-line next action; the 0-D hardcode is a one-line fix that has now survived four session states.

**Cost:** The 0-D utility silently disagrees with the coupled solver's derived clamp (2644 eV vs 100 eV on the LMEA tables).


### THE BC FORK IS LIVE AND DELIBERATELY FROZEN. In the four positiveStreamer beds `etc/changeDictionary` runs AFTER the generators and overrides them, and the two paths disagree on PHYSICS: charged species at the electrodes are `electron/ionDDWallFluxMixed` (ABSORBING) under layer 1 vs `zeroGradient` (NON-absorbing) under changeDictionary; `n_e` at `far` is `zeroGradient` vs `inletOutlet inletValue 1e13`. changeDictionary was left authoritative on 2026-09-10 because it produced the validated 2 ns LFA/LMEA results.
**Evidence:** memory deferred-action-items.md, 'THE FORK THAT MUST BE DECIDED'

**Rule:** PROGRESS.md must name this as a KNOWN FORK with an owner and a decision point (decide, then re-run) — and must repeat the instruction that the benchmark does NOT resolve it.

**Cost:** Two physics models coexist in the published benchmark beds; switching electrodes to absorbing walls would silently invalidate the 2 ns comparison that anchors accuracy.


### REFUTED — PRECONDITIONER AND SCALING (never re-propose): ILU(1) on the transport split (2k: -33% steps; 20k: +3.5% — a 2,000-cell artefact); `selfp` Schur approximation (2k: -22%; 20k: KSP median 11→12); ILU(1)+`selfp` together (they INTERFERE — exactly baseline); per-cell residual/state scaling (2k: -46% steps; 20k: +1.4% mean dt, KSP median 13 vs 13); `sourceAwareScaling` (rms(chemP)/(rms(n)/dt) = 0.833 for all seven source-dominated species and 1.6e-4 for the electron — the source never exceeds ddt, so `max(ddt,src)` is a NO-OP; REMOVED 2026-09-11); more `-fieldsplit_phi_ksp_max_it` as the cure for the 1.15M wall (+14% wall clock for -5% SNES); `extrapolateGuess` as a cure for block imbalance (38 vs 40 converged solves over 11 steps); ion mobility as the cause of Newton's conditioning trouble (mobile and immobile arms identical digit for digit once the chemJacobian sign bug was fixed).
**Evidence:** docs/CAPABILITIES.md:252-275 (table rows, each with its measurement and date)

**Rule:** PROGRESS.md gets a verbatim 'DO NOT RE-ATTEMPT' table with these rows and their numbers. Add the general rule they share: a numerics result measured at 2,000 cells does not transfer — validate at >1 size (SoEEDF rule R43).

**Cost:** Four separate 2k-cell 'wins' that all evaporate at 20k were each implemented before being re-measured.


### REFUTED — GRUBERT DC-GLOW ROUTES (CAPABILITIES 4b, settled 2026-09-11): voltage+ballast CANNOT reach the operating point (at the CVC minimum dR/dI = 0 so the load line is TANGENT and selects nothing — Almeida & Benilov 2017, plus ~10 of our arms at seriesResistor 1e8 and seriesRC 1e6/1e8/1e9/5e9, every one dying at 150/150 correctors with dt collapsing ~25,000x); current control TIME-MARCHED is sound but insufficient (all six `grubert2009_iset*` arms stalled at t~9.5e-07, `grubert2009_iset` after 80,983 steps with 23,575 discards); a gentler ignition ramp helps but does not solve (2.04e7 V/s ignited quasi-statically at -137.7 V and cut rejections 30,547→1,727, then still overshot to 718x setpoint); `ddtSchemes steadyState` as a case setting DOES NOT EXIST — the solver refuses it because ddt IS the diagonal; `relaxationFactors` added by hand BREAKS the working recipe (9472/9472 converged → 0/10) because SoPLASMA already runs adaptive Aitken outer relaxation.
**Evidence:** docs/CAPABILITIES.md:276-308; commits 787c912, 94689f1, 2e1b1d0 (2026-09-11 13:48-14:13)

**Rule:** Reproduce this table in PROGRESS.md. The ONE route still open is the PSEUDO-TRANSIENT recipe (BDF2 + `adjustTimeStep false` + fixed dt + currentSource + NO relaxationFactors), which ran 9472 consecutive converged steps with 0 failures and 4 correctors/step — but only to 9.5 ns of a us-ms problem, at ~5 h per us. STABLE BUT UNPROVEN, not closed.

**Cost:** Recorded cost of not knowing the Almeida result: a full day. The ballast route was then re-proposed and re-run ~10 times.


### REFUTED — NUMERICS AND DIAGNOSTICS: `fvMatrix::residual()` standalone in parallel is BROKEN (misreports by ~21 orders of magnitude — always compute the outer residual by explicit `fvc::`); the `relativeChange` outer criterion was REMOVED 2026-09-06 and is now FATAL (it divides by a nearly-uniform field's deviation about its own mean; and its raw `break` skipped the loop contract — `updateChargeDensity()` ran on 151 of 364,670 timesteps, 0.04%, voiding a day of results); GMRES (not FGMRES) on the transport split gives DIVERGED_BREAKDOWN from the varying operator; `ROUND*01` schemes clamp a field to [0,1] — catastrophic for a density of 1e16 m^-3; tightening the Poisson tolerance to cure the lateral asymmetry (GAMG at 1e-16 costs 2000 iterations, 286x, and leaves `Ey` unchanged); a mesh-symmetrisation utility or symmetry check (rejected with the user — refinement already cures it, `Ey ~ NY^-1.7`); the ngspice bridge (decided against); `chemJacobian`'s P/n applied to ALL species (exact only for the electron; it inverted the diagonal's SIGN); and the negative-`L` extrapolation mechanism (RETIRED 2026-09-11: `plasmaRateTable.C` fits a POWER LAW to the last interval, which cannot cross zero).
**Evidence:** docs/CAPABILITIES.md:252-275; memory raw-break-skips-loop-contract, session-state-2026-09-06, deferred-action-items item 4 'RETIRED 2026-09-11'

**Rule:** Keep each of these as a one-line row with its number. Note the retirement nuance: `max(L,0)` is still cheap insurance against a negative TABULATED value from a mechanism file — a different, unobserved route.

**Cost:** Re-deriving a refutation costs more than reading one; the `relativeChange` regression alone voided an entire day's Grubert results, twice over (45-50% of timesteps then ran ZERO correctors).


### PERMANENTLY OFF THE LIST BY USER INSTRUCTION: solve for log(n) instead of n. 'remove log(ne) from the deferred list or demote it to last and only if i say so' (user, 2026-09-10). Do not propose it, schedule it, or treat it as a prerequisite. The interim workaround — restart Newton from an established Picard state where the clamp is inactive — WORKS.
**Evidence:** docs/CAPABILITIES.md:650-655; memory deferred-action-items.md final section

**Rule:** PROGRESS.md needs a 'NOT ON THE LIST BY USER DECISION' section distinct from 'REFUTED', so a sound idea the user declined is not re-raised as if it were merely unexplored.

**Cost:** Re-proposing something the user explicitly demoted reads as not listening; it also touches the species residual, the preconditioner AND the Picard path, so it is a substantial change, not a patch.


### NO PROGRESS/TODO FILE EXISTS in either tree. `find -maxdepth 2` for *progress*/*status*/*todo*/*INDEX* returns only: `/home/kkourtza/Projects/SoEEDF/docs/INDEX.md` (a documentation index, 149 lines, which points at `docs/framework-state.md` as 'the current snapshot' — that file is dated 2026-09-04 and predates the entire Newton solver), `/home/kkourtza/soplasma-scratch/validation/status.py` (a live-run status tool, not a plan), and `/home/kkourtza/soplasma-scratch/smoke/STATUS` + `STATUS_DT` (dated 2026-08-12, reading `native exit=1 / cantera exit=1 / native-stiff exit=1 / ALL DONE` — STALE, not current health). The repo root `TODO.md` is gitignored and does not exist.
**Evidence:** the find above; `ls -l docs/framework-state.md` → Sep 5 18:44, header '# State of the framework — 2026-09-04'; `cat smoke/STATUS`; `.gitignore` line '/TODO.md'

**Rule:** PROGRESS.md is genuinely new: it must supersede `docs/framework-state.md` for CURRENT STATE (linking to it for the module/default audit), must NOT restate `docs/CAPABILITIES.md` (§4/4b/4c are the refuted/already-run inventory — link to them by section), and must say plainly that smoke/STATUS is a 2026-08-12 artefact.

**Cost:** A fourth overlapping status document. SoEEDF's INDEX.md states the set is deliberately non-overlapping: 'If you find the same explanation in two places, one of them is a bug.'


### THE MEMORY INDEX MUST BE READ IN A SPECIFIC ORDER, and filename dates mislead. For the Grubert thread: session-state-2026-09-07c supersedes -09-07b supersedes -09-07 supersedes -09-06, and -09-08 CONTINUES -09-07c. But the newest state of the CURRENT (Newton) thread is not a session-state file at all — it is `newton-vs-picard-benchmark-state.md` (2026-09-11 16:57) and `ksp-cap-was-the-bug.md` (2026-09-11 17:52), the two most recently written memories.
**Evidence:** `ls -t /home/kkourtza/.claude/projects/-home-kkourtza-Projects-SoEEDF/memory/*.md | head -5` → MEMORY.md, ksp-cap-was-the-bug.md, newton-vs-picard-benchmark-state.md, deferred-action-items.md, session-state-2026-09-08.md; each file's own 'Supersedes [[…]]' line

**Rule:** PROGRESS.md's 'where to read next' block must list, in order: docs/CAPABILITIES.md §4/4b/4c/6 → ksp-cap-was-the-bug → newton-vs-picard-benchmark-state → deferred-action-items → session-state-2026-09-08 → -09-07c. Sort by mtime (`ls -t`), never by the date in the filename.

**Cost:** Reading session-state-2026-09-08 as 'latest' gives the Grubert/CFS thread and misses the ksp-cap invalidation entirely, which is the single fact that changes what to do next.


### THE GRUBERT PHYSICS THREAD'S LAST LIVE STATE (session-state-2026-09-07c/-09-08, superseded for PRIORITY by the Newton thread but not resolved): the model gets INTO the Townsend→glow transition and fails at SHEATH FORMATION — `I_cond` reverses sign while the electrode is still at -123/-282/-255 V, at 0.018-1.9% of I_sc, with Te,max peaking at 37.3 eV and collapsing to 12.1 eV. EXCLUDED by measurement: the circuit, the density floor (1e9 vs 1e7 agree to 0.19%), the energy tables (Joule = Loss + nu_i*U holds to 0.997-1.020 over five decades), the meanE clamp, and resolution/grading (uniform 2000-cell reverses at 1.1241 us vs the graded mesh's 1.121 us — mesh refinement is RETRACTED as an explanation per the case's own pre-committed discipline). REMAINING suspects: the wall-flux closure at a forming sheath, and the n_e/nEps coupling (implicated three times).
**Evidence:** memory session-state-2026-09-07c.md; session-state-2026-09-08.md (CFS delays but does not fix: -200 V bought 0.18 us; SG control ignited 1.7x EARLIER than CFS though SG is provably positivity-preserving, so positivity is NOT the limiting factor)

**Rule:** Keep this as a PARKED thread in PROGRESS.md with its exclusion table intact, so nobody re-litigates the floor, the tables, the clamp or the mesh. Its cheapest next test is a UNIT BED on the n_e/nEps coupling, not another CFD arm.

**Cost:** Five candidate causes each cost multiple multi-hour arms to exclude; re-running any of them buys nothing.


### THE OUTER-COUPLING CONDITIONING ITEM is the largest performance lever measured and has been PROMOTED from optimisation to blocker for any steady solver. On `validation/grubert2009_ballast_clean`: 73,567 PIMPLE iterations over 3,306 steps = 22 correctors/step (tail 141, cap 150), Aitken omega = 0.148, and the dt governor keyed on that SAME omega clamped dt to 1/83 of what accuracy allowed (9.19e-10 vs 7.65e-8) — compounded, ~1000x more work. The temporal error names `nEps_e` as the worst field, so suspect the ENERGY coupling first. Note the semi-implicit Poisson scheme, which exists precisely to condition this, is ALREADY ON.
**Evidence:** memory deferred-action-items.md 'OUTER-COUPLING CONDITIONING'; session-state-2026-09-07b 'Promote from optimisation to blocker'

**Rule:** Carry the exact numbers and the two caveats (weight omega, NOT `rho [contraction]` which reported 18,504 here and is established as unreliable; semi-implicit Poisson is already on, so 'turn it on' is not the fix).

**Cost:** If Newton removes this clamp the win is 80x, which dwarfs its 4-5x per-step cost — the whole economic case for the benchmark rests on this number.


### THE BUILD IS AN ALL-OR-NOTHING REBUILD, and the build guard interacts with running cases. Adding a virtual to a shared header is an ABI change that presents as a startup SEGV looking like a physics bug; `./build-all.sh` to BUILD-COMPLETE is the only safe rebuild. `libplasmaNewtonSolverPETSc` needs BOTH bashrcs — without the project's own `etc/bashrc`, PETSC_DIR is unset, the dlopen fails SILENTLY and SNES reports 'Registered types: 0()'.
**Evidence:** docs/CAPABILITIES.md:178-199 (§3b); memory rebuild-with-allwmake-not-piecemeal, openfoam-abi-partial-rebuild, allwmake-can-silently-skip-a-changed-directory

**Rule:** PROGRESS.md's 'how to resume' block: kill or finish running arms (the build guard refuses while they run), `./build-all.sh` → BUILD-COMPLETE, then verify the handover banner in the first Newton arm.

**Cost:** Recorded: a silent dlopen failure went unnoticed in a reference run; and 21,174 steps once ran Picard while believed to be Newton (commit 68688f9).


## Traps
- `git status` (no flags) takes 181.97 s in /home/kkourtza/soplasma-scratch and will blow past a 120 s tool timeout, looking like a hung command. DETECTOR/FIX: always `git status --porcelain -uno` (0.58 s) or scope to paths (0.00 s).
- 1,310 dirty lines of which 1,109 are untracked field files mean real content is invisible: 14 modified COMPARE.md and 26 modified configuration/ files are buried. DETECTOR: `git status --porcelain -uno | grep -E 'COMPARE|configuration'`.
- `git add -A` or a broad `git add tutorials/plasma` in this repo pulls in case RUN OUTPUT — measured: 436 paths of time directories/logs/postProcessing alongside 10 real dictionary changes. DETECTOR: always `git diff --cached --name-only` after staging and unstage what does not belong.
- An ignore rule using bare `[0-9]*` for time directories also matches `0.orig` and silently untracks the initial-field TEMPLATES. It has already destroyed 0.orig twice here, surfacing much later as 'cannot find file .../0/ePotential'. DETECTOR: use the paired `[0-9]*.[0-9]*/` + `[0-9]*e-[0-9]*/` + explicit `0/` forms the smoke/ block already uses.
- Reading a Newton arm's failure COUNT without the reason breakdown. 10.6% vs 1.0% looked like a solver quality difference and was a 100-iteration budget expiring 4-15 iterations short. DETECTOR: `grep -o 'due to [A-Z_]*' <log> | sort | uniq -c | sort -rn` — if DIVERGED_LINEAR_SOLVE dominates, the cap is the story.
- A CAPPED distribution cannot be read for its tail. The KSP histogram looked 'bimodal — a preconditioner cliff' purely because the cap truncated it at 100; the true one is median 4 / p90 14 / max 115, unremarkable. DETECTOR: check whether the max equals the cap exactly.
- Quoting the 449k COLD head-to-head ladder as a JFNK verdict. By construction it reaches only t=2e-10 with peak n_e at the 1.3e13 seed, dt/tau ~ 2e-7 — seven orders below the stiffness regime, where Picard needs only 2-4 correctors and Aitken is INACTIVE. DETECTOR: check peak n_e and whether Aitken relaxation ever engaged before believing any ladder.
- Assuming an arm ran Newton. The handover gates on every species being off its density floor and LATCHES; a cold-started arm at a dt where Picard's first step dies never hands over at all. DETECTOR: `grep 'Picard warm-up COMPLETE' <log>` — its absence once meant 21,174 steps silently ran Picard.
- BASELINE CONTAMINATION: four of seven ladder arms ran against a soPlasmaFoam binary replaced at 15:48 and were compared against 16:2x arms. DETECTOR: `ls -l` the solver binary's mtime against each arm's log mtime before comparing arms.
- The three `lad_*_k1000` arms were cut off MID-RUN (7, 3 and 1 steps; logs still growing at 18:47 with no process alive now) and have no COMPARE.md. DETECTOR: compare the log mtime against the last `Time = ` line and against `ps`; a log whose last line is mid-iteration is truncated, not a result.
- `smoke/STATUS` and `smoke/STATUS_DT` read 'exit=1' for all five beds and are dated 2026-08-12 — they are stale artefacts, not current health. DETECTOR: `ls -l smoke/STATUS` before reading it.
- Memory files sort misleadingly: session-state-2026-09-07c is NEWER than -09-07b, and -09-08 continues -09-07c, but the newest state of the CURRENT thread is in two non-session-state files written 2026-09-11. DETECTOR: `ls -t .../memory/*.md | head -6`, never the filename dates.
- `docs/framework-state.md` (SoEEDF) calls itself 'the current snapshot' but is dated 2026-09-04 and predates the entire Newton/JFNK solver, the flux-scheme work and the ksp-cap finding. DETECTOR: check its header date against `git log -1` in the live tree.
- A dielectric case run under `outerSolver newton` does NOT refuse — it silently drops the surface charge, because the guard is per-patch work nobody has done. DETECTOR: if `constant/regionProperties` declares a dielectric region or any role is a charging surface, do not trust a Newton run's surface charge.
- Fields populated inside `plasmaTransport::solve()` go STALE under Newton, silently, because Newton REPLACES that solve. Two were caught in a row: `convectiveFlux_` (left `limitSpeciesCo` protecting nothing — `Co_conv (e)` read exactly 0 for 1300+ steps) and `particleFlux_` (which `plasmaExternalCircuit` regulates on, so the current-driven electrode was inoperative for ~3900 steps). DETECTOR: if a field is filled inside solve(), assume Newton never updates it; check it reads non-zero.
- PETSc swallows OpenFOAM's stack trace (`PetscInitialize` installs its own handler; `-no_signal_handler` is read too late). DETECTOR/FIX: `reconstructPar`, then run serially under `gdb --batch -ex run -ex 'bt 40'`.

## Open questions
- Were the three `lad_*_k1000` arms killed deliberately or did the session simply end? None has a COMPARE.md, and `lad_newton_dt1e-10_k1000` died at step 1 with MPI_ABORT errorcode 59 — cause unknown and not recorded anywhere.
- Should `validation/constant/plasmaTables` (1,650 of the 3,228 tracked validation files) stay tracked? deferred-action-items already flags the same ship-vs-sweep inconsistency in the tutorials (`_LMEA_fast`/`_LMEA_minimal` ship 142 files each with `generateTables no`; `_fixedMesh`/`_AMR` sweep at runtime) and notes that shipped tables are exactly what went stale.
- Is `ThirdParty/petsc-3.24.0/` (untracked) meant to be vendored, ignored, or a local build? It is the single largest untracked path and no .gitignore rule covers it.
- Are the 197 unpushed SoPlasma commits / 52 unpushed SoEEDF commits deliberate (fork, no push rights, or a policy), or simply unpushed? `remotes/upstream/main` exists alongside `origin`, so the push target is not obvious from the tree alone.
- Has `needle_mrgate` been re-run with `kspMaxIt 1000`? No arm on disk postdates the 17:49 commit for that case, so the multi-region Newton verdict is currently based on a run made with the known-bad cap.
- Who ran `pcfix_maxit_dt1e-9` / `_dt1e-8`, and was the intent the ceiling re-test the memory calls for? The DIVERGED_MAX_IT 125 count against SNES `maxIt 50` suggests the nonlinear iteration is now the binding limit, but nothing on disk states the question the arms were asked.
- The owed item 'LFA frozen D is 672.2 and should be 264.6 m^2/s' appears in three consecutive session states but neither number occurs anywhere under docs/, tools/ or src/ — where does that value actually live (a case dictionary? a plot script?), and is it still wrong?
- '2450 Td as a cathode-fall reference is UNSOURCED' — it still appears 6 times in docs/design/electron-energy-balance.md and once in steady-and-stability-design.md:466. Whether the two COMPARE.md files and the plot script named in session-state-2026-09-07c still carry it was not verified (the greps that would check it are slow in this tree).
- Does the `verification/` tree (fluxScheme1D with its 504-run sweep and REPORT.md, fluxScheme2Dnonortho, testSnesJFNK, testSnesJFNK2Field) belong in git? It is entirely untracked and contains the only record of the CFS/ROUND scheme ranking.
- Neither tree has `.claude/skills` or `.claude/agents` yet — only `.claude/settings.local.json` with three `Bash(kill…)` permissions — so nothing constrains where the new environment's files should live, or whether both trees get one.
