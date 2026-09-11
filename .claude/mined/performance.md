# performance

## Summary
Performance, cost models and profiling for SoPlasma. The tree ships one profiler, `src/profilers/plasmaSimulationProfiler` — always compiled in, no dictionary switch, cumulative over the whole run, printed only at normal `End`. Measured here from three existing sweep logs: under `outerSolver newton` it accounts for only 2.4-3.7% of the solver's CPU time, because every one of its `start`/`stop` pairs sits in the Picard branch or in `singleRegionPoisson`/`plasmaTransport` methods that JFNK never calls. The real performance record lives in `docs/CAPABILITIES.md` sections 4/4b/4c and in ~95 memory files, not in the profiler. The governing numbers are a JFNK arithmetic cost model (1 Krylov iteration = 1 full residual assembly = ~1 Picard step; 10,000 evals x 4.16 s = 11.6 h/step at 1.15M against 10.06 h observed), a 10x per-step Newton penalty at 449k (142 vs 14.6 s/step), and a hard 1.15M convergence wall proved to be mesh size. Rule A2 (absorbing old rule 43) makes every numerics/scaling claim require two problem sizes and mesh-independence of the ITERATION COUNT, not wall clock — and this tree contains the instance that earned it (refining at fixed dt raised Co_conv 6.50 -> 20.62, exactly the 3.17x refinement, inverting the conclusion). Five optimisations have been measured and refuted at two sizes each; re-proposing them is the single most expensive mistake a performance agent can make.

## Facts

### The profiler instruments almost nothing under Newton: it accounts for 2.4-3.7% of the solver's CPU time. Every `plasmaSimulationProfiler::start/stop` pair at the solver's top level is inside the `else` (Picard) branch, and the rest are inside `singleRegionPoisson::solve()` / `plasmaTransport::solve()` — methods `snesNewtonSolver` never calls. `grep -rn plasmaSimulationProfiler src/numerics/newtonSolverPETSc/` returns ZERO hits.
**Evidence:** soPlasmaFoam.C:878-898 (all six top-level start/stop calls inside the `else` of the newton branch, opening at :876). Measured from the three mesh_sweep logs: mesh_2000 profiler total 2.34 CPU s vs ExecutionTime 69.52 s = 3.4%; mesh_20240 426.87 vs 11469.62 = 3.7%; mesh_200000 716.38 vs 29574.66 = 2.4%. The reports carry only 6 sub-sections (Calc E / Calc Emag / Calc phiE / Calc reducedE / Correct E BCs / chemistry ODE) — no `Build ePotentialEqn`, no `Solve ePotentialEqn`, no `buildEquations`, no `solveEquations`.

**Rule:** NEVER attribute Newton cost from the PROFILING REPORT. For a Newton run, profile from PETSc's own output (`-snes_monitor`, `-ksp_converged_reason`, `-log_view` via case `petscOptions`) and from `ExecutionTime` deltas per step. Treat the profiler as a Picard-path tool only.

**Cost:** Attributing 97% of a Newton run's cost to the 3% the profiler happens to see — the exact class of error that produced the refuted 'Schur complement is starved' hypothesis (-5% SNES for +14% wall).


### `ExecutionTime` is CPU time and `ClockTime` is wall, and their ratio here is 2.3x to 26x and VARIES BETWEEN ARMS of the same sweep. The runs are implicitly multi-threaded (hypre/BLAS OpenMP on a 32-core box) even when launched serially.
**Evidence:** Final lines of the three mesh_sweep logs: mesh_2000 `ExecutionTime = 69.52 s  ClockTime = 30 s` (2.3x); mesh_20240 `11469.62 s / 441 s` (26.0x); mesh_200000 `29574.66 s / 5109 s` (5.8x). `nproc` = 32. `/usr/bin/time -f "wall %e s maxRSS %M kB"` in validation/mesh_sweep.sh recorded wall 29.83 / 441.12 / 5111.73 s.

**Rule:** State WHICH clock every cost number uses. Compare arms in ExecutionTime (CPU) or pin threads (`OMP_NUM_THREADS=1`) before comparing wall. Never quote a wall-clock ratio between two arms without printing both ExecutionTime and ClockTime for each.

**Cost:** A 26x vs 2.3x thread-count difference between two arms of the same sweep makes any wall-clock 'scaling' number meaningless — and it is invisible unless both clocks are read.


### THE JFNK COST MODEL, and it is arithmetic not mystery: one JFNK Krylov iteration = one full nonlinear residual assembly = the cost of a whole Picard step. At 1.15M the budget 50 SNES x 200 KSP = 10,000 residual evaluations x 4.16 s = 11.6 h per timestep, against 10.06 h observed.
**Evidence:** docs/CAPABILITIES.md:335 and memory `newton-vs-picard-benchmark-state.md` (dated 2026-09-11). Bed: `positiveStreamer_fixedMesh`, 1.15M cells.

**Rule:** Before running any Newton arm, compute nSNES x nKSP x t_residual and compare with the budget you can afford. If the product exceeds the wall-clock budget, shrink the bed, not the tolerance.

**Cost:** Two full 1.15M runs (10 h CPU/rank warm, 3 h CPU/rank cold) that never completed step 1.


### Newton costs ~10x per step against Picard at 449k on the streamer bed: Newton 142 s/step, Picard 14.6 s/step. At the same time, warm-started Picard at 1.15M survives dt = 5e-11, so Newton must hold ~5e-10 merely to BREAK EVEN. Nothing measured to date shows Newton taking a larger step at all.
**Evidence:** docs/CAPABILITIES.md:327 (`scale_r4`, dt 1e-11, 2026-09-11); Picard dt ladder table (1e-11/2e-11/5e-11 all exit 0 warm-started from t=1e-9 on 1.15M; 1e-10 COLD = SIGFPE step 1 in `GaussSeidelSmoother::smooth`).

**Rule:** Quote speed as WALL CLOCK PER ns OF SIMULATED TIME AT EQUAL ACCURACY, never per step — a per-step comparison is rigged for Picard by the 4-5x/10x constant factor.

**Cost:** memory `deferred-action-items.md:139-141` states the rule explicitly because a per-step comparison would have declared JFNK dead on a metric it cannot win by construction.


### Picard's COLD ceiling at 449k is 2e-11, not 5e-11 — the 5e-11 figure is WARM-started 1.15M, a different bed and start. And Newton is NOT less robust in dt: at 5e-11 Picard SIGFPEs on step 1 while Newton completes 4 steps. Newton's problem is cost, not stability.
**Evidence:** memory `newton-vs-picard-benchmark-state.md` 'CORRECTED 2026-09-11b'; docs/CAPABILITIES.md head-to-head ladder table (dt 1e-11: Picard 20 steps 2-4 correctors / Newton 11 steps 0 failures; 5e-11: Picard SIGFPE 1 step / Newton 4 steps 3-of-48 failed).

**Rule:** When restating a break-even bar, name the BED and the START (cold vs warm) in the same sentence as the dt. A bar measured warm at 1.15M cannot be applied to a cold 449k arm.

**Cost:** A whole seven-arm ladder read as a negative JFNK verdict when it 'answers nothing' — it reached t=2e-10 with peak n_e still at the 1.3e13 seed, dt/tau ~ 2e-7, seven orders below the stiffness Newton exists to step over.


### `-ksp_max_it 100` was hardcoded and WAS the bug. Raising it to 1000 cut Newton failures 10.5x (10.6% -> 1.0%) for +1.7% residual evaluations and +2.5% wall clock. 38 solves hit EXACTLY 100 and were killed; with the cap raised, 39 solves needed 104 to 115 iterations and ALL converged — 0.92% of 4,235 solves, each 4 to 15 iterations short.
**Evidence:** memory `ksp-cap-was-the-bug.md` and docs/CAPABILITIES.md:405-440. Bed `grubert2009_pseudo` (2000 cells, dt=1e-10, 400 steps, 398 SNES solves, identical in every other respect). Arms: 42/398 @ 326.1 s vs 4/398 @ 334.2 s. Now a dict key: `snesNewtonSolver.C:1134 kspMaxIt_(dict.getOrDefault<label>("kspMaxIt", 100))`, used at :1642. Commits f17a489, aed5dbe.

**Rule:** Raise `kspMaxIt` to 1000 in the `newtonSolver` dict before quoting ANY Newton failure rate. The shipped default is still 100, so nothing moved silently — it must be raised deliberately.

**Cost:** Every negative JFNK verdict recorded on 2026-09-11 (the 449k ladder, the pseudo-transient dt=1e-11 ceiling, the needleDBD multi-region failure) traces to this one number and must be re-tested before being quoted.


### A capped distribution cannot be read for its tail. The 'bimodal preconditioner cliff' diagnosis was a measurement artefact of the cap itself: the true KSP distribution is median 4, p90 14, max 115 — unremarkable. A median of 4 Krylov iterations is the same order as Picard's 2-4 correctors, i.e. one PC application is already worth about one Picard sweep.
**Evidence:** memory `ksp-cap-was-the-bug.md` ('It looked bimodal because the cap had TRUNCATED it'); docs/CAPABILITIES.md:430.

**Rule:** Before reading any iteration-count histogram, check whether its maximum equals a configured cap. If max == cap, the distribution is censored and its tail is not data.

**Cost:** A full day of preconditioner hypotheses chasing an artefact.


### The physics-based PCSHELL (`assembledPmat false`) is the STRONGER preconditioner by every linear-algebra measure — median 3 vs 4, p90 6 vs 14, max 22 vs 115 (5x tighter tail), 8% FEWER residual evaluations (113,208 vs 123,234) — and a first reading of the same data said the opposite. Its remaining 17 failures are ALL line search (-6), a different defect.
**Evidence:** docs/CAPABILITIES.md reason-breakdown table (base / maxit / shell) and memory `ksp-cap-was-the-bug.md`. `assembledPmat_` is `snesNewtonSolver.C:1135`, default true; it was unreachable in production because `snesBridge` prefers PCFIELDSPLIT whenever a Pmat is supplied.

**Rule:** Rank preconditioners by the REASON BREAKDOWN (-3 linear solve / -5 maxIt / -6 line search) and by residual-evaluation count — the real JFNK cost — never by headline failure rate.

**Cost:** The best preconditioner was written off as 'solves a problem that was not broken'.


### RULE A2 (absorbs old rule 43): a numerics choice must state its SCALING and be validated at MORE THAN ONE problem size, and the test is mesh-independence of the ITERATION COUNT, not wall-clock. The recorded inverting instance: refining a mesh at fixed dt RAISES the Courant number, so '10x cells doubled the iterations' measured the CFL, not the mesh — holding Co fixed reversed the conclusion.
**Evidence:** SoEEDF/CLAUDE.md:200-210. The measurement is in validation/mesh_sweep.sh's header: `2000 cells (NX=400): Co_conv(e) = 6.50` vs `20240 cells (NX=1265): Co_conv(e) = 20.62`, ratio 3.17 — 'exactly the 3.16x refinement. Both are already far above 1.'

**Rule:** In any mesh study, scale dt with dx (dt/NX = const) and set `adjustTimeStep false`. Report median KSP iterations per step and SNES iterations per step at >= 2 sizes; wall clock is a secondary column.

**Cost:** The first design of mesh_sweep.sh used a single fixed dt on every mesh and 'the confound was not small — it was the whole effect'.


### The mesh-independence sweep's HEADLINE RESULT WAS NEVER WRITTEN DOWN: validation/mesh_sweep.log ends mid-way through the 200k arm, so the script's final analysis block never printed. Computed here from the existing logs (2026-09-11, read-only): median KSP iterations 5 (2k) -> 6 (20,240) -> 9 (200,000), i.e. ratio 1.8 for 100x cells; SNES iterations/step 4.0 -> 4.0 -> 6.3.
**Evidence:** Computed by awk over validation/mesh_{2000,20240,200000}/logs/log.run, discarding the first 50 relaxation steps as the script specifies. The 200k arm ran 79 steps (29 measured) against the other two's 70 (20 measured), which is a mild confound. Wall/maxRSS from the same logs: 29.83 s / 104,852 kB; 441.12 s / 296,452 kB; 5111.73 s / 2,059,916 kB.

**Rule:** Re-run `NSTEPS=20 NRELAX=50 validation/mesh_sweep.sh` to completion (or re-run only its final python block over the existing mesh_* logs) and record the table in docs/CAPABILITIES.md. Until then, quote 5/6/9 as INDICATIVE, not established.

**Cost:** The one measurement that would settle whether the FIELDSPLIT/Schur preconditioner scales is sitting unread in three log files.


### MEMORY MODEL, measured: ~10.3 kB per cell, serial. 2,000 cells -> 104,852 kB maxRSS; 200,000 cells -> 2,059,916 kB. Extrapolating, a serial 1.15M-cell run needs ~11.8 GB on a 30 GB box.
**Evidence:** `/usr/bin/time -f "wall %e s maxRSS %M kB"` lines in validation/mesh_{2000,200000}/logs/log.time. `free -g` reports 30 GB total.

**Rule:** Forecast maxRSS as 10.3 kB x nCells / nRanks before launching. A serial 1.15M Newton run is within RAM; two concurrent ones are not.

**Cost:** OOM on a 10 h run.


### THE 1.15M WALL IS MESH SIZE, not dt and not the initial state. Newton converges at 40k, 211k and 449k and does not at 1.15M — proved by running 1.15M BOTH warm (10 h CPU/rank, still step 1) and COLD (3 h CPU/rank, stuck right after handover at t=3e-12). Both at 99.9% CPU: computing, not deadlocked.
**Evidence:** docs/CAPABILITIES.md:320-334; memory `newton-vs-picard-benchmark-state.md` 'The wall, and what it is NOT'.

**Rule:** DO NOT re-run 1.15M under Newton to see if it works — answered twice, both ways. The open question is the PRECONDITIONER, not the case.

**Cost:** 13 CPU-hours per repetition, for an answer already in the file.


### THE BEDS, verified on disk with their exact cell counts. scale_r3 = 211,120 cells; scale_r4 = 449,413 cells; streamer_base and positiveStreamer_LMEA_fast and positiveStreamer_LMEA_minimal = 16,900 cells (bare 130x130 block, NO refinement); positiveStreamer_fixedMesh = the same block plus FIVE refineMesh levels = 1.15M cells at ~70 s per time step; grubert2009_pseudo = 2,000 cells; streamer-smoke (make-smoke-case.sh default) = 1,600 cells.
**Evidence:** `nCells:` notes in each constant/polyMesh/owner header. positiveStreamer_LMEA_fast/Allrun-serial:7 'No topoSet/refineMesh chain: this is the bare 130x130 block.' positiveStreamer_LMEA_fast/system/blockMeshDict blocks comment: 'The production streamer bed is this same 130x130 block plus FIVE refineMesh levels -- 1.15M cells, ~70 s per time step.'

**Rule:** Use grubert2009_pseudo (2k) or the 16.9k LMEA_fast bed for anything diagnostic; scale_r3/scale_r4 for scaling claims; positiveStreamer_fixedMesh ONLY for answers about the discharge itself. LMEA_fast at 96 um does not resolve the streamer head — its field magnitudes and propagation speed must never be quoted as physics.

**Cost:** docs/CAPABILITIES.md lists positiveStreamer_LMEA_fast as '~40k, ~2 s/step'; the checked-in tutorial mesh is 16,900 cells, so '~40k' refers to a make-smoke-case.sh derivative (NCELL=130 NREFINE=1), not the tutorial. Quoting either number for the other misstates the bed by 2.4x.


### A6, and the bench size is part of the rule: four Newton hypotheses were guessed and all four refuted on one day; the three INSTRUMENTS each paid off immediately — and all of it was found on the COARSE bed at ~2 s per run, after hours on the 1.15M bed at minutes per attempt. A NaN locator paid off by printing NOTHING through six DIVERGED_FNORM_NAN failures, proving the residual's components were finite and the NORM had overflowed.
**Evidence:** SoEEDF/CLAUDE.md:250-285 (rule A6, stated by the user 2026-09-11).

**Rule:** When something fails for an unknown reason, add the print/counter/locator and run it on the SMALLEST bed that reproduces the failure. A guess costs a full build-run-read cycle and usually returns nothing.

**Cost:** Four refuted hypotheses, one cycle each: ion immobility, sourceAwareScaling, the rate-table edge (table runs to 55,000 Td, state sat at 1,786), negative extrapolated rate coefficients.


### THE INNER SOLVES were at PETSc defaults until 2026-09-10 and were costing an order of magnitude. Setting `-fieldsplit_transport_ksp_rtol 1e-2` (not 1e-5), `-fieldsplit_transport_ksp_max_it 200` (not 10000) and `-fieldsplit_transport_ksp_gmres_restart 100` (not 30) gave: inner iterations/solve median 230 -> 53, mean 260 -> 59; worst single solve 7226 -> 166; solves hitting the cap 377 -> 2; total inner Krylov work 9.3x LESS — normalised per unit SIMULATED time, not per step.
**Evidence:** snesNewtonSolver.C:1712-1739 (the options string and its comment block).

**Rule:** Found only with `-ksp_view`: 'nothing in this solver's own output reports the inner KSP configuration'. Run one step with `petscOptions "-ksp_view"` before tuning any split.

**Cost:** 3.8 million wasted inner Krylov iterations across 377 capped solves.


### The A00 (phi) inner solve must be a FIXED LINEAR OPERATOR, and fixedness is a property of the ITERATION, not of exactness. The shipped default is `richardson` + `max_it 2` + `convergence_test skip` + `hypre` — deliberately not `preonly`+`lu` (affordable at 2000 rows, useless at 1e6: sparse LU is ~O(n^1.5) in 2-D). A convergence-TESTED inner solve makes the Schur operator non-linear in b and raised a SIGFPE.
**Evidence:** snesNewtonSolver.C:1692-1711. Iterative-inner measurements 55, 65, 71, 71, 82, 85 on successive applications. FGMRES (not GMRES) is mandatory for the same varying-PC reason: DIVERGED_BREAKDOWN at iteration 30 (GMRES's default restart) once dt grew to ~6.5e-10.

**Rule:** Any change to the fieldsplit inner solves must keep the iteration count FIXED (`convergence_test skip`) and the method linear in b — richardson or chebyshev, never gmres/cg.

**Cost:** SIGFPE in grubert_steady (docs/design/newton-ignition-experiments.md section 25c).


### REFUTED AT TWO SIZES — do not re-propose. ILU(1) on the transport split: 2k -33% steps, 20k +3.5%. `selfp` Schur approximation: 2k -22% steps, 20k KSP median 11 -> 12. Per-cell residual/state scaling: 2k -46% steps and SNES 1887->1224, 20k +1.4% mean dt with KSP median 13 vs 13 and SNES/step 5.8 vs 5.9. ILU(1)+selfp together are exactly baseline (they interfere).
**Evidence:** docs/CAPABILITIES.md:256-259 (per-cell scaling dated 2026-09-10). Beds exist on disk: validation/pc_ilu1, pc_ilu1_20k, pc_selfp, selfp_20k, pcs_20k, pcs_200k, adt20k_base, adt20k_ilu1.

**Rule:** Every one of these looked like a large win at 2,000 cells and vanished at 20,000. Any preconditioner claim measured only at 2k is an ARTEFACT until repeated at >= 20k.

**Cost:** Three separate 'wins' shipped as defaults would each have been a no-op or a regression on production meshes.


### REFUTED: more `-fieldsplit_phi_ksp_max_it` as the cure for the 1.15M wall. At 449k, 10 steps each at dt=1e-11: 2 cycles = 77 SNES / 387 Krylov / 121 s per step; 8 cycles = 73 / 373 / 138 s. Eight buys -5% SNES and -4% Krylov for +14% wall clock — a net LOSS, and far too small an effect to explain a total failure.
**Evidence:** docs/CAPABILITIES.md:271 and memory `newton-vs-picard-benchmark-state.md`. Tested at 449k deliberately: minutes per run instead of hours (rule A6).

**Rule:** The Schur complement is NOT starved of A_phi^-1. The cause of the 1.15M wall remains unknown; do not spend another arm on phi-cycle count.

**Cost:** A day chasing a -4% effect as the explanation for a total failure.


### REFUTED: `sourceAwareScaling` (scale a species' residual by its chemistry source instead of n/dt) — REMOVED 2026-09-11. `rms(chemP)/(rms(n)/dt)` measured per species is 0.833 for all seven source-dominated species and 1.6e-4 for the electron. The source NEVER exceeds ddt, so `max(ddt, src)` is a NO-OP. And REFUTED: `extrapolateGuess` as a cure for block imbalance — 38 vs 40 converged solves over 11 steps, no effect.
**Evidence:** docs/CAPABILITIES.md:265-266. A per-species scale print disproved sourceAwareScaling in ONE run (SoEEDF/CLAUDE.md:271).

**Rule:** Note the distinction the tree makes explicitly: sourceAwareScaling scaled per CELL within a field; per-FIELD block scaling (phi ~ 1e4 V vs n ~ 1e13-1e20, `MatCreateSNESMF` picks ONE h from ||u||) is a DIFFERENT and still-untested candidate.

**Cost:** Re-proposing a measured no-op, and conflating it with the untested idea that remains open.


### REFUTED: tightening the Poisson tolerance to cure the lateral asymmetry — GAMG at 1e-16 costs 2000 iterations (286x) for nothing, it improves phi/Ex but leaves Ey (the driving component) unchanged, and 1e-16 absolute on a 225 V field is below double precision anyway. Separately, `fvMatrix::residual()` standalone in parallel is BROKEN — misreports by ~21 orders of magnitude; always compute the outer residual by explicit `fvc::`.
**Evidence:** docs/CAPABILITIES.md:261 and :263; snesNewtonSolver.C:355-357 cites memory `fvmatrix-residual-broken-in-parallel-use-fvc-instead`.

**Rule:** Tolerance is never the lever it looks like: the order study measured linear 1e-13 / outer-gate 1e-11 as +34% wall clock (487 vs 363 s at dt=4e-13) for a 6-7 s.f. IDENTICAL answer. Baseline is linear 1e-10 / gate 1e-8.

**Cost:** 286x on a Poisson solve, or 34% on a whole study, for a bit-identical result.


### Counter-intuitive and measured: the FINE-dt arm is CHEAPER PER STEP — 0.21 s/step at dt=1e-13 vs 0.73 s/step at dt=4e-13 — because stiff-ODE work scales with step size. So the 2000-step arm is not the expensive one, and shortening `endTime` is the biggest lever, helping all arms about equally.
**Evidence:** memory `order-study-tolerance-policy.md` (2026-08-21). Also: `reportInterval` 25 not 1 — the LMEA per-term dump is Info-only and verified inert; at 1 it bloats logs and costs time for nothing.

**Rule:** Do not assume a smaller dt is proportionally more expensive. Measure s/step at each dt before sizing a sweep, and cut endTime rather than step count.

**Cost:** Sizing an order study around the wrong axis.


### THE BIGGEST RECORDED PERFORMANCE LEVER is the outer-coupling contraction, and we pay for it TWICE. Picard on ballast_clean needed 73,567 PIMPLE iterations over 3,306 steps = 22 correctors/step (tail to 141, cap 150); Aitken damped omega to 0.148; and the `coupling margin` dt governor is keyed on that SAME omega, clamping dt to 1/83 of what accuracy allowed. Compounded: ~1000x more work than an accuracy-limited solve with a healthy loop.
**Evidence:** memory `deferred-action-items.md:221-233`. Not the timestep, not the Courant caps, not the chemistry — and it would make ANY case slow.

**Rule:** When a case is slow, read the corrector count and the Aitken omega FIRST. If omega is damped and the dt limiter says `coupling margin`, the cost is the coupling, not the physics — and 80x dwarfs Newton's 4-5x per-step penalty.

**Cost:** Optimising the 4-5x while the 80x sits untouched.


### THE dt LIMITER NAMES ITSELF in the log, and counting them is the cheapest performance diagnostic in the tree. needleDBD, 322 steps, LMEA: 197 energy relaxation, 99 growth cap (1.2x), 19 rejection memory, 4 Co_chem, 2 coupling margin (backoff); `maxSpeciesCo 100` NEVER BOUND. adt20k_base: 853 `Co_conv (energy)`, 1 `growth cap (errorMaxGrow)`.
**Evidence:** `Info<< "    deltaT set by:            " << dtLimiterName_` at plasmaTimeControl.C:1996; names bound at :1509 ("energy relaxation"), :1530 ("Co_diff"), :1591 ("Co_chem"), :1644/:1649 ("coupling margin (backoff)/(hold)"), :1692-3 ("growth cap"), :1706 ("rejection memory"). memory `needledbd-first-run-lessons.md`.

**Rule:** Run `grep -o "deltaT set by: *.*" <log> | sort | uniq -c | sort -rn` before proposing any speed change. Optimise the limiter that actually binds; a cap that never binds is a rail, not a control.

**Cost:** Tuning maxSpeciesCo on a case where it never bound once in 322 steps.


### `floor hits` counts TIMESTEP-floor hits (minDeltaT), NOT density clipping. Misreading it produced a whole false hypothesis and a proposed research direction.
**Evidence:** `if (newDeltaT < minDeltaT_) { ++floorHits_; ... }` at plasmaTimeControl.C:1747 (comment '(3) ABSOLUTE FLOOR on deltaT'); printed at :1777, :2086, :2625. memory `floor-hits-is-the-timestep-floor.md` (2026-09-08: 'floor hits 930-1026' read as 100x more density clipping).

**Rule:** A large `floor hits` means the run sat pinned at minDeltaT grinding — a performance symptom, nothing more. Never infer a mechanism from a counter without reading the code that increments it.

**Cost:** A bounded-CFS research direction proposed on a misread counter, and the comparison was ALSO confounded (the two runs were at different states, 2.6x apart in time).


### DECOMPOSITION: every case uses `method scotch` with `numberOfSubdomains` supplied from `configuration/config` — 8 for scale_r3, scale_r4 and positiveStreamer_fixedMesh; 4 for positiveStreamer_LMEA_fast, streamer_base and grubert2009_pseudo. NEVER hand-edit system/decomposeParDict; it reads `$numberOfSubdomains`.
**Evidence:** `grep numberOfSubdomains <case>/configuration/config` on each bed; `method scotch;` in each system/decomposeParDict. A streamer run is '~12 minutes of 8-core compute' (memory `pre-run-checklist.md`).

**Rule:** Change rank count in `configuration/config` only, then prove the key is not dangling (docs/CAPABILITIES.md:87). Editing the generated side is the defect G1 exists to prevent.

**Cost:** A silently ignored edit and a run at the wrong rank count.


### TEST AT LEAST TWO DECOMPOSITIONS: the same collective-behind-a-local-guard bug DEADLOCKS under `simple` (whole electrode on one rank) and raises MPI_ERR_TRUNCATE under `scotch` (patch split). One cause, two symptoms; either alone would have misled. Separately, `reduce` downstream of `if (p.size()==0) return` ran on a SUBSET of ranks and returned garbage that read as a physics failure ((D/delta+uEff)/(D/delta) = -1e+300) — two hours went into the wall-flux closure before a per-rank trace showed the reduction itself was broken.
**Evidence:** memory `collective-behind-a-local-guard.md` (now rule 31, commit 7b47378); docs/design/parallel-blockers.md (decomposePar scotch, numberOfSubdomains 4, Open MPI 4.1.6).

**Rule:** Per-rank `Pout` probes bracketing each construction found it in ONE build after reading the collectives got nowhere three times. And: verify any parallel arm against a serial reference on the same case — serial, np2 and np4 agreed to every printed digit on n_e,max, Te,max, Te,min.

**Cost:** A full session (2026-09-07). 'A parallel run that merely RUNS is not evidence it is right' — a fix applied to one of two sibling models once left electrostatics dead for three weeks.


### THE PROFILER ITSELF IS A LATENT RULE-31 HAZARD. `report()` loops over each rank's OWN `cumulativeTimes_` map and calls the collective `reduce()` inside that loop, so if any rank's key set differs the collective count mismatches. Two sections are already behind runtime conditionals — `if (chem_)` and `if (finalIter && gasHeating_)`.
**Evidence:** plasmaSimulationProfiler.C:55-68 (`for (auto const& [key, time] : cumulativeTimes_) { ... reduce(maxT, maxOp<double>()); reduce(avgT, sumOp<double>()); }`); plasmaTransport.C:929-931 and :1576-1579.

**Rule:** Any new profiled section must be entered by EVERY rank or not at all. Never put `plasmaSimulationProfiler::start/stop` behind a guard that depends on a local patch/cell count, a region membership, or a rank-local branch.

**Cost:** A deadlock or MPI_ERR_TRUNCATE at the very END of a long run, after all the compute is done — the most expensive possible place to hang.


### The profiler is ALWAYS ON (no dictionary switch anywhere: `grep -rn profil` over docs/reference and configuration/config returns nothing), is CUMULATIVE over the whole run (`reset()` is defined at plasmaSimulationProfiler.C:102 and never called), and `report()` runs only at normal termination (soPlasmaFoam.C:1051, after the time loop). A killed, timed-out, crashed or still-running case produces NO profiling report.
**Evidence:** Confirmed: validation/grubert2009_ballast_clean/logs/log.soPlasmaFoam is 221 MB and contains no PROFILING REPORT at its tail — the run never ended normally. The library is `libplasmaProfilers` (src/profilers/Make/files), built at position 25 in build-all.sh.

**Rule:** To get a profile out of a long case, set a short `endTime` so the run REACHES `End`. Do not expect a profile from `run-guarded.sh`'s `timeout` path or from a `pkill`ed run.

**Cost:** Waiting 10 hours for a profile that structurally cannot appear.


### DISK: a runaway warning loop once wrote 9.7 GB / 91M lines before ONE timestep finished and took the machine down. Three tools exist because of it: `run-guarded.sh [maxLogMB] [maxSeconds]` (default 50 MB / 600 s, pipes through `head -c` so the solver dies of SIGPIPE), `run-long.sh` (200 MB cap), and `watchdog.sh <log> [maxMB]` (polls every 3 s, `pkill -9 -x soPlasmaFoam` on overrun).
**Evidence:** docs/CAPABILITIES.md:126; run-guarded.sh header and its `( timeout "$MAXSEC" soPlasmaFoam 2>&1 | head -c $((MAXMB*1024*1024)) > log.solver )`; watchdog.sh.

**Rule:** Nothing unattended runs without a log cap. Use run-guarded.sh/run-long.sh, or pair a bare launch with watchdog.sh.

**Cost:** The machine went down; 9.7 GB written for zero completed timesteps.


### THE WSL DISK TRAP: `df -h /home/kkourtza` is MISLEADING — it reported 713 GB free while the real constraint, Windows C:, was at 98% (13 GB left), and WSL terminates hard when it cannot extend the .vhdx. Deleting inside WSL returns ZERO space to C: until an explicit compact. Current state: C: 465G total, 246G free (48%).
**Evidence:** memory `wsl-disk-space-trap.md` (WSL died mid-session 2026-08-31). Verified now: `df -h /mnt/c` -> `C:\ 465G 219G 246G 48% /mnt/c`; `/dev/sdd 1007G 192G 765G 21% /`.

**Rule:** The pre-run storage check is `df -h /mnt/c`, NOT the WSL path. Forecast the case's WRITE VOLUME against /mnt/c free space. To reclaim: `wsl --manage <distro> --set-sparse true`, or `wsl --shutdown` + `Optimize-VHD -Mode Full`.

**Cost:** A hard WSL termination mid-session, losing every running case.


### SNAPSHOT POLICY (rule B1): `writeInterval <= endTime/10` — at least ten snapshots, set when the case is configured; restart from a snapshot rather than from t=0 whenever it cannot change the conclusion. A SIGFPE diagnosis took 14 ns of simulated time from the 1.940229e-06 snapshot instead of 1.94 us from zero — four orders of magnitude cheaper.
**Evidence:** SoEEDF/CLAUDE.md rule B1; memory `restart-from-snapshots-and-preserve-them.md` (user, 2026-09-10).

**Rule:** A run REWRITES the snapshot it restarts from, `uniform/time`'s deltaT included — measured first steps of 7.71e-12 vs 1.11e-11, each exactly 2x its own stored deltaT. If a restart point will be reused, `cp -a` it aside and `chmod -R a-w` it FIRST. The FIELDS stay byte-identical through that rewrite, so the corruption is invisible unless the time state is checked specifically.

**Cost:** Two 'identical' restarts that diverge at step one, with no visible cause.


### DISK SHARING: sweep cases use `cp -al` so the 320 MB mesh is shared, and `du -sm` per directory then over-counts (217 GB reported against 135 GB real). But `constant/plasmaTables` is WRITTEN at start-up — two fresh cases launched together race on one set of inodes and the loser dies on 'Mechanism hash mismatch'. And `sed -i` is safe on a hardlink while Python `open(w)` and shell `>` truncate IN PLACE, corrupting every sibling.
**Evidence:** memory `hardlinked-case-trees.md`. Cost on 2026-08-18: co0.25/co1.0/co1.5/co2.5 AND ~/streamer-semi-heat all silently acquired `solver adaptiveError` + `chemChangeFloor 1e15`; `rm -rf` on the scrapped case did not undo it.

**Rule:** Share what the run READS (polyMesh only), never what it WRITES. After `cp -al`, give every case a private plasmaTables and verify with `stat -c '%h %n' <case>/constant/plasmaTables/DLN_vs_meanE` — must print 1. Archive a cp -al family into ONE tar (tar preserves hardlinks only within one archive; separate tars inflated the archive ~60 GB).

**Cost:** Silent cross-contamination of an entire sweep, and ~60 GB of archive bloat.


### A HUNG RUN AND A SLOW RUN LOOK IDENTICAL. Both show the same step count and both sit at 99.9% CPU, because OpenMPI busy-polls inside PMPI_Waitall — deadlocked ranks show state `Rl` (running), not D or S. A 28-minute hang was reported as 'slow' and found only by comparing the log's MTIME against the wall clock.
**Evidence:** memory `watching-a-run-log.md` (user, 2026-09-01: 'check running cases regularly ... each 5-10 minutes???!!'). The stack came out in one command: `gdb -p <pid> -batch -ex "bt 8"` identified `gradScheme::grad()` -> `UPstream::waitRequests` -> `PMPI_Waitall`, i.e. a halo-exchange mismatch.

**Rule:** Every wait on a run gets a STALENESS check, not just a progress check: `[ $(( $(date +%s) - $(stat -c %Y "$LOG") )) -gt 600 ]`. Poll every 60 s, report every 5-10 min. On a suspected hang, get the stack BEFORE killing.

**Cost:** 28 minutes of a hung run reported as progress, plus a second wait set up to do the same thing.


### MULTI-REGION UNDER NEWTON IS DONE (commit 10b7b2c, 2026-09-11) and the memory that says otherwise is STALE. `outerSolver newton` now accepts `multiRegionPoisson`; the dielectric and far-field regions are packed into the RAGGED TAIL of the DOF layout and no `lduPrimitiveMeshAssembly` is involved — the monolithic assembly is a LINEAR-SOLVER ACCELERATION, not a correctness requirement.
**Evidence:** docs/CAPABILITIES.md:549-570 (gate table: single-region bit-identical PASS, 62 fields identical; `Poisson tail: 1 extra region(s), 5976 rows` matching cellOffsets 0/71745/77721; first linear solve 2 iterations, ||F|| 1037 -> 245). Perturbation measurement on needleDBD: interface cells (160) d(res) rms 1.549 vs interior (71585) 0.00354, ~440x. CONTRADICTED by memory `newton-vs-picard-benchmark-state.md` 'DO NOT RE-RUN ... needleDBD under Newton -- it refuses multiRegionPoisson'.

**Rule:** CAPABILITIES.md is newer than the memory here. Under rule A3/20, re-read the source and reconcile before quoting either — and note the trap the new path creates: `correctBoundaryConditions()` must be called on EVERY region before each trial-state residual evaluation, or the mapped refValue/refGrad reflect the last PICARD state.

**Cost:** Refusing to run a benchmark arm that now works, or running one whose interface data is silently stale.


### GRUBERT IS NOT A VIABLE BENCHMARK VEHICLE and every time-marched route to its DC-glow operating point is closed. ~10 ballast arms (seriesResistor 1e8, seriesRC 1e6/1e8/1e9/5e9) all died at 150/150 correctors with dt collapsing ~25,000x. Six `iset` current-control arms stalled at t ~ 9.5e-07, one after 80,983 steps with 23,575 discards. `ddtSchemes steadyState` is REFUSED by the solver.
**Evidence:** docs/CAPABILITIES.md section 4b (settled 2026-09-11). The one route still open — pseudo-transient (BDF2 + `adjustTimeStep false` + fixed dt + currentSource, NO relaxationFactors) — measured 9472 consecutive converged steps, 0 failures, 4 correctors/step against a 150 cap, at ~5 h per us of simulated time; only ever run to 9.5 ns.

**Rule:** Benchmark JFNK on the STREAMER, not on Grubert. And never add `relaxationFactors` by hand to the pseudo-transient recipe: it takes 9472/9472 converged to 0/10 — SoPLASMA already runs adaptive Aitken outer relaxation (omega ~0.371) and fixed factors fight it.

**Cost:** Section 4b exists 'so the effort is not restarted from the top'; the ballast route alone cost a full day before Almeida 2017 explained why dR/dI = 0 makes the load line tangent.


### BDF2 is 5.5x faster than Euler, measured, and is part of the pseudo-transient recipe (`ddtSchemes backward`).
**Evidence:** memory `steady-mode-ddt-is-the-diagonal.md:31` ('BDF2; measured 5.5x faster than Euler'), with dt steady at 1e-12 and 4 correctors/step against a 150 cap.

**Rule:** Check `ddtSchemes` before proposing any other speed change; Euler on a case that could take BDF2 is a free 5.5x left on the table.

**Cost:** 5.5x.


### THE UNIT BED IS OFTEN THE STRONGER TEST, NOT MERELY THE FASTER ONE. A ~50 minute streamer run (70 s/step) was staged to ask whether the contraction ratio rho could replace the Aitken omega as the deltaT governor; moving it into `testAitken` answered it in ONE SECOND — and better, because the unit bed has GROUND TRUTH (did the step converge in budget?) where the CFD had only two diagnostics agreeing. It also surfaced a finding the CFD could not: rho >= 1 is the WRONG threshold, since the loop stops converging in budget at rho ~ 0.38.
**Evidence:** memory `shortest-sufficient-test.md` (2026-08-30). Unit beds in the tree: testPoissonSymmetry, testSnesJFNK, testSnesJFNK2Field, testFluxScheme, testWallFlux, testWallLoss, testDischargeCurrent, testEmission, testAitken, testCoulombHeating, testVibRelax, plasmaChemistry0D.

**Rule:** Target UNDER FIVE MINUTES. Before reaching for a CFD case, ask: is there a unit test that has ground truth the CFD case does not? Validate acceleration work (predictor, Anderson, relaxation governors) in testAitken FIRST.

**Cost:** 50 minutes and a weaker answer, against 1 second and a stronger one.


### A one-line fvSolution change was worth 2.75x: adding an `"n_.*(Final)?"` entry took 2000-cap solves from 169 to 0.
**Evidence:** memory `fvsolution-regenerated-and-ion-solver.md:39`.

**Rule:** Check for solves hitting their iteration cap (`grep -c "solution singularity\|2000 iterations"`) before assuming the cost is physics. A capped linear solve is wasted work, not hard work.

**Cost:** 2.75x on every step.


### needleDBD costs ~2.8 s/step and its shipped `endTime 2e-08` is ~3 HOURS — a first run reached only t = 1.71e-9 in 322 steps. The config's own comment says the interesting window opens at ~2.1 ns.
**Evidence:** memory `needledbd-first-run-lessons.md` (2026-09-02).

**Rule:** State the discriminating observable and the time it becomes visible BEFORE setting endTime. A 3-4 ns endTime crosses the observable and finishes in well under an hour.

**Cost:** 3 hours for a window that closed at 2.1 ns.


### An uncommitted fix in the working tree matters for benchmark arms: needleDBD's `Allrun-serial` hardcoded `SoPLASMA_TOOLS=../../../../tools`, which resolves only at that file's original depth — a COPY of the case (a benchmark arm) died at msh2Dto3D with 'No such file or directory'.
**Evidence:** `git diff -- tutorials/plasma/soPlasmaFoam/needleDBD/Allrun-serial` in /home/kkourtza/soplasma-scratch (modified, uncommitted): now `SoPLASMA_TOOLS=${SoPLASMA:+$SoPLASMA/tools}` with the relative path as fallback.

**Rule:** When copying a tutorial case out of the tree to make a benchmark arm, export `SoPLASMA` first, or the setup chain dies on a relative tools path.

**Cost:** A benchmark arm that fails at mesh generation for a reason unrelated to the experiment.


### BUILD AND RUN HYGIENE THAT COSTS WALL CLOCK: `wmake` overwrites the .so a running solver has mmap'd, and `build-all.sh` first runs `check-no-running-solvers.sh` and REFUSES to build while a solver is running — exiting 0 with no BUILD-COMPLETE. A stalled case from twenty minutes earlier silently blocked a rebuild, and the symbols found afterwards were from an earlier manual wmake.
**Evidence:** memory `kill-runs-before-launching.md` (user, 2026-09-05) and `pre-run-checklist.md`. Traps: `pkill -f soPlasmaFoam` does NOT match (the processes are `soPlasmaFoam -parallel`) AND can kill the calling shell (its own command line contains the string); `pgrep -x` fails on names longer than 15 chars (plasmaChemistry0D); `pgrep -c` prints 0 AND exits non-zero, so `|| echo 0` yields "0\n0" and breaks every numeric test.

**Rule:** Kill by NAME as its own simple command: `pkill -9 soPlasmaFoam; pkill -9 mpirun`, then `pgrep -c soPlasmaFoam` before building. Never `pkill -f`.

**Cost:** Measured 2026-08-20: a `pkill -9 -f "run-case.sh"` killed itself between launching two arms, so one arm of an A/B silently never started.


### REUSE THE MESH. The streamer Allrun scripts start with Allclean then blockMesh + 5x(topoSet+refineMesh) + extrudeMesh — many minutes on the 1.15M case, reproducing a bit-identical mesh, on every restart.
**Evidence:** memory `reuse-the-mesh.md`; positiveStreamer_fixedMesh/Allrun-parallel (the full chain). `tools/make_mesh_arm.sh` reuses the mesh-INDEPENDENT Boltzmann/ion tables deliberately; mesh_sweep.sh copies constant/plasmaTables rather than regenerating ('genMechTables is slow and would add nothing').

**Rule:** Restart with a Rerun-nomesh pattern: delete reconstructed and decomposed time dirs (keeping 0.orig and constant/), `cp -r 0.orig/* 0/`, plasmaCreateSpeciesFields, changeDictionary, seed, `decomposePar -force`, run. Deleting stale time dirs is not optional — controlDict has `startFrom latestTime` and a half-written time kills the run.

**Cost:** Many minutes per debugging iteration, dominating turnaround.


## Traps
- THE PROFILER LOOKS LIKE COVERAGE AND IS NOT. Under `outerSolver newton` it reports 6 sub-sections totalling 2.4-3.7% of CPU and prints a clean-looking table, giving no sign that 96%+ is unmeasured. DETECTOR: sum the Max[s] column and divide by the run's final `ExecutionTime`. If it is under ~50%, the profile is not an accounting. One command: `grep -A12 'PROFILING REPORT' <log>` next to `grep ExecutionTime <log> | tail -1`.
- WALL-CLOCK COMPARISONS ACROSS ARMS ARE SILENTLY INVALID because the CPU/wall ratio varies 2.3x-26x between arms of the same sweep (hidden OpenMP threading on 32 cores). DETECTOR: `grep ExecutionTime <log> | tail -1` on every arm and compute ExecutionTime/ClockTime. If the ratios differ, compare CPU seconds or set OMP_NUM_THREADS=1.
- A CAPPED ITERATION HISTOGRAM READS AS A BIMODAL FAILURE MODE. DETECTOR: if the observed maximum equals a configured cap (`-ksp_max_it`, `-fieldsplit_*_ksp_max_it`, nOuterCorrectors), the distribution is censored. `grep -o 'iterations [0-9]*' <log> | sort -n | tail -1` and compare against the cap.
- A PRECONDITIONER WIN AT 2,000 CELLS IS AN ARTEFACT. ILU(1), selfp and per-cell scaling all showed -22% to -46% steps at 2k and evaporated or reversed at 20k. DETECTOR: refuse to record any numerics verdict measured at one problem size (rule A2). The 20k beds already exist: validation/pc_ilu1_20k, selfp_20k, pcs_20k, adt20k_base.
- REFINING A MESH AT FIXED dt RAISES THE COURANT NUMBER AND INVERTS SCALING CONCLUSIONS. Measured: Co_conv(e) 6.50 at 2000 cells -> 20.62 at 20,240, ratio 3.17 = exactly the 3.16x refinement. DETECTOR: print Co on every arm, or scale dt with dx (dt/NX = const) and set `adjustTimeStep false`, as validation/mesh_sweep.sh does.
- A COLD-STARTED NEWTON LADDER CANNOT TEST NEWTON AT ALL. Handover needs one completed Picard step, so at any dt where Picard's first step dies the Newton arm dies inside it (measured at dt=1e-10: both arms, same SIGFPE, 1 step, handover count 0). And a cold arm reaching t=2e-10 sits at dt/tau ~ 2e-7, seven orders below the stiffness regime. DETECTOR: `grep -c 'outerSolver newton (SNES)' <log>` — 0 means the run measured Picard and the arm is VACUOUS. mesh_sweep.sh prints exactly this check.
- A HUNG PARALLEL RUN LOOKS PERFECTLY HEALTHY: 99.9% CPU on every rank, state Rl not D/S, because OpenMPI busy-polls in PMPI_Waitall. DETECTOR: log MTIME staleness, `echo "idle $(( $(date +%s) - $(stat -c %Y $LOG) )) s"`, not step count and not top.
- A STALE LOG READS AS AN INSTANT FAILURE. `Rerun-nomesh` deletes log.soPlasmaFoam but the solver does not create the new one until preprocessing (fields, seed, decomposePar) finishes — a minute or more at 1.15M — so an `until grep ...` loop started right after launch matches the PREVIOUS run's FOAM FATAL ERROR. DETECTOR: `rm -f "$d/log.soPlasmaFoam"` yourself, then `until [ -s "$d/log.soPlasmaFoam" ]; do sleep 5; done` before grepping.
- `pgrep -c soPlasmaFoam == 0` DOES NOT MEAN THE CASE IS DEAD — preprocessing runs for minutes with no solver process, longer when EEDF tables rebuild. Concluding death and relaunching produced TWO mpirun instances writing the same processor* dirs. DETECTOR: check the WRAPPER (`pgrep -f Rerun-nomesh`) or treat alive as `pgrep -c -f 'soPlasmaFoam -parallel'` OR `pgrep -c -x mpirun`.
- `df -h` INSIDE WSL IS THE WRONG FILESYSTEM: it reported 713 GB free while Windows C: had 13 GB and WSL then died hard. DETECTOR: `df -h /mnt/c`. And deleting inside WSL returns nothing to C: until an explicit vhdx compact.
- A COLLECTIVE BEHIND A RANK-LOCAL GUARD FAILS DIFFERENTLY UNDER DIFFERENT DECOMPOSITIONS — deadlock under `simple`, MPI_ERR_TRUNCATE under `scotch` — and a reduce after an early return returns garbage that reads as a physics failure (-1e+300). DETECTOR: per-rank `Pout` probes bracketing each construction; test at least two decomposition methods. The profiler's own report() has this shape (reduce inside a loop over a per-rank map).
- A PROFILE THAT NEVER APPEARS. `plasmaSimulationProfiler::report()` runs only after the time loop (soPlasmaFoam.C:1051), so every killed, timed-out or crashed run yields nothing — confirmed on a 221 MB ballast_clean log. DETECTOR: `tail -c 4000 <log> | grep -c 'PROFILING REPORT'`; if 0, the run did not reach End.
- SHELL `>` AND PYTHON `open(w)` WRITE THROUGH A HARDLINK in a `cp -al` sweep tree, corrupting every sibling case and the source; `sed -i` is safe. Also `constant/plasmaTables` is WRITTEN at start-up, so two fresh cases launched together race and the loser dies on 'Mechanism hash mismatch'. DETECTOR: `stat -c '%h %n' <case>/constant/plasmaTables/DLN_vs_meanE` must print 1.
- A RESTART REWRITES THE SNAPSHOT IT STARTS FROM, including uniform/time's deltaT, so a second restart from 'the same' snapshot begins at a different dt and diverges at step one — while the FIELDS stay byte-identical, making it invisible. DETECTOR: compare `uniform/time` before and after, not n_e. Prevention: `cp -a` the restart point aside and `chmod -R a-w` it.
- `floor hits` IS THE TIMESTEP FLOOR, NOT DENSITY CLIPPING (plasmaTimeControl.C:1747). A high value is a run pinned at minDeltaT grinding. DETECTOR: read the incrementing line before inferring any mechanism from any counter.
- BASELINE CONTAMINATION: four of seven arms in the 449k ladder ran against a soPlasmaFoam replaced at 15:48 and were compared against 16:2x arms. DETECTOR: `ls -la --time-style=full-iso $FOAM_USER_APPBIN/soPlasmaFoam $FOAM_USER_LIBBIN/*.so` and compare against each arm's log start time; re-run every arm on one binary before reading any ladder.
- BUILD-COMPLETE IS NOT PROOF OF A BUILD: build-all.sh printed it with zero `error:` lines while four components had FAIL and libplasmaTools.so did not exist. And it exits 0 with NO BUILD-COMPLETE when a solver is still running. DETECTOR: check all three — `grep -c FAIL`, `grep -c BUILD-COMPLETE`, and `ls -la` the .so with a fresh timestamp.

## Open questions
- The 20,240-cell arm of the mesh sweep is anomalous and nobody has explained it: ExecutionTime 11,469.62 CPU s against ClockTime 441 s (26x — i.e. ~26 threads busy), while the 2,000-cell arm ran at 2.3x and the 200,000-cell arm at 5.8x. Its CPU cost per step (163.9 s) is 165x the 2k arm's (0.99 s) for only 10.1x the cells, whereas 200k costs only 2.3x more per step than 20,240 for another 9.9x cells. Something other than problem size changed between those arms (hypre OpenMP thread count is the first suspect). Until this is resolved, no per-step scaling law should be quoted from that sweep.
- docs/CAPABILITIES.md:325 lists positiveStreamer_LMEA_fast as '~40k cells, ~2 s/step, THE debugging bed', but the checked-in tutorial mesh is 16,900 cells (bare 130x130, Allrun-serial:7 says 'No topoSet/refineMesh chain'). Which artefact is the '~40k, ~2 s/step' bed — a make-smoke-case.sh derivative with NCELL=130 NREFINE=1, or something else? The ~2 s/step figure is not reproducible from anything in the tree as read.
- The mesh_sweep.sh final analysis block NEVER RAN: validation/mesh_sweep.log ends part-way through the 200,000-cell arm, so the project's own mesh-independence table was never printed or recorded. The numbers I computed here (median KSP 5/6/9, SNES/step 4.0/4.0/6.3) are fresh and unreviewed, and the 200k arm's measured window differs (29 steps vs 20) because it ran 79 steps rather than 70.
- Is a ratio of 1.8 in median KSP iterations over 100x cells 'mesh-independent'? mesh_sweep.sh's own criterion is 'MESH-INDEPENDENT if the ratio is ~1', with no tolerance band stated. The project has not decided what ratio constitutes a failing preconditioner.
- memory `newton-vs-picard-benchmark-state.md` says 'needleDBD under Newton -- it refuses multiRegionPoisson' while docs/CAPABILITIES.md records multi-region under Newton as DONE at commit 10b7b2c with a five-row passing gate. Both are dated 2026-09-11. Which is current, and does the memory need a SUPERSEDED marker (rule 19)?
- The PCSHELL's 17 line-search failures (-6) are 'NOT yet diagnosed'. It is the stronger preconditioner on every linear-algebra measure, so this is the single cheapest remaining performance lever and nobody has looked at it.
- The cause of the 1.15M convergence wall remains unknown. The Schur-starvation hypothesis is refuted (-4% Krylov for +14% wall), the cap was a separate bug, and the per-FIELD scaling of the matrix-free differencing (phi ~ 1e4 V vs n ~ 1e13-1e20; MatCreateSNESMF picks ONE h from ||u||) is listed as candidate 3 and is still untested. Note also that the whole 1.15M verdict predates the kspMaxIt fix and should be re-tested with kspMaxIt 1000 before anything else is tried.
- Whether the pseudo-transient dt ceiling really moves ~10x with kspMaxIt 1000 (making the ~12 h steady run ~1 h) is INFERRED from the 21% -> 1.0% failure-rate change at dt=1e-10, not measured. dt=1e-9's earlier 74% failure rate is also listed as suspect and untested.
- Parallel defect 2 (MPI_ERR_TRUNCATE past the boundary-role globalPath fix) was never diagnosed; docs/design/parallel-blockers.md lists three candidate sources. Until it is fixed, 'every result in this thread is from a SERIAL run, and parallel execution of this solver is unverified' — yet scale_r3/scale_r4/positiveStreamer_fixedMesh all carry numberOfSubdomains 8 and the Newton cost numbers (142 s/step at 449k) do not state whether they were serial or 8-rank. That is a material ambiguity in the headline cost model.
- Whether `plasmaSimulationProfiler::report()`'s reduce-inside-a-map-loop has ever actually run in parallel with divergent key sets. The two conditional sections (`if (chem_)`, `if (finalIter && gasHeating_)`) are gated on model-level flags that should be rank-uniform, so it may be safe today — but no test exercises it, and rule 31 exists because this exact shape cost a full session elsewhere.
