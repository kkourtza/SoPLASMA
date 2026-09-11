---
name: perf-compare
description: Rule C3 made executable — run or judge a SoPlasma performance/scaling comparison at TWO problem sizes at matched Courant, metric = iteration count, cost = per unit of simulated time. Use before quoting any speed, cost, preconditioner or scaling number, and before proposing an optimisation (step 0 is the refuted list).
allowed-tools: Bash Read Grep
---

# perf-compare — C3 made executable

**The question this answers:** does this change help *at the size the solver exists for*,
and is the number I am about to quote a measurement of the thing I think it is?

**Rule C3:** a performance claim states its SCALING and is measured at TWO problem sizes
at MATCHED COURANT. The metric is mesh-independence of the **ITERATION COUNT**, not wall
clock. Cost is reported **per unit of SIMULATED TIME at equal accuracy, never per step** —
Newton costs ~4-5x per step (10x measured at 449k), so per-step is rigged for Picard.

**Do NOT use this skill for:** a correctness/regression check (that is the V&V bed), a
crash or NaN (that is a diagnosis — smallest bed, add the print, rule A6), or a claim you
are about to make from ONE problem size. One size is not a weaker version of this
procedure; it is not evidence (rule A2).

---

## 0. BEFORE proposing an optimisation, read the refuted list (bottom of this file)

Each entry below was refuted by exactly this two-size protocol, after looking like a large
win at 2,000 cells. Re-proposing one is the most expensive mistake available here.

```bash
grep -n "REFUTED\|refuted" /home/kkourtza/soplasma-scratch/.claude/mined/performance.md
```
PASS: your idea is not in the list, or you can say what is different about it.

## 1. Name the two sizes, the bed, and the START (cold vs warm) in one sentence

Cell counts verified on disk: grubert2009_pseudo 2,000 · mesh_20240 20,240 ·
mesh_200000 200,000 · scale_r3 211,120 · scale_r4 449,413 · streamer_base /
positiveStreamer_LMEA_fast 16,900 · positiveStreamer_fixedMesh 1.15M (~70 s/step).

```bash
grep -m1 note /home/kkourtza/soplasma-scratch/validation/<arm>/constant/polyMesh/owner
```
PASS: prints `nCells`. A bar measured WARM at 1.15M cannot be applied to a COLD 449k arm —
name the bed and the start in the same sentence as the dt.

## 2. Hold the Courant number fixed — scale dt with dx, `adjustTimeStep false`

Refining at fixed dt RAISES Co: measured Co_conv(e) 6.50 at 2,000 cells -> 20.62 at
20,240 — ratio 3.17, exactly the 3.16x refinement. "10x cells doubled the iterations"
measured the CFL, not the mesh; holding Co fixed reversed the conclusion.

```bash
NSTEPS=20 NRELAX=50 DT=1e-12 /home/kkourtza/soplasma-scratch/validation/mesh_sweep.sh
```
It does the whole protocol: `dt = DT*400/NX` (mesh_sweep.sh:85), `adjustTimeStep false`,
a developed state mapped from t=1.954497e-06, and NRELAX steps discarded. cwd is
irrelevant (`cd "${0%/*}"` at line 42). Runtime ~95 min for three arms (29.83 / 441.12 /
5111.73 s measured).
**--fast (B6):** `NSTEPS=5 NRELAX=10` and edit `GRIDS=("400 5" "1265 16" "4000 50")`
(line 72) to drop the 4000x50 arm -> ~8 min. Two sizes is the MINIMUM C3 accepts.

Then PROVE Co matched, per arm:
```bash
LOG=$(ls -t logs/log.run logs/log.soPlasmaFoam log.soPlasmaFoam 2>/dev/null | head -1)
grep -oE 'Co_conv \(e\): +[0-9.e+-]+' "$LOG" | tail -3
```
PASS: the tail values agree between arms to ~1%. If they do not, the sweep measured CFL.

## 3. Prove both arms ran the SAME binary (four of seven arms once did not)

```bash
ls -la --time-style=full-iso $FOAM_USER_APPBIN/soPlasmaFoam $FOAM_USER_LIBBIN/*.so
head -20 <each arm's log> | grep -i 'Exec\|Date\|Time'
```
PASS: every .so and the solver are OLDER than every arm's start time. The 449k Newton
ladder had four arms against a binary replaced at 15:48 compared with 16:2x arms — the
comparison is void, re-run every arm on one binary before reading the ladder.

## 4. Diff the RESOLVED settings and the START-UP BANNERS of both arms

A **model-dependent DEFAULT is an uncontrolled variable the moment the model differs**
(A2). Diff resolved configuration, not the keys you changed.
```bash
diff <(sed -n '1,120p' armA/logs/log.run) <(sed -n '1,120p' armB/logs/log.run)
grep -n 'newtonSolver' validation/<arm>/system/plasmaSimulationControls
```
PASS: the only differences are the swept variable. Defaults are `kspMaxIt 100`,
`assembledPmat true` (snesNewtonSolver.C:1134-1135) — raise `kspMaxIt` to 1000 in the
`newtonSolver` dict before quoting ANY Newton failure rate.

## 5. Prove the arm is not vacuous — the knob was in force (A1)

```bash
grep -ac "outerSolver newton (SNES)" "$LOG"
```
PASS: > 0. Zero means Picard ran and the arm says nothing about Newton (a cold ladder at
a dt where Picard's first step dies is vacuous by construction). Keep `-a`: a truncated
log contains NUL bytes. Output identical to baseline in every digit is proof the knob is
not wired — go read the code, do not re-run.

## 6. Check the iteration histogram is not CENSORED before reading its tail

```bash
grep -o 'Linear solve .* iterations [0-9]*' "$LOG" | awk '{print $NF}' | sort -n \
  | awk '{a[NR]=$1} END{print "n="NR, "med="a[int(NR/2)], "p90="a[int(NR*0.9)], "max="a[NR]}'
```
PASS: `max` != the configured cap. If max == cap the distribution is truncated and its
tail is not data — the "bimodal preconditioner cliff" was this artefact. Reference
(grubert2009_pseudo, 400 steps): cap 100 -> med 4 / p90 14 / max 89; kspMaxIt 1000 ->
med 4 / p90 14 / max 115, failures 10.6% -> 1.0% for +2.5% wall.

## 7. THE METRIC: median KSP iterations per step and SNES per step, at both sizes

```bash
cd /home/kkourtza/soplasma-scratch/validation
for d in mesh_2000 mesh_20240 mesh_200000; do echo -n "$d: "; awk '/^Time = /{s++} /Linear solve .* iterations [0-9]+/{if(s>50){n=$NF;a[++k]=n;t+=n;if(n>mx)mx=n}} /\(SNES\): reason .* iterations [0-9]+/{if(s>50){sn+=$NF;sc++}} END{m=s-50; asort(a); printf "steps=%d measured=%d KSP med=%d mean=%.1f max=%d | SNESits/step=%.2f\n", s, m, a[int((k+1)/2)], t/k, mx, sn/m}' "$d/logs/log.run"; done
```
PASS: a FLAT median across arms. Measured 2026-09-11 (INDICATIVE — mesh_sweep.log ends
mid-200k arm, so the project's own table was never printed): med KSP 5 / 6 / 9 and
SNES/step 4.00 / 4.00 / 6.28, i.e. ratio 1.8 over 100x cells. The project has NOT decided
what ratio counts as failing; state the ratio, do not label it.

## 8. Cost, if you must quote one: per ns of SIMULATED TIME, and say which clock

```bash
grep ExecutionTime "$LOG" | tail -1      # ExecutionTime = CPU, ClockTime = wall
```
The CPU/wall ratio varies **2.3x to 26.0x between arms of the same sweep** (hidden
hypre/BLAS OpenMP on 32 cores): 69.52/30, 11469.62/441, 29574.66/5109. If the ratios
differ, compare CPU seconds or set `OMP_NUM_THREADS=1`. Exact wall + peak RSS:
`/usr/bin/time -f "    wall %e s  maxRSS %M kB" soPlasmaFoam > logs/log.run 2> logs/log.time`
(~10.3 kB/cell: 104,852 kB at 2k, 2,059,916 kB at 200k).

Reference cost models to compare against — both arithmetic, not mystery:
* **JFNK:** 1 Krylov iteration = 1 full residual assembly = ~1 whole Picard step. At
  1.15M: 50 SNES x 200 KSP = 10,000 evals x 4.16 s = 11.6 h/step vs 10.06 h observed.
  Compute `nSNES x nKSP x t_residual` BEFORE launching; if it exceeds your budget, shrink
  the bed, not the tolerance.
* **Newton vs Picard, 449k, dt 1e-11:** 142 s/step vs 14.6 s/step. Warm Picard at 1.15M
  survives dt 5e-11, so Newton must hold ~5e-10 merely to BREAK EVEN.

## 9. Do NOT attribute Newton cost from the PROFILING REPORT

`src/profilers/plasmaSimulationProfiler` is always on, cumulative, and prints only at a
normal `End` (soPlasmaFoam.C:1051) — a killed/timed-out/crashed run yields nothing. Every
top-level `start`/`stop` sits in the Picard branch (soPlasmaFoam.C:878-898, inside the
`else` opening at :876); `grep -rn plasmaSimulationProfiler src/numerics/newtonSolverPETSc/`
returns ZERO. Coverage check:
```bash
grep -A12 'PROFILING REPORT' "$LOG" | awk '/->/ {s+=$(NF-2)} END{print "instrumented CPU s =", s}'
grep ExecutionTime "$LOG" | tail -1
```
PASS: ratio > 0.5 means the profile is an accounting. Measured under Newton: 2.4-3.7%.
For Newton, profile from `-snes_monitor` / `-ksp_converged_reason` / `-log_view` via the
case's `petscOptions`, and from `ExecutionTime` deltas per step.

## 10. First, though: is the cost even where you are looking?

```bash
grep -oE 'deltaT set by: +.*' "$LOG" | sort | uniq -c | sort -rn | head
```
The cheapest diagnostic in the tree — optimise the limiter that BINDS. (Reference,
adt20k_base: `853 Co_conv (energy)` / `1 growth cap (errorMaxGrow)`; needleDBD's
`maxSpeciesCo 100` never bound once in 322 steps.) The largest recorded lever is the
outer-coupling contraction: 22 correctors/step with Aitken omega damped to 0.148 clamped
dt to 1/83 of accuracy — ~1000x, dwarfing Newton's 4-5x. Also check `ddtSchemes`: BDF2 is
a measured 5.5x over Euler.

## 11. Record it (D1/D2)

A COMPARE.md in the case dir (question, ABSOLUTE baseline paths, varies, matches,
reference numbers, extraction command), the table in `docs/CAPABILITIES.md`, and the
dated line in `PROGRESS.md`. Every measurement carries its date.

---

## REFUTED AT TWO SIZES — do not re-propose (check this BEFORE step 1)

| proposal | 2,000 cells | 20,000+ cells |
|---|---|---|
| ILU(1) on the transport split | -33% steps | +3.5% |
| `selfp` Schur approximation | -22% steps | KSP med 11 -> 12 |
| per-cell residual/state scaling | -46% steps, SNES 1887->1224 | +1.4% mean dt, KSP med 13 vs 13 |
| ILU(1) + `selfp` together | — | exactly baseline (they interfere) |
| more `-fieldsplit_phi_ksp_max_it` | — | 449k: 8 cycles = -5% SNES, -4% Krylov, **+14% wall** |
| `sourceAwareScaling` | removed 2026-09-11 | `max(ddt, src)` is a NO-OP: src/ddt = 0.833 and 1.6e-4 |
| `extrapolateGuess` for block imbalance | — | 38 vs 40 converged over 11 steps |
| tightening Poisson tol to 1e-16 | 286x iterations | Ey unchanged; below double precision |

Beds already on disk: `validation/pc_ilu1{,_20k}`, `pc_selfp`, `selfp_20k`, `pcs_20k`,
`pcs_200k`, `adt20k_base`, `adt20k_ilu1`. **Do not re-run 1.15M under Newton to see if it
works** — answered twice (warm 10 h, cold 3 h, both still step 1); the wall is mesh size
and the open question is the PRECONDITIONER. Still open, untested: per-FIELD block scaling
(phi ~ 1e4 V vs n ~ 1e13-1e20; `MatCreateSNESMF` picks ONE h from ||u||) — a DIFFERENT
thing from the refuted per-cell scaling.
