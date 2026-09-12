---
name: validate
description: Launch or read a SoPlasma validation/comparison case under B3 + A1 + A2 + D3 — full parameter set stated, unit-bed check first, COMPARE.md on disk before the run, disk forecast, backgrounded launch with the watcher armed in the same action, then A1's six checks before any number is believed.
---

# /validate — launch and read a comparison case

**Answers:** what does this case say, and is that number believable?
**Not for:** a numerics question with a closed form (step 0 sends you to a unit bed); a build
(`/build`); a verification-baseline diff (`/regression-gate`); creating a case (`/new-case`).

    source /usr/lib/openfoam/openfoam2412/etc/bashrc && source $HOME/soplasma-scratch/etc/bashrc

## 0. Would a UNIT bed answer it better? (B3 / rule 17; A7)

A ~50-min streamer run asking whether the contraction ratio could replace Aitken's omega as the dt
governor was moved into `testAitken` and answered in **one second** — and was the STRONGER test: the
CFD could only ask "does rho agree with omega", while the bed sets the loop gain and has GROUND
TRUTH. It also found what the CFD could not: rho >= 1 is the wrong threshold, the loop stops
converging at rho ~ 0.38. Target under five minutes per validation test.

    testWallFlux && testWallLoss && testVibRelax && testCoulombHeating && testEmission

Seconds, mesh-free. `&&` not `;`, or a failing bed is silently skipped. Pass = the literal
`ALL PASS`/`PASSED` line AND exit status, never absence of errors (A1). **`testAitken` has no
pass/fail machinery and ALWAYS exits 0.** A bed is ground truth only if its GEOMETRY matches the
formula (plate2D bundles left+right+bottom into one `walls` patch → 2D Laplace, no analytic answer).

## 1. STATE THE FULL PARAMETER SET — every launch AND relaunch (B3, rule 36)

Print, before launching: circuit type and values (`ballastedElectrode`/`currentDrivenElectrode`,
R, C, V_src, compliance) · applied voltage/current and its **profile** (ramp table, rise time) ·
`simulationType` · `electronEnergyModel` (LFA/LMEA/LMEAC) · gas, pressure, temperature · mesh
(`grep -m1 nCells constant/polyMesh/owner`) · `endTime` · dt/Courant strategy (`adjustTimeStep`,
`maxCo`, floor/cap) · `minNumberDensity` · `gammaSEE`, `electronReflection` · the discriminating
observable and the **earliest time it is visible** — that is the kill point, decided at launch.
Sources: `configuration/config`, `configuration/boundaries`, `system/plasmaSimulationControls`.
Never hand-edit a generated dictionary (B4).

## 2. Write COMPARE.md into the case dir BEFORE running (D1; A1 "baseline READ, not guessed")

The 2026-08-30 failure compared against a guessed directory-name family instead of the declared
baseline: every number correct, conclusion backwards by 9000x (1.564e18 vs 1.69e14).

````
```compare
question:  <one sentence, answerable yes/no>
baseline:  /home/kkourtza/soplasma-scratch/validation/<control>   # ABSOLUTE, repeatable
varies:    <the one key>
matches:   minNumberDensity, driftDiffusionFluxScheme, electronEnergyModel, mesh, chemistry, gammaSEE, electronReflection
time:      2e-5
field:     n_e
region:    domain
```
````

Then prose: reference numbers with citation and date (D2), and the extraction command. Only 52 of 75
contracts are machine-readable; `validation/grubert2009` and `validation/mesh_M1` are not:

    for f in /home/kkourtza/soplasma-scratch/validation/*/COMPARE.md; do
      grep -q '```compare' "$f" && echo "BLOCK  $f" || echo "PROSE  $f"; done

`PROSE` → quote its "Baseline — ABSOLUTE PATHS" section by hand and **say it was not machine-checked**.
Diff the block against siblings: `coulomb_compare_withCoulomb/COMPARE.md` carries
`grubert2009_ps_cfs`'s verbatim — a copied contract is a wrong control read as a right one.

## 3. Disk, then kill what you no longer need (B3 / rule 9; class 11)

    df -h /mnt/c    # THE real limit. /mnt/d (1.9T) is archive only — never RUN from drvfs

`df -h /home/kkourtza` once said 713 GB free while C: was at 98% (13 GB) and WSL died mid-session;
deleting inside WSL frees ZERO on C: without a vhdx compact. Then, each as its **own simple command**
(inside a compound command `pkill` can kill the harness's shell and return 144):

    pgrep -ax soPlasmaFoam
    pkill -x soPlasmaFoam
    pkill -x mpirun

**Never `pkill -f`** — it matches the calling shell, and `pgrep -c -f "soPlasmaFoam -parallel"`
returned **3** with no solver running. `pgrep -x` silently fails for names >15 chars
(`plasmaChemistry0D`, 17) — use `pgrep -af "[p]lasmaChemistry0D"`; `check-no-running-solvers.sh` has
that same hole, so it can only ever see soPlasmaFoam, testWallLoss and testVibRelax.

## 4. Launch backgrounded AND arm the watcher in the same tool-call block (B3 / rule 10; class 11)

A promise with no mechanism is worse than no watcher — the user then waits instead of asking
(*"your watchers kind of suck. Nothing happened until I typed news!!!"*). `nohup ... & disown` is not
enough; **setsid** is required. Logs go OUTSIDE the case tree (a `cp -al` sibling shares the inode).
**`./run-guarded.sh` is broken — do not use it:** it hardcodes `CASE=$HOME/streamer-case`, which does
not exist (exits 1 at once), and its `RUN-EXIT=` reports `head`'s status, not the solver's.

    LOGDIR=<scratchpad>          # outside every case tree
    (cd "$CASE" && setsid nohup ./Allrun-serial > "$LOGDIR/$(basename "$CASE").setup.log" 2>&1 </dev/null &)
    nohup /home/kkourtza/soplasma-scratch/watchdog.sh "$CASE/logs/log.soPlasmaFoam" 50 \
          > "$LOGDIR/watchdog.log" 2>&1 &

`watchdog.sh <log> [maxMB=50]` prints `WATCHDOG: <log> exceeded 50MB -- killing solver` and exits 1,
else exits 0 when the solver stops (a 9.7 GB / 91M-line runaway earned it); it kills only
`soPlasmaFoam`, so `mpirun` survives. `./run-long.sh <caseDir>` also takes a case dir but runs
`./Allclean` first — never point it at a case whose time dirs you need.

Then as a **separate** call (background it, or `Monitor`) the wake-up — all four endings plus liveness:

    until grep -qE "^End|FOAM FATAL|signal [0-9]+|Segmentation" "$f"; do
      [ -f "$f" ] && [ $(( $(date +%s) - $(stat -c %Y "$f") )) -gt 600 ] && { echo STALLED; break; }
      a=$(pgrep -c -x soPlasmaFoam 2>/dev/null); a=${a:-0}
      m=$(pgrep -c -x mpirun 2>/dev/null);       m=${m:-0}
      [ "$a" -eq 0 ] && [ "$m" -eq 0 ] && [ -s "$f" ] && { echo "NO PROCESSES"; break; }
      sleep 60
    done

### 4b. WHEN the watcher first fires — set it from the earliest informative moment, not a round number

Arming a watcher is only half the rule. The other half is **the first check goes at the
earliest time the discriminating observable could appear** — which B3 already makes you
state when you design the run, and which is then routinely ignored when deciding when to
look. Get it wrong in either direction and the run teaches you nothing on schedule: poll
every few seconds and you drown in noise; wait for completion and you are blind to a
failure that was visible in the first second.

| what the arm can tell you | when it can tell you | so the first check is |
|---|---|---|
| a restart, a new dict key, a config or mesh change, a rebuilt library — **anything that can die on STEP 1** | seconds after the solver starts | **immediately** (~60–90 s, enough for start-up) |
| a failure RATE, a corrector count, an iteration histogram | once there are enough steps to BE a rate | ~10% of the planned steps |
| a trajectory, an accuracy verdict, "did it reach the common endpoint" | only at the endpoint (A2) | arm the waiter and do NOT poll |

**Measured 2026-09-12, which is why this exists.** A warm-restart ladder arm was launched
and a completion-waiter armed. The user asked why it was not checked at once; the immediate
check then returned the answer in 90 s — `picard dt=1e-11` running clean at 50 steps while
every `dt >= 5e-11` arm had SIGFPE'd at step 1, so the restart was sound and the crashes
were a real dt limit. A waiter would have delivered that ~20 minutes later. The same check
also reported **7 GB of 30 GB free with four solvers running**, which is the quantitative
reason ten concurrent 449k arms had frozen WSL twice that day.

**The first check is also the resource check.** Report free memory alongside the step count
whenever more than two arms of a large bed are running: this bed is ~1.5 GB resident per
solver, and the failure mode is not a slow run, it is the machine going down.

**The liveness fallback is the important half** — it catches endings whose log signature you did not
anticipate, which is by definition the ones that matter. `[ -s "$f" ]` stops it firing during
preprocessing (mesh, fields, seed, `decomposePar`, EEDF table rebuild run for MANY MINUTES with no
solver process; concluding death there once produced two `mpirun` instances writing the same
`processor*` dirs). Never `n=$(pgrep -c X || echo 0)` — yields `0\n0`, and every later `[ ]` dies
with "integer expression expected"; that recurred three times in one day. `rm -f` the log before
watching or you grep the PREVIOUS run's `FOAM FATAL`. The `[ -f "$f" ]` guard is required or the
staleness arithmetic aborts while the log does not yet exist.

**Fast bed (B6)** — ~2 s/step, where all four Newton defects were found after hours lost at
minutes-per-attempt on 1.15M cells; never quotable for physics:

    NCELL=130 NREFINE=1 CASE=$HOME/streamer-fast ./make-smoke-case.sh   # does `rm -rf $CASE`, no prompt
    CASE=$HOME/streamer-fast TIMEOUT=300 ./rerun.sh                     # keeps polyMesh, 0/ from .snapshot0

**To stop a run** (B1): `stopAt writeNow` in `system/controlDict`, then `touch system/controlDict` —
it writes the state AND exits in one action, it IS the kill. Confirm the written directory's name
matches the log's last reported time.

## 5. Read it (D3)

    python3 /home/kkourtza/soplasma-scratch/validation/status.py <case> [<case> ...]

(mode 0644 — invoking by path gives Permission denied; stdlib only.) One line per case:
`t= dt= steps= rej= crash= lim= V= j=`; safe on a live log (complete-float regex + `v < 1.0`
torn-line filter). It hardcodes `<case>/logs/log.soPlasmaFoam`, so it prints `no log` for every
`adt*`/`mesh_*` case (`logs/log.run`) and every tutorial (root `log.soPlasmaFoam`):

    LOG=$(ls -t logs/log.run logs/log.soPlasmaFoam log.soPlasmaFoam 2>/dev/null | head -1)
    grep -oE 'deltaT set by: +.*' "$LOG" | sort | uniq -c | sort -rn   # WHICH limiter binds
    awk '/current deltaT/{s+=$NF} END{print s}' "$LOG"                 # = final Time on retry-free runs
    echo "idle $(( $(date +%s) - $(stat -c %Y "$LOG") )) s"            # >600 s at 99.9% CPU = HUNG
    gdb -p <pid> -batch -ex "bt 8"                                     # stack BEFORE killing

Hung and slow look identical — OpenMPI busy-polls in `PMPI_Waitall`, ranks show `Rl`; a 28-minute
hang was reported as "slow". **If diagnosing needs the log parsed, the monitoring is the defect (D3)**:
the case must write probe time series and `postProcessing/dischargeCurrent/current.csv`.

## 6. Before believing ANY number

**Discharge pre-checks, before the plateau test.** The first Grubert LFA run held current constant to
0.00% over 36 us — what a converged steady state looks like — while `tableKey meanE` made `reducedE`
(~2e-18 SI) index eV-keyed tables: every lookup hit the floor, `S_iz` = 1.9e-210 against 1.1e+20,
`n_e` = 1.53e11 against a 1e11 floor. Assert (1) `n_e` >> `minNumberDensity`, (2) `S_iz` ~ `k_iz*n_e*N`.
Reject any headline equal to a clamp, seed or default (Te,max 1762.9 = (2/3)*2644.46 is the
mean-energy clamp pinned — a reading from non-converged steps).

**A1 — six checks, each of which has caught a wrong result alone:**
* **the tool RAN** — `runApplication` does NOT propagate the solver's status, so Allrun exit 0 means
  only THE SCRIPT finished: `grep -nE 'FOAM FATAL|sigFpe|sigSegv|Floating point exception' "$LOG"`.
* **the knob MOVED** — output identical to baseline in every digit proves the knob is not wired.
  **Read the code; do not re-run.**
* **the answer COULD have differed** — every no-op check needs a paired liveness check that MUST
  fail, at the worst point found by scan. `testWallFlux` had 32 passing checks while the defect
  shipped (all at Dd = 0, where right and wrong coincide); its first CONVICTS check sat at Pe=100
  where the forms differ by 1% — the real worst point is 18% at Pe=1.
* **a DEFAULT is what the resolver returns** — the start-up BANNER, not an initialiser or tutorial.
* **the control is NAMED** — from `COMPARE.md`; report every magnitude against its physical scale.
* **the baseline was READ, not guessed** — never from a directory name.

**A2 — equal samples, equal conditions, on outputs:**
* **equal samples** — a common SIMULATED time, never a common step count (five misreadings came from
  this). `reconstructPar -case <case> -newTimes` if only `processor*` hold it.
* **equal conditions** — refining at fixed dt raises Co: Co_conv(e) 6.50 at 2000 cells vs 20.62 at
  20240, ratio 3.17 = exactly the 3.16x refinement, i.e. the whole effect. Scale dt with dx, set
  `adjustTimeStep false`. Diff RESOLVED config **and the banners** — `adaptiveRelaxation` defaults to
  `hasLMEA`, so an LFA arm silently ran with no relaxation. The diff must show ONLY `varies:`:

      diff <(grep -vE '^\s*(//|$)' A/system/plasmaSimulationControls) \
           <(grep -vE '^\s*(//|$)' B/system/plasmaSimulationControls)

* **on OUTPUTS** — not the START snapshot, not 6-s.f. log output (`writePrecision 6` puts a ~1e-6
  relative floor under any field value, so that compares your own print precision).

**Enforce the contract mechanically.** `tools/compare_cases.py` does NOT exist in this tree — run it
from SoEEDF by absolute path, or state plainly "the contract was read by hand, not enforced":

    ~/ct-env/bin/python /home/kkourtza/Projects/SoEEDF/tools/compare_cases.py <caseDir>

Exit 0 with every reference line `OK`. Exit 1 prints `A DECLARED REFERENCE DID NOT REPRODUCE` → A3's
contradiction protocol: reproduce the prior number in the same pass; if it does not reproduce, say
"I cannot reproduce X" and **STOP**. Needs `~/ct-env/bin/python` (numpy + fluidfoam); plain `python3`
lacks them. Its messages cite OLD rule numbers (12→D1, 13→A1, 14→A3). Defaults `--tol 0.02
--max-offset 2e-11`.

## 7. Report

State the **tier** (A7), never quote above it — 1 analytic unit bed, 2 verification/order,
3 validation against experiment or a published benchmark, **4 = it ran and nothing looked wrong,
which is not evidence; say so.** Date every measurement (D2). A Grubert run is an **open reproduction
attempt**, never "validation": 58 of 60 validation cases are Grubert variants and NONE has a known
answer; measured, we are 1300x high on j and 25300x on n_e. The one genuine prediction is the GAP
VOLTAGE at I_set = I_op = 1.022e-6 A — the current is an input, so matching j proves nothing. If a
baseline moved and you cannot say why, **stop and ask** (B5); then `/save-state` (D5) and report
uncommitted work plainly (D4).
