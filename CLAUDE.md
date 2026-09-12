# SoPhy — SoPlasma + SoEEDF + SoGlobal. THE governing file.

**ONE file governs the whole multi-module project, and this is it.** `Projects/SoEEDF/CLAUDE.md`
is a stub that imports it. Never copy it — a second copy WILL go stale, and did: a session once
ran with 10 rules when the file had 44.

**Read `PROGRESS.md` first. Every session. It is the only orientation file.**

## What this is

A commercial-grade plasma-multiphysics framework, offered eventually as a web/cloud SaaS (the
SimScale model: the case is set up in the browser, the solver runs in the backend, results are
visualised online). More modules will follow (CFD, thermal, PIC, high-frequency EM), so **every
module shares one vocabulary and one infrastructure** — that is what makes a web front-end
possible at all, and why G3 is a directive rather than a preference.

| module | what it is |
|---|---|
| **SoPlasma** | the self-consistent plasma-fluid solver, OpenFOAM-based. Initial version: arXiv:2607.05137 (Kourtzanidis & Pasolari); much has advanced since. |
| **SoEEDF** | a two-term electron Boltzmann solver (BOLSIG+ style), embeddable C++17 — EEDF and transport/rate coefficients from LXCat cross-sections. Builds standalone, **no OpenFOAM dependency**. |
| **SoGlobal** | **NOT IMPLEMENTED — planned.** A 0-D global solver on the same core. No tree exists. |

| tree | what it is |
|---|---|
| `/home/kkourtza/soplasma-scratch` | **THE LIVE SoPlasma TREE.** All SoPlasma work happens here; this file lives here. |
| `/home/kkourtza/Projects/SoEEDF` | the Boltzmann repo. Also holds `docs/rules-postmortems.md`, `Literature/`, `tools/compare_cases.py`. |
| `/home/kkourtza/Projects/SoPLASMA` | a STALE clone of the same remote — **DO NOT EDIT** |
| `/home/kkourtza/square_NS_Diel`, `..._2` | a separate FORTRAN solver — **NEVER TOUCH** |

**"BEST" has four axes**, so a trade-off can be argued rather than asserted: **accuracy** (vs
experiment, benchmarks, analytic truth), **robustness** on hard cases without per-case tuning,
**CPU cost** per unit of resolved physics, and **ease of use** (G1). When two conflict, say which
is being traded. Do not silently optimise one.

## GOVERNING DIRECTIVES

Not preferences to weigh. When a design choice conflicts with one, the choice is wrong.

**G1. Fewest files, fewest hand-written choices.** Nearly every file OpenFOAM needs must be
GENERATED. Prefer, in order: **DERIVE** it from something already stated; **DEFAULT** it from
documented best practice; **GENERATE** the file. Require an input only for genuine physics the
user alone can know. **A setting stated in two places is a defect even when both copies agree**,
because they will not agree later.

**G2. Two-layer case architecture.** The user writes a small SEMANTIC layer; every dictionary is
generated from it. For boundaries the user says WHAT a surface IS — electrode, dielectric barrier,
floating conductor, open boundary — and each module decides HOW it treats that surface. Anything
with literature values becomes a LIBRARY with cited numbers. **Every shipped value is either
PHYSICS-RIGID** (computed from a stated relation) **or LITERATURE-DERIVED with a citation in the
entry.** Where the literature spread exceeds the difference between entries, say so instead of
shipping a point value. Never ship a value measured for one system as general — argon-on-copper is
not air-on-acrylic, and shipping it as a default is worse than a constant because it looks
authoritative.

**G3. One vocabulary across every module.** A concept has ONE name everywhere — dictionary key,
field name, CSV column, log line, doc. Search for the existing name before inventing one; where two
modules disagree, both are wrong until reconciled. Corollary: genuinely DIFFERENT quantities must
not share a name either — `I_collected` was renamed away from `I_plasma` because it invited
confusion with Sato's `I_total`.

## THE RULES — ORGANISED BY WHEN THEY FIRE

**Not by subject.** Rules attached to a recognisable MOMENT fire; rules attached to a TOPIC get
violated. That is the whole point of the 2026-09-10 consolidation (44 → 19), and rule 41 was broken
hours after being written because nothing about "the log is too big" prompts a lookup under
"snapshots".

**OLD NUMBERS REMAIN THE STABLE IDENTIFIERS** — cited in commits, code and memories; the lookup
table maps all 44. **The measured evidence is `../Projects/SoEEDF/docs/rules-postmortems.md`**, keyed
by original number. Read it when a rule seems arbitrary; it is why these are trusted.

**Procedures live in skills, not here** — `/orient`, `/build`, `/regression-gate`, `/validate`,
`/new-case`, `/doc-sweep`, `/commit`, `/perf-compare`, `/save-state`. This file carries only what
must fire unprompted.

### 0. IF YOUR CONTEXT WAS JUST REBUILT

**R0) READ `PROGRESS.md`. IT IS THE ONE FILE.** Triggers: a compaction/summary resume, a
`--continue`/`--resume`, a fresh session, or your own sense that you are reaching for something and
cannot recall whether it exists. It holds the live state and points onward to `docs/CAPABILITIES.md`
(what exists, the tooling, what was refuted). `/orient` does it for you.

A compaction summary preserves the NARRATIVE and LOSES THE INVENTORY. On 2026-09-10 that produced a
proposal for a capability already implemented in three places, and a hand-rolled cold start around a
nonexistent utility while the generator sat in the tree. The user: *"its like you dont remember
anything of what weve done anymore... i cant follow and continue like this to work with you."* The
cost is not the minutes; it is that the user cannot rely on continuity. **Keep both files true (D5).**

### A. BEFORE YOU BELIEVE A RESULT

**A1) A MEASUREMENT IS NOT EVIDENCE UNTIL YOU HAVE SHOWN IT MEASURED WHAT YOU THINK.** *(2, 13, 15,
21, 23, 24, 42)* Six checks, each of which has caught a wrong result alone:
* **the tool RAN** (21) — a positive success marker plus exit status, never the absence of errors.
* **the knob MOVED** (42) — and the proof goes in the report. **Output identical to baseline in
  every digit is proof the knob is not wired: go read the code, do not re-run.**
* **the answer COULD have differed** (23) — exercise it where it MUST come out different.
* **a DEFAULT is what the resolver returns** (24) — not an initialiser, comment or tutorial.
* **the control is NAMED** (15) — state what differs and what must match, or say "no control"; and
  report every magnitude against the physical scale it should be judged by.
* **the baseline was READ, not guessed** (13) — from `COMPARE.md`, never from directory names.
* and results get VALIDATED against data (2), asking for it when needed.

**A2) COMPARE ONLY AT EQUAL SAMPLES, EQUAL CONDITIONS, AND ON OUTPUTS.** *(43)* Five misreadings came
from unequal step counts — wait for a common endpoint. Refining a mesh at fixed dt RAISES the Courant
number, so "10x cells doubled the iterations" measured the CFL, not the mesh; holding Co fixed
reversed it. **A MODEL-DEPENDENT DEFAULT is an uncontrolled variable the moment the model differs** —
diff both arms' resolved settings and start-up banners. And read the FINAL state at full precision: a
comparison off the start snapshot, or off 6-significant-figure log output, compares its own inputs or
its own print precision.

**A3) DO NOT CONCLUDE AHEAD OF THE EVIDENCE; RECONCILE BEFORE EXPLAINING.** *(11, 14, 20, 22)* No
conclusions before the test is run or the code read (11). **Contradiction protocol** (14): reproduce
the prior number in the same pass FIRST; if it does not reproduce, say "I cannot reproduce X" and
STOP — never state a conclusion alongside an unreproduced prior. Never contradict a documented
conclusion from RECALL (20); that a parameter still exists is not that the conclusion holds. Never
explain away a discrepancy before its diagnostic is reconciled (22) — the diagnostic is the first
suspect, not the physics.

**A4) READ THE SOURCE BEFORE CITING IT.** *(32, 35)* Never conclude from a paper not read in full —
an abstract or another paper's citation is a lossy secondhand account; ask for the PDF. When unsure
of a detail INSIDE a paper the user has, **ASK THEM TO LOOK** rather than inferring. → `literature-analyst`

**A5) MATCH THE DIAGNOSTIC TO THE FAILURE CLASS; INVARIANTS BEFORE PHYSICS.** *(27, 39)* Memory
corruption → `-malloc_debug`/valgrind/ASan, never a `free()` backtrace; SEGV → `gdb -batch -ex run
-ex "bt full"`; SIGFPE → suspect a huge sentinel before a zero denominator; N candidates → change ONE
and measure. **READ THE TOOL'S OWN HINT FIRST.** When a case misbehaves check IN ORDER: (a) did every
per-step update actually RUN — a counter, not inspection; (b) is every DERIVED field consistent with
its sources; (c) is any field on its TRIVIAL solution (vacuum field, floor, clamp, initial value).
Only then question the physics. → `diagnostician`

**A6) INSTRUMENT, DO NOT GUESS — AND DEBUG ON THE SMALLEST BED.** *(user, 2026-09-11)* An instrument
returns a fact every time, **including when it stays silent** — a NaN locator that printed NOTHING
through six `DIVERGED_FNORM_NAN` failures proved the residual's COMPONENTS were finite and the NORM
had overflowed. Measured that day: four guesses at Newton's streamer failures, all four refuted, one
build-run-read cycle each; three instruments, all three decisive on first use. **Bed size is part of
the rule** — all of it was found on the ~2 s coarse bed after hours lost at minutes-per-attempt on
1.15M cells. Ask whether a smaller bed reproduces it before debugging anything.

**A7) STATE WHICH TIER OF EVIDENCE SUPPORTS A CLAIM; NEVER QUOTE A RESULT ABOVE ITS TIER.** Tiers are
not interchangeable: **(1) analytic unit bed** — closed-form ground truth (`verification/fluxScheme*`,
`testWallFlux`, `check_series_stack.py`); **(2) verification** — order of accuracy, mesh/time
convergence, symmetry, conservation; **(3) validation** — experiment or a published benchmark
(Grubert 2009, the arXiv:2607.05137 streamer, Biagi/Phelps); **(4) plausibility** — it ran and nothing
looked wrong, which is **not evidence**. A numerical claim needs tier 1–2, a physical claim tier 3;
**a converged run is tier 4.** State the observed order `p`, the grids it was fitted on, and what was
held fixed (A2). A unit bed with ground truth beats a CFD case for a numerics question — prefer it,
and say when no valid analytic reference exists rather than implying one.
→ `numerical-analyst`, `physics-validator`

### B. BEFORE YOU TOUCH A RUNNING CASE OR WORKING CODE

**B1) SNAPSHOTS: WRITE TEN, RESTART FROM THEM, AND TO STOP A RUN FORCE ONE.** *(41)*
`writeInterval <= endTime/10`, set when the case is configured. Restart from a snapshot rather than
re-running from t=0 whenever it cannot change the conclusion — sound for a defect that recurs from any
nearby state, NOT when the conclusion depends on trajectory history; and a restart REWRITES the
snapshot it starts from, so copy it read-only if reusing. To stop: `stopAt writeNow` in `controlDict`
(`runTimeModifiable`; then `touch system/controlDict`) **writes the state AND exits in one action — it
IS the kill.** Confirm the written directory matches the log's last reported time. To change something
and keep running, lower `writeInterval` live. **A run's LATEST TIME DIRECTORY is the only thing that
survives being stopped.** And check whether stopping is needed at all — a log eating the disk is fixed
by `: > the.log`.

**B2) WHEN A PROPOSAL TURNS OUT LESS PROMISING THAN WHEN YOU PITCHED IT, STOP AND ASK.** *(44)*
Trigger: the expected payoff changed after the work was agreed. Mandatory when it risks the SOLVER (a
wrong change corrupts silently) or is a huge effort relative to what was described. Say what changed,
the revised cost/benefit, and the alternatives including doing nothing — then wait.

**B3) BEFORE LAUNCHING: STATE THE PARAMETERS, CHECK THE DISK, BACKGROUND ANYTHING OVER 1–2 MINUTES.**
*(9, 10, 17, 36)* State the FULL parameter set first — circuit type and values, applied
voltage/current and profile, `simulationType`, gas/pressure/temperature, `endTime`, dt strategy — for
every case, arm and **relaunch**. Check storage (in WSL `df` is misleading; the real limit is the
Windows C: drive). Prefer the SHORTEST run that proves the point, and ask whether a UNIT bed answers
it better; state the discriminating observable and the earliest time it is visible — **and that
time is when the FIRST CHECK GOES, not a round interval and not completion.** An arm that can die
on STEP 1 (a restart, a new dict key, a rebuilt library) is checked in ~90 s; a failure RATE at
~10% of the steps; a trajectory or endpoint verdict only by an armed waiter. Measured 2026-09-12:
an immediate check answered a restart-vs-dt-limit question in 90 s that a completion waiter would
have answered 20 minutes later — and reported 7 GB of 30 GB free, the quantitative reason ten
concurrent 449k arms had frozen WSL twice that day. **The first check is also the resource check.**
→ `/validate` §4b

**B4) READ WHAT EXISTS BEFORE PROPOSING; PROVE THE NEED BEFORE ADDING; AUDIT EVERYWHERE AFTER
FIXING.** *(30, 37)*
* **Before proposing any feature or declaring a gap, read what the solver already does**:
  `docs/CAPABILITIES.md`, `README.md`, `etc/boundaryRoles`, `docs/reference/`, `docs/design/*.md`
  (**check the `Status:` line — several are DESIGN with nothing implemented**). → `inventory-scout`
* **This applies to TOOLING, and before DOING the work, not only before proposing it.** Before writing
  any script or setting up any case, read CAPABILITIES' tooling section, `tools/`, and the repo-root
  `*.sh`. **Each script's header records the trap it was written to avoid** — re-deriving the script
  means re-earning the traps.
* **Never hand-edit a generated file.** The value goes in `configuration/config` as a `$key`; the
  boundary in `configuration/boundaries`. Then prove the key is not dangling (command card). A hook
  blocks what it can identify, but it only knows what `.gitignore` records.
* Before adding a mechanism the framework already has, prove BY MEASUREMENT that the framework's own
  fails in the configuration you would ship — the answer to a mechanism failing on bad input is a
  GUARD, not a parallel mechanism (30).
* **A bug is a property of a MECHANISM, not of the call site that exposed it** (37): fix every other
  use of the same mechanism in the same change.

**B5) NO SILENT REGRESSIONS.** Before any commit touching solver source, run `/regression-gate`. A
changed number is one of three things and they are not interchangeable: a **REGRESSION** (revert or
fix), an **INTENDED IMPROVEMENT** (update the baseline IN THE SAME COMMIT, with reason and date), or a
**STALE BASELINE** (say so). **If you cannot tell which, you do not get to guess — ask.** Never update
a baseline to make a test pass: that converts a failure into a silent permanent wrong answer, and a
baseline refreshed without a stated reason has stopped being a test.

**B6) EVERY TEST HAS A `--fast` MODE, AND THE FAST MODE IS WHERE YOU DEBUG.** Every bed, sweep and
suite needs a documented seconds-scale variant, carried in the command card beside the full form. Add
it when you add the test. **A fast mode is for DEBUGGING and must never be quoted for physics** — the
coarse streamer bed does not resolve a streamer.

### C. CODE STRUCTURE AND THE FRAMEWORK

**C1) IN SoPLASMA, WORK INSIDE OPENFOAM'S OWN STRUCTURES — NEVER BUILD A PARALLEL WORLD.** *(29, 31,
40)* **Scoped to SoPlasma; does NOT apply to SoEEDF**, which must keep building with no OpenFOAM
dependency.
* Express physics with `fvm::`/`fvc::`, carry data in `GeometricField`s with real `dimensionSet`s,
  read settings through `dictionary` and runtime selection, follow existing naming (G3). Confine a
  third-party library behind a thin MECHANICAL boundary — `snesBridge.H` is the worked example.
* **A control object's exit path is part of its CONTRACT** (29): never leave its loop with
  `break`/`return`/`goto`; set a flag and DRAIN the loop. `pimple.loop()` promises a `finalIter`, and a
  raw `break` voids it — that froze Poisson's source for 99.96% of steps. Verify per-step INVARIANT
  COUNTS afterwards, not that it runs.
* **A collective must never sit downstream of a guard on a LOCAL quantity** (31): `reduce`/`gSum`/
  `gMax`/`returnReduce` must be reached by every rank equally. Hoist it out rather than guarding it
  harder, and **test two decompositions** — the same bug deadlocks under one and raises
  `MPI_ERR_TRUNCATE` under another, and one variant faked a PHYSICS failure.
* **`fvMatrix::residual()` is BROKEN in parallel** (~21 orders of magnitude, independent of PETSc) —
  compute any outer/Newton residual by explicit `fvc::`.
* **All solution state must be idempotent within a step**; a new field MUST join the discard path.
  → `code-quality-reviewer`, `openfoam-implementer`

**C2) ADVANCED NUMERICS, CPU EFFICIENCY, CORRECT PHYSICS.** *(4, 7)* Restated from the objectives above.

**C3) A PERFORMANCE CLAIM STATES ITS SCALING AND IS MEASURED AT TWO PROBLEM SIZES AT MATCHED
COURANT.** The test is mesh-independence of the **ITERATION COUNT**, not wall clock. Measure **per unit
of simulated time at equal accuracy, never per step** — Newton costs ~4–5x per step, so a per-step
comparison is rigged for Picard. An optimisation that helps at 2,000 cells and not at 20,000 is not an
optimisation; three were refuted exactly that way. → `performance-analyst`

### D. RECORDS THAT OUTLIVE THE CONVERSATION

**D1) DOCUMENTATION IS PART OF "DONE" — FOR THE USER.** *(1, 12)* The whole surface: module
`README.md`, `docs/reference/` for any dictionary key, the `docs/design/` note, any report the result
bears on. A key rename needs a SWEEP in the same change; anything left behind is not stale, it is an
INSTRUCTION a later session will follow. Every comparison case carries a `COMPARE.md` in its own
directory, written WHEN THE CASE IS CREATED: the question, the baseline's ABSOLUTE path, the reference
numbers, the extraction command. → `/doc-sweep`, `docs-sweeper`

**D2) EVERY MEASUREMENT CARRIES ITS DATE; SUPERSEDED TEXT IS MARKED, NEVER DELETED.** *(18, 19)* A
measurement without a date will be read as a live conclusion. Mark a superseded result `SUPERSEDED BY
<what>` in place, stating what changed and what still holds — a finding that quietly disappears leaves
no way to recognise its stale copies elsewhere.

**D3) EVERY CASE WRITES A MACHINE-READABLE TIME SERIES, AT PROBES.** *(6, 26)* One to three probes,
each versus time: every species density, `chargeDensity`, `Emag`, the reduced field, electron energy,
currents; plus run-level `dt` **AND WHICH LIMITER SET IT**, outer iterations, rejected steps. Probes
rather than fields (a field history is unwritable in 3-D); keep min/max too, but a reduction throws
away WHERE and where is usually the answer. A DERIVED quantity sits NEXT TO its sources. **If
diagnosing a case requires parsing the log, the monitoring is the defect.** Validation cases also
produce publication-ready figures AND the python that generates them.

**D4) TELL ME WHEN THERE IS UNCOMMITTED WORK — DO NOT WAIT TO BE ASKED.** *(38, verbatim)* Check `git
status` and say so PLAINLY as its own statement: after any new file, after a logically complete unit,
before anything long-running, before any destructive git operation. A PROMPT, not licence to commit
unasked. **NEW UNTRACKED FILES ARE THE URGENT CASE** — they appear in no diff and are destroyed by a
clean. **STAGE EXPLICITLY: never `git add -A`**, never `git add <dir>` on a tree holding case output;
re-read `git diff --cached --name-only` and unstage what does not belong.

**D5) DOCUMENTATION FOR THE NEXT SESSION — UPDATE `PROGRESS.md` AFTER EVERY MEANINGFUL UNIT OF WORK.**
D1 is for the eventual USER; this is for the next session, and the register differs — what was tried,
what the number was, what is still unknown. Never marketing prose. The contract: **check off completed
items WITH THE DATE**; **note what worked, what didn't, and what is blocked and on what**; **record
failed approaches** as *"Tried X — didn't work because Y. Switched to Z."*; **add new tasks the moment
you discover them**; and **when stuck, keep a dated running log of attempts** so the next session
inherits the search rather than repeating it. `/save-state` does this and prompts D4.

**D6) SMALL, TESTABLE COMMITS.** One logical unit, explicitly staged, with a message saying WHAT
changed and WHY — and the measured number where there is one. Never mix a refactor with a behaviour
change: when the behaviour then moves, nothing says which half did it. If it cannot be described in one
sentence without "and", it is two commits. → `/commit`

### E. WORKING WITH THE USER

**E1) SAY WHAT YOU ARE ABOUT TO DO, BEFORE DOING IT.** *(16, 28)* State plainly and up front what you
are about to do, what question it answers, and the steps. Not narrated as you go, not inferred from
tool calls, not explained afterwards. **It applies EVEN WHEN the work was requested** — "verify Part A"
authorises the work, it does not say what verifying consists of. At minimum, one line naming the
question before any stretch of more than ~3 tool calls whose purpose is not obvious.

**E2) PLAN FIRST; EXECUTE IN SMALL VALIDATED INCREMENTS AFTER APPROVAL.** *(3, 5)* **One change at a
time.** Changing the timestep AND the table range together made neither result interpretable.

**E3) NEVER SWITCH THREADS WITHOUT ASKING.** *(25)* Finish the open thread, or say which you are
leaving and which you are starting, and WAIT. An open-ended "ok", "continue" or "next?" means CONTINUE
THE CURRENT WORK. The test: can I name the thread I am on, and did the user's last message point at
THAT thread?

**E4) RESEARCH IS A VALID OUTCOME, NOT A LOWER-PRIORITY ONE.** *(33)* The user is a professor;
publishing is an equally important track. Never treat "this needs genuine research effort" as grounds
to defer. Assess an idea on whether the physics is sound, whether a gap genuinely exists (checked, not
assumed), and whether it is tractable.

**E5) DELEGATE TO THE SPECIALIST AGENTS.** They carry the traps for their area, so delegating is also
how those traps get applied rather than remembered — and they run at the effort their work needs.
`numerical-analyst` (discretisation, preconditioner, integrator, order study) · `physics-validator` (a
physical claim, units, "is this result physical") · `diagnostician` (a crash, NaN, stall, wrong number)
· `literature-analyst` (a paper must actually be read) · `openfoam-implementer` (solver/BC/model code)
· `test-quality-auditor` (before believing a test) · `code-quality-reviewer` (before committing solver
source) · `performance-analyst` (any speed claim) · `docs-sweeper` (a rename needs the D1 surface swept)
· `inventory-scout` (**always**, before proposing: does the solver already do X?).

**E6) A QUESTION YOU ASK AND THEN BURY IS WORSE THAN NOT ASKING.** *(user, 2026-09-12)* Fires the
moment you want a decision. **First test: do you actually need one?** If the change is small,
reversible, and plainly what the user just asked for, MAKE IT AND REPORT IT — asking turns a
two-minute task into a lost thread. Measured twice in one session: a question was asked, the original
work continued in the same message, and the user found the question only by scrolling back —
*"your question is LOST in the chat"*.

If you genuinely need the decision, then:
* **Ask and STOP.** Never ask and keep working in the same message — the run status buries it.
* **The question is the LAST thing in the message**, on its own line, so it cannot be scrolled past.
* **Carry it forward.** Repeat it as a one-line `STILL OPEN:` in EVERY later message until it is
  answered. A question asked once has been asked zero times.
* **If it outlives the turn, it is a BLOCKER** — put it in `PROGRESS.md` §3 BLOCKED with what it is
  waiting on. That is what §3 is for, and `/orient` then resurfaces it at the next session start, so
  the question survives a compaction, a crash and a reboot rather than only the scrollback.

### Old-number lookup — all 44 resolve

1→D1 · 2→A1 · 3→E2 · 4→C2 · 5→E2 · 6→D3 · 7→C2 · 8→*dropped* · 9→B3 · 10→B3 · 11→A3 · 12→D1 · 13→A1 ·
14→A3 · 15→A1 · 16→E1 · 17→B3 · 18→D2 · 19→D2 · 20→A3 · 21→A1 · 22→A3 · 23→A1 · 24→A1 · 25→E3 · 26→D3 ·
27→A5 · 28→E1 · 29→C1 · 30→B4 · 31→C1 · 32→A4 · 33→E4 · 34→*ref doc* · 35→A4 · 36→B3 · 37→B4 · 38→D4 ·
39→A5 · 40→C1 · 41→B1 · 42→A1 · 43→A2/C3 · 44→B2

**Dropped deliberately:** rule 8 (save sessions for a vscode restart) — a harness concern that governed
no decision. Rule 34 (consult JC-PIC before deriving plasma theory) — a POINTER, moved to
`Literature/reference-codes-and-manuals/JC-PIC_Boeuf/README.md` where it will be found. Rule 17's
5-minute target — routinely impossible, and an unmet rule teaches that rules are optional.

## THE COMMAND CARD

Every command was executed or existence-checked on 2026-09-11; **one that cannot be run is not listed**,
because a wrong command here is an instruction a later session will follow. Each skill carries the traps
and `--fast` variants for its area.

```bash
# ENVIRONMENT — always both, project second so PETSC_DIR wins.
source /usr/lib/openfoam/openfoam2412/etc/bashrc && source $HOME/soplasma-scratch/etc/bashrc
# etc/bashrc does NOT put binaries on PATH; OpenFOAM's does. Python is ~/ct-env/bin/python.

# BUILD                                                                      -> /build
cd $HOME/soplasma-scratch && ./build-all.sh 2>&1 | tee /tmp/build.log
grep -E 'BUILD-COMPLETE|BUILD-FAILED|BUILD-INCOMPLETE|UNCOVERED|^FAIL' /tmp/build.log
#   PASS: BUILD-COMPLETE present AND no FAIL / BUILD-INCOMPLETE / UNCOVERED.
#   Exit 0 does NOT prove a library rebuilt:
find <srcdir> -name '*.[CH]' -newer <that .so>      # must return nothing

# TEST                                                             -> /regression-gate
./check-no-running-solvers.sh                       # exit 0 required (FORCE=1 overrides)
tools/run_electrostatics_tests.sh                   # "N ok, 0 failed"; exit 2 = env unsourced
testWallFlux && testWallLoss && testVibRelax && testCoulombHeating && testEmission
#   Seconds, mesh-free. `&&` not `;`, or a failing bed is silently skipped.
#   testAitken has NO pass/fail machinery — it ALWAYS exits 0. Never treat it as a pass.
verification/fluxScheme1D/Allrun                    # results.txt must have EXACTLY 504 RESULT lines
#   Compare to results.baseline.txt by RELATIVE TOLERANCE, never `diff`.

# RUN                                                                     -> /validate
./run-guarded.sh [maxLogMB] [maxSeconds]            # hard log + wall-clock cap
~/ct-env/bin/python validation/status.py <case>...  # safe on a live log
pkill -9 soPlasmaFoam; pgrep -c soPlasmaFoam        # never `pkill -f`: misses -parallel and
                                                    # can kill the calling shell
NCELL=130 NREFINE=1 CASE=$HOME/streamer-fast ./make-smoke-case.sh   # the ~2 s coarse bed
CASE=$HOME/streamer-fast ./rerun.sh                                # (A6/B6; never for physics)

# CASE SETUP                                                              -> /new-case
./Allrun-serial                                     # PASS: log contains '=== done'
# Prove no configuration/config key is dangling — an inert $key once invalidated a whole arm
# that exited 0. Use grep -rF '$'"$k"; inside double quotes "\$$k" expands to the shell PID.
for k in $(grep -oE '^[a-zA-Z][a-zA-Z0-9_]*' configuration/config | sort -u); do
  grep -rqF '$'"$k" system/ constant/ configuration/ Allrun-serial 2>/dev/null || echo "DANGLING: $k"
done

# SoEEDF — standalone, no OpenFOAM dependency (C1 does not apply there).
cd $HOME/Projects/SoEEDF && cmake -S . -B build && cmake --build build -j
```

**Three traps that bite outside a skill's moment.** `grep FAIL buildlogs/*` CANNOT work — `OK`/`FAIL`
go to build-all.sh's *stdout*, and the repo's `build.log` is a stale Aug-10 artefact that reports a
month-old build as current. `wmake src/numerics` silently does NOT build `libplasmaNewtonSolverPETSc`;
it is a separate target. And **`./Allwmake` is NOT a full build** — it builds 4 of 8 model directories
and prints "Build complete." anyway.
