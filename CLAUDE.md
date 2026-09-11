# SoPhy — SoPlasma + SoEEDF + SoGlobal. THE governing file.

**ONE file governs the whole multi-module project, and this is it.** It lives in
the LIVE SoPlasma tree beside the config it refers to; `../Projects/SoEEDF/CLAUDE.md`
is a stub that imports this. Never copy it — a second copy WILL go stale, and
has: a session once ran with 10 rules when the file had 44.

**Read `PROGRESS.md` first. Every session. It is the only orientation file.**

---

## What this is

A commercial-grade plasma-multiphysics framework, offered eventually as a
web/cloud SaaS (the SimScale model: the user sets the case up online, the solver
runs in the backend, results are visualised in the browser). Three modules today,
more later (CFD, thermal, PIC, high-frequency EM), so **every module must share
one vocabulary and one infrastructure** — that is what makes a web front-end
possible at all, and it is why G3 below is a directive and not a preference.

| module | what it is |
|---|---|
| **SoPlasma** | the self-consistent plasma-fluid solver, OpenFOAM-based. Initial version: arXiv:2607.05137 (Kourtzanidis & Pasolari); much has advanced since. |
| **SoEEDF** | a two-term electron Boltzmann solver (BOLSIG+ style) as an embeddable C++17 library — EEDF and transport/rate coefficients from LXCat cross-sections. Builds standalone, **no OpenFOAM dependency**. |
| **SoGlobal** | **NOT IMPLEMENTED — planned.** A 0-D global solver on the same core. No tree exists (checked 2026-09-04). |

### Where the code lives — READ BEFORE EDITING ANYTHING

| tree | what it is |
|---|---|
| `/home/kkourtza/soplasma-scratch` | **THE LIVE SoPlasma TREE.** All SoPlasma work happens here. This file lives here. |
| `/home/kkourtza/Projects/SoEEDF` | the Boltzmann solver's repo |
| `/home/kkourtza/Projects/SoPLASMA` | a STALE clone of the same GitHub remote — DO NOT EDIT |
| `/home/kkourtza/square_NS_Diel`, `..._2` | a separate FORTRAN solver. **NEVER TOUCH.** |

Both SoPlasma trees clone `git@github.com:kkourtza/SoPLASMA.git`, which is exactly
why the wrong one is easy to edit.

### "BEST" has four axes, so a trade-off can be argued rather than asserted

1. **Accuracy** — against experiment, published benchmarks, and analytic truth.
2. **Robustness** on hard cases (sharp electrodes, thin sheaths, stiff chemistry)
   without per-case hand-tuning.
3. **CPU cost** per unit of resolved physics.
4. **Ease of use** — how little the user must write (G1).

When two conflict, say which is being traded and why. Do not silently optimise one.

---

## GOVERNING DIRECTIVES

Not preferences to weigh. When a design choice conflicts with one, the directive
wins and the choice is wrong.

### G1. Fewest files, fewest hand-written choices

The user edits the LEAST number of files and hand-writes the LEAST number of
choices. Nearly every file OpenFOAM needs must be GENERATED from a few major
inputs. Prefer, in order: **DERIVE** it from something already stated; **DEFAULT**
it from documented best practice; **GENERATE** the file. Require an input only for
genuine physics the user alone can know.

**A setting stated in two places is a defect even when both copies agree**,
because they will not agree later.

### G2. Two-layer case architecture

The user writes a small SEMANTIC layer; every OpenFOAM dictionary is generated
from it. For boundaries the user says WHAT a surface IS — an electrode, a
dielectric barrier, a floating conductor, an open boundary — and each physics
module decides HOW it treats that surface. Anything with literature values
(materials, gas mixtures, mechanisms, boundary roles) becomes a LIBRARY with
cited numbers, not a number typed into a case.

**Every value in a shipped library is either PHYSICS-RIGID** — computed from a
stated relation — **or LITERATURE-DERIVED with a citation** carried in the entry.
A number without one is a number somebody remembered. Where the literature spread
is LARGER than the difference between entries, say so rather than shipping a point
value. And never ship a value measured for one system as general: argon-on-copper
is not air-on-acrylic, and shipping it as a default is worse than a constant
because it looks authoritative.

### G3. One vocabulary across every module

A concept has ONE name everywhere — dictionary key, field name, CSV column, log
line, doc. Search for the existing name before inventing one; where two modules
name the same thing differently, both are wrong until reconciled. A new module
inherits the vocabulary rather than bringing its own.

**Corollary:** genuinely DIFFERENT quantities must not share a name either.
`I_collected` (charge arriving on one conductor) was renamed away from `I_plasma`
on 2026-09-04 precisely because it invited confusion with Sato's `I_total`.

---

## THE RULES — ORGANISED BY WHEN THEY FIRE

**Not by subject.** That is the whole point of the 2026-09-10 consolidation
(44 → 19): rules attached to a recognisable MOMENT fire — rule 42 caught three
inert knobs in one day because it triggers on *I am about to believe a comparison*
— while rules attached to a TOPIC get violated. Rule 41 was broken hours after
being written, because nothing about "the log is too big" prompts a lookup under
"snapshots".

**OLD RULE NUMBERS REMAIN THE STABLE IDENTIFIERS**, cited in commits, code
comments and memories. The lookup table below maps all 44. The measured evidence
— the post-mortems — is in `../Projects/SoEEDF/docs/rules-postmortems.md`, keyed
by original number. Read it when a rule seems arbitrary; it is why they are trusted.

**Procedures live in skills, not here.** A skill is a moment-triggered rule with
the procedure attached, so `/validate`, `/build`, `/regression-gate`,
`/new-case`, `/doc-sweep`, `/commit`, `/perf-compare`, `/orient` and `/save-state`
carry the steps. This file carries only what must fire unprompted.

### 0. BEFORE ANYTHING ELSE, IF YOUR CONTEXT WAS JUST REBUILT

**R0) READ `PROGRESS.md`. IT IS THE ONE FILE.** Triggers: a compaction/summary
resume, `--continue`/`--resume`, a fresh session, or your own sense that you are
reaching for something and cannot recall whether it exists.

`PROGRESS.md` holds the live state — next action, what is in flight, what is
blocked, and the failed approaches not to retry — and points onward to
`docs/CAPABILITIES.md` (the inventory: what exists, the tooling, what was
refuted). `/orient` does this for you.

**Why, measured.** A compaction summary preserves the NARRATIVE of what was done and
LOSES THE INVENTORY of what exists. On 2026-09-10, after one compaction, I proposed a
capability already implemented and documented in three places (B4), and hand-rolled a
cold start around a nonexistent utility (`plasmaCreateFields`) while the generator and
`tools/make_mesh_arm.sh` sat in the tree with every trap already encoded. The user:
*"its like you dont remember anything of what weve done anymore... i cant follow and
continue like this to work with you."* The cost is not the minutes; it is that the user
cannot rely on continuity. **Keep both files TRUE** (D5) — they are only worth
re-reading if they are current.

### A. BEFORE YOU BELIEVE A RESULT

Every failure of the 2026-09-10 session lived here: three vacuous tests, three
inert knobs, five unequal-sample misreadings.

**A1) A MEASUREMENT IS NOT EVIDENCE UNTIL YOU HAVE SHOWN IT MEASURED WHAT YOU
THINK.** *(absorbs 2, 13, 15, 21, 23, 24, 42)* Six checks, each of which has
caught a wrong result on its own:

* **the tool RAN** (21) — a positive success marker plus exit status, never the
  absence of visible errors.
* **the knob MOVED** (42) — a sweep must prove the swept parameter actually
  changed, and the proof goes in the report. **Output identical to baseline in
  every digit is proof the knob is not wired — go read the code, do not re-run.**
* **the answer COULD have differed** (23) — exercise a count or metric on a case
  where it MUST come out different, or it is untested.
* **a DEFAULT is what the resolver returns** (24) — not a member initialiser, a
  header comment, or a tutorial.
* **the control is NAMED** (15) — every comparison row states what differs and
  what must match. If there is no valid control, say "no control". And report
  every magnitude against the physical scale it should be judged by.
* **the baseline was READ, not guessed** (13) — from `COMPARE.md`, never from
  directory names.
* and results get VALIDATED against data (2), asking for it when needed.

**A2) COMPARE ONLY AT EQUAL SAMPLES, EQUAL CONDITIONS, AND ON OUTPUTS.**
*(absorbs 43)*

* **equal samples** — five misreadings came from comparing runs at different step
  counts. Wait for a common endpoint; a partial window is not a verdict.
* **equal conditions** — refining a mesh at fixed dt RAISES the Courant number, so
  "10x cells doubled the iterations" measured the CFL, not the mesh; holding Co
  fixed reversed the conclusion. A **MODEL-DEPENDENT DEFAULT is an uncontrolled
  variable the moment the model differs** — diff both arms' resolved settings and
  start-up banners, not just the keys you changed.
* **on OUTPUTS** — a comparison reading the START snapshot, or 6-significant-figure
  log output, compares its own inputs or its own print precision.

**A3) DO NOT CONCLUDE AHEAD OF THE EVIDENCE; WHEN EVIDENCE CONFLICTS, RECONCILE
BEFORE EXPLAINING.** *(absorbs 11, 14, 20, 22)* No conclusions before the test is
run or the code read (11). **Contradiction protocol** (14): when a new result
contradicts an established one, FIRST reproduce the prior number in the same pass;
if it does not reproduce, say "I cannot reproduce X" and STOP — never state a
conclusion in the same message as an unreproduced prior. Never contradict a
documented conclusion from RECALL (20) — re-read and reconcile; that a parameter
still exists is not that the conclusion about it still holds. Never explain away a
discrepancy before its diagnostic is reconciled (22) — the diagnostic is the first
suspect, not the physics.

**A4) READ THE SOURCE BEFORE CITING IT.** *(absorbs 32, 35)* Never conclude from a
paper not read in full — an abstract or another paper's citation is a lossy
secondhand account; say so and ask for the PDF (32). When unsure of a single detail
INSIDE a paper the user has — an equation, a coefficient, an exponent — **ASK THEM
TO LOOK** rather than inferring it (35). → `literature-analyst`

**A5) MATCH THE DIAGNOSTIC TO THE FAILURE CLASS, AND CHECK INVARIANTS BEFORE
PHYSICS.** *(absorbs 27, 39)* Memory corruption → `-malloc_debug`/valgrind/ASan
(never guess from a `free()` backtrace); SEGV → `gdb -batch -ex run -ex "bt full"`;
SIGFPE → suspect a huge sentinel before a zero denominator; "which of N candidates"
→ change ONE and measure. **READ THE TOOL'S OWN HINT FIRST** (39). When a case
misbehaves check IN ORDER (27): (a) did every per-step update actually RUN — a
counter, not inspection; (b) is every DERIVED field consistent with its sources;
(c) is any field sitting on its TRIVIAL solution (vacuum field, floor, clamp,
initial value). Only then question the physics. → `diagnostician`

**A6) INSTRUMENT, DO NOT GUESS — AND DEBUG ON THE SMALLEST BED.** *(user, 2026-09-11)*
Measure where it fails, do not reason about where it might. A guess costs a full
build-run-read cycle and usually returns nothing; an instrument returns a fact every
time, **including when it stays silent** — a NaN locator that printed NOTHING through
six `DIVERGED_FNORM_NAN` failures proved the residual's COMPONENTS were finite and the
NORM had overflowed, which no amount of guessing produces. Measured on 2026-09-11:
four guesses at Newton's streamer failures, all four refuted, one cycle each; three
instruments, all three decisive on first use (gdb located a `0/0` in `fvc::interpolate`
in one attempt). **Bed size is part of the rule** — all of it was found on the ~2 s
coarse bed after hours lost at minutes-per-attempt on 1.15M cells. Before debugging
anything, ask whether a smaller bed reproduces it. → `diagnostician`, B6

**A7) STATE WHICH TIER OF EVIDENCE SUPPORTS A CLAIM, AND NEVER QUOTE A RESULT
ABOVE ITS TIER.** *(new)* Mathematical rigour here is a hierarchy, and the tiers are
not interchangeable:

| tier | what it proves | where |
|---|---|---|
| **1. Analytic unit bed** | the discretisation is right, against a closed form | `verification/fluxScheme*`, `testWallFlux`, `check_series_stack.py` |
| **2. Verification** | order of accuracy, mesh/time convergence, symmetry, conservation | order studies; `verification/` |
| **3. Validation** | the physics matches experiment or a published benchmark | Grubert 2009, the arXiv:2607.05137 streamer, Biagi/Phelps |
| **4. Plausibility** | it ran and nothing looked wrong | **not evidence.** Say so. |

A numerical claim needs tier 1 or 2; a physical claim needs tier 3. **A converged
run is tier 4.** State the observed order `p` with the grids it was fitted on, and
what was held fixed (A2). A unit bed with ground truth beats a CFD case for
answering a numerics question — prefer it, and say when no valid analytic reference
exists rather than implying one. → `numerical-analyst`, `physics-validator`

### B. BEFORE YOU TOUCH A RUNNING CASE OR WORKING CODE

**B1) SNAPSHOTS: WRITE AT LEAST TEN, RESTART FROM THEM, AND TO STOP A RUN FORCE A
SNAPSHOT RATHER THAN WAITING FOR ONE.** *(absorbs 41)* WRITE: `writeInterval <=
endTime/10`, set when the case is configured. USE: restart from a snapshot rather
than re-running from t=0 whenever it cannot change the conclusion — sound for a
defect that recurs from any nearby state, NOT when the conclusion depends on the
trajectory's history; and a restart REWRITES the snapshot it starts from, so copy
it read-only if reusing. TO STOP: `stopAt writeNow` in `controlDict` (it is
`runTimeModifiable`; then `touch system/controlDict`) **writes the state AND exits
in one action — it IS the kill.** Then confirm the written directory's name matches
the log's last reported time. To change something and KEEP RUNNING, lower
`writeInterval` live. A run's LATEST TIME DIRECTORY is the only thing that survives
being stopped. And check whether stopping is needed at all: `controlDict` is re-read
live, and a log eating the disk is fixed by `: > the.log`.

**B2) WHEN A PROPOSAL TURNS OUT LESS PROMISING THAN WHEN YOU PITCHED IT, STOP AND
ASK.** *(absorbs 44)* Trigger: the expected payoff changed after the work was
agreed. Mandatory when it risks the SOLVER (a wrong change corrupts silently) or is
a HUGE effort relative to what was described. Say what changed, the revised
cost/benefit, and the alternatives including doing nothing — then wait. Not a fait
accompli with a caveat attached.

**B3) BEFORE LAUNCHING: STATE THE PARAMETERS, CHECK THE DISK, BACKGROUND ANYTHING
OVER 1–2 MINUTES.** *(absorbs 9, 10, 17, 36)* STATE THE FULL PARAMETER SET first —
every case, arm and relaunch: circuit type and values, applied voltage/current and
its profile, `simulationType`, gas/pressure/temperature, `endTime`, timestep/Courant
strategy. A relaunch after a redesign needs it again (36). Check storage and forecast
it (9) — in WSL `df` is misleading, the real limit is the Windows C: drive.
Background anything over 1–2 min (10). Prefer the SHORTEST run that proves the point,
and ask whether a UNIT bed answers it better (17); state the discriminating
observable and the earliest time it is visible. → `/validate`

**B4) READ WHAT EXISTS BEFORE PROPOSING WORK; PROVE THE NEED BEFORE ADDING; AUDIT
EVERYWHERE AFTER FIXING.** *(absorbs 30, 37)*

* **BEFORE PROPOSING ANY FEATURE OR DECLARING A CAPABILITY GAP, READ WHAT THE SOLVER
  ALREADY DOES**: `docs/CAPABILITIES.md` (start here), `README.md`,
  `etc/boundaryRoles`, `docs/reference/`, `docs/design/*.md` (**check the `Status:`
  line — several are DESIGN with nothing implemented**). Never say "that would need
  building" without having looked. → `inventory-scout`
* **THIS APPLIES TO TOOLING, AND BEFORE DOING THE WORK, NOT ONLY BEFORE PROPOSING
  IT.** Before writing any script, launching any run, or setting up any case, read
  CAPABILITIES' tooling section, `tools/`, and the repo-root `*.sh`. **Each script's
  header records the trap it was written to avoid** — re-deriving the script means
  re-earning the traps.
* **NEVER HAND-EDIT A GENERATED FILE.** The value goes in `configuration/config` and
  the dictionary references `$key`; the boundary goes in `configuration/boundaries`.
  Then prove the key is not dangling (command card). A hook blocks the generated
  files it can identify, but it only knows what `.gitignore` records.
* Before adding a mechanism the framework already has, prove BY MEASUREMENT that the
  framework's own fails in the configuration you would ship — the correct response to
  a framework mechanism failing on bad input is a GUARD, not a parallel mechanism (30).
* **A bug is a property of a MECHANISM, not of the call site that exposed it** (37):
  search every other use of the same mechanism and fix those too before calling the
  fix complete.

**B5) NO SILENT REGRESSIONS. RUN THE GATE, AND IF A BASELINE MOVED AND YOU CANNOT
SAY WHY, STOP AND ASK.** *(new)* Before any commit touching solver source, run
`/regression-gate`. A changed number is one of three things and they are not
interchangeable: a **REGRESSION** (revert or fix), an **INTENDED IMPROVEMENT** (update
the baseline IN THE SAME COMMIT, with the reason and the date), or a **STALE
BASELINE** (say so explicitly). **If you cannot tell which, you do not get to guess —
ask.** Never update a baseline to make a test pass; that converts a failure into a
silent, permanent wrong answer. A test whose baseline is refreshed without a stated
reason has stopped being a test.

**B6) EVERY TEST HAS A `--fast` MODE, AND THE FAST MODE IS WHERE YOU DEBUG.** *(new)*
Every bed, sweep and suite must have a documented seconds-scale variant, and the
command card carries it alongside the full form. This generalises A6's measured
finding: every Newton defect of 2026-09-11 was found on the coarse bed at ~2 s/run
after hours wasted at minutes/attempt. Add the fast mode when you add the test, not
later. **A fast mode is for DEBUGGING and must never be quoted for physics** — the
coarse streamer bed does not resolve a streamer.

### C. CODE STRUCTURE AND THE FRAMEWORK

**C1) IN SoPLASMA, WORK INSIDE OPENFOAM'S OWN STRUCTURES — NEVER BUILD A PARALLEL
WORLD BESIDE THEM.** *(absorbs 29, 31, 40)* **Scoped to SoPlasma; does NOT apply to
SoEEDF**, which must keep building with no OpenFOAM dependency.

* Express physics with `fvm::`/`fvc::`, carry data in `GeometricField`s with real
  `dimensionSet`s, read settings through `dictionary` and runtime selection, use
  `Info`/`FatalErrorInFunction`, follow existing naming (G3). Confine a third-party
  library behind a thin MECHANICAL boundary that only translates — `snesBridge.H` is
  the worked example (40).
* **A control object's exit path is part of its CONTRACT** (29): never leave its loop
  with `break`/`return`/`goto`; set a flag and DRAIN the loop. `pimple.loop()` promises
  a `finalIter`, and a raw `break` voids it — that froze Poisson's source for 99.96% of
  steps. After any loop-control change, verify per-step INVARIANT COUNTS, not that it runs.
* **A collective must never sit downstream of a guard on a LOCAL quantity** (31):
  `reduce`/`gSum`/`gMax`/`returnReduce` must be reached by every rank the same number of
  times. Hoist the collective out rather than guarding it harder, and **test at least TWO
  decompositions** — the same bug deadlocks under one and raises `MPI_ERR_TRUNCATE` under
  another, and one variant faked a PHYSICS failure.
* **`fvMatrix::residual()` is BROKEN in parallel** — misreports by ~21 orders of
  magnitude, independent of PETSc. Compute any outer/Newton residual by explicit `fvc::`.
* **All solution state must be idempotent within a step**; a new field MUST join the
  discard path. → `code-quality-reviewer`, `openfoam-implementer`

**C2) ADVANCED NUMERICS, CPU EFFICIENCY, CORRECT PHYSICS.** *(absorbs 4, 7)* Restated
from the objectives above, which is where aims belong.

**C3) A PERFORMANCE CLAIM STATES ITS SCALING AND IS MEASURED AT TWO PROBLEM SIZES AT
MATCHED COURANT.** *(new; promoted out of A2, where it was easy to miss)* The test is
mesh-independence of the **ITERATION COUNT**, not wall clock. Measure **per unit of
simulated time at equal accuracy, never per step** — Newton costs ~4–5x per step, so a
per-step comparison is rigged for Picard. An optimisation that helps at 2,000 cells and
not at 20,000 is not an optimisation; three were refuted exactly that way (ILU(1),
`selfp`, per-cell scaling). → `performance-analyst`

### D. RECORDS THAT OUTLIVE THE CONVERSATION

Nothing important may exist only in chat: it does not survive compaction.

**D1) DOCUMENTATION IS PART OF "DONE" — FOR THE USER.** *(absorbs 1, 12)* The WHOLE
surface: module `README.md`, `docs/reference/` for any dictionary key, the
`docs/design/` note, any report the result bears on. A key rename or API change needs a
SWEEP in the same change; anything left behind is not stale, it is an INSTRUCTION a
later session will follow. Every comparison case carries a `COMPARE.md` in its own
directory, written WHEN THE CASE IS CREATED: the question, the baseline's ABSOLUTE
path, the reference numbers, the extraction command. → `/doc-sweep`, `docs-sweeper`

**D2) EVERY MEASUREMENT CARRIES ITS DATE; SUPERSEDED TEXT IS MARKED, NEVER DELETED.**
*(absorbs 18, 19)* A measurement without a date cannot be told from a live conclusion
and will be read as one. A superseded result is marked `SUPERSEDED BY <what>` in place,
stating what changed and what still holds — a finding that quietly disappears leaves no
way to recognise its stale copies elsewhere.

**D3) EVERY CASE WRITES A MACHINE-READABLE TIME SERIES, AT PROBES.** *(absorbs 6, 26)*
One to three probe points per case, each recording versus time: every species density,
`chargeDensity`, `Emag`, the reduced field, electron energy, currents; plus run-level
`dt` **AND WHICH LIMITER SET IT**, outer iterations, rejected steps. Probes rather than
fields (a field history is unwritable in 3-D); keep min/max too, but a reduction throws
away WHERE and where is usually the answer. A DERIVED quantity sits NEXT TO its sources
in the same file. **If diagnosing a case requires parsing the log, the monitoring is the
defect.** For validation cases, produce publication-ready figures AND the python that
generates them from the data files.

**D4) TELL ME WHEN THERE IS UNCOMMITTED WORK — DO NOT WAIT TO BE ASKED.** *(rule 38,
verbatim in force)* Check `git status` and say so PLAINLY as its own statement: after
any new file, after a logically complete unit of work, before anything long-running,
before any destructive git operation. A PROMPT, not licence to commit unasked. **NEW
UNTRACKED FILES ARE THE URGENT CASE** — they appear in no diff and are destroyed by a
clean. **STAGE EXPLICITLY: never `git add -A`**, never `git add <dir>` on a tree holding
case output; re-read `git diff --cached --name-only` and unstage what does not belong.

**D5) DOCUMENTATION FOR THE NEXT SESSION — UPDATE `PROGRESS.md` AFTER EVERY MEANINGFUL
UNIT OF WORK.** *(new)* D1 is for the eventual USER; this is for the next session, and
the register is different — what was tried, what the number was, what is still unknown.
Never marketing prose. The contract:

* **Check off completed items WITH THE DATE.**
* **Note what worked, what didn't, and what is blocked** — and on what.
* **Record failed approaches so they are not re-attempted**, in the form
  *"Tried X — didn't work because Y. Switched to Z."* This section is the one that pays
  for the file.
* **Add new tasks discovered during implementation**, at the moment you discover them.
* **When stuck, keep a running log of attempts** in PROGRESS.md — dated, so the next
  session inherits the search rather than repeating it.

`/save-state` does this and prompts D4. **"Save the state" before ending a session.**

**D6) SMALL, TESTABLE COMMITS.** *(new)* One logical unit per commit, explicitly staged,
with a message saying WHAT CHANGED and WHY — and the measured number where there is one.
Do not mix a refactor with a behaviour change: when the behaviour then moves, nothing
says which half did it. If a change cannot be described in one sentence without "and",
it is two commits. → `/commit`

### E. WORKING WITH THE USER

**E1) SAY WHAT YOU ARE ABOUT TO DO, BEFORE DOING IT.** *(absorbs 16, 28)* Before any
piece of work state plainly and up front: what you are about to do, what question it
answers, what the steps are. Not narrated as you go, not inferred from tool calls, not
explained afterwards. **It applies EVEN WHEN the work was requested** — "verify Part A"
authorises the work, it does not say what verifying consists of. If the plan is more
than one step, list the steps first. At minimum, one line naming the question before any
stretch of more than ~3 tool calls whose purpose is not obvious.

**E2) PLAN FIRST; EXECUTE IN SMALL VALIDATED INCREMENTS AFTER APPROVAL.** *(absorbs 3, 5)*
**One change at a time.** Changing the timestep AND the table range together made
neither result interpretable.

**E3) NEVER SWITCH THREADS WITHOUT ASKING.** *(rule 25)* Finish the open thread, or say
which you are leaving and which you are starting, and WAIT. An open-ended "ok",
"continue" or "next?" means CONTINUE THE CURRENT WORK — never licence to pick a different
task, however important elsewhere. The test: can I name the thread I am on, and did the
user's last message point at THAT thread?

**E4) RESEARCH IS A VALID OUTCOME, NOT A LOWER-PRIORITY ONE.** *(rule 33)* The user is a
professor; publishing is an equally important track. Never treat "this needs genuine
research effort" as grounds to defer. Assess an idea on whether the physics is sound,
whether a genuine gap exists (checked, not assumed), and whether it is tractable.

**E5) DELEGATE TO THE SPECIALIST AGENTS.** *(new)* They exist so the expensive judgement
runs at the right effort and the mechanical work runs cheap. Each carries the traps for
its area, so delegating is also how those traps get applied rather than remembered.

| agent | delegate when |
|---|---|
| `numerical-analyst` | designing or judging a discretisation, preconditioner, time integrator or order study |
| `physics-validator` | a physical claim, a literature comparison, units/vocabulary, "is this result physical" |
| `diagnostician` | a crash, a NaN, a stall, a wrong number with no obvious cause (A5/A6) |
| `literature-analyst` | a paper must actually be read before being cited (A4) |
| `openfoam-implementer` | writing solver/BC/model code inside OpenFOAM's structures |
| `test-quality-auditor` | before believing a test — does it hunt vacuity (A1) |
| `code-quality-reviewer` | before committing solver source — C1 structural violations |
| `performance-analyst` | any speed/cost claim (C3) |
| `docs-sweeper` | a key rename or API change needs the whole D1 surface swept |
| `inventory-scout` | "does the solver already do X?" — always, before proposing (B4) |

### Old-number lookup — all 44 still resolve

| old | now | old | now | old | now | old | now |
|---|---|---|---|---|---|---|---|
| 1 | D1 | 12 | D1 | 23 | A1 | 34 | *ref doc* |
| 2 | A1 | 13 | A1 | 24 | A1 | 35 | A4 |
| 3 | E2 | 14 | A3 | 25 | E3 | 36 | B3 |
| 4 | C2 | 15 | A1 | 26 | D3 | 37 | B4 |
| 5 | E2 | 16 | E1 | 27 | A5 | 38 | D4 |
| 6 | D3 | 17 | B3 | 28 | E1 | 39 | A5 |
| 7 | C2 | 18 | D2 | 29 | C1 | 40 | C1 |
| 8 | *dropped* | 19 | D2 | 30 | B4 | 41 | B1 |
| 9 | B3 | 20 | A3 | 31 | C1 | 42 | A1 |
| 10 | B3 | 21 | A1 | 32 | A4 | 43 | A2/C3 |
| 11 | A3 | 22 | A3 | 33 | E4 | 44 | B2 |

**Dropped, deliberately:** rule 8 (save sessions for a vscode restart) — a harness
concern, governed no decision. Rule 34 (consult JC-PIC before deriving plasma theory)
— a POINTER, moved to `Literature/reference-codes-and-manuals/JC-PIC_Boeuf/README.md`
where it will be found. Rule 17's 5-minute target — routinely impossible, and an unmet
rule teaches that rules are optional.

---

## THE COMMAND CARD

Every command was executed or existence-checked on 2026-09-11; **a command that cannot
be run is not listed**, because a wrong command here is an instruction a later session
will follow. Each skill carries the traps and the `--fast` variants for its area — this
is the short form.

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
pkill -9 soPlasmaFoam; pgrep -c soPlasmaFoam        # never `pkill -f`: misses -parallel,
                                                    # and can kill the calling shell
NCELL=130 NREFINE=1 CASE=$HOME/streamer-fast ./make-smoke-case.sh   # the ~2 s coarse bed
CASE=$HOME/streamer-fast ./rerun.sh                                # (A6/B6; never for physics)

# CASE SETUP                                                              -> /new-case
./Allrun-serial                                     # PASS: log contains '=== done'
# Prove no configuration/config key is dangling — an inert $key once invalidated a whole
# arm that exited 0. Use grep -rF '$'"$k"; inside double quotes "\$$k" expands to the PID.
for k in $(grep -oE '^[a-zA-Z][a-zA-Z0-9_]*' configuration/config | sort -u); do
  grep -rqF '$'"$k" system/ constant/ configuration/ Allrun-serial 2>/dev/null || echo "DANGLING: $k"
done

# SoEEDF — standalone, no OpenFOAM dependency (C1 does not apply there).
cd $HOME/Projects/SoEEDF && cmake -S . -B build && cmake --build build -j
```

**Two build traps that bite outside `/build`'s moment.** `grep FAIL buildlogs/*` CANNOT
work — `OK`/`FAIL` go to build-all.sh's *stdout*, and the repo's `build.log` is a stale
Aug-10 artefact that will report a month-old build as current. And `wmake src/numerics`
silently does NOT build `libplasmaNewtonSolverPETSc`; it is a separate target.
