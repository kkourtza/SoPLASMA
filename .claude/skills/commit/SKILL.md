---
name: commit
description: Commit work in this tree the way rules D6 + D4 require — one logical unit, staged explicitly by pathspec (never `git add -A`, a hook blocks it), the staged list re-read and cleaned of case output, the regression gate run when solver source is touched, and a message saying what changed and why with the measured number. Also reports the unpushed backlog.
---

# /commit — one logical unit, explicitly staged

**The question this answers:** what is uncommitted, which ONE unit goes in *this* commit,
and is exactly that unit staged — no run output, no second unrelated change.

**Do NOT use this** to survey the tree at session start (`/orient`), to write the narrative
state (`/save-state`, run it after), or as licence to commit unasked — D4 is a PROMPT to
tell the user, not permission. Pushing is a separate decision: step 8, ask first.

## Repo facts (measured 2026-09-11, re-check with the commands below)

* Live tree `/home/kkourtza/soplasma-scratch`, branch
  `feature/chemistry_integration_and_BoltzmannSolver`, remote `origin` =
  `git@github.com:kkourtza/SoPLASMA.git`. **205 commits unpushed** — the whole Newton/JFNK
  solver exists only on this machine plus those local commits. `upstream` is rpasolari's
  fork and its pushurl is literally `DISABLED`; never push there.
* SoEEDF `/home/kkourtza/Projects/SoEEDF`, branch `fix/mechc-audit-unitarity-bound`,
  **NO UPSTREAM AT ALL** — never pushed, 52 commits ahead of `origin/master`.
* `validation/`, `verification/` and `tutorials/` run output became gitignored on
  2026-09-11: `git status --porcelain` went **1310 lines / 182 s → 3 lines / 0.00 s**
  (`core.untrackedCache` + `feature.manyFiles` are set). Status is short and trustworthy
  again. If it ever prints hundreds of lines, a NEW output family escaped `.gitignore` —
  fix `.gitignore`, do not stage around it.
* `.claude/hooks/guard-bash.sh` **blocks** `git add -A`, `git add --all`, `git add .`
  (PreToolUse, not advisory).

## Procedure

**1. Look at the whole tree.**
```bash
cd /home/kkourtza/soplasma-scratch && git status --porcelain
```
PASS: returns in under a second with a handful of lines. **Read the `??` lines FIRST and
say them to the user (D4): NEW UNTRACKED FILES ARE THE URGENT CASE** — they appear in no
diff and a `git clean` destroys them. (Measured: `ROUNDW.C/.H` were built into
`libROUNDSchemes` for four days while untracked; `verification/` — the only analytic
ground truth in the tree — was untracked *in its entirety*.)
If it exceeds ~50 lines or ~10 s, stop and diagnose: `git status --porcelain -uno` shows
tracked changes only (0.58 s even in the old 1310-line state).

**2. Name the ONE unit — from the diff, never from filenames.**
```bash
git diff --stat && git diff -- <path>            # unstaged, tracked
git status --porcelain | grep '^??'              # untracked: decide TRACK / IGNORE / DELETE
git check-ignore -v <path>                       # why is this path (not) ignored
```
The D6 test: **can you state what changed and why in one sentence with no "and"?** If not,
it is two commits. **Never mix a refactor with a behaviour change** — when the behaviour
later moves, nothing says which half did it. Generated output that is untracked belongs in
`.gitignore` (mirror the `smoke/` stanza's reasoning), not in the commit.

**3. Stage EXPLICITLY, by pathspec (D4).**
```bash
git add <path> [<path> ...]
```
Never `git add -A` / `--all` / `.` — the hook refuses them. Never `git add <dir>` on a tree
holding case output: measured, a broad `git add tutorials/plasma` pulled **436 paths** of
time directories, logs and postProcessing alongside 10 real dictionary changes. A directory
pathspec is acceptable only for pure-source trees (`src/`, `applications/`, `tools/`,
`docs/`) — and step 4 still runs.

**4. Re-read the staged list and unstage what does not belong. THIS IS THE STEP THAT GETS SKIPPED.**
```bash
git diff --cached --name-only
git diff --cached --name-only | git check-ignore --stdin          # PASS: prints nothing
git diff --cached --name-only | grep -E '(^|/)(logs|postProcessing|polyMesh|plasmaTables|ionTables|VTK|processor[0-9]+)/|\.log$|(^|/)[0-9]+(\.[0-9]+)?(e-?[0-9]+)?/|\.snapshot'
git diff --cached --name-only | wc -l
```
PASS: check-ignore prints nothing (exit 1); the count matches the unit you named in step 2;
and every hit of the output-shape grep falls in one of the three families below. **That grep
already matches 488 TRACKED paths** (measured 2026-09-11 over `git ls-files`), so a hit is a
question, not a verdict — unstaging on sight deletes ground truth from the commit:
* `tutorials/**/constant/plasmaTables/` — 410 paths, three cases (`positiveStreamer_LMEA_fast`,
  `positiveStreamer_LMEA_minimal`, `positiveStreamer_fixedMesh`). Shipped deliberately: a
  tutorial `Allrun-serial` does not regenerate the tables, so a fresh clone could not start them.
* `**/constant/polyMesh/` — 54 paths: `verification/testSnesJFNK{,2Field}/` plus the vendored
  `ThirdParty/{blastAMR,libROUNDSchemes}` meshes. (`verification/fluxScheme{1D,2Dnonortho}/`
  polyMesh IS ignored, `.gitignore:190-191` — those two beds generate their mesh.)
* `<case>/0/` INITIAL CONDITIONS, not time output — 24 paths:
  `verification/{fluxScheme1D,fluxScheme2Dnonortho,testSnesJFNK,testSnesJFNK2Field}/0/` plus
  `ThirdParty/{libROUNDSchemes,petsc4Foam}` tutorials. `verification/fluxScheme1D/0/` is an
  INPUT to the step-5 gate, and `verification/` is the tree's only analytic ground truth — it
  was lost once already to being untracked.
Any hit outside those three is output. Unstage it with `git restore --staged <path>` (git 2.43 here).

**5. Solver source staged → run the gate (B5).**
Trigger, mechanically: `git diff --cached --name-only | grep -E '^(src|applications)/'`
returns anything. Then run `/regression-gate`. Acceptance (the command card's strings):
```bash
./check-no-running-solvers.sh                    # exit 0 required (FORCE=1 overrides)
tools/run_electrostatics_tests.sh                # "N ok, 0 failed"; exit 2 = env unsourced
testWallFlux && testWallLoss && testVibRelax && testCoulombHeating && testEmission
verification/fluxScheme1D/Allrun                 # results.txt: EXACTLY 504 RESULT lines
```
`&&`, not `;` — otherwise a failing bed is silently skipped. Compare `results.txt` to
`results.baseline.txt` by RELATIVE TOLERANCE, never `diff`. `testAitken` has no pass/fail
machinery and ALWAYS exits 0: it is never a pass.
**--fast (B6):** the five mesh-free beds above are the seconds-scale variant and are where
you iterate; the full gate still runs before the commit lands.
**A moved baseline is one of three things** and you do not get to guess: REGRESSION (fix or
revert), INTENDED IMPROVEMENT (update the baseline *in this same commit*, with the reason
and the date), or STALE BASELINE (say so explicitly). If you cannot tell which — **stop and
ask**. Never update a baseline to make a test pass.

**6. Write the message: WHAT changed, WHY, and the measured number.**
Repo style (see `6505eb4`, `0f88693`, `f17a489`): a subject that states the fact, an
optional `scope:` prefix, ~72 chars; then the reason, the number with its date (D2), the
rule ID the change implements, and what verification was actually run. `-F -` reads the
message from stdin, so the blank lines and the footer survive.
```bash
git commit -F - <<'EOF'
<scope>: <what changed, one sentence, no "and">

<why it was wrong / what it buys, with the measured number and the date it was
measured. Cite the rule ID (D6, B5, C1, ...) where the change implements one.
Name the verification run and its result.>

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```
If the change renames a key or alters an API, the D1 documentation sweep (`/doc-sweep`)
belongs in the SAME commit — text left behind is not stale, it is an instruction a later
session will follow.

**7. Verify what landed.**
```bash
git show --stat --oneline HEAD | head -40
git status --porcelain
```
PASS: every path of the unit is in the commit, and nothing belonging to it is still dirty.

**8. State the backup position; offer the push.**
```bash
git rev-list --left-right --count origin/feature/chemistry_integration_and_BoltzmannSolver...HEAD
#   prints "<behind>\t<ahead>"; 0<TAB>205 on 2026-09-11
git push origin feature/chemistry_integration_and_BoltzmannSolver      # origin ONLY
# SoEEDF, which has no upstream at all:
cd /home/kkourtza/Projects/SoEEDF && git rev-list --count HEAD ^origin/master   # 52
git push -u origin fix/mechc-audit-unitarity-bound                     # -u is required
```
Say the unpushed count plainly whenever a unit of work completes. **A commit is not a
backup:** measured 2026-09-09, two days without one left ~6,900 lines of new solver source
untracked with no backup of any kind, and the user had to catch it. Ask before pushing.

**9. Then `/save-state` (D5)** — PROGRESS.md §4 gets the ledger tick with the date, §5 gets
any approach that failed, in the form *"Tried X — didn't work because Y. Switched to Z."*
