---
name: orient
description: Rule R0 made executable. Run at session start, after a compaction/resume or --continue, or whenever you are reaching for something and cannot recall whether it exists. Reports in ~6 lines where the open thread is, whether anything is running, whether run output is NEWER than the last commit (a measurement nobody has read), the git and backup state, and the next concrete action.
allowed-tools: Read, Grep, Bash(git log:*), Bash(git status:*), Bash(git rev-list:*), Bash(pgrep:*), Bash(ls:*), Bash(stat:*), Bash(date:*), Bash(grep:*), Bash(tail:*), Bash(wc:*), Bash(sed:*)
---

# orient — R0 as a procedure

**The question:** *what thread am I on, is anything running, is there a result nobody has
read, and what is the next concrete action?*

**It REPORTS, it does not summarise.** The output is the ~6-line block in step 7. Never
paraphrase PROGRESS.md at length and never restate `docs/CAPABILITIES.md` — a third copy
of their contents is the G1 defect ("a setting stated in two places is a defect even when
both copies agree"). **NOT for:** mid-thread with intact context (use `--fast`, step 8);
as a substitute for `/save-state` (this skill only READS — it never writes PROGRESS.md);
as a build or regression check (`/build`, `/regression-gate`). All commands run in
`/home/kkourtza/soplasma-scratch` unless stated.

## 1. Read the live state (R0)

Read `PROGRESS.md` §1 WHERE WE ARE RIGHT NOW, §2 IN FLIGHT, §3 BLOCKED. Stop there. §5
FAILED APPROACHES is read only when about to propose something (B4); §4 TASK LEDGER only
once the thread in §1 is settled.

PASS: you can name the open thread in one sentence and quote §1's **NEXT CONCRETE
ACTION** verbatim. If not, PROGRESS.md is stale — say so, do not reconstruct it from
memory.

## 2. Is anything running?

```bash
pgrep -x soPlasmaFoam      # rc=1, no output = nothing running
pgrep -x mpirun            # rc=1, no output = nothing running
```

Each as its OWN simple command. Measured in this harness 2026-09-11: `pgrep -af
soPlasmaFoam` returned rc=0 with three hits and nothing running — it matched agent shells
whose command lines merely contain the string, and `ps -eo … | grep soPlasmaFoam` fails
the same way. `-x` matches the process NAME and is clean (verified rc=1), but is blind to
names over procps' 15-character `comm` limit: for `plasmaChemistry0D` (17),
`plasmaCreateSpeciesFields` (25), `testChemistryBackends` (21) use `pgrep -af
"[p]lasmaChemistry0D"`. Never `pkill -f` — it kills the calling shell. If something IS
running, reconcile it against its snapshot before any kill (B1): the latest time
directory is all that survives.

## 3. Is there a result nobody has read?

```bash
C=$(git log -1 --format=%ct)
N=$(ls -t validation/*.log validation/*/log.* 2>/dev/null | head -1)
M=$(stat -c %Y "$N")
echo "last commit $(date -d @$C '+%F %H:%M')   newest output $(date -d @$M '+%F %H:%M')  $N"
[ "$M" -gt "$C" ] && echo "UNREAD RESULT: $N is newer than the last commit"
```

Runtime 0.01 s. **Output newer than the last commit means a measurement exists with no
commit, no COMPARE.md verdict and no memory** — the exact state at 2026-09-11 18:47
against commit 04f815c at 17:50, where ~1 h of Newton measurement sat unread and two arms
had been cut off mid-run. Read it before launching anything: a completed measurement
silently redone is the cost R0 exists to prevent.

## 4. If step 3 fired, read that arm properly — three checks, not one

```bash
f=validation/<arm>.log
echo "steps: $(grep -c '^Time = ' $f)"
grep -o 'due to [A-Z_]*' $f | sort | uniq -c | sort -rn   # A1: never the bare failure count
tail -2 $f                                                # truncation detector
grep -c 'Picard warm-up COMPLETE' $f                      # A1: did it actually run Newton?
```

* **Reason histogram, never the bare count.** `DIVERGED_LINEAR_SOLVE` dominating means
  the KSP budget expired, not a solver failure — that misreading made 42/398 vs 4/398
  look like solver quality when it was `-ksp_max_it 100` expiring 4–15 iterations short.
* **Truncation.** A log whose last line is mid-iteration (e.g. `Discretizing transport
  with standard schemes...`) with no live process is TRUNCATED, not a result. Verified on
  `validation/lad_newton_dt2e-11_k1000.log`: 7 steps, ends mid-step.
* **Newton actually ran.** `0` hits on the handover banner means the arm ran Picard and
  every number in it is void — its absence once meant 21,174 steps silently ran Picard.

## 5. Git and backup state (D4)

```bash
git status --porcelain -uno                                  # 0.58 s / ~193 lines (2026-09-11)
git status --porcelain src docs tools etc build-all.sh tutorials   # 0.00 s; 0 lines = source clean
git rev-list --left-right --count origin/feature/chemistry_integration_and_BoltzmannSolver...HEAD
git -C /home/kkourtza/Projects/SoEEDF status -sb | head -1   # no '...origin/' = NO UPSTREAM
```

**Bare `git status` is now cheap and trustworthy again: 1 line, 0.6 s.** It used to cost
**181.97 s and print 1310 lines** — exceeding the 120 s Bash timeout, so it read as a hung
tool — because `validation/` run output was tracked and 1109 of those lines were experiment
artefacts. Fixed 2026-09-11 (`.gitignore` + 2742 files removed from the index +
`core.untrackedCache`). **If it is ever slow or long again, that is a REGRESSION in the
ignore rules, not something to work around** — the scoped forms above then tell you which
path reintroduced it.

`rev-list` measured `0	205` on 2026-09-11; SoEEDF's branch
`fix/mechc-audit-unitarity-bound` prints no upstream at all, so that work exists only on
this machine. Report uncommitted work plainly as its own statement (D4) — do not commit
unasked.

## 6. The memory index — ordering is NOT filename ordering

```bash
ls -t /home/kkourtza/.claude/projects/-home-kkourtza-Projects-SoEEDF/memory/*.md | head -6
```

Sort by mtime, never by the date in the filename. The newest state of the CURRENT
(Newton) thread is not a `session-state-*` file at all — it is `ksp-cap-was-the-bug.md`
and `newton-vs-picard-benchmark-state.md`. Reading `session-state-2026-09-08.md` as
"latest" gives the Grubert thread and misses the ksp-cap invalidation, the single fact
that changes what to do next. (Within Grubert: `-09-07c` supersedes `-09-07b` supersedes
`-09-07` supersedes `-09-06`; `-09-08` continues `-09-07c`.) Read at most the top TWO, and
only if PROGRESS.md left something unexplained; reading ~129 memories is not orientation.

**There is only ONE memory store.** Since 2026-09-11 the `soplasma-scratch` project memory
directory is a SYMLINK to the `Projects-SoEEDF` one, so both roots see the same ~129 files
and the base cannot fork by starting a session from the other directory. The old pointer
note that used to live on the soplasma side is gone — if you ever see two different
`MEMORY.md` files, the symlink has been broken and that is the bug.

## 7. Report, then stop

Emit exactly this shape and nothing more:

```
RUNNING: nothing | <pid> <case> <elapsed>
UNREAD:  none | <arm>.log <hh:mm> > commit <hash> <hh:mm>
GIT:     <n> dirty (-uno), <m> in src/docs/tools/etc; <k> unpushed
THREAD:  <PROGRESS §1, one sentence>
NEXT:    <PROGRESS §1 NEXT CONCRETE ACTION, verbatim>
BLOCKED: <only the §3 rows that touch NEXT>
```

Then WAIT (E1). An "ok" / "continue" / "next?" after this means the thread in THREAD —
never licence to pick a different task (E3).

For what EXISTS — inventory, tooling, the full refuted list — point at
`docs/CAPABILITIES.md` (§4/4b/4c refuted, §6 open); read it when about to propose a
feature or write a script (B4), not during orientation. Never orient from
`../Projects/SoEEDF/docs/framework-state.md`: dated 2026-09-04, zero mentions of
newton/snes/jfnk against 88 in CAPABILITIES.md — it predates the whole Newton solver.

## 8. `--fast` (B6)

Steps 2, 3 and the first two commands of 5: four commands, under 1 s, no file reads. Use
mid-session for "is anything running, and is anything unread". **The full form, including
step 1, is mandatory after a compaction, a `--resume`/`--continue`, or a fresh session** —
the fast form answers the machine's state, not the thread's.
