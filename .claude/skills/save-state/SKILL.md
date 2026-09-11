---
name: save-state
description: Rule D5 made executable — update PROGRESS.md after a meaningful unit of work, then report uncommitted work (D4), what the hooks caught, and whether a memory is owed. Run before ending a session, before any long run, and whenever a thread changes state. It never stages and never commits.
allowed-tools: Read, Edit, Write, Bash
---

# /save-state — write down what the next session needs

**The question this answers:** if this session ended now, what would the next one have
to re-derive? D5: `PROGRESS.md` is documentation *for the next session*, not for the
user (that is D1). The register is *what was tried, what the number was, what is still
unknown*. Never marketing prose. Every entry carries its date (D2). **Run it** before
ending a session, after a logically complete unit of work, before launching anything
long (B3), and the moment a thread changes state (blocked → unblocked, refuted).
**Do NOT use it** to commit — it stages nothing (that is `/commit`, D6); to record a
hypothesis you have not measured — a §5 entry needs the *because Y*; or to mark an
in-flight arm done — an unread run belongs in §2 with its stopping state, never in §4.
**`--fast`** (B6): steps 1, 2a, 2d and 3 only, between units of work mid-session; the
full form is mandatory before a session ends.

## 1. Get the facts first — do not write from recall (A3)

```bash
cd /home/kkourtza/soplasma-scratch
date +%Y-%m-%d
git log -1 --format='last commit: %h %ad' --date=format:'%Y-%m-%d %H:%M'
ps -eo pid,etime,cmd | grep -E 'soPlasmaFoam|mpirun' | grep -v grep   # empty = nothing running
find validation -maxdepth 2 -newermt "@$(git log -1 --format=%ct)" \
     \( -name '*.log' -o -name 'log.*' \)                             # 0.01–0.2 s
```

PASS: you can state the date, the last commit, what is running, and **every file the
last line printed**. Each of those is a measurement nobody has written up — §2 exists
because ~1 h of it sat unread on 2026-09-11 (five arms, no commit, no `COMPARE.md`).
Do not widen that `find`: `find validation tutorials -maxdepth 3` measured **72 s**.

## 2. Update every PROGRESS.md section that changed — edit in place, never rewrite it

* **a. The `Last updated:` stamp and §1 WHERE WE ARE.** `.claude/hooks/session-start.sh`
  greps that stamp and echoes §1's **first 6 non-blank lines** into the next session's
  context, so the **NEXT CONCRETE ACTION** line must fall inside those 6; anything below
  is silently truncated at the moment it matters most.
* **b. §2 IN FLIGHT.** From step 1's `ps` and `find`. Per arm: the case path and *the
  state it stopped in* — `t=`, steps, the SNES reason, whether it was cut off mid-run.
  A completed-but-unread measurement is the expensive thing here.
* **c. §3 BLOCKED.** Three columns: item | blocked on | **what unblocks it**. If you
  cannot fill the third, the row is a complaint, not a blocker.
* **d. §4 TASK LEDGER — check off completed items WITH THE DATE.** `- [x] 2026-09-11 —
  <what> (<commit hash>)`. Carry the measured number into the line itself (`100 →
  42/398 failures (10.6%); 1000 → 4/398 (1.0%), +2.5% wall clock`), not a pointer to
  it. Add tasks discovered during the work *at the moment you discover them*.

* **e. §5 FAILED APPROACHES — the section that pays for the file.** Exact form:
  **`Tried X — didn't work because Y. Switched to Z.`** Y must be a measurement with
  its bed and size attached — "helped at 2k (−33% steps), *hurt* at 20k (+3.5%)" is an
  entry, "seemed slower" is not (C3: a result at one problem size does not transfer).
  An idea the user declined goes in the **NOT ON THE LIST BY USER DECISION** block, not
  among the refuted; re-raising the first as merely unexplored is an E3 failure.
* **f. §6 DISCOVERED TASKS** — one line each, the date found and the next action.
* **g. §7 RUNNING NOTES**, the dated attempt log kept while stuck. **Clear resolved
  entries into §4 (it worked) or §5 (it didn't) and delete them from §7**; one left
  there reads as still-open. If §7 is non-empty the thread is open and §1 must say so.

Anywhere in the file, superseded text is marked `SUPERSEDED BY <what>` in place, never
deleted (D2) — a finding that quietly vanishes leaves no way to spot its stale copies.

## 3. Report uncommitted work plainly — UNTRACKED FIRST (D4)

Its own statement, not an aside. A prompt, not licence to commit.
```bash
cd /home/kkourtza/soplasma-scratch
git status --porcelain | grep '^??'            # UNTRACKED — the urgent case
git status --porcelain -uno                    # tracked changes (0.08 s, 1 line today)
git rev-list --count @{u}..HEAD 2>/dev/null || echo 'NO UPSTREAM — never pushed'
cd /home/kkourtza/Projects/SoEEDF && git status --porcelain && \
git rev-list --count @{u}..HEAD 2>/dev/null || echo 'NO UPSTREAM — never pushed'
```

* **Untracked first** — they appear in no diff and a `git clean` destroys them. Measured
  2026-09-11 22:00: `?? PROGRESS.md` and `?? .claude/`, i.e. the orientation file and
  the whole hook/skill set with **no backup of any kind**; the same shape cost two days
  and ~6,900 lines of untracked Newton solver on 2026-09-09.
* **Backup state**, stated every time: 205 unpushed commits here (2026-09-11); SoEEDF's
  branch `fix/mechc-audit-unitarity-bound` has **no upstream at all**.
* **If bare `git status` takes minutes again, something got un-ignored** — it cost
  181.97 s / 1310 lines until `validation/` output was ignored on 2026-09-11 and is
  0.01 s / 3 lines now. The timing *is* the detector. Never `git add -A` or
  `git add <dir>` here; a hook blocks both (D4).

## 4. Report what the hooks caught, since last time

```bash
cd /home/kkourtza/soplasma-scratch
seen=$(cat .claude/hook-events.seen 2>/dev/null || echo 0)
tail -n +$((seen+1)) .claude/hook-events.log | cut -f2,3 | sort | uniq -c
wc -l < .claude/hook-events.log > .claude/hook-events.seen      # mark as reported
```

Format is `ISO-8601 ⟨tab⟩ TIER ⟨tab⟩ RULE ⟨tab⟩ DETAIL`; tiers `BLOCK` and `WARN`; rules
today are `generated-fvSolution`, `generated-initial-field`,
`generated-changeDictionaryDict`, `git-add-all`, `wmake-while-solver-running`.

Report the counts and, for each, whether the block was RIGHT. **A wrong BLOCK gets
demoted to WARN, not argued with**; a WARN that fired on a real defect gets promoted —
the hook set is tuned from evidence, not irritation. Say so too if `.claude/hooks-off`
exists: the kill switch was left on and nothing is guarding.

## 5. Write or refresh a memory — only for a DURABLE fact

A memory is for what outlives this thread (a mechanism, a trap, a refuted class). Live
state goes in PROGRESS.md; a rule that must never be missed goes in `CLAUDE.md`, loaded
verbatim every session — a memory body is only *recalled when judged relevant*, and two
GOVERNING directives once lived only in memory and were simply absent. Write into
**`/home/kkourtza/.claude/projects/-home-kkourtza-Projects-SoEEDF/memory/`** — never the
`-home-kkourtza-soplasma-scratch` one, or the base forks (129 files there on
2026-09-11). Frontmatter matches the existing files: `name: <slug>`, `description:
<YYYY-MM-DD> — <the finding, with its number>`, then `metadata:` / `  type: project`.

Add one line to that directory's `MEMORY.md` (`- [Title](slug.md) — summary`) and mark
what it replaces `Supersedes [[…]]`. **Order by `ls -t`, never by filename** —
`session-state-2026-09-07c` supersedes `-09-07b`, and the newest state of the Newton
thread is not a `session-state-*` file at all.

## 6. Say whether this should be committed — and stop there

One sentence: what the logical unit is, and whether it is committable now (D6: one unit
per commit; if describing it needs an "and", it is two). **If solver source was touched,
say `/regression-gate` must run first (B5)** — a moved baseline is a REGRESSION, an
INTENDED IMPROVEMENT updated in the same commit with its reason, or a STALE BASELINE;
if you cannot tell which, ask. `/save-state` does not stage, commit or push. Hand over.

**Self-check before reporting done** — facts, not reminders:

```bash
cd /home/kkourtza/soplasma-scratch
grep -m1 -i 'Last updated:' PROGRESS.md                            # today's date
awk '/^## 1\./{f=1;next} /^## 2\./{f=0} f' PROGRESS.md |
  grep -v '^[[:space:]]*$' | head -6 | grep -c 'NEXT CONCRETE ACTION'   # must be 1
grep -c '^- \[x\] [0-9]\{4\}-' PROGRESS.md                         # every done item dated
awk '/^## 7\./{f=1} f' PROGRESS.md | grep -c 'Tried'               # resolved notes in §7 = 0
```
