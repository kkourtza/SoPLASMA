---
name: doc-sweep
description: Rule D1 made executable. After a dictionary-key rename, an option rename or an API change, sweep the WHOLE documentation and case surface in the SAME commit — both trees, docs/, tutorials, smoke beds, validation COMPARE.md, etc/boundaryRoles — then prove the old name survives nowhere except a dated SUPERSEDED marker (D2).
allowed-tools: Bash, Read, Edit, Grep
---

# doc-sweep — the D1 rename/API sweep

**Question it answers:** after I renamed `<old>` to `<new>`, what is still telling a later
session to use `<old>`? Anything left behind is not stale — D1: *it is an INSTRUCTION a
later session will follow.*

**Use when:** a dictionary key, option, utility name, field name, boundary `kind` or public
API signature changed; or a key was deleted or split.
**Do NOT use for:** a source-internal variable rename no dictionary or doc mentions; a new
measurement (D2/D5 → `PROGRESS.md`); a brand-new key (skip to steps 4-5).
The sweep ships in the SAME commit as the rename (D6); deferring it is the failure this
skill exists to stop.

## 1. Fix the two strings, then PRINT the surface

Write `<old>` and `<new>` down verbatim first, then print the surface instead of recalling
it — counts measured 2026-09-11 in `/home/kkourtza/soplasma-scratch`:
```bash
find docs -name '*.md' | wc -l                                   # 45
git ls-files -- '*COMPARE.md' | wc -l                            # 75  (69 under validation/)
git ls-files -- smoke | wc -l                                    # 167 across 5 beds
git ls-files -- '*/configuration/config' | grep -vc validation   # 16
find /home/kkourtza/Projects/SoEEDF/docs -name '*.md' | wc -l    # 37
```
Surface set (D1): both `README.md` (987 l / 47204 B SoPlasma, 1158 l / 79657 B SoEEDF);
`docs/**/*.md` in BOTH trees — `docs/CAPABILITIES.md`, `docs/reference/` (6), `docs/models/**`
(10), `docs/design/*.md` (22), `docs/theory/`, `docs/getting_started/`,
`docs/simulationManuals/`, SoEEDF `docs/INDEX.md` + 36; `SoEEDF/mechanisms/FORMAT.md` and
`SoEEDF/{data,meshes,Literature}/README.md`; `CLAUDE.md`, `PROGRESS.md`, `PR-DESCRIPTION.md`,
`validation/TUNE_*.md`, `validation/RESOLUTION_TEST.md`; 8 `tutorials/**/tutorial_info.md`
+ 5 smoke ones; **every tutorial, smoke and validation dictionary INCLUDING its COMMENT
text**; every `configuration/{config,boundaries}`; `etc/boundaryRoles` (445 lines).

## 2. The five passes — run ALL five, grepping BOTH names
```bash
cd /home/kkourtza/soplasma-scratch
# A  tracked live tree                                            0.065 s warm
git grep -nE '<old>|<new>' -- ':!validation'
# B  the places D1 names explicitly — comments included            0.54 s
git grep -nE '<old>|<new>' -- '*.md' etc tutorials smoke
# C  validation beds, two cheap slices                             0.63 s / 15 s
git grep -nE '<old>|<new>' -- 'validation/**/COMPARE.md'                        # 69 files
git grep -nE '<old>|<new>' -- 'validation/*/configuration/*' 'validation/*/system/*' \
                              'validation/*/constant/*' 'validation/*.md'       # 3925 files
# D  second tree — 759 tracked files, no pathspec needed           ~5 s
cd /home/kkourtza/Projects/SoEEDF && git grep -nE '<old>|<new>'
# E  untracked / generated / out-of-tree cases git cannot see
grep -rn '<old>' ~/Projects/SoEEDF ~/soplasma-scratch --include='*.md' --include='*Properties'
```
* Grep the NEW name too: its hits are how you see which files were already updated.
* **Never** `grep -rn ... .` from the SoPlasma tree root — measured killed at 120 s, exit
  143, no output. The `':!validation'` pathspec is what makes git grep finish; unscoped it
  takes 1m37s of pure I/O stall (WSL2, deleted `validation/` index entries). The *unfiltered*
  recursion is the anti-pattern — pass E's filtered form is cheap and is the only one that
  reaches untracked trees.
* Pass B PASS: each of `smoke/{native,native-dt2,native-stiff,native-stiff-dt2,cantera}`
  appears or is explicitly stated clean. These five are the documented worst case — ignored
  wholesale, so no sweep reached them, and by 2026-09-02 all five carried FOUR rejected
  settings and **none could start**: `limitVoltageRiseRate`+`maxVoltageRiseRate`, a dangling
  `appliedVoltage`, a missing `electronEnergyModel`, an unsplit `chemistrySolver`
  (`.gitignore:24-45` carries that post-mortem inline).
* A COMPARE.md carries the extraction command and the field names; a renamed field silently
  invalidates it. G2: generated layer-2 dictionaries are gitignored, so no git pass sees
  them — **regenerate them, never hand-edit a hit inside one.**

## 3. Classify every hit: historical record, or live instruction?

The sweep's one judgement. D2: superseded text is MARKED IN PLACE with its date, **never
deleted** — a finding that quietly disappears leaves no way to recognise its stale copies
elsewhere. The in-tree pattern to copy is the smoke beds':
```
// MIGRATED 2026-09-02 from `limitVoltageRiseRate false` + ...   (system/plasmaSimulationControls:37)
// SPLIT 2026-09-02 from the single `chemistrySolver adaptive;`  (constant/plasmaTransportProperties:39)
```

## 4. A DICTIONARY key changed → `docs/reference/` changes in the same commit

Ground-truth default and file:line from the source, never from memory:
```bash
python3 tools/checkConfigReference.py --list-source | grep -w '<new>'
```
Write the entry in the only form the checker parses — ``### `<new>` `` followed, before the
next heading, by ``**Default:** `value` `` or `**REQUIRED -- no default.**`. Without that
shape the claim is invisible to the gate: only 7 level-3 headings in `docs/reference/` carry
a backticked key today, and the two largest files have zero.

## 5. Run the honesty checker (baselines, so regression ≠ baseline)
```bash
python3 tools/checkConfigReference.py                                        # ~8 s
python3 tools/checkConfigReference.py --doc docs/reference/<one>.md          # --fast (B6)
python3 tools/checkConfigReference.py --case tutorials/plasma/soPlasmaFoam/needleDBD  # ~2 s
```
Reproduced 2026-09-11 — the gate is ALREADY RED, so PASS means *nothing worse than*:
* default mode **exit 1**: MISSING READER `uniformValue` (`docs/reference/changeDictionary.md:33`);
  STALE DEFAULT `defaultSEEC` (doc `0.001`, source `<derived; see the reader>`); footer
  `3 documented, 1 verified against source; 301 options read by src/`.
* `--case` **exit 1**: DANGLING `simulationType`; 5 UNRESOLVED (`ePotentialLaplacianScheme`,
  `electronDiffusionLaplacianScheme`, `electronDiffusivityInterpolationScheme`,
  `electronDriftDivScheme`, `electronMobilityInterpolationScheme` — nested in compound
  values, NOT failures); footer `14 variables declared; 8 resolve to a key the source reads;
  5 unresolved; 1 dangling`. Only DANGLING makes it non-zero. Run it once per case whose
  `configuration/config` you touched.
* FAIL = `<old>` appears anywhere in the output, `<new>` appears as MISSING READER or STALE
  DEFAULT, or a new DANGLING name appears.
* Trap: a MISSING READER can be a false positive from the script's own allowlists
  (`uniformValue` is OpenFOAM-owned and simply absent from `OPENFOAM_OWNED`). Before
  deleting a doc line, `git grep -n '<key>' -- src` yourself.

## 6. Links and inventory
```bash
cd /home/kkourtza/soplasma-scratch/docs/reference && grep -o '](\([^)]*\))' README.md \
  | sed 's/](//;s/)$//' | while read l; do case "$l" in http*) continue;; esac; \
    [ -e "$l" ] || echo "BROKEN: $l"; done      # baseline: BROKEN: regions-and-materials.md
```
Two further dead rows there are plain text (`*pending*`) and this check cannot see them.
If the change alters what the framework can DO, `docs/CAPABILITIES.md` is the live inventory
and moves in the same commit; SoEEDF's `framework-state.md` and `capability-overview.md` are
already stale by an entire solver — do not copy into them.

## 7. VERIFICATION — the gate

Re-run every pass on the OLD name ALONE. Required: zero hits, or every surviving hit on a
line that also carries a dated supersession marker.
```bash
cd /home/kkourtza/soplasma-scratch
git grep -nE '<old>' -- ':!validation'
git grep -nE '<old>' -- '*.md' etc tutorials smoke
git grep -nE '<old>' -- 'validation/**/COMPARE.md'
git grep -nE '<old>' -- 'validation/*/configuration/*' 'validation/*/system/*' 'validation/*/constant/*'
cd /home/kkourtza/Projects/SoEEDF && git grep -nE '<old>'
```
`git grep` exits 1 on zero hits — **exit 1 here IS the pass.**

Then prove D2 was obeyed rather than dodged — the diff must ADD markers, not remove records:
```bash
cd /home/kkourtza/soplasma-scratch && git diff -U0 -- smoke tutorials validation docs \
  | grep '^-' | grep -iE '<old>'
```
Any removed line that was a historical record (a dated comment, a post-mortem sentence, a
COMPARE.md reference number) is a D2 violation — restore it and mark it in place instead.

D4: state plainly what is uncommitted (`git status --short`), stage explicitly — never
`git add -A` on a tree holding case output — and commit the rename with its sweep as one
unit (D6).
