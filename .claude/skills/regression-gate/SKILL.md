---
name: regression-gate
description: Rule B5 made executable — run the three fast test tiers before any commit touching solver source, compare every verification bed to its baseline by RELATIVE tolerance, and classify each changed number as REGRESSION, INTENDED IMPROVEMENT or STALE BASELINE. Use before committing solver/library source, or whenever a baseline has moved and you must say why.
allowed-tools: Bash, Read, Grep
---

# regression-gate

**The question:** did this change move a number that was right before — and if so, which of the
three things is it? This is the artifact rule **B5** names; the baselines exist and until now
nothing compared them (`PROGRESS.md:179`).

**Use it:** before any commit touching `src/`, and whenever a baseline has moved.
**Do NOT use it** to decide physics. Every bed here is tier 1–2 evidence (A7): analytic ground
truth and order of accuracy. It says nothing about whether a discharge is right.

---

## 0. Preconditions (30 s)

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc && source /home/kkourtza/soplasma-scratch/etc/bashrc
cd /home/kkourtza/soplasma-scratch && ./check-no-running-solvers.sh        # exit 0 required
find src/<the component you edited> -name '*.[CH]' -newer $FOAM_USER_LIBBIN/<its .so>
```

PASS: `check-no-running-solvers.sh` exits 0, and the `find` prints **nothing**. Order of the two
`source` lines matters — the project's is second so `PETSC_DIR` wins; `etc/bashrc` does NOT put
binaries on PATH, OpenFOAM's does. A non-empty `find` means a partial `wmake` left an
ABI-mismatched library, and the gate then tests a build that does not exist: run `./build-all.sh`
(background it, 3–20 min) and require `BUILD-COMPLETE` with no `FAIL`/`UNCOVERED`.

## 1. Tier 1 — unit beds with analytic ground truth (seconds)

```bash
testWallFlux && testWallLoss && testVibRelax && testCoulombHeating \
  && testEmission -case /home/kkourtza/soplasma-scratch/verification/fluxScheme1D
```

PASS: exit 0 **and** the literal lines, all measured 2026-09-11 —
`ALL PASS: 35 checks, 0 failed` / `all checks passed` / `all checks passed` /
`5/5 checks passed.` / `PASSED`. Never accept absence-of-errors as a pass (A1).

* `&&`, **never `;`** — with `;` a failing bed is stepped over and the shell reports only the
  last exit status.
* The first four are mesh-free (no `setRootCase.H`/`createMesh.H`); `testEmission` needs a case
  with a mesh, hence the `-case`. It also takes `-patch <name>`.
* **This IS the `--fast` mode (B6)** for any wall-flux, wall-loss, vib-relaxation,
  Coulomb-heating or emission question. Warm cost ~0.03 s each; the FIRST run after a build costs
  24–45 s of cold page-in (measured: testWallFlux 44.62 s cold, 0.03 s warm). Run twice, quote
  the second.
* **`testAitken` is excluded from this list and from every pass/fail list.** `testAitken.C:1803-4`
  ends `Info<< nl << "=== end ===" << nl << endl; return 0;` and the file contains no
  `check()`/`nFail` machinery at all: it ALWAYS exits 0 and prints no verdict. It is an
  exploratory contraction-sweep bed a human reads. Putting it in a `&&` chain creates a check
  that can never fail. Same for `testFluxScheme` — it always returns 0; its verdict is produced
  by the caller (step 3).
* A green bed only covers what it VARIES. `testWallFlux` had 32 passing checks while the defect
  shipped, because every `gRatio` check used `Dd = 0` where right and wrong forms coincide.
  Before believing it, check the parameter you changed is actually swept (A1, "the answer COULD
  have differed").

## 2. Tier 2 — the electrostatics suite (seconds)

```bash
/home/kkourtza/soplasma-scratch/tools/run_electrostatics_tests.sh
```

PASS: exit 0 and a final line `  N ok, 0 failed`. **Exit 2 means `WM_PROJECT_DIR` is unset — the
environment was not sourced, NOT a test failure.** Run it after any change to the Poisson path,
multi-region coupling, `surfCharge`, or a boundary condition. It exists because `fe827ae`
(2026-08-11) killed every electrostatics case with SIGFPE and nobody noticed for three weeks
(found by accident 2026-09-03, fixed in `82fc587`).

* A row reading `RAN   no analytic reference -- completion only` is a completion check.
  **Never upgrade a RAN to a PASS in a report.**
* **It is not read-only.** It `Allclean`s and `Allrun`s five tutorial directories, writes
  `log.suite` into each, and `Allclean`s every case that passes — so a standalone
  `tools/check_series_stack.py` run must come BEFORE it, not after (afterwards there is no time
  directory and it exits 1).

## 3. Tier 3 — verification beds vs their baselines, by RELATIVE tolerance

Run the bed, then compare. **Never `diff`**: the `shear=0.0` rows of the 2D bed carry genuine
~1e-9 run-to-run drift, so `diff` fires on a healthy run.

```bash
cd /home/kkourtza/soplasma-scratch/verification/fluxScheme1D && ./Allrun   # full: 504 RESULT rows
cd ../fluxScheme2Dnonortho && SHEARS="0.0 0.20" ./Allrun                   # full: 48 RESULT rows
cd ../fluxScheme2Dnonortho && SHEARS=0.20 ./Allrun                         # --fast: 24 rows (B6)
```

The 1D `Allrun` has **no env knob** for its `PE`/`N` loops (unlike the 2D bed's `SHEARS`), so its
`--fast` form is one point, the loop body verbatim (`Allrun:30-48`):

```bash
cd /home/kkourtza/soplasma-scratch/verification/fluxScheme1D
sed 's/NCELLS/40/' system/blockMeshDict.in > system/blockMeshDict
sed -e 's/AVAL/1.0/' -e 's/BVAL/2.0/' 0/n.in > 0/n
blockMesh > logs/log.blockMesh 2>&1 || echo "blockMesh FAILED"
sed 's|div(phi,n) .*;|div(phi,n) Gauss ROUNDF;|' system/fvSchemes.in > system/fvSchemes
testFluxScheme -v 100 -D 1.0 -S 1.0 -a 1.0 -b 2.0 -scheme standard | grep '^RESULT'
grep 'scheme=std:ROUNDF ' results.baseline.txt | grep 'N=   40 ' | grep 'Pe=1.000000e+02'
# sub-second. It rewrites three TRACKED files — restore them afterwards:
git checkout -- system/blockMeshDict system/fvSchemes 0/n
```

Then the comparator. It handles both beds (the 1D rows have no `nCorr` field, the 2D rows do),
gates L2 and Linf at 1e-6 relative and `nCorr` at exact equality, and **prints the row counts
first** — a silently short `results.txt` (the per-run `testFluxScheme` calls are unguarded;
only `blockMesh` is) shows up as `in baseline only: N`, not as a pass:

```bash
~/ct-env/bin/python - results.txt results.baseline.txt 1e-6 <<'PY'
import re,sys
tol=float(sys.argv[3])
def load(p):
    d={}
    for ln in open(p):
        if not ln.startswith('RESULT'): continue
        kv=dict(re.findall(r'(\w+)=\s*(\S+)',ln))
        d[(kv.get('shear','-'),kv['scheme'],kv['N'],kv['Pe'])]=(
            float(kv['L2']),float(kv['Linf']),int(kv.get('nCorr',-1)))
    return d
a,b=load(sys.argv[1]),load(sys.argv[2])
print('rows: new',len(a),'baseline',len(b),
      '| in baseline only:',len(set(b)-set(a)),'| new only:',len(set(a)-set(b)))
bad=[]
for k in sorted(set(a)&set(b)):
    r=max(abs(a[k][i]-b[k][i])/max(abs(b[k][i]),1e-300) for i in (0,1))
    if r>tol or a[k][2]!=b[k][2]: bad.append((r,k,a[k][2]-b[k][2]))
print(len(bad),'of',len(set(a)&set(b)),'rows deviate >%g or changed nCorr'%tol)
for r,k,dn in sorted(bad,reverse=True)[:8]:
    print('  %-46s rel=%.3e dnCorr=%+d'%(str(k),r,dn))
sys.exit(1 if bad else 0)
PY
```

PASS: `0 of N rows deviate`, with `in baseline only: 0`. Measured 2026-09-11: the 1D bed gives
`0 of 504` (its two files are bit-identical, md5 `efe4eff01b11cd0ae4af242772d626ae`); the 2D bed
gives `24 of 48`, max `rel=9.477e-01`, `dnCorr=+14`.

Two standing facts about these beds, so you do not re-derive them: the 2D `report.py` matches
**0 of 48** rows (its regex predates the `shear=` prefix) and the 2D `Allrun` never calls it —
never read 2D orders from it. The 1D `report.py` silently drops 72 of 504 rows (all `ROUNDW`,
all `CompleteFlux`) via a hardcoded scheme list at line 29.

## 4. CLASSIFY every changed number — this is the point of the skill

Three buckets (B5), not interchangeable:

| verdict | what you do |
|---|---|
| **REGRESSION** | revert or fix. **Do not touch the baseline.** |
| **INTENDED IMPROVEMENT** | `cp results.txt results.baseline.txt` and commit it IN THE SAME COMMIT as the source change, message stating the reason and the date (D2, D6). Stage explicitly — never `git add -A` (D4). |
| **STALE BASELINE** | say so explicitly, in writing, with the evidence. Refreshing it is still a separate, stated act — not a side effect. |

**If you cannot tell which, you do not get to guess — STOP AND ASK (B5).** Never update a
baseline to make a test pass: that converts a failure into a silent, permanent wrong answer, and
a baseline refreshed without a stated reason has stopped being a test.

Mechanical evidence for the classification (`verification/` has been tracked since `bec57c9`,
2026-09-11, so git works here now — but the history starts there, and **mtimes are the real
provenance**):

```bash
ls -la --time-style=+%F_%T results.txt results.baseline.txt $FOAM_USER_APPBIN/testFluxScheme
git log -1 --format='%h %ad %s' --date=short -- results.baseline.txt
git diff --stat verification/          # both results.txt and results.baseline.txt are tracked
```

An **`nCorr` change, or an `nCorr` field appearing where the baseline has none, is a HARNESS or
BINARY change, not a scheme change** — treat it as the first suspect, not the physics (A3/A5).

### Worked example — `fluxScheme2Dnonortho` is a STALE BASELINE (2026-09-11)

1. Comparator: **24 of 48** rows deviate, every one of them `shear=0.20`; max `rel=9.477e-01`.
   The 24 `shear=0.0` rows agree to ~1e-9 — that is the round-off floor, not a signal.
2. Every deviating row changed `nCorr` 2 → 15/16 (`dnCorr` +13/+14). A harness change.
3. That harness change is the bed's own documented bug #2, *"NON-ORTHOGONAL CORRECTORS MUST
   RE-ASSEMBLE THE EQUATION"* — so `nCorr=2` is the **pre-fix** data.
4. Mtimes: `results.baseline.txt` and `README.md` are both `2026-09-08 00:37:20`; `results.txt`
   is `2026-09-08 00:40:48` — the new data was written **three minutes after** the baseline.
5. The errors are ~2x **lower** and the fitted orders move from 0.00 to 0.97–1.00. A regression
   does not simultaneously improve the answer and fix a documented harness bug.

→ **STALE BASELINE.** Consequence beyond the gate: the bed's `README.md` headline ("SG AND CFS DO
NOT CONVERGE on a non-orthogonal mesh", order 0.00) quotes the stale rows and is refuted by the
`results.txt` beside it. Report the contradiction under A3 rather than picking a side — and note
`ScharfetterGummel.H:116` still uses orthogonal `mesh.deltaCoeffs()` (no `nonOrthDeltaCoeffs`,
no `nonOrthCorrectionVectors`), so the operator-level defect is real even if the convergence
claim is not.

### Second live instance — the 1D baseline predates the current binary

`results.txt`/`results.baseline.txt` are 2026-09-08 00:11:34/00:11:37 (the baseline is a copy of
the results, 3 s later — which is why they are bit-identical and why the 1D gate reads clean
today). `testFluxScheme` was rebuilt at 00:40:46, **29 minutes later**. Measured at one point
today (ROUNDF, N=40, Pe=100): the current binary gives `L2=8.47077435e-02` against the baseline's
`8.48360063e-02` (1.51e-3 relative) and emits an `nCorr=50` field the baseline rows lack. So the
first full 1D re-run will flag all 504 rows, and that is a **STALE BASELINE**, not a regression —
but verify it the same way (mtimes, `nCorr`, direction of the error) before saying so.

## 5. Never

* Never `diff` a results file against its baseline. Relative tolerance, or nothing.
* Never put `testAitken` or `testFluxScheme` in a `&&` gate — both always exit 0.
* Never judge a step by a pipeline's exit status: `./build-all.sh | grep BUILD-COMPLETE` returns
  *grep's* status. Redirect to a file and grep the file.
* Never upgrade a `RAN   no analytic reference` to a `PASS`.
* Never quote a `--fast` result as physics (B6). The coarse bed does not resolve a streamer.
