---
name: new-case
description: Create a SoPlasma case, or derive an arm from an existing one, through the layer-1 generators (plasmaSetupBoundaries + plasmaCreateSpeciesFields) instead of hand-writing OpenFOAM dictionaries. Use when setting up any new case, converting a case to layer 1, or copying a case for a new mesh/voltage/geometry arm. Carries the dangling-$key proof, the multi-region checks and the half-converted trap.
allowed-tools: Bash, Read, Edit, Write, Grep, Glob
---

# new-case — build a case through the generator, not by hand

**Question this answers:** what do I actually write, and what does the framework write for me?
**Answer (G2):** you write TWO files — `configuration/config` (values) and
`configuration/boundaries` (what each surface IS). Every `0/` field is generated from those
plus `etc/boundaryRoles` and the mesh topology. `needleDBD` is the reference: 16 tracked
files, of which the layer-1 pair is two, and `0/` is 17 generated files with no `0.orig`.

**Do NOT use this skill for:** re-running an existing case after a code change
(`CASE=<dir> TIMEOUT=300 ./rerun.sh` — Allrun-serial rebuilds the mesh and costs minutes);
a seconds-scale debug bed off the streamer family (`NCELL=130 NREFINE=1 CASE=$HOME/streamer-fast
./make-smoke-case.sh`, B6); or a unit bed under `verification/`, whose `0/` and `system/` are
deliberately hand-authored and outside layer 1.

---

## 0. Environment — both, project second

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc && source $HOME/soplasma-scratch/etc/bashrc
[ -f "$SoPLASMA_ETC/boundaryRoles" ] || echo "FAILED: SoPLASMA_ETC wrong: $SoPLASMA_ETC"
```
PASS: no `FAILED` line, and `which plasmaSetupBoundaries` resolves. Both generators install
under names that differ from their source directories: `foamPlasmaSetupBoundaries` →
`plasmaSetupBoundaries`, `foamPlasmaCreateSpeciesFields` → `plasmaCreateSpeciesFields`.

## 1. If you are deriving from an existing case: check for the half-converted trap FIRST

```bash
grep -c plasmaSetupBoundaries <case>/Allrun*; ls -d <case>/0.orig <case>/etc 2>/dev/null
```
PASS (fully generated, e.g. `needleDBD`): count >= 1 **and** neither `0.orig` nor `etc` exists.

Anything else is **half-converted**: `etc/changeDictionary.<region>` runs AFTER the generators
and **WINS**, and `0.orig/` is a second hand-authored source of boundary truth. Measured
2026-09-11: `positiveStreamer_fixedMesh` returns `Allrun-serial:0` with both directories
present, and ~30 `validation/` arms are identical clones. On those beds the two paths disagree
on PHYSICS — charged species at the electrodes are `electron/ionDDWallFluxMixed` (ABSORBING)
under layer 1 versus `zeroGradient` (NON-absorbing) under changeDictionary. **That fork is
frozen deliberately (PROGRESS.md §4); do not "fix" it as a side effect of making a new case.**
If you copy such a case, you inherit both sources — say which one you intend to own.

## 2. Write layer 1 — and only layer 1

* `configuration/boundaries` — one block per patch, named exactly as the mesh names it.
  **Never declare** `front`/`back` (mechanical, read from the mesh), `<gas>_to_<dielectric>`
  and its partner (created by `splitMeshRegions`, DERIVED), or `procBoundary*`.
* `configuration/config` — a value only if some dictionary references it as `$key`. A value
  used in exactly one dictionary belongs in that dictionary (measured on needleDBD: 22 of 37
  variables were used once; config went 37 → 15).

The nine shipped kinds (`etc/boundaryRoles`), with what each REQUIRES and what species do there:

| kind | requires | species surface |
|---|---|---|
| `drivenElectrode` | `waveform` (any Function1) | conductor |
| `groundedElectrode` | — | conductor |
| `floatingElectrode` | — | conductor |
| `thinDielectricSurface` | `material` | chargingSurface |
| `thinDielectricOnElectrode` | `material thickness backingPotential` | chargingSurface |
| `openBoundary` | — | open |
| `ballastedElectrode` | `circuit` | conductor |
| `currentDrivenElectrode` | `circuit` | conductor |
| `insulatingWall` | — | conductor |

Full descriptions, optional parameters (`material`, `gammaSEE`, `electronReflection`) and the
emission mechanisms: `plasmaSetupBoundaries -listKinds`, run **from inside a case directory**
(the listing sits after `setRootCase.H`, so it needs a valid case root, not just `SoPLASMA_ETC`).
PASS: `Boundary kinds in <path>/etc/boundaryRoles`, nine kinds, ending
`Region INTERFACES and MECHANICAL patches are DERIVED and must not be declared.`

## 3. NEVER hand-edit `system/*` or `0/*` (G1, B4)

A hook (`.claude/hooks/guard-generated-files.sh`) blocks what it can identify, and it only
knows what `.gitignore` records — so it is a backstop, not the rule. It blocks:
`0/<field>` where the case has `configuration/` or `0.orig/`; `system/fvSolution` only where
`fvSolution-foam`/`-petsc` exists beside it; `system/changeDictionaryDict`; and
`constant/<region>/{plasmaSpeciesProperties,plasmaTransportProperties,photoionizationProperties}`.
It WARNS (never blocks) on a bare numeric literal entering `system/plasmaSimulationControls`
or `system/controlDict` — the 2026-09-10 defect was a literal `maxDeltaT` typed there while
`plasmaTimeControl` read its own dict and all three arms ran `dt=1.29e-9` byte-identical.

## 4. Mesh, and the multi-region path

Single region needs neither step: `plasmaSetupBoundaries` handles a case with no
`constant/regionProperties` ("No constant/regionProperties: treating this as a...").

```bash
# system/fvSchemes AND system/fvSolution must exist BEFORE the split -- fvMesh reads them
# MUST_READ. splitMeshRegions then OVERWRITES system/<region>/ with an EMPTY dummy.
splitMeshRegions -cellZones -overwrite
grep -A6 "Sizes of interfaces" log.splitMeshRegions | grep -E "^[0-9]"
```
PASS: at least one numeric row (needleDBD's own guard is stricter — `grep -qE "^0[[:space:]]"`).
**`-A6`, not `-A3`.** The data row is the FOURTH line after the header (blank, column names,
dashes, row), so `-A3` prints an empty-looking table and makes a correct split look exactly
like the failure this check exists to detect. An empty table means the regions never coupled
and **nothing else errors**.

```bash
tools/plasmaSetupRegions.sh > logs/log.plasmaSetupRegions 2>&1   # IMMEDIATELY after the split
grep -E -- "->|NEEDS A VALUE|exists" logs/log.plasmaSetupRegions
grep -q 'ACTION REQUIRED' logs/log.plasmaSetupRegions && echo "FAILED: uncomment epsilonR"
```
PASS: one `system/<r>/{fvSchemes,fvSolution}  ->  #include ../{fvSchemes,fvSolution}` per region.
The `ACTION REQUIRED` banner is a bare `echo` with **no exit** (`plasmaSetupRegions.sh:440-455`),
so `|| die` cannot catch it — grep the log, as above.

## 5. Generate `0/` — both generators, order irrelevant

```bash
rm -rf 0 && mkdir -p 0/gas 0/dielectric logs
plasmaSetupBoundaries && plasmaCreateSpeciesFields
```
PASS: `plasmaSetupBoundaries: N field(s) generated.` with N = regions + 1 (`surfCharge` is
written for the gas region), and `Boundary conditions DERIVED per patch (wall-flux family …):`
followed by one `Creating field: n_<sp>` per species. **Do not grep for `Mixed`** — the family
is read from `wallFluxFamily` and may be `Implicit`; grep the stable prefix
`"Boundary conditions DERIVED per patch"`. Order does not matter: verified 2026-09-04, all 17
generated files on needleDBD are identical whichever runs first (13 `n_<sp>` + `nEps_e` +
`gas/{ePotential,surfCharge}` + `dielectric/ePotential`).

## 6. Prove no `$key` is dangling — this is a gate, not a reminder

```bash
for k in $(grep -oE '^[a-zA-Z][a-zA-Z0-9_]*' configuration/config | sort -u); do
  grep -rlF -- '$'"$k" system constant configuration etc 0.orig Allrun-serial Allrun-parallel 2>/dev/null \
    | grep -q . || echo "DANGLING: $k"
done
```
PASS: no output. **Two traps, both measured.** `grep -rq "\$$k"` inside DOUBLE quotes expands
`$$` to the shell PID and falsely flags every key (49 of 49 on `positiveStreamer_LMEA_fast`);
use `'$'"$k"`. And `grep -rq` returns exit 2 when any named path is missing **even after a
match**, so `||` fires for everything — hence `grep -rlF … | grep -q .` above.

On `needleDBD` today this prints exactly one line, `DANGLING: simulationType` (reproduced
2026-09-11, 14 keys) — **a real, live defect already logged in PROGRESS.md §6**, not a false
positive. An inert `$key` is silent: a dangling `appliedVoltage` once invalidated a whole
"low-field" arm that ran at the original 62.1 Td and exited 0 (A1: the knob must have MOVED).

Deeper check: `tools/checkConfigReference.py --case <dir>`. Its current baseline on needleDBD
is **exit 1** with `simulationType` DANGLING plus 5 UNRESOLVED names (schemes nested in
compound values — not failures). Anything worse than that is a regression.

## 7. Run it

```bash
./Allrun-serial            # from the case directory
grep -q '=== done' log.* logs/log.* 2>/dev/null
```
PASS: **`grep '=== done'`, not "the last line"** — `Allrun-serial:213` prints
`=== done. Worth reading:` and three further echo lines follow. And an `EXIT=0` from any Allrun
means only that the SCRIPT finished; `runApplication` does not propagate the solver's status:

```bash
grep -c '^Time = ' log.soPlasmaFoam; grep -E 'FOAM FATAL|sigFpe|sigSegv' log.soPlasmaFoam
```

**`--fast` (B6):** comment out step 8 (`runApplication "$(getApplication)"`,
`Allrun-serial:210`). Steps 1-7 build mesh, regions, ion tables and all 17 generated fields in
**well under a minute** and exercise the entire layer-1 chain. Use this for every setup change;
never quote it for physics. needleDBD's full run is ~2.8 s/step (tens of minutes at
`endTime 3e-09`, ~3 h at `2e-08`).

## 8. Deriving a mesh arm — and what `make_mesh_arm.sh` is NOT

`tools/make_mesh_arm.sh <parent_case> <target_case> <NX> <BUMP>` (~1-2 min).
PASS: final line `ARM-OK <target>  NX=<n> Bump=<b>  cells=<n>  fields=<n>`.

**It is not a general copy-a-case recipe.** It is hardcoded to ONE geometry family and dies on
anything else: `gap1cm.geo` / `gap1cm.msh` / `gap1cm_3D.msh`, a `NX = <n>;` and `Using Bump <b>;`
sed against that `.geo`, thickness `-t 2e-4`, an inline rename of `side_lo`/`side_hi` to
`symmetryPlane`, and `rm -rf 0 && mkdir -p 0` — **single region, no `splitMeshRegions`, no
`plasmaSetupRegions.sh`**. `find . -name gap1cm.geo` matches only `validation/` arms. It also
does not export `BOLTZMANN_DIR` or `PETSC_DIR`. For any other geometry, follow steps 4-7.

Its one transferable decision: it copies the parent's `constant/plasmaTables` and
`constant/ionTables` rather than regenerating them — table generation is mesh-INDEPENDENT and
costs minutes. Keep them in SEPARATE directories; a shared one is destroyed by `Allclean` and
the Boltzmann sweep does not put `ionTables` back.

## 9. Before you call the case done

* **B1** — set `writeInterval <= endTime/10` NOW. The project's first genuine ignition wrote
  ZERO time directories because `writeInterval 5e-6` was set against `endTime 20e-6`.
* **D3** — one to three probes, every species density, `chargeDensity`, `Emag`, reduced field,
  electron energy, currents, plus `dt` AND which limiter set it. If diagnosing the case needs
  the log parsed, the monitoring is the defect.
* **D1** — if this case is a comparison, write `COMPARE.md` in its own directory **now**:
  question, the baseline's ABSOLUTE path, the reference numbers, the extraction command.
  Check it with `python3 /home/kkourtza/Projects/SoEEDF/tools/compare_cases.py <caseDir>`
  (that tool lives in the OTHER tree).
* **B3** — state the full parameter set before launching: circuit type and values, applied
  voltage/current and its profile, `simulationType`, gas/pressure/temperature, `endTime`,
  timestep strategy. Background anything over 1-2 minutes.
* **D4** — `git status`, then say plainly what is uncommitted. Stage explicitly by path:
  `configuration/`, `constant/` dictionaries and `system/` dictionaries are the work; `0/`,
  time directories and `logs/` are run output. Never `git add -A`, never `git add <dir>` on a
  tree holding case output (once staged 436 paths of output for 10 real dictionary changes).
