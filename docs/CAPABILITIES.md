# CAPABILITIES, TOOLING, AND WHAT HAS BEEN REFUTED

**THIS IS THE FILE TO RE-READ WHENEVER MY CONTEXT WAS REBUILT** -- after a
compaction/summary resume, a `--continue`/`--resume`, or a fresh session. Rule
B5 makes that mandatory. It exists because a compaction summary preserves the
NARRATIVE of what we did and loses the INVENTORY of what we have, and twice on
2026-09-10 that cost real work:

* I proposed current-source electrode control as "a genuine piece of work". It
  was already implemented and documented in three places. (Rule B4.)
* I hand-rolled a case cold-start with `nohup`, two `exit 127`s and the wrong
  utility name, when `tools/make_mesh_arm.sh` already encoded every step
  including the traps. The user: *"its like you dont remember anything of what
  weve done anymore"*.

This file POINTS. It does not duplicate. If a fact lives in `README.md` or
`docs/reference/`, the pointer is the entry.

---

## 1. Where a capability question is answered

| question | authority |
|---|---|
| what can the solver do at all | `README.md` -- the de facto capability list |
| every boundary role and circuit type | `etc/boundaryRoles` |
| every dictionary key, with meaning | `docs/reference/` (`configuration-config.md`, `plasmaSimulationControls.md`, `controlDict.md`, `fvSchemes-fvSolution.md`) |
| a model's physics and options | `docs/models/**` |
| design notes, status, and what was refuted | `docs/design/*.md` -- **check the `Status:` line; several are DESIGN, nothing implemented** |
| the rules and their measured evidence | `../Projects/SoEEDF/CLAUDE.md`, `../Projects/SoEEDF/docs/rules-postmortems.md` |

**ONE DOCUMENTATION FOLDER: `docs/`.** Decided by the user 2026-09-10 after I
had created a second top-level `doc/` alongside the pre-existing `docs/` -- two
folders differing by one letter, which is exactly the ambiguity G1 exists to
prevent. `doc/` is GONE; its 22 design notes are now `docs/design/`, and this
file sits at `docs/CAPABILITIES.md`. 94 references across 63 files were
rewritten in the same pass. SoEEDF already had only `docs/`, so both trees now
agree.

    docs/CAPABILITIES.md   this file -- the inventory
    docs/design/           internal design and status notes (was doc/)
    docs/reference/        every dictionary key
    docs/models/           per-model physics and options
    docs/theory/           derivations
    docs/getting_started/  install and dependencies
    docs/simulationManuals/

**Do not create another top-level documentation folder.** A new category is a
subfolder of `docs/`.

**DESIGN ONLY, nothing implemented** (do not cite these as capabilities):
`numerics-generator-plan.md`, `adaptive-dt-design.md`, `steady-and-stability-design.md`,
`case-monitoring-plan.md`, `schur-semiimplicit-poisson-preconditioner.md`.
`external-circuit-plan.md` is HALF done.

---

## 2. TOOLING THAT ALREADY EXISTS -- check here before writing a script

This is the section the compaction summary never carries, and the one that has
cost the most. **Read the header comment of the script before using it**: each
one records the trap it was written to avoid.

### THE GENERATOR -- this is the point of the whole architecture

**The user edits a SMALL SEMANTIC LAYER and everything else is generated.** This
is G1/G2 made real, and it is what makes the solver usable at all. `needleDBD`
is **two hand-written files**: `configuration/config` and
`configuration/boundaries` (4 declarations). No `etc/`, no `0.orig/`. Every
field in `0/` is generated. (Layer 1 landed 2026-09-04, commits `7ef9a2e`,
`bebd70f`.)

| piece | role |
|---|---|
| `etc/boundaryRoles` | **the LIBRARY.** Boundary kinds -- `drivenElectrode`, `groundedElectrode`, `floatingElectrode`, `currentDrivenElectrode`, `thinDielectricSurface`, `thinDielectricOnElectrode`, `openBoundary`, `insulatingWall`. Each gives description, required/optional params, the electrostatics condition verbatim, and ONE word for transport (`surface` = mechanical \| open \| conductor \| chargingSurface). |
| `configuration/boundaries` | **LAYER 1.** Per-patch: the user says WHAT a surface IS, plus its `circuit {}` and `material`. |
| `configuration/config` | **the variable layer.** Values live here; the real dictionaries reference them as `$key`. |
| **`plasmaSetupBoundaries`** | **THE GENERATOR** (`-listKinds` lists the kinds). Generates every potential/surfCharge BC. |
| `plasmaCreateSpeciesFields` | generates every species/energy field, with BCs DERIVED from the potential's own |

**NEVER hand-edit `system/*` or `0/*`.** A value goes in `configuration/config`
and the dictionary references `$key`. A boundary goes in
`configuration/boundaries`. Editing the generated side is the defect G1 exists
to prevent -- and on 2026-09-10 I put `maxDeltaT` as a literal in
`system/plasmaSimulationControls` and had to move it.

**After editing `configuration/config`, PROVE the key is not dangling:**
`grep -rlF '$'"$key" system/ constant/`. A variable nothing references is
silently inert -- `appliedVoltage` was dangling in four beds while the voltage
came from a hardcoded table, and a whole "low-field" arm ran at the original
62.1 Td reporting exit 0. **Use `grep -rF '$'"$k"`, not `"\$$k"`** -- inside
double quotes `$$` expands to the shell PID and reports nearly everything as
dangling. ([[config-variable-may-be-dangling]])

DERIVED and never declared: mechanical patches (from the mesh), region
interfaces (from topology -- sigma is owned by the GAS side only, and ownership
INVERTS between sides with nothing stating it), and `epsilonR`/`gamma` behind a
named material. Generator order does NOT matter (fixed in `db77095`); both
generators derive from `Foam::boundaryRoleLibrary`.

### Mesh and derived-case construction

| tool | what it does |
|---|---|
| `tools/make_mesh_arm.sh <parent> <target> <NX> <BUMP>` | Builds one Grubert mesh arm from a parent: gmsh -> `msh2Dto3D` -> `gmshToFoam` -> boundary fix -> `symmetryPlane` on side_lo/side_hi -> `checkMesh` must say "Mesh OK" -> `rm -rf 0 && mkdir -p 0` -> the two generators. Reuses the mesh-INDEPENDENT Boltzmann/ion tables deliberately. |
| `<case>/Allrun-serial` | the per-case setup chain: mesh, `fvSolution` chosen by `matrixSolver`, `0.orig -> 0`, generators, `changeDictionary`, seed, run |
| `tools/plasmaSetupRegions.sh` | per-region files for a multi-region mesh |
| `make-smoke-case.sh` | small copy of the streamer case for "does the code still run" (NOT for answers) |
| `tools/msh2Dto3D.py` | extrude a 2-D gmsh mesh; also `--fix-boundary` |

Three traps `make_mesh_arm.sh` encodes, all measured:
1. **NOT `set -e`** -- OpenFOAM's `etc/bashrc` returns non-zero, and with `-e`
   the script dies silently with no log and no message.
2. **`export` on SEPARATE lines** -- in `export A=x B=$A` bash expands `$A`
   before assigning it, which set `SoPLASMA_ETC=/etc` and killed
   `plasmaSetupBoundaries` on `Cannot read ... "/etc/boundaryRoles"`.
3. **`rm -rf 0 && mkdir -p 0`** before the field utilities.

The utilities are `plasmaSetupBoundaries` and **`plasmaCreateSpeciesFields`**
(source dir `foamPlasmaCreateSpeciesFields`; there is no `plasmaCreateFields`).

### Running

| tool | what it does |
|---|---|
| `run-guarded.sh [maxLogMB] [maxSeconds]` | hard log cap + wall clock. A runaway warning loop once wrote **9.7 GB / 91M lines** before one timestep finished and took the machine down. |
| `run-long.sh` | long unattended run with that same cap |
| `watchdog.sh` | kills the solver if its log exceeds a cap |
| `rerun.sh` | re-run only the solver on an already-meshed case |
| `check-no-running-solvers.sh` | refuse to build while a solver is running |
| `build-all.sh` | serial build in dependency order, per-component logs |

### Monitoring and analysis

| tool | what it does |
|---|---|
| `validation/status.py` | **robust** status across cases; only reads complete lines (a live log can catch a half-written one) |
| `tools/news.py` | evolution of peak `n_e` and peak `T_e` |
| `check-run.sh` | cheap health check: mechanism, missing species, finiteness, rate |
| `tools/settle_check.py` | avalanche rate, and whether a run has SETTLED |
| `tools/cathode_fall.py` | normal/abnormal glow similarity parameters from a snapshot |
| `tools/glow_stability.py` | small-signal stability of the glow against a lumped ballast |
| `tools/compare_at_voltage.py` | compare two cases at MATCHED source voltage |
| `tools/flux_profile.py` | where the conduction current reverses first |
| `tools/plot_ne_Te.py` | spatial `n_e`/`T_e` from a snapshot |
| `validation/extract_co.py` | dt and limiter at a COMMON simulated time (rule A2) |
| `validation/mesh_sweep.sh` | the mesh-independence sweep |
| `tools/checkConfigReference.py` | keeps the option reference honest against the code |
| `tools/check_series_stack.py` | two-region series dielectric vs analytic |

Python is **`~/ct-env/bin/python`**, not the system python.

### Test utilities built into the tree

`testPoissonSymmetry`, `testSnesJFNK`, `testSnesJFNK2Field`, `testFluxScheme`,
`testWallFlux`, `testWallLoss`, `testDischargeCurrent`, `testEmission`,
`testAitken`, `testCoulombHeating`, `testVibRelax`, `plasmaChemistry0D`.

---

## 3. What WORKS and is verified

| capability | verification |
|---|---|
| **Newton/JFNK outer solver** (`outerSolver newton`, PETSc SNES, matrix-free) | replaces the segregated Picard sweep; runs on the real 5-field system |
| **Newton in PARALLEL** | verified to 1.96e-09, below SNES rtol 1e-8, by fixed-dt cell-by-cell comparison of final fields. The bug was `ISCreateStride(...,0,...)` using rank 0's rows on every rank; fixed with `MatGetOwnershipRange`. |
| **Newton transport models** | `driftDiffusion` and, since 2026-09-10, `immobile` -- the latter returns ZERO mu()/D() so the existing assembly reduces to ddt(n) == sources with no branching. `diffusion` is still REFUSED (it exposes neither), deliberately: assembling it without its transport would silently drop physics. The guard asks `providesTransportCoefficients()`, not a type. |
| **Newton chemistry sources** | Needs per-species sources. `ode`/`adaptive`/`adaptiveError`/`implicitRate` populate chemP/chemL; `explicitSource` records a net source instead (added 2026-09-10, into its OWN field -- writing chemP_/chemL_ would have switched on the Co_chem limiter via maxChemStateRate(), changing dt for unrelated cases). Ask `chemistrySourcesAvailable()` / `chemNetSourceAvailable()` BEFORE reading: the accessors index the lists directly and segfault when unsized. |

| **`PCFIELDSPLIT` + Schur** | splits named `phi`/`transport` via `PCFieldSplitSetIS` (**names, not `0`/`1`**) |
| **Mesh independence of the Newton solve** | KSP median 8/11/9 over **100x** cells at MATCHED Courant; SNES its/step identical; exponent ~0.05 |
| **Current-driven electrodes** | `currentDrivenElectrode` + `circuit { type currentSource; setCurrent; compliance; capacitance; }` -- `etc/boundaryRoles:377`, `README.md`, `docs/reference/plasmaSimulationControls.md:90` |
| **Adaptive dt** | temporal-error PI controller (`plasmaTimeControl`), `maxDeltaT`/`minDeltaT`/`maxInitialDeltaT` supported. `deltaT` in `configuration/config` is only the INITIAL step. |
| **Semi-implicit Poisson** | `poissonScheme semiImplicit` (auto-switched to `explicit` under `outerSolver newton`) |
| **Flux schemes** | `standard` (div+laplacian, so `fvSchemes` decides) and `ScharfetterGummel` (bypasses `divSchemes`/`laplacianSchemes` **entirely and silently**) |

## 3b. Newton gotchas that cost a session (2026-09-10)

* **Newton needs a successful PICARD step before it hands over.** The handover
  gates on every species being off its density floor, and it LATCHES. So a
  COLD-STARTED arm cannot test Newton at any dt where Picard's own first step
  dies -- measured: at dt=1e-10 on positiveStreamer both arms died with the
  same SIGFPE after 1 step, and the Newton arm never handed over at all.
  **Warm-start from an established state instead**; then handover happens on
  step 1 and the two solvers start from a bit-identical field.
* **Verify the handover MESSAGE, never assume it.** "outerSolver newton:
  Picard warm-up COMPLETE ... this handover is PERMANENT" is the only proof an
  arm is actually Newton. Its absence once meant 21,174 steps silently ran
  Picard (commit b29c8a0).
* **`libplasmaNewtonSolverPETSc` needs BOTH bashrcs.** Without the project's
  own `etc/bashrc`, PETSC_DIR is unset, `libpetsc.so.3.24` is not found, the
  dlopen fails silently and SNES reports "Registered types: 0()".
  `libpetscFoam.so` fails the same way -- and did, unnoticed, in a reference run.
* **Adding a VIRTUAL to a shared header is an ABI change.** It shifts every
  later vtable slot; rebuilding only the changed library leaves the others
  calling the wrong slots, which presents as a startup SEGV that looks like a
  physics bug. Rebuild everything (`./build-all.sh`).

## 4. What has been REFUTED or DECIDED AGAINST -- do not re-propose

| thing | verdict, and the measurement |
|---|---|
| **ILU(1) on the transport split** | 2,000-cell artefact. 2k: -33% steps. 20k: **+3.5%**. Not a default. |
| **`selfp` Schur approximation** | same. 2k: -22% steps. 20k: KSP median 11 -> 12. Not a default. |
| **Per-cell residual/state scaling** | same. 2k: -46% steps, SNES 1887->1224. 20k: **+1.4% mean dt, KSP median 13 vs 13, SNES/step 5.8 vs 5.9**. Not a default. (2026-09-10) |
| ILU(1) + `selfp` together | they INTERFERE -- the combination is exactly baseline |
| **Voltage + ballast control for the DC glow** | CANNOT reach the operating point. The CVC has a minimum, so `dR/dI = 0` and no load line selects a point (Almeida & Benilov 2017). Use `currentSource`. Cost of not knowing this: a full day. |
| **Tightening the Poisson tolerance to cure the lateral asymmetry** | REFUTED. It improves `phi`/`Ex` but leaves `Ey` -- the driving component -- unchanged, and GAMG at 1e-16 costs 2000 iterations (286x) for nothing. `1e-16` absolute on a 225 V field is below double precision anyway. |
| **A mesh-symmetrisation utility, or a symmetry CHECK** | rejected with the user. Refinement already cures it (`Ey ~ NY^-1.7`); point-snapping still left 82/1592 volume pairs differing; an unstructured triangular mesh has no mirror symmetry at all, so a check would fire constantly and mean nothing. |
| **`fvMatrix::residual()` standalone in parallel** | BROKEN -- misreports by ~21 orders of magnitude. Always compute the outer residual by explicit `fvc::`. Independent of PETSc. |
| **`relativeChange` outer criterion** | REMOVED 2026-09-06. It divides by the field's deviation about its own mean, which collapses for a nearly-uniform field -- exactly `nEps_e` before ignition. |
| **ngspice bridge** | decided against |
| **GMRES (not FGMRES) on the transport split** | `DIVERGED_BREAKDOWN` from a varying operator. FGMRES tolerates it. |
| **`ROUND*01` schemes for number densities** | `libROUNDSchemes.so` registers EIGHT names; the `01` variants clamp the field to **[0,1]** -- catastrophic for a density of 1e16 m^-3, and "bounded" is exactly what a user reaches for. |

## 5. Diagnosed, with the action already chosen

* **Grubert lateral asymmetry** (`docs/design/grubert-lateral-asymmetry.md`). Mesh not
  bit-symmetric -> gradient makes a tolerance-immune `Ey` of 3.0e-07 -> `rho`'s
  6 ppm quasi-neutral cancellation amplifies it ~1000x -> `alpha(E)` feedback
  e-folds it every 1.46 ns (1.59x the mean). **Action: `NY = 1` for 1-D,
  normal refinement for 2-D/3-D. No utility, no check, no tolerance change.**
* **Grubert runaway** -> the current-source case, section 4 above.
* **SIGFPE in the Schur inner solve** -> SOLVED by a fixed inner solve
  (`richardson`, `max_it 2`, `convergence_test skip`, `hypre`) + FGMRES outside.

## 6. Open threads -- IN THE USER'S ORDER (agreed 2026-09-10)

1. **Newton-vs-Picard benchmark** on `positiveStreamer_fixedMesh` -- the
   Pasolari redo, with a `COMPARE.md`. Already parametrised; no new code. Its
   KEY measurement is the largest dt each solver survives at fixed accuracy
   with the Picard-era limiters off -- which is also the AP claim's central
   experiment, so it must precede the paper. **Leave `changeDictionary`
   authoritative** (it produced the validated 2 ns LFA/LMEA results that anchor
   accuracy; both arms share it, so the comparison stays fair) and do NOT
   resolve the BC fork here. Measure speed **per ns of simulated time at equal
   accuracy**, never per step -- Newton costs ~4-5x per step, so a per-step
   comparison is rigged for Picard.
2. **Surface charge under Newton**, in
   `tutorials/plasma/soPlasmaFoam/needleDBD` -- **dielectric cases under Newton
   are currently unguarded.** A correctness hole, not an enhancement.
3. **The AP proof / semi-implicit-Poisson paper**, with the benchmark's dt data
   in hand.

Also open, unordered: `grubert_1d_I` (running -- and `I_cond` sat at exactly
-4.5787e-10 A for 3371 steps with `g_dIdV = 0` while V ramped past -281 V,
which needs one measurement); the negative-`L` chemistry guard (cheap, fold in
when next touching `plasmaChemistryODE`); the validation-suite regeneration and
its live BC fork; the second memory-consolidation pass.

**NOT on the list: solve for `log(n)` instead of `n`.** Demoted by the user
2026-09-10 -- *"remove log(ne) from the deferred list or demote it to last and
only if i say so."* Do not propose it, schedule it, or treat it as a
prerequisite. The diagnosis is preserved at the end of
[[deferred-action-items]]; the interim workaround (restart Newton from an
established Picard state, where the clamp is inactive) works.

---

Keep this file current: it is only worth re-reading if it is true. Rule D2 --
every measurement carries its date, superseded text is marked, never deleted.
