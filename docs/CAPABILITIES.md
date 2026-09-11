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
| **Newton on a STREAMER** | `positiveStreamer_LMEA_fast` and `_fixedMesh` run under `outerSolver newton` as of 2026-09-11 (neither did before). See 3c for the four obstacles. |
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

## 3c. Newton on a STREAMER -- what it took (2026-09-11)

`positiveStreamer_LMEA_fast` and `_fixedMesh` now run under `outerSolver
newton`. Neither did before. Four obstacles, each a CLASS of bug:

* **A COEFFICIENT THAT IS EXACTLY ZERO IS NOT SAFE IF IT IS HARMONICALLY
  INTERPOLATED.** `immobile` returns zero mu/D so the ordinary drift-diffusion
  assembly reduces to ddt(n) == sources with no branching anywhere. But the
  streamer beds interpolate mobility with `harmonic` (2ab/(a+b)) and
  diffusivity with `Gauss harmonic corrected`, so an identically-zero
  coefficient is 0/0 on every face. Now 1e-30: harmonic-safe, a normal double,
  26 orders below a real ion mobility, and invisible to Picard (which reaches
  mu()/D() only through a caller that assembles the equation itself).
  **Scharfetter-Gummel divides by D as well** (its Bernoulli argument), so an
  immobile species must never be assembled with SG.
* **THE CHEMISTRY JACOBIAN CAN INVERT THE DIAGONAL'S SIGN.** `chemJacobian`
  approximates `dP_s/dn_s` as `P_s/n_s`: EXACT for the electron
  (`dP_e/dn_e = k_iz*n_gas` IS the ionisation frequency), the WRONG QUANTITY
  for a species produced by electron impact on something else, whose true
  derivative is ~0. Species O at dt = 2.5e-13: `1/dt = +4.0e12` against
  `-P/n = -1.5e13`, net **-1.1e13**. The transport block goes indefinite and
  the outer Krylov STAGNATES (5.03e2 -> 4.04e2 over 28 iterations) rather than
  diverging -- while BOTH sub-solves report themselves converged, which is the
  tell. Now clamped so the sign cannot invert, and restricted to the electron
  by default (`chemJacobianElectronOnly`): halves Krylov work on the streamer,
  measured neutral on grubert.
* **FIELDS FILLED INSIDE `plasmaTransport::solve()` GO STALE UNDER NEWTON,
  SILENTLY** -- Newton REPLACES that solve. Two were caught one after the
  other: `convectiveFlux_` (left `limitSpeciesCo` protecting nothing --
  `Co_conv (e)` read EXACTLY 0 for 1300+ steps) and `particleFlux_` (feeds
  `I_cond`, which `plasmaExternalCircuit` REGULATES ON, so the current-driven
  electrode was inoperative under Newton: `I_cond` fossilised at one Picard
  value for ~3900 steps while the gas broke down unnoticed).
  **If a field is populated inside `solve()`, assume Newton never updates it.**
* **EVERY CASE NEEDS THE `fvSchemes` CATCH-ALLS.** Newton assembles each
  species equation itself, so it interpolates mu/D and forms a flux for EVERY
  transported species. A missing entry is a FATAL lookup, not a fallback.
  `_fixedMesh` had them; `_AMR`, `_LMEA_fast`, `_LMEA_minimal` and `needleDBD`
  did not and died at the first Newton step. All fixed.

**THE COARSE BED IS THE DEBUGGING TOOL.** Every one of the above was found on
`positiveStreamer_LMEA_fast` at ~2 s per run, after hours of chasing the same
bugs on the 1.15M-cell bed at minutes per attempt. That is what it exists for
-- "COARSE BY DESIGN ... exercises the outer-loop coupling cheaply". It does
NOT resolve the streamer and must never be quoted for physics.

**PETSc SWALLOWS OpenFOAM'S STACK TRACE.** `PetscInitialize` installs its own
signal handler and `-no_signal_handler` is read too late to stop it. To locate
an FPE/SEGV inside a Newton run: `reconstructPar`, then drive it serially under
`gdb --batch -ex run -ex 'bt 40'`. That found the harmonic division in one
attempt, after two wrong guesses from reading code.

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
| **`sourceAwareScaling`** (scale a species' residual by its chemistry source instead of n/dt) | REFUTED and REMOVED 2026-09-11. `rms(chemP)/(rms(n)/dt)` measured per species: **0.833** for all seven source-dominated species, 1.6e-4 for the electron. The source NEVER exceeds ddt -- a species produced from nothing has n ~ P*dt by construction, so n/dt IS the source scale, and `max(ddt, src)` is a NO-OP. |
| **`extrapolateGuess` as a cure for block imbalance** | no effect here: 38 vs 40 converged solves over 11 steps. Left at its default; still right for the case it was built for. |
| **Ion mobility as the cause of Newton's conditioning trouble** | REFUTED. With the chemJacobian sign bug fixed, mobile and immobile arms are IDENTICAL digit for digit (18 steps, 76 KSP solves, 5.1 SNES its/step either way). The earlier apparent advantage of mobile ions WAS the bug. |
| **Ballast (voltage + series R) to reach the DC-glow operating point** | REFUTED by Almeida et al. 2017 AND by ~10 of our own arms (`seriesResistor` 1e8, `seriesRC` 1e6/1e8/1e9/5e9): at the CVC minimum dR/dI = 0, so the load line is TANGENT and selects nothing. "The fix is current control, NOT a ballast." Every arm died identically: 150/150 correctors, dt collapsing ~25,000x. Do not re-propose a ballast, including "ballast-limited ignition". |
| **Time-marching to the DC-glow operating point at all** | The same paper: COMSOL's time-dependent module fails in the same place ("convergence was lost shortly before the minimum of the CVC"), and Grubert computed his result AT STEADY STATE by FEM -- so the nanosecond avalanche that kills our runs is a transient he never traversed. Measured again 2026-09-11: a quasi-static ramp ignited cleanly at -137.7 V and still overshot to 718x setpoint, after which Newton could not step (line-search + maxIt, NO NaN). The recorded path is `ddtSchemes steadyState` + `relaxationFactors` FIRST (case settings, no new code), then a stationary driver with j-continuation -- `docs/design/stationary-solver-plan.md`. |
| **`chemJacobian`'s P/n applied to ALL species** | Exact only for the electron. For species produced by electron impact on something else, dP/dn_s ~ 0 while P/n_s is enormous -- it inverted the diagonal's SIGN and stalled the outer Krylov. Now clamped AND electron-only by default. |
| **More `-fieldsplit_phi_ksp_max_it` as the cure for the 1.15M wall** | REFUTED at 449k, 2026-09-11, 10 steps each at dt=1e-11: **2 cycles = 77 SNES / 387 Krylov / 121 s per step; 8 cycles = 73 / 373 / 138 s**. Eight buys -5% SNES and -4% Krylov for +14% wall clock -- a net LOSS, and an effect far too small to explain a total failure. The Schur complement is not starved of `A_phi^-1`. (Tested at 449k deliberately: minutes per run instead of hours -- A6.) |
| **ngspice bridge** | decided against |
| **GMRES (not FGMRES) on the transport split** | `DIVERGED_BREAKDOWN` from a varying operator. FGMRES tolerates it. |
| **`ROUND*01` schemes for number densities** | `libROUNDSchemes.so` registers EIGHT names; the `01` variants clamp the field to **[0,1]** -- catastrophic for a density of 1e16 m^-3, and "bounded" is exactly what a user reaches for. |

## 4b. GRUBERT'S DC-GLOW OPERATING POINT: SETTLED 2026-09-11

**Every route to it by time-marching is closed, and the remaining one is a
build.** Recorded here so the effort is not restarted from the top.

| route | verdict |
|---|---|
| Voltage + ballast | REFUTED. Almeida 2017: at the CVC minimum dR/dI = 0, the load line is TANGENT and selects nothing. ~10 of our arms (`seriesResistor` 1e8, `seriesRC` 1e6/1e8/1e9/5e9) all died at 150/150 correctors with dt collapsing ~25,000x. Includes "ballast-limited ignition". |
| Current control, time-marched | SOUND BUT INSUFFICIENT. All six `grubert2009_iset*` arms stalled at t ~ 9.5e-07 -- `grubert2009_iset` after **80,983 steps** with 23,575 discards. Checked 2026-09-11 per the plan's own instruction; they did not fail for a fixable reason. |
| A gentler ignition ramp | HELPS, DOES NOT SOLVE. 2.04e7 V/s (20 V per ion transit) ignited quasi-statically at -137.7 V instead of over-volting to -245 V, and cut rejections 30,547 -> 1,727. It STILL overshot to 718x setpoint, after which the source sat at its 0 V rail and Newton could not step (line-search + maxIt; NO NaN -- that was my own misreading of PETSc reason -5). |
| `ddtSchemes steadyState` as a case setting | **DOES NOT EXIST.** The solver refuses it: "the ddt term IS the diagonal... removing it leaves rows with no diagonal... SIGFPE in GaussSeidelSmoother". Added 2026-09-09, one day AFTER the stationary plan proposed it. |
| **PSEUDO-TRANSIENT** (BDF2 + `adjustTimeStep false` + fixed dt + currentSource, and NO relaxationFactors) | **WORKS TODAY, NO CODE.** Measured: **9472 consecutive converged steps, 0 failures, 4 correctors/step** against a 150 cap. Cost ~5 h per us of simulated time, and the run reported was only 9.5 ns -- far too early to say anything about the physics. So this route is STABLE BUT UNPROVEN, not closed. See `docs/design/steady-mode-spec.md`. |
| **`relaxationFactors` added by hand** | **BREAKS IT.** 9472/9472 converged -> **0/10** converged. SoPLASMA already runs ADAPTIVE AITKEN outer relaxation (`plasmaOuterRelaxation`, omega ~0.371); fixed factors fight it. Recorded in the spec as its author's own error -- and repeated by me on 2026-09-11 before reading it. The generator must never emit them. |
| A true stationary driver with j-continuation | A SOLVER PROJECT, not a case setting. And NOT reachable via relaxation: "OpenFOAM's equation relaxation cannot rescue this either -- it SCALES a diagonal that must already exist." So something must supply the species diagonal that `ddt` currently provides. |

**THE UNDERLYING REASON, which is not ours:** Almeida et al. 2017 report that
COMSOL's time-dependent module fails in the same place ("convergence was lost
shortly before the minimum of the CVC"), and that Grubert computed his result
AT STEADY STATE by FEM -- so the nanosecond avalanche that kills these runs is
**a transient he never traversed**. Time-marching is the wrong tool for this
operating point, whatever the circuit.

**NOT EVERY ROUTE IS CLOSED -- corrected 2026-09-11.** The PSEUDO-TRANSIENT
recipe above runs stably and was simply never run far enough (9.5 ns of a
problem that needs us-ms). "Settled" means the VOLTAGE/BALLAST and
case-setting-`steadyState` routes are closed, not that the operating point is
unreachable.

**CONSEQUENCE FOR BENCHMARKING.** Grubert is NOT a viable vehicle for the
Newton-vs-Picard comparison: it is blocked behind a stationary solver, which is
separate work. The streamer IS viable and is the published benchmark -- see
section 6.

## 4c. RUNS ALREADY DONE, AND WHAT THEY ANSWERED -- DO NOT RE-RUN THESE

The point of this section is that a case which has already answered its
question must never be run again to answer it a second time. Check here first.

### Newton vs mesh size, on the streamer bed (2026-09-11)

Same case family, only resolution changing. Beds built with
`make-smoke-case.sh` (`NCELL=130`, `NREFINE=n`).

| cells | bed | dt 1e-12 | dt 1e-11 | cost |
|---|---|---|---|---|
| ~40k | `positiveStreamer_LMEA_fast` | converges | -- | ~2 s/step, THE debugging bed |
| 211k | `scale_r3` | converges | converges | -- |
| 449k | `scale_r4` | converges | converges | Newton **142 s/step** vs Picard **14.6 s/step** (~10x) |
| **1.15M** | `positiveStreamer_fixedMesh` | **STUCK** | **STUCK** | burns the full 50 SNES x 200 KSP budget |

**THE 1.15M WALL IS MESH SIZE, NOT dt AND NOT THE INITIAL STATE.** Proven by
running it BOTH ways: warm-started from an established streamer (10 h CPU/rank,
still on step 1) and COLD (3 h CPU/rank, handover at t=3e-12, then stuck on the
next step). Both at 99.9% CPU -- computing, not deadlocked. 211k and 449k
converge under the identical configuration. Do not re-run 1.15M under Newton
expecting a different answer; the open question is the PRECONDITIONER, not the
case.

**The cost model matches exactly**, so this is arithmetic rather than mystery:
each JFNK Krylov iteration costs ONE full nonlinear residual assembly, as
expensive as a whole Picard step. 50 x 200 = 10,000 residual evaluations x
4.16 s = 11.6 h per timestep, against 10.06 h observed.

### The Picard dt ladder (warm-started from t=1e-9 on the 1.15M bed)

| dt | result |
|---|---|
| 1e-11 | 100 steps to t=2e-9, exit 0 |
| 2e-11 | 50 steps to t=2e-9, exit 0 |
| 5e-11 | 20 steps to t=2e-9, exit 0 |
| 1e-10 (COLD) | SIGFPE on step 1, in `GaussSeidelSmoother::smooth` |

So warm-started Picard survives to at least 5e-11. **Newton must therefore
permit >10x larger dt merely to BREAK EVEN at 449k** -- i.e. ~5e-10. That
number is the benchmark's whole question; nothing measured so far shows Newton
taking a larger step at all.

### Newton's handover, measured

Cold-started arms cannot test Newton above PICARD's own limit: handover needs
one completed Picard step, so at any dt where Picard's first step dies, the
Newton arm dies in it and Newton never runs (measured at dt=1e-10: both arms,
same SIGFPE, 1 step, handover count 0). **Warm-start from an established state**
and handover fires on step 1.

### THE HEAD-TO-HEAD dt LADDER AT 449k (2026-09-11) -- AND WHY IT ANSWERS NOTHING

Run cold on `scale_r4`, limiters ALL off (`limitSpeciesCo/limitChemistryCo/
limitDielectricRelaxationRatio false`, `adjustTimeStep false`), fixed dt, both
solvers on the same bed. Handover message VERIFIED present in every Newton arm.

| dt | Picard | Newton |
|---|---|---|
| 1e-11 | 20 steps, 2-4 correctors, exit 0 | 11 steps, 0 SNES failures |
| 2e-11 | 10 steps, 2-4 correctors, exit 0 | 10 steps, 5 of 38 SNES failed |
| 5e-11 | **SIGFPE, 1 step** | 4 steps, 3 of 48 failed, NO crash |
| 1e-10 | **SIGFPE, 1 step** | crash |

**THE LADDER CANNOT ANSWER THE BENCHMARK'S QUESTION, BY CONSTRUCTION.** Cold
start means it reaches only t = 2e-10 -- 10% of the validated 2 ns window --
with peak `n_e` still at the **1.3e13 seed**. A streamer head is ~1e20. The
dielectric relaxation time tau = eps0/(e*mu_e*n_e) is then **~1e-4 s**, so
dt/tau ~ 2e-7: the ladder ran SEVEN ORDERS OF MAGNITUDE BELOW the stiffness
constraint Newton exists to step over. Picard's own numbers confirm it --
**2 to 4 correctors** against a cap of 20, with Aitken relaxation INACTIVE.
That is a nearly-uncoupled fixed-point map, and nothing can beat 2-4 cheap
sweeps of one. At streamer-head density tau is ~1e-11 s, i.e. dt/tau ~ 2, which
is the regime the experiment was specified for.

**So the benchmark still requires the WARM-STARTED ladder from a developed
streamer** (deferred-action-items: "a metric must be tested where its answer
MUST differ"; CAPABILITIES 3b: cold arms cannot test Newton above Picard's own
first-step limit). Do not quote the table above as a verdict on JFNK.

**What it DOES establish, and it corrects an earlier claim:**
* **Picard's cold ceiling at 449k is 2e-11, not 5e-11.** At 5e-11 and 1e-10 it
  SIGFPEs on step 1. The "Newton must reach ~5e-10 to break even" bar below
  rests on WARM-started Picard at 1.15M; that is a different bed and a
  different start, and the bar must be restated against a like-for-like control.
* **Newton is NOT less robust in dt here.** At 5e-11 Picard dies outright while
  Newton completes 4 steps. Newton's cost is the problem, not its stability.
* **BASELINE CONTAMINATION, caught:** four of the seven original arms ran
  against a `soPlasmaFoam` replaced at 15:48 and were compared against 16:2x
  arms. Re-run on one binary before reading any ladder. ([[baseline-contamination]])

### NEWTON'S FAILURE MODE IS THE LINEAR SOLVE, NOT THE NONLINEAR ONE (2026-09-11)

Every single Newton failure in the ladder reads:

    Linear solve did not converge due to DIVERGED_ITS iterations 100
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0

`iterations 0` -- **SNES never took a step.** The outer FGMRES hit its
`-ksp_max_it 100` cap (hardcoded, `snesNewtonSolver.C:1530`) on the FIRST
linear solve of the step. ZERO line-search failures, ZERO NaN, ZERO nonlinear
divergence, in any arm. Successful solves converge in 3-41 iterations, most
under 12 -- so it is BIMODAL: normally trivial, occasionally straight through
the cap. That is a PRECONDITIONER cliff. Newton's nonlinear machinery was never
exercised on the failing steps, so no ladder result to date says anything about
Newton as a nonlinear method.

Three candidates, untested, cheapest first -- and (3) is NOT what
`sourceAwareScaling` refuted (that scaled per CELL within a field; this is
non-dimensionalising the field BLOCKS against each other):
1. `-ksp_max_it 100` simply too low; raising it converts a hard failure into an
   expensive success and separates "stalling" from "slow". One flag.
2. The transport split is a VARIABLE preconditioner (`rtol 1e-2`, `max_it 200`
   -- which is why FGMRES is mandatory). If that inner solve hits its own cap
   the preconditioner is noise for that iteration. NOT instrumented. Same shape
   as the chemJacobian bug: both sub-solves report converged while the outer
   stagnates.
3. No per-FIELD scaling of the matrix-free differencing. phi ~ 1e4 V,
   n ~ 1e13-1e20. `MatCreateSNESMF` picks ONE h from ||u||.

### MULTI-REGION UNDER NEWTON NEEDS NO MONOLITHIC ASSEMBLY -- MEASURED (2026-09-11)

`coupledElectricPotential::updateCoeffs()` sets `valueFraction`, `refValue`
(the mapped neighbour potential) and `refGrad` (the surface-charge jump)
**unconditionally**; `useImplicit` only ADDITIONALLY fills `source()` for
`manipulateMatrix`. So the Robin condition carries the interface physics
explicitly, and the monolithic assembly is a LINEAR-SOLVER ACCELERATION, not a
correctness requirement.

Measured on needleDBD (`scheme explicit`, which is what Newton forces), by
bumping the dielectric potential +1 V and re-evaluating the explicit gas
residual:

    interface cells (160):    d(res) rms 1.549      [base res rms 0]
    interior  cells (71585):  d(res) rms 0.00354

O(1) response at the interface; ~440x smaller one cell in, which is the
non-orthogonal correction's stencil spread. **The explicit residual SEES the
dielectric.**

**Consequence:** the Newton path evaluates each region's phi residual on its
own mesh with `fvc::` and packs them into the ragged tail --
`lduPrimitiveMeshAssembly` is not needed at all, and Newton's own outer
iteration converges the region coupling (exactly what the assembly does inside
the linear solve, but at the NONLINEAR level, which also picks up the
surface-charge and permittivity nonlinearity for free). This is generic by
construction: indifferent to the dielectric's cell fraction, and `farField`
Poisson-only regions are just another region with `laplacian(eps,phi) == 0`.

**THE TRAP THIS CREATES:** `correctBoundaryConditions()` must be called on
EVERY region before each trial-state residual evaluation, or the mapped
`refValue`/`refGrad` reflect the last PICARD state. That is the same class as
the `particleFlux_`/`convectiveFlux_` staleness (3c) -- a field populated in
the Picard path going stale under Newton.

Assembled ordering, from the live probe (gas FIRST, confirmed not inferred):

    cellOffsets = 0  71745  77721     gas 71745 cells, dielectric 5976

**A diagnostic note worth keeping.** The FIRST version of this probe compared
interface-cell residual against interior-cell residual and read EXACTLY 0 at
the interface -- equally consistent with "coupling perfect" and "nothing
happening there". It separated no causes. The perturbation test has only one
possible reading. ([[diagnostics-must-separate-causes]],
[[silent-diagnostic-trap]])

### MULTI-REGION UNDER NEWTON: DONE 2026-09-11 (commit 7df1cbe)

`outerSolver newton` now accepts `multiRegionPoisson`. The dielectric and
far-field regions are packed into the RAGGED TAIL of the DOF layout -- phi
lives in every region, the species only in the gas.

**No `lduPrimitiveMeshAssembly` is involved.** See the measurement above: the
explicit `fvc::` residual already sees the neighbouring region through
`coupledElectricPotential`'s Robin condition, and Newton's own outer iteration
converges the coupling -- which is what the monolithic assembly does inside the
LINEAR solve, but at the nonlinear level. So it is GENERIC: the dielectric's
cell fraction never enters, and a `farField` Poisson-only region is just
another region with `laplacian(eps,phi) == 0`.

| gate | result |
|---|---|
| single-region bit-identical | **PASS** -- 62 fields identical, 0 differing, stripped log identical |
| needleDBD past the guard | **PASS** (it used to die at "supports only singleRegionPoisson") |
| handover on a dielectric case | **PASS** -- t=2e-12 |
| ragged tail wired | **PASS** -- `Poisson tail: 1 extra region(s), 5976 rows`, matching cellOffsets 0/71745/77721 |
| a real Newton step | **PASS** -- first linear solve 2 iterations, \|\|F\|\| 1037 -> 245 |

**THREE DEFECTS THE GATE CAUGHT, none visible without running it:**

* **Empty Pmat tail rows.** Extending the phi split over the tail without
  assembling anything there left 5976 rows with a ZERO DIAGONAL, and the first
  `PCApply` raised an FPE. The residual had already evaluated cleanly to
  `0 SNES Function norm 1.037388548231e+03`, so nothing upstream flagged it.
  Fixed by `blockMatrixCOO::addFvMatrixTail()`.
* **A FIFTH `fvSchemes` catch-all: `"snGrad\(n_.*\)"`.** Newton refreshes
  `particleFlux_` with an explicit `fvc::snGrad(n_s)`, so EVERY species needs a
  snGrad scheme. It only shows up on a case whose species set is bigger than
  the electron: needleDBD died after 545 residual evaluations on
  `Entry 'snGrad(n_N2p)' not found`. The streamer beds needed only three
  (`interpolate(mu_*)`, `div(phi_*,n_*)`, `laplacian(D_*,n_*)`).
* **PARALLEL IS REFUSED, loudly.** The tail's processor-interface coupling is
  not assembled; a FatalError says so rather than building a preconditioner
  from disconnected subdomains.

**A HANDOVER PROBLEM THIS EXPOSED, still open.** `anySpeciesOnFloor()` gates
the Picard->Newton handover on EVERY species having max > its floor. In
needleDBD all 13 species start exactly AT their floor (n_e 1e11 = the default
floor; N2p 7.9e10, O2p 2.1e10, Om/O2m 1e5 = theirs; and eight more at
`uniform 0`, which the clamp lifts to the floor). So Newton is unreachable
until ignition -- thousands of steps at ~3 s each. The gate above seeds all 13
above their floors, which is MACHINERY, NOT PHYSICS and must never be quoted
for DBD results. This matters for the automatic Picard/Newton switch: the
handover is also ONE-WAY and PERMANENT by construction.

### Surface charge under Newton -- ANSWERED

Works. `thinDielectricOnElectrode` (pmma, 100 um, Vb=0) on the coarse streamer
bed ran clean under `outerSolver newton`: surfCharge field, `chargingSurface`
wall flux and the Robin potential condition all exercised, no solver work
needed. **The only blocker for DBD cases is `multiRegionPoisson`**, which Newton
refuses outright ("supports only singleRegionPoisson so far"). needleDBD hits it
because its `regionProperties` declares a `dielectric` MESH REGION; the thin
roles themselves are single-region.

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

Also open, unordered:

* **`grubert_1d_I`, relaunched COLD 2026-09-11** -- and this is now the FIRST
  genuine test of current control, because the regulator has never actually
  worked under Newton (see 3c, `particleFlux_`). Its earlier state was produced
  by ~3900 steps of UNREGULATED ramping and is not trustworthy as physics.
* The negative-`L` chemistry guard (cheap; fold in when next touching
  `plasmaChemistryODE`).
* The validation-suite regeneration and its live BC fork
  (`changeDictionary` vs layer 1 disagree on electrode BCs and the far-field
  `n_e` seed).
* The second memory-consolidation pass.

**NOT on the list: solve for `log(n)` instead of `n`.** Demoted by the user
2026-09-10 -- *"remove log(ne) from the deferred list or demote it to last and
only if i say so."* Do not propose it, schedule it, or treat it as a
prerequisite. The diagnosis is preserved at the end of
[[deferred-action-items]]; the interim workaround (restart Newton from an
established Picard state, where the clamp is inactive) works.

---

Keep this file current: it is only worth re-reading if it is true. Rule D2 --
every measurement carries its date, superseded text is marked, never deleted.
