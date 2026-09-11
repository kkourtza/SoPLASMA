# build-and-test

## Summary
The exact build/test command surface for both modules, verified by reading every script and by measuring what was safe to run read-only. SoPlasma (OpenFOAM) has TWO build entry points that are NOT equivalent: `/home/kkourtza/soplasma-scratch/build-all.sh` (serial, dependency-ordered, per-component logs in `buildlogs/`, a coverage gate that fails on any un-listed application, and a `BUILD-FAILED` exit) versus the older `./Allwmake`, which builds only 3 solvers + 1 utility and NONE of the 11 test beds. The project's own `etc/bashrc` must be sourced AFTER OpenFOAM's, because it is the only place `PETSC_DIR`/`PETSC_ARCH` are set; without it `libplasmaNewtonSolverPETSc` silently goes stale. Six guard/run scripts (`check-no-running-solvers.sh`, `run-guarded.sh`, `run-long.sh`, `watchdog.sh`, `rerun.sh`, `check-run.sh`, `make-smoke-case.sh`) each encode one measured disaster (a 9.7 GB runaway log, a SIGBUS rebuild-over-a-running-solver, a mesh rebuilt for no reason). SoEEDF is plain CMake: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)` producing 10 binaries (there is NO `validateBiagi`). Fast modes already exist for roughly half the surface — `make-smoke-case.sh` (NCELL/NREFINE), `positiveStreamer_LMEA_minimal` (75 s vs ~1 h), `plasmaChemistry0D -endTime`, `genMechTables --EN/--grid`, `mesh_sweep.sh NSTEPS/NRELAX` — and are missing entirely for `build-all.sh`, `testPhysics`, `validateAir*`, `verification/fluxScheme1D` and the three unit beds (which need none: they are already 0.03 s warm).

## Facts

### `./Allwmake` and `./build-all.sh` are NOT interchangeable: Allwmake's application step builds only singleRegionElectrostaticFoam, multiRegionElectrostaticFoam, soPlasmaFoam and foamPlasmaCreateSpeciesFields. None of the 11 testXxx beds, and not plasmaChemistry0D, are built by it.
**Evidence:** /home/kkourtza/soplasma-scratch/src/applications/Allwmake lines 16-22 (four `wmake` lines only) vs build-all.sh BUILD_DIRS which lists 12 utilities

**Rule:** Always rebuild with `/home/kkourtza/soplasma-scratch/build-all.sh`. Treat `./Allwmake` as legacy; if you use it, the test beds you then run are STALE binaries.


### build-all.sh carries a coverage gate: every directory under `src/applications` with a `Make/` must appear in BUILD_DIRS or SKIP_DIRS, otherwise it prints `UNCOVERED <dir>` and `BUILD-INCOMPLETE` and exits 1. It was added because testAitken was absent for nine days: binary built 2026-08-21 19:12, libplasmaNumerics.so rebuilt 20:32 against a changed aitkenRelaxation.H, and the ABI mismatch segfaulted the unit test mid-run (free() on 0x3ff0000000000000, the bit pattern of 1.0). Both electrostatic solvers had NEVER been built at all.
**Evidence:** /home/kkourtza/soplasma-scratch/build-all.sh, the `---- coverage --` block and its comment; find pattern `find src/applications -mindepth 2 -maxdepth 3 -type d -name Make`

**Rule:** When adding any new application, add its directory to BUILD_DIRS (or to SKIP_DIRS with a written reason) in build-all.sh in the same commit.


### The coverage gate scans ONLY `src/applications`, so `/home/kkourtza/soplasma-scratch/applications/test/testChemistryBackends` — a second, separate applications tree — is invisible to it and is in no build list. Its installed binary is dated Aug 12 12:36 while libplasmaChemistry.so is Sep 5 17:55 and libplasmaTransport.so is Sep 11 17:16: exactly the 30-day ABI-staleness pattern that produced the testAitken segfault.
**Evidence:** ls -la $FOAM_USER_APPBIN: `testChemistryBackends  61368  Aug 12 12:36`; ls -la $FOAM_USER_LIBBIN: `libplasmaChemistry.so Sep 5 17:55`, `libplasmaTransport.so Sep 11 17:16`; build-all.sh find root is `src/applications`

**Rule:** Do not trust a `testChemistryBackends` result until you rebuild it: `source etc/bashrc && wmake applications/test/testChemistryBackends`. Better: add it to BUILD_DIRS and widen the coverage find to include the top-level `applications/` tree.


### build-all.sh sources OpenFOAM's bashrc then the project's own `$SoPLASMA_ETC/bashrc`, which is the ONLY place PETSC_DIR=$SoPLASMA/ThirdParty/petsc-3.24.0 and PETSC_ARCH=DPInt32 are set (OpenFOAM's own bashrc points PETSC_DIR at a nonexistent system path on this apt install, and /usr/lib/openfoam/openfoam2412/ThirdParty is a stub FILE, not a directory). Without it libplasmaNewtonSolverPETSc cannot build and goes stale silently — which segfaulted against a changed plasmaTransportModel ABI on 2026-09-10.
**Evidence:** /home/kkourtza/soplasma-scratch/build-all.sh lines 4-13; /home/kkourtza/soplasma-scratch/etc/bashrc last 4 lines; ThirdParty/petsc-3.24.0/DPInt32/lib exists (libHYPRE, libfblas, libflapack present)

**Rule:** In any shell that will build or run SoPlasma: `source /usr/lib/openfoam/openfoam2412/etc/bashrc` FIRST, then `source /home/kkourtza/soplasma-scratch/etc/bashrc`. Order matters — the second overrides PETSC_DIR.


### `wmake src/numerics` builds libplasmaNumerics and NOT libplasmaNewtonSolverPETSc, exiting 0 with the latter untouched. It is a separate target that must be named explicitly: `src/numerics/newtonSolverPETSc`.
**Evidence:** /home/kkourtza/soplasma-scratch/build-all.sh, comment block immediately above the `src/numerics/newtonSolverPETSc` BUILD_DIRS entry

**Rule:** Never conclude 'the numerics are rebuilt' from `wmake src/numerics`. Check `ls -la $FOAM_USER_LIBBIN/libplasmaNewtonSolverPETSc.so`.


### `Allwmake` can exit 0 having silently NOT rebuilt a directory whose .C mtime is newer than its .o. Measured 2026-09-08: after editing localEnergyEnergyModel.H/.C, `FORCE=1 ./Allwmake` exited 0 with no output and libplasmaEnergy.so's timestamp was UNCHANGED, while libplasmaProperties.so and libplasmaTransport.so in the same pass did rebuild. Root cause never identified.
**Evidence:** memory /home/kkourtza/.claude/projects/-home-kkourtza-Projects-SoEEDF/memory/allwmake-can-silently-skip-a-changed-directory.md

**Rule:** After any build, compare `.so` mtimes under $FOAM_USER_LIBBIN against the source .C mtimes you edited. If one is stale: `cd <that dir> && wclean && wmake libso .`, then one more full pass.


### A partial `wmake` of only the libraries you edited produces a phantom SIGFPE. Measured 2026-09-07: plasmaSpecies.C and plasmaEnergy were rebuilt, plasmaTransport (which LINKS plasmaSpecies) was not; the stale libplasmaTransport.so raised SIGFPE in mechanismSourceTerms on the first chemistry call, was misdiagnosed as a pre-existing Scharfetter-Gummel defect, and the 'control experiment' ran against the same stale build so it confirmed nothing. Rebuilding plasmaTransport made the crash vanish and SG then agreed with `standard` to 2.6% in n_e.
**Evidence:** memory rebuild-with-allwmake-not-piecemeal.md; and openfoam-abi-partial-rebuild.md for the header-layout variant (bad_alloc in soPlasmaFoam, segfault in plasmaTimeControl::adjustDeltaT inside a library never edited)

**Rule:** Partial `wmake` is legitimate ONLY for a .C-only edit inside a leaf library you are still iterating on, and its result must never be used to validate a behaviour change. Any .H member-layout edit, and any conclusion drawn from a run, requires a full ./build-all.sh first.


### `BUILD-COMPLETE` used to be printed unconditionally even when components reported FAIL (measured 2026-08-21: BUILD-COMPLETE with zero `error:` lines while four components FAILed and libplasmaTools.so did not exist). build-all.sh now exits 1 with `BUILD-FAILED: at least one component did not build. See buildlogs/.` — but a component's own log can contain `Error 1` or `cannot find -l...` which never appear as `error:`.
**Evidence:** /home/kkourtza/soplasma-scratch/build-all.sh trailing `if [ "$failed" -ne 0 ]` block and its comment; memory pre-run-checklist.md item 2

**Rule:** Judge a build by three things, never one: last line is `BUILD-COMPLETE`, `grep -c '^FAIL' ` of the output is 0, and `ls -la $FOAM_USER_LIBBIN/<the .so you changed>` is newly stamped. NEVER judge it by the exit code of a grep pipeline — `./build-all.sh | grep BUILD-COMPLETE` returns grep's status.


### Both ./Allwmake and ./build-all.sh call ./check-no-running-solvers.sh first and REFUSE to build (exit 1, no BUILD-COMPLETE) while any of soPlasmaFoam, plasmaChemistry0D, plasmaCreateSpeciesFields, foamPlasmaCreateSpeciesFields, testWallLoss, testVibRelax or testChemistryBackends is running. It matches by NAME with `pgrep -x` because `pkill -f soPlasmaFoam` does not reliably match `soPlasmaFoam -parallel` and DOES kill the calling shell.
**Evidence:** /home/kkourtza/soplasma-scratch/check-no-running-solvers.sh, `_apps` list and the `pgrep -x` loop; memory kill-runs-before-launching.md (a stale needleDBD run made build-all return 1 with no compile error at all)

**Rule:** Before any build: `pgrep -af soPlasmaFoam` and kill what you no longer need, as its own simple command. Override only when you accept losing those runs: `FORCE=1 ./build-all.sh`. Note `pgrep -x` fails on names longer than 15 chars (plasmaChemistry0D) — use `pgrep -f` to LOOK, `pkill -x` to kill.


### On this WSL2 box the FIRST invocation of an OpenFOAM/SoPlasma binary after a build costs 24-45 s of cold dynamic-loader page-in; the same binary re-run immediately costs 0.00-0.03 s. Measured today: testWallLoss 24.20 s cold / 0.00 s warm; testWallFlux 44.62 s cold / 0.03 s warm (libplasmaTransport.so alone is 18 MB).
**Evidence:** /usr/bin/time -f 'WALL %e s' testWallLoss → WALL 24.20 s then WALL 0.00 s; testWallFlux → WALL 44.62 s then WALL 0.03 s; ls -la $FOAM_USER_LIBBIN/libplasmaTransport.so = 18018880 bytes

**Rule:** Do not time a unit bed on its first run after a build and do not conclude it is slow. Run it twice; quote the second. Budget ~45 s once per build for cache warm-up.


### Five unit beds need NO case directory, no mesh and no dictionary (they include neither fvCFD.H nor setRootCase*.H): testAitken, testCoulombHeating, testVibRelax, testWallFlux, testWallLoss. Six DO need a case directory: testDischargeCurrent (createMeshes, multi-region), testEmission (setRootCase+createTime, 'Needs a case with a mesh'), testFluxScheme, testPoissonSymmetry, testSnesJFNK, testSnesJFNK2Field. plasmaChemistry0D needs none — it constructs `argList args(argc, argv, false, false, false)`.
**Evidence:** grep of `#include "(fvCFD.H|argList.H|setRootCase*|createMesh*)"` across src/applications/utilities/*/*.C; plasmaChemistry0D.C:236

**Rule:** Run the five case-free beds from anywhere (e.g. /tmp) with only the OpenFOAM+SoPLASMA env sourced. For the other six, cd into a real case first — verification/fluxScheme1D, verification/testSnesJFNK, verification/testSnesJFNK2Field are the shipped ones.


### testAitken ALWAYS returns 0 and prints no PASS/FAIL — its output ends `=== end ===`. It is a reporting bed whose numbers a human reads (separation ratios, `(LEAD)` / `(no useful lead)`), not a gate. Every other unit bed returns nonzero on failure: testWallFlux prints `ALL PASS: 35 checks, 0 failed`, testWallLoss/testVibRelax print `all checks passed`, testCoulombHeating prints `5/5 checks passed.`, testDischargeCurrent prints `PASSED`.
**Evidence:** tail of testAitken.C (`Info<< nl << "=== end ==="; return 0;`); testWallFlux.C:658 printf; testWallLoss.C:252; testVibRelax.C:173; testCoulombHeating.C final printf; testDischargeCurrent.C:185

**Rule:** Never put testAitken in a `&&` chain as a pass gate. For the others `testWallFlux && testWallLoss && testVibRelax && testCoulombHeating` is a valid gate.


### `testFluxScheme` is likewise not self-judging: it prints one `RESULT scheme=... N=... h=... Pe=... PeGrid=... L2=... Linf=... L2rel=... nCorr=...` line and returns 0. The verdict is produced by the caller — verification/fluxScheme1D/Allrun collects 504 RESULT rows (6 domain-Peclet x 6 meshes x 14 scheme variants) into results.txt and runs `python3 report.py`; a baseline sits in results.baseline.txt.
**Evidence:** tail of testFluxScheme.C printf; verification/fluxScheme1D/Allrun loops `for PE in 0.1 1 10 100 1000 10000` x `for N in 20 40 80 160 320 640` x 12 fvSchemes entries + 2 CFS schemes; `wc -l results.txt` = 504; `ls logs | wc -l` = 541

**Rule:** To check a flux-scheme change, diff results.txt against results.baseline.txt. To debug, run one point by hand: `testFluxScheme -v 100 -D 1 -S 1 -a 1 -b 2 -scheme ScharfetterGummel` in verification/fluxScheme1D (sub-second).


### The canonical Python interpreter is `~/ct-env/bin/python` (Python 3.12.3, numpy 2.5.1, matplotlib 3.11.1, cantera 3.2.0, fluidfoam 0.2.9, pytest 9.1.1, pyyaml 6.0.3, gmsh 4.15.2). But bare `python3` on PATH resolves to `/home/kkourtza/ct-plasma-env/bin/python3`, a DIFFERENT venv that has numpy/scipy/matplotlib/fluidfoam but NO cantera, NO pytest and NO pyyaml. Several shipped scripts call bare `python3` (verification/fluxScheme1D/Allrun's `python3 report.py`, validation/cut200k.sh, validation/mesh_sweep.sh's inline `python3 -c`, tutorials' `python3 initGaussianSeed.py`) and two call `/usr/bin/python3` explicitly (the electrostatics Allrun-sweep scripts).
**Evidence:** `which python3` → /home/kkourtza/ct-plasma-env/bin/python3; ls of both site-packages dirs; docs/CAPABILITIES.md:151 'Python is ~/ct-env/bin/python, not the system python'; grep of ct-env in tools/make_mesh_arm.sh:30,34,36 and validation/mesh_sweep.sh:106,110,112

**Rule:** Anything touching cantera, mechc, pyyaml or pytest MUST use the absolute path `~/ct-env/bin/python`. Do not 'fix' a script's bare `python3` to `~/ct-env/bin/python` without checking it — the seeding/plot scripts work under ct-plasma-env and the electrostatics sweeps deliberately use /usr/bin/python3 (stdlib only).


### SoEEDF builds with plain CMake and produces exactly 10 targets: libSoEEDF.so, testBoltzmann, testPhysics, testMechanism, validateAir, validateAirSST, validateArgon, genValidationData, genMechTables, lxcatdump. There is NO validateBiagi anywhere in the tree.
**Evidence:** /home/kkourtza/Projects/SoEEDF/CMakeLists.txt (10 add_library/add_executable blocks); ls -la build/ shows all 10 artefacts dated 2026-08-20 17:47; `grep -rn validateBiagi` returns nothing

**Rule:** Use the real list. If a doc or agent mentions validateBiagi, it is wrong.


### The existing SoEEDF build/ cache was configured from `/home/kkourtza/Projects/BoltzmannSolver`, which is a symlink to Projects/SoEEDF created 2026-08-12. Deleting that symlink breaks incremental builds in build/.
**Evidence:** build/CMakeCache.txt: `CMAKE_HOME_DIRECTORY:INTERNAL=/home/kkourtza/Projects/BoltzmannSolver`, `SoEEDF_SOURCE_DIR:STATIC=/home/kkourtza/Projects/BoltzmannSolver`; `ls -ld` shows BoltzmannSolver -> /home/kkourtza/Projects/SoEEDF

**Rule:** Detector: `grep CMAKE_HOME_DIRECTORY /home/kkourtza/Projects/SoEEDF/build/CMakeCache.txt`. If you ever need a clean configure, `rm -rf build` first — reconfiguring in place from the new path will error on the source-dir mismatch.


### SoEEDF's build is currently UP TO DATE: every binary is stamped 2026-08-20 17:47 and the newest source file (src/MechTables.C) is 2026-08-20 17:47. So a `cmake --build build -j` today is a no-op, not a rebuild.
**Evidence:** ls -la build/ vs `find src include test tools -name '*.C' -o -name '*.H' | xargs ls -la | sort -k6,7 | tail`

**Rule:** Do not report 'rebuilt SoEEDF' from a no-op build. Check a binary's mtime advanced.


### Every SoEEDF test/validate binary loads cross sections by the RELATIVE path `../data/lxcat/cross-sections/...` and validateAir reads `../test/data/Bolsig_N2_O2_noCoulomb.dat`, so they MUST be run with cwd = the build directory. testPhysics says so itself on failure: `(run from the build/ directory)`.
**Evidence:** testPhysics.C:42,50,154,220,800,1084; validateAir.C:127,148; validateArgon.C:63,114; testBoltzmann.C:52

**Rule:** Always `cd /home/kkourtza/Projects/SoEEDF/build` before running any of them; `./build/testPhysics` from the repo root fails on data paths.


### Three SoEEDF validators WRITE a CSV into the current directory: validateAir → air_validation.csv, validateAirSST → air_sst_validation.csv, validateArgon → ar_kion_vs_eps.csv; testPhysics writes a temporary modified LXCat file. genValidationData writes the full validation data set.
**Evidence:** grep of ofstream: validateAir.C:168, testPhysics.C:1088, validateArgon.C:148, validateAirSST.C:193

**Rule:** Expect and ignore these artefacts in build/. Do not run the validators inside a case directory you care about.


### `positiveStreamer_LMEA_minimal` and `positiveStreamer_LMEA_fast` are the project's designed fast beds: 75 s wall, 49 steps to 5e-10 s on the bare 130x130 block (16 900 cells, ~1.4 s/step), measured 75.64 s vs 75.68 s for the two. The production streamer is 1 146 466 cells at ~70 s PER TIME STEP — about 45x slower per step. The two fast beds produce byte-identical fields at the final time across all 45 written fields (verified 2026-09-01), which is itself a regression test of the LMEA defaults.
**Evidence:** tutorials/plasma/soPlasmaFoam/positiveStreamer/positiveStreamer_LMEA_fast/tutorial_info.md (wall-clock table: 96 um/16 900 cells/75 s; 24 um/81 640/435 s; 3 um/1 146 466/~1 h) and positiveStreamer_LMEA_minimal/tutorial_info.md

**Rule:** Ask outer-loop/relaxation/step-rejection questions on positiveStreamer_LMEA_fast. NEVER quote its field magnitudes, propagation speed or head structure as physics — meanE range is 1.398-1.478 eV (ratio 1.1) against production's 0.411-7.759 eV (ratio 18.9), so nothing threshold-dependent is testable there, and the 24 um intermediate bed does not help.


### `make-smoke-case.sh` is the generic seconds-scale bed maker: it copies positiveStreamer_fixedMesh, drops all five refineMesh passes, coarsens blockMesh to NCELL^2 and sets endTime/writeInterval to 2e-12 (two timesteps). It advertises ~68x fewer cells than the 1.1M-cell case which takes ~90 s per check. Knobs are environment variables: NCELL (default 40), NREFINE (default 0, up to 5), CASE (default ~/streamer-smoke).
**Evidence:** /home/kkourtza/soplasma-scratch/make-smoke-case.sh header and the three `sed -i` blocks

**Rule:** For a code change that is not about the answer: `NCELL=40 NREFINE=0 CASE=$HOME/streamer-smoke ./make-smoke-case.sh` then `CASE=$HOME/streamer-smoke ./rerun.sh`. Its own header states what it is NOT for: streamer velocity, propagation, grid convergence.


### `rerun.sh` exists because Allrun-serial rebuilds the mesh through five refinement passes and re-seeds every time, which costs minutes and cannot change between solver rebuilds. It restores 0/ from `.snapshot0`, deletes all time dirs via `foamListTimes -rm -withZero`, and runs soPlasmaFoam alone under `timeout ${TIMEOUT:-1800}`.
**Evidence:** /home/kkourtza/soplasma-scratch/rerun.sh header and body

**Rule:** Use `CASE=<dir> TIMEOUT=300 ./rerun.sh` for every code-change re-test; use Allrun-serial only when the MESH or CASE SETUP (species list, seed, BCs) changed. It refuses with `no .snapshot0 -- run Allrun-serial once first`.


### `run-guarded.sh`, `run-long.sh` and `watchdog.sh` all exist because a runaway warning loop once wrote 9.7 GB and 91M lines before a single timestep completed and took the machine down. run-guarded caps via `head -c $((MAXMB*1024*1024))` (solver then dies of SIGPIPE) plus `timeout $MAXSEC`; defaults 50 MB and 600 s. It prints `RUN-EXIT=<rc>  log=<size>  lines=<n>`.
**Evidence:** /home/kkourtza/soplasma-scratch/run-guarded.sh header + last two lines; run-long.sh (200 MB cap); watchdog.sh (`pkill -9 -x soPlasmaFoam` above the cap)

**Rule:** Never launch an unattended solver bare. `./run-guarded.sh 20 300` for a bounded debug run; `./watchdog.sh <log> 50 &` alongside anything long. Note run-guarded.sh has the case hardcoded to $HOME/streamer-case — edit-free use requires that path.


### `tools/run_electrostatics_tests.sh` runs the whole electrostatics-only suite (2 analytic sweeps + 3 series-stack cases + 1 completion-only case) and 'the whole suite is seconds'. Exit status is 0 only if EVERY case passed; it prints `  N ok, M failed`. It exists because every one of these cases was DEAD for three weeks: fe827ae (2026-08-11) made backgroundDensity fatal in the Poisson coeffs but left `reducedE = Emag/backgroundDensity` unguarded in singleRegionPoisson, whose only setter lives in plasmaSpecies which an electrostatics-only solver does not have, so every case died with SIGFPE. Found by accident 2026-09-03, fixed in 82fc587.
**Evidence:** /home/kkourtza/soplasma-scratch/tools/run_electrostatics_tests.sh header (WHY THIS EXISTS) and the `exit $((fail > 0))` tail

**Rule:** Run `tools/run_electrostatics_tests.sh` after ANY change to the Poisson/electrostatics path or to a boundary condition. It requires WM_PROJECT_DIR set (it exits 2 with 'ERROR: OpenFOAM is not sourced' otherwise).


### The mesh-independence sweep `validation/mesh_sweep.sh` has measured per-arm cost: 2000 cells = 29.83 s wall / 104 MB maxRSS, 20 240 cells = 441.12 s / 296 MB, 200 000 cells had to be cut short by a separate script. Its knobs are env vars NSTEPS (default 20 measured), NRELAX (default 50 discarded) and DT (default 1e-12); dt is scaled as DT*400/NX so the Courant number is constant across arms — a fixed dt instead was measured to be the whole effect (Co_conv 6.50 at 2000 cells vs 20.62 at 20 240, ratio 3.17 = exactly the 3.16x refinement).
**Evidence:** validation/mesh_sweep.log lines 5-16; validation/mesh_sweep.sh header and GRIDS=("400 5" "1265 16" "4000 50")

**Rule:** Fast mode: `NSTEPS=5 NRELAX=10 ./mesh_sweep.sh`. The GRIDS array is hardcoded with no env override — to drop the 200k arm you must edit it, or use the validation/cut200k.sh pattern (never edit a RUNNING bash script; bash reads it incrementally).


### `plasmaChemistry0D` is the documented fastest end-to-end installation check: 'it runs in seconds and its self-tests are printed rather than inferred'. A healthy run prints `plasmaChemistry: 48 reactions (22 electron-impact, 26 heavy), 15 species, ODE solver rodas23`, `charge residual of the RHS at t=0: 3.50314e-16` and `Jacobian vs finite difference: worst relative error 5.00114e-07  at d(O2m)/d(O2)`. Adding `-manifest <mech.json>` also prints `dynamic EEDF: N Boltzmann solves over the run, 0 unconverged`.
**Evidence:** /home/kkourtza/Projects/SoEEDF/docs/installation.md lines 258-296

**Rule:** Use `-endTime` as the fast knob (20e-9 in the doc). The two residuals are the pass criterion — both must be machine-scale; they are computed, not asserted.


### `testChemistryBackends` runs WITHOUT a mesh, costs about a second, and prints `Backends agree.` with differences at the 1e-15 level. It needs a directory containing `constant/*.foam` and `constant/*.heavy.yaml`, and accepts `[-mechanism <file>] [-Tgas 300] [-dt 1e-9]`.
**Evidence:** applications/test/testChemistryBackends/testChemistryBackends.C Description+Usage; SoEEDF docs/installation.md Step 5 'Expect `Backends agree.` ... It runs in about a second and needs no mesh.'

**Rule:** Use `-dt` as the fast knob. Rebuild it first (see the staleness fact) — its binary is a month older than libplasmaChemistry.so.


### An `EXIT=0` from any Allrun-serial means only that THE SCRIPT finished: OpenFOAM's `runApplication` does not propagate the solver's status. The newer scripts fix this individually — positiveStreamer_LMEA_minimal's last lines capture `rc=$?` from soPlasmaFoam, print `soPlasmaFoam exit $rc, <n> steps` and `exit $rc`; needleDBD uses `die()` after every step.
**Evidence:** memory unit-bed-beats-cfd-case.md item 6; positiveStreamer_LMEA_minimal/Allrun-serial tail; needleDBD/Allrun-serial die() pattern

**Rule:** After any Allrun, check the LOG: `grep -c '^Time = ' log.soPlasmaFoam` and `grep -E 'FOAM FATAL|sigFpe|sigSegv' log.soPlasmaFoam`, not the script's exit code.


### Seeding failure used to be silent and is now fatal: `python3 initGaussianSeed.py` was wrapped in `|| true`, so a missing fluidfoam left a quiescent gas in a Laplace field — 2000 steps, 0 fatal errors, n_e pinned at its 1e13 background, no streamer anywhere. Measured 2026-09-11. The script now exits with 'refusing to run a streamer case with no streamer'.
**Evidence:** tutorials/.../positiveStreamer_LMEA_minimal/Allrun-serial, the initGaussianSeed block comment

**Rule:** If you copy a streamer case, keep the seeding step fatal. Detector on any streamer log: peak n_e still equal to the 1e13 background means the seed never ran.


### The installed binary name differs from the source directory name for two utilities: `src/applications/utilities/foamPlasmaCreateSpeciesFields` installs as `plasmaCreateSpeciesFields`, and `foamPlasmaSetupBoundaries` installs as `plasmaSetupBoundaries`.
**Evidence:** Make/files of each: `EXE = $(FOAM_USER_APPBIN)/plasmaCreateSpeciesFields`, `.../plasmaSetupBoundaries`; ls of $FOAM_USER_APPBIN confirms both names

**Rule:** Grep for both spellings when tracing a utility; check-no-running-solvers.sh already lists both.


### A NEW header is invisible to the build until `wmakeLnInclude -update src/<component>`, and the consumer's Make/options may still lack the include path. Cost four build iterations on 2026-09-02; the error names the including header, not the missing symlink (`fatal error: pressureUnits.H: No such file or directory`).
**Evidence:** memory new-source-file-lninclude.md

**Rule:** After adding a source file: `wmakeLnInclude -update src/<component>` and `grep lnInclude <consumer>/Make/options` BEFORE building.


### The two JFNK utilities are deliberately in build-all.sh's SKIP_DIRS, with the reason recorded: they need PETSC_DIR/PETSC_ARCH from the project's etc/bashrc, which build-all.sh does source but which the SKIP comment predates. The documented way to build them is `source etc/bashrc && wmake src/applications/utilities/testSnesJFNK`. Their design constraint is a hard split: testSnesJFNK.C includes fvCFD.H and NEVER a PETSc header; snesBridge.C includes petscsnes.h and NEVER an OpenFOAM header; they communicate only through int/double/void* in snesBridge.H — PETSc's petscmath.h and OpenFOAM's operator overloads are ambiguous in one translation unit.
**Evidence:** build-all.sh SKIP_DIRS comment block; docs/design/newton-outer-solver-design.md lines ~242-285

**Rule:** Build the JFNK beds by hand as above. Any future PETSc-touching SoPlasma code must follow the same two-translation-unit split. PetscInitialize does not strip its own flags, so pass PETSc options via the PETSC_OPTIONS environment variable, not the command line (argList rejects -snes_monitor).


### smoke/{native,native-dt2,native-stiff,native-stiff-dt2,cantera} are NOT tracked by git, so they never appear in a diff and rot behind every key rename — five broken beds were found in one sitting on 2026-09-02, each blocking the next. Their current recorded state is failure: smoke/STATUS reads `native exit=1 / cantera exit=1 / native-stiff exit=1 / ALL DONE`.
**Evidence:** memory smoke-beds-lag-refactors.md; cat /home/kkourtza/soplasma-scratch/smoke/STATUS

**Rule:** Do not treat smoke/ as a working smoke suite; it is currently all-failing and untracked. After any dictionary-key rename, sweep smoke/ explicitly by path — `git status` will not show it.


### A full build-all.sh is a 15-20 minute job per the project's own rule 10 ('background anything over 1-2 min'); the last incremental pass measured 3 minutes end to end (buildlogs timestamps 17:15 → 17:18 on 2026-09-11, all components OK, no `Error 1` or `cannot find -l` in any log).
**Evidence:** memory run-long-cases-in-background.md ('blocking on a 15-20 minute build-all.sh'); ls -la /home/kkourtza/soplasma-scratch/buildlogs/ (constants.log 17:15 … tools.log 17:18); grep -lE 'Error 1|error:|cannot find -l' buildlogs/*.log returns nothing

**Rule:** Always launch build-all.sh in the background and do other work; poll buildlogs/ rather than blocking.


### `tools/checkConfigReference.py` is a read-only doc/code consistency gate with three verdicts (MISSING READER, STALE DEFAULT, UNDOCUMENTED) and real exit semantics: 0 clean, 1 on MISSING READER or STALE DEFAULT, and with --strict also on UNDOCUMENTED. It found three drifts that had survived review, including `ePotentialControls { nonCoupledResidualControl }` which NOTHING has ever read.
**Evidence:** tools/checkConfigReference.py docstring 'HOW TO RUN IT' block and argparse at lines 340-349

**Rule:** Run `tools/checkConfigReference.py` before committing any docs/reference change and after adding/renaming/re-defaulting any option in src/. Use `--case <dir>` for configuration/config variables — checking those directly reports every one as a false MISSING READER.


### `genMechTables` carries the only real fast knobs on the SoEEDF table-generation path: `--EN min:max:n` (default 0.1:2000:200), `--grid N` (solver energy grid, default 200), `--eedf-points N`, `--no-eedf`. The Cantera cross-validation deliberately uses `--grid 3200` because at 200 points Cantera's own Simpson quadrature is still 16% from converged for a steep-onset cross section (1.3% at 1600).
**Evidence:** tools/genMechTables.C usage() block; validation/mech/validate_cantera_rates.py docstring

**Rule:** Debug-speed tables: `--EN 10:1000:12 --grid 200`. Never publish a rate comparison generated below --grid 1600.


### Table generation is mesh-INDEPENDENT and is deliberately not re-run when only the mesh changes: tools/make_mesh_arm.sh and validation/mesh_sweep.sh both copy constant/plasmaTables and constant/ionTables from the parent case because genMechTables/ionmob 'cost minutes and would be identical'.
**Evidence:** tools/make_mesh_arm.sh header lines 5-7; validation/mesh_sweep.sh lines 40-42

**Rule:** When building a mesh-refinement arm, copy constant/plasmaTables and constant/ionTables rather than regenerating — and keep them in SEPARATE directories, since a re-sweep or Allclean destroys a shared one and the sweep does not put ionTables back.


### `set -e` and `set -u` are both actively harmful in these scripts and are deliberately absent, each with a measured failure: OpenFOAM's etc/bashrc returns non-zero (with `set -e` tools/make_mesh_arm.sh died SILENTLY with an empty logs/ directory, 2026-09-07) and references unset variables (with `set -u` the redirected source line killed mesh_sweep.sh, producing one apparently-launched sweep that never ran and a 0-byte log). needleDBD/Allrun-serial documents the same for RunFunctions line 29 (FOAM_LD_LIBRARY_PATH unbound).
**Evidence:** tools/make_mesh_arm.sh lines 6-9; validation/mesh_sweep.sh lines 48-52; validation/cut200k.sh line 10; needleDBD/Allrun-serial lines 17-20

**Rule:** Do not add `set -eu` to any script that sources OpenFOAM's bashrc or RunFunctions. Guard each step with an explicit `|| die "<step>"` instead.


## Traps
- `./Allwmake` looks like the build but does not build a single test bed. Detector: `grep -c wmake /home/kkourtza/soplasma-scratch/src/applications/Allwmake` returns 4, against 12 utility entries in build-all.sh's BUILD_DIRS.
- build-all.sh's coverage gate scans only `src/applications`, so the SEPARATE tree `/home/kkourtza/soplasma-scratch/applications/test/testChemistryBackends` is never built and never reported. Detector: `ls -la $FOAM_USER_APPBIN/testChemistryBackends $FOAM_USER_LIBBIN/libplasmaChemistry.so` — the binary is Aug 12, the library Sep 5.
- Allwmake/build-all.sh can exit 0 without rebuilding a specific directory (measured for plasmaEnergy, 2026-09-08). Detector: `ls -la $FOAM_USER_LIBBIN/libplasma*.so` and compare against `ls -la src/models/plasmaModels/*/ *.C` mtimes for what you edited; fix with `wclean && wmake libso .` in that one directory.
- A partial `wmake` leaves an ABI-mismatched library and the crash lands somewhere you never touched (SIGFPE in mechanismSourceTerms; std::bad_alloc in soPlasmaFoam; segfault in plasmaTimeControl::adjustDeltaT inside libplasmaTools). Detector: when a crash appears right after a build, run `ls -la $FOAM_USER_LIBBIN/libplasma*.so` BEFORE forming any hypothesis — a control experiment against a stale binary is not a control.
- Forgetting `source etc/bashrc` (or sourcing it BEFORE OpenFOAM's) leaves PETSC_DIR pointing at a nonexistent system path; libplasmaNewtonSolverPETSc then silently fails to rebuild and segfaults later against a changed ABI. Detector: `echo $PETSC_DIR` must be /home/kkourtza/soplasma-scratch/ThirdParty/petsc-3.24.0 and `[ -d $PETSC_DIR/$PETSC_ARCH/lib ]` must be true.
- `wmake src/numerics` exits 0 having left libplasmaNewtonSolverPETSc untouched — it is a separate target. Detector: `ls -la $FOAM_USER_LIBBIN/libplasmaNewtonSolverPETSc.so` mtime.
- A running solver makes build-all.sh exit 1 with NO compile error at all, and the symbols you then find came from an earlier manual wmake. Detector: run `./check-no-running-solvers.sh` first, or check that the output ended with `BUILD-COMPLETE` rather than the REFUSING banner.
- `build-all.sh | grep -E 'error:|BUILD-COMPLETE'` returns GREP's exit status, so a clean build reads as 'Exit code 1'. Detector: never judge a build by a pipeline's status; redirect to a file and grep the file.
- A component log can fail with `Error 1` or `cannot find -l<lib>` and contain zero `error:` lines. Detector: `grep -lE 'Error 1|cannot find -l' /home/kkourtza/soplasma-scratch/buildlogs/*.log`.
- A new header is invisible until `wmakeLnInclude -update`, and the error names the including header instead. Detector: `ls src/<component>/lnInclude/<newHeader>.H` and `grep lnInclude <consumer>/Make/options`.
- `pkill -f soPlasmaFoam` matches — and kills — the calling shell, and `pgrep -x` silently fails on names longer than 15 characters (plasmaChemistry0D). Detector: use `pgrep -f` to look and `pkill -x` to kill, each as its own simple command; verify with a follow-up `pgrep`.
- The first run of any SoPlasma binary after a build takes 24-45 s of cold page-in and looks like a hang or a slow test. Detector: run it twice — the second is 0.03 s.
- testAitken and testFluxScheme ALWAYS exit 0. Putting either in a `&&` gate creates a check that can never fail. Detector: `grep -n 'return 0;' src/applications/utilities/testAitken/testAitken.C` (unconditional) — pair them with an explicit comparison against results.baseline.txt or a human read.
- `Allrun-serial` exiting 0 says only that the SCRIPT finished; OpenFOAM's `runApplication` does not propagate the solver's status. Detector: `grep -c '^Time = ' log.soPlasmaFoam` and `grep -E 'FOAM FATAL|sigFpe|sigSegv' log.soPlasmaFoam`.
- A streamer case whose seed never ran completes 2000 steps with 0 fatal errors and no streamer anywhere. Detector: peak n_e still sitting at the 1e13 background, and `test -s log.initGaussianSeed`.
- Bare `python3` in a script resolves to ~/ct-plasma-env, which has no cantera, no pytest and no pyyaml — a mechc or Cantera step there dies on import. Detector: `which python3` and `~/ct-env/bin/python -c 'import cantera'`.
- SoEEDF validators run from the repo root fail on `../data/...` paths. Detector: testPhysics prints `FATAL: ... (run from the build/ directory)`; the others just throw on the parse.
- The SoEEDF build/ cache is keyed to /home/kkourtza/Projects/BoltzmannSolver, a symlink. Deleting the symlink, or reconfiguring in place from the SoEEDF path, breaks the cache. Detector: `grep CMAKE_HOME_DIRECTORY build/CMakeCache.txt`.
- `set -e` / `set -u` in any script that sources OpenFOAM's bashrc or RunFunctions kills it SILENTLY (empty logs/ directory, 0-byte log, no message). Detector: an 'apparently launched' job with a zero-byte log; guard each step with `|| die` instead.
- `smoke/*` is untracked by git so it never appears in a diff and rots behind every key rename; its recorded state is `native exit=1 / cantera exit=1 / native-stiff exit=1`. Detector: `cat /home/kkourtza/soplasma-scratch/smoke/STATUS` — and sweep smoke/ by path after any dictionary-key rename.
- Editing a RUNNING bash script is unsafe (bash reads it incrementally) — the documented workaround is a second script that stops the specific arm. Detector/pattern: validation/cut200k.sh, which also shows how to kill ONE solver by its /proc/<pid>/cwd rather than by a pattern (a pattern kill cost three self-kills and a set of killed monitor tasks on 2026-09-10).
- testWallFlux had 32 passing checks while the defect shipped, because every gRatio check used Dd = 0, where the right and wrong forms coincide. A unit bed only covers what it VARIES. Detector: before trusting a green bed, check that the parameter the change affects is actually swept in it.
- `plasmaChemistry0D -mechanism` handed the master YAML instead of the compiled .foam writes 0 tables and EXITS 0 (it reads a `speciesCharge` block that only the .foam has). Detector: `grep -c 'wrote .* table' logs/log.ionmob` / count the files in the output directory.

## Open questions
- Exact wall-clock of the SoEEDF validators (testPhysics, validateAir, validateAirSST, validateArgon, testBoltzmann, testMechanism) was NOT measured here: each writes a CSV into its cwd and I was constrained to read-only. Docs assert 'seconds'. Measure once with `cd build && for t in testPhysics validateAir validateAirSST validateArgon testBoltzmann testMechanism; do /usr/bin/time -f "$t %e s" ./$t >/dev/null; done`.
- Exact wall-clock of a from-clean `build-all.sh` is inferred (15-20 min, from memory run-long-cases-in-background.md) rather than measured; the only measured figure is the ~3 min incremental pass of 2026-09-11 17:15-17:18.
- Is `./Allwmake` still maintained, or is it dead weight that should be deleted or made a thin wrapper around build-all.sh? It builds a strict subset and has no coverage gate, so keeping both is itself the trap.
- Why Allwmake silently skipped src/models/plasmaModels/plasmaEnergy on 2026-09-08 was never diagnosed (stale .dep files? a src/models/Allwmake dispatch quirk?). Until it is, the .so-timestamp check after every build is mandatory rather than belt-and-braces.
- `verification/fluxScheme1D/Allrun` has no environment knob for its PE/N loops (unlike fluxScheme2Dnonortho's SHEARS), so there is no fast mode without editing the script. Adding `PES=${PES:-...}` / `NS=${NS:-...}` would be a one-line change.
- `validation/mesh_sweep.sh`'s GRIDS array is hardcoded; NSTEPS/NRELAX/DT are overridable but the arm list is not. A `GRIDS` env override would remove the need for the cut200k.sh workaround entirely.
- testPhysics has no subset flag — its 17 test groups are called unconditionally from `int main()`. If a per-group `--only <name>` is wanted, that is a small, self-contained addition.
- `verification/testSnesJFNK` and `testSnesJFNK2Field` ship as cases with processor0-3 directories but no Allrun script was found; the exact serial and 4-rank invocations (and how PETSC_OPTIONS is set for them) are documented only in prose in docs/design/newton-outer-solver-design.md.
- The three `smoke/` beds currently record exit=1. Whether they are broken by another key rename or by something else has not been re-checked since 2026-09-02.
