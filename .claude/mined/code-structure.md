# code-structure

## Summary
The two-layer case architecture (G1/G2) is the product vision: a case is TWO hand-written files -- `configuration/config` (every number, referenced as `$key`) and `configuration/boundaries` (one block per patch saying WHAT that surface IS) -- and every OpenFOAM dictionary and every field in `0/` is GENERATED from them plus the mesh topology plus three shipped libraries under `$SoPLASMA_ETC`. The generators are two OpenFOAM utilities, `plasmaSetupBoundaries` (writes `0/<region>/ePotential` and `0/<region>/surfCharge`) and `plasmaCreateSpeciesFields` (writes every `n_<species>` and `nEps_e`), both reading the SAME `Foam::boundaryRoleLibrary` so their order does not matter -- an ordering bug that existed until commit `db77095` (2026-09-04). The role library `/home/kkourtza/soplasma-scratch/etc/boundaryRoles` defines 9 kinds, each reducing to ONE word for species transport (`mechanical|open|conductor|chargingSurface`). Three things are DERIVED and a case may not declare them: mechanical patches (from the mesh), region interfaces `<a>_to_<b>` (from topology, with surface charge owned by the GAS side only), and epsilonR/gammaSEE behind a named material. The canonical minimal case is `tutorials/plasma/soPlasmaFoam/needleDBD` (two hand-written config files, no `etc/`, no `0.orig/`); the positiveStreamer beds are only HALF converted and are the live trap. The architecture's sharp edge is the `$key` indirection: a config variable nothing references is silently inert, and the naive detector for it is broken by `$$` shell expansion.

## Facts

### The complete boundary-role vocabulary is 9 kinds in /home/kkourtza/soplasma-scratch/etc/boundaryRoles, each with exactly one `plasmaTransport { surface ... }` word. drivenElectrode=conductor (requires `waveform`); groundedElectrode=conductor (requires nothing); floatingElectrode=conductor (optional `initialCharge 0`); thinDielectricSurface=chargingSurface (requires `material`); thinDielectricOnElectrode=chargingSurface (requires `material thickness backingPotential`); openBoundary=open (requires nothing); ballastedElectrode=conductor (requires `circuit`); currentDrivenElectrode=conductor (requires `circuit`); insulatingWall=conductor (requires nothing). Every conductor/dielectric kind also takes optional `material none; gammaSEE -1; electronReflection -1` (floatingElectrode swaps `material` for `initialCharge`).
**Evidence:** etc/boundaryRoles lines 119,153,182,225,264,305,336,377,415 (kind headers); `awk '/^[a-zA-Z][a-zA-Z0-9_]*$/{k=$1} /plasmaTransport *\{ *surface/{print k": "$0}' etc/boundaryRoles` prints all nine

**Rule:** Never invent a `kind`. Run `plasmaSetupBoundaries -listKinds` and pick one of the nine; if a surface is not expressible, add a SUB-DICTIONARY to an existing kind (G3), not a new name.


### `insulatingWall` declares `surface conductor`, NOT chargingSurface. It absorbs the wall flux but accumulates no sigma, because its electrostatics is `zeroGradient` -- the same equation with sigma forced to zero. A wall that is meant to trap charge and is declared `insulatingWall` will silently lose every charge that reaches it; the correct kind is `thinDielectricSurface`.
**Evidence:** etc/boundaryRoles:415-445, comment block above `electrostatics { type zeroGradient; }` in the insulatingWall entry

**Rule:** If charge must stay on a non-electrode surface, use `thinDielectricSurface` (requires a `material`). `insulatingWall` is only for a wall you have decided does NOT charge.

**Cost:** A converged, plausible, wrong DBD: the surface charge that makes a barrier discharge self-limiting is accumulated and discarded with nothing in the log.


### `plasmaSetupBoundaries` (binary at $FOAM_USER_APPBIN/plasmaSetupBoundaries) has exactly ONE CLI option: `-listKinds` (an addBoolOption). `plasmaCreateSpeciesFields` has exactly one: `-dict <file>` (alternative plasmaSpecies dictionary), plus noBanner/noJobInfo/noFunctionObjects. Neither takes a `-region` flag: plasmaSetupBoundaries loops every region in constant/regionProperties itself, and plasmaCreateSpeciesFields finds the single `gas` region itself.
**Evidence:** src/applications/utilities/foamPlasmaSetupBoundaries/plasmaSetupBoundaries.C:111-126 (addNote, addBoolOption "listKinds", args.found("listKinds")); src/applications/utilities/foamPlasmaCreateSpeciesFields/plasmaCreateSpeciesFields.C:666-682 (addOption "dict"); Make/files EXE lines

**Rule:** Do not pass `-region` or guess flags. `plasmaSetupBoundaries -listKinds` is the only discovery command; run both generators bare from the case directory.


### The two generators are ORDER-INDEPENDENT, verified 2026-09-04: all 17 generated files are byte-identical whichever runs first. It briefly was NOT: plasmaCreateSpeciesFields used to read 0/gas/ePotential and map the electrostatics BC type back to a semantic class, so running it first classified EVERY patch as `open` and emitted no wall flux anywhere. Both now derive from configuration/boundaries through the shared Foam::boundaryRoleLibrary.
**Evidence:** commit db77095 (2026-09-04) "Species BCs derive from the DESCRIPTION, not from the electrostatics BC type"; memory boundary-role-library-layer1.md; plasmaCreateSpeciesFields.C:95-112 comment; needleDBD/Allrun-serial step 6 comment

**Rule:** Run `plasmaSetupBoundaries` and `plasmaCreateSpeciesFields` in either order, but ALWAYS after `rm -rf 0 && mkdir -p 0/<each region>`. Never reintroduce a dependency where one module reads another module's generated output.


### `Foam::boundaryRoleLibrary` and `Foam::materialLibrary` are the ONE-OWNER resolvers. The role library is anchored on $SoPLASMA_ETC and is NOT searched for, so which library answered cannot depend on cwd; the case declaration is read from `runTime.globalPath()/configuration/boundaries` -- globalPath(), not path(), because in parallel Time::path() returns <case>/processorN and every rank died looking for its own copy (measured 2026-09-07, blocked parallel execution for every G2 case).
**Evidence:** src/models/electromagnetics/boundaryRoleLibrary/boundaryRoleLibrary.C:33-54 (getEnv("SoPLASMA_ETC")) and :144-169 (globalPath comment); boundaryRoleLibrary.H:16-26 (ONE OWNER rationale)

**Rule:** Before running any generator, assert `[ -f "$SoPLASMA_ETC/boundaryRoles" ]`. If SoPLASMA_ETC is unset the utility aborts with "SoPLASMA_ETC is not set, so the boundary-role library cannot be..."; if it points at the wrong place it aborts with `Cannot read the boundary-role library: "<path>/boundaryRoles"`.


### `export A=x B=$A` in ONE bash statement expands $A BEFORE assigning it. A one-liner setting SoPLASMA and SoPLASMA_ETC together produced SoPLASMA_ETC=/etc and plasmaSetupBoundaries died on `Cannot read the boundary-role library: "/etc/boundaryRoles"`. Measured 2026-09-07.
**Evidence:** tools/make_mesh_arm.sh lines 19-24, comment "SEPARATE lines, deliberately"

**Rule:** Always set SoPLASMA, SoPLASMA_ETC and SoPLASMA_SRC on SEPARATE lines, then guard with `[ -f "$SoPLASMA_ETC/boundaryRoles" ] || exit 1`.


### `configuration/config` is a CASE CONVENTION, not an OpenFOAM feature: dictionaries `#include "../configuration/config"` and reference `$key`. A key nothing references does nothing, silently. MEASURED NOW on the shipped needleDBD case: of its 14 config keys, exactly one is dangling -- `simulationType`, because system/plasmaSimulationControls hardcodes `simulationType transient;` instead of `$simulationType`. Editing it to `pseudoSteady` would change nothing.
**Evidence:** /home/kkourtza/soplasma-scratch/tutorials/plasma/soPlasmaFoam/needleDBD/system/plasmaSimulationControls:29 `simulationType transient;` vs configuration/config:317 `simulationType transient;`; detector run 2026-09-11 reports exactly `DANGLING: simulationType`

**Rule:** After editing any configuration/config value, prove something reads it with the detector below, or verify the effect in the OUTPUT (a field range, a start-up line). Never assume the edit landed.

**Cost:** A whole test arm run at the wrong condition: setting `appliedVoltage 3019` for a ~10 Td arm ran clean, exit 0, 46 steps -- at the original 62.1 Td (memory config-variable-may-be-dangling.md, 2026-09-01).


### THE $$ TRAP in the dangling detector. `grep -rq "\$$k"` inside DOUBLE quotes expands `$$` to the shell PID, so it searches for `12345k` and reports nearly every variable as dangling. MEASURED NOW: the naive form reports 14 of needleDBD's 14 keys dangling (100% false positive); the correct form reports 1. The correct spelling is `grep -rqF '$'"$k"`.
**Evidence:** Both forms run in /home/kkourtza/soplasma-scratch/tutorials/plasma/soPlasmaFoam/needleDBD on 2026-09-11: naive -> 14 DANGLING lines, correct -> 1; memory config-variable-may-be-dangling.md "I 'found' 47-48 per case before noticing"

**Rule:** Use `grep -rqF '$'"$k"` (single-quoted dollar concatenated with the double-quoted variable), never `"\$$k"` inside double quotes. A detector that flags EVERY variable is the signature of this bug -- not of a broken case.


### THE INVERSE DEFECT, live in needleDBD today: `configuration/config` documents a variable that does not exist. Its comment says the library list is held there "as ONE WHOLE VALUE so the list can be swapped from here" and instructs "add libpetscFoam.so to this list" -- but no such key is declared in config, and system/controlDict:18 hardcodes `libs ( "libROUNDSchemes.so" );`. Following the documented PETSc switch edits nothing.
**Evidence:** `grep -n "solverLibs|libs" configuration/config system/controlDict` in the needleDBD case returns only system/controlDict:18; configuration/config has the prose block but no key

**Rule:** To switch needleDBD to PETSc you must edit BOTH `configuration/config` (`ePotentialSolver petsc;`) AND `system/controlDict` (add "libpetscFoam.so" to `libs`). Treat the config comment as stale. tools/checkConfigReference.py --case is the tool for this class.


### NEVER hand-edit: everything under `0/` (every field is generated by the two generators), `system/<region>/{fvSchemes,fvSolution}` and `constant/<region>/{plasmaSpeciesProperties,plasmaTransportProperties,photoionizationProperties}` (ten-line `#include` stubs written by tools/plasmaSetupRegions.sh, and OVERWRITTEN with empty dummies by splitMeshRegions), and `system/fvSolution` in any bed that ships `system/fvSolution-foam`/`-petsc` (Allrun copies one over it, so edits vanish). The value belongs, respectively, in configuration/boundaries, in the case-wide system/fvSchemes+fvSolution, in constant/<dict>, and -- for needleDBD -- directly in the single system/fvSolution.
**Evidence:** needleDBD/Allclean echoes the exact hand-authored list; tools/plasmaSetupRegions.sh stub headers say "GENERATED ... DO NOT EDIT"; tutorials/plasma/soPlasmaFoam/positiveStreamer/positiveStreamer_fixedMesh/Allrun-serial:42-46 `cp system/fvSolution-foam system/fvSolution`

**Rule:** Before editing any dictionary, `head -8` it: a GENERATED banner or an `#include "../<file>"` body means edit the target, not the file. Before editing system/fvSolution, `ls system/fvSolution-*` -- if a `-foam`/`-petsc` pair exists, that file is a build artefact.


### needleDBD deliberately has NO fvSolution-foam/-petsc pair: one `ePotential` block carries BOTH backends and `$ePotentialSolver` picks one. This works because MEASURED 2026-09-02 an OpenFOAM solver dictionary SILENTLY IGNORES entries it does not use -- a PCG block carrying `agglomerator`, `nCellsInCoarsestLevel` and a `petsc { options { ... } }` sub-dictionary runs clean.
**Evidence:** tutorials/plasma/soPlasmaFoam/needleDBD/system/fvSolution:13-30 header comment and the single `ePotential` block carrying agglomerator + petsc{options{}} together

**Rule:** When converting a bed off the -foam/-petsc pair, merge both blocks into one and add `$ePotentialSolver` to configuration/config; do not add a second dictionary.


### Under monolithic region coupling (`useImplicit true`) GAMG MUST use `agglomerator assembledFaceAreaPair`. The default `faceAreaPair` does `refCast<const fvMesh>` on an `lduPrimitiveMeshAssembly` and aborts with "Attempt to cast type lduPrimitiveMeshAssembly to type fvMesh" -- which reads like an internal error and is a solver-configuration error. MEASURED 2026-09-01 on the two-region plate2D case (20000 cells, per timestep): smoothSolver/GaussSeidel 2000 iters (cap), resid 3e-4, 7.3% WRONG; PCG/DIC 114; GAMG default ABORT; GAMG+assembledFaceAreaPair 11.
**Evidence:** tutorials/plasma/soPlasmaFoam/needleDBD/system/fvSolution, ePotential block comment (the four-row measured table)

**Rule:** Any generated/edited fvSolution for a multi-region case must carry `agglomerator assembledFaceAreaPair;` under GAMG. Reaching maxIter is NOT an error in OpenFOAM -- grep the log for solves at the cap.


### THE ONE HAND-WRITTEN FILE NOTHING GENERATES is `constant/regionProperties`, and its `kind` keywords are a closed set of three: `gas` (exactly one, epsilonR defaults 1.0), `dielectric` (any number, epsilonR REQUIRED with no default), `farField` (any number, epsilonR 1.0 by definition). Using a wrong keyword -- e.g. `plexiglass (plexiglass)` instead of `dielectric (plexiglass)` -- used to be SILENT: the dielectric list came back empty, the model was derived as singleRegionPoisson, and the potential was solved in the gas alone. Guarded since 2026-09-02 in electromagneticsModel::validateRegionKinds. `foamListRegions` is no check -- it flattens every list regardless of kind.
**Evidence:** tutorials/plasma/soPlasmaFoam/needleDBD/constant/regionProperties, the FORMAT and "USING A WRONG KIND IS FATAL" comment blocks

**Rule:** Write regionProperties by hand as `regions ( gas (<name>) dielectric (<name> ...) );` and verify the names agree in FIVE places: the gmsh Physical Volume, the cellZone names gmshToFoam prints, the parentheses here, the constant/<region>/ directory names, and any etc/changeDictionary.<region> filename.


### The material library is /home/kkourtza/soplasma-scratch/etc/materials/dielectrics with 12 entries -- vacuum(1.0), air(1.0), ptfe(2.1), pmma(3.0), polyimide(3.4), fusedSilica(3.8), borosilicateGlass(4.6), sodaLimeGlass(7.0), mica(5.4), alumina96(9.0), alumina99(9.8), magnesia(9.7) -- each carrying epsilonR, gammaSEE and a mandatory `reference` string (G2: an entry without a reference is a number somebody remembered). Two sibling libraries exist: etc/materials/workFunctions and etc/materials/ionisationEnergies. Precedence per G1: explicit number in the case WINS (reported as an override) -> material -> region-kind default; a misspelt material is FATAL and lists the twelve.
**Evidence:** etc/materials/dielectrics (entry names at lines 80,87,99,112,126,141,155,167,180,196,213,226); materialLibrary.H:14-33 (PRECEDENCE and ONE OWNER); `ls etc/materials/` -> dielectrics, ionisationEnergies, workFunctions

**Rule:** Never type epsilonR or gammaSEE into a case. Write `material pmma;` in configuration/boundaries (thin-dielectric kinds) or in constant/<region>/electricalProperties. To add a material, add BOTH numbers AND a `reference` string.


### gammaSEE is 0.001 for 11 of the 12 dielectrics (magnesia=0.3 is the deliberate exception), matching the `defaultSEEC` default changed to 0.001 on 2026-09-02 so gamma cannot depend on whether a case named a material. The real material dependence now comes through `emission { ionInducedSEE { yield hagstrum; } }`, verified 2026-09-04: on copper N2p=0.05063424 and O2p=0.00569216, an 8.9x spread between two ions on the SAME material.
**Evidence:** etc/materials/dielectrics header "GAMMA ... IS THE LEAST TRANSFERABLE NUMBER IN THIS FILE"; memory material-library.md (commit afb7913 and 7464bd9)

**Rule:** Testing that a material reaches a BC requires a material whose gamma DIFFERS -- use `magnesia` (0.3), never `pmma` (0.001 = the default, which proves nothing).


### EVERY shipped kind ships with NO emission, deliberately: emission decides whether a discharge is self-sustaining rather than a single avalanche, so a default would silently change the physics of every existing case. Turning it on is either `gammaSEE <n>;` on the patch (the numeric route -- it WINS over `material`, because a validation case reproducing a published gamma must use THAT number) or an `emission { ... }` block on the patch (which overrides the kind's). -1 means NOT DECLARED and is distinct from a declared 0 (a perfectly absorbing, non-emitting wall).
**Evidence:** etc/boundaryRoles header "EVERY KIND SHIPS WITH NO EMISSION, DELIBERATELY"; plasmaCreateSpeciesFields.C:600-640 (gammaSEE/eReflection initialised to -1, patch overrides role) and :495-525 (enableSEE/defaultSEEC emission)

**Rule:** To add SEE, put `gammaSEE 0.06;` on the patch block in configuration/boundaries -- never reach into a changeDictionary. `electronReflection` is written on BOTH the electron and the energy condition automatically (the Hagelaar closure uses it in the energy weight); do not write it on one.


### `speciesSurface` is the ONLY sanctioned override of what species do at a surface, and exists solely for benchmark reproduction (the Bagheri / Pasolari & Kourtzanidis arXiv:2607.05137 sec 6.1 positive-streamer problem prescribes the electrode potential while imposing homogeneous Neumann on the densities -- a non-absorbing electrode, which does not exist). It is NOT a new kind, it is validated against `open|conductor|chargingSurface`, it is FATAL on a mechanical patch, and it is ANNOUNCED on every run because "the electrodes do not absorb" silently changes an answer.
**Evidence:** plasmaCreateSpeciesFields.C:230-300 (the override block, its two FatalErrorInFunction guards and the Info announcement); etc/boundaryRoles "OVERRIDING `surface` ON ONE PATCH" section; used at tutorials/.../positiveStreamer_fixedMesh/configuration/boundaries:81,88

**Rule:** Grep any case for `speciesSurface` before trusting its wall physics: `grep -rn speciesSurface */configuration/boundaries`. Its presence means an electrode absorbs nothing -- no drift flux, no thermal flux, no secondary emission.


### THE HALF-CONVERTED CASE IS THE LIVE TRAP. positiveStreamer_fixedMesh has a `configuration/boundaries` (added 2026-09-10 only because plasmaCreateSpeciesFields refuses to start without layer 1) but its Allrun-serial NEVER runs plasmaSetupBoundaries: it copies hand-written 0.orig/{ePotential,surfCharge}, runs plasmaCreateSpeciesFields, then `changeDictionary`, which WINS on every field it names. The two paths DISAGREE on physics, measured 2026-09-10: at the electrodes layer 1 gives electron/ionDDWallFluxMixed (ABSORBING) while changeDictionary gives zeroGradient (NON-absorbing); at `far`, layer 1 gives zeroGradient while changeDictionary gives inletOutlet with inletValue 1e13.
**Evidence:** tutorials/plasma/soPlasmaFoam/positiveStreamer/positiveStreamer_fixedMesh/Allrun-serial:59-67 (cp 0.orig, plasmaCreateSpeciesFields, cp etc/changeDictionary, changeDictionary -- no plasmaSetupBoundaries) and configuration/boundaries:14-46 (the explicit disagreement table)

**Rule:** Before editing any case's configuration/boundaries, check `grep -c plasmaSetupBoundaries <case>/Allrun*` and `ls <case>/0.orig <case>/etc 2>/dev/null`. If 0.orig or etc/changeDictionary exist, layer 1 is NOT the only source of truth and your edit may be overwritten by a later changeDictionary pass.

**Cost:** Editing configuration/boundaries there changes nEps_e and the excited/neutral species but NOT ePotential, n_e, n_N2p, n_O2p, n_Om, n_O2m -- a partial, invisible change.


### Region interfaces and mechanical patches are DERIVED by plasmaSetupBoundaries and a case is REFUSED if it declares them. A gas/dielectric pair emits `coupledElectricPotential` with `useImplicit true` and the surface charge owned by the GAS side only (`surfCharge surfCharge; surfChargeNbr none;` on gas, inverted on the dielectric); a gas/farField pair emits plain `zeroGradient` with no charging, because a farField region is fictitious air and not a surface. Exactly one owner, or the charge is double-counted.
**Evidence:** plasmaSetupBoundaries.C:283-315 (the interface branch writing both halves) and :11-30 (the THREE THINGS ARE DERIVED header)

**Rule:** Do not write `gas_to_dielectric` blocks in configuration/boundaries. If a patch you did not expect is reported undeclared, check whether msh2Dto3D left an interface edge as a boundary patch -- that is a mesh bug, not a declaration gap.


### tools/plasmaSetupRegions.sh MUST run AFTER splitMeshRegions and is idempotent. OpenFOAM reads fvSchemes/fvSolution PER REGION from system/<region>/ with NO fallback to the case-wide copy (measured 2026-09-02: absent system/gas/fvSchemes -> "cannot find file .../system/gas/fvSchemes"), and splitMeshRegions writes an EMPTY dummy, which is WORSE than missing because checkMesh is happy and the first scheme lookup then dies with "Entry 'grad(...)' not found". The script writes ten-line #include stubs instead of copies, so there is nothing to drift and no include-depth surgery.
**Evidence:** tools/plasmaSetupRegions.sh header "WHY THIS EXISTS" and "MUST RUN AFTER THE SPLIT ... measured, the stub's md5 changes"

**Rule:** Sequence is fixed: `splitMeshRegions -cellZones -overwrite` THEN `tools/plasmaSetupRegions.sh` (run from the case directory; it errors out if system/controlDict is absent). Re-run it after every re-split.


### plasmaSetupRegions.sh writes `constant/<region>/electricalProperties` per region KIND and NEVER overwrites an existing one. For `gas` and `farField` it writes a complete `epsilonR 1.0;`. For a `dielectric` it writes the entry COMMENTED OUT and lists the file under ACTION REQUIRED -- deliberately, because epsilonR=1 for a barrier is not a placeholder but VACUUM: the case would run, converge and give a plausible WRONG answer with the barrier electrically invisible, and nothing in the log would say so. The solver then stops with an error naming that exact path.
**Evidence:** tools/plasmaSetupRegions.sh, writeElectricalProperties() -- the three branches and the closing "ACTION REQUIRED" banner

**Rule:** After plasmaSetupRegions.sh, grep its log for `NEEDS A VALUE` and uncomment the epsilonR in every file it names (the file itself lists common barrier permittivities: quartz 3.8, borosilicate 4.6, soda-lime 7.0, alumina96 9.0, PTFE 2.1, PMMA 3.0, Kapton 3.4, mica 6.0).


### tools/make_mesh_arm.sh <parent_case> <target_case> <NX> <BUMP> is the canonical 'copy a case and regenerate' recipe: it rebuilds the mesh and `0/` ONLY, reusing the parent's Boltzmann and ion tables because they are mesh-INDEPENDENT. It deliberately does NOT use `set -e`: OpenFOAM's etc/bashrc returns non-zero, and with -e the script dies SILENTLY with no log and no message (measured 2026-09-07 -- the first test produced an empty logs/ directory and no output at all). Each step is guarded by die() instead.
**Evidence:** tools/make_mesh_arm.sh lines 1-12 (header) and the per-step `|| die <name>` chain

**Rule:** In any SoPlasma shell script that sources OpenFOAM: use `die()` per step, never `set -e`. Guard the environment with `command -v gmshToFoam >/dev/null || exit 1` and `[ -f "$SoPLASMA_ETC/boundaryRoles" ] || exit 1`.


### THE SINGLE MOST IMPORTANT CHECK IN A MULTI-REGION BUILD is that splitMeshRegions actually built an interface -- if it did not, the two regions never couple and NOTHING ELSE ERRORS. The data row is the FOURTH line after the "Sizes of interfaces between regions" header (blank, column names, dashes, then the row), so `grep -A3` prints an empty-looking table and makes a correct split look like the exact failure the check exists to detect. Use `-A6`.
**Evidence:** tutorials/plasma/soPlasmaFoam/needleDBD/Allrun-serial step 3 guard; TUTORIAL.md:408-410 "NOTE -A6, not -A3"

**Rule:** After splitMeshRegions run `grep -A6 "Sizes of interfaces" log.splitMeshRegions | grep -E "^[0-9]"` and require at least one row.


### msh2Dto3D.py exists because extruding the CAD instead of the MESH generates two measured traps: under SetFactory("OpenCASCADE") extruding two surfaces that SHARE a curve produces coincident-but-separate lateral faces, so the volumes are non-conformal and splitMeshRegions cannot build the interface; and `Coherence;` renumbers entities and invalidated the Extrude indices the physical groups were built from -- face counts collapsed from 103/160/338/352 to 3/1/7/184. It also AUTOMATICALLY detects and excludes region-interface edges: in the shipped needle mesh the group `air_dielectric` mixes 160 interface edges with 178 genuine wall edges, so excluding the whole group would also be wrong.
**Evidence:** tools/msh2Dto3D.py docstring lines 15-55; needleDBD/configuration/boundaries comment "178 faces left over after splitMeshRegions took the 160 that became the interface"

**Rule:** Always go 2-D .msh -> msh2Dto3D.py -> gmshToFoam -> msh2Dto3D.py --fix-boundary. Never re-author the .geo to extrude the CAD, and never name interface edges as a boundary patch.


### Allrun-serial has just been made location-independent (uncommitted change in the working tree): `SoPLASMA_TOOLS=../../../../tools` is replaced by `SoPLASMA_TOOLS=${SoPLASMA:+$SoPLASMA/tools}` with the bare relative path as fallback. The bare path resolves only at the file's original depth, so a COPY of the case (a benchmark arm) died at msh2Dto3D with "../../../../tools/msh2Dto3D.py: No such file or directory".
**Evidence:** `cd /home/kkourtza/soplasma-scratch && git diff tutorials/plasma/soPlasmaFoam/needleDBD/Allrun-serial` -- the only modified hunk, lines 29-37

**Rule:** When templating a new case from needleDBD, derive every tool path from `$SoPLASMA`, never from a relative `../../../..`. Set `export SoPLASMA=$HOME/soplasma-scratch` before running a copied case.


### needleDBD's Allrun-serial has NO `set -u`, deliberately: OpenFOAM's own bin/tools/RunFunctions references FOAM_LD_LIBRARY_PATH unbound at line 29, so `set -u` makes sourcing it fatal with "RunFunctions: line 29: FOAM_LD_LIBRARY_PATH: unbound variable". The variables the script actually depends on are guarded individually.
**Evidence:** tutorials/plasma/soPlasmaFoam/needleDBD/Allrun-serial lines 19-23

**Rule:** Never add `set -u` to a script that sources $WM_PROJECT_DIR/bin/tools/RunFunctions. Guard with `[ -z "${VAR:-}" ] && exit 1` per variable.


### needleDBD's ion transport tables come from SoEEDF's tools/ionmob.py, NOT from the Boltzmann sweep (which solves the ELECTRON energy distribution; an ion drifting through a neutral gas is a different problem). `--mechanism` wants the COMPILED .foam, not the master YAML -- it reads a `speciesCharge` block, and passing the YAML writes 0 tables and EXITS 0. constant/ionTables and constant/plasmaTables are deliberately SEPARATE directories, because when they shared one, any re-sweep or Allclean destroyed the ion tables while the sweep believed its own set complete.
**Evidence:** tutorials/plasma/soPlasmaFoam/needleDBD/Allrun-serial step 5 comments; Allclean `rm -rf constant/plasmaTables constant/ionTables` comment

**Rule:** Pass `--mechanism constant/<name>.foam` to ionmob.py and then assert `ls constant/ionTables | wc -l` is non-zero -- exit 0 does not mean tables were written. Read the NOTES block: `grep -A20 NOTES logs/log.ionmob` -- those are modelling DECISIONS (substituted gas pairs, stitched datasets, Einstein-derived diffusion), not diagnostics.


### If a case DOES keep an etc/changeDictionary.<region> override, changeDictionary STRIPS the `#include` from every field it rewrites, so any `$var` it copied through verbatim is then undefined and the solver dies reading the field. The fix must re-insert the include AFTER the FoamFile block, or the reader fails with "problem while reading header". needleDBD's Allrun-serial carries an inline python3 heredoc that does exactly this after EVERY changeDictionary pass.
**Evidence:** tutorials/plasma/soPlasmaFoam/needleDBD/Allrun-serial step 7, the FIXPY heredoc and its comment; memory derived-species-bcs.md "the include fixup must run after EVERY pass"

**Rule:** changeDictionary is now OPTIONAL and normally absent (needleDBD ships none). If you add one, copy the FIXPY heredoc verbatim and run it after every pass. Verify with `foamDictionary 0/<region>/ePotential -entry boundaryField/<patch>/uniformValue`.


### The ambient/background density used for `inletOutlet` inflow at an openBoundary, and for every generated internalField, is DERIVED from each species' own `minNumberDensity` (G1) rather than restated. openBoundary means "the gas continues outside at ambient", not zeroGradient -- changed in commit 438b3a7 on 2026-09-10 -- because in a positive streamer the outer boundary IS an inflow boundary for electrons (they drift against E, hence inward), and under zeroGradient a net-attaching bulk drains the background out of the domain.
**Evidence:** etc/boundaryRoles openBoundary block; plasmaCreateSpeciesFields.C:341-370 (the inletOutlet branch) and :1035-1050 ("THE AMBIENT (BACKGROUND) DENSITY IS THE SPECIES' OWN DECLARED FLOOR")

**Rule:** Never restate a background density in a boundary entry. Set `minNumberDensity` once per species in plasmaSpeciesProperties; the generator uses it for both the initial field and the inflow value.


### plasmaSetupBoundaries writes `0/<region>/surfCharge` on the GAS region ONLY, as `calculated` everywhere plus the mesh's own type on mechanical patches. The VALUES are written by plasmaTransport::updateSurfaceCharge from the species wall fluxes; a dielectric region's surfCharge is READ_IF_PRESENT and its interface names the gas side's field as `surfChargeNbr`. Accumulation happens only where a species patch field is a plasmaWallBC AND has `enableSurfaceCharging true` (plasmaTransport::updateSurfaceCharge does a dynamic_cast and skips everything else).
**Evidence:** plasmaSetupBoundaries.C:465-520 (the surfCharge block, `if (isGas)`); plasmaSetupRegions.sh farField note "Surface charge accumulates ONLY where ... dynamic_cast"

**Rule:** Never hand-write a surfCharge field. To check who owns sigma: `grep -l enableSurfaceCharging.*true 0/gas/*` -- it should be exactly the chargingSurface patches.


### The template placeholders in etc/boundaryRoles are quoted-and-bracketed `"<param>"`, NOT `$param`, because OpenFOAM's own dictionary reader resolves `$name` at PARSE time and aborts on an unknown one -- so a `$` placeholder never reaches the generator. `"<material:epsilonR>"` resolves through materialLibrary; `"<internalField>"` becomes `uniform 0`; a plain `"<name>"` takes the patch's own entry verbatim (so a whole Function1 like `table ((0 0) (1e-7 8e3))` survives), falling back to the kind's `optional` defaults, else FATAL.
**Evidence:** plasmaSetupBoundaries.C:388-440 (the substitution branch and its comment); memory boundary-role-library-layer1.md trap 1

**Rule:** When adding a kind to etc/boundaryRoles, write `"<param>"` never `$param`. Also: an ITstream returned by lookup() must be INDEXED as a tokenList, not driven with `is >> t` -- the stream form returns EMPTY and generated `type ;` for every patch.


### needleDBD is intentionally sized to finish: endTime 3e-09 (not 2e-08), because MEASURED 2026-09-02 the case runs at ~2.8 s/step, making 2e-08 roughly a THREE-HOUR run. 3 ns crosses discharge initiation (with ~15x tip enhancement, E/N reaches the ionisation range at V ~ 167 V, i.e. t ~ 2.1 ns on the 8 kV / 100 ns ramp) and finishes in well under an hour. Raise to 2e-08 for the barrier-charging physics that makes a DBD self-limiting.
**Evidence:** tutorials/plasma/soPlasmaFoam/needleDBD/system/controlDict:33 `endTime 3e-09;` and configuration/config TIME CONTROL block

**Rule:** For a smoke test of the generator chain, do NOT run the solver at all -- steps 1-7 of Allrun-serial (through both generators) are seconds. Only step 8 costs minutes-to-hours.


### needleDBD_floatingStrip predates the role library and still post-edits GENERATED fields with an inline python regex (rewriting `air_dielectric` in 0/gas/ePotential to floatingElectrodePotential and flipping enableSurfaceCharging to false across 0/gas/*). Layer 1 can now express this in two lines (`kind floatingElectrode; initialCharge 0;`), so the python surgery is obsolete machinery that will break the moment the generator's output format changes.
**Evidence:** tutorials/plasma/soPlasmaFoam/needleDBD_floatingStrip/Allrun-derive lines 30-60 (the `PY` heredoc)

**Rule:** Treat any python/sed post-edit of a file under 0/ as a conversion debt item: it should become a kind + parameters in configuration/boundaries.


### Two design requirements from the agreed architecture are still NOT implemented: an `advanced { }` passthrough in layer 1 (a generator with no bypass is a CEILING -- the first user needing an unmodelled option hand-edits a generated file and loses it), and REFUSAL TO CLOBBER (the generator must refuse to overwrite a layer-2 file whose checksum changed, the rule plasmaSetupRegions.sh already follows). Today Allrun-serial does `rm -rf 0` before generating, so any hand edit under 0/ is destroyed with no warning.
**Evidence:** memory two-layer-case-architecture.md "TWO THINGS THE DESIGN MUST HAVE"; needleDBD/Allrun-serial step 6 `rm -rf 0 && mkdir -p 0/gas 0/dielectric logs`; plasmaSetupBoundaries.C writes with OFstream unconditionally

**Rule:** Assume anything under 0/ is destroyed on the next Allrun. If you must experiment, use etc/changeDictionary.<region> (which runs AFTER both generators and wins), not a direct edit.


### plasmaSetupBoundaries and plasmaCreateSpeciesFields now BOTH handle single-region cases (no constant/regionProperties). plasmaSetupBoundaries' read was MUST_READ until 2026-09-05, so a plain gas gap between two boundary electrodes aborted with "cannot find file constant/regionProperties" and could not use the semantic layer at all -- while plasmaCreateSpeciesFields already fell back correctly. The two halves of one generator disagreed about the same guard.
**Evidence:** plasmaSetupBoundaries.C:160-218 (the isFile guard and its comment, "Measured 2026-09-05, the first time either was run on one"); the fallback treats the default region as kind `gas`

**Rule:** A single-region case needs NO regionProperties and no constant/<region>/ move. If a generator demands one, the guard has regressed -- check both utilities agree.


## Traps
- THE $$ TRAP: `grep -rq "\$$k"` inside double quotes expands `$$` to the shell PID, so the dangling-variable detector searches for `12345k` and reports nearly every variable as dangling. DETECTOR: if the check flags EVERY key (14/14 on needleDBD), the detector is broken, not the case. Fix: `grep -rqF '$'"$k"`.
- A dangling configuration/config variable is silently inert. DETECTOR: the for-loop above; on needleDBD it correctly reports `simulationType`. Cost precedent: an `appliedVoltage 3019` low-field arm ran clean at the ORIGINAL 62.1 Td.
- THE INVERSE: a config COMMENT documenting a variable that does not exist. needleDBD's config explains how to swap the solver library list "from here" but no key exists and system/controlDict:18 hardcodes `libs ( "libROUNDSchemes.so" );`. DETECTOR: `tools/checkConfigReference.py --case <dir>`, or grep the key name in config itself before believing the prose.
- THE HALF-CONVERTED CASE: a case with configuration/boundaries whose Allrun never runs plasmaSetupBoundaries and still applies etc/changeDictionary afterwards -- changeDictionary WINS on every field it names. DETECTOR: `grep -c plasmaSetupBoundaries <case>/Allrun*` and `ls <case>/0.orig <case>/etc`. Live example: positiveStreamer_fixedMesh, where the two paths disagree on whether the electrodes absorb.
- splitMeshRegions OVERWRITES system/<region>/ with an EMPTY dummy and says "Writing dummy". Empty is worse than missing: checkMesh is happy and the first scheme lookup dies with "Entry 'grad(...)' not found in dictionary system/<region>/fvSchemes/gradSchemes". DETECTOR: `wc -l system/*/fvSchemes` -- a real stub is ~25 lines and contains `#include "../fvSchemes"`; a dummy is near-empty.
- A dielectric region with epsilonR left at 1.0 runs, converges, and gives a plausible WRONG answer with the barrier electrically invisible -- and nothing in the log says so. DETECTOR: `grep -rn '^ *epsilonR' constant/*/electricalProperties` -- a `dielectric`-kind region with the line still commented is caught by the solver, but one edited to 1.0 is not.
- A wrong region KIND keyword in constant/regionProperties (e.g. `plexiglass (plexiglass)` instead of `dielectric (plexiglass)`) used to silently derive singleRegionPoisson and solve the potential in the gas alone. Guarded since 2026-09-02, but `foamListRegions` is NO check -- it flattens every list regardless of kind. DETECTOR: `grep -A6 'MONOLITHICALLY' log.soPlasmaFoam` / confirm the model line names multiRegionPoisson.
- splitMeshRegions builds the coupling from INTERNAL faces between cellZones. If an interface edge was named as a boundary patch in the .msh, the regions never couple and NOTHING ERRORS. DETECTOR: `grep -A6 "Sizes of interfaces" log.splitMeshRegions | grep -E "^[0-9]"` must return a row -- and use -A6, never -A3.
- `insulatingWall` is `surface conductor`: it absorbs but accumulates NO sigma. A surface meant to charge and declared insulatingWall loses every charge silently. DETECTOR: `grep -l 'enableSurfaceCharging *true' 0/gas/*` should list exactly the patches you expect to charge.
- `speciesSurface open` on an electrode switches the wall flux OFF entirely -- no drift flux, no thermal flux, no secondary emission. It is announced at run time, but only if you read the log. DETECTOR: `grep -n speciesSurface configuration/boundaries` before trusting any wall result; and `grep 'species surface OVERRIDDEN' log.*`.
- changeDictionary STRIPS the `#include` from every field it rewrites, so a `$var` it copied through verbatim becomes undefined and the solver dies reading the field ("Illegal dictionary entry or environment variable name ..."). Re-inserting it before the FoamFile block instead of after gives "problem while reading header". DETECTOR: `grep -L '#include "../../configuration/config"' 0/*/[a-z]*` on any field containing a `$`.
- `set -e` in a script that sources OpenFOAM's etc/bashrc kills it SILENTLY -- no log, no message (measured 2026-09-07: empty logs/ directory, no output at all). DETECTOR: an empty logs/ directory plus exit status != 0 and zero stdout. Fix: per-step `|| die <name>`.
- `set -u` in a script that sources $WM_PROJECT_DIR/bin/tools/RunFunctions is fatal at RunFunctions line 29 (FOAM_LD_LIBRARY_PATH unbound). DETECTOR: the error names FOAM_LD_LIBRARY_PATH, which has nothing to do with your script.
- `export A=x B=$A` on ONE line expands $A before assigning it. SoPLASMA_ETC came out as `/etc` and plasmaSetupBoundaries died on `Cannot read the boundary-role library: "/etc/boundaryRoles"`. DETECTOR: `[ -f "$SoPLASMA_ETC/boundaryRoles" ] || exit 1` immediately after setting it.
- ionmob.py handed the master YAML instead of the compiled .foam writes 0 tables and EXITS 0. DETECTOR: `ls constant/ionTables | wc -l` must be non-zero -- the exit status proves nothing.
- Comparing a SOLVER-SERIALISED field against fresh generator output produces false differences (the BC's own write() adds refValue/refGradient/valueFraction/source/defaultSEEC) or, worse, a FALSE PASS when a backgrounded run has already overwritten the baseline with the generator's own broken output -- which happened while the generator was emitting `type ;` for every patch. DETECTOR: PRINT the generated file, do not only diff it; snapshot the baseline before running anything that could overwrite it; compare generator output to generator output with no solver in between.
- GAMG's default agglomerator aborts under monolithic coupling with "Attempt to cast type lduPrimitiveMeshAssembly to type fvMesh", which reads like an internal error. DETECTOR: if the run dies on that string, add `agglomerator assembledFaceAreaPair;`. Separately, reaching maxIter is NOT an error in OpenFOAM -- ten consecutive unconverged solves were logged as normal and were 7.3% wrong.
- Everything under 0/ is destroyed by `rm -rf 0` at the top of the generation step -- there is no refuse-to-clobber guard yet. DETECTOR: none; the file is simply gone. Put deliberate overrides in etc/changeDictionary.<region>, which runs after both generators and wins.

## Open questions
- Does `plasmaSetupBoundaries -listKinds` actually print all 9 kinds, or only the subset the installed binary knows? The binary at /home/kkourtza/OpenFOAM/kkourtza-v2412/platforms/linux64GccDPInt32Opt/bin/plasmaSetupBoundaries is dated 2026-09-05 19:32 while etc/boundaryRoles was last edited 2026-09-10 12:31. The library is DATA read at runtime so this should be fine, but I could not run the binary (no OpenFOAM env in this session) to confirm. Verify with a single `plasmaSetupBoundaries -listKinds` in a sourced shell.
- `configuration/case` -- the OTHER half of layer 1, named in the agreed architecture (mesh, regions+materials, gas, chemistry, run control) -- does not exist in the tree. Today `configuration/config` plus hand-written constant/regionProperties, constant/plasmaSpeciesProperties, constant/plasmaTransportProperties and system/* fill that role. Is `configuration/case` still planned, or has `configuration/config` absorbed it?
- The `advanced { }` passthrough and the refuse-to-clobber checksum guard are both recorded as design requirements in memory two-layer-case-architecture.md but I found no implementation in either generator. Are they deferred, or dropped?
- The needleDBD `Allrun-serial` change (SoPLASMA_TOOLS derived from $SoPLASMA) is UNCOMMITTED in the working tree. Is it intended to be committed, and should the same fix be applied to the streamer beds' Allrun scripts, which still use bare relative tool paths?
- needleDBD's `simulationType` is dangling: system/plasmaSimulationControls:29 hardcodes `transient` instead of `$simulationType`. Is this a real defect to fix (wire the $var), or was the config key deliberately left as documentation of the choice? Given the file's own rule -- "if it is here, something reads it" -- it reads as a defect.
- tutorials/.../positiveStreamer_fixedMesh/configuration/boundaries records that layer 1 and etc/changeDictionary DISAGREE on electrode absorption (absorbing wall vs zeroGradient) and on the `far` inflow value, and explicitly defers which is the better model of the benchmark. That is an open physics question the architecture cannot settle by itself; a new-case skill must not silently pick one.
- I could not verify the claimed "17 generated files" count for needleDBD directly (the case is Allclean'd -- no 0/ directory present). It comes from memory boundary-role-library-layer1.md dated 2026-09-04, before the 2026-09-10 openBoundary/speciesSurface changes.
