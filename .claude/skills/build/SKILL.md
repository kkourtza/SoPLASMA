---
name: build
description: Build SoPlasma and PROVE the libraries actually rebuilt. Exit 0, "Build complete." and even "BUILD-COMPLETE" do not prove a library recompiled. Carries the only correct build command, the four-part pass criterion, the mtime/ABI staleness proof, the lnInclude rule and the single-library --fast form. Use before any run whose result you will believe, and after every .H edit.
allowed-tools: Bash, Read
---
# /build — build SoPlasma, then prove it

**Question this answers:** *is the binary a run will load actually built from the source I
just edited?* Nothing else answers it. Seven recorded failures (taxonomy CLASS 4) were
phantom physics bugs — `bad_alloc` in `soPlasmaFoam`, SEGV in
`plasmaTimeControl::adjustDeltaT` inside a library never edited, SIGFPE in
`mechanismSourceTerms` — every one a stale `.so`. A control run against a stale binary is
not a control (A1). **Not for:** SoEEDF-only work (step 9); asking what the build *does* —
read `build-all.sh` (A4). **Never** build with `./Allwmake`, a subdirectory `wmake`, or
`wmake libso src/numerics` (step 6).

## 1. Preflight — no running solver, and BOTH bashrcs with the project SECOND
```bash
/home/kkourtza/soplasma-scratch/check-no-running-solvers.sh    # PASS: exit 0, no output
pgrep -x soPlasmaFoam ; pgrep -af "[p]lasmaChemistry0D" ; pgrep -af "[E]lectrostaticFoam"
source /usr/lib/openfoam/openfoam2412/etc/bashrc && source $HOME/soplasma-scratch/etc/bashrc
echo "$PETSC_DIR | $PETSC_ARCH" ; [ -d "$PETSC_DIR/$PETSC_ARCH/lib" ] && echo PETSC-OK
```
PASS on the last line: `/home/kkourtza/soplasma-scratch/ThirdParty/petsc-3.24.0 | DPInt32`
and `PETSC-OK`. `build-all.sh` runs the guard itself and **exits 1 with no compile error at
all** when it refuses (`REFUSING TO BUILD: a SoPLASMA executable is running.`); it cost a
scare on 2026-08-19, a rebuild over a 45-min streamer at step 230. It is **partly
inoperative**, hence line 2: `pgrep -x` is capped at procps' 15-char comm name, so of its
seven executables only `soPlasmaFoam`, `testWallLoss`, `testVibRelax` can ever match
(`plasmaChemistry0D` 17, `testChemistryBackends` 21, `plasmaCreateSpeciesFields` 25 never
do), and the electrostatic solvers and other beds are not listed at all; `pgrep -af
soPlasmaFoam` matches the agent's own shell. Kill with `pkill -9 -x soPlasmaFoam` **as its
own simple command** (`pkill -f` kills the calling shell). `FORCE=1 ./build-all.sh`
overrides, printing nothing extra when nothing runs — not a pass string. Bashrc order
matters: the project's is the ONLY place `PETSC_DIR`/`PETSC_ARCH` are set and must override
OpenFOAM's, which points at a nonexistent path here
(`/usr/lib/openfoam/openfoam2412/ThirdParty` is a stub FILE). Wrong order and
`libplasmaNewtonSolverPETSc` goes stale **silently** — nothing links it, it is `dlopen`'d
from the case's `controlDict libs`, and the run reports `Registered types: 0()`; that is
what segfaulted on 2026-09-10 after a new virtual reached `plasmaTransportModel`.
`etc/bashrc` also prepends `$PETSC_DIR/$PETSC_ARCH/lib` to `LD_LIBRARY_PATH` (needed at RUN
time too); it does NOT put binaries on PATH — OpenFOAM's does.

## 2. The build
```bash
cd $HOME/soplasma-scratch && ./build-all.sh 2>&1 | tee /tmp/build.log
```
~3 min incremental (measured 2026-09-11, buildlogs 17:15→17:18), 15–20 min clean — background
it and poll `buildlogs/` (B3). **Pipe through `tee`**: `OK`/`FAIL` go to build-all.sh's
*stdout*, so `grep FAIL buildlogs/*` can only ever return 0 (verified: 0 in all 34 logs), and
the repo's `build.log` is a stale **Aug-10** artefact reporting a month-old build as current.

## 3. Pass criterion — four greps on the captured log, never an exit code
```bash
grep -cE '^FAIL |BUILD-FAILED|BUILD-INCOMPLETE|^UNCOVERED' /tmp/build.log  # must be 0
grep -c '^BUILD-COMPLETE' /tmp/build.log                                   # must be 1
grep -c '^OK ' /tmp/build.log                                              # must be 31
grep -lE 'Error 1|cannot find -l' $HOME/soplasma-scratch/buildlogs/*.log    # must be empty
```
* `BUILD_DIRS` has exactly 31 entries (verified). Line 78 is `if [ -d "$d" ]` **with no else
  branch**: a renamed or typo'd entry is skipped in total silence and `BUILD-COMPLETE` still
  prints, and the `OK` count is the only check that catches it.
* A component log can fail with `Error 1` / `cannot find -l<lib>` and contain zero `error:`
  lines — hence the fourth grep.
* `./build-all.sh | grep BUILD-COMPLETE` returns **grep's** status and
  `./build-all.sh > log; echo $?` **echo's**; a `BUILD-FAILED` run was called passing on that.

## 4. THE PROOF: exit 0 does not prove a library rebuilt (A1)
Per component you edited (`<srcdir>` e.g. `src/models/plasmaModels/plasmaEnergy`):
```bash
LIBBIN=/home/kkourtza/OpenFOAM/kkourtza-v2412/platforms/linux64GccDPInt32Opt/lib
find <srcdir> -name '*.[CH]' -newer $LIBBIN/lib<X>.so           # PASS: returns NOTHING
nm -DC $LIBBIN/lib<X>.so | grep '<a symbol you just added>'     # PASS: present
```
Measured 2026-09-08: `FORCE=1 ./Allwmake` exited 0 printing nothing while
`libplasmaEnergy.so`'s mtime was UNCHANGED and two siblings rebuilt in the same pass. A shared
library also tolerates undefined symbols until load — a definition that failed to insert still
gave `BUILD-COMPLETE`, then `symbol lookup error: undefined symbol: ...` at runtime; hence
`nm`. Stale library: `cd <srcdir> && wclean && wmake libso .`, then re-run step 2. Directory
name ≠ library name for four components (`genericPlasmaProperties`→`libplasmaProperties`,
`constants`→`libSoPLASMAConstants`, `bcs`→`libplasmaBcs`, `tools`→`libplasmaTools`):
```bash
cd $HOME/soplasma-scratch && find src ThirdParty -maxdepth 4 -path '*/Make/files' \
 -exec grep -l '^LIB' {} + | while read f; do echo "$(dirname $(dirname $f)) -> $(grep -h '^LIB' $f)"; done
```
`libpetscFoam.so` is NOT in `$FOAM_USER_LIBBIN` — petsc4Foam installs with
`-prefix=openfoam` to `/usr/lib/openfoam/openfoam2412/platforms/linux64GccDPInt32Opt/lib/`,
and it has failed silently in a reference run before. Check it there.

## 5. New header or new source file (cost: four build iterations, 2026-09-02)
Both halves are required; the error names the **including** header, not the missing symlink:
```bash
wmakeLnInclude -update $HOME/soplasma-scratch/src/<component>
ls $HOME/soplasma-scratch/src/<component>/lnInclude/<new>.H           # must resolve
grep -n lnInclude $HOME/soplasma-scratch/src/<consumer>/Make/options  # -I$(SoPLASMA_SRC)/<component>/lnInclude
```
A new `.C` also goes in that component's `Make/files` above the `LIB =` line; a new
application goes in `BUILD_DIRS` (or `SKIP_DIRS` with a written reason) **in the same
commit**, or the coverage gate prints `UNCOVERED` and fails the build.

## 6. When a FULL build is mandatory — the ABI rule (C1)
A **new virtual** in a shared header shifts every later vtable slot; a **new member** changes
object size and every later offset. Libraries not rebuilt call the wrong slot or read garbage,
and the crash lands in a library you never touched. `free()` on `0x3ff0000000000000` (the bit
pattern of 1.0) is the ABI-mismatch signature, not heap corruption in your new code; for a
diagnostic counter prefer a file-local `static` in the `.C` to a member on a widely-included
class. `./Allwmake` is NOT a full build and prints `Build complete.` anyway:
`src/models/plasmaModels/Allwmake` has **four** `wmake` lines for **eight** model dirs
(plasmaBoltzmann, plasmaChemistry, plasmaEnergy, plasmaReactionRates never built) and
`src/applications/Allwmake` builds 3 solvers + 1 utility and **no test bed**. `wmake
src/numerics` builds `libplasmaNumerics` and NOT `libplasmaNewtonSolverPETSc` — separate
target, exits 0 having done nothing.

## 7. Outside `build-all.sh` — build by hand (after step 1)
```bash
wmake $HOME/soplasma-scratch/src/applications/utilities/testSnesJFNK      # SKIP_DIRS, deliberate
wmake $HOME/soplasma-scratch/src/applications/utilities/testSnesJFNK2Field
wmake $HOME/soplasma-scratch/applications/test/testChemistryBackends      # separate tree, in NO list
```
The coverage gate walks only `find src/applications -mindepth 2 -maxdepth 3 -type d -name
Make`, so the top-level `applications/` tree is invisible to it: `testChemistryBackends`'s
binary is Aug 12 against a Sep 11 `libplasmaTransport.so` — the staleness pattern that
segfaulted testAitken, so do not trust its `Backends agree.` until you rebuild it. And
`src/tools/decompositionMethods` (`libplasmaDecompositionMethods`) is in no list, never built.

## 8. --fast (B6) — one library, and its hard limit
```bash
wmake libso $HOME/soplasma-scratch/src/<component>        # seconds to a minute, after step 1
```
Legitimate ONLY for a `.C`-only edit inside a leaf library you are still iterating on. **It
may NEVER be used to validate a behaviour change, and no conclusion from a run on a partial
build may be carried forward** (A7). Any `.H` layout or virtual change, and any run you will
quote, requires step 2 in full first.

## 9. SoEEDF (standalone; C1 does not apply there)
```bash
cd $HOME/Projects/SoEEDF && cmake -S . -B build && cmake --build build -j$(nproc)
```
10 targets (no `validateBiagi` — if a doc names one, the doc is wrong). Every binary and the
newest source are stamped 2026-08-20 17:47, so a build today is a **no-op**: check a
binary's mtime advanced before reporting "rebuilt". `build/`'s cache is keyed to the
`/home/kkourtza/Projects/BoltzmannSolver` symlink (`grep CMAKE_HOME_DIRECTORY
build/CMakeCache.txt`), so a clean configure needs `rm -rf build` first; run the binaries
from `build/` — they load cross sections by `../data/...`.
