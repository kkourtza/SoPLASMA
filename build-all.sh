#!/usr/bin/env bash
# Serial build of the whole tree in dependency order. Logs per component so a
# failure is attributable without re-running everything.
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1
export SoPLASMA=$HOME/soplasma-scratch
export SoPLASMA_SRC=$HOME/soplasma-scratch/src
export SoPLASMA_ETC=$HOME/soplasma-scratch/etc
cd "$SoPLASMA" || exit 1

# Never rebuild over a library a running solver has mapped. See the script for
# what that costs and how to override it.
./check-no-running-solvers.sh || exit 1

mkdir -p buildlogs
failed=0

# Dependency order matters for the libraries; the applications come last.
BUILD_DIRS=(
  src/profilers
  src/constants
  src/numerics
  src/models/electromagnetics
  src/models/plasmaModels/genericPlasmaProperties
  src/models/plasmaModels/plasmaSpecies
  src/models/plasmaModels/photoionization
  src/models/plasmaModels/plasmaBoltzmann
  src/models/plasmaModels/plasmaEnergy
  src/models/plasmaModels/plasmaChemistry
  src/models/plasmaModels/plasmaReactionRates
  src/models/plasmaModels/plasmaTransport
  src/bcs
  src/tools
  src/applications/utilities/foamPlasmaCreateSpeciesFields
  src/applications/utilities/plasmaChemistry0D
  src/applications/utilities/testWallLoss
  src/applications/utilities/testVibRelax
  src/applications/utilities/testDischargeCurrent
  src/applications/utilities/testAitken
  src/applications/solvers/soPlasmaFoam
  src/applications/solvers/singleRegionElectrostaticFoam
  src/applications/solvers/multiRegionElectrostaticFoam
  ThirdParty/libROUNDSchemes
)

# Applications deliberately NOT built, each with the reason. An entry here is a
# decision; an application in neither list is an OVERSIGHT, and the coverage
# check below refuses to let one pass silently.
SKIP_DIRS=()

for d in "${BUILD_DIRS[@]}"; do
  n=$(basename "$d")
  if [ -d "$d" ]; then
    ( cd "$d" && { [ -x ./Allwmake ] && ./Allwmake || wmake libso || wmake; } ) \
      > "$SoPLASMA/buildlogs/$n.log" 2>&1 \
      && echo "OK   $n" || { echo "FAIL $n"; failed=1; }
  fi
done

# ---------------------------------------------------------------- coverage --
#
# AN APPLICATION THAT IS NOT BUILT IS NOT TESTED, AND IT ROTS SILENTLY.
#
# testAitken was absent from this list and went stale for nine days: the binary
# was built 2026-08-21 19:12, libplasmaNumerics.so was rebuilt at 20:32 against
# a changed aitkenRelaxation.H, and the ABI mismatch SEGFAULTED the unit test
# mid-run -- free() on 0x3ff0000000000000, the bit pattern of 1.0. Nothing
# reported it, because nothing was looking. Both electrostatic solvers had
# NEVER been built at all. All three compile clean; they were simply forgotten
# when this list was last hand-edited.
#
# So the list is no longer allowed to be the only record. Every directory under
# src/applications that carries a Make/ must be in BUILD_DIRS or, deliberately,
# in SKIP_DIRS. Anything in neither is reported and fails the build.
missing=0
while IFS= read -r mk; do
  d=$(dirname "$mk")
  case " ${BUILD_DIRS[*]} ${SKIP_DIRS[*]} " in
    *" $d "*) ;;
    *)
      echo "UNCOVERED  $d  -- in neither BUILD_DIRS nor SKIP_DIRS"
      missing=1
      ;;
  esac
done < <(find src/applications -mindepth 2 -maxdepth 3 -type d -name Make | sort)

if [ "$missing" -ne 0 ]; then
  echo "BUILD-INCOMPLETE: an application is in neither list. Add it to"
  echo "  BUILD_DIRS, or to SKIP_DIRS with the reason. A binary that is never"
  echo "  rebuilt will segfault against a changed library ABI without warning."
  exit 1
fi

# A COMPONENT THAT FAILED MUST NOT LOOK LIKE A COMPLETE BUILD.
#
# This script printed "BUILD-COMPLETE" even when a component reported FAIL, and
# a whole session was spent grepping for that string as the build check. That is
# how you end up running a solver against a stale library and debugging the
# wrong thing.
if [ "$failed" -ne 0 ]; then
  echo "BUILD-FAILED: at least one component did not build. See buildlogs/."
  exit 1
fi

echo "BUILD-COMPLETE"
