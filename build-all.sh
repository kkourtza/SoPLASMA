#!/usr/bin/env bash
# Serial build of the whole tree in dependency order. Logs per component so a
# failure is attributable without re-running everything.
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1
export SoPLASMA=$HOME/soplasma-scratch
export SoPLASMA_SRC=$HOME/soplasma-scratch/src
export SoPLASMA_ETC=$HOME/soplasma-scratch/etc
cd "$SoPLASMA" || exit 1
mkdir -p buildlogs
for d in src/profilers src/constants src/numerics src/models/electromagnetics \
         src/models/plasmaModels/genericPlasmaProperties \
         src/models/plasmaModels/plasmaSpecies \
         src/models/plasmaModels/photoionization \
         src/models/plasmaModels/plasmaBoltzmann \
         src/models/plasmaModels/plasmaChemistry \
         src/models/plasmaModels/plasmaReactionRates \
         src/models/plasmaModels/plasmaTransport \
         src/bcs src/tools \
         src/applications/utilities/foamPlasmaCreateSpeciesFields \
         src/applications/utilities/plasmaChemistry0D \
         src/applications/solvers/soPlasmaFoam \
         ThirdParty/libROUNDSchemes; do
  n=$(basename "$d")
  if [ -d "$d" ]; then
    ( cd "$d" && { [ -x ./Allwmake ] && ./Allwmake || wmake libso || wmake; } ) \
      > "$SoPLASMA/buildlogs/$n.log" 2>&1 \
      && echo "OK   $n" || echo "FAIL $n"
  fi
done
echo "BUILD-COMPLETE"
