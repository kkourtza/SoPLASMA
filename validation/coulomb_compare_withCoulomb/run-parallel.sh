#!/usr/bin/env bash
cd "${0%/*}" || exit 1
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1 || true
export SoPLASMA=$HOME/soplasma-scratch
export SoPLASMA_ETC=$SoPLASMA/etc
export SoPLASMA_SRC=$SoPLASMA/src
NP=$(grep -oE 'numberOfSubdomains +[0-9]+' configuration/config | awk '{print $2}')
mpirun -np "$NP" soPlasmaFoam -parallel 2>&1 \
  | grep -vE "^(DILUPBiCGStab|smoothSolver|GAMG|PBiCGStab|PCG):|^Discretizing transport|^Transport models corrected|^Species. number densities clamped|^Charge density updated|^Surface charge updated|^PIMPLE: iteration" \
  > logs/log.soPlasmaFoam
