#!/usr/bin/env bash
# PARALLEL launch. Verified against serial 2026-09-07: n_e,max, Te,max and Te,min
# agree to 7 digits at matched simulated times (np2, np4, np8). Note that
# adaptive dt makes the SAMPLE POINTS drift apart between runs even when the
# trajectories agree -- always compare at matched time, never at matched step.
cd "${0%/*}" || exit 1
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1 || true
export SoPLASMA=$HOME/soplasma-scratch
export SoPLASMA_ETC=$SoPLASMA/etc
export SoPLASMA_SRC=$SoPLASMA/src
NP=$(grep -oE 'numberOfSubdomains +[0-9]+' configuration/config | awk '{print $2}')
mpirun -np "$NP" soPlasmaFoam -parallel 2>&1 \
  | grep -vE "^(DILUPBiCGStab|smoothSolver|GAMG|PBiCGStab|PCG):|^Discretizing transport|^Transport models corrected|^Species. number densities clamped|^PIMPLE: iteration" \
  > logs/log.soPlasmaFoam
