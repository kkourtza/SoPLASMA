#!/usr/bin/env bash
# Filtered launch. The per-solve residual lines run ~12 kB/step, tens of GB
# over a run this long, and carry nothing the diagnostics block and
# postProcessing/dischargeCurrent/current.csv do not already record.
cd "${0%/*}" || exit 1
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1
export SoPLASMA=$HOME/soplasma-scratch SoPLASMA_ETC=$SoPLASMA/etc
soPlasmaFoam 2>&1 \
  | grep -vE "^(DILUPBiCGStab|smoothSolver|GAMG|PBiCGStab|PCG):|^Discretizing transport|^Transport models corrected|^Species. number densities clamped|^PIMPLE: iteration" \
  > logs/log.soPlasmaFoam
