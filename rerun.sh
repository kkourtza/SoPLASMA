#!/usr/bin/env bash
# Re-run ONLY the solver on the already-meshed case.
#
# Allrun-serial rebuilds the mesh from blockMesh through five refinement passes
# and re-seeds the fields every time, which costs minutes and cannot change
# between solver rebuilds. This restores 0/ from the snapshot taken after
# seeding, drops the time directories, and runs soPlasmaFoam alone.
#
# Use ./setup-case.sh + ./Allrun-serial only when the MESH or the CASE SETUP
# changes (species list, seed, BCs); use this for anything that is a code change.
set -e
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1 || true
export SoPLASMA=$HOME/soplasma-scratch SoPLASMA_SRC=$HOME/soplasma-scratch/src SoPLASMA_ETC=$HOME/soplasma-scratch/etc
CASE=$HOME/streamer-case
cd "$CASE"
[ -d .snapshot0 ] || { echo "no .snapshot0 -- run Allrun-serial once first"; exit 1; }
foamListTimes -rm -withZero 2>/dev/null || true
rm -rf 0 && cp -r .snapshot0 0
rm -f log.soPlasmaFoam
timeout "${TIMEOUT:-1800}" soPlasmaFoam > log.soPlasmaFoam 2>&1 || true
echo "steps: $(grep -c '^Time' log.soPlasmaFoam)   last: $(tail -1 log.soPlasmaFoam)"
