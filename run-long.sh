#!/usr/bin/env bash
# Long unattended run, with the log cap that a 9.7 GB runaway earned us.
#   ./run-long.sh <caseDir>
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1 || true
export SoPLASMA=$HOME/soplasma-scratch SoPLASMA_SRC=$HOME/soplasma-scratch/src SoPLASMA_ETC=$HOME/soplasma-scratch/etc
CASE=${1:?usage: run-long.sh <caseDir>}
cd "$CASE" || exit 1
NAME=$(basename "$CASE")

# Mesh + fields first; its own log so a setup failure is not buried in the run.
./Allclean > /dev/null 2>&1 || true
{
  ./Allrun-serial 2>&1 | head -c $((200*1024*1024))
} > "setup.log" 2>&1
echo "[$NAME] setup done, $(grep -c '^Time' log.soPlasmaFoam 2>/dev/null || echo 0) steps so far"
