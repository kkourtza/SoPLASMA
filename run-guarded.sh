#!/usr/bin/env bash
# Run the solver with a HARD log cap. A runaway warning loop once wrote 9.7 GB
# and 91M lines before a single timestep completed, which took the machine
# down; nothing unattended should be able to do that again.
#
#   ./run-guarded.sh [maxLogMB] [maxSeconds]
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1 || true
export SoPLASMA=$HOME/soplasma-scratch SoPLASMA_SRC=$HOME/soplasma-scratch/src SoPLASMA_ETC=$HOME/soplasma-scratch/etc
CASE=$HOME/streamer-case
MAXMB=${1:-50}; MAXSEC=${2:-600}
cd "$CASE" || exit 1

# head -c caps the log at MAXMB and closes the pipe, so the solver dies of
# SIGPIPE rather than filling the disk. timeout caps wall clock independently.
( timeout "$MAXSEC" soPlasmaFoam 2>&1 | head -c $((MAXMB*1024*1024)) > log.solver ) 
rc=$?
echo "RUN-EXIT=$rc  log=$(du -h log.solver 2>/dev/null | cut -f1)  lines=$(wc -l < log.solver)"
