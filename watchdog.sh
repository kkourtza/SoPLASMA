#!/usr/bin/env bash
# Kill the solver if its log exceeds a size cap. Safety net for unattended runs:
# a runaway warning loop previously wrote 9.7 GB before one timestep completed.
LOG=${1:?log}; MAXMB=${2:-50}
while true; do
  [ -f "$LOG" ] || { sleep 2; continue; }
  sz=$(stat -c%s "$LOG" 2>/dev/null || echo 0)
  if [ "$sz" -gt $((MAXMB*1024*1024)) ]; then
    echo "WATCHDOG: $LOG exceeded ${MAXMB}MB -- killing solver"
    pkill -9 -x soPlasmaFoam
    exit 1
  fi
  pgrep -x soPlasmaFoam >/dev/null || exit 0
  sleep 3
done
