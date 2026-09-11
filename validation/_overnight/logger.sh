#!/usr/bin/env bash
# Overnight status log. One line per arm every 10 min -> _overnight/STATUS.txt
cd /home/kkourtza/soplasma-scratch/validation
OUT=_overnight/STATUS.txt
printf "%-19s %-22s %-9s %-10s %-7s %-7s %-12s %s\n" \
  DATE CASE t_us n_e_max Te_max Te_min dt flags >> $OUT
while true; do
  for c in grubert2009_cfs_r036 grubert2009_cfs_250V grubert2009_cfs_200V; do
    L=$c/logs/log.soPlasmaFoam
    [ -f "$L" ] || continue
    t=$(grep -aE '^Time = ' "$L" | tail -1 | awk '{print $3}')
    dt=$(grep -aE 'current deltaT' "$L" | tail -1 | awk '{print $3}')
    nc=$(grep -ac 'PIMPLE: not converged' "$L")
    ab=$(grep -ac 'mpi-abort\|FOAM FATAL' "$L")
    cl=$(grep -ac 'pinned at the' "$L")
    n=0
    for p in $(pgrep soPlasmaFoam 2>/dev/null); do
      cw=$(readlink /proc/$p/cwd 2>/dev/null)
      case "$cw" in *$c*) n=$((n+1));; esac
    done
    line=$(python3 ../tools/news.py "$c" 2>/dev/null | tail -1)
    ne=$(echo "$line" | awk '{print $3}')
    te=$(echo "$line" | awk '{print $4}')
    ti=$(echo "$line" | awk '{print $5}')
    tus=$(python3 -c "print('%.4f'%(float('$t')*1e6))" 2>/dev/null)
    printf "%-19s %-22s %-9s %-10s %-7s %-7s %-12s nc=%s ranks=%s ab=%s cl=%s\n" \
      "$(date '+%Y-%m-%d %H:%M:%S')" "$c" "$tus" "$ne" "$te" "$ti" "$dt" \
      "$nc" "$n" "$ab" "$cl" >> $OUT
  done
  echo "" >> $OUT
  sleep 600
done
