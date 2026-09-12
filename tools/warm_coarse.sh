#!/usr/bin/env bash
# warm_coarse.sh -- warm state on the COARSE streamer bed, for NUMERICS questions.
#
# WHY COARSE: "does the outer Krylov converge at dt/tau > 1 on a developed
# streamer" is a NUMERICS question and does not need 449k cells. A6/B6: debug on
# the smallest bed that reproduces it. All four Newton streamer defects were found
# on a coarse bed after hours were lost at minutes-per-attempt on 1.15M.
# NEVER quote this bed for physics -- it does not resolve the streamer.
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1 || true
export SoPLASMA=$HOME/soplasma-scratch
export SoPLASMA_ETC=$SoPLASMA/etc
export SoPLASMA_SRC=$SoPLASMA/src
source $SoPLASMA/etc/bashrc >/dev/null 2>&1 || true
C=$HOME/streamer-warm
cd $C || exit 1
rm -rf [0-9]*e-* log.warm
sed -i -e 's|^deltaT .*|deltaT                              1e-11;|' \
       -e 's|^endTime .*|endTime                             1e-09;|' \
       -e 's|^writeControl .*|writeControl                        timeStep;|' \
       -e 's|^writeInterval .*|writeInterval                       10;|' configuration/config
# writeControl MUST be timeStep: under `runTime` an interval of 10 means TEN
# SECONDS of simulated time and this run ends at 1e-09 s, so NOTHING is written
# while every grep-the-key check still passes (measured 2026-09-12).
grep -q '^writeControl  *timeStep;' configuration/config || { echo "FAILED: writeControl"; exit 1; }
grep -q '^deltaT  *1e-11;'          configuration/config || { echo "FAILED: deltaT"; exit 1; }
grep -q '^endTime  *1e-09;'         configuration/config || { echo "FAILED: endTime"; exit 1; }
if grep -q '^outerSolver' configuration/config; then
  sed -i 's|^outerSolver .*|outerSolver                         picard;|' configuration/config
else printf '\nouterSolver                         picard;\n' >> configuration/config; fi
grep -q '^outerSolver  *picard;' configuration/config || { echo "FAILED: outerSolver"; exit 1; }
for k in limitSpeciesCo limitChemistryCo limitDielectricRelaxationRatio adjustTimeStep; do
  grep -q "^$k  *false;" configuration/config || { echo "FAILED: $k not false"; exit 1; }
done
echo "  warm_coarse VERIFIED: picard dt=1e-11 -> 1e-09, writeControl timeStep, limiters OFF"
timeout 3600 soPlasmaFoam > $C/log.warm 2>&1
n=$(ls -d $C/[0-9]*e-* 2>/dev/null | wc -l)
echo "WARMC-EXIT=$? steps=$(grep -c '^Time = ' $C/log.warm) last=$(grep '^Time = ' $C/log.warm | tail -1) snapshots=$n"
[ "$n" -ge 5 ] || echo "  *** FAILED: expected >=5 snapshots, got $n ***"
