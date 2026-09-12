#!/usr/bin/env bash
# warm_dt_ladder.sh <solver> <dt> [kspMaxIt]
#
# THE WARM-STARTED dt LADDER — the experiment the Newton-vs-Picard benchmark
# exists for (PROGRESS.md Task 1).
#
# IN tools/ AND NOT THE SCRATCHPAD ON PURPOSE: the scratchpad lives under /tmp,
# which a reboot wipes. This harness was lost exactly that way on 2026-09-12
# after a WSL freeze, along with every other ladder script.
#
# PARAMETERS (B3), restated for every arm:
#   bed        validation/warm449 @ t=1e-09 — 449,413 cells, DEVELOPED streamer
#              (peak n_e 1.14e19 m^-3 against a 1e13 seed)
#   start      t=1e-09, the IDENTICAL snapshot for every arm
#   end        t=2e-09 — a COMMON endpoint (A2) and the validated benchmark window
#   dt         fixed, adjustTimeStep false
#   limiters   limitSpeciesCo / limitChemistryCo / limitDielectricRelaxationRatio
#              ALL false — ASSERTED below, never assumed. Removing them is the
#              point: they exist because a SEGREGATED scheme cannot step over the
#              timescales they bound, and Newton's claim is that it need not.
#   newton     kspMaxIt (default 1000; the hardcoded 100 produced every prior
#              negative verdict), assembledPmat true, rtol 1e-8, maxIt 50
#
# WHY WARM: at t=1e-09 the dielectric relaxation time tau = eps0/(e mu_e n_e) is
# ~1.2e-10 s, so dt/tau = 1 near dt=1.2e-10 — THE STIFF REGIME. Every cold ladder
# ran at tau ~ 1e-4 s, seven orders away, and could not test the claim at all.
# Cold also cannot hand over to Newton: handover needs a completed Picard step.
#
# RESOURCE NOTE: this bed is ~1.5 GB resident per solver on a 30 GB machine. Ten
# concurrent arms froze WSL twice on 2026-09-12. Run a FEW at a time.
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1 || true
export SoPLASMA=$HOME/soplasma-scratch
export SoPLASMA_ETC=$SoPLASMA/etc
export SoPLASMA_SRC=$SoPLASMA/src
source $SoPLASMA/etc/bashrc >/dev/null 2>&1 || true
V=$HOME/soplasma-scratch/validation
s=$1; dt=$2; KMI=${3:-1000}
# BED selects the warm bed. Default is the 449k physics bed; the COARSE bed
# (81,640 cells, $HOME/streamer-warm) reproduces the SAME stiff regime -- peak
# n_e 1.19e19 vs 1.14e19, tau 1.16e-10 vs 1.2e-10 -- at 5.5x less cost per step
# and far less memory, so NUMERICS questions belong there (A6/B6). It does NOT
# resolve the streamer and must never be quoted for physics.
W=${BED:-$V/warm449}
tag=$(basename $W)
name=wlad_${tag}_${s}_dt${dt}
[ -d "$W/1e-09" ] || { echo "FAILED: no warm state at $W/1e-09"; exit 1; }

# REFUSE to clobber a live arm. `rm -rf` on a case directory whose solver is
# still running deletes the tree UNDER the process: the old solver goes on
# writing deleted file handles while the new one writes fresh ones, both under
# the same name, and the log reads as empty. Measured 2026-09-12 -- two
# soPlasmaFoam processes were found in wlad_streamer-warm_newton_dt1e-10. Same
# class as the documented "two mpirun instances writing the same processor* dirs".
for pid in $(pgrep -x soPlasmaFoam 2>/dev/null); do
  if [ "$(readlink /proc/$pid/cwd 2>/dev/null)" = "$V/$name" ]; then
    echo "FAILED: a solver (pid $pid) is ALREADY running in $V/$name."
    echo "        Kill it first, or pick another dt. Refusing to rm -rf under a live process."
    exit 1
  fi
done

# TRIMMED copy: the warm bed is 1.7 GB because it holds ten snapshots and an arm
# needs exactly one. NOT `cp -al` — a hardlinked tree is corrupted by Python
# open(w), which the newtonSolver edit below uses, silently rewriting every sibling.
rm -rf $V/$name; mkdir -p $V/$name
for d in constant system configuration 0.orig etc; do
  [ -e "$W/$d" ] && cp -r "$W/$d" "$V/$name/"
done
cp -r "$W/1e-09" "$V/$name/1e-09"
cd $V/$name || exit 1

sed -i -e "s|^deltaT .*|deltaT                              $dt;|" \
       -e "s|^endTime .*|endTime                             2e-09;|" \
       -e "s|^writeControl .*|writeControl                        timeStep;|" \
       -e "s|^writeInterval .*|writeInterval                       1000000;|" configuration/config
grep -q "^deltaT  *$dt;"    configuration/config || { echo "FAILED: deltaT";  exit 1; }
grep -q '^endTime  *2e-09;' configuration/config || { echo "FAILED: endTime"; exit 1; }
sed -i 's|^startFrom .*|startFrom       latestTime;|' system/controlDict
grep -q '^startFrom       latestTime;' system/controlDict || { echo "FAILED: startFrom"; exit 1; }

for k in limitSpeciesCo limitChemistryCo limitDielectricRelaxationRatio adjustTimeStep; do
  grep -q "^$k  *false;" configuration/config || { echo "FAILED: $k is not false"; exit 1; }
done

if grep -q '^outerSolver' configuration/config; then
  sed -i "s|^outerSolver .*|outerSolver                         $s;|" configuration/config
else printf '\nouterSolver                         %s;\n' "$s" >> configuration/config; fi
grep -q "^outerSolver  *$s;" configuration/config || { echo "FAILED: outerSolver"; exit 1; }

if [ "$s" = newton ]; then
  ~/ct-env/bin/python - "$KMI" <<'PYEOF'
import re, sys
kmi = sys.argv[1]
p='system/plasmaSimulationControls'; s=open(p).read()
blk = ('    newtonSolver      { type SNES; rtol 1e-8; maxIt 50; bounded false; '
       'kspMaxIt %s; assembledPmat true; petscOptions "-ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor"; }' % kmi)
if 'newtonSolver' in s:
    s = re.sub(r'^ *newtonSolver .*$', blk, s, count=1, flags=re.M)
else:
    s = s.replace('    target            converged;',
                  '    target            converged;\n    outerSolver       $outerSolver;\n' + blk, 1)
open(p,'w').write(s)
PYEOF
  grep -q "kspMaxIt $KMI;" system/plasmaSimulationControls || { echo "FAILED: kspMaxIt"; exit 1; }
  if ! grep -q 'libplasmaNewtonSolverPETSc' system/controlDict; then
    ~/ct-env/bin/python -c "
import re
p='system/controlDict'; s=open(p).read()
s=re.sub(r'(^libs\s*\()', r'\1 \"libplasmaNewtonSolverPETSc.so\"', s, count=1, flags=re.M)
open(p,'w').write(s)"
  fi
  grep -q 'libplasmaNewtonSolverPETSc' system/controlDict || { echo "FAILED: newton lib"; exit 1; }
fi

echo "  $name VERIFIED: bed=$tag $s dt=$dt start=1e-09 end=2e-09 limiters OFF kspMaxIt=$([ "$s" = newton ] && echo $KMI || echo n/a)"
timeout 14400 soPlasmaFoam > $V/$name.log 2>&1
echo "WLAD-EXIT=$? $s dt=$dt steps=$(grep -c '^Time = ' $V/$name.log) last=$(grep '^Time = ' $V/$name.log | tail -1)"
grep -q 'sigFpe\|sigSegv' $V/$name.log && echo "  *** CRASHED (SIGFPE/SEGV) ***"
