#!/usr/bin/env bash
# Cut the 200k arm of the mesh sweep down to 20 relax + 10 measured.
#
# The sweep script uses one NRELAX/NSTEPS for every arm, and editing a RUNNING
# bash script is unsafe (bash reads it incrementally). So instead: let the sweep
# build the 200k mesh and map the fields -- that work is needed either way --
# then stop only its solver, shorten endTime, and rerun. Nothing is wasted
# except the few steps it managed before being stopped.
#
# No `set -u`: OpenFOAM's bashrc trips it (see mesh_sweep.sh).
cd "${0%/*}" || exit 1
SOP=$HOME/soplasma-scratch
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1
source "$SOP"/etc/bashrc >/dev/null 2>&1
export SoPLASMA=$SOP SoPLASMA_SRC=$SOP/src SoPLASMA_ETC=$SOP/etc

SRCT=1.954497e-06
DT=1e-12
NRELAX=20
NSTEPS=10

echo "waiting for the sweep to reach the 200k arm..."
while [ ! -f mesh_200000/logs/log.run ]; do
  # if the sweep driver is gone and 200k never started, give up rather than spin
  pgrep -x soPlasmaFoam >/dev/null || { [ -d mesh_200000 ] || { echo "sweep ended without a 200k arm"; exit 1; }; }
  sleep 20
done
echo "200k arm started; letting the mesh + mapFields settle"

# kill ONLY the 200k solver, found by working directory -- never a pattern kill
# (pkill -f matches the killer's own command line; that mistake cost three
#  self-kills and a set of killed monitor tasks on 2026-09-10)
for p in $(pgrep -x soPlasmaFoam); do
  d=$(basename "$(readlink /proc/$p/cwd 2>/dev/null)" 2>/dev/null)
  [ "$d" = "mesh_200000" ] && { echo "stopping the long 200k run (pid $p)"; kill "$p"; }
done
# and stop the sweep driver so it does not run its own analysis on a half arm
for p in $(pgrep -u "$USER" -f mesh_sweep.sh); do
  cl=$(tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null)
  case "$cl" in *pgrep*|*eval*|*cut200k*) continue ;; esac
  case "$cl" in *mesh_sweep.sh*) echo "stopping sweep driver (pid $p)"; kill "$p" ;; esac
done
sleep 5

cd mesh_200000 || exit 1
NEW=$(python3 -c "print(f'{float('$SRCT')+($NRELAX+$NSTEPS)*$DT:.9e}')")
sed -i "s/^endTime  *.*/endTime                             $NEW;/" configuration/config
echo "endTime cut to $NEW  ($NRELAX relax + $NSTEPS measured)"
rm -f logs/log.run
/usr/bin/time -f "    wall %e s  maxRSS %M kB" \
  soPlasmaFoam > logs/log.run 2>logs/log.time
echo "exit $?"
tail -2 logs/log.time
echo "NEWTON ENGAGED: $(grep -ac 'outerSolver newton (SNES)' logs/log.run) SNES solves"
cd ..

echo
echo "############ MESH-INDEPENDENCE RESULT ############"
NRELAX=$NRELAX python3 - <<'PY'
import re, glob, os, statistics as s
rows=[]
for d in sorted(glob.glob("mesh_*"), key=lambda x:int(x.split('_')[1]) if x.split('_')[1].isdigit() else 0):
    if not d.split('_')[1].isdigit(): continue
    L=os.path.join(d,"logs","log.run")
    if not os.path.exists(L): continue
    txt=open(L,errors='ignore').read().split('\n')
    # 200k used a shorter relaxation window than the other arms
    nr = 20 if d=="mesh_200000" else 50
    it=[]; step=0; steps=0
    for l in txt:
        if l.startswith('Time = '): step+=1; steps+=1; continue
        m=re.search(r'^ *Linear solve (?:converged|did not converge) due to [A-Z_]+ iterations (\d+)',l)
        if m and step>nr: it.append(int(m.group(1)))
    snes=sum(1 for l in txt if 'SNES Function norm' in l)
    rows.append((int(d.split('_')[1]), steps, it, snes, nr))
print(f"{'cells':>9} {'steps':>6} {'relaxSkip':>10} {'KSPmed':>7} {'mean':>7} {'max':>5} {'n':>6}")
for c,steps,it,snes,nr in rows:
    if it: print(f"{c:>9} {steps:>6} {nr:>10} {s.median(it):>7.0f} {s.mean(it):>7.1f} {max(it):>5} {len(it):>6}")
    else:  print(f"{c:>9} {steps:>6} {nr:>10} {'-':>7} {'-':>7} {'-':>5} {0:>6}")
print()
good=[r for r in rows if r[2]]
if len(good)>1:
    import math
    print("  cells ratio -> KSP-iteration ratio (mesh-independent PC would give ~1.0):")
    a=good[0]
    for r in good[1:]:
        cr=r[0]/a[0]; kr=s.median(r[2])/s.median(a[2])
        expo=math.log(kr)/math.log(cr) if cr>1 else float('nan')
        print(f"    x{cr:>6.1f} cells  ->  x{kr:.2f} iterations   (exponent {expo:+.2f})")
    print()
    print("  exponent ~0.00 : mesh-independent, the PC scales")
    print("  exponent ~0.25 : ~n^1/4, weak growth (2-D sqrt of a 1-D-ish scaling)")
    print("  exponent ~0.50 : ~sqrt(n), the signature of an UNpreconditioned elliptic operator")
PY
