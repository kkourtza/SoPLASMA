#!/usr/bin/env bash
# MESH-INDEPENDENCE SWEEP (project rule 43).
#
# THE QUESTION: do outer Krylov iterations per Newton step stay constant as the
# mesh refines? That, not wall-clock at one size, is what decides whether a
# preconditioner scales. Everything measured on 2026-09-10 was at 2000 cells
# in serial, so none of it is yet evidence about the 1e6-cell cases this solver
# exists for.
#
# METHOD, and each choice removes a confound:
#   * dt SCALED WITH THE CELL SIZE (dt ~ dx, i.e. dt/NX = const), with
#     adjustTimeStep false so the controller cannot interfere.
#
#     NOT a single fixed dt on every mesh. That was the first design and it was
#     WRONG, for the reason the user pointed out: refining tightens the Courant
#     limit, so at fixed dt a finer mesh runs at a HIGHER Courant number and its
#     linear system is harder for a reason that has nothing to do with how the
#     preconditioner scales with cell count. Measured 2026-09-10, and the
#     confound was not small -- it was the whole effect:
#         2000  cells (NX=400):  Co_conv(e) =  6.50
#         20240 cells (NX=1265): Co_conv(e) = 20.62   <- ratio 3.17
#     and 3.17 is exactly the 3.16x refinement. Both are already far above 1.
#     Scaling dt with dx holds Co constant, so any remaining iteration growth
#     is the preconditioner's own mesh dependence.
#   * from a DEVELOPED state at t = 1.954497e-06, mapped onto each mesh with
#     mapFields. NOT from t = 0: at a cold start every species sits exactly on
#     its floor and the field is 0.1 V/m, so the solver's start-up guard keeps
#     PICARD engaged and Newton never runs at all. That mistake cost one whole
#     vacuous sweep on 2026-09-10 -- 20 steps, 0 SNES solves, 26 Picard
#     correctors. Reaching a developed discharge from cold takes ~2 us of
#     physics, i.e. ~2e6 steps at this dt, which is not a preconditioner test.
#   * the mapped state keeps its ORIGINAL TIME, so the voltage ramp and every
#     time-dependent control see the same instant on every mesh.
#   * a fixed, small step count -- this measures the PRECONDITIONER, not the
#     physics, so it does not need to reach ignition.
#   * both mesh directions refined, keeping the cell aspect ratio comparable, so
#     what changes is resolution and not cell shape.
#
# Tables (constant/plasmaTables, constant/ionTables) are MESH-INDEPENDENT and
# are copied from the base case rather than regenerated -- genMechTables is slow
# and would add nothing.
cd "${0%/*}" || exit 1

BASE=/home/kkourtza/soplasma-scratch/validation/grubert2009_pseudo
# The developed state to map onto every mesh, and its time.
SRC=/home/kkourtza/soplasma-scratch/validation/grubert_fix
SRCT=1.954497e-06
SOP=$HOME/soplasma-scratch
# NOTE: no `set -u` anywhere in this script. OpenFOAM's bashrc
# references unset variables, so `set -u` kills the script on the
# source line -- silently, because that source is redirected. That
# cost one apparently-launched sweep that never ran (0-byte log).
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1
source "$SOP"/etc/bashrc >/dev/null 2>&1
export SoPLASMA=$SOP SoPLASMA_SRC=$SOP/src SoPLASMA_ETC=$SOP/etc

NSTEPS=${NSTEPS:-20}      # steps MEASURED
NRELAX=${NRELAX:-50}      # steps run first and DISCARDED
DT=${DT:-1e-12}

# WHY A RELAXATION WINDOW. mapFields RESAMPLES the developed state onto each
# mesh, so a refined mesh does not start from a state converged on its OWN
# grid -- the first steps carry an interpolation transient, and it is largest
# on the finest mesh, i.e. biasing exactly the arm the sweep is about. The run
# therefore takes NRELAX steps that are thrown away and only the following
# NSTEPS are measured. Single solver invocation; the split is done in the
# analysis, so there is no second restart to go wrong.

# NX NY  -> cells = NX*NY. 400x5=2000 is the reference case.
# DT is scaled as DTREF*400/NX so the Courant number is the SAME on
# every arm -- see the note above.
GRIDS=("400 5" "1265 16" "4000 50")
DTREF=${DT:-1e-12}
NXREF=400

echo "############ MESH-INDEPENDENCE SWEEP ############"
echo "  dt SCALED as $DTREF*400/NX (constant Courant), $NRELAX relax + $NSTEPS measured"
echo "  developed state mapped from $SRC at $SRCT"
echo

for g in "${GRIDS[@]}"; do
  set -- $g; NX=$1; NY=$2
  CELLS=$((NX*NY))
  # dt ~ dx so Co is constant across arms
  DT=$(python3 -c "print(f'{$DTREF*$NXREF/$NX:.6e}')")
  D=mesh_${CELLS}
  echo "=== ${CELLS} cells (NX=$NX NY=$NY, dt=$DT) -> $D ==="

  rm -rf "$D"; mkdir -p "$D/logs"
  # only what a case needs; NOT the base's run output
  cp -a "$BASE"/configuration "$BASE"/system "$BASE"/etc "$D"/ 2>/dev/null
  cp -a "$BASE"/gap1cm.geo "$D"/
  mkdir -p "$D/constant"
  cp -a "$BASE"/constant/argon_plasma.foam "$BASE"/constant/argon_plasma.mech.json \
        "$BASE"/constant/plasmaSpeciesProperties "$BASE"/constant/plasmaTransportProperties \
        "$BASE"/constant/plasmaTables "$BASE"/constant/ionTables "$D"/constant/ 2>/dev/null

  ( cd "$D" || exit 1
    sed -i "s/^NX = .*/NX = $NX;/; s/^NY = .*/NY = $NY;/" gap1cm.geo
    # fixed dt, and stop after NSTEPS
    sed -i "s/^adjustTimeStep .*/adjustTimeStep                      false;/" configuration/config
    sed -i "s/^deltaT  *.*/deltaT                              $DT;/" configuration/config
    sed -i "s/^endTime  *.*/endTime                             $(python3 -c "print(f'{$NSTEPS*$DT:.6e}')");/" configuration/config

    gmsh -2 -format msh2 gap1cm.geo -o gap1cm.msh > logs/log.gmsh 2>&1 || { echo "  gmsh FAILED"; exit 1; }
    ~/ct-env/bin/python "$SoPLASMA"/tools/msh2Dto3D.py gap1cm.msh -o gap1cm_3D.msh -t 2e-4 \
        > logs/log.msh2Dto3D 2>&1 || { echo "  msh2Dto3D FAILED"; exit 1; }
    cp system/fvSolution-foam system/fvSolution 2>/dev/null
    gmshToFoam gap1cm_3D.msh > logs/log.gmshToFoam 2>&1 || { echo "  gmshToFoam FAILED"; exit 1; }
    ~/ct-env/bin/python "$SoPLASMA"/tools/msh2Dto3D.py --fix-boundary constant/polyMesh/boundary \
        >> logs/log.msh2Dto3D 2>&1 || { echo "  fix-boundary FAILED"; exit 1; }
    ~/ct-env/bin/python - <<'PY'
import re
p="constant/polyMesh/boundary"; s=open(p).read()
for n in ("side_lo","side_hi"):
    s=re.sub(r"(\b%s\s*\n\s*\{\s*\n\s*type\s+)patch;" % n, r"\1symmetryPlane;", s)
open(p,"w").write(s)
PY
    checkMesh > logs/log.checkMesh 2>&1
    grep -q "Mesh OK" logs/log.checkMesh || echo "  WARNING: checkMesh not OK"
    grep -oE "cells: *[0-9]+" logs/log.checkMesh | head -1 | sed 's/^/    actual /'

    rm -rf 0 && mkdir -p 0
    plasmaSetupBoundaries     > logs/log.setupBoundaries 2>&1 || { echo "  setupBoundaries FAILED"; exit 1; }
    plasmaCreateSpeciesFields > logs/log.createFields 2>&1 || { echo "  createFields FAILED"; exit 1; }

    # Put the generated fields AT THE SOURCE TIME, then overwrite them with the
    # mapped developed state. Keeping the time (rather than mapping to 0) means
    # the voltage ramp and every time-dependent control see the same instant on
    # every mesh -- and it is what lets the start-up guard release, since the
    # mapped state has no species on its floor.
    mv 0 "$SRCT"
    sed -i "s/^startFrom .*/startFrom       startTime;/" system/controlDict
    sed -i "s/^startTime .*/startTime       $SRCT;/" system/controlDict
    # 7 digits: OpenFOAM auto-bumps precision when WRITING a time name and then
    # cannot find its own directory on restart at the default 6 (rule 41).
    sed -i "s/^timePrecision .*/timePrecision   7;/" system/controlDict
    sed -i "s/^endTime  *.*/endTime                             $(python3 -c "print(f'{float(\"$SRCT\")+($NRELAX+$NSTEPS)*$DT:.9e}')");/" configuration/config

    mapFields "$SRC" -sourceTime "$SRCT" -consistent > logs/log.mapFields 2>&1 \
        || { echo "  mapFields FAILED"; tail -4 logs/log.mapFields; exit 1; }

    /usr/bin/time -f "    wall %e s  maxRSS %M kB" \
      soPlasmaFoam > logs/log.run 2>logs/log.time || echo "  run exited nonzero"
    cat logs/log.time | tail -2
    ns=$(grep -ac "outerSolver newton (SNES)" logs/log.run)
    echo "    NEWTON ENGAGED: $ns SNES solves  (0 => VACUOUS, the run measured Picard)"
  )
  echo
done

echo "############ RESULT: outer Krylov iterations per Newton step ############"
python3 - <<'PY'
import re, glob, os, statistics as s
rows=[]
for d in sorted(glob.glob("mesh_*"), key=lambda x:int(x.split('_')[1])):
    L=os.path.join(d,"logs","log.run")
    if not os.path.exists(L): continue
    txt=open(L,errors='ignore').read().split('\n')
    # Attribute each KSP line to the step it belongs to, and DISCARD the
    # relaxation window -- the interpolation transient lives there.
    NRELAX=int(os.environ.get('NRELAX','50'))
    it=[]; step=0; steps=0
    for l in txt:
        if l.startswith('Time = '):
            step+=1; steps+=1; continue
        m=re.search(r'^ *Linear solve (?:converged|did not converge) due to [A-Z_]+ iterations (\d+)',l)
        if m and step>NRELAX: it.append(int(m.group(1)))
    snes=[int(m.group(1)) for m in
        (re.search(r'\(SNES\): reason [0-9-]+ \(positive = converged\), iterations (\d+)',l) for l in txt) if m]
    fatal=sum(1 for l in txt if 'FOAM FATAL' in l or 'trapped floating' in l)
    rows.append((int(d.split('_')[1]), steps, it, snes, fatal))
print(f"  (KSP stats EXCLUDE the first {os.environ.get('NRELAX','50')} relaxation steps)")
print(f"{'cells':>8} {'steps':>6} {'KSPits/med':>11} {'mean':>7} {'max':>6} {'SNESits/step':>13} {'fatal':>6}")
for c,steps,it,snes,fatal in rows:
    km = f"{s.median(it):.0f}" if it else "-"
    kmn= f"{s.mean(it):.1f}" if it else "-"
    kmx= f"{max(it)}" if it else "-"
    sn = f"{sum(snes)/max(steps,1):.1f}" if snes else "-"
    print(f"{c:>8} {steps:>6} {km:>11} {kmn:>7} {kmx:>6} {sn:>13} {fatal:>6}")
print()
if len(rows)>1 and all(r[2] for r in rows):
    a=s.median(rows[0][2]); b=s.median(rows[-1][2])
    f=(rows[-1][0]/rows[0][0])
    print(f"  cells x{f:.0f}  ->  median KSP its {a:.0f} -> {b:.0f}  (ratio {b/a:.2f})")
    print("  MESH-INDEPENDENT if the ratio is ~1. A ratio growing with the cell")
    print("  count means the preconditioner does NOT scale, whatever it costs at 2000.")
PY
