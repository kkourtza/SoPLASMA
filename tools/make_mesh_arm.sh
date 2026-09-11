#!/usr/bin/env bash
# Build one mesh-convergence arm from a chosen-floor parent.
#   usage: make_mesh_arm.sh <parent_case> <target_case> <NX> <BUMP>
# Regenerates the mesh and the 0/ fields ONLY. The Boltzmann and ion tables are
# mesh-INDEPENDENT and are reused from the parent, so genMechTables/ionmob are
# deliberately NOT re-run (they cost minutes and would be identical).
# NOT `set -e`: OpenFOAM's etc/bashrc returns non-zero, and with -e the script
# then dies SILENTLY with no log and no message. Measured 2026-09-07 -- the
# first test produced an empty logs/ directory and no output at all. Each step
# is guarded by die() instead, which says which one failed.
PARENT=$1; TGT=$2; NX=$3; BUMP=$4
[ -d "$PARENT" ] || { echo "no parent $PARENT"; exit 1; }
rm -rf "$TGT"; cp -r "$PARENT" "$TGT"
cd "$TGT"
rm -rf postProcessing logs COMPARE.md; mkdir -p logs
for d in $(ls -d [0-9]*e-0* 2>/dev/null); do rm -rf "$d"; done
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1 || true
command -v gmshToFoam >/dev/null || { echo "FAILED: OpenFOAM env not loaded"; exit 1; }
# SEPARATE lines, deliberately. In `export A=x B=$A` bash expands $A BEFORE
# assigning it, so a one-liner set SoPLASMA_ETC=/etc and plasmaSetupBoundaries
# died on `Cannot read the boundary-role library: "/etc/boundaryRoles"`.
# Measured 2026-09-07. Allrun-setup uses separate lines for this reason.
export SoPLASMA=$HOME/soplasma-scratch
export SoPLASMA_ETC=$SoPLASMA/etc
export SoPLASMA_SRC=$SoPLASMA/src
[ -f "$SoPLASMA_ETC/boundaryRoles" ] || { echo "FAILED: SoPLASMA_ETC wrong: $SoPLASMA_ETC"; exit 1; }
sed -i -e "s/^NX = [0-9]*;/NX = $NX;/" -e "s/Using Bump [0-9.]*;/Using Bump $BUMP;/" gap1cm.geo
die(){ echo "FAILED: $1"; exit 1; }
gmsh -2 -format msh2 gap1cm.geo -o gap1cm.msh > logs/log.gmsh 2>&1 || die gmsh
~/ct-env/bin/python "$SoPLASMA"/tools/msh2Dto3D.py gap1cm.msh -o gap1cm_3D.msh -t 2e-4 \
    > logs/log.msh2Dto3D 2>&1 || die msh2Dto3D
rm -rf constant/polyMesh
gmshToFoam gap1cm_3D.msh > logs/log.gmshToFoam 2>&1 || die gmshToFoam
~/ct-env/bin/python "$SoPLASMA"/tools/msh2Dto3D.py --fix-boundary constant/polyMesh/boundary \
    >> logs/log.msh2Dto3D 2>&1 || die fix-boundary
~/ct-env/bin/python - <<'PY'
import re
p="constant/polyMesh/boundary"; s=open(p).read()
for n in ("side_lo","side_hi"):
    s=re.sub(r"(\b%s\s*\n\s*\{\s*\n\s*type\s+)patch;" % n, r"\1symmetryPlane;", s)
open(p,"w").write(s)
PY
checkMesh > logs/log.checkMesh 2>&1 || die checkMesh
grep -q "Mesh OK" logs/log.checkMesh || die "checkMesh did not report Mesh OK"
rm -rf 0 && mkdir -p 0
plasmaSetupBoundaries     > logs/log.plasmaSetupBoundaries 2>&1 || die plasmaSetupBoundaries
plasmaCreateSpeciesFields > logs/log.plasmaCreateSpeciesFields 2>&1 || die plasmaCreateSpeciesFields
echo "ARM-OK $TGT  NX=$NX Bump=$BUMP  cells=$(grep -m1 'cells:' logs/log.checkMesh | awk '{print $2}')  fields=$(ls 0 | wc -l)"
