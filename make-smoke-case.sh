#!/usr/bin/env bash
# Build a SMALL copy of the streamer case, for checking that a code change runs
# and conserves what it should.
#
# The full case is 1.1M cells (130x130 refined five times) and takes ~90 s per
# check. Nearly all of that is mesh, and a dictionary change or a source-term
# refactor does not care how fine the mesh is. This drops the refinement passes
# and runs two timesteps: ~68x fewer cells, seconds instead of minutes.
#
# What it is NOT for: anything where the ANSWER matters -- streamer velocity,
# propagation, grid convergence. Those need the real mesh. This checks that the
# code runs, that charge balances, and that fields stay finite.
set -e
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1 || true
export SoPLASMA=$HOME/soplasma-scratch SoPLASMA_SRC=$HOME/soplasma-scratch/src SoPLASMA_ETC=$HOME/soplasma-scratch/etc
SRC=$HOME/soplasma-scratch/tutorials/plasma/soPlasmaFoam/positiveStreamer/positiveStreamer_fixedMesh
CASE=${CASE:-$HOME/streamer-smoke}

rm -rf "$CASE"; cp -r "$SRC" "$CASE"; cd "$CASE"
rm -rf .snapshot0 processor* [0-9]*e-* logs log.* 2>/dev/null || true
sed -i 's|^python initGaussianSeed.py|~/ct-env/bin/python initGaussianSeed.py|' Allrun-serial

# Allclean does `rm -r 0`, and Allrun then does `cp -r 0.orig/* 0/`, which
# fails when 0/ no longer exists -- so ePotential and surfCharge never arrive
# and, worse, never get processed by changeDictionary. Copying them in
# afterwards is NOT equivalent: the raw 0.orig files carry a `.*` catch-all
# patch entry that changeDictionary exists to specialise, so a wedge patch ends
# up with a fixedValue field. Create the directory instead.
sed -i 's|^cp -r 0.orig/\* 0/|mkdir -p 0 \&\& cp -r 0.orig/* 0/|' Allrun-serial

# Coarse base mesh, and no refinement passes at all.
sed -i "s|    (130 130 1)|    (${NCELL:-40} ${NCELL:-40} 1)|" system/blockMeshDict
if [ "${NREFINE:-0}" -lt 5 ]; then
  for lvl in $(seq $(( ${NREFINE:-0} + 1 )) 5); do
    sed -i "/topoSetDict$lvl/d;/refineMeshDict$lvl/d" Allrun-serial
  done
fi

# Two timesteps is enough to exercise every path: sources, the charge
# diagnostic, transport, and the second-step reuse of the generated tables.
sed -i "s|^endTime .*|endTime                             2e-12;|" configuration/config
sed -i "s|^writeInterval .*|writeInterval                       2e-12;|" configuration/config

./Allrun-serial > /tmp/smoke-setup.log 2>&1 || true

rm -rf .snapshot0 && cp -r 0 .snapshot0
echo "smoke case at $CASE  ($(grep -m1 '^[0-9]*$' constant/polyMesh/owner) faces)"
echo "run it with:  CASE=$CASE ./rerun.sh"
