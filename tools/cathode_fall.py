#!/usr/bin/env python3
"""Extract the NORMAL/ABNORMAL GLOW SIMILARITY PARAMETERS from a snapshot.

  d_c    cathode-fall thickness -- distance from the cathode to where E_x
         first crosses zero (the cathode fall ends where the field reverses)
  phi_c  cathode-fall voltage   -- potential drop across that distance
  j      conduction current density at the cathode

and reports the similarity groups p*d_c and j/p^2, which is what Grubert 2009
Fig. 5 plots against phi_c and compares with Engel & Steenbeck.

CALIBRATION, from Grubert 2009 sec. III B: Engel & Steenbeck expect
1-10 mA/cm^2 for this discharge; Grubert's own LMEA reaches ~0.5 mA/cm^2, and
they state ALL hydrodynamic descriptions fall below the expected value. So
0.5 mA/cm^2 is the MODEL target, not the physical one -- do not read a factor
of a few against Engel & Steenbeck as a failure of this code alone.

Usage: cathode_fall.py <case> <time> [pressure_Pa]
"""
import re, sys, os

def rd(p, vec=False):
    s = open(p, errors='ignore').read()
    m = re.search(r'internalField\s+nonuniform\s+List<(?:scalar|vector)>\s*\d+\s*'
                  r'\((.*?)\n\)\s*;', s, re.S)
    if not m:
        u = re.search(r'internalField\s+uniform\s+\(?([-\d.eE+ ]+)\)?\s*;', s)
        return None
    b = m.group(1)
    if vec:
        return [tuple(map(float, t.split())) for t in re.findall(r'\(([^)]*)\)', b)]
    return [float(x) for x in b.split()]

case = sys.argv[1]
t    = sys.argv[2]
pPa  = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0

d = os.path.join(case, t)
cx = rd(os.path.join(d, 'Cx'))
if cx is None:
    sys.exit(f"no Cx in {d} -- run postProcess -func writeCellCentres")

seen = {}
for i, x in enumerate(cx):
    k = round(x, 12)
    if k not in seen: seen[k] = i
cols = [seen[k] for k in sorted(seen)]          # ascending x: cathode -> anode

E   = rd(os.path.join(d, 'E'), vec=True)
phi = rd(os.path.join(d, 'ePotential'))
ne  = rd(os.path.join(d, 'n_e'))
mue = rd(os.path.join(d, 'mu_e'))

# cathode is at x = 0 (phi most negative there)
xs = [cx[i] for i in cols]

# d_c : first zero crossing of E_x moving away from the cathode
dc = None
for k in range(1, len(cols)):
    a, b = E[cols[k-1]][0], E[cols[k]][0]
    if a*b < 0:
        # linear interpolation of the crossing
        f = abs(a)/(abs(a)+abs(b))
        dc = xs[k-1] + f*(xs[k]-xs[k-1])
        kc = k
        break

qe = 1.602176634e-19
jc = qe*ne[cols[0]]*mue[cols[0]]*abs(E[cols[0]][0])   # A/m^2 at the cathode

print(f"case {case}  t = {float(t)*1e6:.4f} us   p = {pPa:g} Pa")
print(f"  phi(cathode)      = {phi[cols[0]]:.3f} V")
if dc is None:
    print("  E_x NEVER CROSSES ZERO -- no cathode fall has formed yet.")
    print("  (a developed glow must have a field reversal at the end of the fall)")
else:
    # phi_c = drop from the cathode to the end of the fall
    phic = abs(phi[cols[kc]] - phi[cols[0]])
    L = max(xs) - min(xs)
    if dc > 0.5*L:
        print(f"  *** THE ZERO CROSSING IS AT {dc/L*100:.0f}% OF THE GAP -- this is the")
        print("      ANODE-SIDE reversal, NOT a cathode fall. A cathode fall at")
        print("      100 Pa cm is a fraction of a mm. The discharge is still")
        print("      Townsend-like, with a near-uniform field across the gap.")
        print("      Treat d_c/phi_c below as MEANINGLESS for similarity purposes.")
    print(f"  d_c               = {dc*1e3:.4f} mm")
    print(f"  phi_c             = {phic:.2f} V")
    print(f"  p*d_c             = {pPa*dc*100:.3f} Pa cm")
    print(f"  j (cathode)       = {jc*0.1:.4f} mA/cm^2   ({jc:.4e} A/m^2)")
    print(f"  j/p^2             = {jc/(pPa**2):.4e} A/(m^2 Pa^2)")
    print()
    print("  Grubert 2009 Fig. 5 plots p*d_c and j/p^2 vs phi_c over")
    print("  pd = 50..1000 Pa cm; this case is p*d = 100 Pa cm.")
    print("  Their LMEA reaches j ~ 0.5 mA/cm^2; Engel & Steenbeck expect 1-10.")
