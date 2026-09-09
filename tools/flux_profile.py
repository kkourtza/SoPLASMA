#!/usr/bin/env python3
"""Where does the conduction current reverse FIRST -- at the wall, or inside?

  usage: flux_profile.py <case> <t1> <t2> ... [--cx <case/timedir with Cx>]

For each requested time it prints the ELECTRON DRIFT FLUX along the axis,

    Gamma_e(x) = -n_e * mu_e * E_x        [m^-2 s^-1]

(electrons drift against E; the sign convention here is that Gamma_e > 0 means
electrons moving toward the ANODE, i.e. the healthy direction.)

THE TEST: if the flux reverses sign at the WALL FACE first and the reversal
then propagates inward, the boundary closure is responsible. If it reverses in
the INTERIOR first and reaches the wall later, the interior discretisation is.
"""
import re, sys, os, argparse
from collections import defaultdict

def read_field(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform List<scalar>\s*\n(\d+)\s*\n\(\s*\n(.*?)\n\)", s, re.S)
    if m: return [float(x) for x in m.group(2).split()]
    m = re.search(r"internalField\s+uniform\s+([0-9.eE+-]+)\s*;", s)
    if m: return None
    raise SystemExit("cannot parse %s" % path)

def read_vec_x(path):
    """x-component of a nonuniform vector field."""
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform List<vector>\s*\n(\d+)\s*\n\(\s*\n(.*?)\n\)", s, re.S)
    if not m: return None
    out = []
    for tok in re.findall(r"\(([^)]*)\)", m.group(2)):
        out.append(float(tok.split()[0]))
    return out

ap = argparse.ArgumentParser()
ap.add_argument("case"); ap.add_argument("times", nargs="+")
ap.add_argument("--cx", required=True)
a = ap.parse_args()

cx = read_field(os.path.join(a.cx, "Cx"))
GAP = 0.01
print("ELECTRON DRIFT FLUX along the axis.  Gamma_e = -n_e*mu_e*E_x  [m^-2 s^-1]")
print("  POSITIVE = electrons moving toward the anode (healthy).")
print("  x/L = 0 is the CATHODE, 1 is the ANODE.\n")
for td in a.times:
    d = os.path.join(a.case, td)
    if not os.path.isdir(d):
        print("  %s : no such snapshot" % td); continue
    try:
        ne = read_field(os.path.join(d, "n_e"))
        mu = read_field(os.path.join(d, "mu_e"))
        Ex = read_vec_x(os.path.join(d, "E"))
    except SystemExit as e:
        print("  %s : %s" % (td, e)); continue
    if ne is None or mu is None or Ex is None:
        print("  %s : a needed field is uniform/absent" % td); continue
    acc = defaultdict(lambda: [0.0, 0])
    for i, x in enumerate(cx):
        g = -ne[i]*mu[i]*Ex[i]
        s = acc[round(x, 12)]; s[0] += g; s[1] += 1
    xs = sorted(acc)
    prof = [(x, acc[x][0]/acc[x][1]) for x in xs]
    # sign changes
    flips = [prof[i][0]/GAP for i in range(1, len(prof))
             if (prof[i-1][1] > 0) != (prof[i][1] > 0)]
    print("  t = %-14s  Gamma_e at x/L = " % td, end="")
    for t in (0.0005, 0.01, 0.1, 0.5, 0.9, 0.99, 0.9995):
        i = min(range(len(prof)), key=lambda k: abs(prof[k][0]/GAP - t))
        print("%.2f:%+9.3g " % (t, prof[i][1]), end="")
    print("\n%s sign changes at x/L = %s" % (" "*20,
          ", ".join("%.4f" % f for f in flips) if flips else "NONE"))
