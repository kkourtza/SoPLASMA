#!/usr/bin/env python3
"""Plot the SPATIAL distribution of n_e and Te from an OpenFOAM snapshot.

  usage: plot_ne_Te.py <case> <timedir> [--cx <case/timedir with Cx>] [-o out.png]

Writes a CSV of the extracted profile next to the figure, and the figure itself,
so the data and the plot are both kept (rule 6). Te = (2/3)*meanE.

Cell centres: OpenFOAM does not write Cx unless asked. Any case whose
constant/polyMesh/{points,owner} are BYTE-IDENTICAL has identical cell ordering
and centres, so Cx may be borrowed from such a case. The caller must pass
--cx and this script VERIFIES the md5 match before using it.
"""
import re, sys, os, csv, hashlib, argparse

def md5(p):
    return hashlib.md5(open(p,"rb").read()).hexdigest()

def read_field(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform List<scalar>\s*\n(\d+)\s*\n\(\s*\n(.*?)\n\)", s, re.S)
    if m: return [float(x) for x in m.group(2).split()]
    m = re.search(r"internalField\s+uniform\s+([0-9.eE+-]+)\s*;", s)
    if m: return None            # uniform -> caller decides
    raise SystemExit("cannot parse %s" % path)

ap = argparse.ArgumentParser()
ap.add_argument("case"); ap.add_argument("timedir")
ap.add_argument("--cx", required=True, help="case/timedir holding Cx")
ap.add_argument("-o", "--out", default=None)
a = ap.parse_args()

cxcase = a.cx.split("/")[0]
for f in ("points", "owner"):
    p1 = os.path.join(a.case, "constant/polyMesh", f)
    p2 = os.path.join(cxcase, "constant/polyMesh", f)
    if md5(p1) != md5(p2):
        raise SystemExit("MESH MISMATCH on %s: cannot borrow Cx from %s" % (f, cxcase))
print("  mesh md5 verified identical to %s -- Cx is transferable" % cxcase)

cx = read_field(os.path.join(a.cx, "Cx"))
ne = read_field(os.path.join(a.case, a.timedir, "n_e"))
nA = read_field(os.path.join(a.case, a.timedir, "n_Arp"))
mE = read_field(os.path.join(a.case, a.timedir, "meanE"))
phi = read_field(os.path.join(a.case, a.timedir, "ePotential"))
if ne is None or mE is None:
    raise SystemExit("n_e or meanE is uniform in this snapshot -- nothing to plot")

GAP = 0.01
# average the 5 lateral cells at equal x -> centreline
from collections import defaultdict
acc = defaultdict(lambda: [0.0,0.0,0.0,0.0,0])
for i, x in enumerate(cx):
    k = round(x, 12)
    s = acc[k]
    s[0] += ne[i]; s[1] += nA[i]; s[2] += mE[i]
    s[3] += phi[i] if phi else 0.0; s[4] += 1
xs = sorted(acc)
prof = [(x, acc[x][0]/acc[x][4], acc[x][1]/acc[x][4],
         acc[x][2]/acc[x][4], acc[x][3]/acc[x][4]) for x in xs]

base = a.out or os.path.join(a.case, "profile_%s" % a.timedir)
csvp = base.replace(".png","") + ".csv"
with open(csvp, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["# case=%s time=%s" % (a.case, a.timedir)])
    w.writerow(["x_m","x_over_L","n_e_m3","n_Arp_m3","meanE_eV","Te_eV","ePotential_V"])
    for x,n,ni,e,pp in prof:
        w.writerow(["%.6e"%x, "%.6f"%(x/GAP), "%.6e"%n, "%.6e"%ni,
                    "%.4f"%e, "%.4f"%(2*e/3), "%.4f"%pp])
print("  wrote %s (%d centreline points)" % (csvp, len(prof)))

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("  matplotlib unavailable -- CSV written, no figure"); sys.exit(0)

xl = [p[0]/GAP for p in prof]
fig, ax = plt.subplots(3, 1, figsize=(7.5, 9), sharex=True)
ax[0].semilogy(xl, [p[1] for p in prof], label="$n_e$")
ax[0].semilogy(xl, [p[2] for p in prof], "--", label="$n_{Ar^+}$")
ax[0].set_ylabel("density  [m$^{-3}$]"); ax[0].legend(); ax[0].grid(alpha=.3)
ax[1].plot(xl, [2*p[3]/3 for p in prof], color="crimson")
ax[1].set_ylabel("$T_e$  [eV]"); ax[1].grid(alpha=.3)
ax[2].plot(xl, [p[4] for p in prof], color="k")
ax[2].set_ylabel("potential  [V]"); ax[2].set_xlabel("$x/L$   (0 = cathode, 1 = anode)")
ax[2].grid(alpha=.3)
ax[0].set_title("%s   t = %s s" % (a.case, a.timedir))
fig.tight_layout()
png = base.replace(".csv","") + ".png" if not a.out else a.out
fig.savefig(png, dpi=130)
print("  wrote %s" % png)
