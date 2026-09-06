#!/usr/bin/env python3
"""Seed n_e and n_Ar+ from Grubert et al. (2009) figure 3 onto this mesh.

WHY: the transient solver cannot carry this discharge THROUGH ignition -- at
ignition the plasma conductance g is small (the plasma is forming) exactly when
the growth rate gamma is large, so the circuit response time C_gap/g is 3-4
decades too slow, for ANY two-terminal circuit. Measured 2026-09-06 with a
current-imposed and a voltage-imposed arm failing identically.

That is a statement about the PATH. This script tests the DESTINATION: start AT
Grubert's profile and ask whether the transient solver holds it. If it does,
Grubert's solution is a stable fixed point of our model -- which is the
validation claim we actually want, and it needs no new solver, because
transient relaxation from a good initial state IS pseudo-timestepping.

It is NOT a production route. A solver that needs the answer as input is not a
solver, and the user rejected seeding for that purpose (2026-09-06).

Reference units, from the digitisation headers: number density in 1e9 cm^-3,
i.e. 1e15 m^-3. `z_over_d` is measured FROM THE CATHODE, and this mesh has the
cathode at x = 0 (gap1cm.geo), so z/d = x/L directly.

Requires Cx: run `postProcess -func writeCellCentres -time 0` first.
"""
import re, os, csv, sys

REF = os.path.expanduser(
    "~/Projects/SoEEDF/validation/dias2025/reference")
GAP = 0.01
SCALE = 1e15                    # 1e9 cm^-3 -> m^-3

def read_ref(name):
    rows = [l for l in open(os.path.join(REF, name)) if not l.startswith("#")]
    d = list(csv.DictReader(rows))
    pts = sorted((float(r["z_over_d"]), float(r["value"])) for r in d)
    return pts

def interp(pts, x):
    if x <= pts[0][0]:  return pts[0][1]
    if x >= pts[-1][0]: return pts[-1][1]
    lo, hi = 0, len(pts)-1
    while hi - lo > 1:
        m = (lo+hi)//2
        if pts[m][0] <= x: lo = m
        else: hi = m
    (x0,y0),(x1,y1) = pts[lo],pts[hi]
    return y0 if x1==x0 else y0 + (y1-y0)*(x-x0)/(x1-x0)

def read_field(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform List<scalar>\s*\n(\d+)\s*\n\(\s*\n(.*?)\n\)", s, re.S)
    if m: return [float(x) for x in m.group(2).split()], s
    m = re.search(r"internalField\s+uniform\s+([0-9.eE+-]+)\s*;", s)
    if m: return None, s
    raise SystemExit(f"cannot parse {path}")

def write_field(path, vals, src):
    body = "\n".join(f"{v:.8e}" for v in vals)
    new = f"internalField   nonuniform List<scalar> \n{len(vals)}\n(\n{body}\n)\n;"
    out = re.sub(r"internalField\s+uniform\s+[0-9.eE+-]+\s*;", new, src, count=1)
    if out == src:
        out = re.sub(r"internalField\s+nonuniform List<scalar>\s*\n\d+\s*\n\(\s*\n.*?\n\)\s*;",
                     new, src, count=1, flags=re.S)
    if out == src:
        raise SystemExit(f"failed to substitute internalField in {path}")
    open(path, "w").write(out)

def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "0"
    cx, _ = read_field(os.path.join(d, "Cx"))
    if cx is None:
        raise SystemExit("Cx is uniform -- run postProcess -func writeCellCentres -time 0")
    n = len(cx)

    # the floors this case declares, so the seed never drops below them
    sp = open("constant/plasmaSpeciesProperties").read()
    def floor(name, default):
        m = re.search(rf"{name}\s*\{{[^}}]*minNumberDensity\s+([0-9.eE+-]+)", sp)
        if m: return float(m.group(1))
        m = re.search(r"minNumberDensity\s+([0-9.eE+-]+)", sp)
        return float(m.group(1)) if m else default
    f_e  = floor("electron", 1e11)
    f_i  = floor("Arp", 1e11)

    for fname, ref, fl in (("n_e",   "grubert_fig3_na_LMEA_electrons.csv", f_e),
                           ("n_Arp", "grubert_fig3_na_LMEA_Arp.csv",       f_i)):
        pts = read_ref(ref)
        _, src = read_field(os.path.join(d, fname))
        vals = []
        for x in cx:
            v = interp(pts, x/GAP)*SCALE
            vals.append(max(v, fl))       # digitisation noise goes slightly negative
        write_field(os.path.join(d, fname), vals, src)
        print(f"  {fname:<7} seeded from {ref}")
        print(f"          min {min(vals):.4g}  max {max(vals):.4g} m^-3  "
              f"(floor {fl:.3g}, peak at x/L = {cx[max(range(n), key=lambda k: vals[k])]/GAP:.4f})")

    # report the net charge the seed implies -- it is the difference of two
    # INDEPENDENTLY digitised curves, so it is noisier than either
    ne, _ = read_field(os.path.join(d, "n_e"))
    ni, _ = read_field(os.path.join(d, "n_Arp"))
    QE = 1.602176634e-19
    rho = [QE*(ni[k]-ne[k]) for k in range(n)]
    print(f"\n  implied space charge: {min(rho):.4g} .. {max(rho):.4g} C/m^3")
    print(f"  net charge sign: {sum(1 for r in rho if r>0)} of {n} cells positive")
    print("  (the ion and electron curves were digitised separately, so their")
    print("   DIFFERENCE is noisier than either -- expect a few tau_diel of settling)")

if __name__ == "__main__":
    main()
