#!/usr/bin/env python3
"""Dump ALL centreline points of every snapshot to profiles/<time>.csv.

Written 2026-09-07 for the frozen-spike diagnosis of grubert2009_spike.

The mesh is 2000 cells = 400 along x x 5 in the extruded direction; the cell
pairs at equal x are averaged, giving 400 centreline points.

reducedE is stored in the snapshots as a `uniform 0` PLACEHOLDER (the solver
uses it as a lookup key internally and never writes real values), so E/N is
reduced HERE from the E vector field, which is populated.  N is the neutral
density at the case's p = 100 Pa, T = 300 K.
"""
import re, os, sys, glob
N_GAS = 100.0/(1.380649e-23*300.0)   # 2.4144e22 m^-3
GAP   = 0.01                          # m

def scalars(d, name, ncells=None):
    p = os.path.join(d, name)
    if not os.path.isfile(p): return None
    s = open(p).read()
    m = re.search(r"internalField\s+nonuniform List<scalar>\s*\n(\d+)\s*\n\(\s*\n(.*?)\n\)", s, re.S)
    if m:
        v = [float(x) for x in m.group(2).split()]
        assert len(v) == int(m.group(1)), f"{name}: {m.group(1)} declared, {len(v)} parsed"
        return v
    m = re.search(r"internalField\s+uniform\s+([0-9.eE+-]+)", s)
    return [float(m.group(1))]*(ncells or 1) if m else None

def vecmag(d, name, ncells):
    s = open(os.path.join(d, name)).read()
    m = re.search(r"internalField\s+nonuniform List<vector>\s*\n(\d+)\s*\n\(\s*\n(.*?)\n\)", s, re.S)
    if not m:
        m0 = re.search(r"internalField\s+uniform\s+\(([^)]*)\)", s)
        if not m0: raise SystemExit(f"cannot parse {name}")
        c = [float(x) for x in m0.group(1).split()]
        return [sum(q*q for q in c)**0.5]*ncells
    out = []
    for line in m.group(2).splitlines():
        c = [float(x) for x in line.strip().strip("()").split()]
        out.append(sum(q*q for q in c)**0.5)
    assert len(out) == ncells, f"{name}: {len(out)} vs {ncells}"
    return out

FIELDS = ("n_e","n_Arp","meanE","chargeDensity","ePotential")

def one(d, outdir):
    cx = scalars(d, "Cx")
    if cx is None: return None
    n = len(cx)
    f = {k: scalars(d, k, n) for k in FIELDS}
    f["Emag"] = vecmag(d, "E", n)
    order = sorted(range(n), key=lambda i: cx[i])
    xs = []; col = {k: [] for k in f}
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and abs(cx[order[j]]-cx[order[i]]) < 1e-12: j += 1
        idx = order[i:j]; xs.append(cx[idx[0]])
        for k, v in f.items():
            col[k].append(sum(v[q] for q in idx)/len(idx) if v else float("nan"))
        i = j
    t = os.path.basename(d)
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, f"{t}.csv")
    with open(p, "w") as fh:
        fh.write(f"# t = {t} s   case = {os.path.abspath('.')}\n")
        fh.write(f"# N_gas = {N_GAS:.6e} m^-3 (100 Pa, 300 K); gap = {GAP} m\n")
        fh.write("# x is measured from the CATHODE (x/L=0) to the ANODE (x/L=1)\n")
        fh.write("x_m,x_over_L,n_e,n_Arp,ni_over_ne,Emag_Vpm,EN_Td,meanE_eV,"
                 "Te_eV,chargeDensity,ePotential_V\n")
        for k in range(len(xs)):
            ne = col["n_e"][k]; ni = col["n_Arp"][k]; U = col["meanE"][k]
            fh.write(f"{xs[k]:.8e},{xs[k]/GAP:.8f},{ne:.8e},{ni:.8e},"
                     f"{ni/max(ne,1e-30):.8e},{col['Emag'][k]:.8e},"
                     f"{col['Emag'][k]/N_GAS/1e-21:.8e},{U:.8e},{2.0*U/3.0:.8e},"
                     f"{col['chargeDensity'][k]:.8e},{col['ePotential'][k]:.8e}\n")
    return p, len(xs)

if __name__ == "__main__":
    outdir = "profiles"
    dirs = sorted([d for d in glob.glob("[0-9]*e-0*") if os.path.isdir(d)], key=float)
    for d in dirs:
        r = one(d, outdir)
        print(f"  {d:<22} {'-> '+r[0]+f'  ({r[1]} pts)' if r else 'SKIPPED (no Cx)'}")
