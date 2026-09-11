#!/usr/bin/env python3
"""Extract the centreline profile from a soPlasmaFoam time directory.

Emits `centreline_<time>.csv` with columns
    z_over_d, x_m, n_e, n_Arp, chargeDensity, Emag, EN_Td, meanE
in SI, plus the reduced field in Townsend, ready to compare against
Grubert et al. (2009) figure 3 (digitised in
SoEEDF/validation/dias2025/reference/grubert_fig3_*.csv, whose densities are
in 1e9 cm^-3 = 1e15 m^-3 and currents in 1e-1 mA/cm^2).

The gap runs along **x** here (gap1cm.geo: cathode at x = 0, anode at x = L),
so `z_over_d` is x/L with z measured FROM THE CATHODE, matching the reference.

Cell centres come from OpenFOAM's own `writeCellCentres`, not from parsing
constant/polyMesh/points: this mesh is a 2-D extrusion, so distinct x values
are fewer than cells and cell order is not sorted by x. Reconstructing
coordinates by hand produced garbage on 2026-09-06.

    postProcess -func writeCellCentres -time <t>
    ./extract_centreline.py <t>
"""
import re, sys, os, csv

NGAS = 100.0/(1.380649e-23*300.0)     # 100 Pa, 300 K
GAP  = 0.01                            # m

def field(d, name, ncells=None):
    path = os.path.join(d, name)
    if not os.path.isfile(path):
        return None
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform List<scalar>\s*\n(\d+)\s*\n\(\s*\n(.*?)\n\)",
                  s, re.S)
    if m:
        v = [float(x) for x in m.group(2).split()]
        if len(v) != int(m.group(1)):
            raise SystemExit(f"{name}: declared {m.group(1)}, parsed {len(v)}")
        return v
    m = re.search(r"internalField\s+uniform\s+([0-9.eE+-]+)", s)
    if m and ncells:
        return [float(m.group(1))]*ncells
    raise SystemExit(f"cannot parse {name}")

def main():
    d = sys.argv[1]
    Cx = field(d, "Cx")
    if Cx is None:
        raise SystemExit("no Cx: run `postProcess -func writeCellCentres -time <t>` first")
    n = len(Cx)
    names = ["n_e", "n_Arp", "chargeDensity", "ePotential", "meanE"]
    f = {k: field(d, k, n) for k in names}

    # collapse the extruded pair: average cells sharing an x, then sort by x
    order = sorted(range(n), key=lambda i: Cx[i])
    rows, i = [], 0
    while i < len(order):
        j = i
        while j < len(order) and abs(Cx[order[j]] - Cx[order[i]]) < 1e-12:
            j += 1
        idx = order[i:j]
        r = {"x_m": Cx[idx[0]]}
        for k in names:
            r[k] = (sum(f[k][t] for t in idx)/len(idx)) if f[k] else float("nan")
        rows.append(r)
        i = j

    # Emag from the potential, central differences on the (graded) 1-D grid
    m = len(rows)
    for k, r in enumerate(rows):
        if 0 < k < m-1:
            dx = rows[k+1]["x_m"] - rows[k-1]["x_m"]
            r["Emag"] = abs((rows[k-1]["ePotential"] - rows[k+1]["ePotential"])/dx)
        else:
            r["Emag"] = float("nan")

    out = f"centreline_{os.path.basename(d)}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["z_over_d", "x_m", "n_e", "n_Arp", "chargeDensity",
                    "Emag", "EN_Td", "meanE"])
        for r in rows:
            en = r["Emag"]/NGAS/1e-21 if r["Emag"] == r["Emag"] else float("nan")
            w.writerow([f"{r['x_m']/GAP:.6f}", f"{r['x_m']:.6e}",
                        f"{r['n_e']:.6e}", f"{r['n_Arp']:.6e}",
                        f"{r['chargeDensity']:.6e}", f"{r['Emag']:.6e}",
                        f"{en:.4f}", f"{r['meanE']:.6f}"])
    print(f"{m} points -> {out}")

    ne = max(r["n_e"] for r in rows); ni = max(r["n_Arp"] for r in rows)
    ke = max(range(m), key=lambda k: rows[k]["n_e"])
    ki = max(range(m), key=lambda k: rows[k]["n_Arp"])
    print(f"  n_e   peak {ne:.4g} m^-3 at z/d = {rows[ke]['x_m']/GAP:.4f}  (ref 2.478e15)")
    print(f"  n_Arp peak {ni:.4g} m^-3 at z/d = {rows[ki]['x_m']/GAP:.4f}  (ref 6.39e15)")
    em = [r["Emag"] for r in rows if r["Emag"] == r["Emag"]]
    if em:
        kE = max(range(1, m-1), key=lambda k: rows[k]["Emag"])
        print(f"  Emag  peak {max(em):.4g} V/m at z/d = {rows[kE]['x_m']/GAP:.4f}"
              f"  -> E/N = {max(em)/NGAS/1e-21:.0f} Td")
        # cathode fall: where |E| has fallen to 10% of its peak
        pk = max(em)
        fall = [rows[k]['x_m'] for k in range(1, m-1) if rows[k]["Emag"] > 0.1*pk]
        if fall:
            print(f"  cathode fall (|E| > 10% of peak) extends to "
                  f"{max(fall)*1e3:.3f} mm   (a normal fall at 100 Pa wants ~4.4 mm)")

if __name__ == "__main__":
    main()
