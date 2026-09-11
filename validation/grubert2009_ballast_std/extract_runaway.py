#!/usr/bin/env python3
"""Extract the runaway's time series from a soPlasmaFoam log.

WHY: the diverged Grubert ballast run showed n_Arp 38x the reference, a sheath
26x too thin, and E/N = 31645 Td (13x the expected cathode fall). Those are a
CLOSED feedback loop, so the end state cannot say what STARTED it. This pulls
the joint history so the ordering can be read: does the FIELD run away first
(field-driven) or the DENSITY (density-driven)?  Written 2026-09-06.

Emits runaway.csv: time, n_e_max, n_Arp_max, Emag_max, EN_Td, sheath_mm.
Plot it with plot_runaway.py, which is deliberately separate so the figure can
be changed without re-parsing a 1.1 GB log.
"""
import re, sys, csv

LOG = sys.argv[1] if len(sys.argv) > 1 else "logs/log.soPlasmaFoam"
OUT = sys.argv[2] if len(sys.argv) > 2 else "runaway.csv"

# 100 Pa, 300 K
N_GAS = 100.0 / (1.380649e-23 * 300.0)
EPS0, QE = 8.8541878128e-12, 1.602176634e-19

t_re  = re.compile(rb"^Time = ([0-9.eE+-]+)")
rng_re = re.compile(rb"^\s+(\S+) \[[^\]]+\]:\s+min = (\S+)\s+max = (\S+)")

rows, cur = [], {}
with open(LOG, "rb") as f:
    for line in f:
        m = t_re.match(line)
        if m:
            if cur.get("t") is not None and "e" in cur and "Emag" in cur:
                rows.append(cur)
            cur = {"t": float(m.group(1))}
            continue
        m = rng_re.match(line)
        if m and cur:
            name = m.group(1).decode()
            if name in ("e", "Arp", "Emag"):
                try:
                    cur[name] = float(m.group(3))
                except ValueError:
                    pass          # a torn line mid-write; skip the field
if cur.get("t") is not None and "e" in cur and "Emag" in cur:
    rows.append(cur)

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["time_s", "n_e_max", "n_Arp_max", "Emag_max",
                "EN_Td", "sheath_mm"])
    for r in rows:
        E  = r.get("Emag", 0.0)
        ni = r.get("Arp", 0.0)
        en = E / N_GAS / 1e-21
        # E = e n_i d / eps0  =>  d = eps0 E /(e n_i): the sheath thickness the
        # measured field and ion density imply, for comparison with the ~4.4 mm
        # a normal cathode fall at 100 Pa should have.
        d  = EPS0 * E / (QE * ni) * 1e3 if ni > 0 else float("nan")
        w.writerow([f"{r['t']:.9e}", f"{r.get('e',0):.6e}", f"{ni:.6e}",
                    f"{E:.6e}", f"{en:.4f}", f"{d:.6f}"])

print(f"{len(rows)} samples -> {OUT}")
if rows:
    print(f"  t: {rows[0]['t']:.4e} .. {rows[-1]['t']:.4e} s")
