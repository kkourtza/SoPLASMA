#!/usr/bin/env python3
"""Report the evolution of peak electron density and peak electron temperature.

  usage: news.py <case_dir> [n_rows]

Te is not logged directly; the solver logs the LMEA MEAN ENERGY range in eV, and
Te = (2/3)*meanE for a Maxwellian, which is the same relation the extraction
scripts use (verified: meanE_eV / Te_eV = 1.5000 in the spike profiles).
"""
import re, sys, os, math

case = sys.argv[1].rstrip("/")
NROWS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
log = os.path.join(case, "logs/log.soPlasmaFoam")

# ramp rate [V/us] inferred from the sourceVoltage table, so V_src is reported
vrate = None
bf = os.path.join(case, "configuration/boundaries")
if os.path.exists(bf):
    m = re.search(r"sourceVoltage\s+table\s*\(\(0 0\)\s*\(([0-9.eE+-]+)\s+(-?[0-9.eE+-]+)\)\)",
                  open(bf).read())
    if m:
        tramp, vmax = float(m.group(1)), float(m.group(2))
        vrate = (vmax / (tramp * 1e6), tramp * 1e6)   # V per us, ramp end in us

t = None; ne = None; rows = []
try:
    for line in open(log, errors="ignore"):
        if line.startswith("Time = "):
            try: t = float(line.split("=", 1)[1])
            except ValueError: t = None
        elif "e [m^-3]:" in line:
            m = re.search(r"min = ([0-9.eE+-]+)\s+max = ([0-9.eE+-]+)", line)
            if m:
                try: ne = (float(m.group(1)), float(m.group(2)))
                except ValueError: ne = None
        elif "LMEA mean energy" in line and t is not None and ne is not None:
            m = re.search(r"\[\s*([0-9.eE+-]+)\s*,\s*([0-9.eE+-]+)\s*\]", line)
            if m:
                try: rows.append((t, ne[0], ne[1], float(m.group(1)), float(m.group(2))))
                except ValueError: pass
except FileNotFoundError:
    print("no log for %s" % case); sys.exit(0)

if len(rows) < 2:
    print("%s: only %d records yet" % (case, len(rows))); sys.exit(0)

print("EVOLUTION -- %s" % case)
print("  t        simulated time [us]")
if vrate: print("  V_src    source voltage [V] = %.0f * t[us], holding at %.0f V after t = %.1f us"
                % (vrate[0], vrate[0]*vrate[1], vrate[1]))
print("  n_e,max  peak electron density [m^-3]   (floor / n_e,min = %.3g)" % rows[-1][1])
print("  Te,max   peak electron temperature [eV] = (2/3) * max mean energy")
print("  Te,min   MINIMUM electron temperature [eV] -- the bulk value; should stay low")
print("  rate     d(ln n_e,max)/dt [1/s] over the preceding row interval")
print()
print("     t[us]   V_src[V]    n_e,max[m^-3]   Te,max[eV]  Te,min[eV]    rate[1/s]")
step = max(1, len(rows)//NROWS)
sel = rows[::step]
if sel[-1] != rows[-1]: sel.append(rows[-1])
prev = None
for r in sel:
    tt, nmin, nmax, emin, emax = r
    v = ("%9.1f" % (vrate[0]*min(tt*1e6, vrate[1]))) if vrate else "     --  "
    rt = "     --   "
    if prev and tt > prev[0] and nmax > 0 and prev[2] > 0:
        rt = "%10.3g" % (math.log(nmax/prev[2])/(tt-prev[0]))
    print("   %7.3f  %s  %14.4g  %11.2f %11.2f  %s"
          % (tt*1e6, v, nmax, (2.0/3.0)*emax, (2.0/3.0)*emin, rt))
    prev = r
