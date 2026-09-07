#!/usr/bin/env python3
"""Report the avalanche rate of a running case, and whether it has SETTLED.

  usage: settle_check.py <case_dir> [threshold_1_per_s]

'Settled' means the peak electron density changes little over ONE ION TRANSIT,
the slowest thing that must equilibrate. At the -200 V plateau tau_i = 7.16 us,
so a threshold of 1e4 1/s is 3% per ion transit. Derived, not chosen.

Prints ONE line, machine-readable, so a Monitor can grep it.
"""
import re, sys, math, os

case = sys.argv[1]
THRESH = float(sys.argv[2]) if len(sys.argv) > 2 else 1e4
WIN = 1.0e-6          # trailing window for the rate fit [s]

log = os.path.join(case, "logs/log.soPlasmaFoam")
t = None
pts = []
try:
    for line in open(log, errors="ignore"):
        if line.startswith("Time = "):
            try: t = float(line.split("=", 1)[1])
            except ValueError: t = None
        elif "e [m^-3]:" in line:
            m = re.search(r"min = ([0-9.eE+-]+)\s+max = ([0-9.eE+-]+)", line)
            if m and t is not None:
                try: pts.append((t, float(m.group(1)), float(m.group(2))))
                except ValueError: pass
except FileNotFoundError:
    print("%s NOLOG" % case); sys.exit(0)

if len(pts) < 5:
    print("%s TOOFEW" % case); sys.exit(0)

tn, nmin, nmax = pts[-1]
# rate over the trailing window
a = min(pts, key=lambda z: abs(z[0] - (tn - WIN)))
rate = math.log(nmax / a[2]) / (tn - a[0]) if (tn > a[0] and a[2] > 0) else float("nan")

# circuit state
vsrc = vel = icond = float("nan")
cs = os.path.join(case, "postProcessing/externalCircuit/circuit.csv")
if os.path.exists(cs):
    for l in open(cs, errors="ignore"):
        if l.startswith("#"): continue
        r = l.strip().split(",")
        if len(r) < 4 or r[0] == "time": continue
        try: _, vsrc, icond, vel = [float(x) for x in r[:4]]
        except ValueError: pass

# RUNAWAY, as measured on the case that actually failed (grubert2009_ballast_low):
#   V_el flipped SIGN to +204 V, and |I_cond| reached 95.6x I_sc = |V_src|/R.
# Either alone is disqualifying: a passive ballast cannot sustain |I_cond| > I_sc,
# and a positive electrode means the "cathode" has become an electron collector.
R_BALLAST = 1e8
isc = abs(vsrc) / R_BALLAST if vsrc == vsrc and vsrc != 0 else float("nan")
iratio = abs(icond) / isc if isc == isc and isc > 0 else float("nan")
runaway = []
if vel == vel and vsrc == vsrc and vsrc < 0 and vel > 0:
    runaway.append("V_el SIGN FLIPPED to %+.1f V" % vel)
if iratio == iratio and iratio > 1.0:
    runaway.append("|I_cond| = %.3gx I_sc" % iratio)

if runaway:
    state = "RUNAWAY"
elif rate == rate and rate < THRESH:
    state = "SETTLED"
else:
    state = "GROWING"
print("%s %s t=%.4gus n_e_max=%.4g rate=%.3g/s V_src=%.1f V_el=%+.1f I_cond=%.3g I/Isc=%.3g%s"
      % (case, state, tn * 1e6, nmax, rate, vsrc, vel, icond, iratio,
         ("  << " + "; ".join(runaway)) if runaway else ""))
