#!/usr/bin/env python3
"""Compare two cases at MATCHED SOURCE VOLTAGE, in n_e and Te.

  usage: compare_at_voltage.py <caseA> <caseB> [nrows]

Matched VOLTAGE, not matched time: cases with different ramp rates sit at
different voltages at any instant, so a matched-time comparison would measure
the ramp rather than the discharge. The ramp rate is read from each case's own
sourceVoltage table.  Te = (2/3)*meanE.
"""
import re, sys, os, math

def ramp(case):
    m = re.search(r"sourceVoltage\s+table\s*\(\(0 0\)\s*\(([0-9.eE+-]+)\s+(-?[0-9.eE+-]+)\)\)",
                  open(os.path.join(case, "configuration/boundaries")).read())
    if not m: raise SystemExit("no sourceVoltage ramp in %s" % case)
    tramp, vmax = float(m.group(1)), float(m.group(2))
    return abs(vmax)/(tramp*1e6), tramp*1e6, abs(vmax)   # |V| per us, ramp end us, |Vmax|

def trace(case):
    t=None; ne=None; out=[]
    for line in open(os.path.join(case,"logs/log.soPlasmaFoam"), errors="ignore"):
        if line.startswith("Time = "):
            try: t=float(line.split("=",1)[1])
            except ValueError: t=None
        elif "e [m^-3]:" in line:
            m=re.search(r"min = ([0-9.eE+-]+)\s+max = ([0-9.eE+-]+)",line)
            if m:
                try: ne=(float(m.group(1)),float(m.group(2)))
                except ValueError: ne=None
        elif "LMEA mean energy" in line and t is not None and ne is not None:
            m=re.search(r"\[\s*([0-9.eE+-]+)\s*,\s*([0-9.eE+-]+)\s*\]",line)
            if m:
                try: out.append((t,ne[1],float(m.group(1)),float(m.group(2))))
                except ValueError: pass
    return out

A,B = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/")
N = int(sys.argv[3]) if len(sys.argv)>3 else 9
rA,rB = ramp(A), ramp(B)
tA,tB = trace(A), trace(B)
tof = lambda r,V: (min(V,r[2])/r[0])*1e-6
at  = lambda tr,tt: min(tr,key=lambda z:abs(z[0]-tt))
vmaxA = min(rA[0]*tA[-1][0]*1e6, rA[2]); vmaxB = min(rB[0]*tB[-1][0]*1e6, rB[2])
VMAX = min(vmaxA, vmaxB)

print("COMPARISON AT MATCHED SOURCE VOLTAGE")
print("  %-24s ramp %.0f V/us, reaches %.0f V at t = %.1f us" % (A, rA[0], rA[2], rA[1]))
print("  %-24s ramp %.0f V/us, reaches %.0f V at t = %.1f us" % (B, rB[0], rB[2], rB[1]))
print()
print("  |V_src|  = source voltage magnitude [V] -- the comparison axis")
print("  t        = simulated time [us] at which each case is at that voltage")
print("  n_e,max  = peak electron density [m^-3]")
print("  Te,max   = peak electron temperature [eV] = (2/3)*max meanE")
print("  Te,min   = MINIMUM Te [eV], i.e. the bulk value")
print("  ratio    = %s / %s at the same voltage" % (os.path.basename(A), os.path.basename(B)))
print()
print("  |V|[V] | %-9s t   n_e,max     Te,max Te,min | %-9s t   n_e,max     Te,max Te,min | n_e ratio Te,max ratio"
      % (os.path.basename(A)[:9], os.path.basename(B)[:9]))
print("  -------|" + "-"*48 + "|" + "-"*48 + "|" + "-"*22)
step = VMAX/N
V = step
while V <= VMAX + 1e-9:
    ca = at(tA, tof(rA,V)); cb = at(tB, tof(rB,V))
    print("  %6.1f | %8.3f %11.4g %7.2f %6.2f | %8.3f %11.4g %7.2f %6.2f | %8.3f %10.3f"
          % (V, ca[0]*1e6, ca[1], 2*ca[3]/3, 2*ca[2]/3,
                cb[0]*1e6, cb[1], 2*cb[3]/3, 2*cb[2]/3,
                ca[1]/cb[1], (2*ca[3]/3)/(2*cb[3]/3)))
    V += step
print("\n  compared up to |V_src| = %.1f V (limited by %s)"
      % (VMAX, A if vmaxA<vmaxB else B))
