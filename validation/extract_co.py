#!/usr/bin/env python3
"""Courant-headroom sweep: dt and limiter AT A COMMON SIMULATED TIME, plus the
breakdown voltage, which is the invariant that decides the test."""
import csv, os, re, sys

CASES = [("grubert2009",15),("tuneCo_50",50),("tuneCo_100",100),
         ("tuneCo_1000",1000),("tuneCo_10000",10000)]
TCOMMON = float(sys.argv[1]) if len(sys.argv)>1 else 1.0e-6
JBD = 0.5   # mA/cm2 -- "broken down"
A = 2.0e-7

def breakdown(case):
    p=f"{case}/postProcessing/dischargeCurrent/current.csv"
    if not os.path.exists(p): return None,None
    rows=[r for r in csv.reader(open(p)) if r and not r[0].startswith('#')]
    h=rows[0]; d=[[float(x) for x in r] for r in rows[1:]]
    i=h.index("I_total"); v=h.index("V_applied")
    for r in d:
        if abs(r[i])/A*1e3/1e4 > JBD:
            return r[v], r[0]
    return None,None

def at_time(case, t):
    """dt and limiter at the step nearest t."""
    p=f"{case}/logs/log.soPlasmaFoam"
    if not os.path.exists(p): return None,None,None
    times=[]; lims=[]
    cur=None
    for line in open(p, errors="replace"):
        if line.startswith("Time = "):
            try: cur=float(line.split("=",1)[1]); times.append(cur)
            except ValueError: pass
        elif "deltaT set by:" in line and cur is not None:
            lims.append((cur, line.split("by:",1)[1].strip()))
    if len(times)<2: return None,None,None
    if times[-1] < t: return "not reached", None, times[-1]
    k=min(range(1,len(times)), key=lambda i: abs(times[i]-t))
    dt=times[k]-times[k-1]
    lim=next((L for (tt,L) in reversed(lims) if tt<=times[k]), "?")
    return dt, lim, times[-1]

print(f"  common time t = {TCOMMON:.2e} s\n")
print(f"  {'case':14s} {'Co':>6s} {'dt @ t':>12s} {'limiter @ t':22s} {'reached':>11s} {'V_breakdown':>12s}")
print("  " + "-"*84)
for c,co in CASES:
    dt,lim,last = at_time(c, TCOMMON)
    vb,tb = breakdown(c)
    dts = f"{dt:.3e}" if isinstance(dt,float) else str(dt or "-")
    print(f"  {c:14s} {co:6d} {dts:>12s} {str(lim or '-'):22s} {last if last else 0:11.3e} "
          f"{(f'{vb:.1f} V' if vb else '-'):>12s}")
print("\n  INVARIANT: V_breakdown must not move. Baseline is grubert2009 (Co=15).")
