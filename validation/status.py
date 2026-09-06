#!/usr/bin/env python3
"""Robust status. Reading a LIVE log can catch a half-written line, so only
lines matching a COMPLETE float to end-of-line are accepted."""
import re, os, csv, sys
NUM = re.compile(r"^Time = ([0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)\s*$")
A = 2.0e-7
for d in sys.argv[1:]:
    L=f"{d}/logs/log.soPlasmaFoam"
    if not os.path.exists(L): print(f"  {d:16s} no log"); continue
    ts=[]; lim="-"; rej=0; crash=0
    for line in open(L, errors="replace"):
        m=NUM.match(line)
        if m: ts.append(float(m.group(1)))
        elif "deltaT set by:" in line: lim=line.split("by:",1)[1].strip()
        elif "DISCARDING this step" in line: rej+=1
        elif "stack trace" in line: crash+=1
    dt=(ts[-1]-ts[-2]) if len(ts)>1 else 0
    j=v="-"
    c=f"{d}/postProcessing/dischargeCurrent/current.csv"
    if os.path.exists(c):
        rows=[r for r in csv.reader(open(c)) if r and not r[0].startswith('#')]
        if len(rows)>1:
            h=rows[0]; last=[float(x) for x in rows[-1]]
            j=f"{abs(last[h.index('I_total')])/A*1e3/1e4:.4f}"
            v=f"{last[h.index('V_applied')]:.1f}"
    print(f"  {d:16s} t={ts[-1] if ts else 0:.4e} dt={dt:.2e} steps={len(ts):<7d} rej={rej:<5d} "
          f"crash={crash} lim={lim:<20s} V={v:>7s} j={j:>10s}")
