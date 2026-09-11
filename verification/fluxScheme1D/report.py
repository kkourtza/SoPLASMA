#!/usr/bin/env python3
"""Fit convergence orders from results.txt and print the table."""
import re, math, collections

rows = []
for ln in open('results.txt'):
    m = re.match(r'RESULT scheme=(\S+)\s+N=\s*(\d+)\s+h=(\S+)\s+Pe=(\S+)\s+'
                 r'PeGrid=(\S+)\s+L2=(\S+)\s+Linf=(\S+)\s+L2rel=(\S+)', ln)
    if m:
        rows.append(dict(sch=m.group(1), N=int(m.group(2)), h=float(m.group(3)),
                         Pe=float(m.group(4)), Pg=float(m.group(5)),
                         L2=float(m.group(6)), Linf=float(m.group(7))))

by = collections.defaultdict(list)
for r in rows:
    by[(round(r['Pe'], 6), r['sch'])].append(r)

print("\nFLUX-SCHEME VERIFICATION -- steady 1D drift-diffusion WITH a source")
print("  d/dx(v n - D dn/dx) = S,  n(0)=1, n(1)=2, S=1, D=1, L=1")
print("  Exact: n = C1 + C2 exp(Pe x) + (S/v) x")
print("  Pe      = v L / D          (domain Peclet)")
print("  PeGrid  = v h / D          (grid Peclet, h = 1/N)")
print("  L2      = RMS error against the exact solution, at cell centres")
print("  order   = log2(L2(h)/L2(h/2)) between successive rows; 2 = second order\n")

for Pe in sorted({k[0] for k in by}):
    print(f"=== domain Pe = {Pe:g} " + "="*54)
    print(f"{'scheme':>14} {'N':>5} {'PeGrid':>10} {'L2':>13} {'Linf':>13} {'order':>7}")
    for sch in ['std:upwind','std:limitedLinear','std:Minmod','std:vanLeer','std:MUSCL','std:SuperBee','std:linear','std:ROUNDF','std:ROUNDA','std:ROUNDAplus','std:ROUNDL','ScharfetterGummel']:
        rs = sorted(by.get((Pe, sch), []), key=lambda r: r['N'])
        prev = None
        for r in rs:
            o = ''
            if prev and r['L2'] > 0 and prev['L2'] > 0:
                o = f"{math.log(prev['L2']/r['L2'])/math.log(2.0):7.2f}"
            print(f"{sch:>14} {r['N']:5d} {r['Pg']:10.4f} {r['L2']:13.5e} "
                  f"{r['Linf']:13.5e} {o:>7}")
            prev = r
        if rs: print()
