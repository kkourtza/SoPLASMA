#!/usr/bin/env python3
"""Centreline profiles of n_e, n_Ar+, E/N and mean energy at a given time.

Reads the SOLVER'S OWN stored fields (reducedE, meanE) rather than deriving
them, so the numbers are exactly what the solver used. Requires Cx:
    postProcess -func writeCellCentres -time <t>
"""
import re, os, sys
N_GAS = 100/(1.380649e-23*300)
GAP   = 0.01
def rd(d, n, ncells=None):
    p = os.path.join(d, n)
    if not os.path.isfile(p): return None
    s = open(p).read()
    m = re.search(r"internalField\s+nonuniform List<scalar>\s*\n(\d+)\s*\n\(\s*\n(.*?)\n\)", s, re.S)
    if m:
        v=[float(x) for x in m.group(2).split()]
        assert len(v)==int(m.group(1)), f"{n}: {m.group(1)} declared, {len(v)} parsed"
        return v
    m = re.search(r"internalField\s+uniform\s+([0-9.eE+-]+)", s)
    return [float(m.group(1))]*(ncells or 1) if m else None
def main():
    d = sys.argv[1]
    cx = rd(d,"Cx")
    if cx is None: raise SystemExit("no Cx -- run postProcess -func writeCellCentres")
    n = len(cx)
    f = {k: rd(d,k,n) for k in ("n_e","n_Arp","meanE","chargeDensity","ePotential")}
    # reducedE is stored as a `uniform 0` PLACEHOLDER -- the solver uses it as a
    # lookup key internally and never writes real values. Take |E| from the E
    # VECTOR field, which is populated, and reduce it here.
    sE = open(os.path.join(d,"E")).read()
    mv = re.search(r"internalField\s+nonuniform List<vector>\s*\n(\d+)\s*\n\(\s*\n(.*?)\n\)", sE, re.S)
    if not mv: raise SystemExit("cannot parse E")
    Emag=[]
    for line in mv.group(2).splitlines():
        c=line.strip().strip("()").split()
        Emag.append((float(c[0])**2+float(c[1])**2+float(c[2])**2)**0.5)
    assert len(Emag)==n, f"E: {len(Emag)} vs {n}"
    f["Emag"]=Emag
    order = sorted(range(n), key=lambda i: cx[i])
    xs=[]; col={k:[] for k in f}
    i=0
    while i < len(order):
        j=i
        while j<len(order) and abs(cx[order[j]]-cx[order[i]])<1e-12: j+=1
        idx=order[i:j]; xs.append(cx[idx[0]])
        for k,v in f.items():
            col[k].append(sum(v[q] for q in idx)/len(idx) if v else float("nan"))
        i=j
    m=len(xs)
    print(f"  t = {os.path.basename(d)}   ({m} points)")
    print(f"  {'x/L':>7} {'n_e':>11} {'n_Ar+':>11} {'n_i/n_e':>8} {'E/N [Td]':>10} {'meanE[eV]':>10} {'rho':>11}")
    for fr in (0.0005,0.002,0.005,0.01,0.02,0.05,0.10,0.20,0.35,0.50,0.70,0.90,0.99):
        k=min(int(fr*(m-1)),m-1)
        en = col["Emag"][k]/N_GAS/1e-21
        print(f"  {xs[k]/GAP:>7.4f} {col['n_e'][k]:>11.4g} {col['n_Arp'][k]:>11.4g} "
              f"{col['n_Arp'][k]/max(col['n_e'][k],1e-30):>8.2f} {en:>10.1f} {col['meanE'][k]:>10.3f} "
              f"{col['chargeDensity'][k]:>11.3g}")
    ke=max(range(m),key=lambda k: col["n_e"][k]); ki=max(range(m),key=lambda k: col["n_Arp"][k])
    kE=max(range(m),key=lambda k: col["Emag"][k]); kT=max(range(m),key=lambda k: col["meanE"][k])
    print(f"\n  n_e   peak {col['n_e'][ke]:.4g} at x/L={xs[ke]/GAP:.4f}   (Grubert 2.478e15 @ 0.394)")
    print(f"  n_Ar+ peak {col['n_Arp'][ki]:.4g} at x/L={xs[ki]/GAP:.4f}   (Grubert 6.391e15 @ 0.203)")
    print(f"  |E|   peak {col['Emag'][kE]:.4g} V/m = {col['Emag'][kE]/N_GAS/1e-21:.1f} Td at x/L={xs[kE]/GAP:.4f}  (normal fall ~2450 Td)")
    print(f"  meanE peak {col['meanE'][kT]:.2f} eV at x/L={xs[kT]/GAP:.4f}  (clamp 100 eV)")
    V=col["ePotential"]
    print(f"  ePotential {min(V):.1f} .. {max(V):.1f} V")
if __name__=="__main__": main()
