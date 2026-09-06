#!/usr/bin/env python3
"""Build an initial state satisfying BOTH current continuity AND Poisson.

WHY. Seeding n_e and n_Ar+ straight from Grubert's digitised figure 3 failed
(2026-09-06): the space charge is their DIFFERENCE, and in the negative glow the
two curves agree to a few percent, so a few percent of digitisation error in each
becomes ~100% error in rho -- and rho sets the field. Result: bulk field 18x too
high, current 9x too high, dt 2000x smaller than the unseeded run.

THE FIX. Use the constraint that SUMS, and DERIVE the one that differences:

    current continuity (exact in steady state):
        n_e mu_e + n_i mu_i = J/(e E)          <- a SUM, robust to digitisation
    Poisson:
        n_i - n_e = (eps0/e) dE/dx             <- the difference, now DERIVED

Two equations, two unknowns at each point:

    n_e = (P - D mu_i)/(mu_e + mu_i),   n_i = n_e + D
    P = J/(e E),  D = (eps0/e) dE/dx

Iterate: densities -> E (from continuity) -> dE/dx -> D -> densities.

The digitised profiles still set the SHAPE through P; only the small difference
is replaced by something physically derived rather than digitised.

Caveat kept in view: this is a DRIFT-ONLY current balance. In the negative glow
E is small and the density gradient steep, so diffusion is comparable to drift
there; the cathode fall, which dominates the voltage integral, is drift-dominated
and unaffected.
"""
import re, os, csv, sys

REF = os.path.expanduser("~/Projects/SoEEDF/validation/dias2025/reference")
T, TI = "constant/plasmaTables", "constant/ionTables"
QE, EPS0 = 1.602176634e-19, 8.8541878128e-12
N   = 100/(1.380649e-23*300)
L, JTOT = 0.01, 5.11          # Grubert: 0.511 mA/cm^2

def tab(d,n):
    s=open(os.path.join(d,n)).read()
    return [(float(a),float(b)) for a,b in re.findall(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)",s)]
def interp(p,x):
    lo=[q for q in p if q[0]<=x]; hi=[q for q in p if q[0]>=x]
    if not lo: return p[0][1]
    if not hi: return p[-1][1]
    (x0,y0),(x1,y1)=lo[-1],hi[0]
    return y0 if x1==x0 else y0+(y1-y0)*(x-x0)/(x1-x0)
def ref(n):
    rows=[l for l in open(os.path.join(REF,n)) if not l.startswith('#')]
    return sorted((float(r["z_over_d"]),float(r["value"])) for r in csv.DictReader(rows))

def main():
    muNe, muNi = tab(T,"muN_vs_reducedE"), tab(TI,"muN_Arp_vs_reducedE")
    NEr, NIr = ref("grubert_fig3_na_LMEA_electrons.csv"), ref("grubert_fig3_na_LMEA_Arp.csv")
    M = 400
    xs = [(k+0.5)*L/M for k in range(M)]
    ne = [max(interp(NEr,x/L),1e-6)*1e15 for x in xs]
    ni = [max(interp(NIr,x/L),1e-6)*1e15 for x in xs]
    FLOOR = 1e11
    E = [1e4]*M

    for it in range(400):
        # E from current continuity, given the densities
        for k in range(M):
            e_ = E[k]
            for _ in range(50):
                en = e_/N
                me, mi = interp(muNe,en)/N, interp(muNi,en)/N
                new = JTOT/(QE*(ne[k]*me + ni[k]*mi))
                if abs(new-e_) < 1e-8*max(e_,1): e_=new; break
                e_ = 0.5*e_ + 0.5*new
            E[k] = e_
        # rho from Poisson, then rebuild the densities keeping the SUM
        dx = L/M
        dEdx = [ (E[min(k+1,M-1)]-E[max(k-1,0)])/((min(k+1,M-1)-max(k-1,0))*dx) for k in range(M) ]
        maxrel = 0.0
        for k in range(M):
            en = E[k]/N
            me, mi = interp(muNe,en)/N, interp(muNi,en)/N
            P = JTOT/(QE*E[k])
            D = -(EPS0/QE)*dEdx[k]        # E decreases with x in the fall -> rho>0
            nen = (P - D*mi)/(me + mi)
            nin = nen + D
            nen, nin = max(nen,FLOOR), max(nin,FLOOR)
            rel = abs(nen-ne[k])/max(ne[k],FLOOR)
            maxrel = max(maxrel, rel)
            ne[k], ni[k] = 0.5*ne[k]+0.5*nen, 0.5*ni[k]+0.5*nin
        if maxrel < 1e-6:
            print(f"  converged in {it+1} iterations"); break
    else:
        print(f"  stopped after 400 iterations, worst relative change {maxrel:.3g}")

    V = sum(E)*(L/M)
    print(f"\n  integrated gap voltage = {V:.1f} V   (Grubert 500 V, ratio {V/500:.3f})")
    pk = max(E); kp = E.index(pk)
    print(f"  |E| peak {pk:.4g} V/m at x/L = {xs[kp]/L:.4f}")
    print(f"\n  {'x/L':>6} {'n_e':>11} {'n_i':>11} {'n_i/n_e':>9} {'E [V/m]':>11} {'E/N [Td]':>9}")
    for z in (0.005,0.05,0.15,0.25,0.40,0.55,0.74,0.90,0.99):
        k=min(int(z*M),M-1)
        print(f"  {xs[k]/L:>6.3f} {ne[k]:>11.4g} {ni[k]:>11.4g} {ni[k]/ne[k]:>9.3f} "
              f"{E[k]:>11.4g} {E[k]/N/1e-21:>9.1f}")
    # how far did the densities move from the digitised ones?
    d0e=[max(interp(NEr,x/L),1e-6)*1e15 for x in xs]
    d0i=[max(interp(NIr,x/L),1e-6)*1e15 for x in xs]
    re_=[abs(ne[k]-d0e[k])/max(d0e[k],FLOOR) for k in range(M)]
    ri_=[abs(ni[k]-d0i[k])/max(d0i[k],FLOOR) for k in range(M)]
    print(f"\n  departure from the digitised profiles: n_e median {sorted(re_)[M//2]*100:.1f}%, "
          f"n_i median {sorted(ri_)[M//2]*100:.1f}%")
    import json
    json.dump({"xs":xs,"ne":ne,"ni":ni,"E":E}, open("consistent_seed.json","w"))
    print("  -> consistent_seed.json")

if __name__ == "__main__":
    main()
