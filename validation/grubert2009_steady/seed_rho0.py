#!/usr/bin/env python3
"""Seed the STEADY-route initial GUESS: quasi-neutral bulk, digitised fall.

WHY THIS CONSTRUCTION. Seeding n_e and n_Ar+ straight from Grubert's digitised
figure 3 fails (2026-09-06): rho = e(n_Arp - n_e) DIFFERENCES two curves that
agree to a few percent in the negative glow, so a few percent of digitisation
error becomes ~100% of the space charge -- and rho sets the field. Measured
consequence: bulk field 18x too high, dt 2000x smaller than the unseeded run.

THE FIX -- use the constraint that is ACCURATE IN EACH REGION:
  * in the FALL (x/L < XLO) the densities differ by DECADES, so their
    difference is well conditioned: use the digitised curves as they are.
  * in the BULK (x/L > XHI) quasi-neutrality holds to (lambda_D/L)^2, so
    n_e = n_Arp is MORE right than differencing two digitised curves.
  * blend smoothly between, so rho carries no artificial discontinuity.

IT IS A GUESS, NOT A SEED. A steady solve needs only a state in the BASIN OF
ATTRACTION, not a self-consistent one -- see doc/steady-and-stability-design.md
section 2d. That is what makes the digitisation usable here after it failed as
a transient seed.

TRANSITION AT x/L ~ 0.40, NOT 0.25. Measured 2026-09-07: 0.25 leaves a 662 Td
bulk field, where 1/nu_i = 5.81 ns against a ~6 ns electron transit -- the bulk
would AVALANCHE. The digitised curves say quasi-neutrality begins at 0.40
(n_e/n_Arp = 0.919 there, 0.306 at 0.30), and the fall-charge field swing
SATURATES beyond it (13,098 / 13,277 / 13,326 Td at XT 0.35/0.40/0.45), so the
transition is fixed by the data rather than chosen.

Requires Cx: run `postProcess -func writeCellCentres -time 0` first.
"""
import re, os, csv, sys

REF   = os.path.expanduser("~/Projects/SoEEDF/validation/dias2025/reference")
GAP   = 0.01
SCALE = 1e15                     # 1e9 cm^-3 -> m^-3
XLO, XHI = 0.35, 0.45            # blend window; quasi-neutral above XHI
VGAP  = -500.0                   # Grubert's gap voltage
QE, EPS0 = 1.602176634e-19, 8.8541878128e-12
N     = 100.0/(1.380649e-23*300.0)

def read_ref(name):
    rows = [l for l in open(os.path.join(REF,name)) if not l.startswith("#")]
    return sorted((float(r["z_over_d"]),float(r["value"]))
                  for r in csv.DictReader(rows))

def read_tab(path):
    s=open(path).read()
    p=[(float(a),float(b)) for a,b in
       re.findall(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)",s)]
    return sorted(p)

def interp(p,x):
    if x<=p[0][0]: return p[0][1]
    if x>=p[-1][0]: return p[-1][1]
    lo,hi=0,len(p)-1
    while hi-lo>1:
        m=(lo+hi)//2
        if p[m][0]<=x: lo=m
        else: hi=m
    (x0,y0),(x1,y1)=p[lo],p[hi]
    return y0 if x1==x0 else y0+(y1-y0)*(x-x0)/(x1-x0)

def read_field(path):
    s=open(path).read()
    m=re.search(r"internalField\s+nonuniform List<scalar>\s*\n(\d+)\s*\n\(\s*\n(.*?)\n\)",s,re.S)
    if m: return [float(x) for x in m.group(2).split()], s
    if re.search(r"internalField\s+uniform\s+([0-9.eE+-]+)\s*;",s): return None, s
    raise SystemExit("cannot parse %s"%path)

def write_field(path, vals, src):
    # COUNT THE SUBSTITUTIONS, do not compare strings. Re-seeding with the same
    # values leaves the text IDENTICAL, and an `out == src` guard then reports
    # "failed to substitute" for a write that in fact succeeded -- a guard that
    # convicts a no-op. Measured 2026-09-07.
    body="\n".join("%.8e"%v for v in vals)
    new="internalField   nonuniform List<scalar> \n%d\n(\n%s\n)\n;"%(len(vals),body)
    out,k=re.subn(r"internalField\s+uniform\s+[0-9.eE+-]+\s*;",new,src,count=1)
    if k==0:
        out,k=re.subn(r"internalField\s+nonuniform List<scalar>\s*\n\d+\s*\n\(\s*\n.*?\n\)\s*;",
                      new,src,count=1,flags=re.S)
    if k==0: raise SystemExit("failed to substitute internalField in %s"%path)
    open(path,"w").write(out)

def smoothstep(t):
    t=max(0.0,min(1.0,t)); return t*t*(3.0-2.0*t)

def main():
    d = sys.argv[1] if len(sys.argv)>1 else "0"
    cx,_ = read_field(os.path.join(d,"Cx"))
    if cx is None:
        raise SystemExit("Cx is uniform -- run postProcess -func writeCellCentres -time 0")
    ncell=len(cx)

    sp=open("constant/plasmaSpeciesProperties").read()
    m=re.search(r"minNumberDensity\s+([0-9.eE+-]+)",sp)
    FLOOR=float(m.group(1)) if m else 1e11

    ne_ref  = read_ref("grubert_fig3_na_LMEA_electrons.csv")
    nAr_ref = read_ref("grubert_fig3_na_LMEA_Arp.csv")

    def dens(xl):
        a=max(interp(nAr_ref,xl)*SCALE, FLOOR)
        e=max(interp(ne_ref, xl)*SCALE, FLOOR)
        w=smoothstep((xl-XLO)/(XHI-XLO))
        return a, (1.0-w)*e + w*a          # quasi-neutral as w -> 1

    # ---- 1-D axis for the Poisson prediction (the mesh repeats x in y,z) ----
    xu=sorted(set(round(v,12) for v in cx))
    nAr_u=[];ne_u=[]
    for x in xu:
        a,e=dens(x/GAP); nAr_u.append(a); ne_u.append(e)
    rho_u=[QE*(a-e) for a,e in zip(nAr_u,ne_u)]

    def cumtrap(f,xs):
        o=[0.0]
        for i in range(1,len(xs)):
            o.append(o[-1]+0.5*(f[i]+f[i-1])*(xs[i]-xs[i-1]))
        return o
    Irho=cumtrap(rho_u,xu); Esh=[v/EPS0 for v in Irho]
    IIrho=cumtrap(Esh,xu)
    E0=(VGAP-IIrho[-1])/GAP
    Eu=[E0+v for v in Esh]
    Vu=[VGAP-c for c in cumtrap(Eu,xu)]
    td=lambda e: abs(e)/N*1e21

    # ---- nEps_e from the LFA equilibrium at THAT field ----
    mE=read_tab("constant/plasmaTables/meanEnergy_vs_reducedE")
    kiz=read_tab("constant/plasmaTables/k_EI_AR_ION_AR_vs_reducedE")
    # Axis is E/N in V m^2, NOT V/m -- read as V/m every lookup silently
    # returns the table MAXIMUM (2644 eV), identical for 12 Td and 13000 Td.
    # CLAMPED to the energy model's own admissible range so the initial state
    # is already inside it: localEnergyEnergyModel defaults meanEnergyMax 100
    # eV and meanEnergyMin 1.5*kB*Tgas. Unclamped the LFA equilibrium at the
    # 13206 Td cathode field is 363 eV, which the solver would clamp on step
    # one anyway -- but those cells sit at the n_e floor (1e11) and so carry
    # no energy content; seeding inside the range keeps t=0 admissible rather
    # than relying on the clamp to repair it.
    UMIN, UMAX = 1.5*8.617333262e-5*300.0, 100.0
    Umap={x:min(max(interp(mE,abs(e)/N),UMIN),UMAX) for x,e in zip(xu,Eu)}

    # ---- write the three fields on the real mesh ----
    out={}
    for fname in ("n_Arp","n_e","nEps_e"):
        _,src=read_field(os.path.join(d,fname))
        vals=[]
        for x in cx:
            a,e=dens(x/GAP)
            if   fname=="n_Arp":  vals.append(a)
            elif fname=="n_e":    vals.append(e)
            else:                 vals.append(e*Umap[round(x,12)])
        write_field(os.path.join(d,fname),vals,src)
        out[fname]=vals
        pk=max(range(ncell),key=lambda k: vals[k])
        print("  %-7s min %.4g  max %.4g  peak at x/L = %.4f"
              %(fname,min(vals),max(vals),cx[pk]/GAP))

    # ---- THE NO-RUN CHECK, printed so a bad seed is refused before launch ----
    ihi=min(range(len(xu)),key=lambda k: abs(xu[k]/GAP-XHI))
    bulk=[td(Eu[k]) for k in range(ihi,len(xu))]
    bmin,bmax=min(bulk),max(bulk)
    nu=lambda t: interp(kiz,t*1e-21)*N
    tr=6e-3/1e6                              # electron transit across the bulk
    print("\n  --- no-run check ---")
    print("  E at cathode      %9.1f Td" % td(Eu[0]))
    print("  fall carries      %9.1f V of %.0f  (%.1f%%)"
          %(Vu[ihi]-Vu[0],abs(VGAP),100*abs(Vu[ihi]-Vu[0])/abs(VGAP)))
    print("  bulk E/N          %9.2f .. %.2f Td" % (bmin,bmax))
    worst=max(bmin,bmax)
    t_iz=1.0/nu(worst) if nu(worst)>0 else float("inf")
    print("  1/nu_i at %6.1f Td  %9.3g s   vs %.3g s transit  -> %.4g x margin"
          %(worst,t_iz,tr,t_iz/tr))
    ok = t_iz > 100*tr
    print("\n  %s: bulk cannot avalanche (need 1/nu_i >> transit)"
          % ("PASS" if ok else "FAIL"))
    if not ok:
        raise SystemExit("SEED REFUSED: bulk would avalanche")
    print("  SEED-OK")

if __name__=="__main__":
    main()
