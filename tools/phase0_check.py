#!/usr/bin/env python3
"""Phase 0 no-run check: does the rho~0 construction give a 12-30 Td bulk
field with the fall carrying most of the 400-500 V?"""
import os, csv

REF = os.path.expanduser("~/Projects/SoEEDF/validation/dias2025/reference")
QE, EPS0 = 1.602176634e-19, 8.8541878128e-12
N   = 100.0/(1.380649e-23*300.0)      # 2.4144e22 m^-3
L   = 0.01
SCALE = 1e15
VCATH = -500.0
XT  = 0.25                            # transition x/L
FLOOR = 1e-6                          # in SCALE units

def read_ref(name):
    rows = [l for l in open(os.path.join(REF,name)) if not l.startswith("#")]
    return sorted((float(r["z_over_d"]), float(r["value"]))
                  for r in csv.DictReader(rows))

def interp(pts, x):
    if x <= pts[0][0]:  return pts[0][1]
    if x >= pts[-1][0]: return pts[-1][1]
    lo, hi = 0, len(pts)-1
    while hi-lo > 1:
        m=(lo+hi)//2
        if pts[m][0] <= x: lo=m
        else: hi=m
    (x0,y0),(x1,y1)=pts[lo],pts[hi]
    return y0 if x1==x0 else y0+(y1-y0)*(x-x0)/(x1-x0)

ne_ref  = read_ref("grubert_fig3_na_LMEA_electrons.csv")
nAr_ref = read_ref("grubert_fig3_na_LMEA_Arp.csv")

print("=== what the digitisation gives AT the transition x/L = %.2f ===" % XT)
for x in (0.20,0.22,0.25,0.28,0.30):
    e,a = interp(ne_ref,x), interp(nAr_ref,x)
    print("  x/L=%.2f  n_e=%.4f  n_Arp=%.4f  (n_Arp-n_e)/n_Arp = %+.3f"
          % (x,e,a,(a-e)/a if a else 0))

# ---- build on a fine uniform grid ----
M = 20001
xs = [L*i/(M-1) for i in range(M)]
nAr, ne = [], []
for x in xs:
    xl = x/L
    a = max(interp(nAr_ref, xl), FLOOR)
    if xl < XT:
        e = max(interp(ne_ref, xl), FLOOR)
    else:
        e = a                          # quasi-neutral by construction
    nAr.append(a*SCALE); ne.append(e*SCALE)

rho = [QE*(a-e) for a,e in zip(nAr,ne)]

# ---- 1D Poisson: dE/dx = rho/eps0, V(0)=VCATH, V(L)=0 ----
dx = L/(M-1)
def cum(f):
    out=[0.0]
    for i in range(1,len(f)): out.append(out[-1]+0.5*(f[i]+f[i-1])*dx)
    return out
Irho = cum(rho)                                   # INT_0^x rho dx'
Eshape = [v/EPS0 for v in Irho]                   # E(x) - E0
IIrho = cum(Eshape)                               # INT_0^x (E-E0) dx'
E0 = (VCATH - IIrho[-1])/L                        # so that INT_0^L E dx = V(0)
E  = [E0+v for v in Eshape]
V  = [VCATH - c for c in cum(E)]

def td(e): return abs(e)/N*1e21

it = int(XT*(M-1))
print("\n=== the check ===")
print("  V(0)=%.1f V   V(L)=%.3f V  (target 0)" % (V[0],V[-1]))
print("  E0 (cathode) = %.4g V/m  =  %.1f Td" % (E0,td(E0)))
print("  voltage across fall  x/L<%.2f : %8.1f V" % (XT, V[it]-V[0]))
print("  voltage across bulk  x/L>%.2f : %8.1f V" % (XT, V[-1]-V[it]))
print("  fall carries %.1f%% of the gap voltage" % (100*abs(V[it]-V[0])/500))
bulk=[td(E[i]) for i in range(it,M)]
print("  bulk E/N: min %.2f  max %.2f  mean %.2f Td   [TARGET 12-30]"
      % (min(bulk),max(bulk),sum(bulk)/len(bulk)))
print("\n  profile:")
for xl in (0.0,0.02,0.05,0.10,0.15,0.20,0.249,0.26,0.40,0.60,0.80,0.99):
    i=min(int(xl*(M-1)),M-1)
    print("   x/L=%.3f  n_e=%9.3e  n_Arp=%9.3e  E/N=%9.2f Td  V=%8.1f V"
          % (xl,ne[i],nAr[i],td(E[i]),V[i]))

# ---------------------------------------------------------------
# The bulk field came out 662 Td, not 12-30. WHY: with rho == 0 the
# bulk field is CONSTANT, so any voltage not carried by the fall is
# forced into it. How much voltage CAN the fall carry?
print("\n=== diagnosis: sweep the imposed gap voltage ===")
print("  Vgap      E0[Td]   bulk[Td]   fall V    fall %")
for Vc in (-500,-450,-420,-400,-394,-390,-380,-370,-350):
    E0s = (Vc - IIrho[-1])/L
    Es  = [E0s+v for v in Eshape]
    Vs  = [Vc - c for c in cum(Es)]
    print("  %6.0f  %9.1f  %9.2f  %8.1f  %6.1f%%"
          % (Vc, td(E0s), td(Es[it]), Vs[it]-Vs[0],
             100*abs(Vs[it]-Vs[0])/abs(Vc)))

# what voltage puts the bulk field exactly at 20 Td?
# E_bulk = (Vc - IIrho[-1])/L + Eshape[it]  ->  solve for Vc
Etarget = -20.0*N/1e21          # negative: field points -x
Vstar = (Etarget - Eshape[it])*L + IIrho[-1]
E0s = (Vstar - IIrho[-1])/L
Es  = [E0s+v for v in Eshape]
Vs  = [Vstar - c for c in cum(Es)]
print("\n  => bulk EXACTLY 20 Td requires Vgap = %.1f V" % Vstar)
print("     then E0 = %.1f Td, fall carries %.1f V (%.1f%%)"
      % (td(E0s), Vs[it]-Vs[0], 100*abs(Vs[it]-Vs[0])/abs(Vstar)))
print("     COMPARE: independent current-continuity no-run result = 393.9 V")
