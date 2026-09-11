#!/usr/bin/env python3
"""Is the Phase 0 answer robust to the transition point XT? That was a CHOICE."""
import os, csv
REF=os.path.expanduser("~/Projects/SoEEDF/validation/dias2025/reference")
QE,EPS0=1.602176634e-19,8.8541878128e-12
N=100.0/(1.380649e-23*300.0); L=0.01; SCALE=1e15; FLOOR=1e-6
def read_ref(n):
    rows=[l for l in open(os.path.join(REF,n)) if not l.startswith("#")]
    return sorted((float(r["z_over_d"]),float(r["value"])) for r in csv.DictReader(rows))
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
ne_ref=read_ref("grubert_fig3_na_LMEA_electrons.csv")
nAr_ref=read_ref("grubert_fig3_na_LMEA_Arp.csv")
M=20001; dx=L/(M-1)
xs=[L*i/(M-1) for i in range(M)]
def cum(f):
    o=[0.0]
    for i in range(1,len(f)): o.append(o[-1]+0.5*(f[i]+f[i-1])*dx)
    return o
print("  where do the digitised curves actually MEET?")
for x in (0.30,0.35,0.38,0.40,0.42,0.45,0.50):
    e,a=interp(ne_ref,x),interp(nAr_ref,x)
    print("    x/L=%.2f  n_e=%.3f  n_Arp=%.3f  ratio n_e/n_Arp = %.3f"%(x,e,a,e/a))
print("\n  XT   fall-charge swing[Td]   Vgap for 20 Td bulk   Vgap for 30 Td")
for XT in (0.20,0.25,0.30,0.35,0.40,0.45):
    nAr=[];ne=[]
    for x in xs:
        xl=x/L; a=max(interp(nAr_ref,xl),FLOOR)
        e=max(interp(ne_ref,xl),FLOOR) if xl<XT else a
        nAr.append(a*SCALE); ne.append(e*SCALE)
    rho=[QE*(a-e) for a,e in zip(nAr,ne)]
    Irho=cum(rho); Esh=[v/EPS0 for v in Irho]; IIrho=cum(Esh)
    it=int(XT*(M-1))
    swing=abs(Esh[it])/N*1e21
    out=[]
    for tgt in (20.0,30.0):
        Et=-tgt*N/1e21
        out.append((Et-Esh[it])*L+IIrho[-1])
    print("  %.2f   %18.1f   %19.1f   %14.1f"%(XT,swing,out[0],out[1]))
