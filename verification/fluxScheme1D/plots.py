#!/usr/bin/env python3
"""Convergence plots for the fluxScheme1D verification bed."""
import re, math, collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rows = []
for ln in open('results.txt'):
    m = re.match(r'RESULT scheme=(\S+)\s+N=\s*(\d+)\s+h=(\S+)\s+Pe=(\S+)\s+'
                 r'PeGrid=(\S+)\s+L2=(\S+)\s+Linf=(\S+)', ln)
    if m:
        rows.append(dict(sch=m.group(1), N=int(m.group(2)), h=float(m.group(3)),
                         Pe=float(m.group(4)), Pg=float(m.group(5)),
                         L2=float(m.group(6)), Linf=float(m.group(7))))

SCH = ['std:upwind','std:linear','std:ROUNDF','std:ROUNDW','std:ROUNDL','ScharfetterGummel','CompleteFlux']
LBL = {'std:upwind':'upwind','std:linear':'linear (central)','std:ROUNDF':'ROUNDF',
       'std:ROUNDA':'ROUNDA','std:ROUNDAplus':'ROUNDAplus','std:ROUNDL':'ROUNDL','std:ROUNDW':'ROUNDW',
       'std:limitedLinear':'limitedLinear 1','std:MUSCL':'MUSCL',
       'std:SuperBee':'SuperBee','std:Minmod':'Minmod','std:vanLeer':'vanLeer',
       'ScharfetterGummel':'Scharfetter-Gummel','CompleteFlux':'Complete Flux (CFS)'}
STY = {'std:upwind':('#999999','o','--'), 'std:linear':('#444444','s',':'),
       'std:ROUNDF':('#d62728','^','-'), 'std:ROUNDA':('#ff7f0e','v','-'),
       'std:ROUNDAplus':('#8c564b','<','-'), 'std:ROUNDL':('#2ca02c','>','-'), 'std:ROUNDW':('#17a2b8','p','-'),
       'std:limitedLinear':('#9467bd','P','-.'), 'std:MUSCL':('#e377c2','X','-.'),
       'std:SuperBee':('#17becf','*','-.'), 'std:Minmod':('#bcbd22','h','-.'),
       'std:vanLeer':('#7f7f7f','d','-.'),
       'ScharfetterGummel':('#1f77b4','D','-'), 'CompleteFlux':('#000000','*','-')}

by = collections.defaultdict(list)
for r in rows: by[(r['sch'], r['Pe'])].append(r)
for k in by: by[k].sort(key=lambda r: r['N'])
Pes = sorted({r['Pe'] for r in rows})

# ---- Fig 1: convergence, one panel per domain Peclet -------------------
fig, axes = plt.subplots(2, 3, figsize=(16, 9.5))
for ax, Pe in zip(axes.ravel(), Pes):
    for s in SCH:
        d = [r for r in by.get((s, Pe), []) if r['L2'] > 0]
        if not d: continue
        c, mk, ls = STY[s]
        ax.loglog([r['h'] for r in d], [r['L2'] for r in d],
                  color=c, marker=mk, ls=ls, ms=5, lw=1.4, label=LBL[s])
    hs = [r['h'] for r in by.get((SCH[0], Pe), [])]
    if hs:
        h0, e0 = max(hs), 3e-2
        ax.loglog([h0, h0/16], [e0, e0/16],     'k-',  lw=0.8, alpha=.45)
        ax.loglog([h0, h0/16], [e0, e0/256],    'k--', lw=0.8, alpha=.45)
        ax.text(h0/5, e0/5,   ' 1st order', fontsize=7, alpha=.7)
        ax.text(h0/12, e0/90, ' 2nd order', fontsize=7, alpha=.7)
    ax.set_title(f"domain Pe = {Pe:g}   (grid Pe = {Pe/20:.3g} .. {Pe/640:.3g})",
                 fontsize=10)
    ax.set_xlabel('h'); ax.set_ylabel(r'$L_2$ error'); ax.grid(alpha=.3, which='both')
axes.ravel()[0].legend(fontsize=6.5, loc='lower right', ncol=2)
fig.suptitle("Flux-scheme convergence: steady 1D drift-diffusion WITH a source\n"
             r"$d/dx(vn - D\,dn/dx)=S$,  $n(0)=1$, $n(1)=2$, $S=1$, $D=1$",
             fontsize=12)
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig('fig1_convergence.png', dpi=140)

# ---- Fig 2: observed local order vs GRID Peclet ------------------------
fig2, ax = plt.subplots(figsize=(9, 6))
for s in SCH:
    xs, ys = [], []
    for Pe in Pes:
        d = [r for r in by.get((s, Pe), []) if r['L2'] > 0]
        for a, b in zip(d[:-1], d[1:]):
            xs.append(math.sqrt(a['Pg']*b['Pg']))
            ys.append(math.log(a['L2']/b['L2'])/math.log(2))
    if not xs: continue
    c, mk, ls = STY[s]
    o = sorted(zip(xs, ys))
    ax.semilogx([p[0] for p in o], [p[1] for p in o],
                color=c, marker=mk, ls=ls, ms=5, lw=1.3, label=LBL[s])
ax.axvline(1.0, color='k', ls=':', lw=1)
ax.text(1.05, -0.55, r'  $|Pe_{grid}|=1$', fontsize=9)
ax.axhline(1, color='k', lw=.6, alpha=.4); ax.axhline(2, color='k', lw=.6, alpha=.4)
ax.set_ylim(-1, 3.2)
ax.set_xlabel(r'grid Péclet  $|Pe_{grid}| = vh/D$')
ax.set_ylabel('observed local order  $p$')
ax.set_title("Observed order vs GRID Péclet\n"
             "SG is 2nd order for $Pe_{grid}\\ll1$ and 1st order for $Pe_{grid}\\gg1$",
             fontsize=11)
ax.grid(alpha=.3, which='both'); ax.legend(fontsize=7, ncol=3)
fig2.tight_layout(); fig2.savefig('fig2_order_vs_gridPeclet.png', dpi=140)

# ---- Fig 3: error at the finest mesh vs domain Peclet -----------------
fig3, ax = plt.subplots(figsize=(9, 6))
for s in SCH:
    xs, ys = [], []
    for Pe in Pes:
        d = [r for r in by.get((s, Pe), []) if r['N'] == 640 and r['L2'] > 0]
        if d: xs.append(Pe); ys.append(d[0]['L2'])
    if not xs: continue
    c, mk, ls = STY[s]
    ax.loglog(xs, ys, color=c, marker=mk, ls=ls, ms=6, lw=1.4, label=LBL[s])
ax.set_xlabel(r'domain Péclet  $Pe = vL/D$')
ax.set_ylabel(r'$L_2$ error at $N=640$')
ax.set_title("Accuracy at fixed resolution (N=640) across the Péclet range",
             fontsize=11)
ax.grid(alpha=.3, which='both'); ax.legend(fontsize=7, ncol=2)
fig3.tight_layout(); fig3.savefig('fig3_error_vs_Peclet.png', dpi=140)
print("wrote fig1_convergence.png fig2_order_vs_gridPeclet.png fig3_error_vs_Peclet.png")
