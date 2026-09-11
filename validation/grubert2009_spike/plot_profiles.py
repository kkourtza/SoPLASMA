#!/usr/bin/env python3
"""Publication-ready centreline profiles through the Grubert current spike.

Reads profiles/*.csv written by dump_profiles.py -- the data files are
independent of this script, so edit freely.

    ~/ct-env/bin/python plot_profiles.py [profiles_dir]

Writes  profiles_logx.png   (log x: cathode layer AND bulk in one view)
        profiles_cathode.png (linear zoom, x/L < 0.05)

Written 2026-09-07.  Grubert et al., Plasma Sources Sci. Technol. 18 (2009)
015006 -- steady normal-glow reference values are marked.
"""
import sys, os, glob, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = sys.argv[1] if len(sys.argv) > 1 else "profiles"
GRU = dict(ne=(0.394, 2.478e15), ni=(0.203, 6.391e15), EN_fall=2450.0)
CLAMP_TE = 100.0*2/3

files = sorted(glob.glob(os.path.join(D, "*.csv")), key=lambda p: float(os.path.basename(p)[:-4]))
data = []
for p in files:
    t = float(os.path.basename(p)[:-4])
    rows = [r for r in csv.DictReader(l for l in open(p) if not l.startswith("#"))]
    data.append((t, {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}))

cmap = plt.get_cmap("plasma")
cols = [cmap(i/max(len(data)-1, 1)*0.88) for i in range(len(data))]

PANELS = [
    ("n_e",           r"$n_e$   [m$^{-3}$]",                          "log"),
    ("n_Arp",         r"$n_{\mathrm{Ar}^+}$   [m$^{-3}$]",            "log"),
    ("EN_Td",         r"$E/N$   [Td]",                                "log"),
    ("Te_eV",         r"$T_e=\frac{2}{3}\langle U\rangle$   [eV]",   "linear"),
    ("chargeDensity", r"$\rho$   [C m$^{-3}$]",                       "symlog"),
    ("ePotential_V",  r"$\phi$   [V]",                                "linear"),
]

def draw(xlog, xlim, out, title):
    fig, axs = plt.subplots(2, 3, figsize=(15.5, 8.4))
    for a, (key, lab, scale) in zip(axs.ravel(), PANELS):
        for (t, d), c in zip(data, cols):
            a.plot(d["x_over_L"], d[key], color=c, lw=1.5, label=f"{t*1e6:.5f}")
        if scale == "symlog":
            a.set_yscale("symlog", linthresh=1e-7)
        else:
            a.set_yscale(scale)
        a.set_xscale("log" if xlog else "linear")
        a.set_xlim(*xlim)
        a.set_xlabel(r"$x/L$    (0 = cathode, 1 = anode)")
        a.set_ylabel(lab)
        a.grid(alpha=0.3, lw=0.5, which="both")
        a.tick_params(labelsize=9)
    axs[0,0].plot([GRU["ne"][0]], [GRU["ne"][1]], "*", ms=17, color="#00b050",
                  mec="k", mew=0.7, zorder=6, label="Grubert steady")
    axs[0,1].plot([GRU["ni"][0]], [GRU["ni"][1]], "*", ms=17, color="#00b050",
                  mec="k", mew=0.7, zorder=6, label="Grubert steady")
    axs[0,2].axhline(GRU["EN_fall"], color="#00b050", ls="--", lw=1.4,
                     label="normal cathode fall, 2450 Td")
    axs[1,0].axhline(CLAMP_TE, color="crimson", ls=":", lw=1.4,
                     label=r"energy clamp (100 eV)")
    axs[1,2].axhline(0.0, color="k", lw=0.6)
    for a in (axs[0,0], axs[0,1], axs[0,2], axs[1,0]):
        a.legend(fontsize=8, loc="best", framealpha=0.92)
    h, l = axs[0,0].get_legend_handles_labels()
    fig.legend(h[:len(data)], l[:len(data)], fontsize=8.5, title=r"$t$  [$\mu$s]",
               title_fontsize=9, loc="center right", bbox_to_anchor=(1.0, 0.5),
               framealpha=0.95)
    fig.suptitle(title, fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 0.915, 0.955))
    fig.savefig(out, dpi=165)
    print("wrote", out)

xmin = min(d["x_over_L"][0] for _, d in data)
draw(True,  (xmin, 1.0),  "profiles_logx.png",
     "Grubert dc glow, 1e8 $\\Omega$ ballast: centreline profiles through the current spike into the "
     "frozen state ($\\Delta t\\to6\\times10^{-14}$ s).  LOG $x$ -- the whole event lives in $x/L<0.03$")
draw(False, (0.0, 0.05), "profiles_cathode.png",
     "Same data, LINEAR zoom on the cathode layer ($x/L<0.05$, i.e. the first 500 $\\mu$m of a 10 mm gap)")
