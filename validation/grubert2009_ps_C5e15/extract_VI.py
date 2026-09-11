#!/usr/bin/env python3
"""Extract the V-I characteristic traced by a currentSource ramp.

WHY THIS IS THE POINT OF A CURRENT SOURCE. Under voltage control the
normal-glow branch is not a function -- V is flat while j varies, so a
voltage-driven case cannot resolve it at all. Ramping the CURRENT slowly makes
the discharge track its own characteristic quasi-statically, so the whole curve
comes out of one run:

    Townsend     V rises steeply with I
    SUBNORMAL    V FALLS with I          <- negative differential resistance
    normal glow  V flat, j pinned at j_n

Emits VI.csv: I_set, I_cond, j_cond, V_gap, dV/dj, branch.

`dV/dj` is the discriminator: its SIGN separates the branches, and the sign
change is a stronger validation than matching any single operating point,
because no fitted coefficient can fake it.

Reads postProcessing/externalCircuit/circuit.csv, written every step, so the
curve is as resolved as the run. Plot with plot_VI.py -- separate, so the
figure can be restyled without re-running (rule 6).
"""
import sys, csv

AREA_CM2 = 2.0e-3          # 5 faces of 200 x 200 um
REF_J, REF_V = -0.511, -500.0      # Grubert et al. 2009 operating point

def main():
    src = sys.argv[1] if len(sys.argv) > 1 \
        else "postProcessing/externalCircuit/circuit.csv"
    out = sys.argv[2] if len(sys.argv) > 2 else "VI.csv"
    nbin = int(sys.argv[3]) if len(sys.argv) > 3 else 200

    rows = []
    for line in open(src):
        if line.startswith("#"):
            continue
        f = line.strip().split(",")
        if len(f) < 5:
            continue
        try:
            rows.append([float(x) for x in f[:5]])
        except ValueError:
            continue          # a torn line mid-write; skip it
    if len(rows) < 3:
        raise SystemExit(f"{src}: only {len(rows)} usable rows")

    # BIN IN I_set, not in time. The ramp is linear in time but the
    # characteristic is a function of current, and binning in the independent
    # variable is what keeps dV/dj from being dominated by step-to-step noise.
    lo, hi = rows[0][1], rows[-1][1]
    if hi == lo:
        raise SystemExit("I_set never changed -- not a ramp")
    bins = {}
    for t, iset, icond, v, g in rows:
        k = min(nbin - 1, int(abs((iset - lo)/(hi - lo))*nbin))
        b = bins.setdefault(k, [0, 0.0, 0.0, 0.0])
        b[0] += 1; b[1] += iset; b[2] += icond; b[3] += v
    keys = sorted(bins)
    pts = [(bins[k][1]/bins[k][0], bins[k][2]/bins[k][0], bins[k][3]/bins[k][0])
           for k in keys]

    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["I_set_A", "I_cond_A", "j_cond_mA_cm2", "V_gap_V",
                    "dV_dj", "branch"])
        for n, (iset, icond, v) in enumerate(pts):
            j = icond/AREA_CM2*1e3
            if 0 < n < len(pts)-1:
                jm = pts[n-1][1]/AREA_CM2*1e3
                jp = pts[n+1][1]/AREA_CM2*1e3
                dvdj = (pts[n+1][2]-pts[n-1][2])/(jp-jm) if jp != jm else float("nan")
            else:
                dvdj = float("nan")
            # branch from the SIGN of dV/dj, in the physical convention
            # (both V and j negative here, so work with magnitudes)
            if dvdj != dvdj:      br = "-"
            elif dvdj > 1e-9:     br = "rising (Townsend)"
            elif dvdj < -1e-9:    br = "FALLING (subnormal, NDR)"
            else:                 br = "flat (normal glow)"
            w.writerow([f"{iset:.6e}", f"{icond:.6e}", f"{j:.6e}",
                        f"{v:.4f}", f"{dvdj:.6e}", br])
    print(f"{len(pts)} binned points -> {out}")

    # the peak in |V| is the Townsend -> subnormal transition
    kmax = max(range(len(pts)), key=lambda k: abs(pts[k][2]))
    print(f"  |V| peak {pts[kmax][2]:.2f} V at j = "
          f"{pts[kmax][1]/AREA_CM2*1e3:.4g} mA/cm2")
    if kmax in (0, len(pts)-1):
        print("  -> the peak is at an ENDPOINT: the run has not yet traversed"
              " the transition, so no branch conclusion can be drawn.")
    else:
        print("  -> an INTERIOR peak: the characteristic turns over here,"
              " which is the subnormal (NDR) branch.")
    # how close does the curve pass to the reference point?
    kref = min(range(len(pts)),
               key=lambda k: abs(pts[k][1]/AREA_CM2*1e3 - REF_J))
    jr = pts[kref][1]/AREA_CM2*1e3
    print(f"  nearest point to Grubert's j = {REF_J} mA/cm2: "
          f"j = {jr:.4g}, V = {pts[kref][2]:.2f} V (ref {REF_V})")
    if abs(jr - REF_J) > 0.3*abs(REF_J):
        print("  -> the ramp has NOT yet reached the reference current;"
              " that comparison is not yet meaningful.")

if __name__ == "__main__":
    main()
