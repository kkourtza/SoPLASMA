#!/usr/bin/env python3
"""Check a two-region series dielectric stack against its analytic solution.

For the plate2D electrostatics tutorials: a gas layer and a dielectric layer in
series, uniform fields, an optional surface charge at the interface. This has a
closed form, so it is a control against ANALYTIC ground truth rather than
against another run of the same code.

WHAT IT DISCRIMINATES

  1. Each region's permittivity really being read from
     constant/<region>/electricalProperties. If one silently defaults to 1.0 the
     run still converges and the field is smooth and plausible -- nothing in the
     log distinguishes it -- but the interface potential moves by a factor of
     order eps2/eps1.

  2. The surface-charge jump condition at the interface, when the case sets a
     non-zero surfCharge.

  3. The linear solve actually converging. An unconverged solve reports as a
     normal solve (reaching `maxIter` is not an error) and shifted the interface
     potential by 7.3% in this very tutorial -- see the measured table in
     system/fvSolution-foam of plate2D_timeVaryingBC_implicitBoundary.

EVERY NUMBER IS READ FROM THE CASE, not hardcoded, so the check stays valid
across the three plate2D variants, which differ in applied voltage, surface
charge and coupling mode. Geometry is the one exception and is asserted below.

Usage, from inside a case directory:
    ../../../../tools/check_series_stack.py [time]
"""

import os
import re
import sys

EPS0 = 8.8541878128e-12         # [F/m]

# Geometry, from configuration/geometry/mesh.py: a unit square split at x_end/2.
# Not read from the script because it is expressed there as gmsh API calls.
D1 = 0.5                        # gas thickness [m]
D2 = 0.5                        # dielectric thickness [m]

GAS = "leftGas"
DIELECTRIC = "rightDielectric"
INTERFACE = "leftGas_to_rightDielectric"

# The floor on any comparison taken from a written field: controlDict has
# `writePrecision 6`, so the file holds 6 significant digits and a value cannot
# be trusted below ~1e-6 RELATIVE regardless of how well the run converged.
WRITE_PRECISION_FLOOR = 5e-6

# Segregated coupling lags the interface -- each region sees the other's
# potential from the previous corrector -- so it converges to the coupling
# tolerance, not to the write precision. Measured on
# plate2D_timeVaryingBC_explicitBoundary 2026-09-01: 8.8e-5 relative.
SEGREGATED_TOL = 5e-4


def read_config(path="configuration/config"):
    """The case's single-source config, as a flat dict of strings."""
    if not os.path.exists(path):
        sys.exit("No %s -- run this from inside a case directory." % path)
    out = {}
    for line in open(path):
        line = line.split("//")[0].strip()
        if not line or ";" not in line:
            continue
        parts = line.rstrip(";").split(None, 1)
        if len(parts) == 2:
            out[parts[0]] = parts[1].strip().strip('"').rstrip(";").strip()
    return out


def latest_time():
    times = [
        (float(d), d)
        for d in os.listdir(".")
        if _isnum(d) and os.path.isdir(os.path.join(d, GAS))
    ]
    if not times:
        sys.exit("No time directory containing %s/ -- run ./Allrun-serial first."
                 % GAS)
    return max(times)[1]


def _isnum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def epsilon_r(region):
    """The permittivity as the SOLVER reads it: from the region's own file."""
    path = os.path.join("constant", region, "electricalProperties")
    if not os.path.exists(path):
        sys.exit(
            "Missing %s.\n"
            "This file is where a region's permittivity lives (one small file\n"
            "per region, as OpenFOAM does for thermophysicalProperties). If it\n"
            "vanished, check whether Allclean removed constant/<region>/\n"
            "wholesale instead of just constant/<region>/polyMesh." % path
        )
    src = open(path).read()
    m = re.search(r"^\s*epsilonR\s+(\S+?)\s*;", src, re.M)
    if not m:
        sys.exit("No `epsilonR` entry in %s." % path)
    val = m.group(1)
    if val.startswith("$"):
        cfg = read_config()
        key = val[1:]
        if key not in cfg:
            sys.exit("%s refers to $%s, which configuration/config does not "
                     "define. A dangling $var is inert and silently leaves the "
                     "permittivity unset." % (path, key))
        val = cfg[key]
    return float(val)


def patch_values(path, patch):
    src = open(path).read()
    m = re.search(r"^\s*" + re.escape(patch) + r"\s*\n\s*\{(.*?)\n\s*\}",
                  src, re.S | re.M)
    if not m:
        sys.exit(
            "Patch `%s` is not in %s.\n"
            "If it is missing entirely, the region interface was never built:\n"
            "splitMeshRegions makes it from INTERNAL faces between cellZones,\n"
            "so naming those faces as a boundary patch in the mesh leaves the\n"
            "two regions uncoupled with no error." % (patch, path))
    blk = m.group(1)

    u = re.search(r"\bvalue\s+uniform\s+([-\d.eE+]+)\s*;", blk)
    if u:
        return [float(u.group(1))]
    n = re.search(r"\bvalue\s+nonuniform[^(]*\((.*?)\)\s*;", blk, re.S)
    if n:
        return [float(x) for x in n.group(1).split()]
    sys.exit("No `value` entry on patch `%s` in %s." % (patch, path))


def v_interface(v0, sigma, e1, e2):
    """Interface potential of a 1-D two-layer stack.

    Continuity of the normal displacement with a free surface charge sigma at
    the interface, D2 - D1 = sigma, plus V(0) = v0 and V(d1+d2) = 0:

        eps0 (e2 Vi/d2 + e1 Vi/d1) = sigma + eps0 e1 v0/d1
    """
    return (sigma / EPS0 + e1 * v0 / D1) / (e2 / D2 + e1 / D1)


def main():
    cfg = read_config()
    t = sys.argv[1] if len(sys.argv) > 1 else latest_time()

    e1 = epsilon_r(GAS)
    e2 = epsilon_r(DIELECTRIC)
    sigma = float(cfg.get("surfCharge", 0.0))
    v0 = float(cfg["dirichletValueLeft"])
    implicit = cfg.get("useImplicit", "true").lower() == "true"

    # The applied value is a ramp reaching dirichletValueLeft at endTime, so the
    # analytic reference is only the full v0 at the final time.
    end = float(cfg["endTime"])
    if abs(float(t) - end) > 1e-12:
        v0 *= float(t) / end
        print("NOTE: t = %s is mid-ramp, so V_applied = %g, not %s"
              % (t, v0, cfg["dirichletValueLeft"]))

    vals = patch_values(os.path.join(t, GAS, "ePotential"), INTERFACE)
    lo, hi = min(vals), max(vals)
    got = 0.5 * (lo + hi)

    expected = v_interface(v0, sigma, e1, e2)
    defaulted = v_interface(v0, sigma, 1.0, 1.0)

    tol = WRITE_PRECISION_FLOOR if implicit else SEGREGATED_TOL

    print("case                    %s" % os.path.basename(os.getcwd()))
    print("time                    %s   (endTime %g)" % (t, end))
    print("coupling                %s"
          % ("MONOLITHIC (useImplicit true)" if implicit
             else "SEGREGATED (useImplicit false)"))
    print("epsilonR  %-14s %g   <- constant/%s/electricalProperties"
          % (GAS, e1, GAS))
    print("epsilonR  %-14s %g   <- constant/%s/electricalProperties"
          % (DIELECTRIC, e2, DIELECTRIC))
    print("V_applied               %g" % v0)
    print("surfCharge              %g" % sigma)
    print()
    print("interface faces         %d" % len(vals))
    print("V_interface  min/max    %.9g / %.9g" % (lo, hi))
    print()
    print("analytic                %.9g" % expected)
    print("measured                %.9g" % got)
    print("if BOTH eps were 1.0    %.9g" % defaulted)
    print()

    if hi > lo and (hi - lo) / max(abs(hi), 1e-300) > 1e-6:
        print("INCONCLUSIVE: the interface potential is not uniform "
              "(relative spread %.2e)." % ((hi - lo) / abs(hi)))
        print("  The 1-D reference assumes zeroGradient on top and bottom, so")
        print("  fix that before reading the verdict below.")
        print()

    rel = abs(got - expected) / max(abs(expected), 1e-300)

    if rel < tol:
        # Report what the measurement IMPLIES, not just a boolean: invert the
        # formula for eps2/eps1 so a pass names the permittivity it recovered.
        print("PASS   relative error %.2e < %.0e" % (rel, tol))
        if sigma == 0.0:
            # eps2/eps1 = (d2/d1)((v0 - Vi)/Vi), valid only at zero surface
            # charge; with sigma != 0 the two cannot be separated from Vi alone.
            ratio = (D2 / D1) * ((v0 - got) / got)
            print("       recovered eps2/eps1 = %.6f (declared %.6f)"
                  % (ratio, e2 / e1))
        else:
            print("       (eps2/eps1 cannot be inverted from V_interface alone"
                  " at non-zero surface charge)")
        if not implicit:
            print("       tolerance is the SEGREGATED one: this coupling lags")
            print("       the interface, so it converges to the coupling")
            print("       tolerance rather than to the write precision.")
        return 0

    print("FAIL   relative error %.2e (tolerance %.0e)" % (rel, tol))
    if abs(got - defaulted) / max(abs(defaulted), 1e-300) < 1e-3:
        print()
        print("       The measured value is the eps = 1/1 value: a")
        print("       permittivity was NOT read. Check both")
        print("       constant/<region>/electricalProperties files exist.")
    else:
        print()
        print("       Not the defaulted-permittivity signature either. Check")
        print("       the linear solve converged -- reaching `maxIter` is")
        print("       reported as a normal solve, and cost 7.3% here:")
        print("         grep 'Solving for ePotential' logs/log.*serial")
    return 1


if __name__ == "__main__":
    sys.exit(main())
