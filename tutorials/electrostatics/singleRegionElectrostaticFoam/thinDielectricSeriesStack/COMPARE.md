# thinDielectricPotential -- analytic series stack (unit bed)

## Question
Does `thinDielectricPotential` reproduce the exact potential of a gas gap in
series with a thin capacitive layer -- including the CAPACITANCE VALUE and the
SURFACE CHARGE with its sign, not merely the limiting cases?

## Bed
1D column, gas gap L = 1 m (gas epsilonR = 1), driven electrode `top` at
V0 = 1 V, `sides` zeroGradient (which is what makes it exactly 1D), `bottom` =
the condition under test.

NOTE: the shipped sibling case plate2D_timeVaryingBC bundles left+right+bottom
into ONE `walls` patch, which makes it 2D Laplace with NO analytic answer. The
patch split here is load-bearing: a bed is only ground truth if its geometry
matches the formula.

## Ground truth (ANALYTIC -- no baseline case required)
With V = A y + B, V(1) = V0, outward normal at `bottom` = -y:
    eps_g dV/dn = sigma - C (V_s - V_b),    C = eps0*epsR/d
    =>  V_s = (sigma + C V_b + eps_g V0) / (eps_g + C)
    and for C = 0 (free-standing):  V_s = V0 + sigma/eps_g

| case | epsR | d           | sigma    | expected V_s |
|------|------|-------------|----------|--------------|
| a | 5 | 0.5            | 0        | 0.0909090909 |
| b | 5 | 5.0            | 0        | 0.5          |
| c | 2 | 0.5            | 0        | 0.2          |
| d | - | free-standing  | 0        | 1.0          |
| e | - | free-standing  | +eps0    | 2.0          |
| f | - | free-standing  | -eps0    | 0.0          |
| g | 5 | 0.5            | +eps0    | 0.1818181818 |

eps0 = 8.854187817620389e-12, i.e. 1/(mu0 c^2) as plasmaConstants.C builds it.

Errors are judged against the VOLTAGE SCALE V0 = 1 V, never against the expected
value: cases (d) and (f) expect exactly zero, where a relative error is
undefined and a converged zero reads as a catastrophic failure.

Plus two NEGATIVE tests -- the guards must fire:
    bad_pair   `thickness` alone                 -> "must be given TOGETHER"
    bad_eps    `epsilonR` with no `thickness`    -> "has NO EFFECT on a free-standing"

The d -> 0 DIRICHLET limit is covered elsewhere: measured on needleDBD
2026-09-03, `thickness 1e-12` + `backingPotential 500` gave 500 V exactly.

## Extraction
    ./Allrun-sweep      # prints measured vs analytic and PASS/FAIL per row;
                        # exit status is nonzero if ANY row fails
    ./Allclean

## Requirement
`singleRegionElectrostaticFoam` does NOT link libplasmaBcs.so, so
system/controlDict loads it via `libs`. Without that the failure reads
`Unknown patchField type thinDielectricPotential`, which looks like a build
failure and is not.
