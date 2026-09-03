# thinDielectricPotential -- analytic series stack (unit bed)

## Question
Does `thinDielectricPotential` reproduce the exact potential of a gas gap in
series with a thin capacitive layer, including the CAPACITANCE VALUE (not just
its limits)?

## Bed
1D column, gas gap L = 1 m (epsilonR = 1), driven electrode `top` at V0 = 1 V,
`sides` zeroGradient (makes it exactly 1D), `bottom` = the condition under test.
sigma = 0 everywhere, so the answer depends only on the capacitance.

## Ground truth (analytic, no baseline case needed)
D-continuity across the layer:  eps0 (V0 - Vs)/L = eps0 epsR (Vs - Vb)/d
    Vs = d V0 / (d + epsR L)          with Vb = 0

| case | epsilonR | thickness d | expected Vs |
|------|----------|-------------|-------------|
| a | 5 | 0.5 | 0.0909090909 |
| b | 5 | 5.0 | 0.5          |
| c | 2 | 0.5 | 0.2          |
| d | 5 | (free-standing, C=0) | 1.0 exactly -- pure Neumann, sigma=0, so V is uniform |

Case (d) is the C -> 0 limit; the d -> 0 Dirichlet limit was verified separately
on needleDBD 2026-09-03 (patch value 500 V from backingPotential 500).

## Extraction
    ./Allrun-sweep        # prints measured vs analytic and a PASS/FAIL per row
