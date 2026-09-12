# fluxScheme2Dnonortho -- non-orthogonal comparison

2D manufactured solution, exact value imposed on EVERY patch at each face's own
centre, so the boundary data is consistent with whatever shape the mesh has:

    n(x,y) = 1 + x + A sin(pi x) sin(pi y)
    s      = v dn/dx - D lap(n)

`shear` offsets the top edge in x. shear = 0 is the orthogonal control;
shear = 0.20 gives 11.3 deg non-orthogonality, skewness 0.2, `Mesh OK`.
Non-orthogonal correctors are ITERATED WITH RE-ASSEMBLY (see below) and the
count used is reported as `nCorr`.

## THE HEADLINE: SG AND CFS CONVERGE ON A NON-ORTHOGONAL MESH

**Measured 2026-09-12.** Two consecutive runs bit-identical; 5.1 s, 42 MB.

    Pe = 1, shear = 0.20 (11.3 deg)
    scheme                      N=100     400       1600      6400     order  nCorr
    standard (linear+corrected) 5.141e-3  2.620e-3  1.335e-3  6.752e-4  0.97  15-16
    ScharfetterGummel           5.282e-3  2.635e-3  1.347e-3  6.845e-4  0.98  15-16
    CompleteFlux                5.745e-3  2.827e-3  1.403e-3  6.994e-4  1.01  15-16

**The control is reproducible again.** `standard` was absent from `Allrun`'s
hardcoded scheme loop, so the bed could not produce the very row its headline
compared against -- a comparison whose control the harness cannot regenerate
(A1). `SCHEMES` now defaults to all three and the control reproduces the row
above to every digit it quotes. SG and CFS bypass `divSchemes` entirely, so the
per-scheme `div(phi,n)` only ever reaches `standard`.

All three converge at first order -- expected for uniform skew with a corrected
treatment -- and SG and CFS now sit within 2-4% of `standard` at every
refinement instead of on a fixed floor.

On the ORTHOGONAL control (shear = 0) both behave as in the 1D bed:

    Pe = 1     SG   order 2.00,  CFS order 2.00, CFS 3.9x more accurate
    Pe = 100   SG   order 1.86,  CFS order 1.99, CFS 4.5x more accurate

## SUPERSEDED BY THE ABOVE -- the original finding, and why it is kept

**This bed's original headline, measured 2026-09-08, was the exact opposite:
"SG AND CFS DO NOT CONVERGE ON A NON-ORTHOGONAL MESH".**

    Pe = 1, shear = 0.20            N=100     400       1600      6400    order  nCorr
    ScharfetterGummel (2026-09-08)  1.009e-2  1.001e-2  1.004e-2  1.004e-2  0.00   2
    CompleteFlux      (2026-09-08)  1.045e-2  1.012e-2  1.006e-2  1.005e-2  0.00   2

That measurement was CORRECT and so was its diagnosis: the operator used the
ORTHOGONAL delta coefficients, so it carried no non-orthogonal treatment at all,
and the two correctors it "converged" in were two passes over a matrix with no
correction term for them to act on.

**THE DEFECT HAS SINCE BEEN FIXED**, in `ScharfetterGummel.H` (see its comment at
lines 117-140, which records this same history at the point of use). The fix is
exact rather than a patch: with the Bernoulli identity `B(-z) = B(z) + z` the SG
face flux factorises as

    Gamma = coeffP n_P - coeffN n_N = phi n_P - Df B(Pe) magSf snGrad(n)

so SG's diffusive part IS an ordinary diffusion flux with effective diffusivity
`D*B(Pe)`. It therefore takes the standard non-orthogonal treatment --
`nonOrthDeltaCoeffs` in the implicit part plus the explicit `snGrad` correction
scaled by `Df*B(Pe)`. CFS inherits it, its homogeneous part being the SG flux.

The rows above are kept, dated and marked, rather than deleted: a finding that
quietly disappears leaves no way to recognise its stale copies elsewhere (D2).
**They are the pre-fix state and must not be quoted as current.**

The README carried the pre-fix headline for four days after the fix landed,
while `results.txt` beside it already showed order ~1. Nothing compared the two
-- which is why `/regression-gate` now exists.

## CONSEQUENCES

  * SG and CFS are verified on orthogonal AND on uniformly skewed meshes.
    The 1D bed's results (`../fluxScheme1D`, all Cartesian) stand unchanged.
  * The Grubert dc-glow cases are Cartesian and were never affected either way.
  * **Still UNVERIFIED, and do not extrapolate to it:** genuinely unstructured
    or strongly graded meshes. This bed tests ONE uniform shear at 11.3 deg.
    First order there is not a claim about a tetrahedral mesh or about skewness
    that varies cell to cell.
  * `nCorr` 15-16 for all three schemes is now the expected signature. **A run
    reporting nCorr 2 for SG or CFS on a skewed mesh has lost the correction
    and is the pre-fix operator** -- treat it as a regression, not a speed-up.

## Two harness bugs found and fixed on the way (do not repeat them)

  1. FIRST ATTEMPT KEPT THE 1D SOLUTION n(x) ON A SHEARED BLOCK. Shearing turns
     the left/right patches into SLANTED planes, so "n = a at the left
     boundary" stops being a condition at constant x. The boundary data then
     disagreed with the geometry by O(shear), and the measured floor tracked the
     shear ~1:1 (0.01 -> 8.0e-3, 0.05 -> 4.7e-2, 0.10 -> 1.07e-1, 0.20 ->
     2.33e-1) for BOTH schemes at EVERY Peclet. The shear = 0 row reproducing
     the 1D answer to all digits is what proved the harness itself was sound.
     Fixed by the 2D manufactured solution above.

  2. NON-ORTHOGONAL CORRECTORS MUST RE-ASSEMBLE THE EQUATION, not just re-solve
     it. `corrected` puts the non-orthogonal part in the matrix SOURCE at
     assembly time; calling solve() again on the same matrix changes nothing --
     measured, identical digits. The loop now rebuilds the equation each pass
     and exits on residual < 1e-13.

  Both were real, and neither was the cause of the original headline: it
  survived both, and was only removed by the operator fix above.
