# fluxScheme2Dnonortho -- non-orthogonal comparison

2D manufactured solution, exact value imposed on EVERY patch at each face's own
centre, so the boundary data is consistent with whatever shape the mesh has:

    n(x,y) = 1 + x + A sin(pi x) sin(pi y)
    s      = v dn/dx - D lap(n)

`shear` offsets the top edge in x. shear = 0 is the orthogonal control;
shear = 0.20 gives 11.3 deg non-orthogonality, skewness 0.2, `Mesh OK`.
Non-orthogonal correctors are ITERATED WITH RE-ASSEMBLY (see below) and the
count used is reported as `nCorr`.

## THE HEADLINE: SG AND CFS DO NOT CONVERGE ON A NON-ORTHOGONAL MESH

    Pe = 1, shear = 0.20 (11.3 deg)
    scheme                       NX=10      20        40        80      order  nCorr
    standard (linear+corrected)  5.141e-3  2.620e-3  1.335e-3  6.752e-4  0.97   15-16
    ScharfetterGummel            1.009e-2  1.001e-2  1.004e-2  1.004e-2  0.00   2
    CompleteFlux                 1.045e-2  1.012e-2  1.006e-2  1.005e-2  0.00   2

The standard scheme converges (first order, expected for uniform skew with
`corrected`). SG and CFS sit on a FIXED ERROR FLOOR that refinement does not
touch.

THE CAUSE IS IN THE OPERATOR, not the bed:

    ScharfetterGummel.H:116
      surfaceScalarField diffCond = Df * mesh.magSf() * mesh.deltaCoeffs();

`mesh.deltaCoeffs()` are the ORTHOGONAL delta coefficients. There is no
`nonOrthDeltaCoeffs`, no `nonOrthCorrectionVectors`, no correction of any kind.
CFS inherits this because its homogeneous part IS the SG flux. That is also why
they converge in 2 correctors while the standard scheme needs 15-16: their
matrix has no explicit correction term for the correctors to iterate on.

On the ORTHOGONAL control (shear = 0) both behave exactly as in the 1D bed:

    Pe = 1     SG   order 2.00,  CFS order 2.00, CFS 3.9x more accurate
    Pe = 100   SG   order 1.86,  CFS order 1.99, CFS 4.5x more accurate

## CONSEQUENCES

  * SG and CFS are sound on ORTHOGONAL / Cartesian meshes only. Every result in
    `../fluxScheme1D` is on such a mesh and stands.
  * The Grubert dc-glow cases are Cartesian, so the CFS result there is not
    affected by this.
  * DO NOT use SG or CFS on an unstructured or graded-skewed mesh until the
    operators carry a non-orthogonal correction. This is a prerequisite for any
    2D/3D application work, and for making either scheme a default.
  * `standard` + a limiter remains the only verified option on skewed meshes.

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

  Neither bug changes the headline: it survived both fixes, and the standard
  scheme converges on the same mesh with the same harness.
