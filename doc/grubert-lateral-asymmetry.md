# The grubert lateral asymmetry: diagnosed 2026-09-10

**Symptom.** The user visualised `grubert_long` in ParaView and found the
solution NOT symmetric about the gap centreline, on a mesh and boundary
conditions that are symmetric. Measured: 25% in `n_e`, 41% in `nEps_e`, 23% in
`chargeDensity`.

**Conclusion.** NUMERICAL, and specific: an imperfectly mirror-symmetric MESH,
whose asymmetry the gradient reconstruction turns into a spurious TRANSVERSE
field, which the ionisation feedback then amplifies exponentially.

**It does not block 2-D/3-D work.** The transverse-field error falls as
`NY^-1.7`, so a properly refined lateral direction cures it. The pathological
case was specifically `NY = 5` -- coarse enough laterally for a large metric
asymmetry, yet with enough cells to admit an antisymmetric mode. Worst of both
worlds.

## The chain, every link measured

| # | link | measurement |
|---|---|---|
| 1 | the mesh is not bit-symmetric | gmsh `Bump 0.02` grading + extrusion put the midplane at `y = 0.000500000000001`; **982 of 1592** cell-volume mirror pairs not bit-identical |
| 2 | the gradient turns that into a transverse field | `Ey` anti-symmetry violation `3.0e-07` of \|E\|, **immune to solver tolerance and solver type** |
| 3 | the Poisson solve adds its own error | `phi` `1.8e-09`, `Ex` `3.1e-08` at `tol 1e-10`; these DO follow tolerance |
| 4 | `rho` amplifies ~1000x | bulk quasi-neutral to **6 ppm** (`net/n_e = 5.96e-06`), so `rho = e(n_Arp + n_Ar2p - n_e)` is a catastrophic cancellation, and `rho` IS the Poisson source |
| 5 | the ionisation feedback grows it | perturbation e-folds every **1.46 ns** against the mean's **2.32 ns** -- **1.59x faster**; rate invariant under 16x dt and 10x cell count |

## The decisive experiment, and why it was needed

Every field in the running case is asymmetric, and `rho -> E -> species -> rho`
is a LOOP, so no measurement of the coupled state can say which link introduces
the asymmetry. The user put it exactly: *"even with measuring E you cant really
understand whats wrong cause the space charge asymmetry also reflects to E --
theres a two way coupling"*.

`src/applications/utilities/testPoissonSymmetry` breaks the loop: it symmetrises
`rho` to BIT equality in memory, performs one Poisson solve with the case's own
`fvSchemes`/`fvSolution`, and reports the mirror asymmetry of `phi` and of
`E = -grad(phi)`. Any asymmetry in the output was manufactured by the solve.

That is what separated link 2 (structural, tolerance-immune) from link 3
(solver error, tolerance-following) -- a distinction invisible in the coupled
run.

## Measurements that decide the actions

    Poisson tolerance and solver, symmetric rho in:
      GAMG 1e-10     7 iters    phi 1.8e-09   Ex 3.1e-08   Ey 3.0003e-07
      GAMG 1e-16  2000 iters    phi 5.2e-12   Ex 1.8e-10   Ey 3.0007e-07
      PCG  1e-10     3 iters    phi 6.0e-09   Ex 2.1e-08   Ey 3.0007e-07
      PCG  1e-16    20 iters    phi 5.2e-12   Ex 1.8e-10   Ey 3.0007e-07

    the MESH, at tol 1e-16:
      original mesh              phi 5.2e-12   Ex 1.8e-10   Ey 3.0007e-07
      points snapped symmetric   phi 2.5e-12   Ex 2.0e-12   Ey 8.6e-12   (35,000x)

    lateral refinement, PCG 1e-14:
      NY = 5  (2,000 cells)      Ey 3.00e-07
      NY = 16 (20,240 cells)     Ey 4.37e-08     -> Ey ~ NY^-1.7

## Actions taken, and the ones deliberately NOT taken

* **1-D validation runs use `NY = 1`.** No lateral degree of freedom exists, so
  no antisymmetric mode can be represented, let alone grow. `validation/grubert_1d`.
* **2-D/3-D uses normal refinement.** `Ey ~ NY^-1.7` means resolution cures it.
* **The Poisson tolerance is NOT tightened.** It improves `phi`/`Ex` but leaves
  `Ey` -- the component that drives the instability -- completely unchanged,
  and GAMG at `1e-16` costs 2000 iterations (maxIter, 286x) for nothing. Note
  also that `1e-16` ABSOLUTE on a 225 V field is `4.4e-19` relative, below
  double precision, hence unreachable by construction.
  (If a tight Poisson solve is ever wanted: PCG reached `8.5e-17` in 20
  iterations where GAMG needed 2000 for `1.2e-15`.)
* **No mesh-symmetrisation utility, and no symmetry CHECK.** Considered and
  rejected with the user. Editing `gmshToFoam` would fork upstream (C1) and
  help only gmsh meshes. A chained utility could DETECT symmetry rather than be
  told -- hash the point set against each candidate mid-plane, snap only if it
  is nearly-but-not-exactly symmetric, so genuinely asymmetric domains fail the
  test and are untouched -- but it is not worth building: refinement already
  cures the problem, point snapping still left **82 of 1592** volume pairs
  differing (face area/centroid sum in FACE-POINT ORDER, which mirroring does
  not preserve), and as the user noted, **an unstructured triangular mesh has no
  mirror symmetry at all**, so a check would fire constantly and mean nothing.

## Two errors of mine, recorded so the pattern is visible

1. **Inferred E from phi without measuring E.** I reported "ePotential is the
   least asymmetric field, so Poisson is not the driver". True and irrelevant:
   E is the GRADIENT, so differentiation amplified a 0.69% `phi` asymmetry into
   a 29% field asymmetry (41x). Worse, my parser silently skipped vector fields,
   so `E` had been in my field list and produced no row -- and I did not notice
   the absence. The user's objection (*"a TINY change in E field distribution
   can cause ... weird stuff that are totally NUMERICAL"*) was correct and I had
   dismissed it on the strength of the wrong quantity.
2. **Verified half the input and called it verified.** My symmetriser rewrote
   only the `internalField` block; my verification measured only the
   `internalField` and reported `0.000e+00`, which I quoted as proof the input
   was symmetric -- while the cathode patch kept an 18% ramp. Four subsequent
   tests (4 tolerances, the snapped mesh, 2 schemes) all returned ~6e-2 because
   they shared that untouched seed. (As it turned out the stored patch values
   are recomputed by the mixed BC anyway, so they were outputs, not inputs --
   but the verification was still incomplete and reported as complete.)

See also `doc/rules-postmortems.md` (A1, A2) -- both are instances of the same
rule: a measurement is not evidence until you have shown it measured what you
think.
