# `system/fvSchemes` and `system/fvSolution` — option reference

Discretisation schemes, linear solvers and their tolerances.

Part of the [option reference](README.md). Defaults are checked against the
source by `tools/checkConfigReference.py`.

## Discretisation schemes

**None of these has a default.** `fvSchemes` entries are named **per species**,
and a missing one is a fatal lookup rather than a fall-through to `default`:

```
Entry 'div(phi_N2p,n_N2p)'     not found in divSchemes
Entry 'laplacian(D_N2p,n_N2p)' not found in laplacianSchemes
Entry 'interpolate(mu_N2p)'    not found in interpolationSchemes
```

They fail one *family* at a time, minutes apart, which is why a case copied from
a bed with immobile ions appears to work and then dies three times.

> The convective flux is `phi_<species>`, **not** `particleFlux_<species>` —
> that is the name of the registered diagnostic *field* the boundary conditions
> read. Guessing it costs a run.

## `ddtSchemes`

Time discretisation for every transported field. `backward` is second-order
implicit.

Second order in *time* is necessary but not sufficient: with a single outer
corrector the scheme is **first** order whatever this says, because no amount of
`ddt` accuracy repairs the Poisson↔species lag. Measured on the shipped
benchmark: p = 0.87 at one corrector, p = 1.94 converged. See
`outerCoupling/target` in `plasmaSimulationControls`.

## `gradSchemes`

Gradient scheme. `leastSquares` rather than `Gauss linear` — it is more accurate
on the non-orthogonal cells a needle mesh has near the tip.

## `driftDiffusionFluxScheme`

How the drift and diffusion fluxes are assembled. Read by
`plasmaSpeciesProperties`, not `fvSchemes`. `standard` is the plain
Scharfetter–Gummel-style form.

## `electronDriftDivScheme`

Divergence scheme for the electron drift flux, `div(phi_e,n_e)`.

`Gauss ROUNDF` is a **third-party** limited scheme and requires
`"libROUNDSchemes.so"` in `solverLibs`. It is chosen for a bounded,
low-dissipation front: an unlimited scheme oscillates across a steep electron
front and a first-order upwind smears it.

## `electronDiffusionLaplacianScheme`

`laplacian(D_e,n_e)`. `Gauss harmonic corrected 0.33` — `harmonic` on the
diffusivity for the reason below, `corrected 0.33` for mesh non-orthogonality.

## `electronMobilityInterpolationScheme` / `electronDiffusivityInterpolationScheme`

Face interpolation of `mu_e` and `D_e`, and of the ion equivalents.

`harmonic`, not `linear`, and this is a correctness point rather than a
preference: these coefficients span **orders of magnitude** across a discharge
front, and linear interpolation of a steep coefficient can put the face value
*above both* adjacent cell values. The harmonic mean cannot.

## `ePotentialLaplacianScheme` / `ePotentialSnGradScheme`

The Poisson operator and its surface-normal gradient. `Gauss linear corrected
0.33` and `corrected 0.33`.

Both `laplacian(epsilon,ePotential)` and
`laplacian((epsilon+(deltaT*electricalConductivity)),ePotential)` must be listed
in `fvSchemes` — the second is the `semiImplicit` operator. Listing only one
works until someone changes `poissonScheme`.

## `snGradSpeciesScheme`

Surface-normal gradient for species, `snGrad(n_e)`. `corrected 0.33` — the
non-orthogonality correction the needle mesh needs.

---

## Linear solvers

## `ePotentialSolver`

The linear solver for `ePotential`. **Default:** none — `fvSolution` requires it.

`GAMG` | `PBiCGStab` | `PCG` | `petsc`. One `ePotential` block carries the
settings for all of them, because **an OpenFOAM solver dictionary silently
ignores entries it does not use** (measured: a `PCG` block carrying
`agglomerator`, `nCellsInCoarsestLevel` and a `petsc { options { ... } }`
sub-dictionary runs clean). So switching backend is this one word.

### Interactions

- **Under `monolithicRegionCoupling true`, GAMG needs
  `agglomerator assembledFaceAreaPair`.** The default agglomerator casts the
  assembled matrix to `fvMesh` and aborts with
  `Attempt to cast type lduPrimitiveMeshAssembly to type fvMesh`. Measured on
  the two-region bed, per timestep: `GAMG + assembledFaceAreaPair` 11
  iterations, `PCG/DIC` 114, `smoothSolver/GaussSeidel` 2000 (its cap, **never
  converging**, and 7.3% wrong).
- **Symmetry follows `poissonScheme`.** Under `explicit` the matrix is a pure
  Laplacian and symmetric, so `PCG`/`DIC` applies. Under `semiImplicit` the
  conductivity enters the operator and it is **not** symmetric — use
  `PBiCGStab`/`DILU`, not `PCG`/`DIC`.
- Species matrices are drift-dominated and strongly asymmetric: `PBiCGStab` with
  **`DILU`**, not `DIC`. Measured: 73 species solves hit the 2000-iteration cap
  without converging in the first 5 steps under Gauss-Seidel smoothing.

## `solverLibs`

Libraries `controlDict` loads, as one whole-value list.

**Default:** none — `controlDict` requires the entry.

`"libROUNDSchemes.so"` is required by `electronDriftDivScheme`'s `Gauss ROUNDF`.
Add `"libpetscFoam.so"` for `ePotentialSolver petsc`.

> **PETSc is installed but not loadable on this machine:**
> `libpetsc.so.3.24: cannot open shared object file`. The wrapper is present,
> PETSc itself is not. It is a warning rather than an error, so an unnecessary
> entry costs a confusing line at every start-up.

## `poissonScheme`

**Default:** `semiImplicit`

`explicit` treats the space charge as a pure source; `semiImplicit` also brings
the conductivity in implicitly, which stabilises the space-charge/Poisson
coupling.

**A scheme, not a linear solver** — that is `ePotentialSolver`. The key was once
`PoissonScheme` fed by a case variable called `poissonSolver`: two names for one
setting, each suggesting the other's meaning.

## `nOuterCorrectors`

PIMPLE outer correctors per timestep. **Default:** none — `fvSolution` requires
it.

This is a *cap*, not the convergence test. What decides whether the outer loop
has converged is `outerCoupling/tolerance` in `plasmaSimulationControls`, and
editing the per-field tolerances to fix an outer-loop failure changes nothing.

---

## Tolerances

Five values, having been ten names (see the history in the case file).

| variable | value | scope |
|---|---|---|
| `linearTol` | `1e-10` | linear tolerance for `ePotential`, `n_e` and every species |
| `linearRelTol` | `0.0` | relative tolerance, all fields — `0` means "converge on the absolute tolerance", the same decision for every field here |
| `TgasTolerance` | `1e-12` | tighter, because the gas energy equation is stiff |
| `ePotentialTolerancePIMPLE` | `1e-8` | the outer gate for `ePotential` |
| `outerTol` | `1e-10` | the outer gate for `n_e` and species |

**Defaults:** none. `fvSolution` requires every one; there is no fall-through.

> **An outer tolerance must stay comfortably above the linear one — at least
> 10×.** The solver refuses otherwise. An outer tolerance sitting at the inner
> solver's own convergence floor can never be met, so the loop silently runs to
> `maxCorrectors` every step: a fixed iteration count wearing the name of a
> convergence test. `ePotentialTolerancePIMPLE 1e-8` against `linearTol 1e-10`
> is the 100× margin that keeps it honest.

