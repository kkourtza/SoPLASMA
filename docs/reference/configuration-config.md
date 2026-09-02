# `configuration/config` — the case variable reference

Every dictionary in a SoPLASMA case begins with

```
#include "../configuration/config"
```

and refers to values as `$name`. **This is a case convention, not an OpenFOAM
feature**: a value defined in `config` that no dictionary references does
nothing, and nothing warns you. Check it after any edit:

```bash
for k in $(grep -oE '^[a-zA-Z][a-zA-Z0-9_]*' configuration/config | sort -u); do
    grep -rqF "\$$k" system/ constant/ etc/ 0.orig/ 2>/dev/null || echo "DANGLING: $k"
done
```

## What is here and what is in the case

This document is the **reference**: what each option is, its type, its default,
and the menu of alternatives. That text is identical in every case, so it lives
in one place.

The **rationale** — why a particular case chose a particular value, and what was
measured there — stays in that case's `configuration/config`. The two are
different things with different lifecycles, and copying the reference into each
case is how it goes stale.

> **Defaults in this document are checked against the source** by
> `tools/checkConfigReference.py`, which fails if a documented option has no
> reader, if a stated default no longer matches the code, or if the code gains
> an option this document does not mention.

---

## Applied voltage

### `voltageRamp`

The driven electrode's potential as a function of time, in volts. Consumed by
`etc/changeDictionary.<gasRegion>` as the `uniformValue` of a
`uniformFixedValue` boundary condition on the driven patch.

**REQUIRED — no default.** There is no sensible one: a discharge with no drive
does nothing, and guessing an amplitude would be inventing the physics.

**Type:** any OpenFOAM `Function1`. The full set in v2412
(`OpenFOAM/primitives/functions/Function1/`):

```
constant  table  tableFile  csv  polynomial  sine  cosine  square  step
linearRamp  quadraticRamp  halfCosineRamp  quarterCosineRamp  quarterSineRamp
ramp  scale  uniform  zero  one  none  coded  inputValueMapper
functionObjectTrigger  functionObjectValue
```

Because the whole waveform is one value, the *type* can be changed from
`config`. An `appliedVoltage` + `riseTime` pair could not do that: `sine {...}`
has no amplitude-and-rise-time to split in two, so such a pair quietly limits a
case to ramps.

#### Examples

A linear ramp to 8 kV over 100 ns, held after. `table` clamps to its last value
beyond the final entry, so this is a ramp-then-hold, not a ramp-then-zero:

```
voltageRamp   table ((0 0) (100e-9 8e3));
```

Constant DC — a Townsend or steady-state study:

```
voltageRamp   constant 8e3;
```

AC, the usual DBD drive. 1 kHz, 20 kV amplitude, i.e. 40 kV peak-to-peak:

```
voltageRamp   sine { frequency 1e3; amplitude 20e3; level 0; };
```

Damped or growing AC — `scale` multiplies one `Function1` by another, so an
envelope is just a second function:

```
voltageRamp   scale
{
    value  sine { frequency 1e3; amplitude 20e3; level 0; };
    scale  linearRamp { start 0; duration 1e-3; };
};
```

Bipolar pulses:

```
voltageRamp   square { frequency 1e4; amplitude 8e3; level 0; markSpace 1; };
```

A measured waveform from disk — how you drive from an experiment:

```
voltageRamp   csv
{
    nHeaderLine       1;
    refColumn         0;        // time
    componentColumns  (1);      // volts
    separator         ",";
    mergeSeparators   no;
    file              "constant/measuredVoltage.csv";
};
```

A ramp with no derivative discontinuity at either end, gentler on the
space-charge/Poisson coupling than `linearRamp`:

```
voltageRamp   halfCosineRamp { start 0; duration 100e-9; };
```

#### Interactions

- **`maxVoltageRisePerStep`** limits `deltaT` so the applied voltage changes by
  at most that much per step. It binds during a ramp and is inert at constant
  voltage, where `dV/step ≈ 0`. It needs `voltagePatchName` to name the driven
  electrode, or it turns itself off and says so.
- The value reaches the field file **unexpanded**. `changeDictionary` expands
  nothing; `$voltageRamp` is resolved when the field under `0/` is parsed, which
  works because that file carries its own `#include`. Check it without running
  the solver:
  ```bash
  foamDictionary 0/<region>/ePotential \
      -entry boundaryField/<patch>/uniformValue
  ```

---

## Gas state

The gas is described by **two** numbers and everything else is derived. There is
no third place to state a pressure and no second unit.

### `gasPressure`

Absolute gas pressure, **in Pascal**. Feeds `backgroundGas/pressure` in
`plasmaSpeciesProperties`.

**Default:** `101325`

Pa is the only input unit anywhere in a case. Anything the solver *prints* shows
Pa first with the conversions after, so a value can be checked against
literature quoted in Torr without converting by hand:

```
100000 Pa (0.986923 atm, 750.062 Torr, 1 bar)
```

#### What derives from it

| derived | where | why it is not a separate setting |
|---|---|---|
| gas number density `N` | `plasmaSpecies`, from `p` and `T` by the ideal gas law | a second definition of `N` is how `E/N` came to be evaluated at 1 atm in a case running at 1 bar |
| the Boltzmann sweep's tabulation density | `plasmaBoltzmann` | three-body processes break `E/N` similarity, so their contribution scales with `N` and the tables are valid only near the density they were solved at. There is therefore exactly **one** correct density for the sweep — the case's own |
| `pGasAtm_` in `plasmaEnergy` and `plasmaTransport` | both | each divides by 101325 itself |

> **`boltzmann/pressureAtm` is REJECTED.** It was a second spelling of this
> quantity, in atm, in a different file. The error prints the derived value. And
> `backgroundGas/energy/pressure` — a nested duplicate of the outer key in the
> same file — is gone; its two readers now follow the owner.

#### Alternative: state the density directly

`backgroundGas` accepts `numberDensity` instead of `pressure` + `T`, in which
case the ideal-gas closure is skipped. Use it to match a reference that quotes
`N` rather than `p`, e.g. the streamer benchmark's `2.414e25 m^-3`.

#### Interactions

- **Changing it invalidates the EEDF tables.** The sweep stamp records
  `pPa=<value>`, so with `generateTables yes` the tables re-solve automatically;
  with `generateTables no` the run **stops** and names the differing field.
- `1e5` is **1 bar, not 1 atm**. The two differ by 1.3% in density, which shifts
  `E/N` and therefore every rate coefficient.

### `gasTemperature`

Gas temperature [K]. Feeds `backgroundGas/energy/T`.

**Default:** `300`

Used for the ideal-gas closure of `N`, and as the sweep's `T_gas`. It is the
*initial* and, unless gas heating is solved, the *fixed* temperature.

#### Interactions

- **Orthogonal to gas heating.** `backgroundGas/energy/solve` (**default:**
  `false`) decides whether a temperature equation is solved at all; this is the
  value it starts from, or holds at when it is not solved.
- Also stamped into the sweep, so changing it re-solves the tables.
- `backgroundGas/energy/kappa` (**default:** `0.026` W/m/K) is the gas thermal
  conductivity, read only when `solve true`.

---

## Region coupling

### `monolithicRegionCoupling`

Whether the gas and every dielectric are assembled into **one** matrix and
solved together. Feeds the `useImplicit` entry of the `coupledElectricPotential`
boundary condition on each interface patch.

**Default:** `false` — i.e. segregated, if the entry is omitted.

`true` satisfies the permittivity jump and the surface charge exactly at every
linear solve. `false` lags the interface: each region sees the other's potential
from the previous corrector.

> **Why the case variable is not called `useImplicit`.** The *dictionary key*
> is, because `useImplicit_` is a member of OpenFOAM's own `fvPatchField` with a
> `useImplicit()` accessor — renaming it would fight the framework. But
> "useImplicit" says nothing about *what* is made implicit, and this solver has
> several independent implicit/explicit choices: the chemistry source
> (`solver implicitRate`), the wall flux (`ddWallFluxImplicit`) and the Poisson
> scheme (`semiImplicit`). The case variable is ours, so it names the thing.

#### Interactions

- **Needs an assembly-aware linear solver.** With `GAMG` you must also set
  `agglomerator assembledFaceAreaPair`, or the run aborts with
  `Attempt to cast type lduPrimitiveMeshAssembly to type fvMesh` — an internal-
  looking message for a configuration error.
- **It changes the achievable outer tolerance.** Measured on the needle-DBD
  geometry: segregated coupling plateaus the outer `ePotential` residual near
  `1e-7`, so at `outerCoupling/tolerance 1e-8` the loop hits its cap and
  `retryStep` discards every step. Monolithic removes the plateau.
- Measured cost of the lag, on the two-region analytic bed: interface potential
  `8.8e-5` relative error segregated versus `2.0e-6` monolithic, the latter
  being the write-precision floor.

---

## Initial condition

### `floorDensity`

Lower bound applied to every species number density [m⁻³]. Feeds
`minNumberDensity` in `plasmaSpeciesProperties`.

**Default:** `0.0` — no floor.

**In a case with no seeded blob this is the initial condition, not a numerical
guard.** The needle-DBD discharge initiates from this background density in the
geometrically enhanced field at the tip, so the value is physics and changing it
changes when the discharge starts.

#### The consequence that is easy to miss

**A clamp is a source.** Electrons are injected wherever attachment would have
removed them. That only conserves charge if the positive-ion floors sum to the
electron floor, split by the mole fractions of the gas that produces them — so
raising one means raising the others. Otherwise the undisturbed background
carries net space charge, which drives a Poisson field out of nothing.

---

## Secondary electron emission

### `seeCoefficient`

> **This variable is being retired.** SEE is a **per-surface** property: the key
> the solver reads is `defaultSEEC`, a boundary-condition entry on the species
> patch fields (`n_e` and `nEps_e`). It belongs in
> `etc/changeDictionary.<region>`, per patch, not as one global.
>
> A single case variable forces every emitting surface to share a value, which
> is wrong for the common geometry: a metal needle and a dielectric barrier are
> different materials with different γ. And **dielectrics emit too** — SEE is
> not an electrode-only phenomenon.

**Default (of `defaultSEEC`):** `0.001`

γ is the probability that an ion striking a surface releases an electron. In a
DBD it is what makes the discharge self-sustaining rather than a single
avalanche, so it is load-bearing rather than a detail.

> **Default changed 2026-09-02, from `0.05` to `0.001`.** γ spans roughly 10⁻³
> for contaminated oxides and dielectric barriers up to 10⁻¹ for clean metals in
> vacuum. `0.05` is a clean-metal figure and the wrong end of that range to
> default to, because the surfaces this solver is aimed at — air-exposed
> electrodes and dielectric barriers — are contaminated. A default fifty times
> too high manufactures self-sustainment the case never asked for, and the run
> looks entirely healthy.
>
> A case that states its own value is unaffected; one that relied on the old
> default should state one. It is still worth recording the material and a
> citation next to whatever value you choose.

---

## Time control

### `adjustTimeStep`

Whether `deltaT` is chosen adaptively from the limiters below.

**Default:** `false`

With `false` the step stays at `deltaT` and every limiter becomes advisory —
they still *report* if their `print*` switch is on, but nothing acts on them. A
plasma case essentially always wants `true`: the stable step spans orders of
magnitude between the quiescent ramp and a propagating front.

### `deltaT`, `endTime`, `writeControl`, `writeInterval`

OpenFOAM's own `controlDict` entries; documented upstream. Two notes specific to
this solver:

- **`deltaT` is the *initial* step only** when `adjustTimeStep true`. It is
  additionally capped on a fresh start by `maxInitialDeltaT` (**default:**
  `1e-12`), which does not apply on a restart.
- **`writeControl adjustableRunTime`, not `runTime`.** With plain `runTime` the
  writes land on whichever step straddles the interval, so two runs of the same
  case write at slightly different times — measured at 7.4e-13 s apart, which
  was 0.28 of an energy relaxation time and enough to make a comparison
  meaningless.

#### Sizing `endTime`

Size it to when the discriminating observable appears, not to a round number.
**Measured 2026-09-02 on this case:** ~2.8 s/step, so `2e-08` is roughly a
**three-hour** run. Initiation appears at ~2.1 ns (E/N reaches the ionisation
range at V ≈ 167 V under ~15× tip enhancement), so `3e-09` crosses it with
margin in well under an hour. `2e-08` is the setting for the full self-limiting
physics, since barrier charging happens over tens of ns.

---

## Timestep limiters

Each limiter has three keys with the same shape: a `limit*` gate, a `print*`
reporter, and a `max*` value. **`print*` works independently of `limit*`**, so a
limiter can be measured before it is switched on — which is the right way to
choose the value.

**Which one actually bound is reported, and you should read it:**

```bash
grep -oE "deltaT set by: +.*" log.soPlasmaFoam | sed 's/deltaT set by: *//' \
    | sort | uniq -c | sort -rn
```

Measured on this case, 322 steps under LMEA:

```
197  energy relaxation      <- the control, and the LMEA signature
 99  growth cap (1.2x)
 19  rejection memory       <- retryStep active, not smooth by luck
  4  Co_chem
  2  coupling margin (backoff)
```

Note what is *absent*: `maxSpeciesCo` never bound once.

### `limitSpeciesCo` / `printSpeciesCo` / `maxSpeciesCo`

Species Courant number — convective and diffusive. `maxSpeciesCo` feeds **both**
`maxSpeciesConvectiveCo` and `maxSpeciesDiffusiveCo`.

**Defaults:** `limitSpeciesCo` `false`, `printSpeciesCo` `false`,
`maxSpeciesConvectiveCo` `1.0`, `maxSpeciesDiffusiveCo` `1.0`

`100` on this case is a **rail, not the control** — a ceiling that stops the
step running away if every other limiter releases, and measured never to bind.
Set the two underlying keys separately in `plasmaSimulationControls` if you need
convective and diffusive caps to differ.

### `limitChemistryCo` / `printChemistryCo` / `maxChemistryCo`

Chemistry Courant number.

**Defaults:** `false`, `false`, `1.0`

**Not built from an effective ionisation coefficient.** It uses the fractional
rate of change of state over all species from the full mechanism's production
and loss, because a *net* rate cancels near the critical field and would report
a long timescale exactly where the chemistry is fastest.

### `limitDielectricRelaxationRatio` / `printDielectricRelaxationRatio` / `maxDielectricRelaxationRatio`

Ratio of the step to the dielectric relaxation time `ε/σ`.

**Defaults:** `false`, `false`, `1.0`

> On this case the gate is `false`, which matches the default and so is
> redundant. With the gate off, `maxDielectricRelaxationRatio` only reaches the
> **printed report** as the reference value the measured ratio is shown against
> (`plasmaTimeControl.C:1433`); the `deltaT` computation that uses it is
> unreachable. Keep it only if you want the report compared against `5` rather
> than `1`.

### `maxVoltageRisePerStep`

Limits `deltaT` so the applied voltage changes by at most this much per step [V].

**Default:** `100.0` — so stating `100` is redundant.

`0` disables. Always active once `voltagePatchName` (**default:** empty) names
the driven electrode; with no patch named it turns itself off and says so, which
is deliberate — it is on by default, and a default that killed every
constant-voltage case would be worse than the problem it solves. At constant
voltage `dV/step ≈ 0` and it never binds.

> **Four keys once described this one quantity.**
> `limitVoltageRiseRate` and `maxVoltageRiseRate` are **rejected** with an error
> carrying the exact translation for your dictionary, because the old pair gated
> on the *switch*: `maxVoltageRiseRate 1` with the switch false was inert, so
> copying the number across would have set 1 V/step — on one measured case a
> ~2700× change in the driving term. `printVoltageRiseRate` prints a notice
> saying it is no longer read; reporting now follows the limiter.

---

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

### `ddtSchemes`

Time discretisation for every transported field. `backward` is second-order
implicit.

Second order in *time* is necessary but not sufficient: with a single outer
corrector the scheme is **first** order whatever this says, because no amount of
`ddt` accuracy repairs the Poisson↔species lag. Measured on the shipped
benchmark: p = 0.87 at one corrector, p = 1.94 converged. See
`outerCoupling/target` in `plasmaSimulationControls`.

### `gradSchemes`

Gradient scheme. `leastSquares` rather than `Gauss linear` — it is more accurate
on the non-orthogonal cells a needle mesh has near the tip.

### `driftDiffusionFluxScheme`

How the drift and diffusion fluxes are assembled. Read by
`plasmaSpeciesProperties`, not `fvSchemes`. `standard` is the plain
Scharfetter–Gummel-style form.

### `electronDriftDivScheme`

Divergence scheme for the electron drift flux, `div(phi_e,n_e)`.

`Gauss ROUNDF` is a **third-party** limited scheme and requires
`"libROUNDSchemes.so"` in `solverLibs`. It is chosen for a bounded,
low-dissipation front: an unlimited scheme oscillates across a steep electron
front and a first-order upwind smears it.

### `electronDiffusionLaplacianScheme`

`laplacian(D_e,n_e)`. `Gauss harmonic corrected 0.33` — `harmonic` on the
diffusivity for the reason below, `corrected 0.33` for mesh non-orthogonality.

### `electronMobilityInterpolationScheme` / `electronDiffusivityInterpolationScheme`

Face interpolation of `mu_e` and `D_e`, and of the ion equivalents.

`harmonic`, not `linear`, and this is a correctness point rather than a
preference: these coefficients span **orders of magnitude** across a discharge
front, and linear interpolation of a steep coefficient can put the face value
*above both* adjacent cell values. The harmonic mean cannot.

### `ePotentialLaplacianScheme` / `ePotentialSnGradScheme`

The Poisson operator and its surface-normal gradient. `Gauss linear corrected
0.33` and `corrected 0.33`.

Both `laplacian(epsilon,ePotential)` and
`laplacian((epsilon+(deltaT*electricalConductivity)),ePotential)` must be listed
in `fvSchemes` — the second is the `semiImplicit` operator. Listing only one
works until someone changes `poissonScheme`.

### `snGradSpeciesScheme`

Surface-normal gradient for species, `snGrad(n_e)`. `corrected 0.33` — the
non-orthogonality correction the needle mesh needs.

---

## Linear solvers

### `ePotentialSolver`

The linear solver for `ePotential`. **Default:** none — `fvSolution` requires it.

`GAMG` | `PBiCGStab` | `PCG` | `petsc`. One `ePotential` block carries the
settings for all of them, because **an OpenFOAM solver dictionary silently
ignores entries it does not use** (measured: a `PCG` block carrying
`agglomerator`, `nCellsInCoarsestLevel` and a `petsc { options { ... } }`
sub-dictionary runs clean). So switching backend is this one word.

#### Interactions

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

### `solverLibs`

Libraries `controlDict` loads, as one whole-value list.

**Default:** none — `controlDict` requires the entry.

`"libROUNDSchemes.so"` is required by `electronDriftDivScheme`'s `Gauss ROUNDF`.
Add `"libpetscFoam.so"` for `ePotentialSolver petsc`.

> **PETSc is installed but not loadable on this machine:**
> `libpetsc.so.3.24: cannot open shared object file`. The wrapper is present,
> PETSc itself is not. It is a warning rather than an error, so an unnecessary
> entry costs a confusing line at every start-up.

### `poissonScheme`

**Default:** `semiImplicit`

`explicit` treats the space charge as a pure source; `semiImplicit` also brings
the conductivity in implicitly, which stabilises the space-charge/Poisson
coupling.

**A scheme, not a linear solver** — that is `ePotentialSolver`. The key was once
`PoissonScheme` fed by a case variable called `poissonSolver`: two names for one
setting, each suggesting the other's meaning.

### `nOuterCorrectors`

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

---

## All 37 variables are documented.

Verify the reference against the source and the case with:

```bash
tools/checkConfigReference.py --case <caseDir>
```
