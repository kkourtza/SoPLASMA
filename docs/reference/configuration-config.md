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

## Moved out of this file

These were `config` variables and are now literals in the dictionary that reads
each one, because each was used exactly **once** — see
[the rule](README.md#what-goes-in-configurationconfig-and-what-does-not).

| was | now lives in | reference |
|---|---|---|
| `endTime` `deltaT` `writeControl` `writeInterval` | `system/controlDict` | [controlDict.md](controlDict.md) |
| the `limit*`/`print*`/`max*` limiter keys | `system/plasmaSimulationControls` | [plasmaSimulationControls.md](plasmaSimulationControls.md) |
| `poissonScheme` `nOuterCorrectors` `solverLibs` `TgasTolerance` `ePotentialTolerancePIMPLE` `gradSchemes` `ddtSchemes` `snGradSpeciesScheme` `ePotentialSnGradScheme` | `system/fvSchemes`, `system/fvSolution` | [fvSchemes-fvSolution.md](fvSchemes-fvSolution.md) |
| `floorDensity` → `minNumberDensity` | `constant/plasmaSpeciesProperties` | pending |
| `seeCoefficient` → `defaultSEEC` | `etc/changeDictionary.<region>` | pending |

## The 13 that remain, and why

**Used more than once** — one decision, N places:

| variable | uses |
|---|---|
| `electronDriftDivScheme` `electronDiffusionLaplacianScheme` `electronMobilityInterpolationScheme` `electronDiffusivityInterpolationScheme` `driftDiffusionFluxScheme` | **5 each** — one per mobile species |
| `linearTol` `linearRelTol` `ePotentialSolver` | 3 each |
| `outerTol` `maxSpeciesCo` `ePotentialLaplacianScheme` | 2 each |

**Headline physics inputs**, kept here even where used once because they are what
changes between runs: `voltageRamp`, `gasPressure`, `gasTemperature`,
`monolithicRegionCoupling`.
