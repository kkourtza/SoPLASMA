# Mechanism-driven plasma chemistry for SoPLASMA

Branch: `feature/chemistry_integration_and_BoltzmannSolver` (27 commits, fast-forward from `main`)

## What this adds

A case declares **one mechanism file**; the species list, their charges and masses, the electron
transport coefficients, the electron-impact rate coefficients and the heavy chemistry are all
derived from it. Nothing is re-typed between codes, and every generated artifact carries a hash of
the master mechanism, so a mismatched pair is refused rather than silently used.

| library | role |
|---|---|
| `plasmaBoltzmann` | solves the EEDF in process at start-up; staleness-checked against the mechanism hash |
| `plasmaChemistry` | the whole reaction set as an `ODESystem`; native and optional Cantera heavy backends |
| `plasmaSpecies` | derives `activeSpecies`, charges and masses from the mechanism |
| `plasmaTransport` | assembles the chemistry source into the species equations |

## Why `ODESystem` and not `chemistryModel`

The `chemistryModel`/`reactingFoam` stack is built on `psiReactionThermo` — mass fractions and one
temperature — where this needs number densities and a separate electron temperature; every
`ReactionRate` in it is `k(T)` Arrhenius where electron-impact rates are `k(E/N)`; and its species
are bound to a thermo package that would have to be fabricated for `N2(A3)`. The valuable part, the
stiff integrators, is separable: `ODESystem` is a bare abstract class with no thermo dependency.

## Numerics

* The source is `dn_s/dt = P_s − L_s·n_s`, loss implicit via `fvm::Sp` — positivity at any timestep,
  and the outer Picard iteration is a contraction on the loss term.
* `solver adaptive` chooses **per cell** between linearising and stiff integration on the
  local `L·dt`, because the stiffness is local: a streamer head is stiff while the bulk is not.
* Temporal order measured by Richardson extrapolation: **p = 1.68**, charge conserved to 2e-16.

## Validation

| check | result |
|---|---|
| 0-D chemistry vs Cantera | 2.8e-03 |
| native vs Cantera heavy rates | 1.4e-15 |
| native vs Cantera Jacobian | 6.7e-16 |
| charge residual of the reaction source | ~1e-15 every step |
| positive streamer benchmark | propagates, 12 species, 48 reactions |

`applications/test/testChemistryBackends` runs mesh-free in a second and cross-validates the two
backends on any mechanism.

## Compatibility

* **Cantera is optional.** Without `CANTERA_DIR` everything builds; selecting `chemistryBackend
  cantera` then fails with a message naming the variable. Cantera 3.x and 4.x both supported.
* **Existing cases are unaffected**: `reactions electronImpact` + `solver explicitSource` is the legacy path and remains available under those names.
* Requires the SoEEDF library (formerly BoltzmannSolver) for the EEDF sweep, located via
  `BOLTZMANN_DIR`.

## Not included

Gas heating — the gas temperature is a parameter. That is the next piece of work.
