# `system/controlDict` — option reference

Run control: when the run starts and stops, and when it writes.

Part of the [option reference](README.md). Defaults are checked against the
source by `tools/checkConfigReference.py`.


## `adjustTimeStep`

Whether `deltaT` is chosen adaptively from the limiters below.

**Default:** `false`

With `false` the step stays at `deltaT` and every limiter becomes advisory —
they still *report* if their `print*` switch is on, but nothing acts on them. A
plasma case essentially always wants `true`: the stable step spans orders of
magnitude between the quiescent ramp and a propagating front.

The exception is a pseudo-steady solve, where a fixed step is the *point* — see
[`simulationType`](plasmaSimulationControls.md#simulationtype), which forces
this off and makes `adjustTimeStep true` a fatal error rather than a silent
override.

## `deltaT`, `endTime`, `writeControl`, `writeInterval`

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

### Sizing `endTime`

Size it to when the discriminating observable appears, not to a round number.
**Measured 2026-09-02 on this case:** ~2.8 s/step, so `2e-08` is roughly a
**three-hour** run. Initiation appears at ~2.1 ns (E/N reaches the ionisation
range at V ≈ 167 V under ~15× tip enhancement), so `3e-09` crosses it with
margin in well under an hour. `2e-08` is the setting for the full self-limiting
physics, since barrier charging happens over tens of ns.

