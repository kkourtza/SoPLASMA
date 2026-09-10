# `system/plasmaSimulationControls` — option reference

Solution control: the timestep limiters, the outer coupling loop, and the global Poisson numerics.

**Read TOP-LEVEL, not per region**, with two exceptions noted in place. `constant/` holds physical properties; this file holds how they are solved.

> **50 keys are readable here; a typical case sets ~15.** The rest default. This is a reference, not a checklist.

Part of the [option reference](README.md). Defaults are checked against the
source by `tools/checkConfigReference.py`.


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

## `simulationType`

Whether this is a real time-marching run or a pseudo-steady solve.

**Default:** `transient`

| value | meaning |
|---|---|
| `transient` | Real time-marching. The time axis is physical. What every shipped tutorial does. |
| `pseudoSteady` | Pseudo-transient continuation toward a stationary state. `deltaT` becomes a **relaxation parameter**, the time axis is **not** physical, and convergence means the fields have stopped changing — not that a transient was resolved. |

Set it in `configuration/config` and pass it through as
`simulationType $simulationType;`, like the other keys here.

### Why one key instead of several

The settings `pseudoSteady` implies interlock across three files, and getting
one wrong either aborts the run or — worse — silently destroys its
convergence. Both were measured (see `docs/design/steady-mode-spec.md`):

- `onNonConvergence retryStep` retries by *shortening* `deltaT`, which is
  meaningless once the step is fixed.
- Manual `relaxationFactors` in `fvSolution` took a converging pseudo-steady
  case from **9472/9472 converged steps to 0/10**, because SoPLASMA already
  applies *adaptive* Aitken outer relaxation (`plasmaOuterRelaxation`, observed
  `omega 0.371`) and fixed factors fight it.

So `pseudoSteady` forces `adjustTimeStep false`, **refuses** manual
`relaxationFactors`, notes the time scheme if it is not `backward` (BDF2 was
measured 5.5x faster to the same pseudo-time), and reports all of it at
start-up. Asking for `adjustTimeStep true` alongside it is a fatal error rather
than a silent override: the two requests contradict, and guessing which was
meant is how a run ends up not being the run its author thinks it is.

### `ddtSchemes steadyState` is refused, deliberately

Not offered under either `simulationType`. The `ddt` term **is** the diagonal
of a segregated species-transport equation. Drop it and any species without an
implicit loss term has a matrix row with no diagonal at all, so the linear
solver divides by zero on its first sweep — observed as a `SIGFPE` inside
`GaussSeidelSmoother::smooth`, a stack that names the linear solver and points
nowhere near the cause. `simpleFoam` can drop `ddt` because SIMPLE's pressure
equation and under-relaxation supply the diagonal dominance; nothing here does.

The supported route keeps `ddt`:

```
simulationType    pseudoSteady;   // configuration/config
ddtSchemes        backward;
```

### If you are chasing a steady DC glow

Drive it with `circuit { type currentSource; }`. Near the current-voltage
characteristic's minimum, voltage is not a monotone function of current, so
voltage control is ill-posed (Almeida et al 2016). Measured on the same case:
**6 outer correctors under current control against 150/150 under voltage
control.**

## `dischargeCurrent` — Sato's discharge current

**ON BY DEFAULT since 2026-09-03, and the sub-dictionary is optional.** The
discharge current is the primary measurable of almost every discharge
simulation — the one number an experiment can be compared against — so a case
does not have to ask for it. It costs **one extra Poisson solve at start-up**
and two surface integrals per write.

```
dischargeCurrent
{
    enabled         true;       // DEFAULT true. `false` is the explicit opt-out
    drivenPatch     <name>;     // DERIVED when absent (see below)
    groundedPatches ( <name> ); // DERIVED when absent
    perSpecies      false;      // per-species current columns
    writeInterval   1;          // steps between CSV rows
    printInterval   0;          // steps between log lines; 0 = never
    crossCheck      false;      // compare against the surface-integral current
}
```

### The electrode patches are derived, not restated

The potential's own boundary conditions already say which patch is driven and
which are grounded, so naming them again is a second source of truth that can
disagree. The rule:

| classification | condition |
|---|---|
| **driven** | a Dirichlet `ePotential` condition that is **time-varying** (`uniformFixedValue` with a table/sine/ramp), or a **non-zero constant** |
| **grounded** | a **non-time-varying** Dirichlet condition equal to **zero** |
| not an electrode | everything else — `zeroGradient`, `empty`, the region interfaces, `thinDielectricPotential`, `processor` |

**The time-varying test is the load-bearing part.** At `t = 0` a ramp reads
*exactly zero*, so a value-only rule would classify the driven electrode of
every ramped case as ground — leaving no drive and a singular weighting-field
problem. The discriminator is therefore the **condition type**, not the present
value.

Derivation is reported at start-up:

```
plasmaDischargeCurrent: electrode patches DERIVED from the ePotential boundary conditions
    driven   left   (imposed, time-varying or non-zero)
    grounded 1(right)   (imposed, zero)
```

If the driven electrode is **ambiguous** (none found, or more than one) the run
aborts and lists the candidates, because which electrode the current is measured
at is a physical choice rather than something to guess. Name it explicitly, or
set `enabled false`.

### The ground may be in another region

In a DBD the ground sits *behind* the barrier, on the dielectric mesh. That is
fine: the weighting field is solved **monolithically across every region**, so a
ground on a dielectric mesh is reached correctly — validated to `1.1e-15`
against the analytic series-stack capacitance.

> **SUPERSEDED 2026-09-03.** An earlier error message claimed the single-region
> weighting solve could not reach another region, and advised disabling the
> diagnostic. That has been false since the monolithic multi-region assembly
> landed. What remains true: a **region interface** can never be named as an
> electrode, because the weighting solve forces grounded patches to `fixedValue`
> and that would destroy the interface coupling.

### What it writes

`postProcessing/dischargeCurrent/current.csv`, with a header recording `C_g` and
the wedge revolution factor, and columns
`time, V_applied, I_total, I_cond, I_disp` (plus per-species columns under
`perSpecies true`).


## `limitSpeciesCo` / `printSpeciesCo` / `maxSpeciesCo`

Species Courant number — convective and diffusive. `maxSpeciesCo` feeds **both**
`maxSpeciesConvectiveCo` and `maxSpeciesDiffusiveCo`.

**Defaults:** `limitSpeciesCo` `false`, `printSpeciesCo` `false`,
`maxSpeciesConvectiveCo` `1.0`, `maxSpeciesDiffusiveCo` `1.0`

`100` on this case is a **rail, not the control** — a ceiling that stops the
step running away if every other limiter releases, and measured never to bind.
Set the two underlying keys separately in `plasmaSimulationControls` if you need
convective and diffusive caps to differ.

## `limitChemistryCo` / `printChemistryCo` / `maxChemistryCo`

Chemistry Courant number.

**Defaults:** `false`, `false`, `1.0`

**Not built from an effective ionisation coefficient.** It uses the fractional
rate of change of state over all species from the full mechanism's production
and loss, because a *net* rate cancels near the critical field and would report
a long timescale exactly where the chemistry is fastest.

## `limitDielectricRelaxationRatio` / `printDielectricRelaxationRatio` / `maxDielectricRelaxationRatio`

Ratio of the step to the dielectric relaxation time `ε/σ`.

**Defaults:** `false`, `false`, `1.0`

> On this case the gate is `false`, which matches the default and so is
> redundant. With the gate off, `maxDielectricRelaxationRatio` only reaches the
> **printed report** as the reference value the measured ratio is shown against
> (`plasmaTimeControl.C:1433`); the `deltaT` computation that uses it is
> unreachable. Keep it only if you want the report compared against `5` rather
> than `1`.

## `maxVoltageRisePerStep`

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

