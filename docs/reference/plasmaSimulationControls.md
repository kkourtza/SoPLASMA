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

