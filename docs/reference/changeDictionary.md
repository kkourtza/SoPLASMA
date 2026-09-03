# `etc/changeDictionary.<region>` — option reference

The boundary conditions, one file per mesh region. `Allrun` copies each to
`system/<region>/changeDictionaryDict` and runs `changeDictionary`, which writes
the entries into the fields under `0/<region>/`.

Part of the [option reference](README.md).

> **Most keys here are OpenFOAM's own** (`type`, `value`, `uniformValue`,
> `inletValue`, `phi`) and are documented upstream. This file covers the ones
> specific to this solver, and the ones whose *use* here is specific.

## Two things that are easy to get wrong

**`changeDictionary` expands nothing.** A `$var` — nested or whole-value — is
copied through verbatim and resolved later, when the field under `0/` is parsed.
That works because the field carries its own `#include`, inherited from
`0.orig/`. So `Illegal dictionary entry or environment variable name <x>` means
a field file that has **lost** its include, not that nesting is unsupported.
Check without running the solver:

```bash
foamDictionary 0/<region>/ePotential \
    -entry boundaryField/<patch>/uniformValue
```

**In a field file the `#include` goes AFTER the `FoamFile` block.** Before it
gives `problem while reading header for object ePotential`.

---
## The driven electrode

### `uniformValue` on a driven electrode

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

## Surface charge, and which side owns it

`coupledElectricPotential` takes `surfCharge` and `surfChargeNbr`. **Exactly one
side of an interface owns the charge, or it is double-counted:**

```
etc/changeDictionary.gas          surfCharge  surfCharge;   surfChargeNbr none;
etc/changeDictionary.dielectric   surfCharge  none;         surfChargeNbr surfCharge;
```

The field name is not a choice — the solver registers the literal `surfCharge` —
and any other value silently means "no surface charge", because the BC does
`if (surfChargeName_ != "none")` and quietly finds nothing.

The gas owns it because the plasma is what deposits it.

## Secondary electron emission

### `defaultSEEC`

γ, the probability that an ion striking a surface releases an electron. Set on
the **species** patch fields (`n_e` and `nEps_e`), per patch.

**Default:** `0.001`

Per surface, not global: a metal needle and a dielectric barrier are different
materials with different γ, and **dielectrics emit too** — SEE is not an
electrode-only phenomenon.

γ spans roughly 10⁻³ for contaminated oxides and dielectric barriers up to 10⁻¹
for clean metals in vacuum. The default is the contaminated end, because that is
what an air-exposed electrode or a barrier actually is. It is load-bearing
rather than a detail: in a barrier discharge γ is what makes the discharge
self-sustaining rather than a single avalanche, so record the material and a
citation beside whatever value you choose.

> **It is checked on a `farField` interface.** A `plasmaWallBC` with
> `enableSurfaceCharging true` there is fatal: a fictitious air region is not a
> barrier, and depositing charge on its interface would invent a dielectric
> surface in the middle of the gas.
