# The applied-voltage boundary condition

## Where the number lives

One place: `configuration/config`.

```
appliedVoltage      18.75e3;     // [V]
voltageRampTime     9.375e-11;   // [s]
```

and the electrode block references them:

```
active_electrode
{
    type            uniformFixedValue;
    uniformValue    table ((0 0) ($voltageRampTime $appliedVoltage));
    value           uniform 0;
}
```

> **`configuration/config` is a case convention, not a framework feature.** Values
> are defined there and pulled into real dictionaries with `$name`. A value that
> no dictionary (or run script) references **does nothing** — it is inert text.
> Before this was wired, `appliedVoltage` was defined in three beds and read by
> none of them: the electrode hardcoded `table ((0 0) (9.375e-11 18750))`, so the
> voltage existed twice and only the literal mattered. Editing `appliedVoltage`
> changed nothing, silently. Fixed 2026-09-01; all four streamer beds now have
> zero dangling variables.

## The waveform is a Function1 — you are not limited to a table

`uniformFixedValue` takes an OpenFOAM `Function1`, so the electrode can follow
any of these. **No solver code is needed for any of them.**

| what you want | `uniformValue` entry |
| ------------- | -------------------- |
| constant | `constant $appliedVoltage;` |
| linear ramp then hold | `table ((0 0) ($voltageRampTime $appliedVoltage));` |
| linear ramp (explicit) | `linearRamp { start 0; duration $voltageRampTime; }` scaled — see `Scale` below |
| smooth ramp | `halfCosineRamp { start 0; duration $voltageRampTime; }` |
| **sinusoid (AC)** | `sine { frequency 1e3; amplitude 20e3; level 0; }` |
| square / pulse train | `square { frequency 1e3; amplitude 20e3; mark 1; space 1; }` |
| cosine | `cosine { frequency 1e3; amplitude 20e3; }` |
| polynomial in t | `polynomial ((0 0) (2e14 1));` |
| waveform from a file | `tableFile { file "constant/waveform"; }` or `csv { ... }` |
| arbitrary expression | `coded { name V; code #{ return 18750*sin(1e3*x); #}; }` |

`Scale` composes them, e.g. an amplitude-scaled ramp:

```
uniformValue    scale
{
    scale   linearRamp { start 0; duration $voltageRampTime; };
    value   constant $appliedVoltage;
}
```

A **step** is simply `type fixedValue; value uniform $appliedVoltage;` — which is
what `positiveStreamer_AMR` does. The other beds ramp instead, because (their own
comment) *"a step change puts a very large field transient on the electrode in
the first timesteps."* Whether the AMR bed should ramp too is an open question,
flagged in its config rather than silently changed.

## Timestep interaction

Two controls touch the waveform, and it is worth knowing which does what.

**`maxVoltageRisePerStep`** (default 100 V, in `plasmaTimeControl`) caps the
per-step change in the electrode voltage. It is an **absolute** criterion, so
what it buys depends on the amplitude: on the shipped 18.75 kV ramp 100 V/step is
~188 points across the ramp, on a 4 kV waveform the same key gives 40.

Measured 2026-09-01, its cost and benefit on the 0.5 ns LMEA bed:

| setting | steps | time | notes |
| ------- | ----- | ---- | ----- |
| 0 (off) | 47 | 69.2 s | `Co_chem` takes over; the run completes normally |
| 100 V/step | 248 | 188.8 s | binds 41% of steps, 2.7x CPU |

An attempt to calibrate it as a points-per-ramp control (N = 5, 10, 20, 40,
187.5) produced a **non-monotonic** error sequence — 1.09%, 2.15%, 2.07%, 1.15%,
ref — so it does not isolate ramp resolution at all: changing this limiter
reshuffles the whole dt history. **Do not quote a resolution number from that
sweep.** At N = 5 the limiter never binds, because the physics limiters already
deliver more than 5 points across this ramp unaided.

For slow waveforms it is essentially inactive. Slew rate decides:

| waveform | max dV/dt | dt allowed at 100 V/step |
| -------- | --------- | ------------------------ |
| 18.75 kV / 93.75 ps ramp | 2.0e14 V/s | 5e-13 s |
| 20 kV / 100 ns ramp | 2.0e11 V/s | 5e-10 s |
| AC 1 kHz, 40 kV pp | 1.26e8 V/s | 8e-7 s |

At 1 kHz any plasma present demands 1e-13–1e-8 s from dielectric relaxation and
chemistry, so the voltage limiter never binds. **This is arithmetic, not a
measured run** — no AC case has been exercised.

**`maxInitialDeltaT`** (default 1e-12, fresh starts only) caps the **first** step.
This is the one case no other limiter can see: every other limiter here is a
backward difference, so on step 1 there is nothing to difference against, and a
case whose `deltaT` exceeds the ramp has already stepped over it. 1e-12 gives ~94
points across the 93.75 ps ramp (1e-10 would give 0.9), and the 1.2x/step growth
cap recovers it to 1e-8 in ~51 steps and 1e-6 in ~76 — so it costs under 80
warm-up steps even for a slow case. Not applied on a restart, which has already
passed the ramp.

A **relative** field-change limiter was also built and removed unshipped: it was
never the binding constraint in four regimes, because the growth cap plus
`Co_chem` (which keys on `maxChemStateRate()`, the fractional rate of change of
the state over all species) already cover it. See `plasmaTimeControl.H`.

## See also

- `docs/models/energy/lmea.md` — the electron energy closure
- SoEEDF `docs/options-reference.md` §2.5 — the full `plasmaTimeControl` reference
