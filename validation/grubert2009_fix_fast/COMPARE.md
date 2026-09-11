# Does the charge-density fix alone make the FAST (near-fixed-voltage) case work?

```compare
question:  With Poisson's source no longer frozen, does a near-fixed-voltage gap at Grubert's -500 V reach a physical cathode fall instead of running away?
baseline:  /home/kkourtza/soplasma-scratch/validation/grubert2009_fix_n11
varies:    resistance, sourceVoltage
matches:   capacitance, minNumberDensity, driftDiffusionFluxScheme, electronEnergyModel, mesh, chemistry
time:      5e-6
field:     n_e
region:    domain
```

## The question, and why it is worth asking FIRST

The R and C values were introduced to tame a breakdown transient that we now
know was an artefact: `updateChargeDensity()` was running on 0.04% of steps, so
the field was the vacuum field and nothing screened (`4618a51`). **If the fix
is sufficient, the elaborate circuit may not be needed at all** -- and Grubert
imposes a fixed -500 V, not a ballast.

User's proposal, 2026-09-06: go back to something like the fast cases run
before the circuit work, with R small so the gap sits at the applied voltage.

## Settings, and why this is near-fixed-voltage

    R = 1e6 Ohm, C = 5e-14 F, V_source = -501 V

* `tau = R*C = 50 ns`, so the gap reaches its -179 V breakdown at **t = 22 ns**
  rather than 1.8 us. This case answers in MINUTES, which is why it goes first
  (rule 17).
* The load line `V_gap = -501 + R|I|` gives **-500.0 V at I_op** and -490.8 V
  even at 10x I_op. Short-circuit current is 490x the operating point, so this
  is a fixed-voltage boundary in all but name -- matching Grubert's own
  condition -- while the tiny ballast keeps the current from being formally
  unbounded.
* `endTime 5e-6`, `writeInterval 5e-7` -> 10 snapshots, per the user's
  standing request for ~10 writes.

## What makes this a valid control

Compared against `grubert2009_fix_n11`, which differs ONLY in the circuit
(R = 1e8, V_src = -602, i.e. a real ballast). `varies` therefore lists
`resistance` and `sourceVoltage` together, because the load line pins the
second to the first -- they are not independent knobs. Capacitance, floors,
flux scheme, energy model, mesh and chemistry are byte-identical.

**Caveat, stated rather than buried:** this is a two-point comparison of
circuit stiffness, not a sweep. If both arms behave, it says the fix was
sufficient and the circuit is optional; if only the ballasted arm behaves, the
circuit is doing real work and the R sweep is worth repeating properly.

## Reference numbers (Grubert et al. 2009, argon, 100 Pa, 1 cm, -500 V)

| quantity | reference |
|---|---|
| j at the cathode | 0.511 mA/cm^2 |
| n_e peak | 2.478e15 m^-3 |
| n_Ar+ peak | 6.39e15 m^-3 |
| ionisation degree | 1.03e-7 |
| cathode fall thickness | ~4.4 mm (d_c*p ~ 0.33 torr cm) |

## THE OBSERVABLE TO CHECK FIRST, before any density

**`Emag` spread.** With the source frozen the field was uniform at the vacuum
value (12943..12989 V/m, 0.4%). A run whose `Emag` min and max stay within a
few percent has NOT screened, whatever its densities do, and is suspect
immediately. Then: does the field concentrate into a cathode fall of order
mm -- not 0.03 mm, which is what the diverged runs produced.

And per rule 27, before reading any of it:

    grep -c "Charge density updated" logs/log.soPlasmaFoam   # vs step count,
                                                             # must be >= 1:1

---

# RESULT 2026-09-06: the fix was NECESSARY but NOT SUFFICIENT here

**Screening was restored.** `Emag` went from a flat vacuum field (0.4% spread)
to 99-100% structured, with a real cathode fall. And the discharge passed
*through* Grubert's operating point on its way up:

| at t = 111 ns | measured | Grubert | ratio |
|---|---|---|---|
| n_e peak | 2.317e15 | 2.478e15 | **0.93x** |
| n_Ar+ peak | 4.201e15 | 6.39e15 | 0.66x |
| ionisation degree | 9.61e-08 | 1.03e-07 | **0.93x** |
| E/N at peak | 1978 Td | ~2450 (normal fall) | 0.81x |

Against the pre-fix runs, which sat at an ionisation degree of 1.6e-4 to
9.5e-4, that is three to four decades of correction.

**But it did not stop there, and that is correct physics, not a bug.** With
R = 1e6 the gap is fixed-voltage in all but name, and *a fixed-voltage gap
above breakdown has no current limit* -- nothing in the equations selects an
operating point. Measured: `dj/dt` accelerated from 3.8e-5 to 1.97
mA/cm^2/ns while V was still ramping (V = -433, j = -12.4 = 24x the
reference), so this is divergence, not a slow approach.

**Conclusions, and they matter for what comes next:**

1. The frozen-charge-density regression was a REAL and dominant bug. Fixing it
   changed the answer by decades and restored screening.
2. It was NOT the whole story. The ORIGINAL motivation for an external circuit
   stands exactly as first argued: above breakdown the current must be set by
   the circuit. `fix_fast` is the demonstration.
3. Therefore the R sweep regains its purpose -- but now on valid physics. The
   earlier sweep was void because the field was frozen; a repeat is a real
   experiment.
4. The discharge DOES traverse the correct operating point, which is the
   strongest evidence yet that the closure and coefficients are right and only
   the current limiting was missing.

Killed at t = 1.29e-7 once `dj/dt` had established acceleration over five
decades; the ballasted arms `fix_n11`/`fix_n8` were kept running, since they
are the configuration that can select an operating point.
