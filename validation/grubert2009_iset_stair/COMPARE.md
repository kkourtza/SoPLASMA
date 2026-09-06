# Trace the V-I characteristic by holding the current, not ramping through it

```compare
question:  Held at constant current until gamma relaxes to 0, does the discharge settle on the V-I characteristic and pass through Grubert's operating point?
baseline:  /home/kkourtza/soplasma-scratch/validation/grubert2009_iset
varies:    setCurrent, maxSpeciesCo
matches:   circuit type, compliance, capacitance, minNumberDensity, driftDiffusionFluxScheme, electronEnergyModel, mesh, chemistry, gammaSEE, electronReflection
time:      7.8e-5
field:     n_e
region:    domain
```

## Why a staircase, from the measurement that killed the ramp

The linear `setCurrent` ramp reached ignition correctly -- pre-ignition
`dV/dt` matched `I_set/C_gap` to **1.33%**, ignition at V = -116 V against a
predicted -121 V (**4%**) -- and then ran away. The cause was measured, and it
was NOT the timestep:

* through ignition the ramp's e-folding time was **945 ns**, against an ION
  TRANSIT of **3.7-6.7 us**. So it drove the current **4-7x faster than the
  ions could follow**; the discharge was never quasi-static and `gamma` never
  returned to 0.
* in that regime the loop response bound `tau_loop = C_gap/g` was **238x too
  slow** to contain the excursion (`g` = 1.86e-09 S at onset against a
  threshold `C_gap/tau_growth` = 4.43e-07 S). `g` was verified as the true
  plant gain against finite differences at 1/10/100/1000-step windows.
* **dt was conclusively excluded:** four arms over six decades of dt (1e-12,
  6e-12, 1.2e-10, 2.4e-10) followed ONE trajectory, agreeing to 0.66% through
  the overshoot, peak and turnover, to the same n_e ~ 1e18.

At `gamma = 0` the stability analysis says EVERY circuit is stable, including
the ballast that ran away. **The destination is stable; the approach was not.**
So hold each current level long enough for `gamma` to relax, and each plateau
becomes a genuine steady state.

## The staircase

    setCurrent table ((0 5e-8) (2e-6 5e-8) (1.4e-5 5e-8)
                      (1.8e-5 2e-7) (3e-5 2e-7)
                      (3.4e-5 5e-7) (4.6e-5 5e-7)
                      (5e-5 1.022e-6) (6.2e-5 1.022e-6)
                      (6.6e-5 2e-6) (7.8e-5 2e-6));

| plateau | I_set | j [mA/cm^2] | window |
|---|---|---|---|
| 1 | 5e-8 | 0.025 | 2-14 us |
| 2 | 2e-7 | 0.100 | 18-30 us |
| 3 | 5e-7 | 0.250 | 34-46 us |
| **4** | **1.022e-6** | **0.511** | **50-62 us** <- GRUBERT |
| 5 | 2e-6 | 1.000 | 66-78 us |

**12 us holds = 2-3 ion transits**, so `gamma -> 0` on each. **Risers RAMPED
over 4 us, not stepped** -- a discontinuous jump in `I_set` would create
exactly the fast excursion this design exists to avoid.

**`maxSpeciesCo = 1500`**, which makes this affordable. Justified by
measurement, not hope: raising all Courant caps was shown to follow the SAME
trajectory (the four-arm dt study above), and the arms that diverged did so for
the controller reason, not the caps. It buys 15-30x.

## What this measures that nothing else can

**The V-I characteristic itself.** Under voltage control the normal-glow branch
is not a function -- V is flat while j varies -- so no voltage-driven case can
resolve it. The branch structure is:

    Townsend     V rises with I
    SUBNORMAL    V FALLS with I     <- negative differential resistance
    normal glow  V flat, j pinned at j_n

**The sign of `dV/dj` is the discriminator, and no fitted coefficient can fake
it.** That is a far stronger validation than matching one operating point.

## Success and failure, stated in advance

* **SUCCESS:** V settles on each plateau (`dV/dt -> 0` with
  `I_cond -> I_set`), and plateau 4 gives **V ~ -500 V at j = 0.511** --
  Grubert's point, with the current imposed so V is the prediction.
* **the plateaus never settle** -> 12 us is not enough; extend the holds. The
  settling time itself is then the result.
* **runaway on a riser** -> 4 us risers are still too fast; the quasi-static
  requirement is tighter than one ion transit.
* **V settles but at the wrong value** -> a real physics discrepancy, and
  finally a clean one: measured against a stable state rather than a transient.
* **the compliance rail (-1500 V)** -> the discharge cannot carry that current
  at any voltage.

## Extraction

    ./extract_VI.py                      # I_set, I_cond, j, V, dV/dj, branch
    ~/ct-env/bin/python ../../../Projects/SoEEDF/validation/dias2025/plot_grubert_profiles.py \
        centreline_<t>.csv --arm LMEA    # profiles vs Grubert fig 3

`extract_VI.py` bins in `I_set` rather than time and refuses to name a branch
while the |V| peak sits at an endpoint.

Reference: Grubert et al., Phys. Rev. E 80 (2009) 036405 -- argon, 100 Pa,
1 cm, j = 0.511 mA/cm^2 at V = -500 V, n_e 2.478e15, n_Ar+ 6.39e15.
