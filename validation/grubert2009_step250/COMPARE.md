# Stepping up in voltage from a stable Townsend discharge

Created 2026-09-07. Three arms, all built from `grubert2009_ramp1us` after it was
stopped cleanly with `stopAt writeNow` at t = 2.02945e-06 s.

## Why this is worth doing NOW rather than after a density plateau

`grubert2009_ramp1us` reached the -200 V plateau and became quasi-steady in
everything that equilibrates fast: bulk `Te` fixed at **2.31 eV**, `Te,max` flat
at ~29.4 eV, structure steady, and the electrode tracking the source to 0.2 V.
What it did NOT do is plateau in DENSITY -- it settled to a CONSTANT FRACTIONAL
growth rate of ~2.1e5 1/s (21% per microsecond, doubling every 3.3 us).

A constant fractional rate is exponential growth, NOT balance. And it cannot
balance yet, because the discharge is too thin to matter to anything:

    n_e,max     2.5e11 m^-3      9.5e-05 of Grubert's 2.478e15
    j           5.4e-4 mA/cm^2   1.1e-03 of Grubert's 0.511
    lambda_D    2.34 cm          2.34x the 1 cm GAP -> CANNOT screen
    |I_cond|/I_sc  ~9e-4         the ballast is not limiting anything
    tau_loop = C_gap/g ~ 42 us   vs 1/rate ~ 4.8 us -> circuit 9x too slow to arrest

So there is no circuit-stabilised steady state to wait for at -200 V. Waiting
9.2 e-folds for the density to reach glow conditions costs ~43 us of simulated
time, about **9.5 hours of wall clock**. Raising the voltage raises the
multiplication and gets there faster -- and it is what a real experiment does.

**THE REGIME REACHED IS A TOWNSEND (DARK) DISCHARGE, NOT A GLOW.** Ionisation is
self-sustaining (multiplication 9.4 > 1 at -200 V) but there is no space charge
structure, no cathode fall. That distinction must not be lost: this is the first
time this gap has held a stable self-sustaining discharge at all, and it is
still 3 decades in current below the operating point being validated against.

## The arms

    grubert2009_step250       RESTART from t = 2.02945e-06 s of ramp1us;
                             sourceVoltage -200 -> -250 V over 0.25 us
                             (the same -200 V/us riser), then held.
    grubert2009_ramp1us_300  FRESH from t = 0; 1 us ramp to -300 V, then held.
    grubert2009_ramp1us_400  FRESH from t = 0; 1 us ramp to -400 V, then held.

All three inherit the density floor of **1e9 m^-3**, the 400-cell graded mesh,
and `seriesResistor` R = 1e8 ohm from `grubert2009_ramp1us`.

Baselines -- ABSOLUTE PATHS:

    /home/kkourtza/soplasma-scratch/validation/grubert2009_ramp1us    (stopped, -200 V, the parent)
    /home/kkourtza/soplasma-scratch/validation/grubert2009_ramp2us    (-200 V, running)
    /home/kkourtza/soplasma-scratch/validation/grubert2009_ramp5us    (-200 V, running)
    /home/kkourtza/soplasma-scratch/validation/grubert2009_ballast_low (floor 1e11, RAN AWAY at -104 V)

**step250 is a CONTINUATION, not a one-variable control.** It changes the
voltage from a state produced by a prior run; it is not comparable to a
from-scratch case at -250 V. The 300 and 400 V arms ARE mutually comparable and
comparable to `ramp1us`, since all three are fresh 1 us ramps differing only in
the final voltage.

## Multiplication is EXPONENTIAL in voltage -- these are not equal steps

Uniform-field Townsend multiplication `M = gamma*(exp(alpha*d)-1)`, gamma = 0.06:

| V [V] | M | vs -200 V |
|---|---|---|
| 200 | 9.41 | 1.0x |
| 250 | 30.6 | 3.3x |
| 300 | 83.6 | 8.9x |
| 400 | 404 | **43x** |

So -400 V is a 43x jump in multiplication over the state we know is stable. That
is the point of the test, and also the reason it may fail.

## Success and failure, stated in advance

* **STABLE at all three** -> the low-floor configuration is robust well above
  breakdown, and voltage can be raised freely to reach glow density fast.
* **STABLE at 250/300, UNSTABLE at 400** -> there is a multiplication ceiling
  between 84 and 404, worth locating.
* **ANY arm RUNS AWAY** -- electrode sign flip, or `|I_cond|` > `I_sc` -> the
  stability is voltage-limited and the staircase must step more finely. This is
  the failure mode `grubert2009_ballast_low` showed at -104 V with the 1e11
  floor, so a recurrence would say the floor fix bought margin, not immunity.
* **dt COLLAPSE without runaway** -> numerically limited rather than physically
  unstable; report the limiter that binds.

## Extraction

    ~/ct-env/bin/python ../../tools/settle_check.py <case>   # state + RUNAWAY flags
    ~/ct-env/bin/python ../../tools/news.py <case>           # n_e and Te evolution
    ~/ct-env/bin/python ../../tools/compare_at_voltage.py <A> <B>
