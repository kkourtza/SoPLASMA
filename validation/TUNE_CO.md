# Courant-headroom test — is the error controller allowed to be in charge?

Written 2026-09-06, BEFORE the runs.

## The question

On the graded mesh the LMEA arm runs at dt = 8.5e-12 s, with

    Co_conv (energy):  15       [ max: 15 ]  <-- BINDING
    Co_conv (e):        9.98    [ max: 15 ]
    Co_diff (e):        8.89    [ max: 15 ]
    Co_chem:            0.0053  [ max: 0.9 ]   -- 170x slack

against an electron transit time of 5.4e-8 s and an ion transit of 4.5e-6 s.
That is ~6400 steps per electron transit in a phase where n_e sits at the 1e11
floor and the field is uniform. Nothing is happening and it is being resolved
to a part in 1e4.

The framework's own start-up banner states the intent being violated:

    "The Courant and stiffness limits still apply as CAPS and should now rarely
     bind -- if one of them is still naming itself in `deltaT set by`, the
     controller is not in charge."

`Co_conv (energy)` is naming itself. The temporal-error PI controller, which
measures the error directly rather than bounding a proxy, is not governing.

## Why the energy/electron RATIO is not the lever

Already measured, in /home/kkourtza/co-sweep/COMPARE_coSweep.md: "Raising the
energy key ALONE saturates at 2.8x because Co_conv (e) then binds -- so the
thing actually under test is one global convective Courant." We see the same
here: the energy Courant is 1.50x the electron one (it is the muEps/mu ratio,
so it moves with mean energy; it was 3.0 on the coarse mesh). Decoupling the
two caps buys at most that ratio and then the electron Courant takes over.
So this test moves ONE global Co, which is what `maxSpeciesCo` already feeds.

## A structural note on the graded mesh

dt is set by (smallest cell)/(largest velocity), and after grading those two
CO-LOCATE at the electrode -- the 1.35 um cell sits exactly where the field is
highest. That is why refinement cost ~6x in dt rather than the ~2x the cell
count suggests, and it is an argument for the error controller over a static
geometric cap.

## Arms

Baseline: `/home/kkourtza/soplasma-scratch/validation/grubert2009` (Co = 15),
the production LMEA arm. Identical mesh, chemistry, BCs and numerics; the ONLY
difference is `maxSpeciesCo`.

| case | maxSpeciesCo |
|---|---|
| `grubert2009` (baseline) | 15 |
| `tuneCo_50` | 50 |
| `tuneCo_100` | 100 |
| `tuneCo_1000` | 1000 |
| `tuneCo_10000` | 10000 |

Co = 1000 and 10000 added 2026-09-06 to find where this actually breaks: at
Co = 100 the limiter was STILL `Co_conv (energy)`, so the error controller had
not taken over even there and the headroom was not yet bounded from above. A
sweep that never reaches the failure it is looking for has not found a limit,
it has only found that it did not look far enough.

endTime 2.5e-6 for the two variants -- just past breakdown, which on the coarse
mesh was V = 185 V at t = 1.88e-6. No variant needs to run further to answer
this.

## THE INVARIANT, and it is not the error norm

**THE BREAKDOWN VOLTAGE MUST NOT MOVE.** It is the one thing in this phase with
a physical value we can check independently (Paschen for argon at
pd = 0.75 Torr cm sits near its minimum, ~140-200 V). If breakdown shifts by
more than a few per cent, the large Co is WRONG no matter what the temporal
error norm reports -- a small measured error on a mis-timed breakdown is a
small error in the wrong solution.

Secondary, and only meaningful if the invariant holds:
* does `deltaT set by` become `temporal error (PI)` rather than `Co_conv`?
* what dt does the controller then choose, and how many steps does it take?

## PREDICTION

* Co = 50: the controller takes over, dt rises several-fold, breakdown voltage
  unchanged.
* Co = 100: dt saturates rather than rising proportionally -- the co-sweep found
  Co = 50 already coupling-margin limited on its case (141 vs 151 steps) -- so
  the gain from 50 to 100 should be small.

**If the breakdown voltage moves, this is retracted and Co stays at 15.** The
speedup is not worth a shifted ignition point, which is the one number here
with an external reference.

## COMPARE dt AT A COMMON SIMULATED TIME, not at a common wall-clock moment

dt here is proportional to (cell size)/(mu_e * E), and E rises with the ramp, so
dt FALLS as the run proceeds. The arms advance at different rates, so at any
wall-clock instant they sit at different voltages and their dt values are not
comparable. First reading, and an example of exactly that trap:

    Co=15   dt 7.93e-12 at t = 8.53e-07 (V = 85 V)
    Co=50   dt 1.02e-10 at t = 1.77e-07 (V = 18 V)
    Co=100  dt 1.52e-10 at t = 2.56e-07 (V = 26 V)

The apparent 13x from Co=15 to Co=50 is mostly the 4.8x lower voltage, NOT the
Courant cap. Every comparison below must be taken at the same t.

## Extraction

    for d in grubert2009 tuneCo_50 tuneCo_100; do
      python - <<'PY'
    # breakdown voltage = V at which |j| first exceeds 0.5 mA/cm2
    PY
    done
