# Is the energy convective Courant an EFFICIENCY knob or an ACCURACY limit?

```compare
question:  Does raising maxEnergyConvectiveCo change the SOLUTION, or only the cost?
baseline:  /home/kkourtza/soplasma-scratch/validation/grubert2009_iset
baseline:  /home/kkourtza/soplasma-scratch/validation/grubert2009_iset_En150
varies:    maxEnergyConvectiveCo
matches:   maxSpeciesCo, circuit, minNumberDensity, driftDiffusionFluxScheme, electronEnergyModel, mesh, chemistry, gammaSEE, electronReflection, endTime
time:      1e-6
field:     n_e
region:    domain
```

## Why this sweep, and why NOW

`grubert2009_iset` is limited almost entirely by ONE thing. Measured over its
first 8315 steps:

    7970 steps :  Co_conv (energy)     <- the electron-ENERGY convective Courant
     321       :  Co_diff
      18       :  PI controller
       6       :  errorMaxGrow
    last 500   :  500/500  Co_conv (energy)

and the run's own temporal-error controller says dt could be **298x larger**:

    temporal err ||e|| = 3.76e-08  [target 1]  ->  dt would be 5.72e-09
    actual dt                                                  1.92e-11

So ~300x of wall clock is being spent satisfying a Courant number on the energy
equation rather than an accuracy requirement. At 6.2e-8 sim-s per wall-minute
that is the difference between a 5-hour run and a 1-minute one.

**This has to be re-measured rather than taken from the earlier sweep.** The
`ballastCo_1000` / `ballastCo_10000` cases were run under the no-op-step defect
(45-50% of timesteps executing zero correctors, `b6a7102`), so they are among
the ten cases whose timings are void.

## What differs

**Only `maxEnergyConvectiveCo`.** The species cap stays at 15 in every arm, so
the sweep isolates the limiter that actually binds. `maxEnergyConvectiveCo`
otherwise DEFAULTS to `maxSpeciesConvectiveCo` (`plasmaTimeControl.C:374`),
which is why the baseline needs no key at all:

| arm | `maxEnergyConvectiveCo` | how |
|---|---|---|
| `grubert2009_iset` | 15 | defaulted from `maxSpeciesCo` |
| `grubert2009_iset_En150` | 150 | set explicitly |
| `grubert2009_iset_En1500` | 1500 | set explicitly |

Both the config variable and the `plasmaSimulationControls` reference were
added, because a `$variable` nothing references is inert and silently changes
nothing -- measured 2026-09-01, and it invalidated a run then.

## THE OBSERVABLE: do the trajectories OVERLAY?

The current is imposed and identical in all three arms, so this is a pure
numerics question:

* **`V_gap(t)` before ignition** is the sharpest test available, because its
  slope is analytic: `dV/dt = I_set/C_gap`, verified on the baseline to 1.33%.
  A Courant number that is only an efficiency knob cannot change it. Any
  deviation is discretisation error, visible immediately and cheaply.
* **the ignition time** (`V_gap` reaching ~-121 V, where `alpha*d` = 2.872).
* **`n_e(t)` and `n_Ar+(t)`** through ignition.
* **cost**: steps, steps/s, and sim-s per wall-minute.

**Judge cost only where the answers agree.** A 100x speedup on a different
trajectory is not a speedup.

## Endpoint, and why it is short (rule 17)

`endTime 1e-6`, just past the ~7e-7 s ignition. The discriminating behaviour --
the analytic pre-ignition ramp and the ignition itself -- is fully contained
there, and it costs minutes rather than the ~5 h a run to 20 us needs.

## Failure modes, stated in advance

* **trajectories diverge** -> the energy Courant is an ACCURACY limit at this
  cap and cannot be relaxed. That is a real answer and it closes the question.
* **rejections appear** (`retryStep`) -> the outer loop cannot converge at the
  larger dt, so the nominal speedup is eaten by discarded steps. Count
  rejections, not just steps.
* **the contraction diagnostic worsens** -> the baseline already reports
  `rho [contraction] max 48.4, residual GREW` at dt = 1.9e-11, which is why
  this is a measurement and not an assumption. If a larger dt makes that worse,
  the outer loop is the real constraint, not the Courant number.
* **no-op steps** -> would mean the `b6a7102` regression has returned. The step
  audit is fatal on it now, so this cannot pass silently.

## Extraction

    for a in iset iset_En150 iset_En1500; do
      L=../grubert2009_$a/logs/log.soPlasmaFoam
      echo "$a: steps=$(grep -ac '^Time = ' $L) rejections=$(grep -ac retryStep $L)"
      grep -a 'deltaT set by:' $L | awk '{print $NF}' | sort | uniq -c | sort -rn | head -3
    done
    tail -1 postProcessing/externalCircuit/circuit.csv   # time,I_set,I_cond,V,g

---

# RESULT 2026-09-06, and it moved the question on

| arm | dt vs base | error in V_gap | discards | binding limiter |
|---|---|---|---|---|
| En150 | 1.5x | 0.054% | 0 | `Co_conv` (species) |
| En1500 | 1.6x | 0.055% | 0 | `Co_conv` (species) |
| All150 (all caps 150) | 15x | 0.105% | 0 | `Co_conv (energy)` |
| All1500 (all caps 1500) | **30x** | **0.267%** | 1 | **`coupling margin (hold)`** |

**The energy cap alone buys almost nothing** -- 1.5x, because the SPECIES cap
takes over immediately. En150 and En1500 are identical to the digit (12742
steps each), which is the clean confirmation.

The question therefore moved to `grubert2009_iset_All150`, which raises ALL the
caps and finds ~30x of dt for ~0.3% in V_gap -- and, at 1500, reveals that the
binding limiter is no longer Courant at all but the omega-based
`coupling margin`. See that case's COMPARE.md and
`docs/design/adaptive-dt-design.md`.

**The answer to this case's own question is therefore YES, unambiguously:** the
energy convective Courant is an EFFICIENCY knob, not an accuracy limit --
0.055% of V_gap for a 100x change in the cap.
