# What do the Courant caps actually cost, and what governs dt when they are lifted?

```compare
question:  Raising ALL Courant caps 10x and 100x -- how much dt is bought, how much accuracy is lost, and what limiter takes over?
baseline:  /home/kkourtza/soplasma-scratch/validation/grubert2009_iset
baseline:  /home/kkourtza/soplasma-scratch/validation/grubert2009_iset_All150
varies:    maxSpeciesCo
matches:   circuit, minNumberDensity, driftDiffusionFluxScheme, electronEnergyModel, mesh, chemistry, gammaSEE, electronReflection, endTime
time:      1e-6
field:     n_e
region:    domain
```

## The question (user, 2026-09-06)

> "why dont we run in parallel another case where ALL Co are not capped to 15
> but higher? 150 or/and 1500?"

`maxSpeciesCo` feeds `maxSpeciesConvectiveCo` AND `maxSpeciesDiffusiveCo`, and
`maxEnergyConvectiveCo` DEFAULTS to the convective one
(`plasmaTimeControl.C:374`) -- so one variable raises all three caps together.

| arm | `maxSpeciesCo` | conv / diff / energy |
|---|---|---|
| `grubert2009_iset` | 15 | 15 / 15 / 15 |
| `grubert2009_iset_En150` | 15 | 15 / 15 / **150** (energy only) |
| `grubert2009_iset_En1500` | 15 | 15 / 15 / **1500** |
| `grubert2009_iset_All150` | **150** | 150 / 150 / 150 |
| `grubert2009_iset_All1500` | **1500** | 1500 / 1500 / 1500 |

## RESULT, measured 2026-09-06 at each arm's own time against the baseline history

| arm | dt vs base | **error in V_gap** | error in I_cond/I_set | discards | binding limiter |
|---|---|---|---|---|---|
| En150 | 1.5x | 0.054% | 1.0% | 0 | `Co_conv` (species) |
| En1500 | 1.6x | 0.055% | 1.0% | 0 | `Co_conv` (species) |
| All150 | **15x** | **0.105%** | 2.7% | 0 | `Co_conv (energy)` |
| All1500 | **30x** | **0.267%** | 3.2% | **1** | **`coupling margin (hold)`** |

Three conclusions:

1. **~30x of dt for ~0.3% in V_gap.** The caps were costing an order of
   magnitude and more for a fraction of a percent.
2. **Raising the ENERGY cap alone buys almost nothing** (1.5x), because the
   SPECIES cap immediately takes over. That is why En150 and En1500 are
   indistinguishable to the digit -- 12742 steps each.
3. **At All1500 the binding limiter is NO LONGER COURANT.** It is
   `coupling margin (hold)`, the omega-based outer-coupling governor. **The
   automatic dt controller already exists; the hand-set caps at 15 were masking
   it.** See `doc/adaptive-dt-design.md`.

`coupling margin` is keyed on the Aitken relaxation factor `omega`, which is
DIMENSIONLESS -- so its threshold means the same thing in any case, satisfying
the universality requirement the caps do not. The code comment states the
reason it is preferred over a non-convergence verdict: *"a diverging step can
grind indefinitely inside the stiff chemistry integrator and never reach the
corrector cap, so that verdict may never arrive. omega is available every
corrector."* That also explains why discards are near-nonexistent -- the margin
backs dt off BEFORE a step fails.

## Caveats, stated rather than buried

* **The comparison is at t ~ 5-6.6e-7 s, BEFORE the cathode fall forms.** The
  errors may grow once the fields stiffen; re-measure at endTime.
* **The error GROWS with dt** (0.054 -> 0.105 -> 0.267%), which is the expected
  first-order behaviour and is the argument for letting the ACCURACY controller
  decide where to stop, not a Courant number.
* **All1500's single discard is one event, not a characterised boundary.** It
  is the first discard in ~1.35 MILLION steps across all cases, so the
  robustness path is only now being exercised at all (rule 23).
* Every timing here post-dates `b54ff16`. Earlier Courant work
  (`ballastCo_1000`, `ballastCo_10000`) is void: 45-50% of its timesteps ran
  zero correctors.

## Extraction

    for a in iset iset_En150 iset_En1500 iset_All150 iset_All1500; do
      L=../grubert2009_$a/logs/log.soPlasmaFoam
      printf "%-14s steps=%-8s discards=%-4s limiter=%s\n" $a \
        "$(grep -ac '^Time = ' $L)" "$(grep -ac 'DISCARDING this step' $L)" \
        "$(grep -a 'deltaT set by:' $L | tail -1 | awk '{$1="";$2="";$3="";print}')"
    done

Accuracy is compared at EACH ARM'S OWN TIME against the baseline's history --
not at the arms' endpoints, which are different instants. Comparing different
instants is not a comparison.
