# grubert2009_n11_clamp{100,Derived} -- what the 100 eV clamp actually cost

Created 2026-09-07 at the user's instruction, after they asked whether the
mean-energy clamp could have destroyed earlier results.

## The question

`meanEnergyMax` defaulted to a hardcoded **100 eV** while the LMEA tables cover
**2644 eV** -- so the clamp was 26.4x tighter than the "range the tables cover"
justification it carried, and it DISCARDED ENERGY in cells inside the tabulated
range (`c131d8c` fixes the default to be derived). It fired on **91,586** outer
iterations of `grubert2009_fix_n11`, at cells with `n_e = 7.2e18` -- dense
plasma, not floor cells -- with raw mean energy reaching 105 eV.

**Does removing it change the answer?**

## The two arms -- A VALID ONE-VARIABLE CONTROL

    /home/kkourtza/soplasma-scratch/validation/grubert2009_n11_clamp100
    /home/kkourtza/soplasma-scratch/validation/grubert2009_n11_clampDerived

Both are `cp -r` of `grubert2009_fix_n11`, restarted from the SAME snapshot
(`startFrom startTime`, `startTime 2.000001e-06`, `timePrecision 12`), to
`endTime 2.26e-6`, `writeControl timeStep`, `writeInterval 200`.

    DIFFERS:  clamp100 sets `energyModelCoeffs { meanEnergyMax 100; }`
              clampDerived omits it -> DERIVED 2644.46 eV
    MATCHES:  everything else. Verified with `diff -rq`: the two trees differ in
              exactly ONE file, constant/plasmaSpeciesProperties.

**Why both arms are re-run rather than comparing against the original log.** A
restart re-initialises from a written snapshot and does not reproduce a
continuous run bit-for-bit, so the original trajectory is NOT a valid baseline
for a restarted arm. Running both from the same snapshot makes the clamp the
only difference. (Rule 15, and [[baseline-contamination]].)

## Why this window

`fix_n11`'s FIRST clamp fired at `t = 2.2418e-06`, and its last step was
`2.2432e-06` -- the clamp was active for the final **1.4 ns** only, long after
the healthy assessment at `t = 1.86 us`. So the window that matters is
2.000001e-06 -> 2.26e-6: it starts clean, crosses where the clamp began, and
runs slightly past where the original stopped.

**Prediction, registered before running** (so it can be wrong): the arms agree
until ~2.2418e-06 and then diverge, with `clampDerived` ionising HARDER --
because `k_ion` rises with mean energy across the whole table, so the clamp
SUPPRESSES ionisation. If `clampDerived` instead runs COOLER, the reasoning
above is wrong and must be reworked, not explained away.

## Discriminating observables

| observable | agreement means | divergence means |
|---|---|---|
| time of first `eV clamp on` in clampDerived | none at all -> 2644 eV is never reached | it clamps too, and the fix is insufficient |
| `n_e` max, `meanE` max vs t | clamp was cosmetic | it was load-bearing |
| `I_cond` vs t | earlier conclusions safe | quantitative results need revisiting |
| step at which each dies / stalls | -- | -- |

## Extraction

    grep -c "eV clamp on" logs/log.soPlasmaFoam        # must be 0 for clampDerived
    grep -aF "e [m^-3]:" logs/log.soPlasmaFoam         # NOTE the -F: without it this matches ZERO lines
    ~/ct-env/bin/python extract_VI.py
    postProcessing/externalCircuit/circuit.csv

## What this does NOT settle

The two arms share `R = 1e8` and an unseeded floor start, so this is a control
for the CLAMP only. It says nothing about whether Grubert's profile is a fixed
point -- that is `../grubert2009_steady`.
