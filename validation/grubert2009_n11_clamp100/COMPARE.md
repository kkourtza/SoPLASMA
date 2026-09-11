# grubert2009_n11_clamp{100,Derived} -- what the 100 eV clamp actually cost

Created 2026-09-07 at the user's instruction, after they asked whether the
mean-energy clamp could have destroyed earlier results.

## The question

`meanEnergyMax` defaulted to a hardcoded **100 eV** while the LMEA tables cover
**2644 eV** -- so the clamp was 26.4x tighter than the "range the tables cover"
justification it carried, and it DISCARDED ENERGY in cells inside the tabulated
range (`8f5bab4` fixes the default to be derived). It fired on **91,586** outer
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

---

# RESULT 2026-09-07: STOPPED EARLY, question answered. The clamp was LOAD-BEARING.

Both arms stopped deliberately at t ~ 2.107e-06 of a 2.26e-06 endTime (~7% of
the window) because they had already answered, and because `clampDerived` was
evolving in a regime the LMEA closure cannot represent -- further wall time
would not have been physics. ~880k steps would have been needed to finish.

## Final state, same physical time, ONE variable apart

| quantity | clamp100 | clampDerived | factor |
|---|---|---|---|
| meanE range [eV] | **[6.97, 100.0]** | **[816.7, 2644.5]** | min 117x |
| n_e max [m^-3] | 1.303e18 | **2.878e19** | **22x** |
| n_e min [m^-3] | 5.489e15 | 9.987e15 | 1.8x |
| Joule/loss, domain | **4.45** | **721** (rising) | **162x** |
| steps | 39372 | 47185 | |

**The clampDerived arm put its ENTIRE DOMAIN above 816 eV.** A glow bulk must be
a few eV, so that arm is unphysical everywhere, not merely at a hot spot.
`clamp100` held a sane bulk at 6.97 eV throughout.

## The verdict, and it inverts the change that prompted it

`8f5bab4` made `meanEnergyMax` DERIVE from the table range (2644 eV) instead of
a hardcoded 100 eV, on the argument that the clamp exists to prevent
extrapolation so its value should be the tabulated range. **That argument is
sound and the outcome is still wrong**, because a table's EXTENT is not its
VALIDITY: the tables are self-consistent to 2644 eV only for a FREELY GROWING
swarm, where `nu_i*U` carries up to 99.4% of the power budget
(`docs/design/electron-energy-balance.md`). The fluid energy equation is robustly
dissipative only below ~1000-2450 Td, i.e. `U` ~ 15-35 eV.

So the old 100 eV was an unsourced constant AND it was holding the solution
inside the regime the closure can represent. **These arms are the measurement
that proves it was load-bearing rather than cosmetic.**

## What NOT to do

Do not simply revert `meanEnergyMax` to 100 eV. That restores the old behaviour
for a reason nobody understood and re-hides this. The three candidate fixes are
in `docs/design/electron-energy-balance.md`; the evidence does not yet choose between
them, and the one independently justified is tightening the `n_e`/`nEps`
coupling, since the growth cancellation `dU/dt = (Joule - Loss) - U*nu_i` arises
only from `n_e` in the DENOMINATOR and is a difference of two large terms at
high `U`.

## Superseded by this result

My earlier statement that the clamp's damage was "bounded -- a few percent, not
orders" was based on the single-ITERATE overshoot (105 vs 100 eV). On the
TRAJECTORY it is 22x in peak `n_e` and 162x in the source balance. The
protective DIRECTION still holds (clamping suppresses ionisation, so clamped
runs under-predict growth rather than inventing it), which is what preserved the
qualitative Grubert conclusions -- but "a few percent" was wrong.
