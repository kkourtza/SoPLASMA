# Is the DENSITY FLOOR bootstrapping its own breakdown?

Created 2026-09-07, from the user's question: is the floor/initial value of n_e
and the ions too high?

## The question, and the arithmetic that motivates it

`minNumberDensity` is not just an initial value -- `0/n_e` is written
`internalField uniform 0` and the per-step clamp RAISES cells to the floor
("clampNumberDensity(), which raises the field rather than capping it"). So the
floor is a CONTINUOUS volumetric electron source, every step, everywhere.

And it is amplified, because ions leave 307x slower than electrons. At 100 V:

    electron transit  tau_e = 3.48e-08 s
    ion drain         tau_i = 1.07e-05 s        ratio 307
    ions made per clamped electron = N*k_iz*tau_i = 691

| n_e floor | n_Arp accumulated | E_sc [Td] | **E_sc/E_applied** |
|---|---|---|---|
| **1e11 (parent)** | 6.91e13 | 518 | **1.251** |
| **1e9 (this arm)** | 6.91e11 | 5.18 | **0.0125** |
| **1e7 (this arm)** | 6.91e9 | 0.052 | **0.000125** |

The applied field at 100 V is 414 Td. **At the 1e11 floor the accumulated ion
space charge alone makes a field 1.25x the applied one** -- enough to distort it
and drive a sub-breakdown gap into local breakdown.

**AND IT PREDICTS THE OBSERVATION.** The parent run showed `n_Arp` = 3.6e13 to
1.6e14 at t = 10 us; the floor arithmetic predicts **6.91e13**, inside that
range. That is why `grubert2009_ballast_low` ran away at `V_src` = -104 V, BELOW
the computed 120.8 V breakdown: the floor bootstrapped it.

## Baseline -- ABSOLUTE PATH, and this IS a valid one-variable control

    /home/kkourtza/soplasma-scratch/validation/grubert2009_ballast_low

`cp -r` of that case with **only `minNumberDensity` changed** (electron and Arp
together, so the initial state stays quasineutral). `diff -rq` confirms the trees
differ in exactly ONE file, `constant/plasmaSpeciesProperties`. Ar2p is left at
1e5 -- it is a token floor on a minor species and changing it would add a second
variable.

TWO ARMS so the floor itself gets a convergence check rather than a single
point: if 1e9 and 1e7 agree, the answer is floor-independent there and the
mechanism is settled; if they differ, the floor is still binding at 1e9.

## What the parent did, for comparison

| quantity | parent (floor 1e11) |
|---|---|
| healthy until | `V_src` = -96.6 V, drawing **1% of I_sc** |
| then, in 0.76 us | `I_cond` -> **95.6x I_sc**, `V_el` -> **+204 V** |
| froze at | t = 1.04291e-05 s, `V_src` only -104.3 V |
| `n_e` max | 2.57e18 |

## Success and failure, stated in advance

* **SUCCESS:** the arm passes t = 1.04291e-05 s (where the parent froze) and
  keeps `|I_cond| < I_sc`, settling into a mild Townsend/subnormal discharge as
  `V_src` continues its ramp. `n_e` stays orders below 2.57e18.
* **PARTIAL:** it survives longer but still eventually flips the electrode
  positive -> the floor delays rather than prevents, and the accumulation is
  driven by the growing plasma rather than the floor.
* **FAILURE (mechanism refuted):** it runs away at the same time and voltage as
  the parent -> the floor is NOT the source, and this hypothesis dies rather
  than being patched with a still-lower floor.
* **NUMERICAL:** `meanE = nEps/n_e` degrades in near-empty cells. The energy
  model's own `electronDensityFloor` defaults to 1.0, so even 1e7 sits 1e7x
  above it -- but watch the reported mean-energy range for noise at the low end.

## Extraction

    ~/ct-env/bin/python -c "..."   # I_cond vs I_sc = |V_src|/R, R = 1e8
    postProcessing/externalCircuit/circuit.csv
    grep -aF "e [m^-3]:" logs/log.soPlasmaFoam
    grep -a "LMEA mean energy" logs/log.soPlasmaFoam

## Note on what this does NOT test

The mesh. No convergence study on the graded mesh has ever been completed --
`../RESOLUTION_TEST.md` designed one and its arms died at t = 1.8e-7 / 6.3e-7 s,
pre-breakdown, so it produced no evidence. That question stays open and is
independent of this one.
