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

## THE PARENT'S OWN TRAJECTORY -- the numbers this must be judged against

Extracted 2026-09-07. **The parent looks just as mild as the low-floor arms up
to ~3 us, so early mildness is NOT evidence of anything.**

| t [us] | n_e min | n_e max | max/min |
|---|---|---|---|
| 0.50 | 1e11 | 1.047e11 | 1.047 |
| 1.00 | 1e11 | 1.229e11 | 1.229 |
| 2.00 | 1e11 | 1.228e11 | 1.228 |
| **2.95** | 1e11 | 1.302e11 | **1.302** |
| 5.00 | 1e11 | 1.954e11 | 1.954 |
| 8.00 | 1e11 | 7.507e11 | 7.507 |
| **9.70** | 1e11 | 7.481e12 | **74.81** |
| 10.40 | 2.918e13 | 5.281e15 | 181 |

At t = 2.95 us the low-floor arms read max/min = 1.30078 (1e7) and 1.30075
(1e9), against the parent's **1.3018** -- identical to 0.08%. All three are in
the LINEAR regime there, where the solution simply scales with the floor. That
agreement is a consistency check, not a result.

**THE DISCRIMINATING WINDOW IS t = 5 to 10.4 us**, where the parent departs from
linearity (1.95 -> 7.5 -> 74.8 -> 181). An arm that is still linear at 9.7 us
has falsified the parent's behaviour; an arm that turns over there has not.

Note also that by t = 10.4 us the parent's n_e MINIMUM has risen to 2.918e13 --
**292x its own floor** -- so the whole domain, not just a peak, had left the
floor. That is the signature to watch.

---

# RESULT 2026-09-07: RUNAWAY at t = 20.77 us. The floor bought MARGIN, not IMMUNITY.

Left running at the user's instruction rather than stopped.

## What fired

Both pre-registered runaway signatures, 0.77 us after reaching the -200 V plateau:

| signature | at plateau entry (t = 20.0 us) | at t = 20.77 us |
|---|---|---|
| `V_el` | -199.1 V | **+21.6 V -- SIGN FLIPPED** |
| `\|I_cond\|/I_sc` | 4.4e-3 | **1.92** -- beyond what the ballast can supply |
| `n_e,max` | 1.601e12 m^-3 | **1.474e17 m^-3** |
| log-rate `d(ln n_e)/dt` | 1.64e6 1/s | 1.19e7 1/s |

Same failure mode as `../grubert2009_ballast_low` (floor 1e11), which flipped to
+204 V and reached 95.6x `I_sc`. **So lowering the floor from 1e11 to 1e9 did
NOT remove the instability -- it delayed it**, from -104 V to the -200 V
plateau. That is the outcome pre-registered in the floor study as the one that
would say the floor is not the whole story.

## THE THRESHOLD IS A DENSITY, NOT A VOLTAGE -- and the other arms prove it

At the SAME -200 V, with the SAME floor, mesh and circuit:

| arm | ramp | n_e,max at plateau entry | state at -200 V |
|---|---|---|---|
| `grubert2009_ramp2us` | -100 V/us | 2.035e11 | STABLE, log-rate ~2.5e5 1/s |
| `grubert2009_ramp5us` | -40 V/us | 2.378e11 | STABLE, log-rate ~2.8e5 1/s |
| **`grubert2009_floor1e9`** | **-10 V/us** | **1.601e12** | **RUNAWAY** |

The slow-ramp arm arrived at -200 V with **9.2x more density** than the fast
arms, because the discharge grows exponentially throughout the ramp and a slower
ramp integrates more of that growth. It crossed the threshold; they have not.

**COROLLARY, and it inverts an earlier reading.** I reported at -74 V that the
ramp rate "barely matters" (spread 17%). That was measured too early: the spread
GROWS with voltage, reaching 9.2x at -200 V, and it is decisive for stability.
A fast ramp is cheaper per VOLT and reaches a given voltage with LESS density;
a slow ramp is the faster route to a given DENSITY. The two questions have
opposite answers and must not be conflated.

Also visible in the same table: the slow arm reached -200 V with `Te,max` 15.3 eV
against 27.6-29.3 eV in the fast arms. Denser AND cooler -- the signature of
screening beginning, which only the slow arm had enough density to show.

## What this does NOT overturn

The floor result stands as measured: at the 1e11 floor the gap ran away at
**-104 V, below its own 120.8 V breakdown**, driven by space charge the floor
itself bootstrapped (691 ions per clamped electron). At 1e9 it survives to
-200 V and 1.6e12 m^-3. That is a real and large improvement, and the mechanism
is unchanged. What is now clear is that it postpones the instability rather than
removing it, so the floor is a necessary fix and not a sufficient one.

## Still open, and now sharper

The question is no longer "why does it run away" but **"what sets the density
threshold, and can the circuit hold the discharge through it?"** `I/Isc` went
from 4.4e-3 to 1.92 in 0.77 us -- the ballast had three decades of headroom and
lost it in under a microsecond, which is the `tau_loop = C_gap/g` bound doing
exactly what it was measured to do.

`../grubert2009_step250` is the arm to watch: at -250 V it has reached
`n_e,max` = 5.46e15 (past Grubert's 2.478e15) with `V_el` pulled back to
-204.6 V and `I/Isc` = 0.275 -- the ballast genuinely loading for the first time,
and not yet flipped.
