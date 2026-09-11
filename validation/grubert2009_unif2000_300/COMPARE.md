# UNIFORM refined mesh as a BENCHMARK -- is the current reversal a mesh artefact?

Created 2026-09-07 at the user's suggestion: run a uniform, well-refined mesh as
a reference, using the fastest ramp that makes sense.

## The question

Four arms have now failed by the SAME route, and it is not a circuit event:
`I_cond` REVERSES SIGN while the electrode is still strongly negative
(-123 V, -282 V, -255 V in three of the four -- `V_el` never flipped at all in
those), at only 0.018-1.9% of the ballast's short-circuit current. What they
share is DENSITY, and therefore SCREENING:

| arm | n_e at reversal [m^-3] | lambda_D / gap |
|---|---|---|
| `grubert2009_ramp1us_300` | 5.33e16 | 0.0073 |
| `grubert2009_ramp1us_400` | 8.83e16 | 0.0057 |
| `grubert2009_step250` | 2.95e17 | 0.0029 |
| `grubert2009_floor1e9` | 2.93e17 | 0.0025 |

So the failure coincides with the plasma becoming able to screen -- the
**Townsend -> glow transition**, where a sheath must form and the electrode has
to collect ions through a thin space-charge layer instead of sitting in a
quasi-uniform field. **Is the model failing at that transition, or is the mesh
failing to resolve the sheath that forms there?**

## Baseline -- ABSOLUTE PATH, and a strict one-variable control

    /home/kkourtza/soplasma-scratch/validation/grubert2009_ramp1us_300

Same 1 us ramp to -300 V, same density floor (1e9 m^-3), same circuit
(`seriesResistor`, R = 1e8 ohm), same endTime. **The mesh is the only
difference:**

    baseline:  400 cells, Bump 0.02 -> GRADED, dx 1.35 um at the electrodes
                                      to 66.8 um at mid-gap, ratio 49
    this arm: 2000 cells, Bump 1.0  -> UNIFORM, dx = 5.0000 um everywhere
                                      (verified max/min = 1.0000)

Uniform removes the GRADING assumption entirely, which is what is in question:
`gap1cm.geo` grades toward the ELECTRODES, but the structures that formed in the
earlier runaway sat 215 um out and at mid-gap.

## Resolution actually delivered, and its limit

Cells per `lambda_D` at Te = 2 eV (the measured bulk value):

| n_e [m^-3] | lambda_D [um] | cells per lambda_D at dx = 5 um |
|---|---|---|
| 5.33e16 (where the baseline reversed) | 45.7 | **9.1** |
| 8.8e16 | 35.4 | 7.1 |
| 2.9e17 | 19.5 | **3.9** |
| 1e18 | 10.5 | 2.1 |

**So this mesh is well resolved AT THE ONSET (9 cells) and only ~4 cells by
3e17.** It is a clean reference for whether the reversal HAPPENS; it is NOT
converged for what follows it. Do not read the post-reversal trajectory as
benchmark-quality.

NX = 4000 (dx 2.5 um) would give 18 cells at onset and 7.8 at 3e17, but costs
~20.8 h against ~5.2 h to reach the reversal, and is only worth spending once
this arm has said whether the reversal survives refinement at all.

## Why the 1 us ramp to -300 V

Fastest route to the physics in question: on the graded mesh this exact
configuration reversed at **t = 1.121 us**, so the interesting regime arrives in
about a microsecond of simulated time. -400 V reversed marginally earlier
(0.875 us) but is a 43x jump in Townsend multiplication over the state known to
be stable, i.e. more extreme than the question needs.

NOTE the ramp-rate caveat established today: a FAST ramp reaches a given voltage
with LESS accumulated density (9.2x less at -200 V than a -10 V/us ramp), because
the discharge grows exponentially throughout the ramp. That is why this arm needs
-300 V rather than -200 V to reach the reversal at all.

## Success and failure, stated in advance

* **NO REVERSAL, or reversal at a much higher density** -> the graded mesh was
  under-resolving the sheath, the earlier failures are partly mesh artefacts, and
  the whole runaway story needs re-measuring on a resolved mesh. This is the
  outcome that would matter most.
* **REVERSAL AT THE SAME DENSITY (~5e16)** -> resolution is EXCLUDED as the
  cause, and the failure is in the physics or the closure at the Townsend->glow
  transition. Per the discipline `../RESOLUTION_TEST.md` pre-committed, mesh
  refinement is then RETRACTED as an explanation rather than pushed further.
* **REVERSAL EARLIER / at LOWER density** -> the graded mesh was accidentally
  stabilising, which would be worse than either.
* **dt COLLAPSE before reaching 5e16** -> the benchmark cannot reach the
  question; report the binding limiter and reconsider.

## Extraction

    ~/ct-env/bin/python ../../tools/settle_check.py .      # state, RUNAWAY flags, I/Isc
    ~/ct-env/bin/python ../../tools/news.py .              # n_e and Te evolution
    postProcessing/externalCircuit/circuit.csv             # I_cond SIGN is the onset signal

---

# RESULT 2026-09-07: the uniform mesh REVERSES TOO. Resolution is EXCLUDED.

Run in PARALLEL on 8 ranks after three parallel defects were fixed (rule 31,
`docs/design/parallel-blockers.md`); parallel verified against serial to 7 digits at
matched simulated times.

## Matched-time comparison against the graded mesh

| t [us] | GRADED n_e,max / Te,mx / I_cond | UNIFORM n_e,max / Te,mx / I_cond |
|---|---|---|
| 0.562 | 5.27e10 / 24.9 / -2.16e-10 | 4.91e10 / 25.1 / -2.00e-10 |
| 0.787 | 5.56e11 / 34.6 / -2.14e-09 | 5.27e11 / 35.3 / -2.02e-09 |
| 0.956 | 4.98e12 / 32.1 / -1.57e-08 | 4.62e12 / 34.6 / -1.49e-08 |
| 1.068 | 6.93e14 / 11.6 / -7.60e-08 | 2.28e14 / 11.6 / -6.85e-08 |
| **1.124** | 2.15e17 / 17.1 / **+1.66e-06** | 7.22e16 / 14.3 / **+9.16e-07** |

`I_cond` NEGATIVE is healthy; POSITIVE is the reversal that has killed every arm.

**Graded reversed at t = 1.121 us. Uniform reversed at t = 1.1241 us.** The two
agree to **5-7%** on `n_e,max` and within 8% on `Te,max` up to t ~ 0.96 us, then
diverge ~3x through the transition -- expected, since exponential growth
amplifies any small difference -- and reverse at the SAME MOMENT.

## The verdict, against the criteria registered before the run

This is the **"REVERSAL AT THE SAME DENSITY"** outcome: resolution is EXCLUDED.
A 5x-refined UNIFORM mesh -- no grading assumption at all, dx = 5.0000 um
everywhere, **9 cells per lambda_D at the reversal density** -- does the same
thing at the same time as the graded 400-cell mesh.

Per the discipline `../RESOLUTION_TEST.md` pre-committed on 2026-09-06,
**mesh refinement is therefore RETRACTED as an explanation, not pushed to a
further level.** NX = 4000 is not worth running for this question.

**The current reversal at the Townsend->glow transition is not a mesh artefact.**
It is in the physics or in the closure.

## What the run also SHOWED, and it is the most interesting part

The temperature profile FLATTENS as the transition begins:

| t [us] | V_src [V] | n_e,max [m^-3] | Te,max [eV] | Te,min [eV] |
|---|---|---|---|---|
| 0.911 | -273.2 | 2.40e12 | **37.31** | 2.53 |
| 0.979 | -293.7 | 6.68e12 | 31.22 | 2.60 |
| 1.103 | -300.0 | 4.82e15 | **12.14** | 2.95 |
| 1.119 | -300.0 | 1.72e16 | 13.26 | **5.49** |

`Te,max` PEAKS at 37.3 eV and COLLAPSES to 12 eV while the bulk `Te,min` rises
2.5 -> 5.5 eV. That is the field collapsing in the fall as space charge finally
screens it -- the Townsend->glow transition actually occurring. **The model gets
INTO the transition and then loses the conduction current's sign.** So the
failure is specifically at the moment the sheath forms, not before it.

## Where this leaves the question

Excluded so far, each by measurement: the CIRCUIT (reversal happens at
0.018-1.9% of I_sc, with V_el still strongly negative), the DENSITY FLOOR
(1e9 and 1e7 agree; it postponed but did not remove the failure), the ENERGY
TABLES (self-consistent to 0.997-1.020 over five decades once the growth term is
included), and now RESOLUTION and GRADING.

What remains is the wall-flux closure at a forming sheath and the `n_e`/`nEps`
coupling -- the latter already implicated twice independently, and now a third
time, since the growth-term cancellation carries 99% of the power budget exactly
where `Te,max` collapses.
