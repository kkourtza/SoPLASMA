# dt should be chosen by the solver, not by the user

Status: DESIGN + step-1 measurement, 2026-09-06. Nothing implemented.

## The complaint, and it is correct

> "we still let Co of meanE to be a user defined value while we dont even know
> what it should be for every case? ... the solver should be able to run
> automatically to the biggest possible dt (while being still robust and
> accurate)"

Three controllers run simultaneously on every step and disagree by FIVE ORDERS
OF MAGNITUDE. Measured on `grubert2009_iset`, ~16000 steps:

| signal | what it wants | vs actual dt | how often it BINDS |
|---|---|---|---|
| accuracy -- temporal error PI, `\|\|e\|\| = 2.35e-09` vs target 1 | 9.85e-09 | **752x larger** | 18 steps |
| **Courant (energy)** -- `maxEnergyConvectiveCo`, HAND-SET to 15 | 1.31e-11 | -- | **15723 steps** |
| robustness -- outer contraction `rho = 6.23` vs target 0.912 | 1.92e-12 | 6.8x smaller | never (`diagnostic; omega still governs`) |

The accuracy controller wants 752x more, the robustness diagnostic wants 6.8x
less, and **the number that decides is the one a user typed.** That is the G1
violation, exactly as stated.

## The universality requirement (user, 2026-09-06)

> "we need to make sure that these fixes are UNIVERSAL and should work for any
> case i.e. they should be based on metrics that scale/translate"

So every signal that governs dt must be a DIMENSIONLESS RATIO whose threshold
means the same thing in any case. Auditing the candidates against that:

| metric | form | dimensionless? | is its THRESHOLD universal? |
|---|---|---|---|
| `\|\|e\|\|` / target | relative temporal error | yes | **yes** -- "1" means "at the accuracy you asked for" |
| `rho` | `rNorm/rNormPrev` | yes (`plasmaOuterRelaxation.C:448`) | **yes in principle** -- "<1" means the loop contracts |
| discard fraction | discarded/total steps | yes | **yes** -- "0" means no step failed |
| Courant number | `u dt/dx` | yes | **NO.** This is the crux |

A Courant number is dimensionless but its ACCEPTABLE VALUE is not universal:
these equations are solved implicitly, so Courant is not a stability bound, and
the "right" cap depends on the scheme, the stiffness and the case. That is why
`maxSpeciesCo = 15` is a number nobody can derive -- and why it must not be the
thing that decides dt.

## STEP 1 RESULT: rho is NOT usable, and the rejection path is UNTESTED

Two findings, both negative, both necessary.

**1. `rho > 1` fires on demonstrably healthy runs.** Correlating median `rho`
against actual discarded steps over 12 cases and ~1.35 MILLION steps:

| case | steps | DISCARDED | median rho | verdict |
|---|---|---|---|---|
| grubert2009 | 6094 | 0 | **378.2** | "residual GREW" |
| grubert2009_iset | 22535 | 0 | **22.05** | "residual GREW" |
| grubert2009_LFA | 9789 | 0 | **4.82** | "residual GREW" |
| grubert2009_fix_n11 | 238819 | 0 | 0.194 | contracting |
| grubert2009_fix_n8 | 383932 | 0 | 0.191 | contracting |
| grubert2009_fix_fast | 26424 | 0 | 0.193 | contracting |
| ballastCo_1000 | 30990 | 0 | 0.327 | contracting |
| ballastCo_10000 | 30820 | 1 | 0.248 | contracting |

`grubert2009` reports `rho = 378` -- a residual growing 378x per corrector --
with ZERO discarded steps. And `iset`, whose pre-ignition trajectory was
verified against the analytic `I_set/C_gap` slope to **1.33%** and which then
ignited within **4%** of the predicted voltage, reports `rho = 22`. **A signal
that screams on a run that is provably accurate cannot govern dt.** `rho` must
be fixed or retired before it is used for anything; it is currently the
silent-diagnostic trap inverted -- it is not silent, it is wrong.

**2. The rejection feedback has never been exercised: 1 discarded step in
~1.35 million.** So there is NO evidence that the discard-and-halve loop works
as a robustness limiter, because dt has always been held ~750x below what
accuracy needs and the outer loop has never been stressed. "Zero rejections" is
not evidence of robustness at large dt -- it is evidence that large dt was
never tried.

    (And the count itself needed fixing: the marker is `DISCARDING this step`,
    not the `retryStep|REJECT` pattern used earlier in the day, which matched
    nothing. The corrected count agrees, so the earlier reports were right by
    luck rather than by measurement.)

## The design that follows

    dt = min( dt_accuracy , dt_robustness , dt_backstop )

* **`dt_accuracy`** -- the existing temporal error controller. Universal: the
  user sets an ERROR TARGET, which is a physics-level choice they can reason
  about, not a Courant number nobody can derive.
* **`dt_robustness`** -- learned from OBSERVED discards, not predicted. dt
  grows while steps succeed and is cut when one is discarded, so the controller
  finds the robustness boundary by walking into it. Universal because "a step
  failed" means the same thing everywhere.
* **`dt_backstop`** -- the Courant caps, demoted, at a deliberately LOOSE value
  whose only job is to stop a pathological first step. Not a tuning knob.

**And the robustness path MUST BE EXERCISED FIRST (rule 23).** A control
mechanism that has fired once in 1.35 million steps is as untested as a guard
that has never fired. Before it can be trusted to bound dt, push the caps up
until discards appear, and confirm the loop recovers and the answer is
unchanged.

## Staged plan

1. **DONE -- step 1:** `rho` is not usable; the rejection path is untested.
2. **Exercise the robustness boundary.** Raise the caps progressively on
   `iset` until discarded steps appear. That is the same walk the automatic
   controller would perform, so it validates the mechanism AND finds the real
   boundary. The `iset_En150` / `iset_En1500` arms are the first two points:
   both run at ~5x the baseline dt with zero discards, so the boundary is above
   5x, and 150 vs 1500 is indistinguishable because the SPECIES cap takes over.
3. **Diagnose `rho`.** Either it is measuring the wrong residual, or its
   normalisation is wrong, or "growing" is not a failure for this loop. Until
   this is settled `rho` stays a printed diagnostic and governs nothing.
4. **Make accuracy primary**, robustness the limiter, Courant the backstop.
5. **Remove `maxSpeciesCo` / `maxEnergyConvectiveCo` from the case layer** --
   the user sets an error target and nothing else. That is the G1 deliverable.

## Provenance

Measured 2026-09-06 from the case logs named above. All timing numbers in this
document post-date `b6a7102`; earlier timings are void because 45-50% of
timesteps ran zero correctors.
