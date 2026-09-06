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

## STEP 2 RESULT (2026-09-06): the caps are NOT merely masking the governor

Ran `iset` (all caps 15) against four arms: energy cap alone at 150 and 1500,
and ALL caps at 150 and 1500. Extended into the ignition transient at the
user's suggestion -- "run the simulations even longer to get error estimates in
more interesting plasma regimes" -- which is what produced the real answer.

**In the quiescent phase everything looks free.** Errors of 0.05-0.27% in
V_gap for 15-30x of dt. That measurement is a trap: nothing is happening.

**In the ignition transient the arms separate sharply:**

| arm | dt vs base | err in V_gap (transient) | outcome |
|---|---|---|---|
| energy cap 150 / 1500 | 1.5x | **0.35%** at 160% overshoot | tracks the baseline |
| ALL caps 150 | 15x | -- | **n_e -> 2.7e17, diverged** |
| ALL caps 1500 | 30x | -- | **n_e -> 1.2e18, diverged** |

`En150` and the baseline agree to **0.35%** at EVERY common time right through
ignition -- 162.9% vs 162.3% of I_set, V -179.52 vs -180.14. They are the same
solution. So **the ENERGY cap is an efficiency knob even in the hard regime**,
and it buys 1.5x because the SPECIES cap immediately takes over.

**But raising ALL the caps diverged.** And the mechanism matters: the
discharge has a REAL current overshoot at ignition -- the baseline reaches
160-200% of I_set too -- which the current source then corrects by pulling
V down. **The coarse arms did not diverge because the caps were preventing an
artefact; they diverged because they could not RESOLVE a real overshoot well
enough to recover from it.**

That reframes the whole question. The Courant caps are not redundant with the
governor. They are what keeps the outer loop resolved enough that the
robustness and accuracy controllers can act in time.

### PROVISIONAL -- THE PREMISE IS NOT YET ESTABLISHED (flagged 2026-09-06)

**The paragraph above depends on the BASELINE RECOVERING from its overshoot,
and that is not yet measured.** At the time of writing `iset` is at
t = 8.4e-7 with `Ic/Is` = 158% and still rising, while the coarse arms diverged
somewhere past t = 9.5e-7 -- a time the baseline has not reached. There is
therefore NO CONTROL at the moment that matters (rule 15).

Two readings remain open:

* **if the baseline recovers** (Ic/Is turns back toward 100%, n_e settles near
  1e15) -- then the coarse arms failed to resolve a recoverable overshoot, the
  reasoning above holds, and the growth-rate proposal follows.
* **if the baseline ALSO diverges** to ~1e18 -- then the coarse arms were
  simply arriving at the same place sooner, the caps prevented nothing, and the
  reasoning above is BACKWARDS. The instability would then be physical and the
  question returns to the ramp rate and the circuit.

**Do not build on this section until `iset` has passed t = 9.5e-7.** The
0.35% agreement between the baseline and the energy-cap arms is measured at
common times and stands independently; the divergence reasoning does not.

**And the accuracy controller reacted TOO LATE.** At All1500 it eventually
clamped dt by 400x and forced 11 discards -- after the solution had left.
`||e|| = 2.35e-09` against a target of 1 reported "everything is fine" while
the trajectory was departing. **That is the defect to fix for a universal
governor: the error measure is a LOCAL step error and is blind to a physical
instability that is about to run away.**

### What a regime-universal governor therefore needs

The user's requirement is that it work for "low pressure glows, sheaths,
cathode layers, streamers (high pressure), surface ionization waves". From this
measurement, three signals are needed, and only two exist:

1. **local truncation error** -- exists (temporal error PI), dimensionless
   against its target. NECESSARY, NOT SUFFICIENT: blind to the instability.
2. **outer-coupling margin** -- exists (`coupling margin`, keyed on the Aitken
   omega, dimensionless). It backs dt off before a step fails, which is why
   discards are near-nonexistent.
3. **A GROWTH-RATE signal, which DOES NOT EXIST.** The instability that
   defeated the coarse arms has a physical rate -- the ionisation e-folding
   `gamma = nu_iz - nu_loss`, measured at 0.400 ns here. A dt that resolves the
   local error to 1e-9 but is comparable to `1/gamma` cannot follow the
   physics. `dt * gamma << 1` is DIMENSIONLESS and is the same requirement for
   a glow, a streamer or a surface wave -- it is exactly the regime-independent
   criterion the Courant caps are a crude proxy for.

**That is the proposal: replace the hand-set Courant caps with `dt*gamma`**,
where gamma is the fastest local growth rate the chemistry itself reports.
`nu_iz` is already computed every step for the source terms, so the signal is
available and costs nothing. Unlike a Courant number its threshold is
physical -- "resolve the fastest growth" -- and therefore transfers between
regimes without retuning.

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
