# What EVERY case must monitor, and the order to debug in

Status: DESIGN, 2026-09-06. Written after a day lost to a bug that one
consistency check would have caught in the first minute.

## The failure this exists to prevent

`updateChargeDensity()` stopped being called (a raw `break` bypassed the
`pimple.finalIter()` guard it sat behind). Poisson's source froze. Measured:
**151 updates in 364670 timesteps -- 0.04%.**

The day then went into: wall-flux boundary forms, Scharfetter-Gummel `f`,
ion mobility tables, the `meanEnergyMax` clamp, table ranges, ballast RC
timescales, `currentSource` design, load-line algebra, the ion wall-flux sign,
Sato per-species currents. **Every one of those was a plausible physics
hypothesis and every one was wrong**, because the field was not a physics
result at all -- it was the vacuum solution of an equation with no source.

Two facts were sitting in plain view the whole time:

1. `E/N = 1104 Td` measured against `V/L = 1106 Td` -- the field equalled the
   VACUUM field to 0.2%, in a run with three decades of accumulated space
   charge. That is impossible, and I read it as physics.
2. `chargeDensity` was flat at ~4e-8 C/m^3 across the whole gap while
   `e(n_Arp - n_e)` varied over three decades. **A single plot of those two
   together ends the investigation immediately.**

Neither was visible without work: the first needed an arithmetic check nobody
was doing, the second needed parsing a **1.15 GB log** with a bespoke script.

## The distinction that organises everything below

**INVARIANTS** must hold BY CONSTRUCTION. A violation is a bug, never a
result, so it is checked automatically and reported loudly.

**OBSERVABLES** are physics. They are written as a time series so a human or
an agent can look at them without parsing a log.

Conflating the two is why this was missed: the min/max ranges in the log are
observables, and no amount of staring at observables tells you the source term
is dead. Only an invariant does.

## Part A -- INVARIANTS (the high-value half)

Checked every `checkInterval` steps (default ~100: all are O(nCells) and cost
nothing next to a Poisson solve). A violation prints ONCE, loudly, with the
cell, the two values, and the ratio -- and says which is the suspect.

| invariant | catches | today? |
|---|---|---|
| `chargeDensity == sum(q_i n_i)` | derived field not tracking its sources | **YES, immediately** |
| `-div(eps grad V) == chargeDensity` | Poisson solved a different problem than posed | yes |
| **per-step update COUNTERS** | any per-step state update that stopped being reached | **YES, immediately** |
| every transported field has its old-time levels | step-discard / error-control breakage | no |
| `sum(species wall flux) == I_collected` | wall-flux vs current-diagnostic disagreement | no |
| `n >= 0` and `n >= floor` per species | a clamp silently doing the physics | no |

**The counters are the most valuable new mechanism**, because they catch the
CLASS rather than the instance. Every per-step state update registers itself;
at the end of each accepted step the solver asserts each ran at least once.
Today's bug becomes:

    *** updateChargeDensity() ran 0 times in step 41273 (expected >= 1) ***
        Poisson's source is STALE. Something exited the corrector loop
        before the update. See doc/case-monitoring-plan.md.

That is a one-line, self-explaining failure instead of a day.

**Design rule for the checks themselves:** each must be exercised on an input
where it MUST fire (rule 23) before it is trusted. A silent invariant is
indistinguishable from an absent one.

## Part B -- OBSERVABLES: one `history.csv` per run

Written every `historyInterval` steps (default: every step up to a cap, then
throttled) to `postProcessing/history/history.csv`. **No case should ever
require parsing the log to get a time series again.**

Columns, using the framework's existing names (G3):

* `time`, `deltaT`, **`deltaT_limiter`** (which limiter set it -- already
  computed, currently only printed), `outerIters`, `rejectedSteps`
* per species: `n_<sp>_min`, `n_<sp>_max`, `n_<sp>_mean`
* `chargeDensity_min/max` **-- next to the densities, deliberately.** This
  adjacency is what makes the flat-vs-growing failure visible at a glance.
* `Emag_min/max`, `reducedE_max` in **Td**
* `meanE_min/max` (LMEA) or `Te` (LFA), plus `meanE_clampedCells`
* `I_total`, `I_cond`, `I_disp` and per-species currents (already exist)
* `V_electrode` per driven/circuit electrode (already exists for circuits)
* `surfCharge_min/max` where a dielectric exists

Plus `tools/plot_history.py` producing the standard six-panel figure, and
emitting the data it plots -- per rule 6, the plotting script is separate so
the figure can be changed without re-running anything.

**A derived quantity goes NEXT TO its sources in the CSV, always.** That is
the cheap structural trick that turns an invisible bug into an obvious one.

## Part C -- THE DEBUGGING ORDER (the actual rule)

When a case misbehaves, in this order, and **do not skip to physics**:

1. **Did every per-step update run?** Counters, or `grep -c` against the step
   count. Cost: seconds.
2. **Are derived fields consistent with their sources?** `chargeDensity` vs
   `sum(q n)`; transport coefficients vs their lookup key; `meanE` vs
   `n_eps/n_e`. Cost: seconds.
3. **Is any field sitting exactly on an analytic no-op?** The vacuum field, a
   floor, a clamp, a uniform initial value. A field that equals its own
   trivial solution is a broken term, not a result. Cost: one arithmetic check
   -- and it is the one I had already done and misread.
4. **Are the limiters and clamps binding?** `deltaT_limiter`,
   `meanE_clampedCells`, floor-hit counts. A guard that binds constantly means
   the run is outside the model's validity.
5. **Only now**: question the physics, the closure, the coefficients.

Steps 1-3 together cost under a minute and would have ended today's
investigation before it started.

## Proposed CLAUDE.md rules, for the user's approval

CLAUDE.md is the authority and its rule numbers are stable identifiers, so
these are proposed rather than added:

> **26) EVERY case writes a machine-readable time series** --
> `postProcessing/history/history.csv` -- carrying at minimum: species
> densities, `chargeDensity`, `Emag`, the electron energy, the currents, `dt`
> and WHICH LIMITER SET IT, outer iterations and rejected steps. A derived
> quantity sits NEXT TO the sources it is derived from. If a diagnosis needs
> the log parsed, the monitoring is the defect. Measured 2026-09-06: a frozen
> Poisson source ran for 364670 steps and needed a 1.15 GB log parsed with a
> bespoke script to see.

> **27) INVARIANTS BEFORE PHYSICS.** A derived field must be checked against
> its own sources, and every per-step update against a call counter, BEFORE
> any physical explanation is entertained. A field equal to its own trivial
> solution (the vacuum field, a floor, a clamp) is a broken term, not a
> result. Measured 2026-09-06: `E/N = 1104 Td` against `V/L = 1106 Td` -- the
> vacuum field to 0.2% -- was read as physics, and a day went into wall
> fluxes, mobility tables, energy clamps and circuit timescales while
> `updateChargeDensity()` had run on 0.04% of steps.

## Cost, honestly

Part A is cheap and is where the value is: O(nCells) arithmetic at a throttled
cadence, plus a counter registry. Part B is a CSV writer over quantities the
solver already computes -- the only real work is agreeing the column names
once, which G3 requires anyway. Part C is free; it is a habit.

The risk to avoid: monitoring so much that nobody reads it. Hence the split --
invariants SHOUT and need no reader, observables are written for plotting and
need no attention until something is wrong.
