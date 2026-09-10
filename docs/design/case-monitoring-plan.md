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
        before the update. See docs/design/case-monitoring-plan.md.

That is a one-line, self-explaining failure instead of a day.

**Design rule for the checks themselves:** each must be exercised on an input
where it MUST fire (rule 23) before it is trusted. A silent invariant is
indistinguishable from an absent one.

## Part B -- OBSERVABLES: PROBES, not fields, not just min/max

**Written at PROBE POINTS.** One to three per case. This is the user's call
(2026-09-06) and it is the right one: a full field history is unwritable in
2-D and unthinkable in 3-D, and domain min/max -- which is what the log already
prints -- throws away WHERE, which is usually the answer. Today the ion peak
sat at x = 9.63 mm and the electron peak at 9.26 mm, and reading their maxima
as a ratio produced a "143x non-neutrality" that does not exist.

**Probe positions are DERIVED by default, per G1.** For a plane-parallel gap
the electrode patches already fix the axis, so the default is three probes at
5%, 50% and 95% of the gap -- near-cathode (in the fall), mid-gap (negative
glow), near-anode. A case overrides with explicit coordinates when it knows
better; it should never be *required* to supply them.

    probes                     // optional; derived if absent
    (
        (0.0005 0.0001 0.0001)   // cathode fall
        (0.005  0.0001 0.0001)   // mid-gap
        (0.0095 0.0001 0.0001)   // near anode
    );

`postProcessing/history/history.csv`, one row per written step:

* run-level: `time`, `deltaT`, **`deltaT_limiter`** (which limiter set it --
  already computed, currently only printed), `outerIters`, `rejectedSteps`
* per probe `k`: `n_<sp>_p<k>` for every species, `chargeDensity_p<k>`,
  `Emag_p<k>`, `reducedE_p<k>` in **Td**, `meanE_p<k>` (LMEA) or `Te_p<k>`
* domain reductions KEPT alongside, because they are what catch a field that
  has gone uniform -- the exact signature of today's bug:
  `Emag_min/max`, `n_<sp>_max`, `chargeDensity_min/max`,
  `meanE_clampedCells`
* currents and electrode potentials, which already exist:
  `I_total`, `I_cond`, `I_disp`, `I_<sp>`, `V_electrode`

**A DERIVED quantity sits NEXT TO its sources**, in the same file and adjacent
columns. `chargeDensity_p1` beside `n_e_p1` and `n_Arp_p1` makes "flat while
they grow three decades" visible at a glance. This adjacency is the cheapest
debugging tool available and it is not an accident of column order.

Plus `tools/plot_history.py`, separate per rule 6, so the figure can be
changed without re-running anything.

Note there IS an OpenFOAM `probes` functionObject, and it is not enough on its
own: it scatters one file per field across time directories, cannot carry
`deltaT_limiter` or the rejection count, and gives no control over adjacency.
It is worth reusing for the interpolation, not for the output.

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

**ADDED to CLAUDE.md on 2026-09-06 as rules 26 and 27** -- 26 under "Records
that outlive the conversation", 27 under "Evidence and claims". CLAUDE.md is
the authority; read it there. Rule 26 as adopted requires PROBES, which is
stronger than the draft this document first carried.

## Cost, honestly

Part A is cheap and is where the value is: O(nCells) arithmetic at a throttled
cadence, plus a counter registry. Part B is a CSV writer over quantities the
solver already computes -- the only real work is agreeing the column names
once, which G3 requires anyway. Part C is free; it is a habit.

The risk to avoid: monitoring so much that nobody reads it. Hence the split --
invariants SHOUT and need no reader, observables are written for plotting and
need no attention until something is wrong.
