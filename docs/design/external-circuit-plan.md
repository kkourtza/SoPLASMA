# External circuit coupling — staged plan

Written 2026-09-06, after the Grubert dc glow showed why the capability is
needed. Stages 1 and 2a are IMPLEMENTED, stage 3 is half done,
stage 4 was decided against; the rest is design.

## Why the framework needs this at all

A gap held at a FIXED applied voltage above breakdown has nothing to limit the
current it draws. Measured on the Grubert case: the gap broke down at -179 V
(correct — Paschen for argon at pd = 0.75 Torr cm sits near its minimum), and
the current then went 0.51 → 6.7 → 40 → 59 mA/cm² and was still climbing when
the timestep collapsed. It passed straight through the 0.511 mA/cm² the
reference reports.

The reason is a property of the discharge, not of the solver. Over that stretch
the current TRIPLED while the applied voltage moved 0.1 V:

    V = -182.29  j = 40.2      V = -182.37  j = 52.2      V = -182.40  j = 58.9

That is the flat (normal-glow) part of the V-I characteristic, where dV/dj ≈ 0
and the current at fixed voltage is essentially undetermined. Grubert's own
operating point is on the ABNORMAL branch at -500 V, where dV/dj > 0 and the
discharge IS stable at fixed voltage — but a transient cannot get there,
because it must first traverse the flat region, and it runs away in it.

Grubert avoid this entirely by solving the TIME-INDEPENDENT equations: their
operating point is computed as a boundary-value problem, never integrated
towards. An unstable equilibrium is still a perfectly good solution of the
steady equations. A transient solver does not have that luxury and needs the
circuit the experiment actually has.

**Related capability, tracked separately: a STEADY-STATE solver.** It reaches
these operating points directly and is the more faithful reproduction of
published fluid-model results. OpenFOAM's steady solvers make it tractable.
That is not a substitute for this plan — a real discharge has a real circuit,
and pulsed/RF cases are transient by nature — but the two together cover the
field properly.

## Stage 1 — series resistor (IMPLEMENTED)

`plasmaExternalCircuit`, `externalCircuit` in `system/plasmaSimulationControls`:

    externalCircuit
    {
        type          seriesResistor;
        electrode     cathode;
        sourceVoltage table ((0 0) (5e-6 -1011));
        resistance    5e8;      // [Ohm]
        relaxation    1.0;      // optional
    }

    V_electrode(t) = V_source(t) - R * I_circuit

`I_circuit` is Sato's `I_total` from `plasmaDischargeCurrent` — by construction
the current an external circuit carries — passed IN rather than looked up, so
the circuit cannot grow a second, independently-derived current.

KNOWN LIMITATION, stated rather than hidden: the coupling is EXPLICIT with one
step of lag (the current is measured after the step that produced it), so a
large R is stiff and `relaxation` exists to damp it. Stage 3 removes this.

## Stage 2 — the lumped models worth owning

Same interface, more topologies. These cover the overwhelming majority of
laboratory plasma sources and are a few hundred lines in total:

| type | for |
|---|---|
| `seriesRC` | ballast with cable/blocking capacitance; DBD driving circuits |
| `seriesRLC` | pulsed sources, ringing, real supply inductance |
| `currentSource` | current-driven operation — pins j directly, which is exactly what a V-I characteristic sweep wants |
| `matchedRF` | 13.56 MHz CCP with a matching network; the standard capacitively-coupled configuration |
| `dcWithBlockingCap` | self-bias, which needs the DC-blocking capacitor to develop |

`currentSource` deserves priority: it makes the V-I characteristic directly
computable, and the characteristic is what tells you whether a discharge is
normal, abnormal or arcing — the very distinction that explains the Grubert
runaway.

## Stage 2a — `currentSource` (IMPLEMENTED 2026-09-06)

### What a current source is, physically

A ballast makes the operating point self-selecting but leaves it to be FOUND.
A current source PINS it: the external circuit forces a prescribed current
through the electrode and lets the voltage be whatever the discharge needs.
Physically it is a supply with output impedance large compared with the
discharge's own |dV/dj| -- the R -> infinity limit of a ballast, which is why a
real lab "current-regulated" supply is just a very stiff ballast with feedback.

The gap voltage is then a PREDICTION of the simulation rather than an input,
and the discharge sits wherever its characteristic says it must for that
current.

### Why it is the better instrument for validating against a paper

Published discharge studies almost always report the current density. Grubert
report j = 0.511 mA/cm^2 at a gap voltage of -500 V. With a current source we
impose the half we are most confident of and PREDICT the other half:

    impose j = 0.511 mA/cm^2   ->   does V_gap come out at -500 V?

That has NO free parameter. The ballast route does not have that property: R
must be chosen, and choosing it so the gap lands at -500 V when j = 0.511
ASSUMES the answer. Measured 2026-09-06, that is exactly what was done here,
and it is circular in a way a validation must not be.

### Why it is also far more ROBUST numerically

This is the argument that decided the priority. A ballast closes a feedback
loop: V sets I, I sets V. That loop is what destabilised the runs --
post-breakdown the plasma current jumped to 8.5x what the circuit could supply
at that instant, because the current the circuit sees is one step old, and the
gap voltage was driven through zero to a positive cathode.

A current source has NO SUCH LOOP. The current is imposed, so there is nothing
to iterate: the instability being fought does not exist. It sidesteps stage 3
rather than depending on it.

### Steady versus non-steady behaviour, which differ and must be documented

* **Steady / quasi-steady.** The natural mode. The discharge relaxes to the
  point on its characteristic corresponding to the imposed j, and the gap
  voltage settles. This is the mode that reproduces a published operating
  point.
* **Ramped current.** Sweeping j slowly traces the V-I CHARACTERISTIC directly
  -- Townsend, subnormal, normal, abnormal -- which is the single most
  informative diagnostic of a dc discharge and is not obtainable at fixed
  voltage at all, because the flat normal-glow branch is not a function there.
  This is worth having for its own sake.
* **Before ignition, a current source is UNPHYSICAL and must be handled.**
  A cold gap cannot carry an imposed current: the voltage required is unbounded
  and the model will chase it to breakdown-and-beyond. Two admissible
  treatments, and the case must state which:
    - `softStart`: ramp the imposed current from ~0, so the gap voltage rises
      with it and ignition happens on the way up. Physically this is a supply
      in voltage-limited mode until the discharge lights.
    - `compliance`: give the source a voltage limit V_max, exactly as a real
      supply has. Below ignition it behaves as a fixed-voltage source at V_max;
      once the discharge can carry the current, it switches to current
      regulation. This is what a laboratory supply actually does and is the
      recommended default.
  A current source WITHOUT a compliance limit is a modelling error before
  ignition, and the implementation must refuse it rather than produce a
  spectacular transient.

  **REFINED BY THE IMPLEMENTATION, 2026-09-06.** The two treatments above were
  framed as alternatives the case must choose between. They are not: keeping
  the `C dV/dt` term in the regulated equation makes the pre-ignition phase
  physical *by itself*, with no mode switch and no `softStart` key. With no
  plasma, `dV = dt I_set / C` -- the source charges the electrode capacitance
  at constant current, which is what a real current-limited supply does into a
  capacitor, and the gap voltage RAMPS to breakdown on its own at a rate the
  user controls through C. `softStart` was therefore NOT implemented: it would
  be a hand-written ramp standing in for one the circuit already produces,
  which is precisely what G1 forbids.

  What DID survive is `compliance`, and it is required as this section says --
  but for a sharper reason than "below ignition it behaves as a fixed-voltage
  source at V_max". It never behaves as a fixed-voltage source; it behaves as a
  constant-current charger, and the rail is what stops the ramp if the gas
  never breaks down at all.
* **Fast transients (pulsed, RF).** A current source is the wrong instrument:
  there the circuit's own dynamics ARE the physics. Use seriesRC/RLC.

### THE DISCRETISATION AS IMPLEMENTED

A current source regulates the TERMINAL current, which is conduction plus the
displacement drawn by the electrode capacitance:

    I_set = I_cond + C dV/dt                                            (1)

Backward Euler on (1), with the plasma linearised by the same secant
conductance `g = dI/dV` the ballast uses, gives ONE update that is correct in
both regimes with no switch between them:

    dV = dt (I_set - I_cond) / (C + |g| dt)                             (2)

**Before ignition** `g -> 0` and (2) becomes `dV = dt I_set / C`: the source
charges the electrode capacitance at constant current, `dV/dt = I_set/C`. That
is exactly what a real current-limited supply does into a capacitor, and it
means **the pre-ignition voltage ramp is a consequence of the circuit rather
than a ramp anybody writes**. Its rate is a design knob set by C alone -- which
is the same lever that made the ballasted case break down gently.

**After ignition** `g` is large and (2) becomes `dV = (I_set - I_cond)/|g|`, a
damped Newton step on `I(V)` whose fixed point is `I_cond = I_set` exactly, for
any `g`.

`|g|` rather than `g`, for the reason documented on the ballast: `dI/dV` is
genuinely negative through a glow's negative differential resistance, and a
signed `g` would inflate the step and invert its direction precisely where the
discharge is stiffest.

The **compliance rail** is clamped on both sides. Above it the supply cannot
go. Below zero it would have to reverse polarity, which a single-quadrant
supply cannot do either -- and an unclamped regulator does try to, because
overshooting the set current asks for a voltage of the other sign.

### IT IS THE LIMIT OF THE BALLAST SWEEP, which is how it was validated first

A ballast's short-circuit current, with the load line pinned through the
operating point (`V_src = V_gap + R I_op`), is

    I_sc/I_op = 1 + V_gap/(R I_op)

so a ballast only pins the current as `R -> infinity` -- and that limit IS this
model. Measured on the Grubert case, 2026-09-06: 5.89x the operating point at
R = 1e8, 1.49x at 1e9, 1.10x at 5e9. A current source is the 1.00x end of the
same sweep, which is why the R sweep in `validation/grubert2009_R*/COMPARE.md`
tests this model's physics before the model is used.

### Binding it in the generator, so the user writes one thing

Same route as the ballast, on the electrode it feeds:

    cathode
    {
        kind        currentDrivenElectrode;
        circuit
        {
            type        currentSource;
            setCurrent  1.022e-6;   // [A] MAGNITUDE, or a Function1 to sweep
            compliance  -1500;      // [V] supply rail, SIGNED, REQUIRED
            capacitance 5e-14;      // [F] optional; sets dV/dt = I_set/C
        }
    }

and nothing else: the electrode is the patch it is written on,
`plasmaCreateSpeciesFields` emits `fixedValue` for it as now, and the solver
reads it through `boundaryRoleLibrary::caseDeclaration`. Current in AMPERES,
not a density -- the electrode area is a property of the mesh and deriving j
from it is the solver's job, not the user's.

**TWO ROLES, not one, and `ballastedElectrode` KEEPS ITS NAME.** This resolves
the open question this section used to carry, and it reverses a rename I made
and then withdrew the same day.

The rename argument was that the role should say WHAT the surface is, and that
`ballasted` named one particular circuit. That is half right and reached the
wrong conclusion. What the surface *is*, physically, includes **which quantity
the supply imposes** -- and that is not a topology detail, it is the same axis
as driven / grounded / floating. A voltage source through a ballast and a
current source are genuinely different instruments:

| role | you set | you find out |
|---|---|---|
| `ballastedElectrode` | supply voltage, R, C | the current, and the gap voltage |
| `currentDrivenElectrode` | the current, and the rail | the gap voltage |

`ballastedElectrode` is the DEFAULT to reach for, because **the operating
current is usually what the case is trying to determine and is not known at the
start.** Collapsing both into one generic role would have hidden that choice
behind a `type` key and taken it away from the user, against G1's "require an
input only for genuine physics the user alone can know" -- this is exactly such
an input.

A `type` that disagrees with the `kind` is FATAL in both directions, and so is
an unknown `kind`: an unrecognised `kind` here falls through to a plain
fixed-voltage electrode, the exact failure this whole class exists to remove.

**The sign convention, stated once.** `sourceVoltage` and `compliance` are
SIGNED (negative for a cathode); `setCurrent` is a MAGNITUDE whose polarity
comes from the sign of `compliance`. This is deliberate: `I_cond` is negative
on a negative electrode -- verified against
`postProcessing/externalCircuit/circuit.csv` on 2026-09-06, not assumed -- and
that is the one convention a user cannot be expected to guess. Both a negative
`setCurrent` and a `resistance` on a `currentSource` are fatal with a message
saying why.

### MEASURED DEFECT 2026-09-06: the regulator is an UNDAMPED INTEGRATOR

First real use of `currentSource` (`validation/grubert2009_iset`, argon glow,
`setCurrent` ramped 1e-8 -> 1.022e-6 A over 10 us, `capacitance 0`). It ignited
correctly -- pre-ignition `dV/dt` matched `I_set/C_gap` to **1.33%** and
ignition occurred at V = -116 V against a predicted -121 V (**4%**) -- and then
went unstable in a way that is entirely the controller's doing.

**THE OSCILLATION, measured:**

| t | Ic/Is | V |
|---|---|---|
| 8.59e-7 | 178% | -175 |
| 9.19e-7 | **229%** (peak) | -139 |
| 9.42e-7 | 69% | -127 |
| 9.44e-7 | **0%, sign reversal** | -127 |
| 9.51e-7 | -1064% | -146 |
| 9.54e-7 | **-21175%** | -214 |

Overshoot to 229% -> the source correctly pulls V from -175 to -127 -> the
conduction current COLLAPSES SMOOTHLY THROUGH ZERO (verified continuous:
-1.19e-9, -8.9e-10, -5.9e-10, -2.9e-10, **+1.3e-11**, +3.2e-10, so not a
diagnostic glitch) -> `|I_cond| < |I_set|` now, so the integrator correctly
drives V back negative -> the second swing re-ignites into a runaway that
never recovers, n_e -> 2e18.

**THE CAUSE: the damping term never engages.** The update is

    dV = dt (I_set - I_cond) / (C_gap + |g| dt)

and `|g| dt` **never exceeded 12% of `C_gap` at any point in the oscillation** --
typically `1e-20` against `C_gap = 1.77e-16`, i.e. four orders of magnitude
short. So the controller is a PURE INTEGRATOR on the gap capacitance, with no
damping at all, driving a plant whose gain is EXPONENTIAL in V. That is
unstable by construction.

**IT IS NOT A TIMESTEP PROBLEM.** Established with a control at two
resolutions: `iset` (dt ~6e-12) and `En150` (dt ~1e-12) track each other to
**0.66% at worst** through the entire overshoot, peak and turnover --
`Ic/Is` 59/59, 124/123, 191/191, 216/216, 229/228 -- and BOTH diverge. Coarser
arms at dt 1.2e-10 and 2.4e-10 diverge to the same 1e18. Six decades of dt, one
trajectory.

**CANDIDATE 1 (the estimator) IS EXCLUDED, measured.** The reported `g`
matches a direct finite-difference `dI_cond/dV` at EVERY window width -- 1, 10,
100 and 1000 steps -- to 1-2%:

    t=8.50e-7  reported 4.24e-9  w=1 4.24e-9  w=10 4.23e-9  w=100 4.23e-9  w=1000 4.29e-9
    t=9.42e-7  reported 1.70e-7  w=1 1.70e-7  w=10 1.73e-7  w=100 1.71e-7  w=1000 1.54e-7

So `g` IS the true plant gain, and the estimator is sound.

**THE REAL BOUND: the loop response time is `tau_loop = C_gap/g`, and it is
PHYSICS, not a design parameter.** `|g|dt/C_gap` is just `dt/tau_loop`, so the
per-step "damping" reading was a red herring. Regulation requires
`tau_loop < tau_growth`, with `tau_growth = 0.400 ns` measured:

| t | Ic/Is | g [S] | tau_loop | tau_loop/tau_growth |
|---|---|---|---|---|
| 7.99e-7 | 84% | 1.86e-09 | 95 ns | **238x too slow** |
| 9.15e-7 | 227% | 1.16e-09 | 153 ns | **383x too slow** |
| 9.45e-7 | -96% | 6.36e-08 | 2.8 ns | 7x too slow |
| 9.535e-7 | -12045% | 4.15e-05 | 4.3 ps | fast enough -- TOO LATE |

    threshold:  g > C_gap/tau_growth = 4.43e-07 S
    at onset:   g = 1.86e-09 S        -> 238x short

The loop only becomes fast enough once `n_e ~ 1e18`, i.e. after the discharge
is gone.

**AND THIS CORRECTS AN EARLIER CLAIM IN THIS PROJECT.** The statement that
"the gap self-discharges in `tau_gap = C_gap/G = 0.28 ns`, faster than the
0.68 ns growth, so the gap wants to quench itself" used the CHORD conductance
`G = I/V = 6.25e-7 S`. **A feedback loop is governed by the DIFFERENTIAL
conductance `dI/dV`, which at the critical moment is 336x SMALLER.** That is
why the earlier reasoning looked encouraging and the outcome was not. The
0.28 ns figure is superseded for any stability argument.

**THE BOUND IS GENERAL for a two-terminal circuit.** The only way to change the
gap voltage is to move charge on or off the electrode, so ANY circuit whose
sole actuator is the electrode potential has response time `>= C_gap/g`. This
is the same character as the ballast bound
`tau_RC*(I_sc/I_op - 1) = C_gap*V_gap/I_op = 86.6 ns` -- set by geometry and
plasma state, not by component values.

Unlike the earlier "no lumped circuit works" claim, which was WITHDRAWN because
it rested on a crude `nu'` with a 1.95x margin, this rests on a MEASURED `g`
(verified at four window widths) with a **238x** margin.

**dt IS CONCLUSIVELY EXCLUDED.** Four arms spanning six decades of dt -- 1e-12,
6e-12, 1.2e-10, 2.4e-10 -- followed ONE trajectory (agreeing to 0.66% through
the overshoot, peak and turnover) to the same endpoint, n_e ~ 1e18.

**AND A SEPARATE CONTRIBUTOR:** the `setCurrent` ramp is too fast through
ignition. `I_cond` overshot to 229% because the discharge's own growth outran
the ramp, so the loop was never in the quasi-static regime the design assumed.
A slower ramp through ignition reduces the excursion the controller has to
handle, independently of the damping fix.

### STAIRCASE TEST 2026-09-06: the "ramp too fast" diagnosis is FALSIFIED

The `tau_loop = C_gap/g` bound was first attributed to the ramp being 4-7x
faster than the ion transit, with the remedy being a STAIRCASE -- hold each
current level for 2-3 ion transits so `gamma` relaxes to 0. Both a CURRENT
staircase (`grubert2009_iset_stair`, I_set held at 5e-8 A) and a VOLTAGE
staircase (`grubert2009_stair_V`, V_src held at -250 V through R = 1e8, no
shunt C) were run.

**BOTH DIVERGED INSIDE THE FIRST PLATEAU**, with the set point CONSTANT:

    iset_stair  t=9.3e-7  I_set = -5e-8 A held since t=0   n_e -> 3.2e18
    stair_V     t=1.4e-7  V_src = -250 V held since t=0    n_e -> 4.4e18

Plateau 1 spans t = 0..14 us in both, so nothing was ramping. **The remedy
addressed the wrong cause.**

**WHY, and this is the durable result.** `tau_loop/tau_growth` through the
ignition of the held-current arm:

| t | Ic/Is | g [S] | tau_loop/tau_growth |
|---|---|---|---|
| 3.6e-7 | 6% | 9.25e-11 | **4790** |
| 5.2e-7 | 23% | 4.13e-10 | 1070 |
| 7.0e-7 | 80% | 4.05e-09 | 109 |
| 9.3e-7 | -674% | 1.33e-07 | 3.3 |
| 9.3e-7 | -4330% | 9.52e-06 | 0.047 (fast enough -- far too late) |

**At ignition `g` is small BECAUSE THE PLASMA IS JUST FORMING, and that is
exactly when `gamma` is large.** So `tau_loop >> tau_growth` at ignition is
INTRINSIC, not a consequence of how the gap is driven. The loop becomes fast
enough only after the current has overshot by ~4000x.

**CONCLUSION: no two-terminal circuit can carry this discharge through
ignition.** Changing the gap voltage requires moving charge on or off the
electrode, so the response time is bounded below by `C_gap/g` for ANY such
circuit; at ignition `g` is 3-4 decades too small and no choice of R, L, C or
set point alters it. The ballast bound
`tau_RC*(I_sc/I_op - 1) = C_gap*V_gap/I_op = 86.6 ns` and the current-source
bound `C_gap/g` are two faces of the same limit.

Confirmed by the two arms having UNRELATED ACTUATORS -- one imposes current,
one imposes voltage -- and failing identically. That is why the voltage arm was
worth running.

**THEREFORE the steady, current-imposed solver is the REQUIRED instrument, not
a fallback:** it never traverses ignition. This now rests on measurement,
unlike the earlier `nu'`-based "no lumped circuit works" argument, which was
withdrawn for resting on a 1.95x margin.

**What is still untested:** whether the steady operating point, once reached by
any means, is stable under the transient solver. The linear analysis says yes
(`gamma = 0` is stable for every circuit) but no run has ever held it.

### WHAT IS NOT YET VERIFIED

Syntax-checked and reasoned, **but not yet run**: the build guard correctly
refused to relink while the R sweep was running, so this has had no full
`build-all.sh` and no case. Before it is used for a result it needs:

1. a full `build-all.sh` reporting BUILD-COMPLETE;
2. a case where the regulated current is reached and held, checked against
   `I_cond` in `circuit.csv` -- the fixed point is exact, so agreement should
   be to solver tolerance, not to a few percent;
3. a case that HITS the compliance rail, to confirm the clamp behaves as a real
   supply rather than as a failure. A guard whose silence has not been tested
   is not a guard, and this one has two branches that only fire off-nominal.

## Stage 3 — IMPLICIT coupling, which is the real work

**STATUS 2026-09-06: HALF DONE, and the half that is missing is the one that
matters.** The circuit's own RC dynamics are now implicit, which removed a
divergence of 885x per step and is verified to 0.3% against the analytic ramp
response. But the PLASMA CURRENT is still explicit -- it is measured after the
step that produced it -- and post-breakdown that is the coupling that fails:
the discharge drew 8.5x the current the circuit could supply at that instant,
the ballast demanded a gap voltage the source could not provide, and the
cathode was driven to POSITIVE potential. Implicit in the circuit and explicit
in the plasma is not implicit coupling.

Stage 1 chases the operating point; it does not solve for it. The fix is to
expose

    V(I)   and   dV/dI

from the circuit, and let the electrode potential be solved WITH the discharge
rather than lagged behind it — a 1-D Newton iteration on the electrode
potential nested in the outer loop, using dV/dI as the Jacobian. For R, RC,
RLC that derivative is analytic and trivial.

This is the stage that removes `relaxation` and makes a stiff ballast usable.
It is also why the interface must expose the DERIVATIVE, not just the value.

## Stage 4 — ngspice bridge: CONSIDERED AND DECIDED AGAINST (2026-09-06)

Recorded so it is not proposed again. The framework OWNS its lumped models.

The technical reason, which is the one that decided it: the hard part of
circuit coupling is IMPLICITNESS, not solving the circuit. Stage 3 needs
`dV/dI`; for a lumped model that is one analytic line, but through a black box
it has to be probed numerically -- slower, and fragile exactly where the
coupling is stiffest. A SPICE back end would solve the easy half of the problem
and make the hard half worse.

Three further reasons, any one of which would have been enough on its own:
licensing (ngspice is mostly New BSD but not uniformly, which matters for a
commercial product), a runtime dependency added to a solver that today needs
only OpenFOAM, and ngspice managing its own internal timestepping, which would
have to be negotiated against ours.

The interface requirement from stage 3 stands regardless: a circuit exposes
`V(I)` AND `dV/dI`. That is what makes implicit coupling possible, and it is
worth stating as the contract even though there is now only one implementer.

## What must be true of every stage

* The electrode potential is an OUTPUT. Anything that silently reverts it to an
  input is a defect.
* One current, from one place. `I_total` is Sato's, and a circuit that derives
  its own would eventually disagree with the diagnostic.
* An unrecognised `type` is FATAL. A case asking for an RC ballast must never
  silently receive a resistive one.

---

## 2026-09-07: "NO TWO-TERMINAL CIRCUIT SURVIVES IGNITION" IS NOT ESTABLISHED

The user asked whether a lower applied voltage should give a stable glow, and
whether voltage/current continuation from that state had been done PROPERLY.
Checking that reopened the conclusion above.

### The breakdown voltage, computed independently

Uniform-field Townsend criterion on our own `alphaN_vs_reducedE`, d = 1 cm,
gamma = 0.06, needing `alpha*d = ln(1+1/gamma) = 2.872`:

    V_b = 120.8 V     at E/N = 500 Td

| V [V] | alpha*d | multiplication |
|---|---|---|
| 98 | 2.189 | 0.476 (sub-critical) |
| **120.8** | **2.872** | **1.000** |
| 200 | 5.061 | 9.4 |
| **500 (Grubert)** | **9.939** | **1243** |

**Grubert's 500 V is 4.14x breakdown.** Applied to an UNSCREENED gap it gives
1243x multiplication -- the whole-gap avalanche, as arithmetic rather than
mystery. It is only self-consistent once a cathode fall exists to absorb it.

**AND THIS VALIDATES THE MODEL'S BREAKDOWN:** the `iset` runs ignited at
**-113 to -116 V** against this independently computed 120.8 V (3-6%). The
ignition physics is right.

### The continuation was DESIGNED correctly -- in CURRENT, which is the right parameter

`validation/grubert2009_iset_stair`: five plateaus, j = 0.025 -> 0.100 -> 0.250
-> **0.511 (Grubert)** -> 1.000 mA/cm^2, **12 us holds = 2-3 ion transits**,
risers RAMPED over 4 us. Rationale recorded there: "hold each current level long
enough for `gamma` to relax, and each plateau becomes a genuine steady state."
Current, not voltage, is correct: the glow branch has `dV/dI < 0`, so V is not a
function of I there and voltage continuation cannot traverse it.

### BUT IT WAS NEVER EXECUTED. Not one run reached its first plateau.

| case | reached | endTime | fraction |
|---|---|---|---|
| `grubert2009_iset` | 9.536e-07 | 20e-6 | **4.8%** |
| `grubert2009_iset_stair` | 9.340e-07 | 90e-6 | **1.0%** |
| `grubert2009_stair_V` | 1.383e-07 | 90e-6 | **0.2%** |
| `grubert2009_iset_All1500` | 9.531e-07 | 20e-6 | 4.8% |
| `grubert2009_iset_En1500` | 9.538e-07 | 20e-6 | 4.8% |

`iset_stair`'s plateau 1 window was **2-14 us**. It died at **0.93 us** -- before
the first plateau began. **The quasi-static plateau hypothesis was never tested**,
and the conclusion drawn from these runs is therefore about the first ignition
only, not about continuation.

Four of the five stalled within 2% of the SAME physical time (0.93-0.95 us),
i.e. they all died at one event.

### AND IT WAS MOMENTARILY WORKING. The trace, at I_set = 5e-8 A:

| t [s] | I_cond [A] | V_el [V] | g [S] |
|---|---|---|---|
| 2e-12 | 6.2e-15 | -0.0006 | 0 |
| 4.12e-07 | **-4.66e-09** | **-113.1** | 1.49e-10 |
| 9.24e-07 | **+8.25e-08** | -112.5 | 8.92e-07 |
| 9.328e-07 | 5.27e-06 | -154 | 2.80e-05 |
| 9.339e-07 | **7.34e-05** | **-266.1** | 1.03e-04 |

V ramped to -113 V at `I_set/C_gap` and then **SAT at breakdown for ~500 ns**
while the avalanche built -- which IS the marginal `gamma ~ 0` state the design
wanted. `|I_cond|` then reached 8.25e-8 against a 5e-8 setpoint (**1.65x** --
essentially regulated). Ten nanoseconds later it was 7.34e-5 A, 890x over, with
the regulator having driven V from -112.5 to **-266 V**.

### THE SUSPECT -- RETRACTED 2026-09-07, checked and it does not hold

**The regulator is Kirchhoff-correct and there is no sign defect.** Equation (1)
in `plasmaExternalCircuit.C` is `I_set = I_cond + C dV/dt` -- the terminal
current is conduction plus the displacement drawn by the electrode capacitance
-- and (2) `dV = dt (I_set - I_cond)/(C + |g| dt)` is its backward-Euler
discretisation. A sign flip in `I_cond` therefore produces a LARGE error term,
and a large voltage excursion is the CORRECT response to it: the surplus current
charges `C_gap`. Nothing is inverted.

**And the step size is tiny, which kills the "amplifier" reading.** With the
observed imbalance `I_set - I_cond = -1.3245e-07` A and the measured timesteps:

| t [s] | dt [s] | dV per step |
|---|---|---|
| 9.24006e-07 | 4.275e-11 | **-0.032 V** |
| 9.328e-07 | 8.43e-13 | -0.0006 V |
| 9.339e-07 | 1.979e-13 | -0.00015 V |

So the drift from -112.5 V to -266 V accumulated over THOUSANDS of steps at
hundredths of a volt each. There is no per-step blow-up in the controller. My
earlier reading came from an arithmetic slip -- I multiplied `dt*dI/C` wrong by
six orders of magnitude and inferred a 740 V per-step jump that does not exist.

WHAT SURVIVES from this section: the runs all died within 2% of the same
physical time, the `I_cond` sign reversal is real and still unexplained
PHYSICALLY (why does the conduction current reverse?), and the recorded
`tau_loop = C_gap/g` bound stands as the explanation -- at ignition `g` is small,
so changing `I_cond` needs a large `dV`, and charging `C_gap` to it takes
`C/g` = 95 ns against a 0.4 ns growth time. The discharge outruns the regulator,
and that is a property of the plant, not a bug in the controller.

WHAT IS STILL OPEN: the continuation was never EXECUTED (the table above), so the
quasi-static plateau hypothesis remains untested regardless of this retraction.

### The original suspect, kept for the record

At t = 9.24e-07 `I_cond` is **+8.25e-08** while `I_set` is **-5e-08** -- the
right magnitude, the WRONG SIGN. (`I_cond` was correctly negative at
t = 4.12e-07, so it reversed.) The regulator's error term is

    dV = dt*(Iset - Icond)/(C + g*dt)

so with opposite signs the error is `|Iset| + |Icond|` = **2.65x the setpoint**,
pointing to drive `|V|` UP at exactly the moment it should have backed off. The
`|g|` rather than `g` choice is documented as protecting against an inverted
step, but it assumes the plant-gain SIGN is consistent.

**This is at the handover from "ramp V at I_set/C_gap" to "regulate on the error"
-- which is precisely where every one of these runs died.**

NOT ESTABLISHED, and the distinction matters: I cannot yet separate (i) a
regulator that mishandles a sign reversal in `I_cond` from (ii) a genuine
physical loss of control. The first is a bug with a one-line test; the second is
the recorded `C_gap/g` bound. **The `I_cond` SIGN CONVENTION must be checked
against its definition first** -- see [[verify-before-claiming]] -- because if
the reversal is a diagnostic artefact the whole reading changes.

### What this means for the standing conclusion

`tau_loop = C_gap/g` being 238x too slow at ignition is still a real measurement.
But it was measured on runs that lost control at the regulator handover, so it
may be describing the consequence rather than the cause. **The staircase design
remains untested and is the cheapest way to find out**: if a plateau at
j = 0.025 mA/cm^2 can be held, continuation is viable and the standing
conclusion is wrong.
