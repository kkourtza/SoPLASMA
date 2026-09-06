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

**TWO CANDIDATE FIXES, not yet chosen between:**

1. **`g` is estimated too small.** It is a secant `dI/dV` from successive
   accepted steps; at dt ~1e-12 the per-step `dV` is tiny, so the estimate may
   be dominated by noise or suppressed by the relative-change guard on it, and
   never reflect the true plant gain. If so, fix the estimator.
2. **The controller FORM is wrong for this plant.** A pure integrator cannot
   stabilise an exponentially-nonlinear plant however well `g` is estimated; it
   needs proportional action, a slew limit on `dV/dt`, or both.

These are distinguishable: instrument `g` against a directly measured
`dI_cond/dV` over a finite perturbation. Do that before changing the form.

**AND A SEPARATE CONTRIBUTOR:** the `setCurrent` ramp is too fast through
ignition. `I_cond` overshot to 229% because the discharge's own growth outran
the ramp, so the loop was never in the quasi-static regime the design assumed.
A slower ramp through ignition reduces the excursion the controller has to
handle, independently of the damping fix.

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
