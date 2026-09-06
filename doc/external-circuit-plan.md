# External circuit coupling — staged plan

Written 2026-09-06, after the Grubert dc glow showed why the capability is
needed. Stage 1 is IMPLEMENTED; stages 2-4 are design.

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

## Stage 2a — `currentSource`, and why it comes FIRST

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
* **Fast transients (pulsed, RF).** A current source is the wrong instrument:
  there the circuit's own dynamics ARE the physics. Use seriesRC/RLC.

### Binding it in the generator, so the user writes one thing

Same route as the ballast, on the electrode it feeds:

    cathode
    {
        kind        ballastedElectrode;      // role name to be revisited:
                                             // "drivenByCircuit" is truer once
                                             // the source need not be a ballast
        circuit
        {
            type        currentSource;
            current     -1.022e-6;           // [A], or a Function1 to sweep
            compliance  -1500;               // [V] supply limit, REQUIRED
        }
    }

and nothing else: the electrode is the patch it is written on, `plasmaCreate-
SpeciesFields` emits `fixedValue` for it as now, and the solver reads it
through `boundaryRoleLibrary::caseDeclaration`. Current in AMPERES, not a
density -- the electrode area is a property of the mesh and deriving j from it
is the solver's job, not the user's.

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
