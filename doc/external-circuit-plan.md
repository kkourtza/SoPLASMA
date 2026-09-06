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

## Stage 3 — IMPLICIT coupling, which is the real work

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
