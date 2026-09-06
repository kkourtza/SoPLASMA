# R sweep: does a larger ballast get through the breakdown transient?

## The question

The ballasted case reaches Grubert's operating point in passing and then
diverges past it -- n_e climbing to 3.8e18 / 2.3e19 m^-3 against a reference
2.478e15, an ionisation degree of 1.6e-4 / 9.5e-4 against 1.03e-7 for a glow.

The physical proposal being tested (user's, 2026-09-06): this is a TRANSIENT
overshoot, not a wrong fixed point. The load line already passes through the
right operating point; what fails is the passage to it. A larger R holds the
gap voltage down harder while the current rises, so the discharge may reach the
fixed point without overshooting first.

## Why larger R pins the current -- and why this contradicts my earlier change

The load line is `V_gap = V_src - R*I`. Fix it through Grubert's point
(I_op = 1.022e-6 A at V_gap = -500 V) by setting `V_src = -500 - R*I_op`. The
ballast's short-circuit current, as a multiple of the operating point, is then

    I_sc/I_op = V_src/(R*I_op) = 1 + 500/(R*I_op)

    R = 1e8   ->  5.89x
    R = 1e9   ->  1.49x
    R = 5e9   ->  1.10x

**Larger R is a TIGHTER current pin, not a looser one.** In the limit
R -> infinity with V_src = 500 + R*I_op, the ballast IS a current source at
I_op. So this sweep and the planned `currentSource` (stage 2a of
`doc/external-circuit-plan.md`) approach the same thing from two directions,
and the sweep tests the physics before the new model is written.

**This cuts against my own earlier change**, 5e8 -> 1e8 (recorded in
`configuration/boundaries`), which lowered R to give the ignition transient
MORE headroom. Both cannot be right, and the difference is what the overshoot
IS:

- if the transient current demand is genuine physics, a source that cannot
  supply it drives `V_src - R*I` positive -- a positive cathode -- and starves
  the discharge. That was the measured 5e8 failure, and headroom is the fix.
- if the demand is a RUNAWAY, headroom feeds it. n_e four decades above both
  the reference and the independent screening-arrest estimate (1.6e16) says
  runaway, and a pin is the fix.

The 5e8 failure also had a confound that is now removed: the 5 us source ramp
had only reached -257 V at breakdown, so the source was starved by the RAMP as
well as by R. Both arms here take full source voltage in 50 ns.

## Arms, and what makes this a valid control

Baselines, absolute paths:

- `/home/kkourtza/soplasma-scratch/validation/grubert2009_ballast`
  R = 1e8, C = 5e-14, V_src = -602 V.  tau = R*C = 5.0 us.
- `/home/kkourtza/soplasma-scratch/validation/grubert2009_R1e9`
  R = 1e9, C = 5e-14, V_src = -1522 V. tau = 50 us.
- `/home/kkourtza/soplasma-scratch/validation/grubert2009_R5e9`
  R = 5e9, C = 5e-14, V_src = -5610 V. tau = 250 us.

**What differs:** R, by 50x across the sweep.
**What must match, and does:** **C is held at 5e-14 F in every arm**, so R is
the only circuit variable. `tau = R*C` therefore does NOT match across arms
(5 us / 50 us / 250 us) -- and holding C rather than tau is a deliberate
correction to an earlier version of this file.

CORRECTED 2026-09-06, before the arms passed breakdown. The first version held
tau = 5 us fixed and claimed that made the voltage rise rate equal. It does
not. The rise rate at breakdown is `dV/dt ~ V_src/tau`, and V_src is pinned by
the load line to `500 + R*I_op`, so at fixed tau it scales with R:

    tau fixed at 5 us:   dV/dt = 1.2e8 -> 3.0e8 -> 1.1e9 V/s   (UP 9.3x)
    C   fixed at 5e-14:  dV/dt = 1.2e8 -> 3.0e7 -> 2.2e7 V/s   (DOWN 5.4x)

The first version therefore made the approach to breakdown 9x HARSHER as R
grew, confounding the thing being tested against the one intervention already
known to matter -- a gentle rise was the whole reason breakdown became well
behaved. At fixed C the rise rate instead falls and converges to
`I_op/C = 2.0e7 V/s`, its current-source limit, so larger R is a tighter
current pin AND a gentler rise: the two effects agree instead of cancelling.

Fixed C is also the physically honest sweep. C is stray capacitance, a property
of the rig; R is the ballast, the knob an experimenter actually turns. Nobody
retunes their coax to keep R*C constant.

Consequence to expect: breakdown arrives LATER at larger R (1.8, 6.3, 8.1 us
from the RC rise), and endTime 45 us covers all three. V_src is not free
either -- it is fixed by the load line through the same operating point. Mesh,
chemistry, `electronEnergyModel LMEA`, Courant settings, and now
`driftDiffusionFluxScheme standard` are identical.

**Caveat on the baseline arm:** `grubert2009_ballast` was run with
`fluxScheme ScharfetterGummel`, which was found on 2026-09-06 to impose the
wrong wall flux and is now fatal. At the refined electrode cell's Pe ~ 0.1-0.2
the two branches agree to ~1%, so its numbers stand, but a re-run on
`standard` is required before it is quoted as a like-for-like control.

## Reference numbers (Grubert et al. 2009, argon, 100 Pa, 1 cm, -500 V)

| quantity | reference |
|---|---|
| j at the cathode | 0.511 mA/cm^2 |
| n_e peak | 2.478e15 m^-3 |
| n_Ar+ peak | 6.39e15 m^-3 |
| ionisation degree | 1.03e-7 |

Independent check, not from the paper: space-charge screening arrests growth
when tau_diel = tau_ionisation, giving n_e ~ 1.6e16 m^-3 -- within a factor 6.4
of the reference, and four decades below what the diverged runs reach.

## Extraction

    cd <case> && ~/ct-env/bin/python status.py

Reports steps, sim time, rejections, n_e range, and |j| in mA/cm^2 with the
reference alongside. The discriminating observable is **whether n_e settles
near 1e15-1e16 or continues past 1e17**; per rule 17 that is visible within a
few us of breakdown, so an arm can be judged long before endTime.

## The machine-readable contract

```compare
question:  Does a tighter current pin get the discharge through breakdown without overshooting into runaway?
baseline:  /home/kkourtza/soplasma-scratch/validation/grubert2009_R1e9
baseline:  /home/kkourtza/soplasma-scratch/validation/grubert2009_R5e9
varies:    resistance
matches:   capacitance, driftDiffusionFluxScheme, electronEnergyModel, mesh, chemistry, Courant settings, endTime
time:      4.5e-5
field:     n_e
region:    domain
```

`varies: resistance` only -- `sourceVoltage` also differs, but it is NOT a
second free variable: the load line pins it to `V_src = -500 - R*I_op` so that
every arm has the SAME correct operating point available. This arm is
R = 1e8 Ohm, V_src = -602 V, short-circuit current 5.89x the operating point.

**`time: 4.5e-5` is endTime, deliberately.** The question is where the
discharge SETTLES, so the contract compares the settled state. An arm that
never reaches endTime has diverged -- and the tool erroring on a missing time
is then the answer, not an obstacle. Use `--time` to look earlier.

**`region: domain`, and the band regions MUST NOT be used on this case.**
`compare_cases.py`'s `cathode-band` / `anode-band` / `bulk` slice along **y**
(`b[0]`, `b[1]` are the y bounds). This case's gap runs along **x** -- see
`gap1cm.geo`, "cathode at x = 0, anode at x = L" -- and y is the 200 um
transverse direction closed by symmetry planes. A band region here would
measure a slab of the transverse direction and mean nothing. Recorded
2026-09-06; the tool needs an axis option before those regions are usable on a
plane-parallel case.

**No `reference:` line yet, and that is not an oversight.** The tool re-measures
`reference` entries against a baseline CASE, so it cannot hold a literature
value. Grubert's numbers are in the table above and are the physical target; a
`reference:` line gets added once a baseline arm has a settled time directory
to reproduce.
