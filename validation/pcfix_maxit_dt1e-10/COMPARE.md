# Can the solver reach Grubert's operating point FROM SCRATCH, under current control?

```compare
question:  Under current control, does the discharge reach Grubert's operating point from scratch -- no seeding -- instead of running away?
baseline:  /home/kkourtza/soplasma-scratch/validation/grubert2009_fix_n11
varies:    circuit
matches:   minNumberDensity, driftDiffusionFluxScheme, electronEnergyModel, mesh, chemistry, gammaSEE, electronReflection
time:      2e-5
field:     n_e
region:    domain
```

## Why this case exists

The user rejected seeding the transient from Grubert's digitised profiles --
"the solver should be able to solve this from scratch" -- which is right: a
solver that needs the answer as input is not a solver, and a SaaS user cannot
supply a seed. This is the from-scratch alternative.

## The argument, from `tools/glow_stability.py` (2026-09-06)

The instability criterion is `c0 = R(I*nu' - gamma*G) - gamma`, where
`gamma = nu_iz - nu_loss` is the distance from equilibrium. Two facts follow:

* **at `gamma = 0` every circuit is STABLE**, including the R = 1e8, C = 5e-14
  ballast that ran away. The operating point is not the problem.
* **at `gamma = 1.47e9` (the measured runaway) no series R, RL or RLC helps**,
  because the blocking condition `nu' > gamma/V` contains no circuit parameter.
  (Caveat: that failed by only 1.95x on a crude `nu'`, so it is
  "none found", not "none exists".)

So the fix is to never let `gamma` leave zero. **A voltage ramp cannot promise
that** -- V keeps rising past the self-sustaining value, `alpha*d` reaches
2-4x equilibrium, and `gamma` runs away. **A current ramp can**, by
construction: the discharge cannot draw more than `I_set`, so it is pinned to
its own characteristic and tracks it quasi-statically.

## What differs from the baseline

`varies: circuit` only. `grubert2009_fix_n11` is
`ballastedElectrode` + `seriesRC` (R = 1e8, C = 5e-14, V_src = -602); this is
`currentDrivenElectrode` + `currentSource`. Floors, flux scheme, energy model,
mesh, chemistry, `gammaSEE 0.06` and the anode `electronReflection 0.36` are
identical.

    setCurrent  table ((0 1e-8) (10e-6 1.022e-6))   ramp to Grubert's I_op
    compliance  -1500                                signed rail, 3x headroom
    capacitance 0                                    LOAD-BEARING, see below

**`capacitance 0` is the load-bearing setting.** With `C_gap` alone
(1.77e-16 F) an excess conduction current empties the gap in **0.286 ns**,
faster than the 0.68 ns ionisation growth, so the source can quench a fast
excursion. The 5e-14 F shunt the ballasted arms carried makes that **80.9 ns --
119x too slow**, which is why they could not hold it. The cost is a fast
pre-ignition ramp (`dV/dt = I_set/C_gap`), and that is harmless here precisely
because the CURRENT is pinned.

## Reference (Grubert et al. 2009, argon, 100 Pa, 1 cm, -500 V)

| quantity | reference |
|---|---|
| j at the cathode | 0.511 mA/cm^2 |
| gap voltage | -500 V |
| n_e peak | 2.478e15 m^-3 |
| n_Ar+ peak | 6.39e15 m^-3 |
| ionisation degree | 1.03e-7 |
| cathode fall thickness | ~4.4 mm |

**The discriminating observable is the GAP VOLTAGE**, and that is new. The
current is imposed, so matching `j` proves nothing -- it is an input. What the
run must produce is `V_gap -> -500 V` at `I_set = I_op`. That is a genuine
prediction and it is the one number Grubert's operating point supplies that
this case does not.

Then the structural discriminators, per
`../grubert2009_fix_n11/COMPARE.md`: density peaks moving to the CATHODE side,
`|E|` localised into a ~4.4 mm fall with a field-free bulk, positive space
charge at the cathode.

## What failure looks like, stated in advance

* **`V_gap` pinned at the -1500 V rail** -> the discharge cannot carry `I_set`
  at any voltage; either `I_set` is too high or the model under-produces.
  Real supply behaviour, not a solver failure (see the model's own banner).
* **`V_gap` collapsing toward 0** -> the discharge carries `I_set` far too
  easily, i.e. over-productive.
* **still runs away** -> the current pin is not sufficient, and the honest
  conclusion is that a transient route needs the `nu'` measurement and a
  proper circuit search, or the steady solver.
* **relaxation oscillation** -> bounded but unsteady. Would mean the source
  quenches and re-ignites; interesting, and NOT the same as a runaway.

## First use of `currentSource`

This model was written on 2026-09-06 and has **never been run**. Its two guards
(`resistance` on a current source; negative `setCurrent`) are exercised, but
the update path is not. Treat the first failure as a defect in the model before
concluding anything about the physics.

## Extraction

    grep -c "Charge density updated" logs/log.soPlasmaFoam   # vs step count, >= 1
    tail postProcessing/externalCircuit/circuit.csv          # time,I_set,I_cond,V,g
    postProcess -func writeCellCentres -time <t> && ./extract_centreline.py <t>
