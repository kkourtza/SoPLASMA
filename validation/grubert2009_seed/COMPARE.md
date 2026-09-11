# Is Grubert's solution a STABLE FIXED POINT of our model?

```compare
question:  Started AT Grubert's profile with the current imposed at I_op, does the transient solver HOLD it?
baseline:  /home/kkourtza/soplasma-scratch/validation/grubert2009_iset
varies:    initial condition, setCurrent
matches:   circuit type, compliance, capacitance, maxSpeciesCo, driftDiffusionFluxScheme, electronEnergyModel, mesh, chemistry, gammaSEE, electronReflection
time:      2e-5
field:     n_e
region:    domain
```

## What this tests, and what it deliberately does NOT

**It tests the DESTINATION, not the path.** Measured 2026-09-06: the transient
solver cannot carry this discharge THROUGH ignition, because at ignition the
plasma conductance `g` is small (the plasma is forming) exactly when the growth
rate `gamma` is large, so the loop response time `C_gap/g` is 3-4 decades too
slow -- for ANY two-terminal circuit. Confirmed by a current-imposed and a
voltage-imposed staircase failing identically with their set points HELD
CONSTANT.

That is a statement about the approach. This case asks the separate and more
important question: **does our model have Grubert's solution as a stable fixed
point?** If yes, the closure, the coefficients, the wall fluxes and the circuit
are all validated at the operating point, and only the approach remains open.

**It is NOT a production route.** A solver that needs the answer as input is
not a solver, and that was rejected (user, 2026-09-06) for exactly the right
reason. Seeding is legitimate as a TEST and illegitimate as a method.

**And no new solver is needed for it:** transient relaxation from a good
initial state IS pseudo-timestepping, which is the only steady mechanism
OpenFOAM offers anyway.

## The seed, and its verified consistency

`n_e` and `n_Ar+` interpolated onto the mesh from
`SoEEDF/validation/dias2025/reference/grubert_fig3_na_LMEA_{electrons,Arp}.csv`
(1e9 cm^-3 = 1e15 m^-3; `z_over_d` from the cathode, and this mesh has the
cathode at x = 0). Clamped at the declared floors -- the digitisation dips
slightly negative at z/d = 0.

    n_e    1e11 .. 2.477e15 m^-3, peak at x/L = 0.398   (Grubert: 2.478e15 @ 0.394)
    n_Ar+  3.2e14 .. 6.391e15 m^-3, peak at x/L = 0.205 (Grubert: 6.391e15 @ 0.203)

**THE CONSISTENCY CHECK THAT MATTERS.** The two curves were digitised
INDEPENDENTLY, so their DIFFERENCE -- which is the space charge, and therefore
the field -- carries the error of both. Integrating Poisson through the seeded
profile with the physical boundary condition `E = 0` at the anode:

| quantity | seeded | expected | ratio |
|---|---|---|---|
| integrated gap voltage | **-601.8 V** | -500 V | **1.20** |
| \|E\| at the cathode | 3.31e5 V/m | ~1.8e5 (2V_c/d_c) | 1.84 |
| cathode fall thickness | 2.95 mm | ~4.4 mm | 0.67 |

**Within 20% on the gap voltage from two independently digitised curves.** That
is a strong independent check on the digitisation, and it means the seed is a
sensible state rather than a guess. `nEps_e` is not seeded -- the solver
derives it from the LFA equilibrium at the solved field.

(An earlier verdict that the seed was "not self-consistent" compared the field
EXCURSION to the MEAN field, which is meaningless in a glow because the field
is deliberately concentrated in the fall. Withdrawn.)

## Structure of the seed -- it is a textbook glow

    x/L 0.000-0.10 : n_Ar+ ~5e15, n_e at the floor   POSITIVE SHEATH (cathode fall)
    x/L 0.74       : n_Ar+ 1.15e15, n_e 1.10e15      QUASINEUTRAL (negative glow)
    x/L 0.96-1.00  : ratio 3-7                       anode region

Note this is the OPPOSITE of every transient run to date, which put both peaks
at z/d ~ 0.97 -- a Townsend avalanche, not a glow.

## Success and failure, in advance

* **SUCCESS:** the profile holds. `n_e` stays ~2.5e15, `n_Ar+` ~6.4e15, the
  peaks stay at x/L ~ 0.4 and ~0.2, `V_gap` relaxes toward -500 V, and
  `I_cond -> I_set`. Grubert's solution is then a stable fixed point of our
  model and the physics is validated at the operating point.
* **it drifts to a NEARBY state:** the fixed point exists but is displaced --
  a real, LOCALISED discrepancy, and finally a diagnosable one, because it is
  measured against a stable state instead of a transient.
* **it collapses to the floor:** the operating point is not sustainable in our
  model at this current, which would contradict the stability analysis and
  point at the wall fluxes or gamma.
* **it runs away as before:** the fixed point is unstable in the full nonlinear
  system despite the linear analysis, which would be the strongest possible
  argument for a genuinely steady (Newton) solver.

**Judge only after a few `tau_diel` of settling** -- the seed is 20% off in
voltage and its space charge is the difference of two noisy curves, so an
initial transient is expected and is not a failure.

## Extraction

    postProcess -func writeCellCentres -time <t> && ./extract_centreline.py <t>
    ~/ct-env/bin/python ../../../Projects/SoEEDF/validation/dias2025/plot_grubert_profiles.py \
        centreline_<t>.csv --arm LMEA

---

# RESULT 2026-09-06: THE SEED IS COMPROMISED -- differencing two digitised curves

**The test could not be run as designed, and the reason is worth more than the
test would have been.**

The space charge is `e(n_Ar+ - n_e)`. In the negative glow those two curves lie
within a few percent of each other, so **a few percent of digitisation error in
each becomes ~100% error in their difference** -- and the difference is what
sets the field. Measured from the seeded state with the case's own mobility
table:

| x/L | n_e | E/N [Td] | j_e = e n_e mu_e E [A/m^2] |
|---|---|---|---|
| 0.035 | 1e11 | 12443 | 0.078 |
| **0.497** | **2.18e15** | **380** | **93.1** |
| 0.832 | 6.82e14 | 284 | 22.6 |

Grubert's TOTAL is **5.11 A/m^2**, so the seeded state carries **18x** that in
the negative glow. For the current to be right there the bulk field would have
to be ~500 V/m (21 Td) rather than 9170 V/m -- i.e. **the negative glow should
be nearly field-free and the seeded field is not.**

**The -601.8 V consistency check MASKED this**, because it is an INTEGRAL over
the whole gap: it averages the distribution error away. Checking an integral
and concluding the profile is right was the mistake -- the quantity that
mattered was the LOCAL field in the bulk.

## What this run still established

* the resume path works: `RESUMED from the field, V_electrode = -601.827239 V`
* the solver HELD the seeded ion profile to **0.1%** (6.401e15 vs 6.391e15) and
  the field to the seeded value
* `dV/dt` matched `(I_set - I_cond)/C_gap` to **0.2%** -- a third independent
  confirmation of `currentSource` after 1.33% (pre-ignition) and 4% (ignition
  voltage)
* **and it explains the tiny timestep.** The primary limiter was
  `temporal error (PI)` on 1319 of ~1865 steps, with 558 outer-loop failures;
  dt collapsed from 4e-12 to 2e-13. For scale, `iset` ran at dt up to 4.5e-10
  with ZERO discards on the same mesh and chemistry. **The 2000x smaller dt is
  a symptom of the inconsistent initial state, not a solver limitation** -- the
  accuracy controller is correctly refusing to take large steps through a
  violent relaxation.

## The fix, for the next attempt

Seed `n_Ar+` as digitised, and **construct** `n_e` rather than digitising it:

* inside the cathode fall, where the two curves differ by ORDERS of magnitude,
  differencing is reliable -- keep the digitised electron profile.
* in the quasineutral bulk, **set `n_e = n_Ar+`** so `rho ~ 0` by construction,
  which is what a negative glow requires, instead of inheriting a spurious
  space charge from two independently digitised curves.

The transition point is where the digitised ratio `n_Ar+/n_e` falls below a few
-- measured at x/L ~ 0.25 in this data.

**And check the LOCAL field, not the integrated voltage, before running.**


---

# THE NO-RUN RESULT, and why the consistent-seed construction failed

## THE RESULT WORTH KEEPING: 394 V against Grubert's 500 V

Using Grubert's digitised densities, OUR transport tables, and CURRENT
CONTINUITY -- exact in steady state -- with no seeding and no run at all:

    j_tot = e(n_e mu_e + n_i mu_i) E = 5.11 A/m^2 everywhere
      =>  E(x) = 5.11/(e(n_e mu_e(E) + n_i mu_i(E)))     [iterate, mu depends on E/N]

| x/L | n_e | n_i | E [V/m] | E/N [Td] |
|---|---|---|---|---|
| 0.006 | 1e9 | 4.73e15 | **3.07e5** | 12473 |
| 0.151 | 1.87e12 | 6.12e15 | 1.43e5 | 5903 |
| 0.401 | 2.48e15 | 2.69e15 | **297** | **12.3** |
| 0.741 | 1.11e15 | 1.15e15 | 718 | 29.7 |
| 0.991 | 6.88e13 | 3.93e14 | 1.71e4 | 708 |

    integrated gap voltage = 393.9 V   vs Grubert 500 V   -> 0.788 (21% low)

**The field-free negative glow emerges by itself** (12-30 Td), which is what a
glow requires and what the Poisson-on-a-difference seed could not produce.

**WHY THIS WORKS AND THE OTHER DID NOT -- the general lesson.** Current
continuity uses `n_e mu_e + n_i mu_i`, a **SUM**. Poisson uses `n_i - n_e`, a
**DIFFERENCE**. Where the two densities agree to a few percent, a few percent of
digitisation error stays a few percent in the sum and becomes ~100% in the
difference. **Use the constraint that sums.**

CAVEAT: drift-only. In the negative glow E is small (12 Td) and the density
gradient steep, so diffusion is comparable to drift there and that region's field
is uncertain. The cathode fall, which dominates the voltage integral, is
drift-dominated and unaffected.

## THE CONSISTENT-SEED CONSTRUCTION FAILED, and the reason is structural

Attempted: impose current continuity AND Poisson simultaneously, solving

    n_e = (P - D mu_i)/(mu_e + mu_i),  n_i = n_e + D,
    P = J/(e E),   D = (eps0/e) dE/dx

by iteration, so the SUM stays pinned to the digitised data while the DIFFERENCE
is derived from the field rather than digitised. Measured outcome:

    integrated gap voltage collapsed to 13.8 V (from 394 V)
    |E| peak moved to x/L = 0.999 -- the ANODE
    the cathode fall was destroyed; n_i at x/L=0.151 fell 6.1e15 -> 1.1e15

**Two of the three steady constraints is UNDERDETERMINED.** The third is the
PARTICLE BALANCE (ionisation = wall loss + recombination). Without it nothing
sustains the cathode fall -- the fall exists because ion flux from the bulk plus
secondary emission maintain it -- so the iteration lowered the field and drifted
to a fall-free state.

**And imposing all three IS the steady solve.** The seeding route is therefore
circular: constructing a consistent state requires the solver it was meant to
avoid. That is a real argument for building it, and it is the third independent
one today, after the `C_gap/g` bound and the staircase falsification.
