# Reaching the Grubert operating point: two routes, and the analysis that decides

Status: DESIGN for review, 2026-09-06. Nothing implemented.
Written in response to two objections from the user, both of which change the
design.

## 0. Part A monitoring: VERIFIED (record, since it changed no code)

Rule 23 needs both directions. Three stages, each with its own rebuild:

| stage | change | expected | result |
|---|---|---|---|
| A | reintroduce the `finalIter` regression | invariant fires | **FIRED**: `chargeDensity = 0` vs `sum q n = 1.602e-14`, cause #1 listed = the actual cause |
| B | keep regression, disable invariant | audit fires | **FIRED** at step 1: `PER-STEP UPDATE DID NOT RUN: updateSurfaceCharge` |
| C | revert both | silence | exit 0, no fire, `3.00 updates/step`, closing report correct |

Stage B is worth noting: it caught `updateSurfaceCharge`, whose staleness has
NO other diagnostic. The audit is not redundant with the invariant.

## 1. The user's objection to a steady solver, and why it is right

> "openfoam's steady state solvers use pseudo-timestepping so not sure if the
> transient can indeed be avoided"

Correct, and it kills the naive version of the plan. A SIMPLE-type or
pseudo-transient steady solver still marches, and a fixed point that is
unstable to time-marching is unstable to pseudo-marching too. Relaxation does
not rescue an unstable mode; it only rescales it.

**But the instability is a fixed-VOLTAGE instability, not a property of the
fixed point.** On the normal-glow branch `dV/dj < 0`, and a negative
differential resistance is unstable under voltage control and STABLE under
current control. That is why every laboratory glow has a ballast, and it is the
lever here:

* **impose the CURRENT, solve for the voltage.** The gap voltage becomes a
  single scalar unknown with the constraint `INT j dA = I_set`. The fixed point
  is then stable, and pseudo-transient marching converges to it.
* the machinery already exists. `floatingElectrode` solves exactly this shape
  of problem -- potential UNKNOWN, an integral quantity KNOWN -- in CLOSED FORM
  by superposition, using the unit-potential field `psi` and `C_self`:
  `V_f = (Q - Q_rho)/C_self`. A current-imposed electrode is the same
  construction with current in place of charge. This is not new
  infrastructure, it is a second use of validated infrastructure.
* Newton on the steady residual `F(x) = 0` is indifferent to time stability in
  any case, and is the fallback if the current-imposed march still struggles.

So the design is **steady + current-imposed**, not steady + voltage-imposed.
Grubert solved at fixed voltage time-independently, which works for Newton;
current-imposed works for both and is more robust.

## 2. The user's circuit suggestion, and the analysis that settles it

> "if the circuit responded faster to the high and fast growth of discharge
> current ... wouldn't it reduce the E field and thus the growth and thus
> regulate the glow? ... what about a more complex circuit?"

**The physics permits it, and this is the most important thing measured today
about the circuit.** The gap self-discharges through its own plasma
conductance:

    tau_gap = C_gap/G_plasma = 1.77e-16/6.25e-7 = 0.28 ns

which is FASTER than the 0.68 ns ionisation growth. **The gap wants to quench
itself.** What defeats that is the ballast holding the voltage up, with
`tau_ext = R*C_gap` = 17.7 ns at R = 1e8. So the runaway is not inherent -- it
is the external circuit preventing a self-quench that would otherwise happen.

And a second clue points the same way: **the classical DC ballast criterion is
already satisfied.** Stability to slow perturbations needs `R > |dV/dI|`, and

    |dV/dI| = 1/|g| = 1/6.25e-7 = 1.7e6 Ohm     vs     R = 1e8 Ohm

a 59x margin. The operating point IS stable to slow perturbations. What we
observe is a FAST mode that this criterion does not cover, living on `C_gap`.

### What must NOT be done: pick L from an envelope

Two plausible criteria for a series inductance disagree by 60x:

| criterion | L | comment |
|---|---|---|
| develop the full gap voltage at the runaway `dI/dt` | 1.09 mH | `L*dI/dt = 200 V` at `dI/dt = 1.84e5 A/s` |
| dominate the impedance at the instability frequency | 68 mH | `omega*L > R`, `omega = 1/tau_g = 1.47e9` |

and at 68 mH the circuit's own LC period is 21.7 ns, still slower than 0.68 ns.
So neither settles it, and shipping either number would be exactly the
"remembered value" G2 forbids. **This needs a small-signal stability analysis,
not arithmetic.**

### THE ANALYSIS (task 2a -- cheap, no run required)

Linearise the discharge about the operating point as a one-pole element with
differential conductance `g < 0` and response time `tau_g`:

    i(s) = g_0/(1 + s*tau_g) * v(s)

in series with the circuit impedance `Z(s)`, and in parallel with `C_gap`. The
loop is stable iff `1 + Z(s)*Y_gap(s)` has no right-half-plane zeros -- a
Nyquist criterion. Evaluate for:

    Z = R                    (seriesResistor)
    Z = R, C across the gap  (seriesRC -- what we ran)
    Z = R + sL               (seriesRL)
    Z = R + sL + 1/(sC)      (seriesRLC)

and produce the STABILITY BOUNDARY in `(R, L)` for the measured `g_0`,
`tau_g = 0.68 ns` and `C_gap = 1.77e-16 F`. That map is the deliverable: it
either exhibits a realisable `(R, L)` that is stable AND has
`I_sc/I_op ~ 1`, or it proves no lumped two-terminal circuit can do both --
and either answer is worth having.

**Why this is the right instrument:** it is a unit analysis with ground truth
(the Nyquist criterion is exact for the linearised system), it costs minutes,
and it tells us WHICH circuit to build before we build one. Per rule 17 that
beats a CFD sweep, which could only show that two arms disagree.

Validate the analysis against what we already measured: it must predict
INSTABILITY for `R = 1e8, C = 5e-14` (fix_n11, observed unstable) and for
`R = 1e6, C = 5e-14` (fix_fast, observed unstable). An analysis that does not
reproduce those two is wrong and must not be used to choose a circuit.

### Then (task 2b) implement whatever 2a selects

`plasmaExternalCircuit` already has the topology hook. Adding `seriesRL` /
`seriesRLC` is one more branch and one more implicit update; the inductive term
is `V = V_src - L dI/dt - R I`, which discretised backward-Euler with the same
secant `g` gives

    x = [V_src - R I^n - V^n - (L/dt)(I^n - I_tot^n)] / (1 + a(L/dt + R)),
    a = g + C/dt,   x = V^{n+1} - V^n

so it costs no new machinery. **It also preserves the transient solver, which
the user wants** -- and it is a real capability regardless of Grubert, since
inductive ballasts are common in dc and pulsed rigs.

## 2c. RESULT of the stability analysis, 2026-09-06 (`tools/glow_stability.py`)

Model: the discharge as `Y_d(s) = G + I*nu'/(s - gamma)` -- conductance plus an
ionisation pole -- in parallel with `C_gap`, in series with `Z(s)`. Clearing the
pole gives a polynomial; Routh-Hurwitz is the criterion. `gamma = nu_iz -
nu_loss` is the LINEARISATION POINT: zero at the operating point by definition,
`1/tau_g = 1.47e9` during the runaway.

**LIVENESS PASSED** -- the model predicts UNSTABLE for both configurations
already observed unstable (`fix_n11` R=1e8 C=5e-14; `fix_fast` R=1e6 C=5e-14),
so it is usable.

**THE OPERATING POINT IS STABLE -- WITH THE CIRCUIT WE ALREADY HAVE.**

| circuit | at gamma = 0 | at gamma = 1.47e9 |
|---|---|---|
| R=1e8, C=5e-14 (as run) | **STABLE** | unstable |
| R=1e8, C=0 | **STABLE** | unstable |
| R=1e6, C=0 | **STABLE** | unstable |
| R=1e8, L=1 mH, C=0 | **STABLE** | unstable |

The blocking coefficient is `c0 = R(I*nu' - gamma*G) - gamma`. **At gamma = 0
this is `R*I*nu' > 0` ALWAYS**, for any circuit. So the instability is a
property of being FAR FROM THE OPERATING POINT, not of the ballast. No circuit
change is needed for the target state, and none rescues the path to it.

### WITHDRAWN: "no lumped circuit can stabilise it"

At `gamma = 1.47e9` the criterion needs `nu' > gamma/V`, which contains no
circuit parameter -- which is why the R x L search found nothing. But that
conclusion is **NOT ROBUST**: it fails by only 1.95x, and `nu'` here is a
two-point estimate that assumed a UNIFORM field (`V = E*d`), while a cathode
fall is anything but. A 1.95x error in `nu'` flips the answer, and that is well
inside the uncertainty. So the honest statement is: *no circuit was found, and
the search is inconclusive because `nu'` is not known well enough.* Tightening
`nu'` would require differentiating `nu_iz` along the ACTUAL profile.

## 2d. THE CONSEQUENCE: seed the transient AT the operating point

Since the target state is stable with the existing circuit, the fix is not a
new circuit or a new solver -- it is **not traversing the unstable path**.
Start the transient near the operating point so `gamma ~ 0` from the first step.

**And the ideal seed already exists: Grubert's own digitised figure 3.**
`SoEEDF/validation/dias2025/reference/grubert_fig3_na_LMEA_{electrons,Arp}.csv`
carry `n_e(z/d)` and `n_Ar+(z/d)` for exactly this case, in 1e9 cm^-3 = 1e15
m^-3. Interpolate them onto the mesh, let Poisson solve the field from the
resulting charge density, and run.

This is a STRONGER test than reaching the state from scratch, and cheaper:

* if the solution STAYS -- Grubert's profile is a fixed point of our model.
  That is a direct validation of the closure, the coefficients and the wall
  fluxes, and it is the claim we actually want to make.
* if it DRIFTS -- we learn WHERE and HOW, on a profile we understand, instead
  of watching a runaway. A localised discrepancy is diagnosable; a runaway is
  not.
* it needs no new solver, no new circuit topology, and no new physics. Only an
  initial-condition utility.

Note this also finally uses the figure-3 digitisation for the purpose it was
done for.

**Risk, stated in advance:** the seed will not be exactly self-consistent --
digitised profiles carry a few percent of error, `nEps_e` must be seeded from
the LFA equilibrium at the seeded field, and the ion and electron profiles were
digitised independently so their difference (the space charge) is noisier than
either. Expect an initial transient of a few `tau_diel`; judge drift only after
it settles, and judge it on the structural discriminators, not on peak values.

## 3. Recommended order

**REVISED after 2c, which changed the answer.**

1. **SEED FROM GRUBERT'S FIGURE 3 and run the existing transient solver**
   (section 2d). Cheapest by far, needs no new solver or circuit, uses the
   circuit already validated as stable at the operating point, and answers the
   question we actually care about: is the reference profile a fixed point of
   our model?
2. **Steady + current-imposed solver** (section 1). Still worth building --
   it is the only route to a V-I characteristic, it removes the dependence on
   a good initial guess, and the `floatingElectrode` superposition machinery
   already does the hard part. But it is no longer the FIRST thing to try.
3. **`seriesRL`/`seriesRLC`** only if 2d drifts in a way that a faster circuit
   would fix, and only after `nu'` is measured along the actual profile so the
   search is conclusive. Worth having as a capability regardless.

Both routes end with the same validation: the centreline profiles against
Grubert figure 3 via `extract_centreline.py` and `plot_grubert_profiles.py`,
judged on the structural discriminators recorded in
`validation/grubert2009_fix_n11/COMPARE.md` -- density-peak location, field
localisation, and the sign structure of `chargeDensity`.

## Provenance and corrections

All numbers measured 2026-09-06 from `grubert2009_fix_n11` unless stated.
Corrections made while writing this, recorded so the wrong versions are
recognisable if they surface elsewhere:

* the ballast's adequacy was first judged by `tau_RC` vs the TIMESTEP. The
  right comparison is `tau_RC` vs the instability e-folding time.
* a series L was first sized at 1.09 mH from `L*dI/dt = V`, which is the wrong
  criterion; the impedance criterion gives 68 mH, and neither is defensible
  without the Nyquist analysis.
* `tau_gap = 0.28 ns` supersedes any earlier suggestion that the gap itself is
  unstable. It is not; the external circuit is what prevents its self-quench.
