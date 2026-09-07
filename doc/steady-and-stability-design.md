# Reaching the Grubert operating point: two routes, and the analysis that decides

Status: DESIGN for review, 2026-09-06. Nothing implemented.
Written in response to two objections from the user, both of which change the
design.

**AMENDED 2026-09-07: sections 1, 2d and 3 are SUPERSEDED.** The route is now
**steady + BALLAST**, not steady + current-imposed. The user asked why the
ballast circuit we already have cannot be used with a steady solver; it can, it
needs no new code, and the argument is in section 4. Sections 0, 2, 2a-2c
(including the Routh-Hurwitz result, which is what PROVES the ballast adequate)
stand unchanged.

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

> **SUPERSEDED 2026-09-07 BY SECTION 4, in its CONCLUSION only.**
> What still holds: the user's objection itself (pseudo-marching cannot rescue
> an unstable fixed point), and that a fixed-VOLTAGE gap on a `dV/dj < 0` branch
> is unstable. What is WRONG: the inference that current-imposed control is
> therefore REQUIRED. That treated "voltage-imposed" and "current-imposed" as
> the only two options, and **a ballast is neither** -- it is the one-parameter
> family between them, stable whenever `R > |dV/dI|`, which section 2c's own
> table confirms for `R = 1e8` at `gamma = 0`. Also OVERSTATED here: that the
> `floatingElectrode` superposition "already does the hard part". It does not --
> its closed form depends on `Q(V_f)` being LINEAR, and `I(V_f)` is not. See
> section 4.

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

> **SUPERSEDED 2026-09-06 (same day) by the attempt itself**, recorded in
> `validation/grubert2009_seed/COMPARE.md`: the construction is CIRCULAR.
> Current continuity plus Poisson is UNDERDETERMINED -- the third steady
> constraint is PARTICLE BALANCE -- and imposing all three IS the steady solve
> this was meant to avoid. Measured: the gap voltage collapsed 394 V -> 13.8 V
> and the cathode fall was destroyed.
>
> **BUT REHABILITATED 2026-09-07 FOR A DIFFERENT USE.** The digitisation fails
> as a TRANSIENT seed because `rho = e(n_Arp - n_e)` differences two curves that
> agree to a few percent, so the field is ~100% wrong and the solver correctly
> refuses to step through the relaxation. A STEADY solve needs no consistent
> state -- only a guess in the BASIN OF ATTRACTION. The same digitisation is
> therefore a perfectly good INITIAL GUESS for section 4, and the objection
> that killed it here does not apply there.
>
> The no-run result in `grubert2009_seed/COMPARE.md` also stands and is worth
> more than the seed was: Grubert's densities with OUR transport tables carry
> his current at 394 V against his 500 V (0.79x), with the field-free negative
> glow (12-30 Td) emerging by itself.

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

> **SUPERSEDED 2026-09-07 BY SECTION 5.** Item 1 died the same day it was
> written (see 2d). Item 2's premise was wrong (see 1 and 4). Item 3 is
> unchanged and still correct. Kept here so stale copies of this order are
> recognisable.

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

---

## 4. WHY THE BALLAST WORKS WITH A STEADY SOLVER (2026-09-07)

Raised by the user: *"why can't we use the ballast circuit with a steady solver
instead? why do we need the current-imposed solve?"* It can, and we do not.
This section supersedes the conclusion of section 1.

### 4a. The ballast's fixed point IS the steady load line, already

`plasmaExternalCircuit::seriesResistor` updates

    V_target = (V_src - R*I_cond + (Rg + a)*V) / (1 + Rg + a),
    Rg = R*|g|,   a = tau/dt,   tau = R*(C_gap + C_ext)

Setting `V_target = V` cancels `Rg` and `a` identically:

    V*(1 + Rg + a) = V_src - R*I_cond + (Rg + a)*V   =>   V_gap = V_src - R*I_cond

**`Rg` and `a` are DAMPING ONLY** -- they set the rate of approach, never the
location of the fixed point. And `a -> 0` as `dt -> inf`. So under
pseudo-transient marching the ballast we already have, UNCHANGED, is a damped
Newton iteration onto the correct steady load line, with contraction factor
`1/(1 + Rg)` in `(0, 1]`.

Consequence: **steady + ballast needs no new circuit code and no new electrode
constraint.** The only new thing is the steady solve of the plasma equations.

### 4b. The section-1 dichotomy was false

Fixed voltage and imposed current are the two ENDS of a one-parameter family,
and the ballast is the interior:

| R | load line | on a `dV/dI < 0` branch |
|---|---|---|
| 0 | horizontal | UNSTABLE -- this is section 1's objection, and it is right |
| finite | slope `-1/R` | **STABLE iff `R > |dV/dI|`** -- the classical criterion |
| inf | vertical | stable; and this limit IS `currentSource` |

Measured 2026-09-06: `|dV/dI| = 1/|g| = 1.7e6 Ohm` against `R = 1e8`, a **59x
margin**, and section 2c's Routh-Hurwitz table independently reports `R = 1e8`
**STABLE at `gamma = 0`** for both `C = 5e-14` and `C = 0`. That table's
instabilities are all at `gamma = 1.47e9`, the RUNAWAY linearisation point,
which a steady solver never visits. Consistent with the recorded unification
`I_sc/I_op = 1 + V_gap/(R*I_op)`.

**Use `R ~ 3e8`, not `1e8`.** A 300 V ballast drop at `I_op = 1.022 uA` needs
`R = 2.94e8`; `1e8` is 3x faster than the operating point demands
(`validation/grubert2009_spike/COMPARE.md`, 2026-09-07).

### 4c. What `floatingElectrode` does and does NOT transfer

Section 1 claimed its superposition "already does the hard part". Checked
against `floatingElectrode.H` on 2026-09-07 -- it does not.

Its closed form exists because Poisson is LINEAR in `V`, so `Q(V_f)` is exactly
linear, `C_self = dQ/dV_f` is exact from one homogeneous solve, and ONE
correction lands on target. Swap the constraint to `INT j.n dA` and that is
gone: `j` depends on the densities, `mu(E/N)` is a table lookup and `alpha(E/N)`
is exponential, so `I(V_f)` is strongly nonlinear -- which is not incidental, it
IS the negative differential resistance. `C_self` is also the WRONG sensitivity
for this: it is `dQ/dV_f` (displacement), not `dI_cond/dV_f` (conduction).

TRANSFERS: (i) the architectural split -- the BC only CARRIES the equipotential
value while a separate class DETERMINES it, which is exactly the shape a
current-imposed electrode needs; (ii) the actuator `V += dV_f*psi`, valid for
ANY `dV_f` because `L(psi) = 0`, so only the CHOICE rule changes; (iii) `psi` /
`C_self` as the exact DISPLACEMENT part of the sensitivity; (iv) the
self-convicting invariant that reports closure in VOLTS and catches `psi` being
built from a different operator than the solve.

The closer relative for the nonlinear part is `currentSource`, whose secant
`g = dI/dV` from the last two accepted steps is the piece `floatingElectrode`
cannot supply.

### 4d. Why `g` being poorly known stops mattering

`g` is the known weak spot: `093ffcb` found the regulator was an undamped
integrator, and section 2c's search was inconclusive because `nu'` was known
only to ~1.95x. **In a steady solve `g` sets the convergence RATE, not the
answer** -- shown by the algebra in 4a, where the fixed point is independent of
`Rg`. In a transient, by contrast, `g` enters the trajectory itself. This
inverts the earlier concern and is the strongest argument for this route.

### 4e. What current-imposed still buys, and why it can wait

Both are one scalar unknown plus one scalar constraint, so it costs the SAME to
build; it is a second constraint row, not a second machine.

1. **Branch uniqueness.** A vertical load line cuts a non-monotonic
   characteristic exactly once; a sloped one can cut the Townsend, subnormal and
   glow branches, so the initial guess selects the solution.
2. **Setting the validation target directly**, and V-I sweeps. Grubert's
   `j = 0.511 mA/cm^2` becomes an input rather than something `(V_src, R)` is
   tuned to hit; and sweeping `I` is far better conditioned than sweeping
   `V_src` at large `R`, where `dV_gap/dV_src` is small.

### 4f. What is NOT made easier

The steady residual must include **PARTICLE BALANCE** -- species and energy
converged simultaneously with Poisson and the circuit relation. That is the
third constraint whose absence made the 2d seeding circular, and it is the
actual work. Choosing the ballast removes the CIRCUIT question from the scope;
it does not shrink the steady-solve question.

## 5. RECOMMENDED ORDER (2026-09-07), superseding section 3

Escalating, cheapest first, each phase with the observable that decides whether
the next is needed. Phases 1-2 need NO new solver code.

**Phase 0 -- fix the initial guess (no run).** Build the `rho ~ 0` variant:
digitised `n_Arp` everywhere, digitised `n_e` only inside the cathode fall,
`n_e = n_Arp` in the bulk, transition at `x/L ~ 0.25`. Judge it on the LOCAL
field, never the integrated voltage -- that is the trap that hid the last one.
Note this is a GUESS, not a seed: per 2d it needs only to be in the basin.

  NEEDS A NEW SCRIPT. Checked 2026-09-07: `grubert2009_seed/build_consistent_seed.py`
  implements the OTHER, FAILED construction (current continuity + Poisson, the
  one that collapsed to 13.8 V), not this one. The `rho ~ 0` variant is far
  simpler -- no iteration at all, just a piecewise assignment.

  WHY IT IS SOUND: in the bulk, quasi-neutrality holds to `(lambda_D/L)^2`, so
  `n_e = n_Arp` there is MORE right than differencing two digitised curves. In
  the fall the densities differ by DECADES, so their difference is well
  conditioned and the digitisation is usable. Each region uses the constraint
  that is accurate there.

  PASS/FAIL, available with no run: the bulk field must come out at
  **12-30 Td** and the fall must carry most of the ~400-500 V. Those are the
  independent no-run numbers from `grubert2009_seed/COMPARE.md`, where the
  field-free negative glow emerged by itself from current continuity. A seed
  that misses them is wrong before it is ever launched.

### PHASE 0 EXECUTED 2026-09-07: PASSES, but `x/L = 0.25` WAS WRONG

The check was run before any case was built, and it caught the transition point.

**`XT = 0.25`, the value in the 2026-09-06 note, FAILS.** It leaves a **662 Td**
bulk field, because with `rho == 0` the bulk field is CONSTANT, so whatever
voltage the fall cannot carry is forced into the bulk uniformly. At 662 Td,
`1/nu_i = 5.81 ns` against a ~6 ns electron transit across the 6 mm bulk --
**the bulk would avalanche**, which is exactly the failure the seed must avoid.

**THE DATA SAYS THE TRANSITION IS AT `x/L ~ 0.40`, not 0.25.** Where the
digitised curves actually meet: `n_e/n_Arp` = 0.306 at 0.30, 0.718 at 0.35,
0.919 at **0.40**, 0.969 at 0.45, 0.986 at 0.50.

**AND THE ANSWER IS ROBUST THERE.** The fall-charge field swing SATURATES --
13,098 Td at XT 0.35, 13,277 at 0.40, 13,326 at 0.45 (a 0.4% spread) -- because
beyond `x/L ~ 0.4` the digitised difference is negligible. A plateau is the
signature of the quasi-neutral region beginning, so XT is no longer a free
choice. Contrast XT 0.20 -> 0.25 -> 0.30, where the swing moves 8511 -> 10808
-> 12416 Td and the implied voltage moves 221 -> 345 -> 451 V. **The earlier
0.25 sat on the steep part, which is why it produced a number at all.**

**THE CONSISTENCY RESULT, and it is better than either previous estimate.** In
the robust range a field-free (20 Td) bulk requires

    Vgap = 519.4 V (XT 0.40) .. 524.3 V (XT 0.45)   vs Grubert's 500 V

i.e. **within 4-5%**, from Grubert's digitised densities + Gauss's law + his
stated voltage, and NOT circular -- the physically-defensible XT range maps onto
a narrow voltage band that contains his value. Compare the two independent
earlier estimates, both LOW: 393.9 V from drift-only current continuity (0.79x)
and 344.9 V from this construction at the wrong XT = 0.25 (0.69x).

**THE IRREDUCIBLE LIMIT, stated so it is not rediscovered.** The bulk field is a
20-60 Td residual of a 13,277 Td swing -- **0.15-0.45%** -- so it is the
DIFFERENCING PROBLEM ONE LEVEL UP: not `n_Arp - n_e` pointwise, but
`E0 - INT rho`. Imposing exactly 500 V at XT 0.40 gives ~60 Td REVERSED rather
than 20 Td forward. The seed's bulk field therefore carries ~+-80 Td of
irreducible uncertainty and CANNOT be controlled by this construction.

**Tolerable, and that is the actual pass criterion.** Across the whole +-80 Td
band `1/nu_i` runs 4.06 us (60 Td) to 0.301 s (12 Td), all >> the 6 ns transit,
so nothing avalanches anywhere in the bulk. The guess is in the basin, which is
all 2d requires of it. Verified against the ionisation table with the reader
cross-checked on two recorded values (14.2 ns at 431 Td vs 14.19 recorded;
0.731 ns at 2450 Td vs 0.73).

RECORDED TRAP: `k_*_vs_reducedE` tables carry the axis in **SI, V m^2**, not Td.
Read as Td every lookup silently returns the table MAXIMUM (2.728e-13 m^3/s),
which looks like a plausible rate. Multiply Td by 1e-21.

**Phase 1 -- existing transient solver, ballast `R = 3e8`, from that guess.**
Zero new code. Question: does it SETTLE, and is what it settles to Grubert's?
Discriminating observables: `n_e` peak magnitude and location, cathode-fall
`E/N`, and `j` against `0.511 mA/cm^2`. FAILURE SIGNATURE, already measured
once: `dt` collapses (the accuracy controller refusing a violent relaxation) ->
the guess is still inconsistent, go to phase 2.

**Phase 2 -- pseudo-transient, only if phase 1 will not settle.** The path is
not wanted, so stop paying for its accuracy: error control off, large fixed
`dt`. A CONFIG change, not code. `a = tau/dt -> 0` makes the ballast update
pure damped Newton (4a).

**Phase 3 -- a real steady solver, only if 1-2 both fail.** Drop the `d/dt`
terms and solve the coupled residual, with the ballast relation as one scalar
row. This is where the work is (4f), and the outer-coupling conditioning item
in memory becomes a blocker rather than an optimisation.

**Phase 4 -- current-imposed as a second constraint row** (4e), for V-I
characteristics and branch uniqueness. Independent of 1-3.

Every phase: a `COMPARE.md` written when the case is created (rule 12), probes
and a machine-readable time series (rule 26), and validation on the structural
discriminators in `validation/grubert2009_fix_n11/COMPARE.md`.
