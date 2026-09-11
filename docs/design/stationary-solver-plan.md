# A stationary (pseudo-steady) solve path for SoPLASMA

**Status: PARTLY SUPERSEDED, 2026-09-11. Section 2's "no new code" first step
DOES NOT EXIST.** The plan proposes `ddtSchemes { default steadyState; }` as a
case setting. It was tried the DAY AFTER this plan was written and the solver
now refuses it outright (`plasmaTimeControl.C`, commit `0ccdea1`, 2026-09-09):

> "The ddt term IS the diagonal of a segregated species transport equation.
> Removing it leaves rows with no diagonal for any species without an implicit
> loss term, and the linear solver then divides by zero on its first sweep
> (observed: SIGFPE in GaussSeidelSmoother::smooth). simpleFoam can drop ddt
> because SIMPLE's pressure equation and under-relaxation supply the diagonal
> dominance; nothing here does."

So under-relaxation would have to supply the diagonal itself (OpenFOAM's
`fvMatrix::relax()` does exactly that -- it divides the diagonal by the factor
and moves the remainder to the source), and the guard fires before relaxation
is ever applied. Making the stationary path real therefore needs SOLVER WORK,
not case settings: at minimum, permit `steadyState` when relaxation is active
AND ensure `relax()` reaches the species equations.

Section 2 item 4 (CURRENT CONTROL) was checked on 2026-09-11 and stands, but is
NOT sufficient on its own: all six `grubert2009_iset*` arms are current
controlled and all stalled at t ~ 9.5e-07, one of them after 80,983 steps. They
did not fail for a fixable reason and were not abandoned early.


Written 2026-09-08 after Almeida, Benilov, Cunha & Gomes, *Computing Different
Modes on Cathodes of DC Glow and High-Pressure Arc Discharges: Time-Dependent
Versus Stationary Solvers*, Plasma Process. Polym. **14**, 1600122 (2017)
[Literature/Almeida_2016_CathodesModelingDCGlowandArc.pdf].

---

## 1. WHY. The Almeida findings, and what they say about our two days of failure

This paper is about the exact problem we have been failing at, in someone
else's code, and it concludes that the tool is wrong.

**(a) The 1D DC glow CVC has a MINIMUM and HYSTERESIS.** In argon at 120 Torr,
0.5 mm gap, they find hysteresis between the subnormal and abnormal discharge
over 3 kA/m^2 <~ j <~ 9 kA/m^2, and remark that such hysteresis "is unusual in a
1D DC discharge."

**(b) A mature commercial time-dependent solver fails there exactly as ours
does.** Of COMSOL Multiphysics' Plasma module:

  > "The convergence was lost shortly before the minimum of the CVC."
  > "the straightforward application of the Plasma module of COMSOL
  >  Multiphysics allows one to readily compute only the ABNORMAL discharge."

**(c) THE FIX IS CURRENT CONTROL, NOT A BALLAST.**

  > "A more robust code can be built by adding the possibility of using the
  >  discharge current density j as a control parameter (instead of the
  >  discharge voltage U) ... without expressly introducing a ballast
  >  resistance."

  WHY, and this is the part that explains everything: at the CVC minimum the
  differential resistance passes through zero. NO voltage-plus-ballast load
  line can select an operating point there -- the load line is tangent to the
  characteristic. Current control can. This is a property of the problem, not a
  numerical preference.

  It also retroactively answers the question asked at the very start of this
  effort -- "why can't we use the ballast circuit with a steady solver instead
  of the current-imposed solve?" Because of the CVC minimum.

**(d) Time-dependent solvers are the wrong tool for these steady states.**

  > "Of all the modes considered, ONLY ONE could be computed in the whole
  >  region of its existence without gaps by means of the time-dependent
  >  solver: the 1D mode of DC glow discharge in xenon in the framework of the
  >  local-field approximation."
  > "In all the other cases, the modeling results are incomplete."
  > "the stationary solver is capable of computing all modes in the whole range
  >  of their existence and is, therefore, A TOOL OF CHOICE."
  > "the use of stationary solvers is generally simpler ... there is more
  >  freedom in the choice of spatial mesh and NO NECESSITY TO RESOLVE
  >  DIFFERENT TIME SCALES."

  And from the abstract, the sentence that reframes our whole effort:

  > "numerical stability of the time-dependent solver was NOT EQUIVALENT to
  >  physical stability."

**(e) WHAT THIS EXPLAINS ABOUT OUR RUNS.** Every arm died the same way --
150/150 correctors, dt collapsing ~25000x to ~4e-14 -- regardless of voltage
(300/400/250/200 V), reflection coefficient (0 / 0.36) or flux scheme (standard
/ SG / CFS). We treated it as a scheme problem, then a voltage problem. It is
neither. Lowering the drive to 1.66x breakdown bought 0.18 us; CFS bought 1.7x.
Neither addressed the obstacle, because the obstacle is that we are
time-marching toward a state that voltage control cannot select.

Note also Grubert 2009 computed "at steady state" with a finite-element
formulation. The nanosecond avalanche that destroyed our runs may be a
transient GRUBERT NEVER TRAVERSED.

---

## 2. WHAT. A simpleFoam-style stationary path

OpenFOAM already carries the pattern. `simpleFoam` is documented as
"Steady-state solver for incompressible, turbulent flows" and is driven by
`simpleControl::loop()`:

    while (simple.loop())      // iterate; on convergence -> writeAndEnd()
    { ...assemble and solve...  }

with `simpleControl::criteriaSatisfied()` testing the per-field
`residualControl` entries in the SIMPLE dictionary. Steady state itself comes
from two case-level settings, not from solver code:

    system/fvSchemes:   ddtSchemes { default steadyState; }   // drops d/dt
    system/fvSolution:  relaxationFactors { ... }             // under-relaxation

### Proposed ingredients for SoPLASMA

  1. **`ddtSchemes steadyState`** on the species and electron-energy equations.
     `fvm::ddt` then contributes nothing and each equation becomes its own
     steady balance. This is a CASE setting; no solver change needed for it.

  2. **Under-relaxation** on n_e, the ion densities, nEps_e and ePotential, via
     `relaxationFactors` in fvSolution. This replaces the timestep as the
     mechanism that limits how far the nonlinear system moves per iteration --
     which is precisely the thing that collapsed to 4e-14 in the transient runs.

  3. **A `simpleControl`-style outer loop** as an alternative to the existing
     PIMPLE loop, with residual control on the plasma fields. SoPLASMA's outer
     loop already exists (`maxCorrectors`, `outerCoupling`, residualControl) --
     the work is to let it run without a time axis and to converge on residuals
     rather than on a temporal-error criterion.

  4. **CURRENT CONTROL, and this is the essential one.** `plasmaExternalCircuit`
     already supports `type currentSource` and six `grubert2009_iset*` cases
     exist from an earlier attempt in this effort. Per Almeida (c), the
     stationary solve MUST be current-controlled to traverse the CVC minimum.
     Re-examining those `iset` cases is the cheapest possible first step and
     should precede any new code: they were the right approach, and it matters
     whether they failed for a fixable reason or were abandoned early.

### What this does NOT give

  * The ignition transient. A stationary solve finds steady states; it says
    nothing about how the discharge gets there. That is fine -- the steady
    state is what Grubert reports and what the Engel & Steenbeck similarity
    data describes.
  * Stability information. Almeida are explicit that a stationary solver
    computes unstable states too, and that "numerical stability of the
    time-dependent solver was not equivalent to physical stability" -- so
    stability must be assessed separately, not inferred from either solver.
  * Continuation. Traversing a hysteretic CVC needs parameter continuation
    (step j, use each converged state as the initial guess for the next), which
    is how Almeida generate their curves. That is a driver-level feature.

---

## 3. Suggested order

  1. Read the `grubert2009_iset*` cases: what was imposed, how far they got,
     why they stopped. NO new code.
  2. If current control is sound, try `ddtSchemes steadyState` + relaxation on
     an existing case with the current-imposed circuit. Still no new code --
     both are case settings.
  3. Only then consider a `simpleControl`-based driver and j-continuation.

Step 2 is the cheap experiment that would tell us whether a pseudo-stationary
solve is viable in SoPLASMA at all.
