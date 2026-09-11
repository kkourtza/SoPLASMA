# Steady / pseudo-steady mode for SoPLASMA -- design spec

2026-09-08. Written after measuring what does and does not work. Motivated by
Almeida et al 2016 (see `docs/design/stationary-solver-plan.md` and the memory entry):
a 1D DC glow cannot be reached reliably by voltage-driven time-marching, and a
stationary solver is "a tool of choice".

**HEADLINE: a TRUE `steadyState` mode is a SOLVER PROJECT, not a case setting.
A PSEUDO-TRANSIENT mode works TODAY with existing settings and is the thing to
expose first.**

---

## 1. What was measured

### 1.1 A working pseudo-transient recipe (no code changes)

    ddtSchemes        backward     # BDF2. Measured 5.5x faster than Euler here
    adjustTimeStep    false        # dt CANNOT collapse -> guaranteed progress
    deltaT            1e-12        # FIXED; a relaxation parameter, not physics
    onNonConvergence  (omit it)    # `retryStep` is REJECTED when adjustTimeStep
                                   # is off -- it retries by shortening deltaT
    relaxationFactors NONE         # see 1.2 -- do NOT add them
    circuit           currentSource  # per Almeida: voltage control is ill-posed
                                   # near the CVC minimum

RESULT: **9472 consecutive converged steps, 0 non-converged, 4 correctors per
step** against a 150 cap, dt rock-steady at 1e-12. Every previous attempt died
because the adaptive controller shrank dt chasing a transient it could not
resolve, then the degraded-step guard aborted. With `adjustTimeStep false` that
cannot happen: the run either converges each step or fails loudly.

COST: 9.47 ns of simulated time in 180 s, i.e. ~5 h per us. Slow, and 9.5 ns is
far too early to say anything about the physics.

### 1.2 MANUAL `relaxationFactors` BREAK IT. Do not add them.

Adding `relaxationFactors` (fields and equations, 0.3) to the above took it from
**9472/9472 converged** to **0/10 converged**. SoPLASMA already runs ADAPTIVE
AITKEN outer relaxation (`plasmaOuterRelaxation`, enrolled automatically,
observed `omega 0.371`); fixed factors layered on top fight it. This was my own
error, introduced as "the steady-solve counterpart of the timestep", and it cost
several wrong conclusions. **The generator must never emit them.**

### 1.3 A TRUE `steadyState` mode SEGFAULTS, and the reason is structural

    ddtSchemes steadyState + deltaT 1 + chemistry solver implicitRate
      -> SIGFPE on iteration 1

Stack: `plasmaTransport::solve -> fvMatrix::solveSegregated -> smoothSolver
-> GaussSeidelSmoother::smooth` -- i.e. **division by a ZERO MATRIX DIAGONAL**.

WHY. `ddt` contributes `V/dt` to the diagonal. With `steadyState` that is gone,
and the log confirms the implicit chemistry loss contributes nothing either
(`plasmaChemistry: implicitRate, 2000 active cells, max(L*dt) = 0`). For a
species with no loss term the row then has NO diagonal at all: the matrix is
singular and Gauss-Seidel divides by zero on the first sweep.

**`ddt` IS the diagonal of a segregated species-transport equation.**
`simpleFoam` can drop it because the SIMPLE algorithm's pressure equation and
its under-relaxation supply diagonal dominance; a species equation has no such
mechanism. OpenFOAM's *equation* relaxation cannot rescue this either -- it
scales a diagonal that must already exist.

### 1.4 Why `deltaT` and `adjustTimeStep` are read at all

In `simpleFoam` neither matters: `deltaT` is a meaningless iteration counter and
`simpleControl` never reads `adjustTimeStep`. In SoPLASMA dt enters in at least
five places besides `ddt`:

  * `computeChemistrySources(dt)` -- ODE integration for `adaptive`/`ode`
  * `solveGasEnergy(dt)`
  * Courant, dielectric-relaxation and temporal-error limiters
  * the outer-loop retry/rejection logic

`chemistry/solver implicitRate` removes the first (VERIFIED: `max(L*dt) = 0`),
but the rest remain. So dt cannot simply be declared meaningless.

---

## 2. What to build

### 2.1 NOW: expose the pseudo-transient recipe (low risk, tested)

A single generator switch, because the settings INTERLOCK and a user assembling
them by hand hits a fatal error (measured, twice):

    simulationType   transient | pseudoSteady

`pseudoSteady` must emit, coherently:

    ddtSchemes        backward
    adjustTimeStep    false
    deltaT            <pseudoTimeStep>      # exposed, documented as a
                                            # RELAXATION PARAMETER not physics
    onNonConvergence  omitted (or an explicit non-dt-shortening option)
    relaxationFactors NOT EMITTED           # section 1.2
    chemistry/solver  implicitRate          # optional but preferred: removes
                                            # dt from the chemistry
    circuit/type      currentSource         # REQUIRED for the CVC minimum

and it should REFUSE, with a message, any combination known to fail:
`retryStep` with `adjustTimeStep false`; manual `relaxationFactors`;
`steadyState` ddt (until 2.2 exists); voltage control with a ballast if the
user asks for a steady glow near the CVC minimum.

### 2.2 LATER: a true steady mode

Needs a diagonal source to replace `ddt`. Two standard options:

  (a) **Pseudo-time diagonal.** Add `V/dTau` to the diagonal with `dTau` an
      explicit relaxation parameter, decoupled from `runTime.deltaT()` and from
      every Courant/dielectric limiter. This is pseudo-transient continuation
      done deliberately rather than by accident -- and it is precisely why the
      recipe in 1.1 converges.
  (b) **Assembly-level under-relaxation**, the usual `diag/alpha` boost with
      the balance carried to the source. Equivalent in effect to (a); differs
      in how the parameter is exposed.

Either way the outer loop must converge on FIELD RESIDUALS, not temporal error.
Note the existing loop is already `pimple.loop()` with `residualControl`, so
that half largely exists.

### 2.3 Also needed for the actual goal

  * **j-continuation.** A hysteretic CVC (Almeida: 3-9 kA/m^2 in Ar) can only
    be traversed by stepping the control parameter and using each converged
    state as the next initial guess. This is a driver feature, not a scheme.
  * **No stability information.** Almeida are explicit that a stationary solver
    computes unstable states too, and that "numerical stability of the
    time-dependent solver was not equivalent to physical stability". Stability
    must be assessed separately and cannot be inferred from either solver.

---

## 3. Honest status

  * Pseudo-transient + current control: **converges reliably, 4 correctors,
    no dt collapse.** First configuration in this whole effort with a
    plausible path to a steady glow. Not yet run long enough to have one.
  * True `steadyState`: **blocked on a singular matrix**, understood, fixable
    only in the solver.
  * NO cathode fall has formed in ANY run to date, so the Engel & Steenbeck
    similarity comparison remains out of reach.
