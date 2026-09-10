# Non-local electron kinetic theory: read carefully, assessed against Grubert -- worth prototyping, NOT worth committing to yet

2026-09-08. All six papers read (Bernstein-Holstein 1954 excluded -- paywalled,
not fetched; the other five plus a bonus Tsendin 2010 review, all now in
`Literature/non-local-kinetics/`). This follows directly from
`docs/design/electron-electron-collisions-gap.md`, which identified non-local kinetic
theory as the documented middle ground between SoPLASMA's local fluid
closure and full PIC for the negative-glow trapping problem.

## The method, in one paragraph

Bernstein & Holstein (1954) and Tsendin (1974) showed that once the electron
energy relaxation length exceeds the system size, the electron distribution
function becomes a function of TOTAL energy `eps = u - Phi(r)` (kinetic minus
the local potential) rather than of position and velocity separately. Averaged
over the spatial volume accessible to electrons of a given total energy, the
kinetic equation collapses to a single ORDINARY differential equation in
`eps` (Kortshagen et al. 1996, eq. 48) -- REGARDLESS of how many spatial
dimensions the actual problem has. The spatially resolved EDF at any point is
then recovered algebraically: `F0(u, r) = F0^(0)(eps = u - Phi(r))` (eq. 58).
Electron-electron collisions enter this equation as an additional term
(eq. 54-57) with no extra structural cost -- they do NOT need the separate,
harder anisotropic treatment SoEEDF's Coulomb model needs, because the
non-local equation is already a 1D energy-space problem, not a two-term
(F0, F1) system.

## Quantitative applicability check against Grubert 2009

Kortshagen et al. (1996) give a hard numerical criterion for argon, derived by
comparing the non-local and local approximations against a numerical solution
of the FULL spatially-resolved kinetic equation in a positive column:

    N0*R <~ 3e21 m^-3 cm   -- nonlocal approximation valid
    N0*R >~ 1e23 m^-3 cm   -- local approximation valid
    in between             -- "the solution of the complete spatially
                              dependent kinetic equation is inevitable,
                              and is a much more complicated problem"

Grubert's case: 100 Pa argon (N0 = 2.414e22 m^-3), 1 cm gap:

    N0 * R = 2.414e22 m^-3 cm

This is **8x above the nonlocal threshold and 4x below the local threshold**
-- squarely in the INTERMEDIATE regime, where neither cheap approximation is
valid by this criterion and the full 2D-in-(space, energy) kinetic PDE is
needed instead (which is NOT cheaper than PIC in any obvious way -- it is
itself an elliptic boundary-value problem with an irregular domain, per
Kortshagen sec. 3.3).

CAVEAT, stated plainly: this specific N0*R threshold was derived for a
POSITIVE COLUMN (bulk, radially confined by ambipolar diffusion), not for a
cathode fall / negative glow. The relevant confinement length for Grubert's
negative glow is more likely its own thickness, L_NG, not the full 1 cm gap.
Using Boeuf's own reported L_NG ~ 4 mm for a comparable He case (not
Grubert's argon case directly): N0*R_NG = 2.414e22 * 0.4 = 9.66e21 m^-3 cm --
still 3.2x above the nonlocal threshold, i.e. still intermediate, just less
severely so. Neither number is a rigorous verdict for Grubert's actual
geometry and gas; both point the same direction: **this is not a clean
nonlocal case, even optimistically.**

## Two further reasons for caution, from the modern reviews

1. **Transient extension is explicitly unsolved.** Kolobov & Godyak's own
   2019 perspective (Phys. Plasmas 26, 060601) lists "nonlocal and transient
   effects" as ONE OF FOUR open research frontiers in the field, not a solved
   problem. Grubert's case -- and every DC-glow case we run -- is a genuine
   IGNITION TRANSIENT, not a steady positive column. The classical theory
   (as validated in Kortshagen 1996) is developed and checked for STATIONARY
   or slowly-varying conditions; the authors of the modern reviews do not
   claim the transient case is solved. Building this for our actual use case
   would mean doing open research, not implementing an established recipe.

2. **The full-generality reference case in Kortshagen 1996 EXCLUDES
   electron-electron collisions** ("the implementation of the
   electron-electron collisions, with the nonlinearity of the collision
   integrals, is beyond the present scope" -- their sec 3.3, describing the
   full 2D kinetic-equation solver used as the accuracy REFERENCE). So even
   the most rigorous available validation of the non-local approximation's
   accuracy was performed WITHOUT the physics (e-e collisions) that is the
   entire reason we are looking at this in the first place. The accuracy of
   the non-local approximation SPECIFICALLY for e-e-collision-dominated
   trapping has less direct validation in what was read than the accuracy of
   the non-local approximation in general.

## A genuinely useful, much cheaper idea that fell out of this reading

Kolobov 2013 (Phys. Plasmas 20, 101610), section D, gives a three-electron-
group picture of the cathode region that matches Boeuf's JC-PIC description
exactly: (i) fast beam/runaway electrons from the sheath, (ii) "intermediate"
electrons carrying the negative-glow/FDS current, (iii) TRAPPED electrons.
About the trapped group specifically:

  "The Trapped Electrons do not participate in the current flow, but are
  responsible for a sharp peak of plasma density at the point xm where the
  electric potential has a local minimum. The EEDF of trapped electrons is
  OFTEN MAXWELLIAN, THEIR DENSITY IS DEFINED BY THE BOLTZMANN RELATION, and
  the temperature is slightly above room temperature."

If that characterisation transfers to Grubert's argon case, it suggests a
MUCH cheaper fluid-closure idea than either full non-local kinetics or PIC:
represent the free/beam population with SoPLASMA's existing LMEA
drift-diffusion-reaction closure (which is reasonably suited to that
non-equilibrium, field-driven population), and add a SEPARATE, SECOND
electron population for the trapped fraction, closed with a simple Boltzmann
relation `n_trapped ~ exp(e*Phi/Te_trapped)` relative to the LOCAL potential
minimum, at a temperature near gas temperature -- NOT tied to the bulk
mean-energy variable at all. This is a "two-group" fluid extension, not a
kinetic module: no energy-space grid, no non-local spatial averaging, no
elliptic PDE.

THIS IS MY OWN SYNTHESIS, NOT A RECIPE STATED IN ANY OF THE SIX PAPERS. None
of them propose this as a closure for a fluid code -- they describe it as a
feature of PIC/kinetic results. It is offered here as a scoped, testable
hypothesis, not a validated method. Before committing effort: check whether
this two-group Boltzmann-trapped-electron ansatz, evaluated as a POST-HOC
diagnostic on the existing PIC/kinetic literature data (Carlsson benchmark,
JC-PIC's own He case), actually reproduces the reported negative-glow
density -- a cheap, decisive test that requires no new SoPLASMA code at all.

## Bottom line

  * The theory is real, well-established (30+ years), and would be a
    genuine capability if built -- NOT snake oil, NOT overhyped.
  * It is NOT free: Grubert's own case sits in the awkward intermediate
    regime by the one quantitative criterion found, and the transient
    extension needed for an ignition problem is open research, not
    engineering.
  * Building it now would be a RESEARCH PROJECT (comparable in scope to,
    though smaller than, a PIC/hybrid capability), not a quick win.
  * A cheaper, testable alternative exists and should be tried FIRST: the
    two-group (LMEA free population + Boltzmann-closed trapped population)
    idea above, checked against literature data before writing any code.
  * Recommendation: do NOT start a non-local kinetics module now. Keep it as
    a documented, scoped option (alongside the PIC/hybrid idea in
    `future-pic-hybrid-capability`) for if/when the negative-glow gap
    becomes an active blocker after the Grubert circuit/floor/current fixes
    already made are given a fair chance to run their course.

## UPDATE 2026-09-08 (same day): the two-group idea is NOT novel -- Boeuf & Pitchford (1995) already did it

Following a direct question from the user on rigor/novelty/whether this is a
"sheath model in disguise", searched for prior art rather than continuing to
speculate. Found via Pinheiro's review (arXiv:physics/0611052, saved to
`Literature/non-local-kinetics/`), which devotes a full section to it:

**J. P. Boeuf and L. C. Pitchford, J. Phys. D: Appl. Phys. 28, 2083 (1995).**
A "simple fluid model" with exactly two electron groups for the negative
glow, giving an analytical formula for the field-reversal location depending
only on cathode sheath length, gap length, and the electron energy
relaxation length lambda_eps.

### The actual formulation (from the review, eq. 8-24)

  * FREE/INTERMEDIATE group: ordinary drift-diffusion continuity equation,
    `dGamma_e/dx = S_T(x)`, with a PRESCRIBED (not self-consistently solved)
    ionization source: zero across the sheath, then exponentially decaying
    into the negative glow with relaxation length lambda_eps, calibrated
    against separate Monte Carlo runs.
  * TRAPPED group: an ORDINARY AMBIPOLAR DIFFUSION EQUATION,
    `d/dx(D_a dn/dx) + I(x) - n/tau = 0` -- isothermal, source I(x), loss
    n/tau. NOT a Boltzmann relation. Structurally the same kind of equation
    SoPLASMA already solves for every quasineutral species.
  * The two groups are LINKED by exchange source/sink terms (the free group's
    exponential-decay ionization deposit standing in for the free->trapped
    rate; the modern picture would replace this with a proper
    trapping-rate model, e.g. via e-e collisions, which this 1995 model does
    not include at all).

### Answering the three questions directly

  1. RIGOR: yes -- it is a well-posed two-species continuity system, not
     hand-waving. The gap is that the free->trapped exchange rate in the
     1995 model is an EMPIRICAL, offline-calibrated ansatz (the exponential
     decay length lambda_eps), not derived from a first-principles local
     trapping criterion.
  2. NOVELTY: the CONCEPT is 30 years old, published, and validated (matches
     a Monte-Carlo/fluid hybrid per the review). NOT publishable as "a
     two-group model exists". What COULD be publishable: replacing the
     empirical exponential-decay exchange rate with a first-principles local
     trapping criterion (e.g. comparing an electron's total energy against
     the running extremum of the SELF-CONSISTENT, evolving 1D potential
     profile between its position and each boundary -- an O(N) array scan
     alongside the existing Poisson solve, cheap, and NOT requiring the full
     non-local energy-space kinetic machinery), embedded in SoPLASMA's
     general TRANSIENT multi-species framework (the 1995 model is 1D,
     steady-state only -- it cannot address an ignition transient at all),
     and validated against modern kinetic benchmarks (Carlsson, JC-PIC) that
     did not exist in 1995. That is a genuine, scoped, checkable research
     contribution: a from-first-principles, transient generalization of a
     known steady-state analytical result.
  3. "SHEATH MODEL IN DISGUISE": NO, not quite. The trapped population in
     the actual precedent is closed by ambipolar diffusion, not a Boltzmann
     relation -- i.e. it is much closer to "a second species with its own
     already-familiar transport equation" than to a repurposed sheath
     closure.

### Revised recommendation

Given the user is weighing this explicitly as a RESEARCH direction (not only
framework engineering) -- research is an equally valid and desired outcome,
not something to be talked out of -- this is now assessed as a WELL-SCOPED,
PRECEDENTED, ATTRACTIVE research candidate: a validated 1995 starting point
to modernize, a concrete idea for the one piece that was empirical
(first-principles local trapping criterion), and two modern benchmarks to
validate against. Distinct from, and considerably cheaper than, both full
non-local kinetics (see the original assessment above -- Grubert's own case
sits in the awkward intermediate N0*R regime for THAT approach) and full
PIC/hybrid. Worth scoping as an actual project if the user wants to pursue
it, rather than being deferred indefinitely.

## UPDATE 2026-09-08 (later same day): all five two-group/hybrid papers now read IN FULL

Per the new rule (32, `full-paper-before-conclusions`), all four requested
papers plus the already-summarized Boeuf & Pitchford were read completely,
not via search summaries. This corrects and completes the picture above.
Full details: `Literature/non-local-kinetics/README.md`.

**The lineage, confirmed by direct reading:**
Kolobov & Tsendin, Phys. Rev. A 46, 7837 (1992) [foundational, three-group
kinetic theory] -> Boeuf & Pitchford, J. Phys. D 28, 2083 (1995) [simpler,
independent, single-population analytic result] -> Rafatov, Bogdanov &
Kudryavtsev, Phys. Plasmas 19, 033502 (2012) [plain LMEA reliability --
cites Grubert 2009 as the same lineage] -> Eylenceoglu, Rafatov &
Kudryavtsev, Phys. Plasmas 22, 013509 (2015) [two-group hybrid, Te still
fitted] -> Eliseev, Bogdanov & Kudryavtsev, Phys. Plasmas 24, 093503 (2017)
[adds the self-consistent slow-electron energy balance, closing the gap].
One continuous, identifiable research school (Tsendin/Kolobov/Kudryavtsev
and students), 1992-2017+.

**Correcting the earlier read of Boeuf & Pitchford**: it is NOT two coupled
electron species with separate transport equations. It is ONE ambipolar
equation with a prescribed exponential ionisation source, and it explicitly
states a self-consistent Te equation is "beyond the scope of this paper."

**The rigorous theoretical foundation (Kolobov & Tsendin 1992)**: three
kinetic groups (fast/intermediate/trapped). The trapped population's EDF is
PROVEN Maxwell-Boltzmann when nu_ee >> delta*nu (a checkable condition, not
assumed), with its temperature from a genuine integral energy-balance
equation. This is the actual rigorous answer to "how do the populations
link" -- confirming the physical picture originally proposed, sourced now
to its real origin.

**The single most actionable result (Eliseev 2017)**: NOT a new species at
all. A single Coulomb-heating SOURCE TERM added to the electron energy
equation SoPLASMA already solves, with a closed-form formula depending only
on local n_e, T_e, and known cross-sections. Validated against real probe
data in ARGON at 57-107 Pa -- Grubert's own 100 Pa almost exactly. Central,
quantitative, falsifiable result: omitting this term forces an artificially
high negative-glow Te (2-4 eV vs the real ~0.2-1 eV) to sustain the
ionisation balance, and since n_e,max*Te is found to be ~pressure-
independent, this DIRECTLY predicts an order-of-magnitude underestimate of
negative-glow density -- a specific, checkable signature to look for in any
future SoPLASMA Grubert run that reaches a steady negative glow.

**A genuine caution (Eylenceoglu 2015)**: their own two-group hybrid model
(run at gamma = 0.06, the same value we use) found NO significant advantage
over plain LMEA for INTEGRAL characteristics (CVC) -- both are dominated by
the fitted gamma value. The advantage of the two-group treatment is
specifically in LOCAL structure (density/field profiles), not integrated
current-voltage behaviour. This tempers how much a two-group (or
Coulomb-term) fix should be expected to change CVC-level agreement on its
own.

**Revised novelty verdict**: the general two-group/hybrid concept is NOT
novel -- confirmed by full reading, not just search results. What remains
undone, as far as these five papers show: implementing the Eliseev 2017
Coulomb-heating term (or an equivalent) inside a general-purpose,
multi-dimensional, TRANSIENT, multi-species OpenFOAM framework, validated
against both this school's classic benchmarks and newer kinetic PIC-MCC
benchmarks (Carlsson, JC-PIC) in the same study. That is the real, scoped,
incremental contribution -- not "inventing two-group models."

Per the standing rule that research is a valid track (not a lower-priority
one), this is assessed as a worthwhile candidate, now with a much more
concrete, minimal, and already-validated starting point than the original
non-local-kinetics or two-species framings suggested: a single source term
in an equation SoPLASMA already solves.
