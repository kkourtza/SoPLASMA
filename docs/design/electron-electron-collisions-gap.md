# Electron-electron (Coulomb) collisions: SoEEDF has them, SoPLASMA doesn't use them -- an integration gap, and a deeper non-local limit underneath it

2026-09-08. Prompted by reading Boeuf's JC-PIC book (section V.A, DC glow
discharge) -- see `Literature/reference-codes-and-manuals/JC-PIC_Boeuf/`, the standing reference (memory
`jcpic-boeuf-standing-reference`). **Corrects an error made earlier the same
day**: the first version of this document claimed SoEEDF had no e-e collision
term. That was wrong -- it only checked SoPLASMA's fluid-solver source, not
SoEEDF, a separate module in the same multi-module project.

## What is actually true

**SoEEDF implements electron-electron collisions**, both the isotropic
("Rockwood") term and the anisotropic term of Hagelaar 2016 (*Coulomb
collisions in the Boltzmann equation for electrons in low-temperature gas
discharge plasmas*, PSST 25, 015015), in `src/CoulombTerms.C`, exposed as the
`coulomb` option (default `None`), documented as mattering "once the
ionisation degree is above ~1e-5" (`docs/options-reference.md`). It is
validated against the paper's own checks: residual convergence under
refinement, the correct direction and onset ionisation degree of the
mobility shift, and the EEDF anisotropy slope converging to the theoretical
+1/2 as Coulomb collisions come to dominate (`docs/anisotropic-ee-plan.md`).

**SoPLASMA's table-generation pipeline never uses it.** Grep for
`coulomb`/`ionizationDegree` across `src/models/plasmaModels/plasmaBoltzmann/`
returns nothing. Tables are generated as a function of E/N (LFA) or mean
energy (LMEA) alone -- never as a function of local ionisation degree n_e/N,
which is what the `coulomb` option needs as its second input. **This is an
integration gap: real, well-scoped engineering work (a second table axis,
n_e/N piped from the local cell state into the lookup, `coulomb` turned on
in the generator call), not a research problem, and not a fundamental
modelling impossibility.**

## Why this correction still leaves a real, deeper limit

The initial write-up also over-generalised the "conserves the fluid moments"
argument into "cannot matter for SoPlasma at all". That conflated two
different effects of electron-electron collisions:

1. **Bulk Maxwellianization at a point.** e-e collisions reshape the EEDF at
   a GIVEN LOCAL (E/N, mean energy), which shifts every tabulated transport
   and rate coefficient as a function of a second local variable, n_e/N. This
   is exactly the standard swarm-vs-plasma distinction (what BOLSIG+-class
   solvers, including SoEEDF, already compute given an ionisation degree).
   **Wiring this into SoPLASMA is tractable and would be a genuine
   improvement**, wherever the local ionisation degree is non-negligible.

2. **Spatially non-local trapping/de-trapping**, the SPECIFIC mechanism
   Boeuf's book names for the negative glow: a cold electron confined in a
   real-space potential well, escaping only via a rare or cumulative
   Coulomb kick after bouncing for many transit times. This depends on the
   SHAPE of the potential over a finite spatial extent and the electron's
   history in it -- a trajectory/non-local quantity. NO local closure, no
   matter how many local parameters it is given (E/N, mean energy, n_e/N,
   anything else evaluated pointwise), can represent this, because a local
   closure is local by construction: it only ever sees the instantaneous
   state at one point.

So fixing (1) is real, useful, and currently missing -- but it does NOT fix
(2), the actual negative-glow problem. Any SoPLASMA case with a trapped
cold-electron population may still plateau at the wrong bulk Te/density even
after (1) is wired up, and this is expected, not a sign the wiring was done
wrong.

## What would actually address (2), short of full PIC

**Non-local electron kinetic theory** (Tsendin's original formulation;
developed further by Kolobov, Godyak, and others) is the established
middle ground. It solves a reduced kinetic equation in TOTAL ENERGY
(kinetic + local potential energy, eps = m v^2/2 - e*phi(x)), coupled
self-consistently to the actual spatial potential profile across the WHOLE
relevant domain -- not per-cell. A trapped population appears naturally as
electrons whose total energy places them below the local potential maximum,
confined between their turning points; collisions (elastic, inelastic, and
e-e) act as energy-space diffusion/friction operators that let them escape.
This captures trapping at a cost far below full velocity-space PIC-MCC,
because it is a reduced (1D-in-energy, coupled to 1D-in-space) kinetic
equation, not a 6D phase-space simulation.

The concept of non-locality itself (electron energy relaxation length
comparable to or larger than the system size) is discussed extensively in
the JC-PIC book, citing Godyak's experimental and theoretical work on the
local-to-non-local transition in RF discharges (search "non-local" in
`Literature/reference-codes-and-manuals/JC-PIC_Boeuf/library_article_book.txt`) -- but the book treats
it via full PIC (which captures non-locality automatically, being kinetic),
not via the reduced Tsendin-style energy-space equation, since that
technique belongs to a different, older line of low-temperature-plasma
modelling literature rather than to a PIC code's own manual.

## Run-time EEDF solving -- the user's proposed alternative, evaluated

Solving the EEDF "at run time" inside SoPlasma is NOT automatically wrong,
but the answer depends entirely on which kind:

  * A LOCAL run-time solve -- re-running SoEEDF's two-term solver every
    timestep, per cell, using only that cell's own local state instead of
    precomputed tables -- is mathematically equivalent to a local closure
    with infinitely dense, always-fresh tables. It inherits the SAME
    limitation: still cannot represent trapping, no matter how often it is
    refreshed, because it is still local.
  * A NON-LOCAL run-time solve -- in total-energy space, coupled
    self-consistently to the real, evolving potential profile across the
    domain, updated as that profile changes -- is the Tsendin/Kolobov-Godyak
    approach, and CAN in principle capture trapping. This would be a
    substantial, specialised development effort inside SoPlasma: more than
    the fluid closure it has today, but well short of a full PIC or hybrid
    (PIC/fluid) capability.

## Consequence for validation work (unchanged from the original finding)

Any SoPLASMA case resembling a negative glow (Grubert 2009, the planned
Carlsson/JC-PIC GD benchmark reproduction -- `docs/design/verification-map.md` item
L) may plateau at a Te/density that is wrong regardless of run duration,
floor settings, or circuit design, UNLESS AND UNTIL either (1) is wired up
(improves bulk coefficients but not trapping) or a non-local kinetic
treatment is added (addresses trapping directly). If the cathode fall
proper reproduces kinetic benchmarks well while the negative-glow bulk does
not, that is the expected signature of exactly this gap, localised where
theory says it should be (the cathode fall is beam-like/non-equilibrium,
where e-e collisions are not the controlling physics per the same source).

## Related, longer-horizon idea (recorded, not scheduled)

The user separately noted a longer-term goal: a PIC and/or hybrid
(PIC/fluid) simulation capability alongside SoPLASMA's fluid solver --
tentatively "SoKinetic" or "SoHybrid" -- explicitly NOT for now. Non-local
kinetic theory (above) is a smaller, more tractable step in the same
direction, worth considering as an intermediate target if resolving the
negative-glow gap becomes an active priority before a full kinetic/hybrid
capability exists. See the memory `future-pic-hybrid-capability`.
