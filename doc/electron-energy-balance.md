# The electron energy balance at high E/N, and what the mean-energy clamp is for

Measured 2026-09-07, prompted by the user asking how the whole domain could sit
at hundreds of eV when the plasma bulk must be below ~10 eV.

## The question

`grubert2009_n11_clampDerived` (mean-energy clamp DERIVED at the table top,
2644 eV) put the **entire domain** above 570 eV and climbing, against its
one-variable control `grubert2009_n11_clamp100` (clamp pinned at the historical
100 eV) which held a sane bulk at **6.95 eV**. Domain-integrated Joule/loss
ratio: **4.35** in the clamped arm, **577 and rising** in the derived arm.

## FIRST ANSWER, AND IT WAS WRONG

The tables all come from one Boltzmann solution, so at the LFA equilibrium
`U_eq(E/N)` the Joule input should equal the tabulated power loss. Measured
`Joule/Loss = (muN/N)E^2 / ((Pel+Pin)N)`:

| E/N [Td] | U_eq [eV] | Joule/Loss |
|---|---|---|
| 100 | 6.80 | 1.055 |
| 2450 | 34.84 | 2.76 |
| 13000 | 353.09 | 21.6 |
| 55000 | 2644.46 | **155** |

`PelecN`, `PgasN` and `PvibN` are IDENTICALLY ZERO for atomic argon, so
`Pel + Pin` is the complete collisional loss -- no missing channel. This looked
like a defect in the sweep's power-loss moments.

## THE DISCRIMINATOR: it is not a defect. The identity was incomplete.

A self-similar GROWING swarm holds `U` constant only if the field also pays to
endow each new electron with the mean energy. The balance is

    Joule  =  Loss  +  nu_i * U ,        nu_i = N * k_iz(U)

| E/N [Td] | U_eq [eV] | Loss | nu_i*U | loss-only ratio | **WITH GROWTH** |
|---|---|---|---|---|---|
| 100 | 6.80 | 1.908e8 | 1.055e7 | 1.055 | **1.0001** |
| 500 | 10.23 | 3.079e9 | 9.917e8 | 1.324 | **1.0014** |
| 1000 | 15.38 | 9.08e9 | 5.725e9 | 1.636 | **1.0032** |
| 2450 | 34.84 | 2.749e10 | 4.766e10 | 2.760 | **1.0094** |
| 5000 | 80.02 | 4.935e10 | 2.131e11 | 5.407 | **1.0168** |
| 13000 | 353.09 | 7.395e10 | 1.49e12 | 21.57 | **1.0199** |
| 25000 | 1038.72 | 8.515e10 | 5.14e12 | 62.09 | **1.0118** |
| 55000 | 2644.46 | 1.127e11 | 1.742e13 | 155.1 | **0.9968** |

**0.997 to 1.020 over five decades of E/N.** The tables are self-consistent and
the sweep is fine. Recorded because the loss-only table above is exactly the
shape of a real bug (see `docs/models/energy/lmea.md`, trap 1, where the loss
term really was 2.4e25 too small) and would have been reported as one.

## THE PHYSICS THIS EXPOSES, which is the actual result

The GROWTH term dominates the high-field budget:

| E/N [Td] | growth share of the power budget |
|---|---|
| 100 | 5% |
| 2450 | 63% |
| 13000 | 95% |
| 55000 | **99.4%** |

So above ~2450 Td essentially all of the field's power goes into CREATING
ELECTRONS, not into collisional dissipation. Two consequences:

1. **`meanEnergy_vs_reducedE` is the equilibrium of a FREELY GROWING swarm.** It
   is only the right closure where the electron population really is
   multiplying at the full `nu_i`. A steady, NON-multiplying plasma at high E/N
   has only `Loss` to balance the field, and `Loss` is 155x too small at
   55,000 Td -- so **no such state exists**. That is real physics, not a code
   defect: a 55,000 Td field in argon cannot support a steady electron
   population, it must avalanche.

2. **The fluid model gets the growth term for free ONLY IF `n_e` and `nEps` are
   advanced consistently.** `U = nEps/n_e`, so
   `dU/dt = (Joule - Loss) - U*nu_i` arises from the `n_e` growth in the
   DENOMINATOR, not from any term written in the energy equation. If the outer
   loop does not converge the `n_e`/`nEps` pair tightly, that cancellation is
   inexact exactly where it carries 99% of the budget.

   **This connects to the deferred outer-coupling item**, which measured Aitken
   damping `omega = 0.148` and named `nEps_e` as the worst field. That is the
   same coupling. At high `U` the cancellation is between two large terms, so
   loose convergence there is not an efficiency question -- it is an accuracy
   question about the dominant balance.

## WHAT THIS MEANS FOR THE CLAMP -- my change was wrong in EFFECT

`c131d8c` made `meanEnergyMax` derive from the table range (2644 eV) on the
argument that the clamp exists to prevent EXTRAPOLATION, so its value should be
the tabulated range. That argument is still correct as far as it goes, and the
old hardcoded 100 eV was still an unsourced constant.

**But a table's EXTENT is not its VALIDITY for this purpose.** The tables are
valid to 2644 eV for a growing swarm; the fluid model's energy equation is only
robustly DISSIPATIVE while `Loss` can balance the field, i.e. below roughly
1000-2450 Td (`U` ~ 15-35 eV). Between there and the table top the equation
depends on an exact cancellation against `n_e` growth.

So the 100 eV clamp was **load-bearing, not cosmetic**: it held the solution in
the regime where the energy equation is dissipative. Removing it did not reveal
a hidden bug in the tables -- it removed a guard against a regime the fluid
closure cannot robustly represent.

MEASURED CONSEQUENCE, the two arms at the same physical time:

| | clamp100 | clampDerived |
|---|---|---|
| meanE range | [6.95, 100.0] eV | **[622.97, 2644.40] eV** |
| Joule/loss, domain | 4.35 | **577** |
| local E/N at the max | 18,261 Td | **123,073 Td** (2.2x past the table) |

The likely mechanism for the spread across the whole domain: `DEpsN` grows with
`U` while `Loss` SATURATES (Pin grows 4.1x while `U` grows 76x), so a hot region
diffuses energy outward faster than it can dissipate it. The clamp broke that
feedback by capping `U`.

## NOT SETTLED, and it decides what the clamp should be

Whether the right fix is (i) a clamp tied to the DISSIPATIVE limit rather than
the table extent, (ii) tightening the `n_e`/`nEps` coupling so the growth
cancellation is exact, or (iii) accepting that the cathode fall of this
discharge is outside what a two-term LMEA closure can represent, and saying so.
These need different work and the evidence does not yet choose between them.

**Do NOT simply revert the clamp to 100 eV.** That restores the old behaviour
for a reason nobody understood and re-hides this.
