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

`8f5bab4` made `meanEnergyMax` derive from the table range (2644 eV) on the
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

---

# IS THE meanE SOLUTION WRONG? No -- the FIELD is. (2026-09-07)

The user judged that even 60 eV in the bulk is far too high for this discharge,
and then proposed that the field may not be properly screened in the bulk. Both
were tested. The second is correct.

## The meanE solution is CONSISTENT with its own local field

`meanE` from the solver against `U_eq(E/N_local)` from `meanEnergy_vs_reducedE`,
at the same cells:

| x/L | EARLY (healthy, -98 V) |  | FROZEN (+209 V) |  |
|---|---|---|---|---|
| | E/N [Td] | solver/U_eq | E/N [Td] | solver/U_eq |
| 0.05 | 702 | 0.671 | 4227 | 1.184 |
| 0.15 | 669 | 0.903 | 5839 | 0.877 |
| 0.30 | 596 | 0.925 | 198 | 1.921 |
| 0.50 | 447 | 0.931 | 198 | 0.944 |
| 0.72 | 223 | 1.000 | 1.8 | 1.683 |
| 0.90 | 57 | 1.120 | 162 | 0.853 |

Ratios cluster on 1. The one large excursion is **8.9x at x/L = 0.020 in the
frozen state**, which is the field-reversal / potential-maximum cell: 81 eV
where the local field alone would give 9 eV. That is energy TRANSPORTED in from
the adjacent high-field region -- exactly the effect LMEA exists to capture, and
the reason `meanE` is a transported field rather than a lookup.

**IN THE EARLY, HEALTHY STATE meanE IS 3.5-10.7 eV (Te 2.3-7.1 eV)** -- i.e.
precisely the range the user expects. The 60-100 eV values occur ONLY in the
already-diverged state, where the local field is 4200-5800 Td. At 5000 Td,
`U_eq = 80 eV`. So the energy equation is returning the right answer to the
field it is given. **The defect is upstream of the energy equation.**

## THE FIELD IS NOT SCREENED, and it CANNOT BE at these densities

Screening requires `lambda_D << ` the bulk scale. Measured:

| state | n_e | gap/lambda_D |
|---|---|---|
| **Grubert's own glow** | 2.478e15, Te~3 eV | **39** |
| early, x/L 0.05 | 1.25e11 | **0.20** |
| early, x/L 0.30 | 7.49e11 | 0.44 |
| early, x/L 0.90 | 3.45e13 | 3.64 |
| frozen, x/L 0.05 | 1.74e18 | 248 |

**In the early state `lambda_D` is 0.27 to 4.9 cm against a 1 cm gap.** The
plasma is 2-4 DECADES too thin to screen anything, so the field is essentially
the applied vacuum field: `E/N_actual / E/N_vacuum` = 1.73, 1.65, 1.47, 1.10,
0.55, 0.14 across the gap, against a 406 Td vacuum value.

### Why that is the whole problem

At **400-700 Td EVERYWHERE**, ionisation happens throughout the entire gap. A
normal glow has 12-30 Td in the bulk and thousands of Td in a THIN cathode fall;
this gap has ~500 Td uniformly. So the discharge never builds a cathode-fall
structure -- it avalanches as a bulk Townsend discharge across the whole gap,
which is the correct behaviour for an unscreened gap and the wrong behaviour for
the glow we are trying to reach.

By the time it CAN screen (frozen state, gap/lambda_D = 10-248) it screens into
the wrong structure: the field piles up at x/L = 0.05-0.15 rather than at the
cathode, and the electrode has already flipped positive.

Ohmic check at the frozen state, using THIS case's own `j = 509 A/m^2`:
`E_actual/E_ohmic` = 1031 at x/L 0.05, 105 at 0.15, then 0.32, 0.12, 0.00, 0.01
beyond x/L 0.3. So the near-electrode field is set by SPACE CHARGE (1000x more
than conduction needs) while the bulk is over-dense and needs almost no field to
carry the current. That is a sheath in the wrong place, not a glow.

CORRECTION recorded: the first pass at this Ohmic comparison used
`I_cond = 0.0269 A` from the DIVERGED `grubert2009_steady` run instead of the
spike's own `1.018e-4 A` -- wrong by 264x. The table above is the corrected one.
Rule 13: never mix cases.

## What this reframes

The chain is `unscreened gap -> uniform ~500 Td -> whole-gap avalanche ->
over-dense plasma -> abrupt screening into a mislocated sheath -> electrode
flip`. The energy equation, the mesh and the tables are all doing their job.

The open question is therefore NOT "why is meanE high" but **"why does this case
start from a state that cannot screen, and is there a path to the glow that does
not pass through a whole-gap avalanche?"** In a real experiment the answer is
that the glow is struck and then the current is raised slowly, so the fall forms
first. That is a statement about the PATH, and it is the same conclusion the
`C_gap/g` bound reached from the circuit side.
