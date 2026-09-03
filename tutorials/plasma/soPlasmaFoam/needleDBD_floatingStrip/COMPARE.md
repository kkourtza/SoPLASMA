# floatingElectrode in a PLASMA -- the charge ledger and Gauss closure

## Question
With a plasma present, does the floating electrode
  (a) integrate its charge ledger Q(t) = Q0 + INT I_plasma dt' correctly, and
  (b) still close Gauss's law exactly, given that psi must be rebuilt every
      step from the SEMI-IMPLICIT operator eps + dt*sigma rather than eps?

(b) is the point. Superposition V += dV_f*psi is exact ONLY if psi solves the
same operator as the field it corrects. semiImplicit is the DEFAULT scheme, so
if that were got wrong every plasma run would be quietly wrong.

## Baseline
DERIVED FROM, and must be compared against:
    /home/kkourtza/soplasma-scratch/tutorials/plasma/soPlasmaFoam/needleDBD
This case is that one with `air_dielectric` (178 faces) REINTERPRETED as an
isolated metal strip -- a floating conductor -- instead of a free-standing thin
dielectric. That patch is used because it already carries real species
wall-flux boundary conditions, so there is genuine charge arriving.

Everything else is identical: same mesh, same Boltzmann tables, same chemistry,
same drive. `./Allrun-derive` copies them rather than regenerating, so the two
cases cannot drift apart.

## Reference numbers, measured 2026-09-03 (endTime 1e-10, 18 steps)

    C_self                  7.03018250107e-16 F
    psi rebuild cadence     every step (semiImplicit operator)
    max |closure|           2.82e-14 V          <- Gauss's law, machine level
    ledger self-consistency 7.5e-12 relative
        worst |Q_k - (Q_k-1 + I_k dt_k)| / |Q_k| over all steps
    I_plasma at t_end       -2.41e-12 A         <- NEGATIVE, electrons arriving
    Q at t_end              -2.373472e-22 C     <- NEGATIVE, charging negative
    V_f at t_end             7.24693 V

HONEST LIMITS OF THIS CASE:
  * At 1e-10 s V_f is INDUCTION-dominated: the charge contributes Q/C_self
    = -3.4e-07 V out of 7.25 V, i.e. 2e7 times smaller. So this validates the
    LEDGER and the CLOSURE, NOT the settling to a floating potential. The
    probe-theory estimate V_f - V_plasma ~ -(kTe/2e) ln(2 pi me/mi) needs a far
    longer run.
  * Only n_e has a wall-flux condition on that patch in needleDBD; nEps_e and
    every ion are zeroGradient there (see the note in the parent case). So
    I_plasma here is ELECTRON-ONLY. That makes the negative sign expected
    rather than surprising, but it is still a real check of the sign
    convention.

## Extraction
    ./Allrun-derive      # derives, runs, and checks against the numbers above
