# Why the Grubert case cannot be reached transiently, measured 2026-09-06

## Summary

The runaway is **real physics**, quantitatively predicted by the framework's
own coefficients, and **no series-resistor circuit can arrest it**. Grubert et
al. solved the TIME-INDEPENDENT equations, and this analysis says that is not
a convenience -- it is what the problem requires.

Nothing here indicts the framework: gamma, alpha, k_iz and the space-charge
screening all check out independently.

## 1. The instability rate is what our own coefficients say

At the measured cathode-fall field:

| at E/N = 2610 Td | value |
|---|---|
| equilibrium mean energy (`meanEnergy_vs_reducedE`) | 37.16 eV |
| `k_iz` at that energy (`k_EI_AR_ION_AR_vs_meanE`) | 6.075e-14 m^3/s |
| `nu_iz = k_iz * N`, N = 2.414e22 | 1.467e9 1/s |
| **predicted e-folding** | **0.682 ns** |
| **measured e-folding** (fix_n11 runaway) | **0.77 ns** |

13% agreement. The growth is not numerical and not an over-productive model.

## 2. gamma and alpha corroborate each other independently

Self-sustaining (Townsend) condition: `alpha*d = ln(1 + 1/gamma)`.

    gamma = 0.06 (boundary condition)  ->  ln(1+1/gamma) = 2.872
    alpha table (Boltzmann solver)     ->  alpha*d = 2.87 at exactly 500 Td

Two unrelated sources, agreeing to three figures. Neither is suspect.

    E/N = 500 Td  -> alpha*d =  2.87  = 1.00x  <- breakdown, by construction
    E/N = 1056 Td -> alpha*d =  6.34  = 2.21x  <- the gap at V = -190 V
    E/N = 2000 Td -> alpha*d =  9.77  = 3.40x
    E/N = 2610 Td -> alpha*d = 10.89  = 3.79x  <- the fall, once formed

## 3. A series resistor CANNOT limit it -- an R-independent bound

    tau_RC = R * C_gap
    I_sc/I_op - 1 = V_gap/(R * I_op)
    =>  tau_RC * (I_sc/I_op - 1) = C_gap * V_gap / I_op = 86.6 ns    [R cancels]

With `C_gap = 1.77e-16 F` (the diagnostic's own value), `V_gap = 500 V`,
`I_op = 1.022e-6 A`. Responding on 0.68 ns therefore demands
`I_sc/I_op ~ 113`, i.e. no current limit at all. Measured, per e-folding:

| config | tau_RC | ballast moves | I_sc/I_op |
|---|---|---|---|
| R=1e8, C=5e-14 (`fix_n11`) | 5018 ns | 0.02% | 5.9 |
| R=1e9, C=5e-14 | 50177 ns | 0.00% | 1.5 |
| R=1e8, C=0 | 17.7 ns | 4.3% | 5.9 |
| R=1e7, C=0 | 1.77 ns | 35% | 49.9 |

**CORRECTION to an earlier statement in this work:** the ballast's adequacy was
first argued by comparing `tau_RC` to the TIMESTEP (`a = tau/dt`). That is the
wrong comparison -- the circuit does not need to act within one step, it needs
to act within one e-folding of the instability. The verdict happens to be the
same here, but the reasoning was wrong and the corrected form is above.

## 4. Why a real glow is stable, and ours is not

In steady state the self-regulation is an INTEGRAL condition over the actual
profile:

    INT alpha(E(z)) dz = ln(1 + 1/gamma)

A real glow satisfies it by dropping nearly all of its voltage across a thin
cathode fall and leaving the bulk field-free, so the integral sits at 2.87
however large the applied voltage is. Our transient has the field spread over
the gap -- measured `|E|` flat from the cathode to z/d ~ 0.35 at t = 2.0 us --
so the integral is 2-4x too large and the discharge grows.

A second mechanism is unavailable in 1-D: a real normal glow holds its CURRENT
DENSITY at j_n and adjusts the AREA it covers on the cathode. A 1-D model
cannot vary area, so it must vary j -- and the normal-glow branch has
`dV/dj < 0`, which is why the current must be set externally.

## 5. What this means for the roadmap

* **The steady-state solver is now required, not preferred.** It reaches the
  regulated state without traversing an instability that no lumped circuit can
  hold. The user flagged it as a wanted capability on 2026-09-06; this analysis
  is the physical argument for it.
* **The external circuit remains correct and necessary** -- it is what selects
  the operating point once a stable solution exists. `fix_fast` demonstrated
  that a fixed-voltage gap has no current limit at all.
* **`currentSource` is not a shortcut past this.** Its regulator responds as
  `dV = dt(I_set - I)/(C + |g| dt)`, so with `C = C_gap` and the measured
  `g ~ 6e-7 S` the capacitive term dominates by orders of magnitude and it is
  as slow as the ballast. It is the right instrument for a STEADY or slowly
  swept case, which is what it was argued for.
* **Not indicted:** gamma (corroborated), alpha and k_iz (they predict the
  measured rate), the screening (tau_diel = 0.018 ns, 43x faster than
  ionisation, and `Emag` min ~ 0 in the bulk), the wall-flux closure, and the
  charge-density fix of `23cd599`, without which none of these numbers would
  have been meaningful.

## Provenance

All numbers measured 2026-09-06 from `validation/grubert2009_fix_n11` and the
case's own generated tables in `constant/plasmaTables` and
`constant/ionTables`. Two unit errors were made and corrected while producing
this: the alpha table is REDUCED (`alpha/N`, m^2) and was first read as alpha,
and `DLN_Arp` (longitudinal diffusion) was first read as a mobility, giving a
259 us ion transit instead of the correct 3.7-6.7 us.
