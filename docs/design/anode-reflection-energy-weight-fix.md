# The anode-reflection defect and its fix

Opened 2026-09-07. Status: mechanism identified, fix drafted, NOT yet applied.

## What is measured, and therefore not in doubt

`validation/grubert2009_unif_probe` (r = 0.36) vs `validation/grubert2009_anode_r0`
(r = 0). Two lines differ: `electronReflection` on the anode in `0/n_e` and
`0/nEps_e`. r = 0.36 runs away at 1.125 us (n_e,max 1.94e18 m^-3, Te,max 1763 eV,
Te,min 817 eV). r = 0 completes 1.3 us with Te,max 11.58 eV, Te,min 2.66 eV.

At t ~ 0.4 us the two runs have the SAME n_e,max (6.53e9 vs 6.64e9, 1.7%) and
Te,max already differs 2.4x (16.43 vs 6.83 eV). The temperature diverges BEFORE
the density does, so the carrier is an energy channel, not a density one.

## The mechanism

`hagelaarClosure()` (electronDDWallFluxMixedFvPatchScalarField.H:347):

    wW = max((A + Dd - Gc/n)/(1+r), Dd);   W = (1-r)*wW;   gRatio = (r*wW + Gc/n)/A

`hagelaarEnergyWeight()` then gives `ew = max(0, 5/3 - (2/3)*gRatio)`, and the
energy condition's loss speed is `W_eps = W*ew`.

`wW` contains the drift and grows proportional to Dd, but `gRatio` divides by
`A`, THE THERMAL SPEED ALONE. So gRatio is unbounded in the field, and the
`max(0, ...)` fires. Unit sweep, Gc/n = 0:

    r = 0.36:  Dd/A = 2     ew 1.137   W_eps/A  1.606   W/A  1.412
               Dd/A = 5     ew 0.467   W_eps/A  1.493   W/A  3.200
               Dd/A = 8.44  ew 0.000   W_eps/A  0.000   W/A  5.402
               Dd/A = 50    ew 0.000   W_eps/A  0.000   W/A 32.000

Past Dd/A = 8.44 THE WALL CARRIES AWAY ZERO ENERGY WHILE THE ELECTRON SINK KEEPS
GROWING. Electrons deliver charge to the anode and deposit no energy, so
everything they took from the field stays in the plasma. Thermal runaway by
construction. W_eps is also NON-MONOTONIC in Dd, which a physical wall sink
cannot be. Cliff at Dd/A = 50 for r = 0.1, 8.44 for r = 0.36, 5.0 for r = 0.5.

## Why it was never seen

gRatio is identically zero unless r > 0 or SEE is on at that wall. The anode has
`enableSEE false`. So for every case this framework ran before, gRatio == 0 and
`hagelaarEnergyWeight` returned a constant 5/3 -- the whole path was dead code.
`configuration/boundaries` records this as framework-state gap 16: "nothing in
this framework had ever exercised a non-zero reflection coefficient before".
The first use of the path is the one that broke.

## The proposed fix -- one normalisation

Divide gRatio by the TOTAL one-way arriving flux `A + Dd`, not by `A`:

    const scalar arriving = A + Dd;
    gRatio = (arriving > VSMALL) ? (r*wW + Gc/n)/arriving : 0.0;

`hagelaarClosure` already takes Dd as a parameter, so the change is contained to
that one static function.

Behaviour: gRatio -> r/(1+r) at Dd = 0 (unchanged), -> r as Dd -> inf (bounded).
ew never reaches 0 for r < 1: 1.428 at r = 0.36 and Dd/A = 200, 1.335 at r = 0.5,
1.070 at r = 0.9. W_eps monotonically increasing in Dd for every r.

WHY IT CANNOT REGRESS ANYTHING ALREADY VALIDATED:
  * At Dd = 0 it is IDENTICAL to the current expression, for every r. That is
    the drift-free limit eq. (6.6)/(6.15) was derived for.
  * At r = 0 with no SEE both give exactly 5/3.
  * At the cathode, electrons are repelled, so `Dd = max(0, uDriftNormal) = 0`
    and it is identical even with gammaSEE 0.06 on.
It changes behaviour ONLY where electrons drift INTO a wall that also reflects
or emits -- exactly the untested combination.

## NOT SOURCED YET (G2)

Eq. (6.15) and the definition of gRatio have NOT been checked against Hagelaar's
text. The case for `A + Dd` above is a physical argument -- boundedness,
monotonicity, correct drift-free limit -- not a citation. THIS MUST BE CHECKED
BEFORE THE FIX IS COMMITTED. It is possible the intended quantity is something
else again and `A` is a transcription of a symbol that means the total flux.

## The ion condition is NOT a defect -- closing an earlier suspicion

`ionDDWallFluxMixedFvPatchScalarField.C:51`: `uEff = A_i + max(0, -uDrift_n)`.
This was raised as a suspect on 2026-09-07 and does not hold up:
  * drift away from wall: total flux = n_p(uDrift + A_i - uDrift) = n_p*A_i.
  * drift into wall:      total flux = n_p(uDrift + A_i).
Both correct, and `uEff > 0` always, so `f = uEff/(D/delta + uEff)` is in (0,1)
and n_p is bounded. The single-cell Ar+ hole seen at r = 0.36 was a CONSEQUENCE
of the E spike, not a cause -- at r = 0 the ion profile is smooth.

## Plan

  A. CONFIRM (running 2026-09-07): split r between the two conditions.
     `grubert2009_r_eps_only` (r on nEps_e only) -- PREDICTED RUNAWAY.
     `grubert2009_r_ne_only`  (r on n_e only)    -- PREDICTED BENIGN.
     If r_ne_only also runs away, this mechanism is WRONG and so is the fix.
  B. SOURCE the energy weight against Hagelaar ch. 6. Blocking for commit.
  C. APPLY the one-line change in hagelaarClosure().
  D. GUARD: `ew == 0` while `W > 0` is an unphysical state (electron sink on,
     energy sink off) and must warn rather than clamp silently. Follow the
     existing `nearSingularReported_` pattern in ddWallFluxMixed.
  E. REGRESS: r = 0 cases must be bit-identical. Then re-run
     `grubert2009_unif_probe` at r = 0.36 and require it to track
     `grubert2009_anode_r0`.
  F. OPEN, SEPARATE: the anode wall-cell Poisson residual is nonzero even at
     r = 0 (sum(grad phi . n) = 6.4 against a required 0.029). 21000x smaller
     than at r = 0.36, but a converged FV Poisson solve should satisfy it.
