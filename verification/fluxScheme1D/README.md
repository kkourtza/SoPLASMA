# fluxScheme1D -- interior flux-scheme verification

`./Allrun` sweeps mesh (N = 20..640) x domain Peclet (0.1, 1, 10, 100) x scheme
{standard+upwind, standard+linear, standard+ROUNDF, ScharfetterGummel} against
an EXACT solution and fits the convergence order. `results.baseline.txt` is the
2026-09-07 reference.

## The problem

    d/dx ( v n - D dn/dx ) = S,   n(0) = a,  n(1) = b
    n(x) = a + K (exp(Pe(x-1)) - exp(-Pe))/(1 - exp(-Pe)) + (S/v) x
    Pe = v L / D,   K = b - a - S L / v

THE SOURCE IS NOT A DECORATION. Scharfetter-Gummel derives its face flux from
this same two-point problem WITH s SET TO ZERO. The Complete Flux Scheme (Liu
et al 2014 PSST 23 015023) keeps s. A source-free bed cannot separate them, so
it would be useless as the CFS acceptance test.

## Results, 2026-09-07 (asymptotic order at N = 640)

    scheme              Pe=0.1   Pe=1    Pe=10   Pe=100
    std:upwind           0.93    2.00     0.98    0.75
    std:linear           2.00   (exact)   2.00    2.00
    std:ROUNDF           1.14    2.00     1.02    1.15
    ScharfetterGummel    2.00    2.00     2.00    2.00

1. ROUNDF -- THE PRODUCTION SCHEME in every SoPLASMA case -- converges at
   ~1.0-1.15, NOT second order, at Pe = 0.1, 10 and 100. Expected of a bounded
   limiter (limiters drop to first order near extrema), but now measured. At
   Pe = 100, N = 640 its error is 1.45e-3 against SG's 3.05e-7: FOUR ORDERS.
   Boundedness is worth paying for on a density that must stay positive; this
   says what it costs.

2. SG is uniformly second order at every Peclet tested, 0.1 to 100. The
   implementation is verified.

3. SG's ERROR IS A UNIFORM OFFSET -- L2 == Linf to all printed digits in every
   SG row -- and at N = 20 it is 3.125e-4 = h^2/8 EXACTLY. This is the source
   term it omits from the flux. The Pe = 1 column proves it: there
   K = b - a - S L/v = 0, so the exact solution is the straight line n = 1 + x;
   `linear` reproduces it to 1e-15 while SG carries the full h^2/8. That
   degeneracy is deliberate -- it isolates the source contribution.

## SG IS FIRST ORDER AT HIGH GRID PECLET -- and the reason is the source

The table above measures the ASYMPTOTIC order, and h -> 0 drives PeGrid -> 0 by
construction, so it only ever samples the second-order regime. Liu's claim is
about GRID Peclet. Sweeping domain Pe = 1000 and 10000 keeps |PeGrid| >> 1
across the whole refinement and shows it plainly:

    domain Pe = 10000, ScharfetterGummel
    N      PeGrid      L2            order
    20     500.0       2.50000e-06
    40     250.0       1.25000e-06   1.00
    80     125.0       6.25000e-07   1.00
    160     62.5       3.12500e-07   1.00
    320     31.2       1.56250e-07   1.00
    640     15.6       7.80618e-08   1.00

Exactly 1.00 over five refinements. At domain Pe = 1000 the order transitions
1.00 -> 1.01 -> 1.12 -> 1.49 -> 1.81 -> 1.95 as PeGrid falls through 1.

THE MECHANISM. The SG error stays a UNIFORM offset (L2 == Linf everywhere) and
equals the source flux SG omits over half a cell:

    low  PeGrid:   S h^2 / (8 D)      -> O(h^2), second order
    high PeGrid:   S h   / (2 v)      -> O(h),   FIRST order

Checked: Pe=1e4, v=1e4, N=20  -> S h/(2v) = 2.50e-6, measured 2.50000e-06.
         Pe=1e3, v=1e3, N=20  ->            2.50e-5, measured 2.50000e-05.

So SG's first-order behaviour at high grid Peclet IS the omitted source term.
Note SG is nevertheless FAR more accurate in absolute terms than ROUNDF there
(2.5e-6 vs 64.9 at Pe=1e4, N=20). Small error and first-order decay are both
true and not in tension.

## What CFS has to do

TWO acceptance criteria, both from the omitted source term:

  1. Drive the SG rows to machine precision at Pe = 1, where the exact solution
     is linear and purely source-driven. Failing this means the inhomogeneous
     flux is not implemented at all.
  2. RESTORE SECOND ORDER at domain Pe = 1000 and 10000, where SG is measured
     at exactly 1.00. This is Liu's "second order uniformly in Peclet" claim and
     is the sharper of the two tests.
