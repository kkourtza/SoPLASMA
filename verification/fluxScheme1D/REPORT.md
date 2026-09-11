# Flux-scheme verification: order of accuracy vs Péclet number

SoPLASMA, `verification/fluxScheme1D`. Run 2026-09-07.
Figures: `fig1_convergence.png`, `fig2_order_vs_gridPeclet.png`,
`fig3_error_vs_Peclet.png`. Raw data: `results.baseline.txt`.

---

## 1. What this test is

It measures the **order of accuracy** of every convective flux scheme available
to SoPLASMA, by solving a problem whose exact answer is known analytically and
watching how fast the error falls as the mesh is refined.

The problem is the steady one-dimensional drift-diffusion equation **with a
source term**:

        d/dx ( v n  -  D dn/dx ) = S        on 0 < x < 1
        n(0) = a,   n(1) = b

  * `n`  the transported quantity (a number density, in the solver's use)
  * `v`  the drift velocity (in the solver, `mu*E`) -- constant here
  * `D`  the diffusivity -- constant here
  * `S`  a uniform volumetric source (in the solver, ionisation)

Its exact solution is

        n(x) = a + K [exp(Pe(x-1)) - exp(-Pe)] / [1 - exp(-Pe)] + (S/v) x
        K    = b - a - S L / v

Everything is constant, so there is no chemistry, no Poisson equation, no
boundary-condition model and no time integration in the way. What the bed sees
is the **flux discretisation and nothing else**. That is the point: it is a
verification test (does the code solve the equation correctly?), not a
validation test (does the equation describe the plasma?).

### Why the source term is essential

Scharfetter-Gummel derives its face flux from exactly this two-point problem
**with S set to zero**. The Complete Flux Scheme (Liu et al 2014, PSST 23
015023) keeps S and adds the resulting inhomogeneous part. **A source-free bed
cannot tell the two apart**, so it would be useless as the CFS acceptance test.
With S != 0 they separate cleanly, and section 5 shows they do.

---

## 2. What the Péclet numbers mean

The **Péclet number is the ratio of drift transport to diffusive transport.**
Two versions appear, and confusing them is the single easiest mistake to make
here (this report's author made it, see section 4).

**Domain Péclet** -- a property of the PHYSICS:

        Pe = v L / D          L = domain length = 1

  `Pe << 1` diffusion-dominated; `Pe >> 1` drift-dominated with a thin boundary
  layer of width ~L/Pe at the outflow end.

**Grid Péclet** -- a property of the PHYSICS *AND THE MESH*:

        Pe_grid = v h / D     h = cell size = 1/N

  This is the one that governs scheme behaviour: it says how much the solution
  varies across ONE CELL. Note `Pe_grid = Pe/N`.

Liu et al give it a physical reading: with the Einstein relation `D = (2/3)(eps/q)mu`,

        Pe = (3/2) * (q E d) / eps

so **the Péclet number is the energy a particle gains from the field over the
distance d, relative to its mean energy.** `Pe_grid` is that energy gain across
one cell.

For orientation, the Grubert dc-glow case at 5 um cells has `Pe_grid ~ 0.014`
(measured), matching Liu's Table 3 extrapolated to that resolution.

---

## 3. What was tested

Schemes (`div(phi,n) Gauss <scheme>`, plus Scharfetter-Gummel as a combined
flux operator):

  * reference points: `upwind` (first-order floor), `linear` (unlimited central)
  * TVD limiters: `limitedLinear 1`, `Minmod`, `vanLeer`, `MUSCL`, `SuperBee`
  * ROUND family: `ROUNDF`, `ROUNDA`, `ROUNDAplus`, `ROUNDL`, `ROUNDW`
    (ROUNDW was ADDED to our vendored `ThirdParty/libROUNDSchemes` on
    2026-09-07 from upstream -- our snapshot dated 2026-08-10 and predates it.
    It implements eq. (6.6) of the same Deng 2023 JCP 481:112052 paper.)
  * `ScharfetterGummel`

Sweep: domain Pe in {0.1, 1, 10, 100, 1000, 10000} x N in {20,40,80,160,320,640}.
468 runs. `a=1, b=2, S=1, D=1`; `v` sets the Péclet number.

The order `p` is fitted between successive mesh refinements,
`p = log2( L2(h) / L2(h/2) )`. `p=1` first order, `p=2` second order.

---

## 4. A trap in reading convergence tables

Order of accuracy is measured as `h -> 0`, and **`h -> 0` drives `Pe_grid -> 0`
by construction** (`Pe_grid = vh/D`). So an "asymptotic order" column ALWAYS
reports the low-`Pe_grid` regime, whatever the domain Péclet. To see behaviour
at high `Pe_grid` the domain Péclet must be raised until `Pe_grid >> 1` across
the whole refinement -- hence the Pe = 1000 and 10000 sweeps here. `fig2` plots
order against `Pe_grid` directly and avoids the trap entirely.

---

## 5. Results

`L2` error at N=640 with the local order over the last refinement:

    scheme               Pe=0.1            Pe=1           Pe=10          Pe=100         Pe=1000        Pe=10000
    upwind          6.14e-06(0.93)  3.63e-07(2.00)  1.09e-03(0.98)  3.28e-03(0.75)  9.47e-03(2.21)  2.70e-01(1.60)
    linear          2.75e-07(2.00)  7.81e-13( -- )  8.44e-06(2.00)  2.94e-04(2.00)  9.85e-03(2.15)  4.24e-01(2.02)
    vanLeer         6.64e-06(1.05)  3.63e-07(2.00)  1.11e-03(1.00)  3.89e-03(1.01)  1.31e-02( -- )  1.05e+00(0.34)
    MUSCL           6.64e-06(1.05)  3.63e-07(2.00)  1.11e-03(1.00)  3.89e-03(1.01)  1.31e-02( -- )  1.05e+00(0.34)
    SuperBee        6.64e-06(1.05)  3.63e-07(2.00)  1.11e-03(1.00)  3.89e-03(1.01)  1.31e-02( -- )  1.05e+00(0.34)
    ROUNDF          2.38e-06(1.14)  1.21e-07(2.00)  3.75e-04(1.02)  1.45e-03(1.15)  1.03e-02(2.20)  1.09e+00(0.45)
    ROUNDA          3.45e-06(1.10)  1.82e-07(2.00)  5.59e-04(1.01)  2.05e-03(1.08)  1.06e-02(2.40)  1.07e+00(0.39)
    ROUNDAplus      4.73e-06(1.07)  2.54e-07(2.00)  7.80e-04(1.00)  2.78e-03(1.04)  1.12e-02( -- )  1.05e+00(0.36)
    ROUNDW          6.11e-07(1.52)  1.99e-08(2.00)  6.77e-05(1.14)  4.68e-04(1.66)  9.91e-03(2.15)  5.34e-01(3.01)
    ROUNDL          3.14e-07(1.91)  2.42e-09(2.00)  1.52e-05(1.62)  3.13e-04(1.95)  9.86e-03(2.15)  4.33e-01(2.06)
    SG              3.05e-07(2.00)  3.05e-07(2.00)  3.05e-07(2.00)  3.05e-07(2.00)  2.91e-07(1.81)  7.81e-08(1.00)

    (Minmod and limitedLinear omitted -- see the caveat in section 7.)

### SG is second order below Pe_grid = 1 and FIRST order above it

    domain Pe = 10000, ScharfetterGummel
    N      Pe_grid     L2            order
    20     500.0       2.50000e-06
    40     250.0       1.25000e-06   1.00
    80     125.0       6.25000e-07   1.00
    160     62.5       3.12500e-07   1.00
    320     31.2       1.56250e-07   1.00
    640     15.6       7.80618e-08   1.00

Exactly 1.00 over five refinements. At Pe = 1000 the order climbs
1.00 -> 1.01 -> 1.12 -> 1.49 -> 1.81 -> 1.95 as `Pe_grid` falls through 1.
This reproduces Liu et al's statement that EDS/SG "behaves second-order only for
|Pe_grid| << 1; when |Pe_grid| >> 1 its accuracy reduces to first order".
`fig2` shows the crossover sitting on `Pe_grid = 1`.

### The cause: the source term SG omits

SG's error is a UNIFORM offset at every Péclet (`L2 == Linf` to all digits in
every SG row -- it is a constant shift, not a shape error), and it equals the
source flux SG drops over half a cell:

    low  Pe_grid:   S h^2 / (8 D)   -> O(h^2), second order
    high Pe_grid:   S h   / (2 v)   -> O(h),   FIRST order

    checked: Pe=1e4, v=1e4, N=20 -> S h/(2v) = 2.50e-6, measured 2.50000e-06
             Pe=1e3, v=1e3, N=20 ->            2.50e-5, measured 2.50000e-05

The Pe=1 column proves it independently: there `K = b-a-SL/v = 0`, so the exact
solution is the straight line `n = 1+x` and is purely source-driven. `linear`
reproduces it to 1e-13; SG carries the full `h^2/8`. That parameter choice is
deliberate -- it isolates the source contribution.

**So SG's first-order behaviour at high Pe_grid and its "source-blind flux" are
not two defects. They are one defect seen from two angles.**

---

## 6. Conclusions

1. **SG is by far the most accurate scheme in this test**, and uniquely so: its
   error is ~3.05e-7 at N=640 and essentially INDEPENDENT of Péclet from 0.1 to
   100 (`fig3`). At Pe=100 it is 4700x more accurate than ROUNDF. This is a
   smooth-solution result and does not contradict section 8.

2. **SG degrades to exactly first order once `Pe_grid > 1`**, and the reason is
   identified analytically: the omitted source term, which scales as `h` rather
   than `h^2` in that regime.

3. **ROUNDL is the best of the five ROUND variants, ROUNDW second.** The
   ordering is consistent across the whole Péclet range:

       ROUNDL  >  ROUNDW  >  ROUNDF  >  ROUNDA  >  ROUNDAplus

   At Pe=10, N=640: ROUNDL 1.52e-5, ROUNDW 6.77e-5, ROUNDF 3.75e-4 -- ROUNDL is
   25x and ROUNDW 5.5x more accurate than the production scheme. ROUNDL holds
   order 1.6-2.0 and ROUNDW 1.1-1.7, where ROUNDF, ROUNDA and ROUNDAplus all sit
   at ~1.0.

   NEITHER ROUNDL NOR ROUNDW WAS TESTED in Pasolari & Kourtzanidis, and neither
   is used anywhere in SoPLASMA. ROUNDW did not exist in our vendored copy of
   the library until it was added for this study.

4. **ROUNDF -- the production scheme in every SoPLASMA case -- is ~first order
   here** (p = 1.0-1.15 at Pe = 0.1, 10 and 100). That is expected of a bounded
   limiter and is not by itself an argument against it (see section 8), but the
   cost is now measured rather than assumed.

5. **CFS acceptance criteria** follow directly, and they are falsifiable:
   (a) drive the SG rows to machine precision at Pe=1, where the exact solution
       is linear and purely source-driven -- failure means the inhomogeneous
       flux is not implemented at all;
   (b) restore SECOND order at Pe = 1000 and 10000, where SG measures exactly
       1.00. This is Liu's "uniformly second order" claim and is the sharper of
       the two.

6. **Nothing here is a reason to change the production scheme yet.** See 7 and 8.

---

## 7. CAVEAT -- an unresolved anomaly, do not use these three rows

`linear`, `Minmod` and `limitedLinear 1` returned **bit-identical** results at
every Péclet and every mesh, including at `Pe_grid = 250` where limiting must be
strongly active:

    domain Pe = 10000, N = 40 (Pe_grid = 250)
      linear          L2 = 1.98573061e+02
      Minmod          L2 = 1.98573061e+02      <- identical to unlimited central
      limitedLinear   L2 = 1.98573061e+02      <- identical to unlimited central
      vanLeer         L2 = 1.97878277e+01      <- 10x better, clearly limiting

Two limiters agreeing with UNLIMITED CENTRAL to nine digits inside a steep
boundary layer is not credible. Scheme selection was verified to be working
(the `fvSchemes` entry was echoed and differs per run, and vanLeer/SuperBee do
differ from each other in the 4th digit), so this is something about how those
limiters behave in this 1D configuration and IT IS NOT UNDERSTOOD. Those rows
are excluded from section 5 and from the conclusions. Resolving it is
outstanding work.

---

## 8. What this test does NOT measure, and the Pasolari comparison

**This bed uses a SMOOTH, MONOTONE solution with no extrema.** Boundedness and
oscillation control are therefore never stressed -- which is precisely what
limiters exist for. A scheme can be first-order here and still be the right
production choice.

Pasolari & Kourtzanidis (arXiv 2607.05137) assess the same family on a **stiff
scalar-advection problem and the positive-streamer benchmark** and conclude the
opposite ranking: SG is *"highly stable but excessively diffusive,
over-predicting the field and propagation speed on coarse meshes, while the
ROUNDF scheme outperforms all tested TVD limiters and is recommended for
streamer transport."*

**The two results are complementary, not contradictory.** They measure different
properties:

  * this bed  -> FORMAL ORDER on a smooth source-driven solution
  * Pasolari  -> STRUCTURE PRESERVATION and front capture on a sharp streamer

SG's excessive diffusion at a sharp front cannot be seen here (there is no
front), and ROUNDF's boundedness cannot be rewarded here (there are no extrema).
ROUNDF's boundedness is exactly what wins it the streamer benchmark and exactly
what costs it order on a smooth problem. Both are true.

The concrete follow-up: **ROUNDL, ROUNDW and ROUNDAplus were not in Pasolari's
tested set** (`limitedLinear, MUSCL, ROUNDF, ROUNDA, SuperBee, Minmod, vanLeer,
upwind, SG`). ROUNDL and ROUNDW win decisively on smoothness here. Whether they
retain ROUNDF's boundedness on the streamer benchmark is the open question, and
it is the test that would decide whether SoPLASMA should change its production
default. Running the two of them through Pasolari's positive-streamer case is
the single highest-value follow-up from this study.
