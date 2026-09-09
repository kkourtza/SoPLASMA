# A verification map for SoPLASMA

Written 2026-09-07 from DeChant 2023 (PSST 32 044006, "Verification methods for
drift-diffusion reaction models for plasma simulations") and Teunissen 2020
(PSST 29 015010, "Improvements for drift-diffusion plasma fluid models with
explicit time integration"), both in Literature/.

## The finding that motivates all of it

`validation/` holds 60 cases. 58 are Grubert variants. NONE is a verification
case in DeChant's sense -- none has a known answer. We have spent this whole
effort debugging a fully coupled 4-field nonlinear system (n_e, n_Arp, n_Ar2p,
nEps_e + Poisson + chemistry) in which no single component has ever been
verified against an exact solution. DeChant's central point is that this is
backwards: verify components against known answers, THEN run the discharge.

Every defect found on 2026-09-07 was found by hand-instrumenting a diverging
discharge. Each one would have been caught faster, or at all, by a test below.

## Tier 0 -- what already exists

`src/applications/utilities/`: testWallFlux, testWallLoss, testEmission,
testDischargeCurrent, testVibRelax, testAitken, plasmaChemistry0D.

testWallFlux is the model to copy: it checks the SHIPPED algebra against
ANALYTIC ground truth from Hagelaar ch. 6, and it now carries the eq. (6.15)
regression. It is also the cautionary tale -- it had 32 passing checks while the
defect shipped, because every gRatio check used Dd = 0, where the right and
wrong forms coincide. A unit bed only covers what it varies.

## Tier 1 -- analytical solutions (DeChant sec. 4.1). Cheap, exact, decoupled.

  A. CHEMISTRY, LINEAR CHAIN.  DeChant eqs (26)-(31). A -> B -> C with
     first-order rates; exact exponentials. Target: the chemistry ODE
     integrator alone. Host: plasmaChemistry0D. Cost: hours.

  B. CHEMISTRY, NONLINEAR WITH AN EQUILIBRIUM.  DeChant eqs (32)-(33): argon
     ionisation plus three-body recombination. Has a KNOWN STEADY STATE.
     Target: whether the integrator finds and HOLDS an equilibrium -- the
     0D analogue of the question that killed the Grubert runs, where the
     LMEA restoring exponent collapses from 20.9 at 6 eV to 0.20 at 649 eV.
     HIGH VALUE. Host: plasmaChemistry0D.

  C. DEBYE SCREENING AT A WALL.  DeChant eq. (34)-(35): Boltzmann electrons,
     Laplace(V) = (e n0/eps0)(exp(V/Te) - 1), exact V = V0 exp(-|x|/lambda_De)
     for V << Te. Target: the Poisson solve and the near-wall field.
     HIGHEST VALUE OF TIER 1. The 2026-09-07 anode failure was a 49x
     single-cell E spike where lambda_De was 683 um = 137 cells. This test
     asks the one question that failure poses -- does the code screen a wall
     correctly -- and it has an exact answer. It ALSO tests the semi-implicit
     potential treatment that Teunissen indicts (see Tier 4).

  D. AMBIPOLAR DIFFUSION.  DeChant eqs (36)-(39): with E = (Di-De)/(mui-mue)
     grad(n)/n the two-species flux collapses to Gamma = Da grad(n), and
     -Da Laplace(n) = G0 has the exact parabola
     n = (G0 l^2/8Da)[1 - (2x/l)^2] + 1.
     Target: the drift-diffusion flux assembly. AND IT IS THE NATURAL TEST
     FOR fluxScheme: `standard` and `ScharfetterGummel` must BOTH reproduce
     the same parabola. That is a two-line test that would immediately expose
     the SIGFPE the SG path hits today, in a case with no chemistry to blame.

  E. GLOBAL (0D) PARTICLE BALANCE.  DeChant eq. (40), Lieberman sec. 10.2:
     ng*deff = uB/Kiz(Te) with uB = sqrt(e Te/M). Given the gas and geometry
     this FIXES Te. Target: the closed loop of ionisation rate against wall
     loss. The closest analytical anchor to a steady discharge, and a direct
     check on whether our tables + wall conditions can hold a steady state at
     all. Note DeChant's own caveat: it holds only for an ionisation-only
     source.

  NOT APPLICABLE: DeChant's electronegative-bulk case (eqs 41-44) is an oxygen
  problem. Argon is electropositive.

## Tier 2 -- method of manufactured solutions (DeChant sec. 3.2, 4.2)

  F. MMS ON THE COUPLED SYSTEM.  Pick a smooth solution (trig functions, so
     every derivative is non-trivial), substitute it into the PDEs, and read
     off the source terms that make it exact. DeChant's rules: the solution
     need not be physical, but it MUST exercise every term, and terms should
     be of COMPARABLE MAGNITUDE -- otherwise the test is dominated by one
     term and says nothing about the others.

     THIS IS THE ONLY METHOD THAT TESTS THE LMEA ENERGY EQUATION COUPLED TO
     EVERYTHING ELSE. It is also where the 2026-09-07 measurement points: at
     the failure the energy equation's div and lap terms were each ~680x the
     Joule source. An MMS with comparable-magnitude terms is exactly the
     instrument for asking whether that assembly is right.

     Caveat worth stating: MMS verifies that the code solves the equations it
     claims to solve. It CANNOT catch a wrong equation. The eq. (6.15) defect
     was a wrong equation faithfully implemented -- MMS would have passed it.
     Analytical solutions and source-checking catch those; MMS does not.

## Tier 3 -- convergence studies (DeChant sec. 3.3)

  G. SPATIAL AND TEMPORAL ORDER.  ||e|| = C h^p and C dt^p, L2 and Linf,
     against an exact or Richardson-extrapolated solution. Expected p = 2 for
     BDF2 + second-order space.

     DIRECTLY RELEVANT: Teunissen measures the SEMI-IMPLICIT field treatment
     -- which SoPLASMA uses -- at roughly FIRST-order convergence, with errors
     larger than the alternative "also for time steps larger than the
     dielectric relaxation time". If our temporal order is 1 and not 2, that
     is measurable in an afternoon and would reframe every adaptive-dt result
     we have.

## Tier 4 -- replicating the DC glow, which is the actual goal

  H. TEUNISSEN'S SEMI-IMPLICIT ARTIFACT TEST.  His sec. 3.4: 1D, 10 mm,
     dx = 20 um, 10 kV, SOURCE TERM ZERO, constant mu_e = 0.03, De = 0.1,
     n_e = n_p = 1e20 m^-3 in 4-6 mm and zero elsewhere. Reference: same grid,
     tiny dt. His result: the semi-implicit method "predicts LARGER PEAKS IN
     THE ELECTRIC FIELD AT THE BOUNDARIES OF THE INITIALLY IONIZED AREA".
     That is our anode signature, in a case with NO chemistry and NO tables to
     blame. If SoPLASMA reproduces the spurious peak here, the 49x spike is a
     semi-implicit artifact and the tables/BCs are exonerated.
     CHEAPEST HIGH-INFORMATION TEST ON THIS LIST. Diagnostic only -- we are
     NOT adopting his current-limited scheme.

  I. TEUNISSEN'S PARALLEL-DIFFUSION DIAGNOSTIC.  His sec. 4: electrons diffuse
     parallel to E out of a screened high-density region into the high-field
     region and ionise there, growing the plasma unphysically; "using a finer
     grid spacing SLOWS DOWN the unphysical growth but DOES NOT PREVENT IT."
     That matches our mesh study exactly (graded -> uniform 5 um moved the
     failure 0.03%). His f_eps = 1 - (Ehat . Gamma_diff)/|Gamma_drift| can be
     computed as a PURE DIAGNOSTIC on stored fields -- where it goes to zero
     or negative is where the model is generating ionisation it should not.
     NOTE: he states this is LFA-specific. Under LMEA the ionisation is keyed
     on meanE, not local E, so the mechanism is damped by energy relaxation --
     which is Grubert's whole argument for LMEA. Whether an LMEA analogue
     exists (electrons carrying nEps into the high-field region) is open.
     Compute the diagnostic before believing either way. DO NOT IMPLEMENT
     f_eps as a model change.

  J. PASCHEN CURVE.  Sweep pd, find breakdown voltage, compare to the
     published argon curve. We already measured V_b = 120.8 V for this gap.
     Cheap, physical, and tests ionisation + secondary emission + the
     Townsend criterion as a closed loop, with a literature answer.

  K. NORMAL CATHODE FALL SIMILARITY (the strongest DC-glow anchor, and NOT
     from either paper). A normal glow has THREE tabulated similarity
     constants for a given gas/cathode pair: the normal cathode fall voltage
     V_n, the normal thickness (pd_c)_n, and the normal current density
     j_n/p^2. For argon these are in von Engel and in Raizer's "Gas Discharge
     Physics" ch. 8. This is a STEADY STATE with published numbers, in exactly
     the regime Grubert 2009 sits in, and it tests the cathode-fall structure
     -- which is where the 2026-09-07 failure actually developed (hot spot at
     x = 1.26 mm, E/N 63 kTd, while the bulk was at 1.6 kTd).
     If SoPLASMA cannot reproduce a normal cathode fall, Grubert is out of
     reach and this says so with far less machinery.

  L. CARLSSON ET AL. GD BENCHMARK (He, 3.5 Torr, 0.62 cm, JC-PIC/Boeuf
     reproduction).  A PUBLISHED, PIC-MCC-VALIDATED DC glow benchmark, more
     tractable than K because it comes with independently reproduced numbers
     from TWO kinetic codes (EDIPIC, LSP) plus JC-PIC's own reproduction, not
     just an analytic similarity law. Conditions: He at 3.5 Torr (300 K), gap
     0.62 cm, cathode voltage -211 V, gamma = 0.28 (SEE, ion bombardment),
     emitted-electron temperature 2 eV. Published targets: field at the
     cathode surface E0 and the cathode-fall width d_c (JC-PIC matches EDIPIC
     and LSP to ~2% on E0 and ~0.02 cm on d_c; all three codes read ~10%
     above the DenHartog/O'Brian/Lawler experimental E0).
     WHY THIS ONE, NOT JUST GRUBERT: gives an independent, non-Grubert DC-glow
     target with numbers from THREE different kinetic treatments already
     agreeing with each other, so a SoPLASMA (fluid) reproduction is compared
     against a converged consensus rather than a single paper's numbers.
     CAVEAT (from the same source): the higher-voltage 600 V case in this same
     benchmark is explicitly UNRESOLVED even in JC-PIC's own kinetic
     reproduction -- Carlsson himself calls it "only a partial success", and
     JC-PIC's attempt did not converge at 150 us of simulated time (their
     estimate: needs ~10x longer, into the ms range, because the cold trapped
     electron population in the negative glow builds up on the ambipolar
     diffusion timescale). Reproduce the 211 V case first; do not expect quick
     convergence at 600 V even in principle.
     Also worth noting: NO electron-electron (Coulomb) collision term exists
     anywhere in SoPLASMA (see doc/electron-electron-collisions-gap.md) --
     the mechanism the source names as dominant for de-trapping the cold
     negative-glow population. A SoPLASMA reproduction may therefore plateau
     at the wrong bulk Te/density regardless of run length; comparing against
     this benchmark is exactly the test that would show it.
     Reference: Literature/reference-codes-and-manuals/JC-PIC_Boeuf/library_article_book.txt, section V.A
     (search "GD Benchmark").

## Suggested order

  1. D (ambipolar) -- exposes the SG SIGFPE with no chemistry in the way.
  2. H (Teunissen semi-implicit) -- decides whether the E spike is an artifact.
  3. C (Debye screening) -- exact near-wall field.
  4. G (convergence order) -- is the semi-implicit path first order here?
  5. K (normal cathode fall) -- the DC-glow anchor.
  6. L (Carlsson/JC-PIC GD benchmark, 211 V case) -- an independent,
     multi-code-consensus DC-glow target; do this alongside or right after K.
  7. B, E, F, J as the ladder fills in.

1-3 are each roughly a day and each can falsify a current hypothesis. K is the
one that answers "can this code do a DC glow at all", independently of Grubert.
