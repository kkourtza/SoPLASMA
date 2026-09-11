# Complete Flux Scheme -- implementation status

2026-09-08. **WORKING AND VERIFIED at the operator level.** Not yet wired into
the solver; see "Remaining work".

## What it is

CFS derives the face flux from the local two-point problem

    d/dx ( v n - D dn/dx ) = s,   x_j < x < x_j+1

WITH THE SOURCE RETAINED. Scharfetter-Gummel solves the same problem with
s = 0. The flux therefore splits (Liu et al 2014, PSST 23 015023, eq. 18) into
a homogeneous part -- which IS the SG flux -- plus an inhomogeneous part
carrying s. Gamma^i depends on s only, never on n, so it is EXPLICIT.

## Verified results (verification/fluxScheme1D)

  1. UNIFORM SOURCE, Pe = 1e4, where Gamma^i is exact and SG's error is
     analytically S h/(2v):

         N     SG              CFS
         20    2.50000000e-06  1.110e-16
         40    1.25000000e-06  7.850e-17
         80    6.25000000e-07  6.568e-17
         160   3.12500000e-07  1.314e-16

     CFS is EXACT TO MACHINE PRECISION.

  2. MANUFACTURED SOLUTION n = 1 + x + 0.5 sin(pi x), source derived, at
     domain Pe = 1e4 (grid Peclet 500 -> 15.6) where SG is FIRST order:

         N     Pe_grid  SG          p      CFS         p
         40    250      1.860e-02   0.99   9.123e-05   2.02
         80    125      9.339e-03   0.99   2.275e-05   2.00
         160    62.5    4.690e-03   0.99   5.683e-06   2.00
         320    31.2    2.364e-03   0.99   1.421e-06   2.00
         640    15.6    1.204e-03   0.97   3.551e-07   2.00

     CFS RESTORES SECOND ORDER where SG is first order -- Liu's "second order
     uniformly in Peclet" claim, reproduced. 3390x more accurate at N = 640
     (147x at Pe = 100). Liu advertise ~10x on a glow discharge; the gap grows
     as h vs h^2, so this is consistent.

## THE TWO CORRECTIONS THAT MADE IT WORK

Both were found by re-deriving Gamma^i from scratch rather than transcribing
Liu's coefficients. Solving the local BVP for a general source gives

    Gamma^i = S(x_face) - [lambda/(1-e^-Pe)] Int_0^h S(x) e^{-lambda x} dx,
    S(x) = Int_0^x s.

  (a) THE SECOND-MOMENT TERM. Liu's eq. (18c), Gamma^i = h(gamma s_j +
      delta s_j+1), is EXACT ONLY FOR A CONSTANT SOURCE. For s linear across
      the interval the exact result is

          Gamma^i = h [ s_L (1/2 - W(Pe)) + (s_R - s_L)(1/8 - V(Pe)/2) ]
          V(z)    = [2 - e^{-z}(2 + 2z + z^2)] / [z^2 (1 - e^{-z})]

      At large Pe, V -> 0 and (18c) drops h(s_R - s_L)/8 ~ h^2 s'. Checked
      against pure diffusion at Pe -> 0, where the exact inhomogeneous flux is
      -h(s_R-s_L)/24 and 1/8 - V(0)/2 = 1/8 - 1/6 = -1/24. Agrees.

  (b) THE BOUNDARY COEFFICIENTS ARE NOT THE INTERIOR ONES, and this was the
      decisive one. For an INTERIOR face the local problem spans the two cell
      centres and the face is at its MIDPOINT, so S(h/2) = s h/2 supplies the
      1/2 in (1/2 - W). At a BOUNDARY face the interval runs from the cell
      centre to the face and the flux is evaluated at its END, where
      S(H) = s H gives 1:

          interior:  h [ s_L (1/2 - W) + (s_R - s_L)(1/8 - V/2) ]
          boundary:  H [ s_C (1 - W)   + (s_F - s_C)(1/2 - V/2) ]

      Using the interior coefficients at the boundary halves that term, and in
      a drift-dominated problem the inflow boundary SETS the error: CFS then
      removed EXACTLY half of SG's error (measured ratio 0.536, 0.518, 0.509,
      0.505 -> 1/2). That clean factor of two is what identified it.

  NOT the cause, contrary to an earlier note here: the vertex- vs cell-centred
  grid arrangement. The interior derivation reproduces (18c) exactly from first
  principles on a cell-centred mesh, and eq. (19e)'s weighted diffusivity
  D = Dtilde*Ptilde/Pbar is a NO-OP for constant coefficients, since
  W(z) + W(-z) = (B(-z) - B(z))/z = 1 by the Bernoulli identity.

  A TEST-DESIGN ERROR worth remembering: the first acceptance test used a
  UNIFORM source, which cannot test CFS at all -- Gamma^i is then identical on
  every interior face, div(Gamma^i) == 0, and CFS collapses onto SG in the
  interior. The uniform source chosen to ISOLATE the source term is the least
  favourable case for the scheme that fixes it. Hence the `-mms` manufactured
  solution, which is a real improvement to the bed regardless of CFS. The
  second error was testing at domain Pe = 1 and 100, where the grid Peclet is
  << 1 at the fine end and SG is ALREADY second order -- there is nothing there
  for CFS to fix. CFS must be judged where SG degrades, Pe_grid >> 1.

## Remaining work before it can be used in a case

  1. Wire it as a THIRD option, `fluxScheme (standard | ScharfetterGummel |
     CompleteFlux)`, in driftDiffusion and in localEnergyEnergyModel. NOT a
     mode of SG: it differs from SG only by Gamma^i, and folding it in would
     hide the one term that distinguishes them.
  2. Feed it the NET chemistry source -- chemP_[s] - chemL_[s]*n plus Sph and
     any explicit source, NOT ionisation alone -- LAGGED one outer iteration,
     since plasmaTransport::solve() builds the transport matrix (~line 1081)
     before the sources are added (~1140-1181). chemP_/chemL_ already persist
     as members, so no new storage, and the lag vanishes at PIMPLE convergence.
     Liu treat the source explicitly for the same reason.
  3. 2D/3D: the derivation above is 1D along the face normal. Non-orthogonal
     meshes need thought before this is trusted off a Cartesian grid.
  4. Boundary conditions other than fixedValue have not been exercised.
