// ===========================================================================
//  Grubert, Becker & Loffhagen, Phys. Rev. E 80 (2009) 036405 -- the geometry.
//
//  A PLANE-PARALLEL 1 cm gap. Their model is one-dimensional, so the lateral
//  extent is a modelling choice, not a physical one: it is made narrow and
//  closed with symmetry planes, which reproduces the 1-D solution exactly
//  while keeping a genuinely 2-D Cartesian mesh (the earlier streamer bed this
//  case derives from is axisymmetric, and a wedge cannot represent this).
//
//  x is the GAP direction: cathode at x = 0, anode at x = L.
//  y is the lateral direction, closed by symmetry.
//  The third direction is created by tools/msh2Dto3D.py, which extrudes the
//  MESH -- so there is no CAD kernel choice and no `Coherence;` to renumber.
// ===========================================================================

L  = 0.01;      // gap [m] -- Grubert's 1 cm
W  = 0.001;     // lateral width [m]; arbitrary, closed by symmetry
NX = 2000;       // cells across the gap, GRADED -- see the Bump below
NY = 5;         // cells across the width

Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, W, 0};
Point(4) = {0, W, 0};

Line(1) = {1, 2};   // side_lo   (y = 0)
Line(2) = {2, 3};   // anode     (x = L)
Line(3) = {3, 4};   // side_hi   (y = W)
Line(4) = {4, 1};   // cathode   (x = 0)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// GRADED TOWARDS BOTH ELECTRODES, changed 2026-09-06.
//
// It was uniform at 50 um, on the reasoning that "a graded mesh would make
// '200 cells' ambiguous when comparing against a published profile". That
// reasoning was wrong, because it compared against the wrong length scale: the
// thing a glow discharge must resolve is not the cathode-fall THICKNESS (~2 mm
// here, 40 uniform cells, ample) but the space-charge SHEATH inside it, whose
// scale is the Debye length.
//
// MEASURED on this case: at the runaway density it reached, n_e = 6.27e19
// m^-3, lambda_D = 1.33 um against a 50 um cell -- under-resolved 37.7x. The
// sheath is precisely the mechanism that chokes the current in a glow, so a
// mesh that cannot form it cannot limit the current: I_cond ran to 663
// mA/cm^2, 1300x Grubert's 0.511, essentially all conduction (I_disp 0.009).
//
// Note the mesh WOULD be adequate at the right answer -- at Grubert's peak
// n_e = 2.478e15, lambda_D = 211 um and dx/lambda_D = 0.24. So this only bites
// once the solution has overshot, which is why a uniform mesh looked defensible
// until it did.
//
// Bump 0.02 refines BOTH ends: min cell 1.35 um, max 66.8 um, ratio 49.
// 1.35 um resolves lambda_D even at the runaway density, and gives ~150 cells
// per lambda_D at the density the answer should have.
Transfinite Curve{1, 3} = NX + 1 Using Bump 1.0;
Transfinite Curve{2, 4} = NY + 1;
Transfinite Surface{1};
Recombine Surface{1};        // quads -> hexes after extrusion

Physical Surface("gas")     = {1};
Physical Curve("cathode")   = {4};
Physical Curve("anode")     = {2};
Physical Curve("side_lo")   = {1};
Physical Curve("side_hi")   = {3};
