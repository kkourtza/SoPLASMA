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
NX = 200;       // cells across the gap: dx = 50 um
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

// Structured: the gap direction must be uniform, because the cathode fall is
// resolved by cell count and a graded mesh would make "200 cells" ambiguous
// when comparing against a published profile.
Transfinite Curve{1, 3} = NX + 1;
Transfinite Curve{2, 4} = NY + 1;
Transfinite Surface{1};
Recombine Surface{1};        // quads -> hexes after extrusion

Physical Surface("gas")     = {1};
Physical Curve("cathode")   = {4};
Physical Curve("anode")     = {2};
Physical Curve("side_lo")   = {1};
Physical Curve("side_hi")   = {3};
