# Mesh convergence for the Grubert dc glow -- PREPARED, NOT LAUNCHED

Written 2026-09-07. Sequenced after the floor study at the user's instruction:
the floor value must be CHOSEN first, or the mesh study is confounded by an
unsettled initial condition.

## Why this is open at all

`RESOLUTION_TEST.md` designed this properly on 2026-09-06 with a pre-registered
prediction, and **it never produced evidence**: the graded-mesh arms
(`grubert2009`, `grubert2009_LFA`) died at t = 1.83e-7 and 6.32e-7 s --
pre-breakdown, `n_e` still at the floor, peak current 0.009 mA/cm^2. So no
mesh convergence study on this case has EVER been completed.

What IS known, from the uniform 50 um mesh: `lambda_D` = 1.33 um at the
runaway density, **under-resolved 37.7x**, and the current reached 663 mA/cm^2
(1300x Grubert). So resolution demonstrably CAN dominate this case's answer.

And there is a specific reason to suspect the graded mesh is refined in the
WRONG PLACE. `gap1cm.geo` grades toward the ELECTRODES on the argument that the
sheath scale is `lambda_D` at the electrode. But the structures that actually
form are not at the electrodes: the double layer sits at x/L = 0.0216 (215 um
out) and the second field reversal at x/L = 0.721 (mid-gap), where the Bump has
already coarsened.

## The four arms, all built by `tools/make_mesh_arm.sh`

| arm | NX | Bump | cells (1-D) | dx min | dx max | ratio | cells/lambda_D at the double layer | at the mid-gap reversal |
|---|---|---|---|---|---|---|---|---|
| **M1** | 400 | 0.02 | 400 | 1.35 um | 66.8 um | 49.3 | **5.05** | 11.8 |
| **M2** | 800 | 0.02 | 800 | 0.67 um | 33.4 um | 49.7 | 10.0 | 23.5 |
| **M4** | 1600 | 0.02 | 1600 | 0.34 um | 16.7 um | 49.8 | 20.1 | 47.1 |
| **MU** | 800 | 0.2 | 800 | 4.05 um | 20.2 um | 5.0 | 6.4 | **37.3** |

`lambda_D` values are those measured on the frozen spike state (34.6 um at the
double layer, 636 um at the reversal). M1 independently reproduces the 5.1
cells/lambda_D measured from the profile data -- a cross-check on both.

**TWO QUESTIONS, deliberately separated:**

* **Convergence** -- M1 -> M2 -> M4, same grading law, 2x and 4x cells. Does the
  answer move?
* **Placement at FIXED COST** -- M2 vs MU, both 800 cells. MU trades electrode
  resolution (10.0 -> 6.4 at the double layer) for mid-gap (23.5 -> 37.3). Does
  it matter WHERE the cells are?

## M1 MUST BE REGENERATED, not taken from an existing run

Measured 2026-09-07: an existing case's `0/` accumulates BC keys the solver
writes back (`refValue`, `refGradient`, `valueFraction`, `source`, `enableSEE`,
`defaultSEEC`) which a freshly generated `0/` does not carry -- see
[[defaults-must-not-be-written-back]]. Using a run's `0/` as the M1 baseline
would confound "mesh refined" with "`0/` regenerated". All four arms are
therefore built by the same script from the same parent, so the ONLY difference
is `NX`/`Bump`.

The generated set is 6 primary fields (`ePotential`, `nEps_e`, `n_e`, `n_Arp`,
`n_Ar2p`, `surfCharge`); the other ~18 in a run's `0/` are derived fields the
solver wrote. 6 is correct.

## Success and failure, stated in advance

* **CONVERGED** -- M1/M2/M4 agree on breakdown time, peak `j`, `n_e` max and
  whether `dt` collapses. The mesh is then settled and the physics stands as
  measured. This is the outcome that would let every earlier conclusion keep its
  numbers.
* **NOT CONVERGED, monotone improvement** -- peak `j` and `n_e` fall with
  refinement toward Grubert's values. Resolution is then a real part of the
  discrepancy and the earlier runaway magnitudes are mesh artefacts.
* **PLACEMENT MATTERS** -- MU differs from M2 at equal cost. The grading law is
  then wrong for this solution and should follow the structures, not the
  electrodes.
* **NO CHANGE ANYWHERE** -- refinement is excluded as an explanation.
  `RESOLUTION_TEST.md` pre-committed the discipline here and it still applies:
  explanation (1) is then RETRACTED, not patched with a further level.

## Cost note

M4 is 4x the cells of M1 and `dt` will fall with the Courant limit, so expect
~8-16x the wall time of an M1 arm. Run M1/M2/MU first and add M4 only once the
first three show whether anything moves.

## Two bugs found while testing the setup script, both recorded in it

1. `set -e` plus OpenFOAM's `etc/bashrc`, which returns non-zero, killed the
   script SILENTLY -- empty `logs/`, no message, no output at all. Replaced with
   per-step `die()` guards plus a positive check that `gmshToFoam` is on PATH.
2. `export SoPLASMA=... SoPLASMA_ETC=$SoPLASMA/etc` in ONE statement expands
   `$SoPLASMA` before assigning it, so `SoPLASMA_ETC` became `/etc` and
   `plasmaSetupBoundaries` died on `Cannot read the boundary-role library:
   "/etc/boundaryRoles"`. Separate lines, plus an assertion that
   `$SoPLASMA_ETC/boundaryRoles` exists. `Allrun-setup` already used separate
   lines for this reason.

The Boltzmann and ion tables are mesh-INDEPENDENT and are reused from the
parent; `genMechTables`/`ionmob` are deliberately not re-run.
