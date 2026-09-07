# Mesh convergence on the Grubert dc glow -- M1 / M2 / MU

Created 2026-09-07. Design and provenance: `../../doc/mesh-convergence-plan.md`.

## The question

**Does mesh resolution change the answer on this case, and does it matter WHERE
the cells are?** No mesh convergence study on this case has ever been completed:
`../RESOLUTION_TEST.md` designed one and its arms died pre-breakdown at
t = 1.8e-7 / 6.3e-7 s, producing no evidence. What IS known is that resolution
CAN dominate -- on a uniform 50 um mesh `lambda_D` was under-resolved 37.7x and
the current reached 663 mA/cm^2, 1300x Grubert.

## Baseline -- ABSOLUTE PATHS

    /home/kkourtza/soplasma-scratch/validation/grubert2009_ballast_low
        the original floor-1e11 run: healthy to V_src = -96.6 V drawing 1% of
        I_sc, then 95.6x I_sc and V_el -> +204 V in 0.76 us, frozen at
        t = 1.04291e-05 s.

    /home/kkourtza/soplasma-scratch/validation/grubert2009_floor1e9
        the CHOSEN floor, and the parent all three arms were built from.

## The floor is FIXED at 1e9, and that choice is measured

Floor 1e9 and 1e7 agreed on `n_e max/min` to **0.19%** at t = 6.35 us, while
both differed from 1e11 (0.813x its ratio there and falling). So the answer is
floor-INDEPENDENT below 1e9 and floor-DEPENDENT above it. 1e9 is picked because
it is converged AND sits 100x further from the `meanE = nEps/n_e` conditioning
risk in near-empty cells than 1e7 would.

## The three arms -- ONE variable

| arm | NX | Bump | cells | dx min | dx max | cells/lambda_D at double layer | at mid-gap reversal |
|---|---|---|---|---|---|---|---|
| **M1** | 400 | 0.02 | 2000 | 1.35 um | 66.8 um | **5.05** | 11.8 |
| **M2** | 800 | 0.02 | 4000 | 0.67 um | 33.4 um | 10.0 | 23.5 |
| **MU** | 800 | 0.2 | 4000 | 4.05 um | 20.2 um | 6.4 | **37.3** |

`lambda_D` reference values are those measured on the frozen spike state
(34.6 um at the double layer at x/L = 0.0216, 636 um at the reversal at
x/L = 0.721). All three built by `tools/make_mesh_arm.sh` from the same parent;
`diff -rq` confirms nothing outside the mesh files differs.

**M1 vs M2** = convergence at fixed grading law.
**M2 vs MU** = PLACEMENT at fixed cost (4000 cells each). MU trades electrode
resolution for mid-gap, which is where the structures actually are: the
`gap1cm.geo` grading refines toward the ELECTRODES, but the double layer sits
215 um out and the second reversal at mid-gap.

M4 (NX 1600) is deliberately HELD until these three show whether anything moves
-- it is 4x the cells with a smaller dt, so ~8-16x an M1 arm.

## M1 also serves a second purpose

M1 reproduces `grubert2009_floor1e9`'s mesh and floor exactly, differing only in
that its `0/` is FRESHLY GENERATED rather than carrying the BC keys a solver
writes back (`refValue`, `refGradient`, `valueFraction`, `source`, `enableSEE`,
`defaultSEEC`). If M1 and `grubert2009_floor1e9` agree, that write-back is
harmless here; if they differ, it is not, and every case built by copying a
completed run is suspect.

## Success and failure, stated in advance

* **CONVERGED** -- M1/M2 agree on breakdown time, peak `j`, `n_e` max and
  whether `dt` collapses. The mesh is settled and the physics stands as
  measured. Earlier conclusions keep their numbers.
* **NOT CONVERGED, MONOTONE** -- peak `j` and `n_e` fall with refinement toward
  Grubert's values. Resolution is then a real part of the discrepancy and the
  recorded runaway magnitudes are partly mesh artefacts. Launch M4.
* **PLACEMENT MATTERS** -- MU differs from M2 at equal cost. The grading law is
  wrong for this solution and should follow the structures, not the electrodes.
* **NO CHANGE ANYWHERE** -- refinement is EXCLUDED. Per the discipline
  `RESOLUTION_TEST.md` pre-committed, explanation (1) is then RETRACTED, not
  patched with a further refinement level.

## Discriminating observables

| observable | where |
|---|---|
| does it pass t = 1.04291e-05 (parent freeze) | `grep -a '^Time = ' logs/log.soPlasmaFoam \| tail -1` |
| `\|I_cond\|` vs `I_sc = \|V_src\|/R`, R = 1e8 | `postProcessing/externalCircuit/circuit.csv` |
| `n_e` max/min vs time | `grep -aF "e [m^-3]:" logs/log.soPlasmaFoam` |
| mean-energy range | `grep -a "LMEA mean energy" logs/log.soPlasmaFoam` |
| dt and WHICH limiter set it | `grep -a "deltaT set by" logs/log.soPlasmaFoam` |

Judge on MATCHED PHYSICAL TIMES, not step counts -- the arms take different
numbers of steps by construction.
