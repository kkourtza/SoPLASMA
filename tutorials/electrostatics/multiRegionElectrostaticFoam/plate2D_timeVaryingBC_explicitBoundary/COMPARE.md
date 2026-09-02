# COMPARE.md — plate2D_timeVaryingBC_explicitBoundary

Written 2026-09-01, when this case was recruited as an acceptance test for the
region/material refactor (permittivity moved from a single global
`constant/electromagneticsProperties` to per-region
`constant/<region>/electricalProperties`).

## The question this case answers

**What does SEGREGATED interface coupling cost, against the same
analytic reference?**

This case differs from `plate2D_timeVaryingBC_implicitBoundary` in
`useImplicit` and nothing else, so the two together isolate the coupling
mode. It also confirms the per-region permittivity read on the
segregated path, which is separate code.

## Why this is a valid control

The reference is **ANALYTIC**, not another run of the same code, so it tests the
physics path itself rather than showing that two diagnostics agree. A 1-D
two-layer stack with uniform layers has a closed form:

```
eps0 (eps2 Vi/d2 + eps1 Vi/d1) = sigma + eps0 eps1 V0/d1
```

from continuity of the normal displacement with a free surface charge,
`D2 - D1 = sigma`, plus `V(0) = V0` and `V(d1+d2) = 0`.

**The conditions that make it exact here.** If any of these changes the numbers
below are void — re-derive them, do not reuse them:

| condition | where it is set | value |
|---|---|---|
| gas thickness `d1` | `configuration/geometry/mesh.py`, `x_mid = x_end/2` | 0.5 m |
| dielectric thickness `d2` | same | 0.5 m |
| `epsilonR` gas | `constant/leftGas/electricalProperties` | 1.0 |
| `epsilonR` dielectric | `constant/rightDielectric/electricalProperties` | 5.0 |
| applied `V0` at `left`, at `endTime` | `configuration/config` `dirichletValueLeft` | 1.0 V |
| `right` patch | `etc/changeDictionary.rightDielectric` | 0 V |
| top / bottom | both `changeDictionary` files | `zeroGradient` ⇒ 1-D |
| space charge | `chargeDensity` `uniform 0` | none |
| interface surface charge | `configuration/config` `surfCharge` | 0 |
| coupling | `configuration/config` `useImplicit` | `false` — SEGREGATED |

## Reference numbers

| quantity | value |
|---|---|
| **`V_interface` at t = 10 (analytic)** | **0.166666667 V  (= 1/6)** |
| measured 2026-09-01 | 0.166652 V |
| relative error | 8.8e-05  — 44x the monolithic error, and this IS the measurement |
| if a permittivity silently defaulted to 1.0 | 0.5 V — a factor of 3 |

Tolerance is **5e-4 relative**, not the 5e-6 used on the monolithic
variant, and the difference is the point. Segregated coupling lags the
interface — each region sees the other's potential from the previous
corrector — so it converges to the coupling tolerance rather than to the
write precision. Judging it at the monolithic tolerance would report a
correct run as a failure.

## Extraction command

```bash
cd /home/kkourtza/soplasma-scratch/tutorials/electrostatics/multiRegionElectrostaticFoam/plate2D_timeVaryingBC_explicitBoundary
./Allrun-serial
/home/kkourtza/soplasma-scratch/tools/check_series_stack.py
```

Also check, in `logs/log.multiRegionElectrostaticFoam.serial`:

```bash
# the model is DERIVED from regionProperties, not read from a file
grep "electromagneticsModel:" logs/log.*serial
#   -> multiRegionPoisson (DERIVED -- 1 dielectric region(s) ...)

# no legacy fallback fired -- this case is migrated, so this must be SILENT
grep -E "DEPRECATED|no longer read" logs/log.*serial

# the linear solve CONVERGED. Reaching `maxIter` is not an error and is
# reported as a normal solve; it cost 7.3% on this tutorial before 2026-09-01.
grep "Solving for ePotential, Initial" logs/log.*serial | tail -3

# segregated, so the log must say NOT coupled
grep "Regions are" logs/log.*serial
```

## Baselines

- **Analytic**, as above. Primary.
- **No run predating 2026-09-01 is a valid baseline for the permittivity
  question**, because every earlier run read the permittivity from the global
  file that no longer exists. Earlier runs remain valid baselines for the fields
  being otherwise unchanged — the refactor was not meant to change the answer,
  only where the number is read from.
- The other two `plate2D_*` cases are valid controls for **one variable each**:
  `explicitBoundary` differs from `timeVaryingBC_implicitBoundary` ONLY in
  `useImplicit`, and `surfCharge_...` differs ONLY in `dirichletValueLeft` and
  `surfCharge`. Everything else — mesh, permittivities, BCs — matches.
