# COMPARE.md — plate2D_timeVaryingBC_implicitBoundary

Written 2026-09-01, when this case was recruited as an acceptance test for the
region/material refactor (permittivity moved from a single global
`constant/electromagneticsProperties` to per-region
`constant/<region>/electricalProperties`).

## The question this case answers

**Is each region's permittivity actually read from its own
`constant/<region>/electricalProperties`, under MONOLITHIC coupling?**

The failure mode of the replacement is silent: if the per-region file is
not found and the legacy fallback does not fire either, `epsilonR`
defaults to 1.0, the run converges, and the field is smooth and
plausible. This case discriminates that by a factor of 3.

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
| coupling | `configuration/config` `useImplicit` | `true` — MONOLITHIC |

## Reference numbers

| quantity | value |
|---|---|
| **`V_interface` at t = 10 (analytic)** | **0.166666667 V  (= 1/6)** |
| measured 2026-09-01 | 0.166667 V |
| relative error | 2.0e-06 |
| if a permittivity silently defaulted to 1.0 | 0.5 V — a factor of 3 |

Tolerance is **5e-6 relative**, which is the `writePrecision 6` floor:
the field file holds 6 significant digits, so 1/6 is written as 0.166667
and cannot agree better than 2e-6 however well the run converged. The
solution is linear in each layer and the interface is a cell FACE, so
there is no discretisation error to allow for — anything at 1e-2 is a
real defect, not a mesh effect.

## Extraction command

```bash
cd /home/kkourtza/soplasma-scratch/tutorials/electrostatics/multiRegionElectrostaticFoam/plate2D_timeVaryingBC_implicitBoundary
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

# monolithic assembly actually happened
grep "MONOLITHICALLY" logs/log.*serial
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
