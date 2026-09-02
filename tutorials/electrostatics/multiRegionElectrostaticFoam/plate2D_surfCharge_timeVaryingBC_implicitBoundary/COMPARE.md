# COMPARE.md — plate2D_surfCharge_timeVaryingBC_implicitBoundary

Written 2026-09-01, when this case was recruited as an acceptance test for the
region/material refactor (permittivity moved from a single global
`constant/electromagneticsProperties` to per-region
`constant/<region>/electricalProperties`).

## The question this case answers

**Is the interface surface-charge jump condition satisfied, at a
realistic field?**

This case differs from `plate2D_timeVaryingBC_implicitBoundary` only in
`dirichletValueLeft` (1e7 V instead of 1 V) and `surfCharge` (1e-2 C/m^2
instead of 0), so it isolates the surface-charge term. That term
DOMINATES here: sigma/eps0 = 1.13e9 against eps1 V0/d1 = 2e7, so the
interface potential is 9.58e7 V rather than the 1.67e6 V it would be at
zero surface charge — a 57x shift. A run that dropped the surface charge
could not be mistaken for a converged one.

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
| applied `V0` at `left`, at `endTime` | `configuration/config` `dirichletValueLeft` | 1e7 V |
| `right` patch | `etc/changeDictionary.rightDielectric` | 0 V |
| top / bottom | both `changeDictionary` files | `zeroGradient` ⇒ 1-D |
| space charge | `chargeDensity` `uniform 0` | none |
| interface surface charge | `configuration/config` `surfCharge` | 1e-2 C/m² |
| coupling | `configuration/config` `useImplicit` | `true` — MONOLITHIC |

## Reference numbers

| quantity | value |
|---|---|
| **`V_interface` at t = 10 (analytic)** | **95784088.9 V** |
| measured 2026-09-01 | 95784100 V |
| relative error | 1.15e-07 |
| if a permittivity silently defaulted to 1.0 | 287352267 V — a factor of 3 |

Tolerance is **5e-6 relative**, the `writePrecision 6` floor. The
measured 1.15e-7 is below it, i.e. exact to everything the written file
can express.

## Extraction command

```bash
cd /home/kkourtza/soplasma-scratch/tutorials/electrostatics/multiRegionElectrostaticFoam/plate2D_surfCharge_timeVaryingBC_implicitBoundary
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
