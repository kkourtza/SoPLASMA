# floatingElectrode -- analytic closed form (unit bed)

## Question
Does the floating-electrode constraint reproduce the exact potential of a
conductor whose CHARGE is known and whose POTENTIAL is the unknown?

## Bed
1D column, gap L = 1 m, gas epsilonR = 1. `top` is a DRIVEN plate at
V0 = 1 V (a non-zero constant Dirichlet). `bottom` is the FLOATING conductor.
`sides` zeroGradient, which is what makes the problem exactly 1D.
Plate area A = 1 m (x) * 0.1 m (z) = 0.1 m^2.

## Ground truth (ANALYTIC)
With V = A y + B on 0 <= y <= 1 and V(1) = V0, the charge on the bottom plate is
    Q = + INT eps snGrad(V) dA,   n_bottom = -y,  snGrad = -A
      => Q = -eps A_area A            (A_area = plate area)
      => A = -Q/(eps A_area)
and therefore
    V_f = V(0) = B = V0 - A = V0 + Q/(eps A_area).

THIS IS NOT A DEGENERATE CASE. The closed form being tested is
    V_f = (Q - Q_rho)/C_self
and here BOTH pieces are non-trivial:
    Q_rho  = -eps*A_area     (the charge induced with the floating patch at 0)
    C_self =  eps*A_area/L
so  V_f = (Q + eps*A_area)/(eps*A_area) = V0 + Q/(eps*A_area)  as above.
An error in the induced-charge term would show up here.

| case | Q0 [C]        | expected V_f |
|------|---------------|--------------|
| q0   | 0             | 1.0  -- zero charge => zero field => the conductor floats to the driving potential |
| qp   | +eps0*0.1     | 2.0          |
| qm   | -eps0*0.1     | 0.0          |

C_self must equal eps0*A_area/L = 8.854187817620389e-13 F.

eps0 = 8.854187817620389e-12, i.e. 1/(mu0 c^2) as plasmaConstants.C builds it.

Errors are judged against the VOLTAGE SCALE V0 = 1 V, never against the expected
value: case `qm` expects exactly 0, where a relative error is undefined.

## Extraction
    ./Allrun-sweep      # measured vs analytic, PASS/FAIL per row; nonzero exit on any failure
    ./Allclean

## Requirement
`singleRegionElectrostaticFoam` does not link libplasmaBcs.so, so
system/controlDict loads it via `libs`.
