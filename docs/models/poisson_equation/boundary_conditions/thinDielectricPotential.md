# `thinDielectricPotential` — an unmeshed dielectric layer as a boundary condition

## What it is for

A dielectric barrier thin enough that meshing it would be wasteful, or one that
is not part of the computational domain at all. Instead of a second mesh region
coupled through `coupledElectricPotential`, the layer is collapsed onto the gas
boundary and represented by its **surface capacitance**.

Two physically distinct situations share one implementation:

| situation | what is behind the layer | `thickness` / `backingPotential` |
|---|---|---|
| **free-standing** | open gas or vacuum — nothing | omitted |
| **conductively backed** | a metal electrode held at a known potential | both given |

The free-standing case is the common one for a wall that merely must not let
charge escape; the backed case is the classical DBD stack
*electrode │ dielectric │ gas*.

## Theory

### The backed layer

Let the gas occupy $z>0$ with permittivity $\varepsilon_g=\varepsilon_0\varepsilon_{r,g}$,
and let a layer of thickness $d$ and permittivity $\varepsilon_0\varepsilon_r$
sit between the gas boundary and a conductor at potential $V_b$. Because the
layer is thin, the field inside it is taken as uniform and normal, so the
potential drop across it is linear. Charge conservation on the interface —
Gauss's law applied to a pillbox straddling it — gives

$$\varepsilon_g \frac{\partial V}{\partial n}\bigg|_{\text{gas}} = \sigma - C\,(V_s - V_b),
\qquad C \equiv \frac{\varepsilon_0\varepsilon_r}{d},$$

where $V_s$ is the potential on the gas side of the surface, $\sigma$ the
accumulated free surface charge, and $C$ the **surface capacitance** per unit
area (F m⁻²). The first term is the charge delivered by the plasma, the second
the displacement charge stored in the layer.

This is a **Robin** condition: a linear relation between $V$ and $\partial V/\partial n$.
OpenFOAM expresses such conditions through `mixedFvPatchScalarField`, which
evaluates

$$V_s = f\,V_{\text{ref}} + (1-f)\big(V_c + \Delta^{-1} g_{\text{ref}}\big),$$

with $f$ the `valueFraction`, $V_c$ the adjacent cell value and $\Delta^{-1}$ the
patch `deltaCoeffs`. Matching the two expressions:

$$f = \frac{C}{C + \varepsilon_g \Delta^{-1}}, \qquad
V_{\text{ref}} = V_b + \frac{\sigma}{C}, \qquad g_{\text{ref}} = 0 .$$

The ratio $C/(\varepsilon_g\Delta^{-1})$ compares the layer's capacitance to that
of the first gas cell, which is what makes the two limits below fall out.

### The free-standing layer

With no conductor behind it there is no return path for displacement charge:
$C\to 0$ and the condition degenerates to pure **Neumann**,

$$\varepsilon_g \frac{\partial V}{\partial n} = \sigma
\quad\Longrightarrow\quad f = 0,\qquad g_{\text{ref}} = \frac{\sigma}{\varepsilon_g}.$$

**This branch must be taken separately, not obtained by putting $C=0$ into the
backed formulae.** $V_{\text{ref}} = V_b + \sigma/C$ diverges as $C\to0$ even
though the product $f\,V_{\text{ref}}$ stays finite, so the limit is correct but
unreachable through that expression. The implementation therefore branches on
`C <= 0` and never forms `refValue` there.

### The two limits of the backed branch

- $d\to 0$: $C\to\infty$, $f\to1$, and the condition becomes **Dirichlet**
  $V_s = V_b + \sigma/C \to V_b$ — a bare electrode, as it should be.
- $d\to\infty$: $C\to0$, recovering the free-standing Neumann condition.

A useful consequence: with $\sigma=0$ a gas gap of thickness $L$ in series with
the layer gives the **series-stack** potential

$$V_s = \frac{d\,V_0}{d + \varepsilon_r L}\bigg|_{\varepsilon_{r,g}=1},$$

which is the exact result used to validate the capacitance *value* (below).

## Usage

```cpp
// free-standing: an insulating wall that accumulates charge
air_dielectric
{
    type            thinDielectricPotential;
    epsilonR        4.6;        // the layer's relative permittivity
    surfCharge      surfCharge; // field supplying sigma; `none` for sigma = 0
    value           uniform 0;
}

// conductively backed: the classical DBD stack
barrier
{
    type                thinDielectricPotential;
    epsilonR            4.6;
    thickness           0.5e-3;  // [m]
    backingPotential    0;       // [V] the electrode behind the layer
    surfCharge          surfCharge;
    value               uniform 0;
}
```

| key | default | meaning |
|---|---|---|
| `epsilonR` | `1.0` | relative permittivity of the *layer* |
| `thickness` | *(absent → free-standing)* | layer thickness $d$ [m] |
| `backingPotential` | *(absent → free-standing)* | conductor potential $V_b$ [V] |
| `surfCharge` | `surfCharge` | name of the $\sigma$ field; `none` ⇒ $\sigma=0$ |

`thickness` and `backingPotential` are validated **as a pair** — one without the
other is fatal. Each alone is meaningless: a thickness with nothing behind it is
a capacitance to nowhere, and a backing potential with no layer has nothing to
drop across. The error message states both cases.

$\varepsilon_g$, the *gas* permittivity, is not an input here. It is read from
the region's registered `dielectricProperties`, falling back to the
`electromagneticsModel` — the same lookup `coupledElectricPotential` uses, so the
two conditions cannot disagree about it.

## Validation (measured 2026-09-03)

Unit bed `tutorials/electrostatics/singleRegionElectrostaticFoam/thinDielectricSeriesStack` — a 1D column, gas gap
$L=1\,$m, $V_0=1\,$V at `top`, `sides` `zeroGradient` so the problem is exactly
one-dimensional, $\sigma=0$ so only the capacitance is under test. Ground truth
is analytic, so this is a stronger test than any CFD comparison: run
`./Allrun-sweep`.

| case | $\varepsilon_r$ | $d$ | analytic $V_s$ | measured | rel. err |
|---|---|---|---|---|---|
| a | 5 | 0.5 | 0.0909090909 | 0.0909091 | 1e-7 † |
| b | 5 | 5.0 | 0.5 | 0.5 | 0 |
| c | 2 | 0.5 | 0.2 | 0.2 | 0 |
| d | 5 | free-standing | 1.0 | 1.0 | 0 |

† write precision (`writePrecision 6`), not solver error.

Cases (a)–(c) test the capacitance **value** across both $\varepsilon_r$ and
$d$; case (d) is the $C\to0$ Neumann limit, where zero flux everywhere but the
driven electrode must make $V$ uniform at exactly $V_0$.

The $d\to0$ **Dirichlet** limit was measured separately on the `needleDBD` bed:
`thickness 1e-12` with `backingPotential 500` gave a patch value of 500 V
(`valueFraction` = 0.9999997), confirming the collapse to a bare electrode.

The paired-key guard was checked by giving `thickness` alone — fatal, as
intended.

## Relation to the meshed alternative

Where the barrier **is** meshed as its own region, this condition is the wrong
tool: use `coupledElectricPotential`, which solves Laplace inside the barrier
and needs no thin-layer assumption. Cost aside, the meshed route is more
accurate whenever the field inside the barrier is not uniform and normal —
near an electrode edge, for instance. `thinDielectricPotential` assumes it is,
which is exactly the assumption that makes a capacitance sufficient.

## References

- Boeuf, J. P. (2003). Plasma display panels: physics, recent developments and
  key issues. *J. Phys. D: Appl. Phys.* **36**, R53 — the thin-dielectric
  capacitive surface condition in DBD/PDP modelling.
- Golubovskii, Yu. B. et al. (2003). Influence of interaction between charged
  particles and dielectric surface over a homogeneous barrier discharge in
  nitrogen. *J. Phys. D: Appl. Phys.* **35**, 751 — surface charge accumulation
  on a barrier and its feedback on the gap field.
- Georghiou, G. E. et al. (2005). Numerical modelling of atmospheric pressure
  gas discharges leading to plasma production. *J. Phys. D: Appl. Phys.* **38**,
  R303 — review, including surface boundary treatments.
