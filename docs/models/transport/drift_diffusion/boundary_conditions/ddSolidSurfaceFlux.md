> # SUPERSEDED 2026-09-04 — read this first
>
> **The class names in this document no longer exist.** `ddSolidSurfaceFlux`,
> `electronDDSolidSurfaceFlux`, `ionDDSolidSurfaceFlux` and
> `neutralDDSolidSurfaceFlux` were replaced by two families:
> `ddWallFluxMixed` (the default: `electron`/`energy`/`ion`/`neutral` members)
> and `ddWallFluxImplicit` (`electron`/`ion`/`neutral` only). Selected with
> `wallFluxFamily Mixed | Implicit`.
>
> **Two physics statements below are also superseded:**
>
> 1. Every `1/4 v_th` in this document — the wall loss speed is now Hagelaar
>    **eq. (6.6)**, `vT/sqrt(pi)`, which is **twice** the `(1/4) sqrt(8kT/pi m)`
>    of eq. (6.3). Changed in `7e6a4b3`. The shifted Maxwellian, not the
>    centred one.
> 2. `Gamma_e,wall = u_wall n_e - Gamma_SEE` — the emission now also enters
>    **inside** the loss speed via eq. (6.8), so the net emission contribution
>    carries a factor `2/(1+r)`. Changed in `70f6674` (Mixed) and 2026-09-04
>    (Implicit).
>
> **Not covered here at all:** electron reflection (`electronReflection`), the
> `max(..., 0)` clamp of eq. (6.6), and the electron-energy weight of
> eq. (6.15).
>
> **What still holds:** the `refValue`/`refGradient`/`valueFraction` mechanics
> of the mixed condition, the `includeDriftFlux` rationale, the SEE summation
> `Gamma_SEE = SUM gamma_i Gamma_i,wall`, and the neutral/ion/electron split.
>
> **The current reference** is the "What the wall-flux conditions actually
> solve" section of the top-level [`README.md`](../../../../../README.md), with
> the derivations in the class headers and the verification in
> `src/applications/utilities/testWallFlux`.

# Solid Surface Flux (drift-diffusion)

The ``ddSolidSurfaceFlux` boundary condition implements a physically consistent wall boundary condition for species transport. It accounts for thermal collection, electric-field-driven drift, and secondary electron emission.

It is designed as a **base class framework (abstract class)** that provides the common numerical logic for mapping physical flux balance to OpenFOAM’s `mixed` boundary condition. This base class is specialized into three **derived classes**:

1.  **`electronDDSolidSurfaceFlux`**: Specifically handles high-mobility electrons and Secondary Electron Emission (SEE).
2.  **`ionDDSolidSurfaceFlux`**: Handles positive and negative ions influenced by the electric field.
3.  **`neutralDDSolidSurfaceFlux`**: Handles non-charged species driven purely by random thermal motion.

---

## 1. Theory

At a solid interface (such as a dielectric wall or electrode), the net species flux ($\mathbf{\Gamma}$) leaving the plasma volume must equal the flux absorbed or emitted by the wall surface ($\Gamma_{wall}$).

### A. Flux Balance Equation
In the Drift-Diffusion approximation, the net flux at the boundary is the sum of the directed drift (driven by the electric field) and the random diffusion (driven by the density gradient). This sum must balance the collection flux at the wall:

$$\Gamma_{plasma} = (Z \mu \mathbf{E})n - D \nabla n = \Gamma_{wall}$$

Where:
* $Z$: Species charge number (e.g., -1 for electrons, +1 for ions).
* $\mu$: Mobility [$m^2/(V \cdot s)$].
* $D$: Diffusion coefficient [$m^2/s$].
* $\Gamma_{wall}$: The physical flux collected by the surface, defined as $u_{wall} \cdot n$.

### B. Wall Velocity ($u_{wall}$)
The "collection velocity" ($u_{wall}$) represents how effectively the wall surface "sinks" the incoming particles. 

#### 1. Neutrals
Neutrals are unaffected by electric fields ($Z=0$). Their flux is determined solely by the random thermal velocity ($v_{th}$):
$$u_{wall, neutral} = \frac{1}{4} v_{th} = \frac{1}{4} \sqrt{\frac{8 k_B T}{\pi m}}$$

#### 2. Ions (Positive and Negative)
Ion collection is enhanced or inhibited by the electric field in the sheath. To maintain consistency with the flux balance derivation, the wall velocity is:
$$u_{wall, ion} = u_{th} + \max(0.0, Z \cdot (\mu \mathbf{E} \cdot \mathbf{n_s}))$$
* **Field pushes ions INTO the wall ($Z \mu E_n > 0$):** The wall velocity increases to $u_{th} + u_{drift}$. This effectively removes the potential barrier, allowing the electric field to drive the ions into the wall.
* **Field pulls ions AWAY from the wall ($Z \mu E_n < 0$):** The $\max$ term becomes zero. $u_{wall}$ remains at the thermal speed ($u_{th}$), representing the random particles that manage to reach the wall despite the repelling field.

#### 3. Electrons
Electrons follow a similar collection logic to ions but are significantly more mobile and susceptible to surface-interaction effects.

* **Drift Flux Switch (`includeDriftFlux`):** While electron collection is often dominated by thermal motion ($u_{th}$), the drift component ($Z \mu E_n$) can be significant in regions with strong electric fields (e.g., near anodes or in high-pressure discharges). This BC provides a toggle to include or exclude the directed drift component in the wall velocity calculation to suit the specific physical regime being modeled.
* **Secondary Electron Emission (SEE):** A critical feature for electrons is the emission flux ($\Gamma_{SEE}$) triggered by ion bombardment. When ions strike the wall, they can knock off secondary electrons, creating a source flux that opposes the primary collection flux:
    $$\Gamma_{e, wall} = u_{wall} n_e - \Gamma_{SEE}$$
    The wall velocity $u_{wall}$ is calculated using the same logic as ions to account for the potential barrier or attraction, while the SEE flux is handled as an independent source term injected via the gradient ($refGrad$).

### C. Secondary Electron Emission (SEE) Logic

In many plasma-surface interactions, the primary source of electrons at a dielectric or cathode wall is **Secondary Electron Emission**. This occurs when heavy positive ions strike the surface with enough energy to overcome the material's work function, ejecting one or more electrons back into the plasma.

#### 1. The Multi-Species Source Term
The total secondary flux ($\Gamma_{SEE}$) is calculated as the sum of all individual ion fluxes reaching the wall, each weighted by its specific secondary emission coefficient ($\gamma_i$):

$$\Gamma_{SEE} = \sum_{i \in \text{ions}^+} \gamma_i \cdot \Gamma_{i, \text{wall}}$$

Where:
* $\gamma_i$: The secondary electron emission coefficient (SEEC) for ion species $i$.
* $\Gamma_{i, \text{wall}}$: The net flux of ion species $i$ into the wall surface.

---

## 2. Numerical Formulation

The `ddSolidSurfaceFlux` suite leverages OpenFOAM’s native `mixed` boundary condition framework to enforce the physical flux balance. This approach provides a stable transition between fixed-value and fixed-gradient behaviors depending on the local plasma conditions.

### A. The OpenFOAM `mixed` BC Mechanism
In OpenFOAM, a `mixed` boundary condition calculates the patch face value ($n_f$) as a weighted linear interpolation:

$$n_f = f \cdot n_{ref} + (1 - f) \cdot \left( n_c + \frac{g}{\Delta} \right)$$

Where:
* $f$ is the `valueFraction` (ranges from 0 to 1).
* $n_{ref}$ is the `refValue` (target Dirichlet value).
* $g$ is the `refGrad` (target Neumann gradient).
* $n_c$ is the internal cell-center value.
* $\Delta$ is the inverse distance between the cell center and the patch face center ($1/d$).

### B. Derivation

To implement the Drift-Diffusion flux at the wall, we begin with the physical balance for a generic particle with charge $Z$:

$$\Gamma_{plasma} = (Z \mu E_n) n_f - D \nabla n = u_{wall} n_f$$

Using a first-order discretization for the gradient at the wall, $\nabla n \approx \frac{n_f - n_p}{d}$, we substitute this into the balance equation:

$$(Z \mu E_n) n_f - D \left( \frac{n_f - n_p}{d} \right) = u_{wall} n_f$$

Rearranging the terms to isolate the wall density $n_f$:

$$n_f \left( \frac{D}{d} + u_{wall} - Z \mu E_n \right) = \frac{D}{d} n_p$$

Dividing by the leading coefficient:

$$n_f = \left( \frac{\frac{D}{d}}{\frac{D}{d} + u_{wall} - Z \mu E_n} \right) n_p$$

By comparing this result to the simplified `mixed` equation (assuming $n_{ref} = 0$ and $g = 0$ for a standard sink), we identify the **Value Fraction ($f$)**:

$$1 - f = \frac{\frac{D}{d}}{\frac{D}{d} + u_{wall} - Z \mu E_n} \implies$$

$$\boxed{f = \frac{u_{wall} - Z \mu E_n}{\frac{D}{d} + u_{wall} - Z \mu E_n}}$$

### C. Species-Specific Formulations

The following sections detail how the generic coefficients are mapped to specific plasma species.

#### 1. Neutrals (`neutralDDSolidSurfaceFlux`)
Neutrals carry no charge ($Z=0$), meaning they are unaffected by the electric field. Their wall interaction is governed entirely by random thermal motion.

* **Wall Velocity ($u_{wall}$):** Set to the one-quarter thermal velocity to represent the physical one-sided flux of a Maxwellian gas hitting a surface:
    $$u_{wall, neutral} = \frac{1}{4} v_{th} = \frac{1}{4} \sqrt{\frac{8 k_B T}{\pi m}}$$
* **Value Fraction ($f$):** Simplifies to the ratio of surface collection speed to total transport speed:
    $$f = \frac{\frac{1}{4} v_{th}}{D/d + \frac{1}{4} v_{th}}$$
* **Reference Gradient ($g$):** $0$.
* **Reference Value ($n_{ref}$):** $0$.

#### 2. Ions (`ionDDSolidSurfaceFlux`)
Ions are influenced by the electric field ($Z \neq 0$). The boundary condition must account for whether the field is pushing ions into the wall or repelling them. To maintain consistency with the flux balance derivation ($Z \mu E_n n - D \nabla n = u_{wall} n$), the formulation is:

* **Wall Velocity ($u_{wall}$):** The sum of the random thermal flux and the directed drift (if pointing toward the wall):
    $$u_{wall, ion} = \frac{1}{4} v_{th} + \max(0.0, Z \mu E_n)$$
* **Value Fraction ($f$):** Accounts for the net collection speed relative to the diffusion rate:
    $$f = \frac{\frac{1}{4} v_{th} + \max(0.0, Z \mu E_n) - Z \mu E_n}{D/d + \frac{1}{4} v_{th} + \max(0.0, Z \mu E_n) - Z \mu E_n}$$
    * **If field pushes ions IN ($Z \mu E_n > 0$):** The drift terms  cancel ($u_{th} + u_D - u_D$), leaving $f = \frac{u_{th}}{D/d + u_{th}}$. This allows the field to drive ions out of the domain without a numerical barrier.
    * **If field pulls ions AWAY ($Z \mu E_n < 0$):** $u_{wall}$ stays at $u_{th}$, and the negative drift increases the "suction" to ensure thermal particles are still collected.
* **Reference Gradient ($g$):** $0$.
* **Reference Value ($n_{ref}$):** $0$.

#### 3. Electrons (`electronDDSolidSurfaceFlux`)
Electrons ($Z = -1$) are highly mobile and follow a similar collection logic to ions, but with the critical addition of **Secondary Electron Emission (SEE)** and a toggle for drift flux inclusion.

* **Wall Velocity ($u_{wall}$):** Calculated to account for thermal collection and, optionally, directed drift when the field opposes electron escape:
    $$u_{wall, e} = \frac{1}{4} v_{th} + \max(0.0, - \mu E_n)$$

* **Value Fraction ($f$):** Defines the balance between the collection speed and the diffusion rate at the wall:
    $$f = \frac{\frac{1}{4} v_{th} + \max(0.0, - \mu E_n) + \mu E_n}{D/d + \frac{1}{4} v_{th} + \max(0.0, - \mu E_n) + \mu E_n}$$
    * When `includeDriftFlux` is enabled, this ensures the effective "suction" at the wall accurately represents the net physical flux.

* **refValue and refGrad:**
    When Secondary Electron Emission (SEE) is enabled, the flux balance must account for the additional source of electrons at the wall. The previous expression for the patch density $n_f$ is modified by the emission flux $\Gamma_{SEE}$:
    
    $$n_f = \left( \frac{\frac{D}{d}}{\frac{D}{d} + u_{wall} - Z \mu E_n} \right) n_p + \frac{\Gamma_{SEE}}{\frac{D}{d} + u_{wall} - Z \mu E_n}$$
    
    

    This extra term represents the contribution of secondary electrons to the surface density. Within the framework of the OpenFOAM **mixed boundary condition**:
    
    $$n_f = f \cdot refValue + (1 - f) \cdot \left( n_p + \frac{refGrad}{\Delta} \right)$$
    
    We can satisfy the physical flux balance by choosing to account for the emission via either the `refGrad` or the `refValue`.

    **Method 1: Reference Gradient (`refGrad`)**
    By setting $refValue = 0$, the wall is treated as a sink for plasma electrons, and the emission term is equated to the gradient $g$:
    
    $$g = refGrad = \frac{\Gamma_{SEE}}{D_e}$$
    
    In this configuration, the boundary condition "injects" a specific slope at the wall. This is the standard convention for Drift-Diffusion solvers as it treats the emission as a pure flux source.

    **Method 2: Reference Value (`refValue`)**
    By setting $refGrad = 0$, the emission term is equated to the reference value $refValue$:
    
    $$refValue = \frac{\Gamma_{SEE}}{u_{wall} - Z \mu E_n} = \frac{\Gamma_{SEE}}{\frac{1}{4} v_{th} + \max(0.0, - \mu E_n) + \mu E_n}$$
    
    

    This method effectively "lifts" the wall density to a non-zero value that is physically consistent with the required emission rate.

---

## 3. Usage

The `ddSolidSurfaceFlux` suite is configured within the species density files located in the `0/` directory. Each species type requires a specific derived boundary condition and its corresponding set of keywords.

### A. Electron Configuration (`0/n_e`)
The electron boundary condition includes additional logic for Secondary Electron Emission (SEE) and a toggle for drift-flux inclusion in the wall velocity.

```cpp
// 0/<regionName>/n_e
plasmaRegion_to_dielectricRegion
{
    type            electronDDSolidSurfaceFlux;
    
    T               20000;         // Temperature [K] (can be a constant or the temperature field name, e.g., T_e)

    // --- Physics Switches ---
    enableSurfaceCharging   true;  // Accumulate species flux to surface charge field
    includeDriftFlux        true;  // Include Z*mu*En in wall velocity calculation
    
    // --- Secondary Electron Emission (SEE) ---
    enableSEE       true;          // Toggle for emission logic
    setRefValue     false;         // true: use Method 2 (refValue)
                                   // false: use Method 1 (refGrad) [Standard]
    
    defaultSEEC     0.1;           // Default emission yield for all positive ions
    speciesSEEC                    // Species-specific yield overrides
    {
        Ar+         0.05;          // Yield specifically for Ar+
        He+         0.15;          // Yield specifically for He+
    }

    value           uniform 0;     // Initial guess/placeholder
}
```

### B. Ion Configuration (`0/n_pIon`)
The ion boundary condition is designed to handle charged species (positive or negative). It focuses on the field-driven collection at the wall and ensuring the accumulation of charge is correctly passed to the potential solver via `enableSurfaceCharging`.

```cpp
// 0/<regionName>/n_pIon
plasmaRegion_to_dielectricRegion
{
    type            ionDDSolidSurfaceFlux;
    
    T               300;           // Temperature [K] (can be a constant or the temperature field name, e.g., T_pIon)

    // --- Physics Switches ---
    enableSurfaceCharging   true;  // Accumulate species flux to the surface charge field

    value           uniform 0;     // Initial guess/placeholder
}
```

### C. Neutral Configuration (`0/n_neutral`)
Neutrals are governed solely by random thermal motion. Because they carry no charge ($Z = 0$), they are unaffected by electric fields and do not contribute to surface charging. The logic defaults to the physical $1/4 v_{th}$ collection flux.

``` cpp
// 0/<regionName>/n_neutral
plasmaRegion_to_dielectricRegion
{
    type            neutralDDSolidSurfaceFlux;
    
    T               300;           // Temperature [K] (can be a constant or the temperature field name, e.g., T_neutral)

    // Neutrals ignore enableSurfaceCharging and includeDriftFlux as they are non-charged

    value           uniform 0;     // Initial guess/placeholder
}

---

# The flux scheme and the wall closure — why they are not independent

Added 2026-09-06.

## How the mixed condition actually imposes a flux

The condition does **not** set a flux. It sets a face VALUE, and lets the
discretisation's own face-flux formula turn that into a flux. With
`refValue = 0` and `refGradient = 0`, OpenFOAM's `mixed` gives

    n_p = (1 - f) n_c

so `f` is the only lever, and it must be chosen so that the flux the
discretisation extracts equals the closure's `n_p * W` — Hagelaar eq. (6.1).

**Therefore `f` is scheme-dependent.** The two flux schemes extract the
boundary flux with different formulae, so they need different `f`, and each is
exact for its own. This is not a defect; it is the condition doing its job
against two different discretisations.

## Standard: `f = uEff/(D/δ + uEff)`

The boundary flux is diffusion plus drift:

    Γ = (D/δ)(n_c - n_p) + uDrift·n_p

and the choice of `f` gives `(D/δ)(n_c - n_p) = n_p·uEff`, hence

    Γ = n_p(uEff + uDrift) = n_p·W                          [eq. (6.1)]

This is an identity, not a fit. It is why `uEff` is defined as `W - uDrift` in
the first place — see `calcEffectiveWallVelocity()`.

## Scharfetter–Gummel: `f = uEff/((D/δ)·Bern(Pe) + W)`

`fvm::ScharfetterGummel` builds boundary coefficients
`pCoeffP = (D/δ)·Bern(-Pe)` and `pCoeffB = (D/δ)·Bern(Pe)` and combines them
with this condition's `valueInternalCoeffs` (see `ScharfetterGummel.H`, the
physical-boundary branch), so the boundary flux is

    Γ = (D/δ)[Bern(-Pe)·n_c - Bern(Pe)·n_p]

Requiring `Γ = n_p·W` with `n_p = (1-f)n_c`:

    (D/δ)Bern(-Pe) = (1-f)[W + (D/δ)Bern(Pe)]

and using the Bernoulli identity **`Bern(-x) - Bern(x) = x`**, which follows
from `Bern(x) = x/(e^x - 1)` in one line:

    f = [W + (D/δ)Bern(Pe) - (D/δ)Bern(-Pe)] / [W + (D/δ)Bern(Pe)]
      = uEff / [(D/δ)Bern(Pe) + W]

and **`uAbs` IS `W`** in every one of the four conditions —
`calcAbsorptionVelocity()` returns the total loss speed, and for ions
`v_th + max(0,u_d)` equals `u_d + uEff` on both drift branches — so the shipped
`D_delta*Bern(Pe) + uAbs` is exactly this.

## Both reduce to each other as Pe → 0

`Bern(0) = 1` and `uDrift → 0`, so `W → uEff` and both give
`uEff/(D/δ + uEff)`. This is the check that catches a wrong SG form, and it is
now in `testWallFlux`.

## WITHDRAWN: the "SG imposes 50x the intended flux" claim

**Earlier on 2026-09-06 this document asserted that the SG branch imposed the
wrong flux — 1.23x at Pe = 1, 5.50x at Pe = 10, 50.5x at Pe = 100, −5.19x at
r = 0.36 — and the combination was briefly made FATAL. That was wrong, and the
claim is withdrawn.** The text is kept here rather than deleted so that stale
copies of it elsewhere can be recognised.

The error was in the diagnostic, not the code: I computed SG's `n_p` and then
pushed it through the **standard** extraction formula. That compares each
branch's `f` against the other branch's discretisation, which means nothing.
The numbers did not even reproduce on recheck — the same wrong quantity gives
1.00 / 1.67 / 1.96, not 1.23 / 5.50 / 50.5 — so the original measurement was
wrong in a way that could not be reconstructed, and per rule 22 a number from
an unreconcilable diagnostic is evidence for nothing.

Two lessons worth more than the retraction:

1. **`fluxScheme` is the INTERIOR scheme** (`driftDiffusion.C:161`), read by
   the BC only to pick its own `f`. Making it fatal punished a legitimate and
   valuable interior discretisation for a supposed defect in the boundary
   treatment. Even had the BC form genuinely been wrong, the right response was
   to fix or collapse the BC — never to refuse the scheme.
2. **Compare each branch against ground truth, not against each other.** The
   two `f` are *supposed* to differ; the invariant they share is
   `Γ = n_p·W`, each under its own extraction.

## Consequence for the singularity

`ddWallFluxMixed` reports a singularity when `D/δ + uEff <= 0`, i.e. when the
drift into the wall exceeds what the thermal term and near-wall diffusion can
carry. With reflection `r` the threshold is

    uDrift > (1-r)/(2r) · A          (only 0.889·A at r = 0.36)

so a reflecting wall reaches it far sooner than a non-reflecting one. The
remedies are physical: resolve the near-wall cell (`D/δ` grows as the cell
shrinks), or accept what it represents — a reflecting wall genuinely piles
density up, and past a point the mixed form cannot represent the wall value
that implies.

**Switching to `ScharfetterGummel` changes the `f` and so can move that
threshold, but it is a change of discretisation, not a fix for
under-resolution.** If the standard branch is singular at your electrode, the
mesh is telling you something.

## What is verified, and what is not

`testWallFlux` verifies:

* the CLOSURE ALGEBRA — eqs. (6.1), (6.2), (6.8), (6.15), the `(1-r)` factors,
  the `Dd` floor, degenerate inputs — to 1.7e-14 over 2430 `(r, Dd, Gc)` points;
* **the `f`-IDENTITY for BOTH schemes** (added 2026-09-06): each branch imposes
  `Γ = n_p·W` under its own extraction formula, over 36 `(Pe, r)` points with
  Pe from 0 to 300 and r from 0 to 0.8, with liveness controls that convict
  each `f` when used with the *other* branch's extraction, plus the Bernoulli
  identity and the `Pe → 0` coincidence.

It does **not** exercise the assembled condition inside a running
discretisation — `valueInternalCoeffs` is reproduced analytically here, not
called. That remains the gap.
