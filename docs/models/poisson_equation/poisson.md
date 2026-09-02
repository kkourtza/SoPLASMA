# Poisson Equation

The Poisson equation is used to calculate the self-consistent electric potential $\phi$ and the resulting electric field $\mathbf{E}$ based on the distribution of charged species (electrons and ions) and fixed charges (dielectrics).

---

## 1. Theory

The electrostatic potential $\phi$ is governed by the Poisson equation:

$$\nabla \cdot ( \epsilon_r \epsilon_0 \nabla \phi ) = -\rho$$

Where:
* $\phi$: Electric potential [$V$]
* $\epsilon_0$: Vacuum permittivity ($\approx 8.854 \times 10^{-12}$ F/m)
* $\epsilon_r$: Relative permittivity (dielectric constant) of the medium
* $\rho$: Net space charge density [$C/m^3$]



### A. Space Charge Density ($\rho$)
The net charge density is calculated by summing the contributions of all charged species $i$:

$$\rho = e \sum_{i} Z_i n_i$$

Where $e$ is the elementary charge, $Z_i$ is the charge state of species $i$, and $n_i$ is the number density.

### B. Electric Field Calculation
Once the potential field is solved, the electric field $\mathbf{E}$ is obtained via the gradient:

$$\mathbf{E} = -\nabla \phi$$


### C. Semi-Implicit Treatment

The equation provided above represents the **explicit form** of the Poisson equation. In this formulation, the solution for the potential at the current time-step, $\phi^{k+1}$, is calculated using the space charge density obtained from the previous time-step, $\rho^k$. 

This explicit treatment imposes a strict constraint on the time-step $\Delta t$, governed by the **dielectric relaxation time** $\tau_d$. The stability of the solution depends on the ratio:

$$\Lambda = \frac{\Delta t}{\tau_d}$$

where $\tau_d$ is the dielectric relaxation time, defined as:

$$\tau_d = \frac{\epsilon}{\sigma}$$

Here, $\epsilon$ is the permittivity of the medium and $\sigma$ is the electrical conductivity, calculated as the sum over all charged species:

$$\sigma = \sum_{i} |q_i| \mu_i n_i$$

where $q_i$ is the charge, $\mu_i$ the mobility, and $n_i$ the number density of species $i$.

For a stable numerical solution, this ratio must satisfy $\Lambda \leq 1.0$. Physically, if the time-step exceeds the dielectric relaxation time ($\Delta t > \tau_d$), the numerical scheme fails to resolve the rapid physical response of the medium to charge imbalances. This leads to **unphysical oscillations**: the electric field overcorrects the motion of charged particles within a single time-step, causing charge densities to swing between extreme positive and negative values rather than relaxing toward equilibrium.

#### Semi-Implicit Formulation
To bypass this restrictive time-step limitation, a **semi-implicit** form is employed. This approach is derived from a first-order Taylor expansion of the space charge density $\rho$ at time $t^{k+1}$ around the current time $t^k$:

$$\rho^{k+1} \approx \rho^k + \Delta t \left( \frac{\partial \rho}{\partial t} \right)^k$$

In a fully consistent scheme, reaction terms would be included in this expansion. However, here we account only for **drift** and **diffusion**, assuming that the sum of the reaction terms for all species equals zero. By substituting the continuity equation ($\frac{\partial \rho}{\partial t} = -\nabla \cdot \mathbf{J}$) and linearizing the drift flux ($\mathbf{J}_{drift} = -\sigma \nabla \phi$), the semi-implicit Poisson equation becomes:

$$\nabla \cdot \left[ \left( \epsilon + \Delta t \sigma^k \right) \nabla \phi^{k+1} \right] = -\rho^k - \Delta t \sum_{i} q_i \nabla \cdot (D_i \nabla n_i^k)$$

By incorporating the conductivity term into the divergence operator, the solver "anticipates" the change in charge distribution due to the electric field at the next time-step, allowing for stable solutions even when $\Delta t$ exceeds the dielectric relaxation time.

---

## 2. Explicit vs Implicit coupling of regions

Most solvers developed in the current toolkit allow for an **implicit coupling** of different regions. In a standard segregated approach, each region's Poisson or Laplace equation is solved separately, and the solver iterates to converge the interface values, which is  computationally expensive.

### Explicit Coupling Formulation

Under explicit coupling, the system matrices contain only intra-region coefficients, while the interface interaction appears only in the right-hand sides through explicit source terms. 

Consider two regions, $A^1$ and $A^2$, sharing a common interface. Region $A^1$ consists of cells $1 \dots N_1$, with the interface at cell $N_1$, and region $A^2$ consists of cells $1 \dots N_2$, with the interface at cell 1. The discretized systems can be written as:

$$\mathbf{A}^1 \phi_{A^1} = \mathbf{b}_{A^1}, \quad \mathbf{A}^2 \phi_{A^2} = \mathbf{b}_{A^2}$$

#### Region $A^1$ system (interface at cell $N_1$)

$$
\begin{bmatrix}
a^1_{1,1} & a^1_{1,2} & \cdots & a^1_{1,N_1-1} & a^1_{1,N_1} \\
a^1_{2,1} & a^1_{2,2} & \cdots & a^1_{2,N_1-1} & a^1_{2,N_1} \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
a^1_{N_1-k-1,1} & a^1_{N_1-k-1,2} & \cdots & a^1_{N_1-k-1,N_1-1} & a^1_{N_1-k-1,N_1} \\
a^1_{N_1-k,1} & a^1_{N_1-k,2} & \cdots & a^1_{N_1-k,N_1-1} & a^1_{N_1-k,N_1} \\
a^1_{N_1-k+1,1} & a^1_{N_1-k+1,2} & \cdots & a^1_{N_1-k+1,N_1-1} & a^1_{N_1-k+1,N_1} \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
a^1_{N_1-1,1} & a^1_{N_1-1,2} & \cdots & a^1_{N_1-1,N_1-1} & a^1_{N_1-1,N_1} \\
a^1_{N_1,1} & a^1_{N_1,2} & \cdots & a^1_{N_1,N_1-1} & a^1_{N_1,N_1}
\end{bmatrix}
\begin{bmatrix}
\phi_{A^1,1} \\
\phi_{A^1,2} \\
\vdots \\
\phi_{A^1,N_1-k-1} \\
\phi_{A^1,N_1-k} \\
\phi_{A^1,N_1-k+1} \\
\vdots \\
\phi_{A^1,N_1-1} \\
\phi_{A^1,N_1}
\end{bmatrix}
=
\begin{bmatrix}
b_{A^1,1} \\
b_{A^1,2} \\
\vdots \\
b_{A^1,N_1-k-1} \\
b_{A^1,N_1-k}^{int} + S_{A^1,N_1-k}^{(A^2)} \\
b_{A^1,N_1-k+1}^{int} + S_{A^1,N_1-k+1}^{(A^2)} \\
\vdots \\
b_{A^1,N_1-1}^{int} + S_{A^1,N_1-1}^{(A^2)} \\
b_{A^1,N_1}^{int} + S_{A^1,N_1}^{(A^2)}
\end{bmatrix}
$$

The source term $S_{A^1,i}$ depends only on the previous-iteration values of the cells of region $A^2$ that face the interface:
$$S_{A^1,i} = \alpha_i^1 \phi_{A^2, \text{nb}(i)}^{\text{old}}, \quad i = N_1 - k, \dots, N_1$$
where $\text{nb}(i)$ denotes the neighbouring cell of region $A^2$ that is geometrically connected to cell $i$ of region $A^1$ across the interface.


#### Region $A^2$ system (interface at cells $1, \dots, k$)

$$
\begin{bmatrix}
a^2_{1,1} & a^2_{1,2} & \cdots & a^2_{1,N_2-1} & a^2_{1,N_2} \\
a^2_{2,1} & a^2_{2,2} & \cdots & a^2_{2,N_2-1} & a^2_{2,N_2} \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
a^2_{k-1,1} & a^2_{k-1,2} & \cdots & a^2_{k-1,N_2-1} & a^2_{k-1,N_2} \\
a^2_{k,1} & a^2_{k,2} & \cdots & a^2_{k,N_2-1} & a^2_{k,N_2} \\
a^2_{k+1,1} & a^2_{k+1,2} & \cdots & a^2_{k+1,N_2-1} & a^2_{k+1,N_2} \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
a^2_{N_2-1,1} & a^2_{N_2-1,2} & \cdots & a^2_{N_2-1,N_2-1} & a^2_{N_2-1,N_2} \\
a^2_{N_2,1} & a^2_{N_2,2} & \cdots & a^2_{N_2,N_2-1} & a^2_{N_2,N_2}
\end{bmatrix}
\begin{bmatrix}
\phi_{A^2,1} \\
\phi_{A^2,2} \\
\vdots \\
\phi_{A^2,k} \\
\phi_{A^2,k+1} \\
\vdots \\
\phi_{A^2,N_2-1} \\
\phi_{A^2,N_2}
\end{bmatrix}
=
\begin{bmatrix}
b_{A^2,1}^{int} + S_{A^2,1}^{(A^1)} \\
b_{A^2,2}^{int} + S_{A^2,2}^{(A^1)} \\
\vdots \\
b_{A^2,k}^{int} + S_{A^2,k}^{(A^1)} \\
b_{A^2,k+1} \\
\vdots \\
b_{A^2,N_2-1} \\
b_{A^2,N_2}
\end{bmatrix}
$$

$$S_{A^2,j} = \alpha_j^2 \phi_{A^1, \text{nb}(j)}^{\text{old}}, \quad j = 1, \dots, k$$

---

### Implicit Coupling Formulation

The interface coupling terms are now placed on the left-hand side, forming off-diagonal blocks that connect the degrees of freedom of $A^1$ and $A^2$. We implemented this formulation for the Poisson equation governing the electric potential. The idea is illustrated by the block matrix system:

$$
\begin{bmatrix}
A^1 & C_{12} \\
C_{21} & A^2
\end{bmatrix}
\begin{bmatrix}
\phi_{A^1} \\
\phi_{A^2}
\end{bmatrix}
=
\begin{bmatrix}
b^1 \\
b^2
\end{bmatrix}
$$

Here $C_{12}$ and $C_{21}$ contain the implicit interface contributions coupling the two regions. Unlike the explicit case, these terms appear directly on the left-hand side, so both regions are solved simultaneously in a single linear solve.

#### How to Use the Implicit Coupling Approach

In order to use the implicit coupling approach for the Poisson equation between the two domains, the user should enable the `useImplicit` switch that appears in the boundary conditions of the electric potential field (`ePotential`) across the coupled domains for the interface patch.

```cpp
// 0/plasmaRegion/ePotential
ePotential
{
    internalField   uniform 0;

    boundaryField
    {
        plasmaRegion_to_dielectricRegion
        {
            type            coupledElectricPotential;
            value           $internalField;
            phiNbr          ePotential;
            surfCharge      surfCharge;
            surfChargeNbr   none;
            useImplicit     true;
        }
    }
}
```

```cpp
// 0/dielectricRegion/ePotential
ePotential
{
    internalField   uniform 0;

    boundaryField
    {
        dielectricRegion_to_plasmaRegion
        {
            type            coupledElectricPotential;
            value           $internalField;
            phiNbr          ePotential;
            surfCharge      none;
            surfChargeNbr   surfCharge;
            useImplicit     true;
        }
    }
}
```

#### Important Notes on Monolithic Coupling

A. **Parallel Decomposition Constraints:**
   For implicit coupling to work in parallel, the coupled patch-pairs must reside on the same processor (this applies to the neighboring cells across the interface). Additionally, they must have a matching number of faces. To ensure this during the decomposition process, you must use a constraint in your `decomposeParDict`:

```cpp
// system/decomposeParDict
constraints
{
    faces
    {
        type    preserveFaceZones;
        zones   (interfaceZone);
        enabled true;
    }
}
```

B. **PETSc Solver Configuration:**
    If the `PETSc` library is used to solve the electric potential, you must explicitly specify `implicitRegionCoupling true` in the `fvSolution` dictionary:

``` cpp
// system/fvSolution
ePotential
{
    solver                  petsc;
    implicitRegionCoupling  true;

    petsc
    {
        options
        {
            ksp_type  cg;
            pc_type   hypre;
            mat_type  aij;
        }
    }
    tolerance   1e-8;
    relTol      0.0;
}
```

C. **Solver Compatibility:**
    With native OpenFOAM solvers, `GAMG` works with monolithic coupling **only if you also select an assembly-aware agglomerator**:

```cpp
ePotential
{
    solver          GAMG;
    agglomerator    assembledFaceAreaPair;   // NOT the default
    ...
}
```

GAMG's *default* agglomerator (`faceAreaPair`) does `refCast<const fvMesh>` on the matrix's mesh, and the monolithic matrix is an `lduPrimitiveMeshAssembly`, not an `fvMesh`. The run aborts with

```
Attempt to cast type lduPrimitiveMeshAssembly to type fvMesh
```

which reads like an internal error but is a solver-configuration error. OpenFOAM ships `assembledFaceAreaPair` (`TypeName` in `src/finiteVolume/lduPrimitiveMeshAssembly/assemblyFaceAreaPairGAMGAgglomeration`) for exactly this case.

> **Corrected 2026-09-01.** *SUPERSEDED BY this paragraph:* both this file and
> `coupledElectricPotential.md` previously stated that GAMG "will not work with
> monolithic coupling... because the geometric grid agglomeration logic cannot
> handle the discontinuous connectivity". The agglomeration limitation is real,
> but it belongs to the *default* agglomerator only, and it is a configuration
> choice rather than a limitation of the coupling. Non-GAMG solvers
> (`PBiCGStab`, `PCG`) and PETSc need no such entry.

**Measured 2026-09-01**, `plate2D_timeVaryingBC_implicitBoundary`, 20000 cells,
per timestep, all four reaching the same analytic interface potential where they
converge at all:

| solver | iterations/step | `V_interface` (analytic 1/6) |
|---|---|---|
| `smoothSolver`/`GaussSeidel`, `maxIter 2000` | 2000 (cap), final residual 3e-4 | 0.154443 — **7.3% wrong** |
| `PCG`/`DIC` | 114 | 0.166667 ✓ |
| `GAMG`, default agglomerator | **abort** | — |
| `GAMG` + `assembledFaceAreaPair` | **11** | 0.166667 ✓ |

The first row is what that tutorial shipped with, and it is the more instructive
failure: reaching `maxIter` is **not** an error, so the log reported ten
consecutive unconverged solves as normal ones. Gauss-Seidel needs O(N) sweeps to
carry information across a Laplace problem, so raising `maxIter` buys wall time
for the same answer rather than fixing it.


D. **Segregated Mode (useImplicit False):**
    If `useImplicit` is set to `false`, a non-monolithic (segregated) coupling is used, and the interface iterations get their own convergence controls:

```cpp
// system/plasmaSimulationControls
poisson
{
    nonCoupledResidualControl
    {
        tolerance       1e-8;
        maxIter         1000;
    }
}
```

> **Corrected 2026-09-01.** *SUPERSEDED BY the block above:* this was documented
> as `ePotentialControls { nonCoupledResidualControl { ... } }` in
> `system/fvSolution`. **Nothing has ever read that.** `ePotentialControls`
> appears in no source file and in no shipped case; the entries were read from
> `<model>Coeffs` in `constant/electromagneticsProperties`, and now from the
> `poisson` block in `system/plasmaSimulationControls`. A user following the old
> text set a tolerance and an iteration cap that had no effect, and got the
> defaults (`1e-6`, `100`) instead.

---

## 3. Numerical Implementation & Controls

### A. Solver Mode Configuration

The global Poisson numerics live in the `poisson` block of
`system/plasmaSimulationControls`. They are uniform over the gas and every
dielectric, because one model owns the whole coupled solve:

``` cpp
// system/plasmaSimulationControls
poisson
{
    scheme                      semiImplicit;   // explicit | semiImplicit
    EScheme                     reconstruct;    // grad | reconstruct
    nNonOrthogonalCorrectors    0;
}
```

`scheme` is a **discretisation scheme**, not a linear solver: the linear solver
for `ePotential` is the `ePotential` block in `fvSolution`. Its default is
`semiImplicit`.

> **Changed 2026-09-01.** *SUPERSEDED BY the block above:* these settings were
> previously `PoissonScheme`, `EScheme` and `nNonOrthogonalCorrectors` inside
> `<model>Coeffs` in `constant/electromagneticsProperties`, and the case variable
> feeding `PoissonScheme` was called `poissonSolver` — two names for one
> setting, each suggesting the other's meaning. `constant/` is for *physical
> properties*, so a discretisation scheme did not belong there; permittivity,
> which is physical, moved the other way, into
> `constant/<region>/electricalProperties`. `constant/electromagneticsProperties`
> no longer exists. An unmigrated case still runs: the old dictionary is read as
> a fallback and prints where each setting has moved to.

### B. Numerical Schemes (`fvSchemes`)

Standard setups typically use second-order linear interpolation for the Laplacian and gradient terms:

``` cpp
// system/fvSchemes
laplacianSchemes
{
    default         Gauss linear corrected; // 'corrected' accounts for non-orthogonality
}
```

### C. Linear System Solvers (`fvSolution`)

The configuration of the fvSolution dictionary follows standard OpenFOAM conventions, with specific requirements depending on whether you are using Implicit (Monolithic) or Explicit (Segregated) region coupling. A typical `fvSolution` file looks like:

``` cpp
// system/fvSolution
solvers
{
    ePotential
    {
        solver                     petsc;
        implicitRegionCoupling     true; // use false for explicit region coupling

        petsc
        {
            options
            {
                ksp_type  cg;
                pc_type  hypre;
                mat_type  aij;
            }
            caching
            {
                matrix
                {
                    update always;  // you can use never for the explicit Poisson
                }

                preconditioner
                {
                    update always;
                }
            }
        }
        tolerance       1e-8;
        relTol          0.0;
    }
    ePotentialFinal
    {
        $ePotential;
        relTol          0;
    }
}

// nonCoupledResidualControl does NOT live here -- see section 2D. It is read
// from the `poisson` block in system/plasmaSimulationControls.

PIMPLE
{
    nOuterCorrectors    3;
    residualControl
    {
        ePotential
        {
            tolerance   1e-6;
            relTol      0;
        }
    }
}
```

### D. Commonly used Boundary Conditions

| Type | Description | Common Use Case |
| :--- | :--- | :--- |
| **`coupledElectricPotential`** | Matches $\phi$ and displacement flux $D$ across regions. | Plasma-Dielectric or Dielectric-Dielectric interfaces. |
| **`fixedValue`** | Dirichlet condition; sets a specific voltage. | Grounded walls (0V), DC electrodes etc. |
| **`zeroGradient`** | Neumann condition; zero electric field flux. | Symmetry planes or ideal (non-charging) insulators. |
| **`fixedGradient`** | Sets the normal component of the electric field ($-\nabla \phi \cdot \mathbf{n}$). | Surfaces with a known, fixed surface charge density. |

For a detailed breakdown of the custom boundary conditions, please refer to the specialized documentation.
