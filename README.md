# SoPLASMA

**SoPLASMA** is a collection of solvers, models, boundary conditions, and utilities developed for plasma applications in **OpenFOAM**. It aims to extend the capabilities of OpenFOAM to simulate plasma-related phenomena, including plasma-fluid interactions, plasma chemistry, and  coupled electromagnetic effects.


##	Description

This toolkit provides an evolving framework for implementing plasma physics features within
OpenFOAM. It is designed for research and development purposes, focusing on modularity and compatibility with the OpenFOAM v2412 (OpenCFD Ltd.).  

##	Status

The project is currently in its **initial development phase**. Functionality, structure, and documentation are expected to change as new components are added.

##	License

This project is distributed under the terms of the **GNU General Public License v3.0 (GPLv3)**.  
For full license details, see the [LICENSE](LICENSE) file.

© 2026 Rention Pasolari

This software is not part of OpenFOAM, but it is developed using the OpenFOAM framework and linked against OpenFOAM libraries (v2412).

## Requirements

### Core Dependency
- **OpenFOAM (OpenCFD Ltd. release)**  
  *(Developed and primarily tested with version **v2412**.)*

This toolkit relies on the OpenCFD Ltd. version of OpenFOAM for compilation and compatibility.

### Additional Dependencies

The toolkit includes a number of additional components that may be required depending on what you intend to use:

- **Crucial dependencies** needed for building the core functionality.
- **Optional dependencies** that enable extended features. 
- **Tutorial-specific tools** required only for certain example cases.  
- **External libraries**, such as **petsc4foam** *(see [`docs/getting_started/petsc4foam.md`](docs/getting_started/petsc4foam.md))*

For a complete and recommended overview of all dependencies, see **[`docs/getting_started/dependencies.md`](docs/getting_started/dependencies.md)**.

## Tested Platforms

The toolkit has been tested on the following configurations:

| Ubuntu Version | OpenFOAM Version (OpenCFD Ltd.) |
|----------------|---------------------------------|
| 24.04 LTS      | v2412                           |
| 22.04 LTS      | v2412                           |

Additional Linux distributions and OpenFOAM (OpenCFD Ltd.) versions may work, but are not officially validated at this time.

---

## Installation / Setup

## 1. Install Required Dependencies

Before installing OpenFOAM and the SoPLASMA, make sure all required dependencies are installed. Refer to **[`docs/getting_started/dependencies.md`](docs/getting_started/dependencies.md)** and install all packages marked as *required*.


## 2. Install OpenFOAM v2412 (required)

This toolkit is designed **for OpenFOAM v2412 (OpenCFD Ltd.)**. It has not been tested in other  versions.

If you already have OpenFOAM installed, skip to **Step 2 (Download and build SoPLASMA)**.


### A. Install OpenFOAM prerequisites

```bash
sudo apt-get update
sudo apt-get install build-essential autoconf autotools-dev cmake gawk gnuplot
sudo apt-get install flex libfl-dev libreadline-dev zlib1g-dev openmpi-bin libopenmpi-dev mpi-default-bin mpi-default-dev
sudo apt-get install libgmp-dev libmpfr-dev libmpc-dev
sudo apt-get install libfftw3-dev libscotch-dev libptscotch-dev libboost-system-dev libboost-thread-dev libcgal-dev
```

### B. Get the source code

Choose the directory where you want OpenFOAM v2412 to be installed.

Common system-wide installation paths are:

- `/usr/lib/openfoam/openfoam2412/`
- `/opt/openfoam2412/`

However, both locations require **sudo** privileges. This can be inconvenient when compiling ThirdParty packages or running tests, as some steps may require elevated permissions.

To avoid this, you can install OpenFOAM in your **home directory**, which allows you to compile and work without using `sudo`.

#### Example: Install OpenFOAM v2412 in your home directory

```bash
cd ~
mkdir OpenFOAM
cd OpenFOAM
```

You can download the official OpenFOAM v2412 source code directly from **openfoam.com**:

```bash
wget -O -  http://dl.openfoam.com/source/v2412/OpenFOAM-v2412.tgz | tar xvz
wget -O -  http://dl.openfoam.com/source/v2412/ThirdParty-v2412.tgz | tar xvz
```

If `wget` fails to establish a connection, you can download the archives manually from the official mirrors and extract them into your chosen OpenFOAM directory.

Manual download links:

- [OpenFOAM-v2412.tgz (source code)](https://sourceforge.net/projects/openfoam/files/v2412/OpenFOAM-v2412.tgz/download)  
- [ThirdParty-v2412.tgz (third-party libraries)](https://sourceforge.net/projects/openfoam/files/v2412/ThirdParty-v2412.tgz/download)


After downloading, extract both archives inside your OpenFOAM installation folder:

```bash
tar -xvzf OpenFOAM-v2412.tgz
tar -xvzf ThirdParty-v2412.tgz
```

### C. Set the environment variables

Open the `.bashrc` file in your home directory with your preferred editor (note the leading dot):

```bash
nano ~/.bashrc
# or
gedit ~/.bashrc
```

Add the following line at the end of the file:

```bash
source $HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc
```

Save the file, then reload your shell configuration:

```bash
source ~/.bashrc
```

### D. Install ThirdPart and OpenFOAM

First, install the ThirdParty packages:

```bash
cd ~/OpenFOAM/ThirdParty-v2412
./Allwmake -j
```

After the ThirdParty build completes, install OpenFOAM itself:

```bash
cd ~/OpenFOAM/OpenFOAM-v2412
./Allwmake -j
```

When the compilation finishes, verify that OpenFOAM is correctly installed:

```bash
foamVersion   # should print: OpenFOAM-v2412
```

## 3. Install SoPLASMA

### A. Get the source code

Choose the directory where you want to install **SoPLASMA**, then download the toolkit by cloning the GitHub repository.

#### Example: Install SoPLASMA in your home directory

```bash
cd ~
git clone https://github.com/rpasolari/SoPLASMA.git
cd SoPLASMA
```

Activate the toolkit environment (so OpenFOAM can detect its extensions):
```bash
source etc/bashrc
```

⚠️ **IMPORTANT:** To avoid having to run this manually every time, add it to your `~/.bashrc` as shown below.

Open the `.bashrc` file in your home directory using any text editor (note the leading dot):

```bash
gedit ~/.bashrc
```

Then add the following line at the end of the file (adjust the path to match where you installed the toolkit):

```bash
source $HOME/SoPLASMA/etc/bashrc
```

Save the file, then reload your shell configuration:

```bash
source ~/.bashrc
```

To verify that the environment was successfully loaded, open a **new terminal** and run:

```bash
echo $SoPLASMA   # should print the toolkit folder
```

### B. Build SoPLASMA
From inside the SoPLASMA directory, compile everything:

```bash
./Allwmake
```
If the compilation finishes without errors, the installation is complete.

⚠️ **IMPORTANT:** Some ThirdParty packages included with the SoPLASMA (but not required as prerequisites) may fail to build (for example, **petsc4foam**). To install these optional components, refer to the corresponding instructions in the `docs/` directory.

⚠️ **The plasma chemistry libraries need SoEEDF.** `plasmaBoltzmann` links the
electron Boltzmann solver, so `./Allwmake` fails on it unless `BOLTZMANN_DIR`
points at a SoEEDF checkout. See the section below.


## Mechanism-driven plasma chemistry

A case declares **one mechanism file**. The transported species, their charges and
masses, the electron transport coefficients, the electron-impact rate coefficients
and the heavy chemistry are all derived from it — nothing is re-typed between
codes, and every generated artifact carries a hash of the master mechanism so a
mismatched pair is refused rather than silently used.

| library | role |
|---|---|
| `plasmaBoltzmann` | solves the EEDF **in process at start-up** and writes the coefficient tables; staleness-checked against the mechanism hash |
| `plasmaChemistry` | the whole reaction set as an OpenFOAM `ODESystem`; native and optional Cantera heavy backends |
| `plasmaSpecies` | derives `activeSpecies`, charges and masses from the mechanism |
| `plasmaTransport` | assembles the chemistry source into the species equations |

The chemistry source is returned as production and a **loss coefficient**
(`dn/dt = P − L·n`) so the sink is implicit via `fvm::Sp`: densities stay positive
at any timestep, and the outer Picard iteration is a contraction on the loss term.
`solver adaptive` chooses **per cell** between linearising and stiff
integration, because the stiffness is local — a streamer head is stiff while the
bulk is not.

### Requirements

**[SoEEDF](https://github.com/kkourtza/SoEEDF)** — the two-term electron Boltzmann
solver that produces the EEDF and the mechanism compiler that generates the case
files. Build it first, then point `BOLTZMANN_DIR` at it (it defaults to
`$HOME/Projects/SoEEDF`):

```bash
git clone git@github.com:kkourtza/SoEEDF.git ~/Projects/SoEEDF
cd ~/Projects/SoEEDF && cmake -S . -B build && cmake --build build -j
```

**Cantera is optional.** Without it everything builds and runs; selecting
`chemistryBackend cantera` then fails with a message naming the variable to set.
Both the 3.x and 4.x series are supported.

### Documentation

Installation, a step-by-step first case, modelling guidelines and the full
validation record live in the SoEEDF repository, indexed at
[`docs/INDEX.md`](https://github.com/kkourtza/SoEEDF/blob/master/docs/INDEX.md).
The chemistry–transport coupling — operator splitting, temporal order, the
dictionary reference and the Cantera backend — is in
[`docs/numerics-chemistry-coupling.md`](https://github.com/kkourtza/SoEEDF/blob/master/docs/numerics-chemistry-coupling.md).

### Backward compatibility

Existing cases are unaffected: `reactions electronImpact` + `solver explicitSource` is the previous behaviour
and remains available. The mechanism-driven path is opt-in per case.


## The electron energy model — one switch

All electron-energy configuration is a single top-level key in
`constant/plasmaSpeciesProperties`:

```
electronEnergyModel  LMEA;      // LFA | LMEA   (REQUIRED - no default)
```

**LFA** (local field approximation) reads every electron coefficient and rate at
the local reduced field `E/N`. **LMEA** (local mean energy approximation) solves
an electron energy-density equation and reads them at the local mean energy
instead, which relaxes the assumption that the electron distribution is in
equilibrium with the local field — the assumption that fails in a streamer head
and at boundaries.

Setting `LMEA` is a *complete* configuration. It derives the electron transport
lookup key, the chemistry rate key, and the whole `energyModelCoeffs` block from
the mechanism's own tables. `tutorials/.../positiveStreamer_LMEA_minimal` is that
one line and nothing else, and produces byte-identical fields to
`positiveStreamer_LMEA_fast`, which keeps the explicit 80-line block as a
regression test of the override path.

Gas heating is a **separate, orthogonal** switch — all four combinations work and
neither implies the other:

```
backgroundGas { energy { solve true; T 300; } }
```

Heavy species have no energy key; they follow `backgroundGas/energy`.

An explicit entry that *contradicts* the closure is fatal, not warned. The
contradiction — transport following one variable while the reaction rates follow
another, the "half-LMEA" — is a measured runaway that produces no error of its
own: plausible fields, a clean run, wrong physics.

Full treatment, including the two measured traps the defaults are shaped to
avoid and the invertibility condition for when LMEA is unavailable, in
[`docs/models/energy/lmea.md`](docs/models/energy/lmea.md).

### Migration note (2026-09-01)

The per-species `energyModel` key is superseded.

- On the **electron** it is **rejected** with an error naming the one-line
  replacement. Keeping it readable would leave two live spellings of one
  setting, which is how the two-vocabulary problem arose in the first place.
- On **heavy species it never had a reader** — it filled five index lists nothing
  consumed — and is also rejected, so a case cannot go on believing it had
  configured something.
- **`electronEnergyModel` is REQUIRED and has no default.** Every case must name
  its closure. This is deliberate: LFA and LMEA are different physics that give
  different answers, so neither can be assumed on a user's behalf — the same
  reason the half-LMEA is fatal rather than warned. The "a required entry whose
  only sensible value is the historical one is a migration tax" argument that
  justifies defaulting `energyModelCoeffs` does not apply, because here there
  are two sensible values. A case without the key stops at start-up with an
  error that explains both closures and what each costs.


## Regions and materials — one place each

A case declares its regions **once**, in `constant/regionProperties`, using
OpenFOAM's own `regionProperties` class and format:

```
regions
(
    gas         (gas)
    dielectric  (barrier1 barrier2)
    farField    (air)                 // fictitious: Poisson only, epsilonR 1
);
```

Three kinds, and the kind is a fixed keyword while the names in parentheses are
free:

| kind | solves | `epsilonR` |
|---|---|---|
| `gas` | species, chemistry, electron energy, Poisson | defaults to 1.0 |
| `dielectric` | Poisson only | **required**, no default |
| `farField` | Poisson only | defaults to 1.0 (that *is* the kind) |

`farField` is a **fictitious air region** whose only job is to extend the
*electrostatic* domain without extending the plasma one. Poisson is long-range;
species transport, chemistry and photoionization are not. Verified exactly
transparent: at εᵣ = 1 the coupled interface reduces to continuity of V and
∂V/∂n, and the two-region analytic bed returns the single-medium answer to a
relative error of 0.00e+00. The solver refuses to run if a species boundary
condition would deposit surface charge on such an interface.

Everything else follows from it:

| what | where | notes |
|---|---|---|
| which regions exist, and their kind | `constant/regionProperties` | the single declaration |
| a region's permittivity | `constant/<region>/electricalProperties` → `epsilonR` | one small file per region |
| the Poisson numerics | `system/plasmaSimulationControls` → `poisson { }` | global, one block |
| which Poisson model runs | **derived** | `multiRegionPoisson` if any `dielectric` is listed, else `singleRegionPoisson` |

This is `chtMultiRegionFoam`'s convention, deliberately: `constant/<region>/` is
that region's physical description, one small file per physics, exactly where
that solver keeps `constant/<region>/thermophysicalProperties`. **Adding a
region means adding a directory**, and adding a physics later (a solid heat solve
in a barrier, say) means adding a file next to `electricalProperties` — neither
edits anything global.

The split follows OpenFOAM's: `constant/` holds *physical properties*, `system/`
holds *how they are solved*. Permittivity is physical and is per region; a
discretisation scheme is neither, since one model owns the whole coupled solve
and forming `E` by a different scheme on either side of an interface would make
the coupled flux inconsistent.

**Any number of dielectrics; exactly one gas.** That is not a limit on the
geometry — a stack `electrode|dielectric|gas|dielectric|gas|dielectric|electrode`
is *one* gas region whose cellZone is geometrically disconnected (an OpenFOAM
mesh region may consist of disjoint parts) plus three dielectrics, and every gap
then shares one species set, one chemistry and one EEDF. What is genuinely
unsupported is two gas regions holding *different mixtures*: the species list,
mechanism and EEDF are global.

### Migration note (2026-09-01)

`constant/electromagneticsProperties` **no longer exists**. It held four
different kinds of thing in one file, three of which were in the wrong place:

- `electromagneticsModel` — **derived** from `regionProperties`. It was never a
  choice: both models solve the same equation and only the topology decides
  which applies. Stating one that contradicts the topology is now fatal; stating
  one that agrees prints a notice and continues.
- `dielectricConstant` and a sub-dictionary per region — **moved** to
  `constant/<region>/electricalProperties` as `epsilonR`. Per-region data in a
  global file meant every region's `Allrun` copy of that file held every other
  region's data; those copies were never read.
- `PoissonScheme`, `EScheme`, `nNonOrthogonalCorrectors`,
  `nonCoupledResidualControl` — **moved** to the `poisson` block in
  `system/plasmaSimulationControls`. `PoissonScheme` is now `scheme`, and its
  default changes from `explicit` to **`semiImplicit`**. The case variable
  feeding it was called `poissonSolver` while the entry was called
  `PoissonScheme`: two names for one setting, each suggesting the other's
  meaning. It selects a *scheme*; the linear solver is in `fvSolution`.
- `backgroundDensity` — already removed, as a second definition of the gas
  density that disagreed with `backgroundGas`.

**Unmigrated cases still run.** The old file is read as a fallback and prints,
per setting, where it has moved to. Two documentation defects were corrected in
the same change: `nonCoupledResidualControl` was documented in a
`system/fvSolution` block named `ePotentialControls` that **no reader has ever
existed for**, and `GAMG` was documented as incompatible with monolithic
coupling when it needs only `agglomerator assembledFaceAreaPair`.


## Dielectric surfaces — meshed, or collapsed onto the boundary

A barrier can be represented two ways, and the choice is about cost, not physics:

- **Meshed** as its own region, coupled by `coupledElectricPotential`. Laplace is
  solved inside the barrier, so no thin-layer assumption is made. Declare the
  region in `constant/regionProperties` and its `epsilonR` in
  `constant/<region>/electricalProperties`; the Poisson model is then *derived*,
  not chosen.
- **Collapsed onto the gas boundary** by `thinDielectricPotential`, which
  replaces the layer with its surface capacitance `C = eps0*epsilonR/d`:

  ```
  eps_g dV/dn = sigma - C (V - Vb)
  ```

  Two situations share the one condition. A **free-standing** sheet — open gas or
  vacuum behind it, nothing to be a capacitor to — is `C = 0`, pure Neumann, and
  needs only `surfCharge`. A **conductively backed** sheet, the classical DBD
  stack, additionally takes `thickness` and `backingPotential`. Those two are
  validated as a pair, and `epsilonR` is *rejected* on a free-standing patch
  because it does nothing there.

  Note that `zeroGradient` is this same condition with `sigma` forced to zero —
  so a wall that is meant to trap charge and is left `zeroGradient` will
  accumulate surface charge and then silently discard it.

Validated against analytic ground truth (the series stack, the surface charge and
its sign, both limits, both guards, and restart round-tripping) — see
[`docs/models/poisson_equation/boundary_conditions/thinDielectricPotential.md`](docs/models/poisson_equation/boundary_conditions/thinDielectricPotential.md).


## Electrodes: the three things a conductor can know

```
                        potential          charge
  driven electrode      known (waveform)   whatever the supply gives
  grounded electrode    known (0)          whatever flows to ground
  floatingElectrode     UNKNOWN            KNOWN (conserved)
```

A **floating electrode** is a conductor connected to nothing — a probe, an
isolated pin, a floating guard ring, one disconnected segment of a segmented
electrode. The usual boundary condition is inverted: the charge is known and the
potential is solved for.

Because Poisson is linear in `V`, this closes in **closed form** rather than by
iteration. With `psi` the unit-potential field (1 on the floating conductor, 0
on every known electrode) and `C_self = ∮ε∇psi·n dA` its self-capacitance,

```
V_f = (Q − Q_rho) / C_self
```

exactly, with no relaxation parameter and no convergence criterion. Declare it
with one entry and no potential — supplying a potential would make it a driven
electrode, so it is rejected:

```
myProbe { type floatingElectrodePotential; initialCharge 0; value uniform 0; }
```

Validated against analytic ground truth, including a 2-D antisymmetric induction
case where `V_f = 0` follows from symmetry alone. In a **plasma**, the charge
ledger `Q(t) = Q0 + ∫I_plasma dt` is fed from the species wall fluxes, and the
electrode charges negative as probe theory requires. Because the default Poisson
scheme is semi-implicit, `psi` is rebuilt every step from the operator the solve
actually used — and the run *verifies* that by checking Gauss's law closes after
each correction (measured 2.8e-14 V). `V_f`, `Q` and `I_plasma` are written per
step to `postProcessing/floatingElectrode/floating.csv`. Details and the
measured tests:
[`docs/models/poisson_equation/floating-electrode.md`](docs/models/poisson_equation/floating-electrode.md).


## One semantic file describes every boundary

You say what each surface **is**; every module decides how to treat it. That is
the whole boundary input for the needle-DBD tutorial:

```
// configuration/boundaries
active_electrode { kind drivenElectrode;       waveform table ((0 0) (100e-9 8e3)); }
ground           { kind groundedElectrode; }
air_dielectric   { kind thinDielectricSurface; material pmma; }
out              { kind openBoundary; }
```

`plasmaSetupBoundaries` turns those four declarations into
`0/<region>/{ePotential,surfCharge}` for every region, reading
[`etc/boundaryRoles`](etc/boundaryRoles) for what each **kind** means to each
module. `plasmaSetupBoundaries -listKinds` prints the kinds and their required
parameters. Seven ship today: `drivenElectrode`, `groundedElectrode`,
`floatingElectrode`, `thinDielectricSurface`, `thinDielectricOnElectrode`,
`openBoundary`, `insulatingWall`.

**Three things are derived and a case may not declare them:**

| derived from | what |
|---|---|
| the **mesh** | `empty`, `wedge`, `symmetry`, `cyclic`, `processor` — a mechanical constraint is not a physical description |
| the **topology** | region interfaces. A gas/dielectric pair gets `coupledElectricPotential` with the surface charge owned by the **gas side only**; a gas/`farField` pair gets `zeroGradient` and never charges |
| the **material library** | `epsilonR` and `gammaSEE` behind a named material |

The interface σ-ownership inverts correctly between the two sides of a barrier
without being stated anywhere — exactly one owner, or the charge is
double-counted.

**needleDBD now ships no `etc/` and no `0.orig/`.** Every field in `0/` is
generated, and the two hand-written files are `configuration/config` and
`configuration/boundaries`. Verified 2026-09-04 against the hand-written
pipeline it replaces: `C_g = 9.03816931478e-17 F` and
`I_disp = 7.23053545181e-06 A`, identical to twelve digits — and `C_g` depends
on which patches are driven and grounded, so the generated conditions are
electrostatically indistinguishable from the ones they replace.


## Materials come from a cited library, not from a case

A dielectric region names what it is *made of*, in one line:

```
// constant/<region>/electricalProperties
material    borosilicateGlass;
```

and both `epsilonR` and `gammaSEE` follow from
[`etc/materials/dielectrics`](etc/materials/dielectrics), each with a
literature reference. Twelve materials ship today: `vacuum`, `air`, `ptfe`,
`pmma`, `polyimide`, `fusedSilica`, `borosilicateGlass`, `sodaLimeGlass`,
`mica`, `alumina96`, `alumina99`, `magnesia`.

**Precedence** — an explicit number always wins, and is reported as an
override; otherwise the material supplies it; otherwise the region kind's
default applies. Naming a material that is not in the library is **fatal** and
lists what is, rather than falling back to a default nobody chose.

`γ` is flagged in the library as **the least transferable number in it**: it
depends on the ion, its energy and above all the surface condition, and the
spread across surface states exceeds any difference between the materials. So
every dielectric carries `0.001`, the contaminated-barrier figure — the same
value and the same reasoning as the wall-flux `defaultSEEC`, so that `γ` cannot
depend on whether a case happened to name a material. `magnesia` is the one
genuine exception and says why. If `γ` matters to your result, it is a number
to measure or fit, not to take from a library.

Verified against analytic ground truth: with `material alumina96` the plate2D
series stack gives `V_interface = 0.1`, exactly `1/(1+εᵣ)` for `εᵣ = 9`, where
the case's own hand-typed 5.0 gives `1/6`.

Conductors are **deliberately absent**: there is no way yet for a case to say
what a metal *patch* is made of, so shipping copper/steel entries would be a
library nothing can reach. Until the patch-material route exists, an electrode's
`γ` is set with `defaultSEEC` on its wall-flux condition.


## Species boundary conditions are derived, not written per case

You declare what a boundary *is* — through the potential's boundary conditions
and `constant/regionProperties` — and `plasmaCreateSpeciesFields` derives what
every species does there:

| patch | charged species | `nEps_e` | surface charge |
|---|---|---|---|
| interface to a **dielectric** region | wall flux | energy wall flux | **accumulates** |
| a `thinDielectricPotential` surface | wall flux | energy wall flux | **accumulates** |
| a metal electrode (driven or grounded) | wall flux | energy wall flux | no |
| a `floatingElectrodePotential` conductor | wall flux | energy wall flux | no — becomes `Q(t)` |
| interface to a **farField** region | zeroGradient | zeroGradient | never |
| `empty` / `wedge` / `symmetry` / `processor` | kept as-is | kept as-is | — |

The condition follows the species **charge the mechanism declares**: the
electron gets one *with* secondary emission, every other charged species one
*without* (an ion arriving is not itself an emission event), and neutrals get
none. `nEps_e` is generated the same way but never accumulates surface charge —
an energy density is not a charge density.

**σ accumulates where the potential consumes it.** A dielectric surface charges
because something reads that charge; a metal electrode does not, because the
arriving charge is conducted away through the circuit, and a floating conductor
redistributes it into a *total* charge rather than a local σ.

**σ on a conductor is not a setting with two valid values — it is refused.**
Charge cannot sit still on a metal: it redistributes over `ε/σ ≈ 10⁻¹⁸ s`, so a
σ distribution on a conductor is not a quantity that exists. The generator
derives `false` there, and a case that overrides it to `true` **aborts**, naming
where the arriving charge actually went instead — `I_cond` in Sato's discharge
current for a driven or grounded electrode, `Q(t)` for a floating one.

Why this replaced hand-written per-species blocks: restating what the potential
already says is a second source of truth, and it drifted. `needleDBD` ended up
with a wall flux on `n_e` and `zeroGradient` on `nEps_e` — the pairing that
leaves the energy equation singular at a sharp electrode — because its
`0.orig/nEps_e` was a copy of the *streamer* case's template, naming `axis`,
`wedge_0` and `grounded_electrode`, patches that mesh does not have.

**The potential's own boundary conditions are generated too**, from
`configuration/boundaries` and [`etc/boundaryRoles`](etc/boundaryRoles) — see
the next section. An explicit block in `etc/changeDictionary.*` still wins over
either generator: it is applied after both, and is normally absent. `wallFluxFamily Mixed | Implicit` in
`system/plasmaSimulationControls` selects the family (`Mixed` by default).

## What the wall-flux conditions actually solve

Every `ddWallFlux*` condition implements chapter 6 of Hagelaar's HDR
(`Literature/`), which is the reference for the whole family. The net flux and
the flux coming back off the wall are

```
(6.1)  Gamma.n  = n*w_w - Gamma_w
(6.2)  Gamma_w  = r*n*w_w + SUM_j gamma_j n_j w_w,j
(6.8)  w_w      = max( vT/sqrt(pi) - Gamma_w/n , 0 ) + max( drift toward wall, 0 )
```

with `r` the reflection probability and the sum over incident species the
secondary emission. Note what (6.2) says: `Gamma_w` is built from the **total**
`w_w`, and (6.8) feeds `Gamma_w` back into `w_w`. **Reflection is therefore
implicit, not a prefactor** — it reflects the drift-driven incident flux too.
The two branches of the `max` partition exactly, so the relation has a closed
form:

```
w_w = max( (A + Dd - Gc/n)/(1 + r), Dd )        A = vT/sqrt(pi)
Gamma.n = (1 - r)*n*w_w - Gc
```

Set `Gc = 0` and `Dd = 0` and this is **eq. (6.7) exactly**, including the
`1/sqrt(pi)` — the one reflection case Hagelaar closes in the text, and the
ground truth `testWallFlux` checks. Combining `r` with (6.8)'s drift term is
not itself quoted from the source; it follows from reading (6.2) literally.

**Two consequences that look like bugs and are not.** The emission carries a
factor **2** at `r = 0` — that is just (6.1) closed with (6.6),
`Gamma = n(A - Gamma_w/n) - Gamma_w = nA - 2 Gamma_w`. And the `max` is
**physics, not a numerical guard**: Hagelaar notes it fires at a cathode under
ion-impact emission, where setting `w_w = 0` is correct *because published
gamma values were deduced neglecting thermal electron loss to the cathode*.

**The thermal base is eq. (6.6), which is TWICE eq. (6.3)** — the shifted
Maxwellian, not the centred one. `(6.3) = vT/(2 sqrt(pi))` is the classical
one-sided flux `(1/4) sqrt(8kT/pi m)`; (6.6) is `(1/2)` of that root, i.e.
`vT/sqrt(pi)`. Adopted 2026-09-04 in `7e6a4b3`.

### The electron energy weight is 5/3, not 4/3

Because the particle base is the shifted Maxwellian, the energy weight must be
its partner **eq. (6.15)**, not the centred **(6.14)**. The consistent pairs
are `(6.3, 6.14)` and `(6.6, 6.15)`; mixing them is what this repo did between
`7e6a4b3` and `2026-09-04`.

```
eps_w = ( 5/2 - (Gamma_w/n_e)/A ) T_e      =>   eps_w/eps = 5/3 - (2/3)(Gamma_w/n_e)/A
```

so the weight **falls** as emission and reflection rise, rather than being the
constant it used to be. Reflected electrons carry their energy back: splitting
(6.13)'s single `eps_s` by (6.2) into the reflected part (returning at the
incident `eps_w`) and the created part (born at the surface with `eps_s`) makes
the energy closure exactly parallel to the particle one.

Three checks, all reproduced exactly and all from the text: `5/3 = (5/4)(4/3)`,
the published `eps_w -> 2 T_e` limit as `r -> 1`, and the HDR's own energy
equation carrying `W = (5/3) n_e w_e`, which a BC weight of 4/3 contradicts.

### `electronReflection` defaults to 0, and the Implicit family has a limit

`r` is set per patch with `electronReflection` on either family. It defaults to
**0**, at which the closure reproduces the pre-2026-09-04 form algebraically
(measured: 3.7e-16 relative over a 372-point sweep — *algebraically*, not
bitwise, since the two forms associate the additions differently).

`ddWallFluxImplicit` carries reflection as three exact constant factors —
`(1-r)/(1+r)` on the thermal and drift terms, `2/(1+r)` on the creation source
— so it stays fully implicit in `n`. What it does **not** carry is (6.6)'s
`max(..., 0)` clamp, which needs `n` at the wall and would make an
intentionally-linear condition lagged and nonlinear. Consequence: under
emission strong enough that `Gamma_w/n` exceeds `A`, that family lets the net
wall flux become an electron *source*. **Use the Mixed family (the default) at
a cathode with significant SEE.**

`ddWallFluxImplicit` also has no electron-energy member, so
`wallFluxFamily Implicit` with LMEA is refused at setup rather than dying later
on an unregistered type.

Verified by [`testWallFlux`](src/applications/utilities/testWallFlux) — 26
checks, no mesh or case needed, including the closed form against an
independent fixed-point solve of (6.2)+(6.8) over 2430 `(r, Dd, Gc)` points
(1.7e-14) and two liveness controls that convict the wrong floor and the wrong
`(1+r)` placement.


## Two currents, and they are not the same thing

| | where | what it is |
|---|---|---|
| **`I_total`** | `postProcessing/dischargeCurrent/current.csv` | the **external-circuit** current from Sato's equation — weighted over the whole domain, and `I_cond + I_disp`, so it *includes* displacement. In a DBD it is dominated by `I_disp`. |
| **`I_collected`** | `postProcessing/floatingElectrode/floating.csv` | the net charge flux onto **one** conductor — a **conduction** current only, local to that patch. It is what charges a floating electrode and what `Q(t)` integrates, and it is broken down per species. |

Neither is a check on the other, and both CSV headers say so.


## Discharge current, on by default

Sato's discharge current is the primary measurable of almost every discharge
simulation — the one number an experiment can be compared against — so since
2026-09-03 it is computed **by default** for every plasma run, with no
configuration. It costs one extra Poisson solve at start-up.

The electrode patches are **derived from the potential's own boundary
conditions** rather than restated: a *driven* electrode is a Dirichlet
`ePotential` condition that is time-varying or a non-zero constant; a *ground*
is a non-time-varying Dirichlet equal to zero. The time-varying test matters,
because at `t = 0` a ramp reads exactly zero and a value-only rule would call
the driven electrode a ground. If the driven electrode is ambiguous the run
aborts and lists the candidates, since which electrode the current is measured
at is a physical choice. `dischargeCurrent { enabled false; }` is the opt-out,
and naming the patches explicitly still wins.

The ground may sit in **another region** — behind the barrier, as in any DBD.
The weighting field is solved monolithically across every region, validated to
`1.1e-15` against the analytic series-stack capacitance. Full reference:
[`docs/reference/plasmaSimulationControls.md`](docs/reference/plasmaSimulationControls.md).


## Running the tests

```bash
tools/run_electrostatics_tests.sh     # seconds; exits nonzero if any case fails
```

Runs every electrostatics-only case and checks each against the reference
declared in its own `COMPARE.md`. It exists because all of these cases were dead
for three weeks and nothing noticed: a commit made `backgroundDensity` fatal but
left `reducedE = Emag/backgroundDensity` unguarded in `singleRegionPoisson`, and
since no electrostatics-only solver has a `plasmaSpecies` to publish a density,
every case died with SIGFPE. It was found by accident, while testing an
unrelated boundary condition.

The plasma tutorials are hours long and cannot be run casually. These are
cheap, so there is no excuse — **a case that is not run is not tested, and it
rots silently.**


## Contributors

**Rention Pasolari** <r.pasolari@gmail.com>  
Software architecture, implementation, testing, validation, documentation, maintenance

**Konstantinos Kourtzanidis** <kkourtza@gmail.com>  
Conceptual design, physical modeling, debugging, testing, validation, scientific supervision
