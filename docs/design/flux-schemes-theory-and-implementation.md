# Flux schemes for drift-diffusion-reaction transport in SoPLASMA

**Theory, cell-centred finite-volume implementation, and verification.**

Written 2026-09-08. Covers the Scharfetter-Gummel (SG) scheme, the Complete
Flux Scheme (CFS) as implemented here, their non-orthogonal treatment, and the
verification beds that established the results. Everything quantitative below
is measured, and the beds that produced it are in `verification/`.

---

## 1. The problem

Each charged species obeys a drift-diffusion-reaction (DDR) equation

$$\partial_t n + \nabla\!\cdot\Gamma = s, \qquad \Gamma = v\,n - D\,\nabla n$$

with drift velocity `v = Z mu E`, diffusivity `D`, and a net volumetric source
`s` from the chemistry. In a discharge `v`, `D` and `s` all vary by orders of
magnitude across a sheath, and the two terms of `Gamma` can each be far larger
than their difference. How the FACE FLUX is built is therefore the dominant
accuracy question, not how the cell values are stored.

The **grid Peclet number** governs which regime a face is in:

$$P_{\rm grid} = \frac{v h}{D}$$

with `h` the owner-to-neighbour distance. Liu et al (2014) give it a physical
reading: with the Einstein relation `D = (2/3)(eps/q) mu`,

$$P = \tfrac{3}{2}\,\frac{qEd}{\varepsilon}$$

so the Peclet number is **the energy a particle gains from the field over the
distance `d`, relative to its mean energy**. `P_grid` is that gain across one
cell. For the Grubert dc-glow case at 5 um cells, `P_grid ~ 0.014` (measured),
consistent with Liu's own table extrapolated to that resolution.

---

## 2. Scharfetter-Gummel

### 2.1 Derivation

SG solves the local two-point boundary value problem across a face

$$\frac{d}{dx}\Big(v n - D \frac{dn}{dx}\Big) = 0,\quad n(0)=n_P,\ n(h)=n_N$$

**with the source set to zero**. The solution is exponential and the resulting
flux is exact for that problem:

$$\Gamma^h = \alpha n_P - \beta n_N,\qquad
\alpha = \frac{D}{h}B(-P),\quad \beta = \frac{D}{h}B(P)$$

with the Bernoulli function `B(z) = z/(e^z - 1)`.

### 2.2 What SG is and is not exact for

MEASURED (`verification/fluxScheme1D`), steady 1D DDR at `Pe = 1e4`:

| source | SG error, N = 20..160 |
|---|---|
| `s = 0` | 5.4e-16, 2.5e-15, 9.7e-15, ... (machine precision) |
| `s = 1` | 2.500e-6, 1.250e-6, 6.250e-7, 3.125e-7 |

**With no source SG is exact. The whole of its error is the source term it
omits from the flux.** That error is a UNIFORM offset (`L2 == Linf` to all
digits) equal to the source flux dropped over half a cell:

$$\text{low }P_{\rm grid}: \frac{S h^2}{8D} \ (\mathcal{O}(h^2))
\qquad
\text{high }P_{\rm grid}: \frac{S h}{2v}\ (\mathcal{O}(h))$$

both confirmed to six digits. **The linear-in-`h` scaling at high `P_grid` is
exactly why SG degrades to first order there** -- measured 1.00 over five
refinements at `P_grid` from 500 down to 15.6. SG's first-order behaviour and
its source-blindness are one defect seen from two angles, not two properties.

---

## 3. The Complete Flux Scheme

### 3.1 The idea

CFS solves the same local problem **with `s` retained**:

$$\frac{d}{dx}\Big(v n - D\frac{dn}{dx}\Big) = s(x)$$

Since `dGamma/dx = s`, the flux is no longer constant across the interval:
`Gamma(x) = Gamma_0 + S(x)` with `S(x) = int_0^x s`. Solving for `Gamma_0` and
evaluating at the face gives the exact result

$$\boxed{\ \Gamma^i = S(x_f) - \frac{\lambda}{1-e^{-P}}\int_0^h S(x)\,e^{-\lambda x}dx\ }
\qquad \lambda = v/D$$

and `Gamma = Gamma^h + Gamma^i`. **`Gamma^i` depends on `s` only, never on
`n`,** so it is an explicit contribution: CFS = the SG matrix plus
`div(Gamma^i)` on the right-hand side.

### 3.2 Constant source -- Liu's published coefficients

With `s` constant and the face at the interval MIDPOINT, the boxed expression
reduces to

$$\Gamma^i = s\,h\left(\tfrac12 - W(P)\right),\qquad
W(z) = \frac{e^z - 1 - z}{z(e^z-1)} = \frac{1-B(z)}{z}$$

which is Liu et al (2014) eq. (18c) with `gamma + delta = 1/2 - W`. Their
split `gamma = max(1/2-W, 0)`, `delta = min(1/2-W, 0)` takes `s` from the
upwind side when `1/2-W > 0` and the downwind side otherwise.

### 3.3 EXTENSION 1 -- the second-moment term

**Eq. (18c) is exact only for a constant source.** With `s` varying linearly
across the interval, `s = s_L + (s_R-s_L)x/h`, the boxed expression gives

$$\Gamma^i = h\Big[\,s_L\big(\tfrac12 - W(P)\big)
   + (s_R-s_L)\big(\tfrac18 - \tfrac12 V(P)\big)\Big]$$

$$V(z) = \frac{2 - e^{-z}(2 + 2z + z^2)}{z^2\,(1-e^{-z})}$$

At large `P`, `V -> 0` and (18c) drops `h(s_R-s_L)/8 ~ h^2 s'`. Checked
independently at `P -> 0` against pure diffusion (`-D n'' = s`), where the
exact inhomogeneous flux is `-h(s_R-s_L)/24` and
`1/8 - V(0)/2 = 1/8 - 1/6 = -1/24`. Agrees.

Series used near `z = 0`, where the direct forms are `0/0`:
`W = 1/2 - z/12 + z^3/720`, `V = 1/3 - z/12 + z^2/360`.

### 3.4 EXTENSION 2 -- boundary faces use different coefficients

This was the decisive correction. For an INTERIOR face the local problem spans
the two cell centres and **the face lies at its midpoint**, so `S(h/2) = s h/2`
supplies the `1/2`. For a BOUNDARY face the interval runs from the cell centre
to the face and **the flux is evaluated at its END**, where `S(H) = s H`
supplies `1`:

$$\text{interior}: \ h\Big[s_L(\tfrac12-W) + (s_R-s_L)(\tfrac18-\tfrac{V}{2})\Big]$$
$$\text{boundary}: \ H\Big[s_C(1-W) + (s_F-s_C)(\tfrac12-\tfrac{V}{2})\Big],
\quad H = 1/\text{deltaCoeffs}$$

Using the interior coefficients at the boundary **halves** that term. In a
drift-dominated problem the inflow boundary sets the error, so CFS then removed
exactly half of SG's error -- measured ratio 0.536, 0.518, 0.509, 0.505,
converging on 1/2. That clean factor of two is what identified it.

### 3.5 Result

With both extensions, on the uniform-source case at `Pe = 1e4` where `Gamma^i`
is exact:

| N | SG | CFS |
|---|---|---|
| 20 | 2.500e-06 | **1.11e-16** |
| 40 | 1.250e-06 | **7.85e-17** |
| 80 | 6.250e-07 | **6.57e-17** |
| 160 | 3.125e-07 | **1.31e-16** |

and on a manufactured solution `n = 1 + x + 0.5 sin(pi x)` at `Pe = 1e4`,
where SG is first order:

| N | `P_grid` | SG | order | CFS | order |
|---|---|---|---|---|---|
| 40 | 250 | 1.860e-02 | 0.99 | 9.123e-05 | **2.02** |
| 80 | 125 | 9.339e-03 | 0.99 | 2.275e-05 | **2.00** |
| 160 | 62.5 | 4.690e-03 | 0.99 | 5.683e-06 | **2.00** |
| 320 | 31.2 | 2.364e-03 | 0.99 | 1.421e-06 | **2.00** |
| 640 | 15.6 | 1.204e-03 | 0.97 | 3.551e-07 | **2.00** |

**CFS restores second order exactly where SG is first order** -- Liu's
"uniformly second order in Peclet" claim, reproduced -- and is 3390x more
accurate at N = 640 (147x at `Pe = 100`).

---

## 4. Non-orthogonal meshes

### 4.1 The defect

As shipped, `fvm::ScharfetterGummel` built its face conductance from
`mesh.deltaCoeffs()` -- the ORTHOGONAL delta coefficients -- with no
non-orthogonal correction of any kind. CFS inherited this, its homogeneous part
being the SG flux. MEASURED on an 11.3 deg sheared mesh
(`verification/fluxScheme2Dnonortho`):

| scheme | NX=10 | 20 | 40 | 80 | order | correctors used |
|---|---|---|---|---|---|---|
| standard (linear + corrected) | 5.141e-3 | 2.620e-3 | 1.335e-3 | 6.752e-4 | 0.97 | 15-16 |
| SG (before) | 1.009e-2 | 1.001e-2 | 1.004e-2 | 1.004e-2 | **0.00** | 2 |
| CFS (before) | 1.045e-2 | 1.012e-2 | 1.006e-2 | 1.005e-2 | **0.00** | 2 |

A fixed error floor untouched by refinement. They "converged" in 2 correctors
against standard's 15-16 because **their matrix had no correction term for the
correctors to act on** -- more correctors could not have helped.

### 4.2 The fix, and the identity that makes it well defined

Using the Bernoulli identity `B(-z) = B(z) + z`, the SG face flux factorises
EXACTLY:

$$\Gamma = \alpha n_P - \beta n_N
        = \underbrace{\phi\,n_P}_{\text{drift}}
        - \underbrace{D_f B(P)\,|S_f|\,\mathrm{snGrad}(n)}_{\text{diffusion}}$$

So **SG's diffusive part is an ordinary diffusion flux with effective
diffusivity `D B(P)`**, and takes the standard non-orthogonal treatment:

* implicit part uses `nonOrthDeltaCoeffs` instead of `deltaCoeffs`;
* explicit correction `-D_f B(P)|S_f|\,\mathrm{snGradCorrection}(n)` added to
  the flux, entering the matrix source as `source -= V\,div(\cdot)`.

AFTER, on the same mesh:

| scheme | NX=10 | 80 | order |
|---|---|---|---|
| SG | 5.282e-3 | 6.845e-4 | **0.97-1.00** |
| CFS | 5.745e-3 | 6.994e-4 | **1.00-1.02** |

and the orthogonal control is **bit-unchanged** at order 2.00 -- the correction
vectors are zero there, so the term vanishes identically. That is why it is
applied unconditionally with no user switch.

**Note on the remaining first order.** All three schemes -- including
`standard` -- sit at ~1.0 on this mesh. That is a property of the mesh and of
OpenFOAM's non-orthogonal treatment, not of SG or CFS: after the fix they are
at parity with the standard scheme, where before they did not converge at all.

### 4.3 Automatic correctors

A corrector count is a property of the MESH, not of the physics, so the user
should not supply it. `electromagneticsModel` now measures the largest face
non-orthogonality angle and derives the count (`<5 deg -> 0`, `<35 -> 2`,
`<60 -> 3`, else 4). The derived value is a **floor**, not merely a default: a
default would never fire because every generated case writes the key, so a
stale `0` would keep silently dropping the correction. Too few correctors
returns a wrong potential with no symptom; too many only costs time. A case may
still ask for more, and is warned when its value is raised.

The reduction over faces is unconditional -- a rank holding no internal faces
must still participate or the collective deadlocks (project rule 31).

**THE TWO CORRECTOR FAMILIES ARE SEPARATE AND MUST NOT BE CONFLATED.**

| | control | treatment |
|---|---|---|
| Poisson | `nNonOrthogonalCorrectors` in the `poisson` block | mesh-derived floor, automatic (above) |
| Transport (DDR) | `maxCorrectors`, the outer PIMPLE loop | no separate knob -- see below |

The drift-diffusion-reaction equations need NO non-orthogonal corrector setting
of their own. `soPlasmaFoam.C:420` calls `transport.solve()` INSIDE
`while (pimple.loop())`, and `plasmaTransport::solve()` calls `nEqn()` for every
species on each pass, so the transport matrices -- and with them the EXPLICIT
non-orthogonal correction now carried in the SG/CFS operator -- are re-assembled
every outer corrector. The correction therefore converges through the loop that
already exists, automatically, with nothing for the user to set. Adding a second
corrector count for transport would duplicate that loop.

This also explains why the standalone `testFluxScheme` utility needs its own
re-assembly loop: it is single-shot and has no outer iteration to borrow.

---

## 5. Implementation in SoPLASMA

**Selection.** `fluxScheme (standard | ScharfetterGummel | CompleteFlux)`,
set once via `driftDiffusionFluxScheme` in the case config. CFS is a THIRD
option, not a mode of SG: it differs from SG only by `Gamma^i`, and folding it
in would hide the one term that distinguishes them.

**Propagation.** The electron's setting propagates automatically to
(a) all derived ions -- `ionFluxScheme` now DEFAULTS to the electron's value
rather than to `standard`, so "use SG" no longer means "use SG for electrons
only"; and (b) the LMEA electron-energy equation, which is itself a
drift-diffusion equation for `n_eps` and inherits `fluxScheme` from the
electron transport model. An explicit override is available at each level.

**The source.** `s` is the NET chemistry source -- production minus loss across
the whole mechanism, `chemP - chemL*n`, not ionisation alone. It is handed to
the transport model through a no-op virtual `setNetSource()` on
`plasmaTransportModel`, so `plasmaTransport` needs no `driftDiffusion` cast.

**It is LAGGED one outer iteration**, and must be: `plasmaTransport::solve()`
builds the transport matrices before it assembles the chemistry sources, so the
current-iterate value does not exist at assembly time. `chemP_`/`chemL_` persist
as members, so no extra storage is needed, and the lag vanishes as the PIMPLE
loop converges. Liu et al treat the source explicitly for the same reason. The
energy equation needs no lag -- its `Psrc`/`Lsrc` are already available in the
same call.

**Wall boundary conditions.** `ddWallFluxMixed` builds its `valueFraction` from
the interior flux formulation, so `CompleteFlux` takes the SG branch -- CFS's
homogeneous part IS the SG flux. Falling through to the `standard` branch would
impose a face value extracted under the wrong formulation; `testWallFlux`
measures that mismatch at 18.4% of the imposed flux one way and 99.3% the
other. The inhomogeneous part is deliberately not represented in
`valueFraction`: it is an explicit face flux the CFS operator already adds on
boundary faces, and doing both would double-count it.

---

## 6. Verification method

Two beds, both with exact solutions:

* `verification/fluxScheme1D` -- steady 1D DDR **with a source**, 504 runs over
  {upwind, linear, limitedLinear, Minmod, vanLeer, MUSCL, SuperBee, ROUNDF,
  ROUNDA, ROUNDAplus, ROUNDL, ROUNDW, SG, CFS} x 6 Peclet x 6 meshes.
* `verification/fluxScheme2Dnonortho` -- 2D manufactured solution with the exact
  value imposed on EVERY patch, so the bed is valid on any mesh shape.

### Three methodological traps, all of which were fallen into

1. **A uniform source cannot test CFS.** `Gamma^i` is then identical on every
   interior face, `div(Gamma^i) == 0`, and CFS collapses onto SG in the interior
   (measured: 0.2% apart). The uniform source chosen to ISOLATE the source term
   is the least favourable case for the scheme that fixes it. Hence `-mms`.

2. **Asymptotic order always reports the low-`P_grid` limit.** `h -> 0` drives
   `P_grid -> 0` by construction, so a convergence table at domain `Pe <= 100`
   never samples the regime where SG degrades. Testing CFS there showed nothing.
   The domain Peclet must be raised until `P_grid >> 1` across the whole
   refinement -- hence the `Pe = 1e3, 1e4` sweeps.

3. **A 1D exact solution is invalid on a sheared mesh.** Shearing turns the
   left/right patches into slanted planes, so "n = a at the left boundary" stops
   being a condition at constant `x`. The boundary data then disagrees with the
   geometry by `O(shear)`, and the measured error floor tracked the shear ~1:1
   (0.01 -> 8.0e-3, 0.10 -> 1.07e-1, 0.20 -> 2.33e-1) for both schemes at every
   Peclet. The `shear = 0` row reproducing the 1D answer to all digits is what
   proved the harness sound. Fixed by the 2D manufactured solution.

Also: **non-orthogonal correctors must RE-ASSEMBLE the equation**, not re-solve
it. The correction enters the matrix source at assembly time; calling `solve()`
again on the same matrix changes nothing (measured: identical digits).

---

## 7. Known limitations

* **CFS and SG are verified on the schemes' own terms only for SMOOTH
  solutions.** Neither bed has extrema, so boundedness and oscillation control
  -- the whole purpose of a limiter -- are never exercised. SG and CFS are
  UNBOUNDED and can undershoot a density negative. Pasolari & Kourtzanidis
  (arXiv 2607.05137) find SG "excessively diffusive" on a positive streamer and
  recommend ROUNDF; that is not a contradiction, it measures a different
  property. **Neither scheme should become a default on this evidence.**
* CFS on non-orthogonal meshes is now at parity with `standard`, but neither is
  second order there.
* Boundary conditions other than `fixedValue` and the plasma wall conditions
  have not been exercised for CFS.
* The automatic corrector floor is verified on an orthogonal mesh only; the
  raise-and-warn path has not been exercised on a skewed production case.
* 2D/3D application use is untested beyond the verification beds.

---

## 8. References

1. D. L. Scharfetter and H. K. Gummel, IEEE Trans. Electron Devices **16**, 64 (1969).
2. J. H. M. ten Thije Boonkkamp and M. J. H. Anthonissen, J. Sci. Comput. **46**, 47 (2011).
3. L. Liu, J. van Dijk, J. H. M. ten Thije Boonkkamp, D. B. Mihailova, J. J. A. M. van der Mullen, Plasma Sources Sci. Technol. **23**, 015023 (2014).
4. L. Liu, J. van Dijk, J. H. M. ten Thije Boonkkamp, D. B. Mihailova, J. J. A. M. van der Mullen, J. Comput. Appl. Math. **250**, 229 (2013).
5. G. J. M. Hagelaar, *Modelling methods for low-temperature plasmas*, HDR, ch. 6.
6. C. DeChant et al, Plasma Sources Sci. Technol. **32**, 044006 (2023).
7. J. Teunissen, Plasma Sources Sci. Technol. **29**, 015010 (2020).
8. X. Deng, J. Comput. Phys. **481**, 112052 (2023).
9. S. Pasolari and K. Kourtzanidis, *SoPlasmaFoam*, arXiv:2607.05137.
