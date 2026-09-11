# The floating electrode

**Status: IMPLEMENTED and VALIDATED, electrostatics and plasma (2026-09-03).**
The design below is what the implementation satisfies; all three analytic tests
pass, and the charge ledger and Gauss closure are verified on a plasma case.

## What it is

A conductor inside the domain connected to **nothing** — no power supply, no
ground. Physically common: a probe, an isolated pin, a metal fragment, a
floating guard ring, or a segmented electrode with one segment disconnected.

It differs from both existing electrode kinds in a way that matters
computationally:

| kind | potential | charge |
|---|---|---|
| `drivenElectrode` | **known** (the waveform) | whatever the circuit supplies |
| `groundedElectrode` | **known** (0) | whatever flows to ground |
| `floatingElectrode` | **unknown — solved for** | **known** (conserved) |

So the usual boundary condition is inverted: we do not know `V`, we know `Q`.

## The two conditions

**1. It is an equipotential.** Charge in a conductor redistributes at the speed
of its own relaxation time, `ε/σ ~ 10⁻¹⁸ s` for a metal, which is
instantaneous against every plasma timescale. So

```
V = V_f   (uniform over the whole surface, value unknown)
```

Note this is *not* a Dirichlet condition with a known value, and it is not
`zeroGradient`. It is a **constraint**: one unknown scalar shared by every face
of the patch.

**2. Its total charge is conserved**, changing only by what the plasma delivers:

```
Q(t) = Q₀ + ∫₀ᵗ I_plasma(t') dt'

I_plasma = ∮ Σ_s q_s Γ_s·n dA          net charge flux onto the surface
```

`Q₀` is the initial charge, normally 0 for an initially neutral electrode. The
flux `Γ_s` is already computed for surface charging, so `I_plasma` costs
nothing new — but note that on a conductor there is **no local σ**: charge
redistributes, so only the *total* is meaningful. That is the essential
difference from `thinDielectricSurface`, where σ is local and stays where it
lands.

Gauss's law closes the system:

```
∮ D·n dA = Q(t)          the surface integral of the displacement equals Q
```

## Why this is exact in one solve, not an iteration

The obvious implementation is a shooting method: guess `V_f`, solve Poisson,
measure `∮D·n`, correct, repeat. That works but costs several solves per step.

It is unnecessary, **because Poisson is linear in V.** Superposition gives, for
a fixed space charge `ρ` and fixed other-electrode potentials,

```
V = V_ρ  +  V_f · V̂
```

where

```
V_ρ   solves ∇·(ε∇V) = -ρ    with the floating patch held at 0
V̂     solves ∇·(ε∇V) = 0     with the floating patch held at 1, every other
                              electrode at 0, and no space charge
```

`V̂` depends **only on geometry and permittivities**, so it is computed **once**
at start-up and reused for the whole run. Then

```
∮D·n dA = Q_ρ + V_f · C_self

  Q_ρ    = ∮ D_ρ·n dA        induced by the space charge and the other
                             electrodes, with the floating patch at 0
  C_self = ∮ D̂·n dA          the electrode's self-capacitance, a CONSTANT
```

Setting that equal to `Q(t)` gives the floating potential in closed form:

```
        Q(t) − Q_ρ
V_f  =  ───────────
          C_self
```

and the field follows by superposition with no second solve:

```
V = V_ρ + V_f · V̂
```

**Cost: one extra Poisson solve at start-up, and two surface integrals per
step.** No iteration, no relaxation parameter, no convergence criterion — the
answer is exact because the relation `Q(V_f)` is exactly linear.

### Why `C_self` is constant, and when it is not

`C_self = ∮ D̂·n dA` involves only `ε` and the mesh. It must be recomputed if
either changes: an **AMR** step, a moving mesh, or a `semiImplicit` Poisson
scheme — where the operator is `ε + Δt·σ_electrical` and therefore depends on
the *conductivity*, which changes every step. Under `semiImplicit` the
superposition still holds within a step, but `V̂` and `C_self` must be rebuilt
whenever the operator does. That is the main implementation subtlety and the
place a naive cache would be silently wrong.

## Secondary emission still applies

A floating electrode is metal. It emits, so its species boundary conditions are
the same wall fluxes as any other electrode, with γ from its material. What it
does *not* have is a local surface charge — the ions arriving contribute to
`I_plasma` and hence to `Q(t)`, not to a `surfCharge` field.

## Physical behaviour worth expecting

A floating electrode in a plasma charges **negative**. Electrons are far more
mobile than ions, so it collects electron current until its potential drops
enough to repel most of them — the floating potential of probe theory. It should
settle near

```
V_f − V_plasma  ≈  −(kT_e/2e) · ln(2π m_e/m_i)      (planar, Maxwellian)
```

which is a few times `−kT_e/e`. That is a **useful sanity check**: a floating
electrode that charges positive, or that keeps drifting without settling, means
either the sign of `I_plasma` is wrong or the electron flux to the wall is not
being resolved.

## Interface

```
floatingElectrode  { material copper;        // -> gamma
                     initialCharge 0; }      // Q₀ [C], default 0
```

`initialCharge` is the only genuinely free parameter and 0 is the honest default
— an electrode that has not been connected to anything starts neutral.

There is deliberately **no** potential entry: providing one would make it a
`drivenElectrode`, and accepting one here would let a case silently ask for a
contradiction.

## Reporting

The floating potential is a **result**, not an input, and it is the single most
informative number about such an electrode. It should be written per timestep
alongside the discharge current, with `Q(t)` and `I_plasma(t)`, so the settling
behaviour above can be checked rather than assumed.

## Implementation status

| piece | state |
|---|---|
| the equipotential constraint | **done** — `floatingElectrodePotential` (a BC that carries the value) plus `Foam::floatingElectrode` (which determines it) |
| `psi` and `C_self` | **done** — via the shared `unitPotentialField`, the same monolithic unit-potential solve Sato's current uses |
| the closed form `V_f = (Q − Q_ρ)/C_self` | **done**, as an exact post-solve correction |
| rebuild of `psi`/`C_self` when the operator moves | **done** — every step under `semiImplicit`, cached under `explicit` on a static mesh |
| `Q(t) = Q₀ + ∫I_plasma dt'` | **done** — fed from `plasmaTransport::updateSurfaceCharge`, which owns the species wall fluxes; verified to 7.5e-12 on a plasma case |
| the step-discard invariant for `Q` | **done** — baselined (`Q = Q_base + I·dt`, never `Q += I·dt`) *and* restored by `discardStep()` |
| more than one floating conductor | **not supported** — fatal. Two of them couple, so the constraint becomes an N×N mutual-capacitance matrix |
| a floating conductor inside a dielectric region | **not supported** — fatal |
| per-step reporting of `V_f`, `Q`, `I_plasma` | **done** — `postProcessing/floatingElectrode/floating.csv`, one row per step |

### How it is applied: a correction, not two solves

The derivation above suggests solving twice — once with the patch at zero for
`V_ρ`, then superposing. The implementation instead runs the *ordinary* solve
with the patch at its current value `V_f_old` and removes the residual charge
error afterwards:

```
dV_f = (Q_target − Q_measured)/C_self
V_f += dV_f
V   += dV_f · psi          in EVERY region
```

This is algebraically identical — `L(psi) = 0`, so adding `dV_f·psi` keeps the
field a solution while shifting the patch value by exactly `dV_f` — and it
leaves the main Poisson solve untouched. It is **exact, not iterative**: one
correction lands on `Q_target` by construction, because `Q(V_f)` is exactly
linear. Applying it to every region matters: `psi` is continuous across the
interfaces by construction of the monolithic solve, so correcting only the gas
would leave a jump at the dielectric face.

### The sign convention, which is where this breaks

`Q = +∮ε·snGrad(V) dA`, with **no** leading minus. The patch normal points out
of the fluid and *into* the metal, while a Gauss surface enclosing the electrode
has its outward normal pointing into the fluid, so `n_enclosing = −n_patch`.
`C_self` uses the identical convention, and it must: `dV_f` divides one by the
other, so a sign slip in either drives the charge *away* from its target and the
run diverges rather than being slightly wrong. The first version of `C_self` did
have it backwards and returned exactly `−C`, caught immediately by the
two-electrode test in `testDischargeCurrent`.

## The three tests it must pass

Each is analytic, cheap, and would catch a sign error immediately.

**All three pass. Measured 2026-09-03**, unit bed
`tutorials/electrostatics/singleRegionElectrostaticFoam/floatingElectrode`
(`./Allrun-sweep`; also run by `tools/run_electrostatics_tests.sh`).

1-D column, gap `L = 1 m`, plate area `A = 0.1 m²`, `top` driven at `V₀ = 1 V`,
`bottom` the floating conductor, `sides` `zeroGradient` so the problem is
exactly 1-D. Analytic `V_f = V₀ + Q/(ε₀A)`:

| case | `Q₀` [C] | analytic `V_f` | measured | err |
|---|---|---|---|---|
| `q0` | 0 | 1.0 | 1.0 | 0 |
| `qp` | `+ε₀A` | 2.0 | 2.0 | 0 |
| `qm` | `−ε₀A` | 0.0 | −2.72e-08 | 2.7e-8 |
| `sym` | 0, antisymmetric 2-D | 0.0 | 4.52e-09 | 4.5e-9 |

`C_self` measured `8.85419e-13 F` against analytic `ε₀A/L =
8.854187817620389e-13 F`.

**`sym` is the real test of the induced-charge term `Q_ρ`.** `xLow = +1 V`,
`xHigh = −1 V`, `top` `zeroGradient`, and the floating conductor spans the whole
bottom edge. The configuration is antisymmetric under `x → 1−x` with `V → −V`,
so `V_f = 0` *by symmetry* — and its `C_self` is `5.4639e-12 F`, a completely
different number from the 1-D cases, which is what shows the result is not an
artefact of the 1-D formula.

Note also that the 1-D cases are **not degenerate**: `Q_ρ = −ε₀A ≠ 0` there, so
the induced-charge term is exercised even in the simplest case.

Two further invariants are checked in the same sweep:

- **The conductor is an equipotential.** The spread of `V` over its faces is
  `0.0e+00` exactly in every case — not a tolerance that was relaxed.
- **The plasma guard fires.** A case with non-zero space charge is refused, so
  the missing charge ledger cannot silently produce a wrong answer.


1. **`Q₀ = 0`, no plasma.** With no space charge and no flux, `Q_ρ` is whatever
   the other electrodes induce and `V_f = −Q_ρ/C_self`. For a symmetric
   two-electrode geometry with ±V applied, a floating electrode on the midplane
   must sit at **0** by symmetry.
2. **`C_self` against a known geometry.** For a concentric-sphere or
   parallel-plate arrangement the self-capacitance is analytic; the computed
   `∮D̂·n` must match.
3. **A charged, plasma-free electrode.** Set `Q₀ ≠ 0` with no plasma: the
   potential must be exactly `Q₀/C_self` above the induced value, and the field
   must match the analytic solution for an isolated charged conductor.

## The plasma case: the charge ledger and Gauss closure

`Q(t) = Q₀ + ∫I_plasma dt'` is fed from `plasmaTransport::updateSurfaceCharge`,
which owns the species wall fluxes — **the same fluxes** that build the local
surface charge on a dielectric, and the same `dt`. Computing the current
anywhere else would let the two accountings disagree about either.

A conductor has **no local σ**: charge redistributes, so only the total is
meaningful. Its wall fluxes therefore feed `Q(t)` instead of a `surfCharge`
field, and they are counted *regardless* of `enableSurfaceCharging` — which
governs the local σ a metal cannot have. Setting `enableSurfaceCharging true` on
a floating conductor is a contradiction and is refused rather than quietly
reinterpreted.

### Why `psi` must be rebuilt every step, and the check that proves it

Under `scheme semiImplicit` — **the default** — the operator is `ε + Δt·σ`, not
`ε`, and both factors change every step. Superposition `V += dV_f·psi` is exact
*only* if `psi` solves the same operator as the field it corrects, so `psi` and
`C_self` are rebuilt every step there.

This is verified rather than asserted. After each correction the conductor must
hold exactly the charge the ledger says it holds, and that residual is reported
as a **closure** in volts (`|Q_measured − Q_target|/C_self`). If `psi` came from
the wrong operator the closure is of order the potential itself; a mismatch
aborts the run.

### Measured 2026-09-03

Case `tutorials/plasma/soPlasmaFoam/needleDBD_floatingStrip` — `needleDBD` with
`air_dielectric` (178 faces) reinterpreted as an isolated metal strip, chosen
because that patch already carries real species wall-flux conditions.
`./Allrun-derive` copies the parent's mesh, tables and fields rather than
regenerating them, so the two cases cannot drift apart.

| quantity | measured | meaning |
|---|---|---|
| `C_self` | 7.03018250107e-16 F | the strip's self-capacitance |
| `psi` rebuild | every step | semiImplicit operator, as required |
| max closure | **2.82e-14 V** | Gauss's law closes at machine level |
| ledger consistency | **7.5e-12** rel. | worst `|Q_k − (Q_{k-1} + I_k Δt_k)|/|Q_k|` |
| `I_plasma` | −2.42e-12 A | **negative** — electrons arriving |
| `Q` | −2.37e-22 C | **negative** — charging negative, as expected |
| `V_f` | 7.24693 V | induction-dominated at this time |

**Honest limits of that case.** At `1e-10 s` the charge contributes
`Q/C_self = −3.4e-07 V` out of 7.25 V — about 2×10⁷ times smaller than the
induced part. So it validates the **ledger** and the **closure**, *not* the
settling to a floating potential; the probe-theory estimate above needs a far
longer run. And only `n_e` has a wall-flux condition on that patch in
`needleDBD` (`nEps_e` and every ion are `zeroGradient` there), so `I_plasma` is
**electron-only** — which makes the negative sign expected rather than
surprising, though it still checks the sign convention.

### The step-discard invariant

`Q` is accumulated solution state, so a step that is thrown away and re-run must
not leave its charge behind. Both mechanisms are used deliberately: `Q` is
**baselined** (`Q = Q_base + I·Δt`, never `Q += I·Δt`, so calling the update
twice in one step recomputes instead of accumulating — the same pattern
`surfCharge` uses) **and restored** by `discardStep()`, called from the solver's
rejection block beside `transport.discardStep()`, so correctness does not depend
on how a retry treats `timeIndex`.

### A diagnostic that could not be reconciled

The CSV first wrote one row per *outer iteration* — 129 rows for 18 steps — and
because `correct()` runs at the start of a step while the ledger updates at the
end, a row carried step *k*'s time beside step *k−1*'s charge. Reconciling `Q`
against `∫I dt` then missed by **11.5%**, which looks exactly like a quadrature
error and is not one. It now emits one row per step with all quantities from the
same step, which is what makes the 7.5e-12 above meaningful.

## Relation to the other kinds

`floatingElectrode` completes a set that is now closed under the two things a
boundary can know:

```
                     potential known    charge known
  drivenElectrode         yes               no
  groundedElectrode       yes               no
  floatingElectrode       no                yes
  thinDielectricSurface   no (Robin)        local sigma, not redistributed
```

Surface charging is **derived, not switched**: charge accumulates on exactly
those patches whose `ePotential` boundary condition consumes it. A dielectric
surface therefore charges because something reads the charge; a metal electrode
does not, because nothing does. `enableSurfaceCharging` survives only as an
explicit override.

`thinDielectricSurface` and `floatingElectrode` are the two "charge is the
input" cases, and they differ in exactly one respect: whether charge
**redistributes** over the surface (conductor) or **stays where it lands**
(insulator). That single physical distinction produces a global constraint in
one case and a local condition in the other.
