# The floating electrode

**Status: SPECIFIED, not implemented (2026-09-03).** This document is the
design; the derivation below is what the implementation must satisfy, and the
three limits at the end are its tests.

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

## The three tests it must pass

Each is analytic, cheap, and would catch a sign error immediately.

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
