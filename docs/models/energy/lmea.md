# The electron energy model: LFA and LMEA

## One switch

All electron-energy configuration is a single top-level key in
`constant/plasmaSpeciesProperties`:

```
electronEnergyModel  LMEA;      // LFA | LMEA   -- REQUIRED, no default
```

Gas heating is a **separate, orthogonal** switch:

```
backgroundGas
{
    energy
    {
        solve   true;   // solve the gas temperature (incl. the vibrational reservoir)
        T       300;    // fixed value when solve = false; initial value when true
    }
}
```

All four combinations are supported and neither switch implies the other:

|                        | `solve false` | `solve true` |
| ---------------------- | ------------- | ------------ |
| `electronEnergyModel LFA`  | ✓ | ✓ |
| `electronEnergyModel LMEA` | ✓ | ✓ |

**The key is required and has no default.** LFA and LMEA are different physics
that give different answers, so neither is assumed for you; a case without it
stops at start-up with an error explaining both. Note this is the opposite of the
policy for `energyModelCoeffs`, and deliberately so: that block has exactly one
correct content given the closure, so requiring it is a migration tax, whereas
the closure itself is a genuine choice.

Heavy species (ions, neutrals) have **no** energy key. They follow
`backgroundGas/energy`. An earlier per-species `energyModel` key existed with a
second, incompatible vocabulary (`isothermal`/`backgroundGas`/`localField`/
`solveEnergy`) that filled five index lists with no readers; it was removed on
2026-09-01. The `localField` energy model went with it: it interpolated a
hand-supplied `T(E/N)` table to publish `T_<species>`, a field nothing reads
(`plasmaEnergy::T(label)` has no callers), and the sweep already writes
`meanEnergy_vs_reducedE`. Under LMEA `T_e = (2/3)⟨ε⟩e/k_B` follows from the
solved `meanE`; under LFA no electron temperature is wanted, and none is even
created — `plasmaEnergy` is not constructed unless gas heating or LMEA asks for
it. The selectable models are now `gasTemperature`, `isothermal`, `localEnergy`. The per-species spelling is **rejected**, on the electron and on
heavy species alike, with an error naming the replacement — two live spellings of
one setting is the problem this change exists to remove.

## What the two closures are

The two-term Boltzmann solution supplies every electron transport coefficient
and rate coefficient as a function of a single scalar. The two closures differ
in **which** scalar.

**LFA — local field approximation.** Coefficients are read at the local reduced
field `E/N`. This assumes the electron energy distribution is in equilibrium
with the local field, i.e. that an electron gains and loses energy over a
distance short compared with the scale on which `E/N` varies. It is the standard
closure for streamer modelling and is accurate where the field gradient is mild
relative to the energy relaxation length.

**LMEA — local mean energy approximation.** An electron energy-density equation
is solved,

```
d(n_eps)/dt + div( n_eps * mu_eps * E  -  D_eps * grad(n_eps) )  =  S_Joule - P_loss
```

with mean energy `eps_bar = n_eps / n_e`, and every electron coefficient and
rate is read at `eps_bar` instead of at `E/N`. This relaxes the equilibrium
assumption: the electron energy is transported, so it can lag the field. That
matters exactly where LFA breaks — in a streamer head, where `E/N` varies over a
length comparable with the relaxation length, and near boundaries.

The distinction and its limits are set out in
Hagelaar & Pitchford, *Plasma Sources Sci. Technol.* **14** (2005) 722, which is
also the reference for the two-term formulation the sweep implements.

## What the switch derives

Setting `LMEA` alone is a complete configuration. It:

1. solves the electron energy-density equation;
2. keys the electron's `driftDiffusionCoeffs` on `meanE`;
3. keys the chemistry rate tables on `meanE`;
4. fills in the whole `energyModelCoeffs` block from the mechanism's tables.

Every one of these is **derived, not read**. An explicit entry that *contradicts*
the closure is fatal, because the contradiction is the "half-LMEA" — transport
responding to the mean energy while ionisation still responds to the field, or
the reverse — which is a measured runaway, not a slightly different model. It is
fatal rather than warned because it produces no error of its own: plausible
fields, a clean run, and wrong physics.

Ion coefficients are never touched. Ion mobilities and diffusivities are
functions of `E/N` or `|E|`, always, and their blocks deliberately omit
`lookupVariable` so the `reducedE` default applies.

## The defaulted `energyModelCoeffs`

Written by the Boltzmann sweep into `constant/plasmaTables` (the directory is
taken from `chemistry/tableDir`). A case that writes any of these leaves keeps
its own version **whole** — see the merge rule below.

| leaf | type | table | keyed on |
| ---- | ---- | ----- | -------- |
| `mobility` | `fromMechanism` | `muN` | `meanE` |
| `diffusivity` | `fromMechanism` | `DLN` | `meanE` |
| `energyMobility` | `fromMechanism` | `muEpsN` | `meanE` |
| `energyDiffusivity` | `fromMechanism` | `DEpsN` | `meanE` |
| `powerLoss/elastic` | `tabulated1D` | `PelasticN_vs_meanE` | `meanE` |
| `powerLoss/inelastic` | `tabulated1D` | `PinelasticN_vs_meanE` | `meanE` |
| `initialMeanEnergy` | `tabulated1D` | `meanEnergy_vs_reducedE` | **`reducedE`** |

`reportInterval` is defaulted to **25 on this path only** (the model's own
default stays 1 for any case that writes its own block). The per-term dump is a
whole-mesh evaluation of `ddt`/`div`/`laplacian` once per *corrector*, so its
cost is mesh-dependent and severe exactly where it matters:

| mesh | reportInterval 1 → 25 |
| ---- | --------------------- |
| 1.15M-cell, 2 ns streamer | ~44 s → ~23 s **per step** (nearly 2×) |
| 130x130 fast bed (2026-09-01) | 76.94 s → 75.68 s total (1.7%) |

At small size the dump is negligible next to the step; at production size it
nearly doubles the run. 1 is therefore the wrong value to hand a new user, who
would reasonably conclude "LMEA is slow" from a diagnostics setting.

### Merge rule

> If you write a property sub-dictionary, it is yours entirely. If you omit it,
> it is derived entirely. There is no half-inherited sub-dictionary.

This is deliberately *not* OpenFOAM's `dictionary::merge`, which recurses into
sub-dictionaries. Under a recursive merge, a user writing

```
mobility { lookupVariable reducedE; }
```

would inherit `type fromMechanism; quantity muN;` from the default and get a
**working half-LMEA**. Under the replace rule the entry is incomplete and the
validator rejects the contradicting key outright.

`powerLoss` is the one exception, and it defaults **per leaf**: a case that
writes only `elastic` still gets `inelastic`. Omitting a loss channel can only
make the equation wrong in the runaway direction, and a case that genuinely
wants none can say so explicitly with `inelastic { type constant; value 0; }`.

## Two traps worth knowing

Both were measured, and both are why the defaults are shaped as they are.

### 1. `fromMechanism` divides by the gas number density

Correct for the *reduced* quantities `muN`, `DLN`, `muEpsN`, `DEpsN`, whose
tables hold `mu*N` and `D*N` — the evaluator divides by `N`, so the values track
density automatically.

Wrong for `PelasticN`/`PinelasticN`, which are **already per density**, and for
`meanEnergy`, which is an absolute energy in eV. Running those through
`fromMechanism` divides a second time: the loss term comes out **2.4e25 times
too small**, the electrons heat with no sink, and the mean energy climbs to the
`meanEnergyMax` clamp. The seed comes out at ~4e-26 eV.

Hence `tabulated1D` with an explicit `file` for those three, and a fatal error if
a case asks for `fromMechanism` there.

### 2. The energy coefficients are not `5/3` times the particle ones

If `energyMobility`/`energyDiffusivity` are absent the model falls back to
`(5/3) x` the particle values, which is the Maxwellian relation. Measured
departure on this air set:

```
E/N [Td]    mu_eps/mu    D_eps/D_T     (Maxwellian: 1.6667 both)
    13        1.557        1.062
   106        1.500        1.505
   428        1.492        1.305
  1723        1.545        1.467
```

i.e. wrong by up to 35% in mobility and 60% in diffusion. The defaults therefore
supply the real EEDF moments `muEpsN`/`DEpsN`, and the start-up report states
which of the two applied. Note these two are **all-or-nothing** in the model:
supplying only one falls back to Maxwellian for both.

## When LMEA is not available

The sweep writes `*_vs_meanE` tables only where the mean electron energy is
**strictly increasing** with `E/N`. The reason is invertibility, not cost: a
table keyed on `eps_bar` is used as a lookup, so the map `E/N -> eps_bar` must be
one-to-one or a single mean energy would correspond to two fields with two
different coefficients.

This fails in **attachment-dominated** mixtures. Attachment removes electrons
energy-selectively — dissociative attachment in O2 peaks near 6-8 eV, taking the
high-energy tail — so over some range of `E/N` the mean energy of the *surviving*
population can fall as the field rises.

The sweep handles this gracefully: it excludes the thermal plateau (below the
thermal floor every point returns the same `1.5 kT`, a flat and equally
non-invertible stretch), trims to the longest monotonic suffix so the high-field
range where a mean-energy key is actually wanted is kept, and reports what it
trimmed. It gives up only if the trimmed range is still non-monotonic.

Where it does give up, **the LMEA closure itself does not apply** — LMEA presumes
`eps_bar` determines the distribution shape, and where the map is not invertible
it does not. Asking for `LMEA` on such a mechanism is therefore an error, not
something to silently substitute: the solver reports the non-invertibility as the
cause and names the remedies, rather than quietly running a different closure
than the one requested or failing on a missing file.

Measured 2026-09-01, the repository's air mechanism is unaffected:
`muN_vs_meanE` and `muN_vs_reducedE` both have 263 rows (nothing trimmed),
strictly increasing from 0.137 eV.

## Start-up report

Provenance is always printed, never silent — a reader must be able to tell
whether a coefficient came from the Boltzmann solution or from an assumption
about it:

```
plasmaSpecies: electronEnergyModel LMEA -- electron coefficients and reaction
    rates are keyed on `meanE`.
plasmaEnergy: energyModelCoeffs for `e` resolved from electronEnergyModel LMEA.
    tables: "constant/plasmaTables"
    DERIVED from the model (the case omitted these): (mobility (muN_vs_meanE)
      diffusivity (DLN_vs_meanE) energyMobility (muEpsN_vs_meanE) ...)
      -- every one is the Boltzmann sweep's own EEDF moment; nothing here is a
      fit or an assumed ratio.
    READ from the case (used whole, the default dropped): (...)
```

## See also

- `docs/models/transport/drift_diffusion/driftDiffusion.md` — the species equations
  the electron closure feeds.
- `tutorials/plasma/soPlasmaFoam/positiveStreamer/positiveStreamer_LMEA_fast` —
  the reference LMEA bed, which keeps the full explicit block as a regression
  test of the override path.
