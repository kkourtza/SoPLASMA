# positiveStreamer_LMEA_minimal — LMEA from one switch

The same physics as `positiveStreamer_LMEA_fast`, configured with **one key**.

## What this case demonstrates

The entire electron energy configuration is:

```
electronEnergyModel  LMEA;
```

at the top of `constant/plasmaSpeciesProperties`. Nothing else. The electron's
species block carries no energy entry, no `energyModelCoeffs`, and no
`lookupVariable` anywhere; `chemistry` in `plasmaTransportProperties` carries no
`lookupVariable` or `tableKey`.

All of it is derived or defaulted from that one key:

| what | derived/defaulted to |
| ---- | -------------------- |
| electron energy-density equation | solved |
| electron `driftDiffusionCoeffs/*/lookupVariable` | `meanE` (derived) |
| `chemistry/lookupVariable`, `tableKey` | `meanE` (derived / defaulted) |
| `energyModelCoeffs` (7 leaves) | from the mechanism's tables |

## Why it is a test and not just a demo

Verified 2026-09-01: this case produces **byte-identical fields at the final
time** (all 45 written fields) to `positiveStreamer_LMEA_fast`, which keeps the
explicit 80-line `energyModelCoeffs` block. The two beds together test both
directions:

- **this case** — the defaults reproduce a correct configuration from nothing;
- **`positiveStreamer_LMEA_fast`** — a hand-written block still overrides them.

If they ever diverge, the defaults have drifted from the block they are meant to
reproduce, and that is a defect rather than a tolerance question.

## Running

```
./Allrun-serial
```

~75 s on the 130x130 bed, 49 steps to 5e-10 s. Same cost as the explicit case
(75.64 s vs 75.68 s measured).

## What to look for in the log

The provenance report, which exists so that a reader can always tell a
table-derived coefficient from an assumed one:

```
plasmaSpecies: electronEnergyModel LMEA (read) -- electron coefficients and
    reaction rates are keyed on `meanE`.
plasmaSpecies: electron driftDiffusionCoeffs (mobility diffusivity) keyed on
    `meanE`, DERIVED from electronEnergyModel.
plasmaEnergy: energyModelCoeffs for `e` resolved from electronEnergyModel LMEA.
    DERIVED from the model (the case omitted these): (mobility (muN_vs_meanE) ...)
```

See `docs/models/energy/lmea.md` for the closures, the defaults table, and the
two measured traps the defaults are shaped to avoid.
