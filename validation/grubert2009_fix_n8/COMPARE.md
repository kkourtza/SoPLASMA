# Does the discharge screen and settle now that Poisson sees the charge?

```compare
question:  With the frozen-charge-density regression fixed, does the discharge screen its own field and settle near Grubert's operating point instead of running away?
baseline:  /home/kkourtza/soplasma-scratch/validation/grubert2009_fix_n11
varies:    minNumberDensity
matches:   resistance, capacitance, sourceVoltage, driftDiffusionFluxScheme, electronEnergyModel, mesh, chemistry, endTime
time:      4.5e-5
field:     n_e
region:    domain
```

## Why this pair exists

The three-arm R sweep of 2026-09-06 is VOID: it ran with Poisson's source
frozen (151 charge-density updates in 364670 steps), so the field was the
vacuum field, nothing screened, and the runaway it was built to study was the
bug. Fixed in `23cd599`.

So the first question is no longer "which R" but **"does it behave at all"**,
and the R sweep is only worth repeating if it does not.

## What differs, and what must match

**`varies: minNumberDensity`** -- the background floor, which IS the initial
condition (there is no seed Gaussian in this case).

- `grubert2009_fix_n11`: 1e11 m^-3 on both `e` and `Arp` (the shipped value)
- `grubert2009_fix_n8` : 1e8  m^-3 on both (the user's proposal, 2026-09-06)

Set on **both** species together: they must stay equal or the case starts
non-neutral, which would inject exactly the space charge whose response is
being measured. `Ar2p` keeps its 1e5 token floor in both.

Everything else is byte-identical, including R = 1e8, C = 5e-14,
V_src = -602 V, the mesh, and `fluxScheme standard`.

## Why a lower floor might matter

The floor is a background reservoir everywhere, so it is also a background
ionisation source. At 1e11 the pre-breakdown gap already carries 1e11 m^-3 of
plasma over the whole 1 cm; at 1e8 the Townsend phase starts three decades
lower and the avalanche has to build the plasma itself, which is closer to how
the reference case is posed.

**The risk, stated in advance:** LMEA recovers the mean energy as
`n_eps/n_e`, so a very low `n_e` floor can make that quotient stiff or noisy.
`nEfloor_` exists for exactly this. If the 1e8 arm becomes stiff where the 1e11
arm does not, that is a result about the floor, not a failure of the fix --
report it as such rather than lowering the floor further.

## Reference numbers (Grubert et al. 2009, argon, 100 Pa, 1 cm, -500 V)

| quantity | reference |
|---|---|
| j at the cathode | 0.511 mA/cm^2 |
| n_e peak | 2.478e15 m^-3 |
| n_Ar+ peak | 6.39e15 m^-3 |
| ionisation degree | 1.03e-7 |

Independent check: screening arrests growth when tau_diel = tau_ionisation,
giving n_e ~ 1.6e16 -- within 6.4x of the reference.

## THE OBSERVABLE THAT MATTERS MOST, and it is new

Not `n_e`. **Whether `Emag` develops a spread.** With the source frozen the
field was uniform at the vacuum value (12943..12989 V/m, 0.4% spread); with the
fix it is structured (8062..14310, 78%). A run whose `Emag` min and max stay
within a few percent of each other has NOT screened, whatever its densities
do, and should be treated as suspect immediately rather than at the end.

Extraction: `~/ct-env/bin/python status.py`, and
`grep -c "Charge density updated"` against the step count -- it must be >= 1
per step, never 0.04%.

## `region: domain` -- band regions are WRONG for this case

`compare_cases.py`'s `cathode-band`/`anode-band`/`bulk` slice along **y**, but
this gap runs along **x** (`gap1cm.geo`: cathode at x = 0, anode at x = L) with
y the 200 um transverse direction closed by symmetry planes. Recorded
2026-09-06.

---

# STRUCTURAL OBSERVABLE, added 2026-09-06 after I misread the field

At t = 2.00 us (V = -190 V, still on the RC ramp) the profile is a **TOWNSEND
AVALANCHE, not a glow**, and the summary statistic I had been quoting hid that:

| z/d | n_e | n_Ar+ | \|E\| [V/m] | E/N [Td] |
|---|---|---|---|---|
| 0.0002 (cathode) | 1.00e11 | 4.14e11 | 2.549e4 | 1056 |
| 0.16 | 3.14e11 | 1.39e13 | 2.529e4 | 1047 |
| 0.50 | 4.49e12 | 8.26e13 | 2.306e4 | 955 |
| **0.98** | **5.90e14** | **1.29e15** | **28** | ~0 |
| 1.0 (anode) | 5.64e14 | 5.23e14 | 2181 | 90 |

The field is UNIFORM from the cathode to z/d ~ 0.35 and collapses only in the
anode-side plasma. **The plasma formed at the ANODE and screens the field
THERE.** I read the 99.9% `Emag` spread as a cathode fall; it was screening at
the wrong end of the gap. A domain-wide min/max cannot tell those apart --
which is precisely why rule 26 requires PROBES.

This is a legitimate intermediate state: electrons multiply from the cathode
seed toward the anode, so the plasma appears there first, and the glow forms
only once ions drift BACK to build the cathode fall -- ~4.5 us of ion transit,
against t = 2 us here.

## THE OBSERVABLES THAT ACTUALLY DISCRIMINATE

Not the spread. In order:

1. **WHERE the density peaks are.** Townsend: both near the anode
   (z/d ~ 0.97, measured). Glow: the negative glow sits just downstream of a
   thin cathode fall, so the peaks move to the CATHODE side.
2. **WHERE the field maximum is, and how localised.** Glow: a fall of order
   4.4 mm at 100 Pa (`d_c*p ~ 0.33 torr cm`), with a near-field-free bulk
   beyond it. Measured here: |E| flat over the cathode-side 35% -- no fall.
3. **The sign structure of `chargeDensity`.** A cathode fall is a positive
   space-charge layer AT THE CATHODE. Measured here: rho > 0 everywhere,
   peaking at z/d = 0.96.

Only when (1) and (2) have flipped is a comparison against Grubert's figure 3
meaningful. Until then the case is mid-formation and quoting its profiles
against the reference would be comparing different physical states.

Extraction: `extract_centreline.py <timeDir>` after
`postProcess -func writeCellCentres -time <t>`, then
`SoEEDF/validation/dias2025/plot_grubert_profiles.py <csv> --arm LMEA`.
