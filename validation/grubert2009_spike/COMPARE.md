# grubert2009_spike -- spatial profiles through the current spike

Created 2026-09-07.

## The question

The ballasted Grubert case (`grubert2009_ballast_clean`) ignites, overshoots, and
then FREEZES: `deltaT` collapses from 6.6e-10 s to 6.0e-14 s and the clock stalls
at t = 1.0429072e-05 s.  This case restarts from that run's last usable snapshot
with frequent writing, to answer:

  WHERE are the high densities / fields / Te, HOW high are they against the
  reference and against physical ceilings, and IS A RECOVERY ABOUT TO HAPPEN?

## Baseline

  /home/kkourtza/soplasma-scratch/validation/grubert2009_ballast_clean

Restarted from its `9.99987e-06` snapshot (`startFrom startTime`,
`startTime 9.99987e-06`, `timePrecision 12`, `writeControl timeStep`,
`writeInterval 400`, `endTime 1.06e-5`).  Everything else is UNCHANGED from the
baseline, so the only difference is write frequency -- this is a diagnostic
re-run, not a physics variant.

Reference: Grubert, Becker & Loffhagen, Phys. Rev. E 80 (2009) 036405 --
steady normal glow, 100 Pa argon, 1 cm gap, gamma = 0.06:

  n_e   peak 2.478e15 m^-3 at x/L = 0.394
  n_Ar+ peak 6.391e15 m^-3 at x/L = 0.203
  cathode fall E/N ~ 2450 Td
  j = 0.511 mA/cm^2 ; V = 500 V

## Extraction

    ~/ct-env/bin/python dump_profiles.py            # -> profiles/<time>.csv, 400 pts each
    ~/ct-env/bin/python plot_profiles.py profiles   # -> profiles_logx.png, profiles_cathode.png
    ~/ct-env/bin/python profile.py <timedir>        # single-snapshot text table

`Cx/Cy/Cz` are copied from the `1.03256868404e-05` snapshot into the later time
dirs.  This is valid ONLY because the mesh is static: 2000 cells in every
snapshot, no `polyMesh` in any time dir, no `dynamicMeshDict`.  Verified
2026-09-07.

NOTE: `reducedE` is written to the snapshots as `internalField uniform 0` -- a
PLACEHOLDER.  The solver uses it as a table lookup key internally and never
writes real values.  E/N here is reduced from the `E` VECTOR field, which IS
populated, using N = p/(k_B T) = 2.4144e22 m^-3 at 100 Pa / 300 K.  Anything
that reads `reducedE` from a snapshot reads zero.

## Measured, at the frozen state t = 1.042906962e-05 s

| quantity | measured | reference | ratio |
|---|---|---|---|
| n_e peak | 2.482e18 m^-3 at x/L = 0.022 | 2.478e15 @ 0.394 | 1002x |
| n_Ar+ peak | 2.529e18 m^-3 at x/L = 0.020 | 6.391e15 @ 0.203 | 396x |
| E/N at cathode | 24835 Td | 2450 Td | 10.1x |
| Te (= 2/3 <U>) at cathode | 66.66 eV | -- | AT the 100 eV clamp |
| j (= I_cond/A) | 50.9 mA/cm^2 | 0.511 | 99.6x |
| V_electrode | +208.3 V | -500 V (cathode) | SIGN REVERSED |

## Why this is not a numerical blow-up (measured 2026-09-07)

- Sheath is RESOLVED.  Mesh is graded: dx = 1.37 um at the cathode.  At the n_e
  peak lambda_D = 34.6 um = 25 cells.
- Growth rate is BELOW the physical ceiling.  Between the last two snapshots
  n_e peak went 1.850e18 -> 2.482e18 over dt = 1.167e-10 s, i.e. 2.52e9 1/s.
  The absolute ceiling N*k_ion,max over the whole mechanism table is 6.59e9 1/s
  (k_ion,max = 2.728e-13 m^3/s at 2644 eV) -> we are at 0.38x.  At the LOCAL
  mean energy (100 eV, k_ion = 1.236e-13) nu_i = 2.98e9 1/s -> we are at 0.84x.
  The avalanche is exactly what the local field and energy prescribe.
- The mean-energy table is NOT extrapolated: `ENmax 55000` gives a table axis
  running to 2644 eV, so 100 eV is well inside it.  (Contrast the needleDBD air
  mechanism, where a 2000 Td sweep topped out at 29.85 eV -- that concern does
  not apply here.)

## Why the recovery cannot happen -- the structural result

`C_g = 1.7708e-16 F` (A = 2.0e-7 m^2, d = 1 cm) and `R = 1e8 Ohm`, so
tau_RC = R*C_g = 17.71 ns.  Note tau_RC = (R*A)*eps0/d is AREA-INDEPENDENT, so
this is not an artefact of the 1D cross-section.

Avalanche e-folding times from `k_EI_AR_ION_AR_vs_reducedE`:

| E/N [Td] | 1/nu_i [ns] | e-folds in one tau_RC |
|---|---|---|
| 431 (breakdown, 104 V / 1 cm) | 14.19 | 1.2 |
| 2450 (normal cathode fall) | 0.73 | **24.2** |
| 25000 (frozen state) | 0.20 | **87.6** |

R is pinned by the operating point: I_op = 0.511 mA/cm^2 * 2e-3 cm^2 = 1.022 uA,
and a 300 V ballast drop needs R = 2.94e8 Ohm.  We used 1e8, i.e. 3x SMALLER
(faster) than the operating point demands.  So:

  **A two-terminal RC ballast that sets this operating point NECESSARILY has
  tau_RC >> tau_ion.  The overshoot is structural, not a bad parameter choice.**

**SUPERSEDED 2026-09-07 AS TO CAUSAL ORDER, by this case's own snapshots (see
"THE ORDERING RESULT" below). The electrode flip is a CONSEQUENCE, not the
trigger: an interior potential MAXIMUM had already formed mid-gap while the
electrode was still at -90 V, the correct polarity.** The mechanics described in
this paragraph are accurate; the implied ordering is not.

Once I_cond reaches 1.02e-4 A, dV/dt = I/C_g ~ 5.7e11 V/s slews the electrode
straight through 0 V into positive: the phi panel shows the family marching
+25 -> +56 -> +100 -> +158 -> +208 V monotonically.  The "cathode" is now an
anode with a 208 V collecting sheath, which ionises hard, so I_cond never turns
around.  There is a POSITIVE POTENTIAL MAXIMUM at x/L = 0.022 (+247.6 V) with
field reversal on its anode side (the E/N notches at x/L ~ 0.02-0.04) -- a
double layer that traps electrons from both sides.

Amplitude check: a slab of half-width L with net charge fraction f gives
dphi = f*(Te/e)*(L/lambda_D)^2/2 = 0.12*53.7*(2.89)^2/2 = 26.9 V, against the
39 V measured from the cathode cell to the hump peak.  So the hump is
physically consistent RELATIVE TO THE CATHODE; the +248 V absolute offset is
entirely the external circuit's overshoot.

SUPERSEDES the earlier reading that "j came within 3% of Grubert": that match
occurred during the spike transit with a TOWNSEND structure (both density peaks
on the anode side, x/L = 0.72-0.93, 1-2 decades BELOW Grubert), not a glow
structure.  See the two earliest snapshots in `profiles/`.

---

# ANALYSIS 2026-09-07: the reversals, the mesh, and the ordering

Prompted by the user asking whether the localised field reversals and the "dive"
in the log E/N plot are a MESH problem.

## The E/N dives ARE field reversals -- and their DEPTH is meaningless

Both dives sit exactly at potential extrema, computed from the signed field
`Ex = -dphi/dx` (the CSV carries `Emag`, a magnitude, so the sign has to come
from the potential):

| dive | x/L | extremum | phi |
|---|---|---|---|
| 1 | 0.0216 | potential **maximum** | **+247.6 V** |
| 2 | 0.7210 | potential **minimum** | -8.4 V |

`|E|` must pass through ZERO at an extremum, so on a log axis every reversal
reads as a bottomless spike. The measured floors -- 23.8 Td and 1.8 Td against
thousands of Td two cells away -- are simply **how close the nearest cell centre
landed to the exact crossing**. Read the dive POSITIONS as physics and ignore
their DEPTHS entirely.

## NOT THE MESH. Measured, per local Debye length

| reversal | dx | lambda_D | cells per lambda_D |
|---|---|---|---|
| x/L 0.0216 | 6.77 um | 34.6 um | **5.1** |
| x/L 0.7210 | 54.3 um | 636 um | **11.7** |

Both resolved. And across the eight snapshots the reversals move SMOOTHLY and
MONOTONICALLY (0.816 -> 0.781 -> 0.752 -> 0.732 -> 0.721) rather than hopping
between cells, which is the opposite of an under-resolved wiggle. The 400 CSV
points ARE the cell centres (400 along x x 5 extruded), so there is no
resampling artefact either.

## On the densities: 1000x, but they are a SYMPTOM, not the fault

`n_e` 2.482e18 against a 2.478e15 reference is 1002x. But `n_e ~ n_Arp` to ~1%
(quasi-neutral) and `j` is 100x the reference, so **the density is the correct
density for the current the circuit pushed**. It is a right answer to a wrong
boundary condition, not a numerical blow-up.

## THE ORDERING RESULT -- the double layer PRECEDES the electrode flip

| t [ns] | electrode | phi_max | reversals x/L | n_e max |
|---|---|---|---|---|
| 9999.9 | -98.13 V | +0.03 @ 0.9732 | 0.9732 | 3.48e13 |
| 10325.7 | **-90.12 V** | **+8.17 @ 0.5827** | 0.5827 | 4.62e14 |
| 10426.6 | **+24.78 V** | +43.9 @ 0.0404 | 0.0416, 0.8158 | 1.11e17 |
| 10429.1 | +208.79 V | +247.6 @ 0.0216 | 0.0216, 0.7210 | 2.48e18 |

At 10325.7 ns the electrode is still at **-90 V, the correct polarity**, and an
interior potential maximum of +8.17 V has ALREADY formed at mid-gap. The
maximum then marches onto the electrode (0.58 -> 0.04 -> 0.022) while the second
reversal migrates in from the anode. Growth rate accelerated **480x** over the
sequence, 5.25e6 -> 2.52e9 1/s (tau 191 ns -> 0.40 ns), ending at 0.38x the
absolute ionisation ceiling.

## REFUTED, and it was MY hypothesis: the early state is NOT over-producing

Proposed cause: the ionisation/loss balance is already tilted toward net
production before anything visibly fails. **It is not.** Townsend integral over
the gap with our own `alphaN_vs_reducedE`, self-sustaining at
`alpha*d = ln(1+1/gamma) = 2.872`:

| t [ns] | electrode | alpha*d | multiplication |
|---|---|---|---|
| 9999.9 | -98.13 V | 2.267 | **0.519** |
| 10049.4 | -98.26 V | 2.286 | **0.530** |
| 10325.7 | -90.12 V | 2.460 | **0.643** |
| 10429.0 | +158.20 V | 2.872 | 0.9999 |
| 10429.1 | +208.79 V | 3.099 | 1.27 |

The healthy, normal-polarity snapshots are **SUB-CRITICAL (0.52-0.64)** -- the
gap is not self-sustaining there and should be DECAYING by gas-phase ionisation
alone. So the runaway is NOT caused by an inherently over-producing coefficient
set, and that hypothesis is dead.

CAVEAT on the late rows: once the electrode is positive the structure is no
longer a cathode fall, so a 1-D `alpha*d` over the gap is not strictly
interpretable there. The EARLY rows are the meaningful ones, and they are the
ones that refute the hypothesis.

## What the evidence now supports

`n_e` GREW 13x (3.48e13 -> 4.62e14) over 10000 -> 10326 ns while multiplication
was only 0.52-0.64. Sub-critical ionisation cannot do that, so the accumulation
is not a local avalanche. The candidate left standing is **slow ION
ACCUMULATION**: the ion transit is 3.7-6.7 us and the run had been going ~10 us,
i.e. ~2 transits, so ions produced by even sub-critical ionisation build
positive space charge faster than they drain. That distorts the field until it
locally exceeds breakdown -- and the load line then has too little current
headroom to hold the electrode.

NOT YET TESTED. The discriminating measurement is the ion continuity balance
per region: production `INT S_iz dx` against wall/drift removal, on the EARLY
snapshots. If accumulation is the cause it shows there, at -98 V, with no
runaway in sight.
