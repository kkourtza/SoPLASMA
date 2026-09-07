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
