# grubert2009_steady -- phase 1 of the steady route

Created 2026-09-07. Plan: `doc/steady-and-stability-design.md` sections 4 and 5.

## The question

**Is Grubert's published profile a fixed point of our model?**

Every earlier attempt asked the opposite question -- can the transient solver
REACH the operating point -- and the measured answer is no, for a reason now
established four independent ways (the `C_gap/g` bound; the current/voltage
staircase falsification; the seeding circularity; and the structural
`tau_RC >> tau_ion` result in `../grubert2009_spike/COMPARE.md`). This case
starts AT the destination and asks only whether it HOLDS.

Phase 1 needs **no new solver and no new circuit code**: setting `V_target = V`
in `seriesResistor`'s update cancels `R*g` and `a = tau/dt` identically, leaving
`V_gap = V_src - R*I_cond` -- the algebraic steady load line. So the ballast we
already have IS a damped Newton iteration onto the steady state, and those
terms only set the rate of approach.

## Baselines -- ABSOLUTE PATHS

    /home/kkourtza/soplasma-scratch/validation/grubert2009_fix_n11
        the healthiest earlier ballast arm (R = 1e8). Structural
        discriminators for judging drift are recorded in its COMPARE.md.

    /home/kkourtza/soplasma-scratch/validation/grubert2009_spike
        the FAILURE mode this case is designed to avoid: electrode slewed
        -500 -> +208 V, densities 400-1000x reference, j 99.6x.

    /home/kkourtza/soplasma-scratch/validation/grubert2009_seed
        the same digitisation used as a TRANSIENT seed, which failed. Its
        no-run result (394 V from current continuity) still stands.

**NOT A VALID CONTROL against any of them for absolute agreement**, and this
must not be forgotten when the numbers are read: this case changes THREE things
at once relative to `fix_n11` -- seeded initial state (vs floor), `R = 3e8`
(vs 1e8), and `meanEnergyMax` now DERIVED at 2644 eV (vs a hardcoded 100 eV).
It is a test of a HYPOTHESIS, not a one-variable comparison. The one-variable
control for the clamp alone is `../grubert2009_fix_n11_clamp`.

## Reference (Grubert, Becker & Loffhagen, Phys. Rev. E 80 (2009) 036405)

Steady normal glow, 100 Pa argon, 1 cm gap, gamma = 0.06, V = 500 V:

| quantity | reference |
|---|---|
| n_e peak | 2.478e15 m^-3 at x/L = 0.394 |
| n_Ar+ peak | 6.391e15 m^-3 at x/L = 0.203 |
| j | 0.511 mA/cm^2 (1.022e-6 A on this 2.0e-3 cm^2 electrode) |
| gap voltage | -500 V |

## The circuit, and why R = 3e8

`sourceVoltage = V_gap + R*I_set = -500 + 3e8*(-1.022e-6) = -806.6 V`, so the
load line passes through Grubert's operating point exactly. `R = 3e8` because a
300 V ballast drop at `I_op` needs 2.94e8; the earlier arms' 1e8 is 3x faster
than the operating point demands. Stability: `R > |dV/dI| = 1.7e6 Ohm`
(measured 2026-09-06) with a 176x margin, and `tools/glow_stability.py` reports
this circuit STABLE at `gamma = 0` independently.

NOT RAMPED. The seeded plasma conducts ~1 uA from step one, so the `R*I` drop
appears immediately and the gap sits near -500 V. Ramping from 0 would hold the
gap near 0 V and let the seeded profile decay -- destroying the initial
condition.

## The seed, and its no-run gate

    ~/ct-env/bin/python seed_rho0.py 0

Quasi-neutral bulk, digitised fall, blended over x/L = 0.35..0.45. The gate is
in the script and it REFUSES to seed if the bulk could avalanche. Measured:

| quantity | value | note |
|---|---|---|
| n_Arp peak | 6.391e15 at x/L = 0.2053 | reference 6.391e15 at 0.203 |
| n_e peak | 2.585e15 at x/L = 0.3979 | reference 2.478e15 at 0.394 (4.3% high, from the blend) |
| E at cathode | 13206 Td | |
| fall carries | 507.2 V of 500 (101.4%) | so the bulk field is slightly REVERSED |
| bulk E/N | 57.5 Td | target was 12-30; see the limit below |
| 1/nu_i at 57.5 Td | 4.84e-6 s vs 6e-9 s transit | **806x margin -- cannot avalanche** |

**THE IRREDUCIBLE LIMIT.** The bulk field is a ~57 Td residual of a 13,277 Td
swing -- 0.4% -- so it is the differencing problem ONE LEVEL UP: not
`n_Arp - n_e` pointwise but `E0 - INT rho`. It carries ~+-80 Td that no care in
the construction can remove. That is ACCEPTABLE because the seed only needs to
be in the BASIN, and across the whole band nothing avalanches. Judge the seed on
the avalanche margin, NOT on hitting 12-30 Td.

## Discriminating observables, and when they become visible

| observable | holds | drifts |
|---|---|---|
| n_e peak magnitude / location | stays ~2.5e15 near x/L 0.39 | moves or grows decades |
| j | settles near 1.022e-6 A | climbs past it, as in every earlier arm |
| V_electrode | sits near -500 V | slews toward 0 and positive (the spike mode) |
| cathode fall | stays localised, x/L < 0.4 | collapses or moves to the anode |

Earliest signal is a few `tau_diel` (sub-ns), so a state that is NOT a fixed
point shows within the first ~100 steps. `endTime 10e-6` is ~1.5-2.7 ion transit
times (3.7-6.7 us), which is what makes a HOLD meaningful rather than merely
slow drift. Stop early if the electrode is clearly slewing.

## Extraction

    ~/ct-env/bin/python extract_centreline.py <timedir>    # profiles vs x
    ~/ct-env/bin/python extract_VI.py                      # V and I vs t
    postProcessing/probes/*/                               # rule 26 time series
    postProcessing/externalCircuit/circuit.csv             # V_src, I_cond, V, g
    postProcessing/dischargeCurrent/current.csv

## Known limitations of this case, stated up front

1. `meanEnergyMax` is DERIVED (2644 eV) here and was a hardcoded 100 eV in every
   earlier run. See `../grubert2009_fix_n11_clamp` for the isolated control.
2. Rule 26 is only PARTLY met: `probes` covers the per-probe field series, but
   dt-with-its-limiter, outer iterations and rejected steps are still log-only,
   because the framework's `plasmaProbes`/`history.csv` is design-only.
3. `driftDiffusionFluxScheme standard`, not ScharfetterGummel -- SG is fatal on
   the Hagelaar wall-flux conditions (`ac1bc6a`). The standard branch has a
   documented singularity mode at the anode; if this run dies there, that is the
   first suspect.
