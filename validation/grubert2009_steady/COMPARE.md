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

---

# RESULT 2026-09-07: DIVERGED in 6.5 ns. Two errors of mine, both findable before launch.

## What happened

| quantity | seeded | at t = 6.46 ns | reference |
|---|---|---|---|
| n_e max | 2.585e15 | **1.378e19** (5330x) | 2.478e15 |
| meanE max | <= 100 eV | **2644 eV -- PINNED at the table top** | -- |
| V_electrode | ~-500 V | **+7041 V** | -500 V |
| I_cond | ~1e-6 A | **-2.69e-2 A** (26 mA) | 1.022e-6 A |
| dt | 1e-12 | 3.65e-15, with retries | -- |

Probes: near-cathode `n_e` 9.31e18, `meanE` 1937 eV, `chargeDensity`
**-0.1177 C/m^3 -- NEGATIVE**, i.e. electron-dominated where a cathode fall must
be ion-dominated. Mid-gap `meanE` 1471 eV. The clamp fired at cells with
`n_e = 7.4e17`, 900x the domain minimum, so this is NOT the empty-cell
pathology -- it is a genuinely hot, dense plasma. Same terminal mode as
`../grubert2009_spike` (electrode driven positive, double layer, no recovery)
reached ~1500x faster.

## ERROR 1 -- the phase 0 gate tested the wrong region

`seed_rho0.py` gates on the avalanche margin at the **BULK** field: 57.5 Td,
`1/nu_i` 4.84e-6 s against a 6e-9 s transit, 806x, "PASS". It never tested the
**FALL**, where the seeded field is **13,206 Td** -- 230x higher and about four
decades faster in ionisation. The gate checked the region where the answer was
reassuring. This is the recurring failure in [[never-read-a-number-without-its-control]]:
a quantity chosen where the error is not.

**The correct gate is the TOWNSEND CRITERION**, and it convicts the seed with no
run. Self-sustaining needs `alpha*d = ln(1 + 1/gamma) = 2.872` at gamma = 0.06.
Measured on the seeded profile with our own `alphaN_vs_reducedE`:

    alpha*d = 3.146   =>  net multiplication gamma*(exp(alpha d) - 1) = 1.33

33% over unity -- a net GROWTH condition. And the excess is in the FALL, not the
adjustable part: of the 3.146, **3.020 is accumulated by x/L = 0.30**, against a
2.872 requirement, so **the fall alone over-satisfies the criterion**.

**No gap voltage fixes it.** `alpha*d` has a MINIMUM near -500 V and rises in
BOTH directions (3.146 at -500, 3.319 at -480, 4.682 at -400), because lowering
|V| drives the bulk residual field more positive and its |E| larger. There is no
self-sustaining root for this profile.

  RETRACTED IN THE SAME BREATH: a bisection here reported "self-sustaining at
  Vgap = 520.0 V, 1.040x Grubert" -- that is a BRACKET ENDPOINT, not a root
  (`alpha*d` there is 3.17, not 2.872), and its apparent agreement with the
  519.4 V from the field-free-bulk construction is COINCIDENCE. Recorded so the
  wrong number is recognisable if it surfaced anywhere. Rule 22.

## ERROR 2 -- R = 3e8 HALVED the current headroom, and the formula was already known

`I_sc/I_op = 1 + V_gap/(R*I_op)` is recorded in memory. Applied:

| R | V_src | I_sc | headroom |
|---|---|---|---|
| 1e8 (the old arms) | -602 | 6.02e-6 | **5.89x** |
| **3e8 (this case)** | -806.6 | 2.69e-6 | **2.63x** |
| 5e8 | -1011 | 2.02e-6 | 1.98x |

Raising R to 3e8 made the ballast stiffer in VOLTAGE and cut the CURRENT
headroom from 5.89x to 2.63x. The electrode flips as soon as `I > I_sc`, so a
33%-over-multiplying seed only needs to reach 2.69 uA -- nanoseconds -- and then
it is the unrecoverable positive-electrode mode. **The trade-off is intrinsic to
a two-terminal ballast: stiffness in voltage IS loss of headroom in current**,
and `currentSource` (R -> inf, headroom 1.0x) is the extreme, not the cure.

I chose R = 3e8 from the ballast-drop argument alone and did not compute
`I_sc/I_op` before launching, though the formula was in memory.

## What this does and does not settle

SETTLED: this seed is not in the basin, and the reason is quantitative and
pre-computable -- 1.33x multiplication concentrated in the fall.
NOT SETTLED: whether Grubert's profile is a fixed point. The construction never
got close enough to ask. A seed must satisfy `alpha*d = 2.872` to be a candidate
at all, and that is now the gate.

## Also flagged: `2450 Td` is an UNSOURCED number

"cathode fall E/N ~ 2450 Td" appears in `../grubert2009_spike/COMPARE.md`,
`../grubert2009_fix_fast/COMPARE.md` (scoring a run at "0.81x") and
`../grubert2009_spike/plot_profiles.py`, always labelled "(normal fall)". It is
in NONE of the digitised Grubert material or its READMEs. G2 forbids exactly
this -- a number somebody remembered, used as a reference. The seeded fall here
needs ~13,000 Td to carry 500 V over this geometry, 5.4x that figure, so the two
cannot both describe this case.
