# V4 part B — the Grubert dc glow, LMEA against LFA

Created 2026-09-05. Relaunched 2026-09-06 after the FIRST RUN PRODUCED AN
ARTEFACT. **Running; no result yet.**

## The first run was invalid — read this before trusting any number here

The LFA arm ran to completion and reported a current density constant to
**0.00% over 36 us**, which is exactly what a converged steady state looks
like. It was not one. `tableKey meanE` had been inherited from the LMEA
streamer bed this case descends from, so under LFA the `reducedE` field
(V m^2, ~2e-18) indexed tables tabulated against mean energy (eV, 0-2644).
Every lookup hit the table floor, `S_iz` was 1.9e-210 against an expected
1.1e+20, and the discharge never ignited: n_e stayed at 1.53e11 m^-3 against
the 1e11 floor, where Grubert's LFA peak is 6.2e13.

**The plateau check specified in this file would have PASSED on that run.** A
flat current is necessary for steady state and nowhere near sufficient. Check
that `n_e` has risen well above `minNumberDensity`, and that `S_iz` is of the
order `k_io * n_e * N`, BEFORE reading a current.

Now guarded: indexing a field into a table of a different physical quantity is
fatal (`plasmaTransport.C`, committed 2026-09-06), and `tableKey` is deleted
from both arms so it is derived.

Reference: Grubert, Becker & Loffhagen, *Phys. Rev. E* **80** (2009) 036405.
Argon, **100 Pa**, cathode **−500 V**, **1 cm** gap, steady-state dc,
γ = 0.06, electron reflection 0.36.

## The question

Does the LFA/LMEA closure choice reproduce the factor-of-5.9 spread in
discharge current that Grubert measured on this case?

Part A (`SoEEDF/validation/dias2025/COMPARE_v4.md`) is the CONTROL for this:
our Boltzmann coefficients match their multiterm ones to better than 7% below
5000 Td, so a difference here is attributable to the **closure** rather than to
different input data. Do not read this case's result without that one.

## The cases

| arm | path |
|---|---|
| LMEA | `/home/kkourtza/soplasma-scratch/validation/grubert2009` |
| LFA | `/home/kkourtza/soplasma-scratch/validation/grubert2009_LFA` |

**THE CONTROL CHANGED 2026-09-06, and this is the important paragraph.**

The two arms no longer differ in exactly one key, and they must not, because
GRUBERT'S TWO CLOSURES DO NOT. Sec. III A: "In the description using the LFA,
the transport coefficients of the electrons at the mean electron energy of
1.5 eV have been used [64]." Their LFA freezes the electron TRANSPORT and lets
only the RATE coefficients follow the local field. So the arms differ by the
closure AS THEY DEFINE IT:

| | LMEA arm | LFA arm |
|---|---|---|
| `electronEnergyModel` | `LMEA` | `LFA` |
| electron transport | tabulated vs mean energy | **constant**, mu_e = 264.55 m^2/Vs, D_e = 672.16 m^2/s |
| rate coefficients | tabulated vs mean energy | tabulated vs `reducedE` |
| everything else | identical | identical |

The two constants are OUR OWN Boltzmann tables evaluated at <U> = 1.5 eV
(muN = 6.387072e+24, DLN = 1.622804e+25 SI, divided by N = 2.4143e22), so the
CLOSURE differs between the arms but the CROSS-SECTIONS do not -- which is what
keeps the comparison about the closure.

This is NOT a licence to differ freely: every other setting is still verified
identical, and the general `electronEnergyModel LFA` still keys transport on
the tabulated `_vs_reducedE` data (the textbook form). Grubert's frozen
transport is a simplification they cite [64] for, expressed per-case, and the
solver announces the override at start-up.

WHY IT MATTERS: with transport keyed on the local field, our LFA arm drove
reducedE to 6.4e7 Td in the sheath -- 1170x past the top of its own tables --
and ran away to I_cond 1.9e7 mA/cm^2 before SIGFPEing inside the chemistry ODE.
Grubert report the opposite failure mode for their LFA, and only at a different
pd: "an ignition of the discharge was not predicted". Ours could not have been
their model, and no mesh or solver tuning would have closed that gap.

**The former control, superseded:** the two trees differ in **exactly one
key** — `electronEnergyModel LMEA` against `LFA` in
`constant/plasmaSpeciesProperties`. Same mesh, same mechanism, same ion
transport tables, same boundary declarations, same numerics. Verified with

```
diff -r --brief grubert2009 grubert2009_LFA | grep -vE "logs|postProcessing|0/|plasmaTables|\.msh|log\."
diff grubert2009/constant/plasmaSpeciesProperties grubert2009_LFA/constant/plasmaSpeciesProperties
```

which must report that one line and nothing else. `chemistry/lookupVariable`
was deleted from both: it is DERIVED from `electronEnergyModel`, and stating
`meanE` under LFA is the half-LMEA the solver rejects.

## Reference numbers to reproduce

Total conduction current density at mid-gap, digitised from their figure 3
(`SoEEDF/validation/dias2025/reference/README_grubert_fig3.md`, all 12 curves,
current-continuity checked to 1.3–2.9%):

| closure | j_tot | relative to LMEA |
|---|---|---|
| LMEA | **0.511 mA/cm²** | 1.00 |
| LMEAC | 0.201 mA/cm² | 0.39 |
| **LFA** | **0.087 mA/cm²** | **0.17** |

We have no LMEAC arm (no excited-state kinetics in `argon_plasma.yaml`), so
LMEAC is context, not a target. Our 3-species mechanism matches their **LMEA**
arm, which also has no excited kinetics.

Profiles to compare against, same directory:
`grubert_fig3_na_{LMEA,LFA}_{Arp,electrons}.csv` (densities, 10⁹ cm⁻³) and
`grubert_fig3_ja_{LMEA,LFA}_{Arp,electrons}.csv` (current densities,
10⁻¹ mA/cm²), both against `z/d`.

## Extraction

```
# current, both arms
tail -1 <case>/postProcessing/dischargeCurrent/current.csv
# columns: time,V_applied,I_total,I_cond,I_disp,I_<species>...

# steady state is judged from the TRACE, not from endTime
awk -F, 'NR>1{print $1","$3}' <case>/postProcessing/dischargeCurrent/current.csv
```

Convert to a current density with the electrode area: the anode/cathode
patches are 5 faces of 1 mm × 0.2 mm = **2.0e-7 m²** total.

## Cost, and why it is a background run

Sized before launch, and then measured:

* ion transit across 1 cm at ~2000 Td: **4.5 µs**; steady state needs several,
  hence `endTime 45e-6`
* **measured 2026-09-06**, after two further changes (below), 60 s of wall
  clock each:
  * **LFA**: 3641 steps to 1.10e-7 s, 40 rejections → ~7 h
  * **LMEA**: 3395 steps to 4.71e-8 s, **2 rejections** → ~16 h

  The LMEA arm previously stalled outright: 31% of steps rejected and 3.8e-9 s
  of progress in 4.5 hours. The cause was physical, not numerical -- a 100 ns
  ramp onto a uniform 1 cm gap at 2071 Td gives alpha*d ~ 930, an explosive
  avalanche the solver had to resolve step by step. Two changes fix it without
  touching the steady state, which is what the reference is:
    1. the ramp is **5 us**, still 9x shorter than endTime;
    2. the initial densities are **n_e = n_Arp = 1e15 m^-3**, quasi-neutral and
       near Grubert's own peaks (LMEA n_e 2.5e15, Ar+ 6.4e15), so the run does
       not have to avalanche four decades. The FLOOR stays at 1e11: this is an
       initial condition, not a source.
  Because the start is deliberately near the expected answer, the steady state
  must be shown to be an ATTRACTOR: rerun the cheap LFA arm from the 1e11 start
  as well and confirm the same plateau, before quoting either.

The LMEA arm is 45x slower in simulated time per wall-second because its outer
loop couples the energy equation; LFA has no energy equation to couple. That
ratio is itself a result worth reporting — it is the CPU cost of the closure,
one of the four axes this project is judged on.

Writes are throttled to **10 snapshots** (`writeInterval 4.5e-6`) because only
the steady state is wanted; the current CSV every 50 steps is what the plateau
is judged from.

## THE INITIAL CONDITION IS PHYSICS HERE, NOT A CONVENIENCE

Measured 2026-09-06, and it cost a day: seeding the domain QUASI-NEUTRAL at
1e15 m^-3 -- done to spare the LMEA arm a violent avalanche under the original
100 ns ramp -- makes the LFA arm RUN AWAY.

A 1e15 plasma shields the bulk and forces the whole applied drop into a sheath
thinner than a cell. E/N spikes there; the LFA keys its ionisation rate on that
local field with no energy relaxation to damp it, and the loop closes: field ->
ionisation -> charge -> sharper field.

| | 1e15 seed | 1e11 seed (used now) |
|---|---|---|
| sim time reached | 1.10e-07 (19045 steps) | **1.22e-06 in 477 steps** |
| rejections | 3022 | **0** |
| Emag max | 8.22e+07 V/m (E/N 3.4e6 Td) | **13183 V/m** (mean is 12200) |
| n_e max | 2.25e+19 from 1e15 | 3.61e+12 from 1e11 |
| dt | ~5e-14 | **~2.6e-09** |
| locality advisories | 3796 | 9 |

The 1e15 arm was 60x past the top of the coefficient tables. So this case runs
from the WEAK seed, and there is no `internalField uniform 1e15` anywhere in
it. The reasoning is recorded beside `minNumberDensity` in
`constant/plasmaSpeciesProperties` so it cannot be "optimised" away again.

**General lesson worth carrying:** starting a transient nearer the expected
answer is a legitimate trick for a steady-state problem, but it changes WHICH
transient is being integrated. Here it replaced a gradual avalanche in a
uniform field with sheath formation in an unresolved one, and only one of the
two closures survives that.

## Numerical settings that had to be changed, each measured

These are NOT free choices: each one was forced by a failure, and each is
identical in both arms so the comparison stays a valid control.

## Known deviations from the reference, each deliberate

1. **Ion reflection 5e-4 is not applied.** The ion wall condition takes no
   reflection coefficient. At 5e-4 it changes the ion wall flux by 0.05%, far
   below the ~5% spread between our Boltzmann coefficients and theirs.
2. **Ar₂⁺ mobility is clamped above 170 Td** (the LXCat set stops there; Ellis
   Parts I–III contain no Ar₂⁺ table). Ar₂⁺ is ~10⁻³ of Ar⁺ here — the 3-body
   conversion rate is ~233 s⁻¹ against a 4.5 µs transit.
3. **The LFA arm uses `wallTeV 1`**, the documented default. The wall flux goes
   as sqrt(Te), so 1 eV against 2 eV is a factor 1.41. If the LFA current is
   off by about that, test this before concluding anything about the closure.
4. **`fluxScheme ScharfetterGummel`, not `standard`.** The standard branch of
   the wall-flux condition went SINGULAR at the anode at t = 1.1e-7 s,
   (D/delta + uEff)/(D/delta) = -9.26. That is not a bug in the closure: with
   reflection r the closure returns W = (1-r)*w_w, so
   uEff = (1-r)/(1+r)*(A + uDrift) - uDrift is negative once
   uDrift > (1-r)/(2r)*A -- only 0.889*A at r = 0.36. A reflecting wall
   accumulates density and the mixed form inverts when the wall value it
   implies is unbounded. SG's denominator, D/delta*Bern(Pe) + uAbs, is
   non-negative for any r, and SG is the right scheme for a drift-dominated
   near-wall cell anyway. `includeDriftFlux true` does NOT help here, despite
   what the error message said before 2026-09-06.
5. **`adaptiveRelaxation true`.** It used to default on only under LMEA, so the
   LFA arm ran undamped: `residual GREW` every step from t = 1.1e-7 s and the
   governor printed "too slow to finish" with every physical limiter at ~3% of
   its cap. 340 rejections became 3. Now the framework default for every energy
   model, so this line is redundant -- kept in the case only until the two arms
   have been re-verified against it.
6. **`maxSpeciesCo` is 15, not 1.5.** Licensed by
   `/home/kkourtza/co-sweep/COMPARE_coSweep.md` (measured 2026-08-31): the
   error SATURATES with Co -- domain maxima within 0.32% (n_e) and 0.51%
   (meanE) for every Co from 3 to 50 -- and that study's own pre-registered
   prediction of degradation beyond Co~5 was RETRACTED. It licenses large Co in
   the inert phase only and only to ~15; both conditions hold here, and the
   guard it names, `Co_chem`, is 0.003 against its 0.9 cap. MEASURED EFFECT on
   the LMEA arm: 4418 steps to reach t = 1.6e-7 where Co=1.5 needed 33546, a
   7.6x speedup against that study's predicted 7.4x, with zero rejections.
7. **`maxCorrectors` is 150, not 60.** Tuned on the LFA arm (`../TUNE_LFA.md`):
   60 rejected 15.9% of steps, 150 rejects 5.1%. Damping harder instead
   (relaxOmegaStart 0.4 / 0.2) reached only 7.6% / 7.3%, so the loop needed
   ITERATIONS rather than more damping -- which is why those knobs were varied
   separately.
8. **`outerCoupling/tolerance` is 1e-6, not the 1e-8 default** — at 1e-8 the
   loop hit its 20-corrector cap and `retryStep` discarded ~24% of steps.
   `maxCorrectors` is 60, not 20. Both are the SAME in both arms.
