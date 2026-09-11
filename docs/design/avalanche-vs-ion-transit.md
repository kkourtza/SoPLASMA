# Why no cathode fall has ever formed: the avalanche outruns the ions

2026-09-08. This is the first mechanistic account of the Grubert 2009 DC-glow
failure that is supported by a measured causal chain rather than by elimination.
Every prior explanation was a candidate ruled out; this one has numbers.

## 1. The measurement

`validation/grubert2009_ps_knee` -- CompleteFlux, BDF2, fixed dt = 1 ns,
`adjustTimeStep false`, no `relaxationFactors`, current-driven cathode with
`I_set` crawling 5.0e-8 -> 5.4e-8 A. It ran 1264 consecutive converged steps
and then failed at t = 1.267 us with two `PIMPLE: not converged within 150
iterations` followed by a SIGFPE in `plasmaEnergy::solveSpeciesEnergy`.

The circuit trace (`postProcessing/externalCircuit/circuit.csv`):

| t [us] | I_set [A]  | I_cond [A] | V [V]  | I_cond/I_set |
|--------|------------|------------|--------|--------------|
| 0.9    | -5.282e-8  | -3.765e-8  | -199.5 | 0.71         |
| 1.0    | -5.323e-8  | -5.893e-8  | -202.3 | 1.11         |
| 1.1    | -5.363e-8  | -7.748e-8  | -193.1 | 1.44         |
| 1.2    | -5.403e-8  | -1.707e-7  | -159.3 | 3.16         |
| 1.264  | -5.429e-8  | -3.155e-6  |   -6.9 | 58.1         |
| 1.266  | -5.430e-8  | -7.180e-6  |    0.0 | 132.2        |

and the diagnostics over the same window:

| t [us] | n_e,max [m^-3] | Te,max [eV] | Te,min [eV] | d(ln n_e)/dt [1/s] |
|--------|----------------|-------------|-------------|--------------------|
| 0.325  | 1.98e11        | 5.57        | 1.58        | 3.8e6              |
| 0.625  | 1.38e12        | 7.22        | 1.74        | 6.6e6              |
| 0.925  | 9.78e12        | 8.51        | 1.85        | 6.7e6              |
| 1.125  | 1.80e14        | 10.20       | 1.92        | 2.1e7              |
| 1.250  | 3.95e15        | 9.62        | 3.42        | 3.0e7              |

## 2. What it says

**The regulator lost authority; it was not out-run by the ramp.** The source
drove the gap voltage to EXACTLY ZERO while the conduction current rose to 132x
its setpoint. A secant-Newton regulator on I(V) cannot control this, and no
retuning would help: lowering V halts FURTHER ionisation but does not remove
the 4e15 m^-3 of electrons already present. The current is n_e*mu*E, so pinning
it once n_e has run away demands E -> 0, which is what the trace shows.

**The ramp rate is exonerated by direct test.** A deliberately 24x slower
current ramp (`grubert2009_ps_slow`) bought 0.05 us: failure moved from 0.962
to 1.013 us. The knee crawl, 25x slower again THROUGH the ignition knee,
reached 1.267 us. Three ramp rates spanning ~600x, three near-identical deaths.

**The timescales do not permit a cathode fall.** The avalanche e-folding time is
1/6.7e6 = 149 ns, and it SHORTENS to 33 ns as the runaway accelerates. The ion
transit time across the 1 cm gap is 4-7 us. So the electrons multiply 30-45x
faster than the ions can drift anywhere. A cathode fall IS a space-charge
structure -- it exists only once ions have accumulated near the cathode -- so
on this separation of timescales it cannot form, in this or any run with the
same separation. **This is why no fall has appeared in any case to date.**

`Te,min` -- the BULK temperature, which should stay near 1 eV -- rising from
2.00 to 3.42 eV in the last 25 ns is the same fact seen from the field side:
the voltage is spread across the whole gap, heating everywhere, rather than
concentrated in a fall.

## 3. The trigger

By t = 0.9 us the gap already held **-199.5 V** with n_e still ~1e13 m^-3, i.e.
with no plasma to speak of. That is not the drive being too strong; it is the
pre-ignition charging ramp:

    C_ext = 0  ->  C = C_gap = 1.77e-16 F
    dV/dt = I_set/C = 5e-8 / 1.77e-16 = 2.8e8 V/s = 282 V/us

so the gap reaches ~200 V in under a microsecond, **uniformly**: 200 V / 0.01 m
against N = 2.4e22 m^-3 is ~83 Td ACROSS THE ENTIRE GAP, which is strongly
ionising everywhere. The whole gap avalanches simultaneously, and a structure
that requires ions to sort themselves out near one electrode never has a
chance. Measured breakdown for this gap is ~-121 V, so 200 V is a ~65%
overshoot -- reached because a current source with no plasma to regulate
against does the only thing it can: charge the electrode capacitance.

## 4. The falsified counter-argument

`capacitance 0` was previously chosen deliberately, and the reasoning was
sound arithmetic:

> With C_gap alone (1.77e-16 F) an excess conduction current discharges the gap
> in 0.286 ns -- FASTER than the 0.68 ns ionisation growth, so the source can
> quench a fast excursion. Adding the 5e-14 F shunt the ballasted arms used
> makes that 80.9 ns, 119x too slow, and is exactly why they could not hold the
> discharge. The cost is a fast pre-ignition ramp, but that is harmless
> BECAUSE THE CURRENT IS PINNED.

The current was **not** pinned: I_cond/I_set reached 132. The arithmetic was
right and the conclusion wrong, because the 0.68 ns growth time is the growth
rate *at ~200 V across a uniform gap* -- the very state the fast ramp creates.
The argument assumed that state and asked whether the source could survive it.
It cannot, and neither can any capacitance, because quenching requires removing
electrons and lowering V does not remove them.

## 5. What follows

Attack the STATE, not the response. `plasmaExternalCircuit.C:208` already says
what a real rig has: "far more capacitance across the gap than the gap".
`capacitance` adds to C_gap and sets the pre-ignition ramp:

| C_ext [F] | x C_gap | dV/dt [V/us] | reaches -121 V at | vs ion transit |
|-----------|---------|--------------|-------------------|----------------|
| 0         | 1       | 282          | 0.43 us           | 10-16x FASTER  |
| 2.5e-15   | 14      | 20           | 6.1 us            | ~1 transit     |
| 5e-15     | 28      | 10           | 12.1 us           | 2-3 transits   |

Ignition then happens NEAR BREAKDOWN, where the net avalanche rate is ~0 by
definition and the growth time is microseconds rather than nanoseconds, so
space charge has time to build a fall. The quench-time objection of section 4
weakens by the same factor, and for the same reason: it is only binding in the
200 V uniform state this avoids.

Running as the control pair `validation/grubert2009_ps_C25e15` and
`validation/grubert2009_ps_C5e15` (one variable, `capacitance`).

**Falsifiable prediction.** A cathode fall forms -- a voltage drop concentrated
within ~1 mm of the cathode, with `Te,min` staying near 1 eV. If both arms
still run away, the drive is exonerated and the defect lies in the wall
emission / energy boundary physics rather than in how the gap is charged.

---

## 6. OUTCOME: the prediction in section 5 FAILED, and the circuit is exonerated

Both capacitance arms ran away and died. So did a ballast arm run afterwards.

| arm | circuit | died at | terminal state |
|-----|---------|---------|----------------|
| `_ps_knee`    | currentSource, C_ext = 0    | 1.267 us  | V -> 0, I_cond/I_set = 132 |
| `_ps_C25e15`  | currentSource, C = 2.5e-15  | 8.192 us  | ratio 0.85 -> 5.15 in 200 ns |
| `_ps_C5e15`   | currentSource, C = 5e-15    | 13.684 us | same runaway |
| `_ps_ballast` | seriesResistor, R = 1e9     | 4.342 us  | **V_gap -> +0.2 V** |

The slow ramp DID do what section 5 said it would, and that part holds:

  * 7841 / 9622 / 3521 consecutive converged steps, 12-14 correctors against a
    150 cap that used to saturate at 150/150.
  * n_e growth 8x slower (7.9e5 vs 6.5e6 /s).
  * Positive space charge of 200x (n_Arp 2.1e13 vs n_e 1e11).
  * Field concentrating at the cathode: E_cath/E_anode 1.24 -> 1.57 -> 2.63 in
    the capacitance arm, and **4.19** under the ballast -- the furthest an
    incipient cathode fall has ever developed here.

But no arm reached a field REVERSAL, and all three died the same way. **Three
independent circuit topologies with the same terminal signature means the
defect is not in the circuit.** The capacitance trade-off is real
(`I_set - I_cond = C dV/dt`, so the C that tames the ramp is what stops V
responding) and the ballast does have algebraic authority -- it held the
current e-folding at ~1200 ns, 11x slower than the capacitance arm's 108 ns --
but it too lost the gap at ignition.

## 7. THE ACTUAL DEFECT: the runs were driven 20x below the operating point

The cathode patch is **2.0e-3 cm^2** (y in 0..1 mm, z in 0..0.2 mm). Grubert's
j = 0.511 mA/cm^2 therefore corresponds to

    I_op = 0.511e-3 A/cm^2 * 2.0e-3 cm^2 = 1.022e-6 A

Every run above was driven at ~5e-8 A -- **20x too low**. This was a sizing
error on my part, not a physics discovery.

It is not merely "too small", and this is the point. The cathode electron
density that gamma-SEE sustains is
n_e = gamma*j / (e*(1+gamma)*v_drift), so at v_d ~ 3e5 m/s:

| I [A]    | j [mA/cm^2] | n_e(cathode) from gamma-SEE | vs the 1e11 floor |
|----------|-------------|------------------------------|-------------------|
| 1.63e-8  | 0.008       | 9.6e10                       | **the floor wins** |
| 1.0e-7   | 0.050       | 5.9e11                       | 6x                 |
| 1.022e-6 | 0.511       | 6.0e12                       | 60x                |

At 5e-8 A the `minNumberDensity 1e11` clamp is THE SAME SIZE as the real
cathode electron source. The clamp is applied cell-by-cell across the whole
domain, so in that regime **the entire gap acts as a weak distributed cathode**
and ionisation proceeds volumetrically everywhere rather than closing the
cathode-fall/gamma loop. Every run has been in exactly that regime.

The profiles say so directly: n_e sits at EXACTLY 1.000e11 through the whole
cathode region at every output time, in every arm -- the solution wants to be
lower there and is being held up. Meanwhile n_Arp is free at 2e13, which is
why space charge and field concentration appear at all.

Note this also retires an earlier hypothesis of mine: the floors are NOT
enforcing quasineutrality (n_Arp/n_e ~ 200 in the cathode region, so space
charge is real). The floor's damage is as a distributed SOURCE, not as a
constraint on charge separation.

### The test

`validation/grubert2009_ps_jop`: ballast R = 5e8 Ohm with V_src ramping to
-1100 V, so I_sc = 2.2e-6 A **brackets** I_op = 1.022e-6 A. The ramp sweeps the
load line across the operating point instead of guessing it -- j-continuation
carried out by the ballast -- and the point where V_gap settles near Grubert's
-500 V is the measurement. Leg 1 crosses to -100 V in 1 us (below the -121 V
breakdown); leg 2 creeps at 25 V/us for 40 us, i.e. 6-10 ion transits.

The floor is left ALONE at 1e11 in that arm, so it is one variable.

## 8. CORRECTION: raising the current is NECESSARY BUT NOT SUFFICIENT

Section 7 claimed the floor is "60x below the physical density and inert" at
I_op. **That was wrong**, and the error was mine: it used an electron drift
velocity of 3e5 m/s. The solution's own mobility is mu_e = 92.8 m^2/V/s, so
mu_e*E is ~1e7 m/s -- two orders higher -- and the sustained density is
correspondingly two orders lower.

The right test is field- and current-INDEPENDENT. At a cathode surface the
gamma-SEE closure is Gamma_e = gamma*Gamma_i with both fluxes pure drift in the
same field, hence

    n_e/n_Arp = gamma * mu_Arp / mu_e

Measured on `_ps_ballast` at t = 3.5 us, using the solution's own mobilities
(mu_e = 92.8, mu_Arp = 0.0695, ratio 1335):

| quantity | value |
|----------|-------|
| predicted `n_e/n_Arp` | 4.50e-5 |
| **observed** `n_e/n_Arp` | **5.44e-3** |
| observed / predicted | **121x** |
| n_e that gamma-SEE alone sustains | 8.27e8 m^-3 |
| n_e actually present | 1.00e11 m^-3 (**the floor**) |

So the cathode electron density is set by the CLAMP, not by the emission
physics, and the Townsend loop is swamped 121-fold.

The floor only goes inert once n_Arp > 1e11 / 4.50e-5 = **2.2e15 m^-3**, which
needs roughly 2e-6 A -- *above* Grubert's operating point of 1.022e-6 A. At
I_op itself gamma-SEE sustains only ~5e10 m^-3, still BELOW a 1e11 floor.
**The floor has to come down as well.**

### The second test

`validation/grubert2009_ps_floor`: same ballast as `_ps_jop`, with

  * `minNumberDensity 1e4` for `e` and `Arp` -- seven orders below the physical
    cathode density, a positivity guard for the unlimited flux scheme and
    nothing more;
  * `0/n_e` and `0/n_Arp` seeded EXPLICITLY at `uniform 1e11`.

The seed had to be made explicit first. `plasmaSpeciesProperties` stated the
conflation outright -- "THE 1e11 FLOOR IS THE INITIAL CONDITION" -- and
`0/n_e` was `uniform 0`, so lowering the clamp without seeding would have
started the case from vacuum. Separating the two roles is what lets the floor
be chosen on positivity grounds alone, and it is a framework fix independent of
this case.

## 9. The inconsistent initial state, and why it violated positivity on step 1

Measured across every arm: **the first negative-density clip occurs at
t = 1e-9 s -- the FIRST timestep** -- with 564 to 3050 cumulative clips before
death.

| arm | died | first clip | clip events | cumulative |
|-----|------|-----------|-------------|------------|
| `_ps_knee`    | 1.267 us | **1e-9 s** | 79 | 1326 |
| `_ps_C25e15`  | 8.192 us | **1e-9 s** |  3 |  564 |
| `_ps_ballast` | 4.342 us | **1e-9 s** | 22 | 3050 |
| `_ps_jop`     | 3.742 us | **1e-9 s** |  3 |  564 |
| `_ps_floor`   | (running) | **none at 6 us** | **0** | **0** |

The cause is the seed/floor conflation. `0/n_e` was `internalField uniform 0`
and `minNumberDensity 1e11` supplied the background, so the run STARTED from a
state the fields and the clamp disagreed about: zero density in the file, 1e11
after clamping. CompleteFlux carries no limiter -- the warning says so
outright, "POSITIVITY IS VIOLATED. The flux scheme in use carries no limiter,
so this is possible by construction" -- and it went negative immediately.

Each clip then injected up to the floor. In `_ps_knee` one event reads
`min = -18.19 in 238 cell(s), floor is 1e11`: 238 cells raised from -18 to
**1e11**, repeated across 79 events. That is a large, repeated, spatially
distributed electron source with no physical origin.

So the chain is:

    0/n_e = 0 with a 1e11 clamp   (inconsistent initial state)
      -> CFS goes negative on step 1  (unlimited scheme)
      -> clips inject 1e11 per affected cell
      -> a distributed volumetric electron source
      -> volumetric avalanche everywhere, no localisation
      -> no cathode fall, and a conductance no circuit can hold

Seeding `0/n_e` and `0/n_Arp` explicitly at `uniform 1e11` and demoting the
clamp to `1e4` breaks it at the first link: **zero clips in 6 us**, against a
violation on step 1 in every previous arm.

This is a framework-level fix, independent of the Grubert case. It also bears
on the pending "should CFS be the default" question: CFS's unboundedness was
being triggered on step 1 by an inconsistent initial state, and the resulting
clips were amplified by a floor set 7 orders too high. The scheme deserves
re-evaluation on a consistent initial state before its boundedness is judged.

## 10. `_ps_floor` died differently: a fixed dt too coarse for the mesh, not physics

`_ps_floor` reached t = 8.871 us with **zero clips throughout** -- the seed/floor
fix held. It then failed the OTHER way: 10 consecutive steps at 150/150
correctors, the outer-loop cap the pseudoSteady recipe relies on to fail loudly
instead of silently degrading.

The diagnostics at failure show why, and it is unrelated to the seed/floor
work:

    Co_conv (e):        2073.86
    Co_conv (energy):   3107.95
    Energy relax. ratio: 0.625      (target 0.2)

`limitSpeciesCo`/`limitEnergyRelaxation` are configured (`maxSpeciesCo 15`,
`maxEnergyRelaxationRatio 0.2`) but, like every limiter, are ADVISORY ONLY
under `adjustTimeStep false` -- which pseudoSteady forces on by design (section
2.1). So the fixed dt = 1e-9 s ran at 138x its own configured Courant limit.

Checking the mesh confirms it is a resolution mismatch, not a scheme problem:
the finest cell is dx = 1.35e-6 m (400 planes across the 1 cm gap, 49x
stretch), and at the pre-ignition field (E ~ 29500 V/m, mu_e ~ 92.8 m^2/V/s)
v_drift = 2.74e6 m/s gives Co(min cell) = 2023 at dt = 1e-9 s. Co = 15 needs
dt = 7.4e-12 s -- **135x smaller** than the value carried over unchanged from
the earlier arms, which used the same mesh at currents ~20x lower and never
drove the field this high before dying for other reasons.

So `deltaT` in the pseudoSteady recipe is not "any small fixed number": it must
be sized against Co_conv on the FINEST cell at the field this specific case
will reach, which for an ignited glow is well past the pre-ignition value the
earlier arms died before ever exercising. This is a case-setting error, not a
finding about the physics -- unlike sections 6-9, which stand.
