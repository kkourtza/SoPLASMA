# HOW FAST CAN THE VOLTAGE RAMP BE? -- a 4-point rate sweep

Created 2026-09-07. The user asked whether the ramp could be even shorter than
5 us, e.g. 1 us. Bracketed with 2 us so a failure at 1 us is immediately
localised instead of needing another round.

## The arms -- ONE variable, the sourceVoltage ramp duration

    grubert2009_floor1e9   table ((0 0) (2e-05 -200))   -10  V/us   (the reference)
    grubert2009_ramp5us    table ((0 0) (5e-06 -200))   -40  V/us
    grubert2009_ramp2us    table ((0 0) (2e-06 -200))  -100  V/us
    grubert2009_ramp1us    table ((0 0) (1e-06 -200))  -200  V/us

All are `cp -r` of `grubert2009_floor1e9`; `diff -rq` confirms the trees differ
in exactly ONE file, `configuration/boundaries`. Same mesh (400 cells, Bump
0.02), same density floor (1e9 m^-3), same circuit (`seriesResistor`, R = 1e8),
same endTime (100 us).

## Pre-flight: nothing else is violated, so this tests ONE thing

| ramp | dV/dt | I_disp = C_gap*dV/dt | I_disp/I_sc | dV per step at dt = 1e-9 |
|---|---|---|---|---|
| 20 us | 1e7 V/s | 1.77e-09 A | 0.0009 | 0.01 V |
| 5 us | 4e7 V/s | 7.08e-09 A | 0.0035 | 0.04 V |
| 2 us | 1e8 V/s | 1.77e-08 A | 0.0089 | 0.10 V |
| 1 us | 2e8 V/s | 3.54e-08 A | **0.0177** | 0.20 V |

Displacement current stays under 2% of `I_sc` = 2e-6 A even at 1 us, and
`maxVoltageRisePerStep` (100 V) is never approached. So the only quantity at
stake is the quasi-static criterion below.

## The criterion under test -- MINE, and already in doubt

`(dV/dt * tau_i)/V`: how much the voltage changes during one ion transit,
relative to itself. I proposed this as the reason not to shorten the ramp.

|  V [V] | tau_i [us] | 20 us | 5 us | 2 us | 1 us |
|---|---|---|---|---|---|
| 50 | 17.12 | 3.42 | 13.69 | 34.23 | **68.46** |
| 100 | 10.70 | 1.07 | 4.28 | 10.70 | **21.39** |
| 200 | 7.16 | 0.36 | 1.43 | 3.58 | **7.16** |

**MEASURED SO FAR (5 us vs 20 us, up to -124 V): the criterion is NOT
predictive.** `Te,max` agreed to within 2% at every voltage and `Te,min` (the
bulk) to two decimals; `n_e,max` agreed to 5-7%, with the FASTER arm slightly
LOWER -- the opposite sign to the prediction. And the 5 us arm reached each
voltage in 3.9x fewer steps.

## Success and failure, stated in advance

* **ALL FOUR AGREE at matched voltage** -> the ramp rate is irrelevant over 20x,
  the quasi-static criterion is RETRACTED outright, and 1 us becomes the default
  (a 20x cost saving over the reference).
* **1 us DEPARTS but 2 us AGREES** -> the limit is bracketed between -100 and
  -200 V/us, and there IS a rate constraint, just far weaker than I claimed.
* **BOTH FASTER ARMS DEPART** -> the criterion has a grain of truth and 5 us is
  near the practical limit.
* **A FASTER ARM RUNS AWAY** (electrode sign flip, `|I_cond|` > `I_sc`) -> the
  rate matters a great deal and the reference ramp was load-bearing.

## COMPARE AT MATCHED VOLTAGE, NOT MATCHED TIME

    V_src = -(200/t_ramp[us]) * t[us]

so 20 us at t = 20, 5 us at t = 5, 2 us at t = 2 and 1 us at t = 1 are all the
SAME operating point (-200 V). A matched-time comparison measures the ramp, not
the discharge.

## Extraction

    ~/ct-env/bin/python ../../tools/compare_at_voltage.py <caseA> <caseB>
    ~/ct-env/bin/python ../../tools/news.py <case>            # n_e and Te evolution
    ~/ct-env/bin/python ../../tools/settle_check.py <case>    # rate + circuit, one line
