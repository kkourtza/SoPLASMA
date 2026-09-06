# LFA numerics tuning — why these cases exist and what they may NOT be used for

Created 2026-09-06. The LFA arm of V4 part B (`grubert2009_LFA`) will not
integrate past t ~ 1.1e-7 s: 1216 rejections, ALL "outer loop did not
converge", the contraction factor rho oscillating 0.63 <-> 1.18, omega pinned
at its 0.8 start, and dt stuck near 5e-14 -- which puts endTime at ~1e9 steps.
Meanwhile every physical limiter has headroom and the temporal error is 0.033
against a target of 1, so nothing physical is asking for a step that small.

There is a structural reason LFA is the weaker arm here, recorded in
plasmaTransport.C: under LFA only ONE field (`n_e`) enrols in the joint Picard
residual, and that single-field configuration is the one measured to be "the
WORST of all -- 191 correctors and omega on its floor, against 91 for the joint
pair". Under LMEA there are two (`n_e` and `nEps_e`). So Aitken has strictly
less to work with under LFA, by construction.

## THESE CASES ARE NOT A COMPARISON AND MUST NOT BE QUOTED AS ONE

They vary NUMERICS only, to find settings under which the LFA arm integrates.
Whatever wins is then applied to **BOTH** arms of V4 part B and both are rerun,
because the validity of that comparison rests on the two arms differing in
exactly one key -- `electronEnergyModel`. Tuning one arm alone would buy a
running LFA case at the cost of the only thing that makes the LFA/LMEA
difference attributable to the closure.

Baseline for all of them:
  /home/kkourtza/soplasma-scratch/validation/grubert2009_LFA
which is the arm as launched, and is the row these must be read against.

## The variants

| case | changed | why |
|---|---|---|
| `tuneLFA_w04` | `relaxOmegaStart 0.4` | damp harder from the first corrector. The source calls 0.8 "A STARTING POINT, NOT A UNIVERSAL CONSTANT, and it is CASE-SPECIFIC" |
| `tuneLFA_w02` | `relaxOmegaStart 0.2` | is the trend monotone, or is there an optimum? One point cannot tell |
| `tuneLFA_c150` | `maxCorrectors 150` | is it too few iterations, or the wrong omega? This separates those two causes |

## The discriminating observable, and the earliest time it is visible

Rejections accumulated by t = 1.5e-7 s, just past the choke at 1.1e-7. The
baseline takes ~1200 by then. A variant that works shows tens, and dt recovers
towards the 1e-9 scale the Courant and dielectric limits allow. No variant needs
to run to endTime to answer this, so none of them do -- endTime is cut to
2e-7 s.

## Extraction

    for d in tuneLFA_*; do
      printf "%-14s steps=%-7s sim=%-12s rej=%s\n" $d \
        "$(grep -cE '^Time = ' $d/logs/log.soPlasmaFoam)" \
        "$(grep -E '^Time = ' $d/logs/log.soPlasmaFoam | tail -1 | awk '{print $3}')" \
        "$(grep -c DISCARDING $d/logs/log.soPlasmaFoam)"
    done
