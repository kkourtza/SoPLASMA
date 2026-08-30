# positiveStreamer_LMEA_fast — a fast bed for the electron-energy coupling

## What this is for

The **only** LMEA case in this repository that runs in seconds. Everything else
that solves an electron-energy equation lives on the 1.15M-cell streamer at
**~70 s per time step**, which is why so many electron-energy questions have
been answered slowly, or on a synthetic bed, or not at all.

This case exists so that a question about the **n_e / nEps_e outer-loop
coupling** — relaxation, acceleration, time-step control, step rejection — can
be asked and answered inside a coffee break.

## What it is NOT

**It does not resolve the streamer.** The production bed is this same 130x130
block plus *five* `refineMesh` levels; here there are none, so the cell is
~96 um — the same size as the streamer head itself.

Field magnitudes, propagation speed and the head structure are therefore **not
physically converged** and must never be quoted as physics. Use
`positiveStreamer_fixedMesh` (LFA) or the production LMEA cases for anything
about the discharge.

What *is* faithful here is the **coupling**, because the Picard gain between
`n_e` and `nEps_e` is set by the time step and the local rates, not by the cell
size. Verified rather than assumed — measured over the full run to 5e-10:

| | |
|---|---|
| wall clock | **75 s**, 49 steps (~1.4 s/step) |
| production bed, same physics | ~70 s **per step** — about 45x slower |
| `omega` | median 0.698, **min 0.158**, below 0.5 on 12% of steps |
| `rho` | median 0.807, max 174 |
| `deltaT GOVERNED` | 0 steps |

`omega` reaching 0.158 means the period-2 mode this machinery exists to damp is
genuinely active — that is the property the bed had to reproduce, and it does.

Two honest limits. The time-step governor **never fires** here (0 governed
steps): `omega` bottoms at 0.158 against a `relaxMarginHold` of 0.15, just
above. So this bed exercises the *relaxation* but not the *governor*, and a
question about the latter needs a stiffer configuration.

And `rho` sits at a median of **0.807** while every step converges — far above
the 0.30 boundary measured on the synthetic bed in `testAitken` CASE 15. **The
synthetic threshold does not transfer.** A real residual carries many more modes
than a two-field model problem, so any production `rho` threshold has to be
calibrated here, not there.

## Running it

```
./Allrun-serial          # serial; the mesh is small enough that MPI costs more
                         # than it saves
```

`endTime 5e-10` reaches the regime where the coupling bites. Shortening it to
2e-11 gives a ~5 second smoke test, but **omega then never leaves its 0.8 seed**
— the coupling is not yet active that early, so a run that short proves the case
builds and nothing more.

## What to watch in the log

```
omega [coupling margin]:  min 0.263779104555  max 1
rho [contraction]:        max 0.196437843183   (diagnostic; omega still governs)
```

* `omega` — the joint Aitken factor. 1 means nothing needed damping; small
  means the coupling is only stable because it is being damped hard, which is
  the time-step controller's stability margin.
* `rho` — max `||r_k||/||r_{k-1}||` over the step. Diagnostic only. It is the
  cross-scheme signal, since Anderson acceleration has no `omega`.
  **Caveat measured in `testAitken` CASE 15:** computed over *all* cells it goes
  blind wherever densities sit on their floor, because a clamped cell's residual
  is frozen and those terms swamp the norm.

## Settings that matter here

| key | value | why |
|---|---|---|
| `energyModel` | `localEnergy` | LMEA. This is the point of the case. |
| `adaptiveRelaxation` | on (LMEA default) | the joint Aitken pair |
| `onNonConvergence` | `retryStep` (default) | never bank a non-converged step |
| `maxSpeciesCo` | 1.5 | as production |
| `numberOfSubdomains` | 4 | only used by a parallel variant |

## Related

* `testAitken` — the synthetic unit bed for the same coupling, with **ground
  truth**: the loop gain is set directly, so "did this converge" is known rather
  than inferred. Prefer it when the question can be posed there; it runs in a
  second. This case is for when the question needs real chemistry, a real
  Poisson solve and spatially varying coefficients.
* `docs/options-reference.md` §2.6 — the outer-loop controls.
