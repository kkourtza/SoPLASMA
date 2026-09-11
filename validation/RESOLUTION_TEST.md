# Sheath-resolution test — PREDICTION RECORDED BEFORE THE RESULT

Written 2026-09-06, both arms relaunched on the graded mesh, no result yet.

## What was measured on the uniform 50 um mesh

Both arms break down at V ~ 185 V, which is CORRECT: Paschen for argon at
pd = 100 Pa cm = 0.75 Torr cm sits near its minimum, ~140-200 V. After
breakdown the LMEA arm went to

    I_cond  = 663 mA/cm^2      (Grubert: 0.511; I_disp = 0.009, so this is
                                real conduction, not a displacement artefact)
    n_e max = 6.27e19 m^-3     (Grubert peak: 2.478e15 -- 25300x)
    lambda_D = 1.33 um         against a 50 um cell: UNDER-RESOLVED 37.7x

and dt collapsed with the current oscillating between 300 and 1170 mA/cm^2.

## The two competing explanations

1. **RESOLUTION.** The space-charge sheath is what chokes the current in a
   glow. A mesh that cannot form a 1.33 um sheath cannot limit the current, so
   the discharge keeps conducting; the density then rises and the sheath thins
   further -- a feedback through the mesh. Note the uniform mesh WOULD be
   adequate at the right answer (at n_e = 2.478e15, lambda_D = 211 um and
   dx/lambda_D = 0.24), so this only bites once the solution has overshot.

2. **SOLUTION BRANCH.** The runaway is real, and Grubert's 0.511 mA/cm^2 is a
   different branch that their STEADY solve finds (they solve the
   time-independent equations by FEM) and a ramped transient at fixed voltage
   does not reach.

## PREDICTION

The mesh is now graded to both electrodes: 400 cells, Bump 0.02, min 1.35 um,
max 66.8 um, ratio 49. Everything else is unchanged and identical in both arms.

* If (1) is right: breakdown still at ~185 V, but the current SATURATES far
  lower -- somewhere in the 0.5-10 mA/cm^2 range -- because the sheath can now
  form. dt should also stop collapsing, since the runaway is what drove it.
* If (2) is right: the current runs away again to hundreds of mA/cm^2 despite
  the resolution, and the discrepancy is not a mesh problem.

**If the current runs away again, explanation (1) is retracted, not patched
with more refinement.** One more refinement level would then be evidence about
convergence, not a fix.

## Extraction

    ~/ct-env/bin/python - <<'PY'
    import csv
    for d in ("grubert2009","grubert2009_LFA"):
        rows=[r for r in csv.reader(open(f"{d}/postProcessing/dischargeCurrent/current.csv"))
              if r and not r[0].startswith('#')]
        h=rows[0]; dd=[[float(x) for x in r] for r in rows[1:]]
        i=h.index("I_total"); A=2e-7
        pk=min(dd,key=lambda r:r[i])
        print(d, "peak", abs(pk[i])/A*1e3/1e4, "mA/cm2 at t =", pk[0])
    PY
