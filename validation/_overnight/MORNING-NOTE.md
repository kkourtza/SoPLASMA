# Overnight state -- read this first (2026-09-08, ~01:30)

## What is running (24 ranks, 8 each)

All three carry **Grubert's own r = 0.36** and the **CompleteFlux** scheme, and
all branch from the same converged state, so voltage is the ONLY variable.
Measured breakdown voltage is 120.8 V.

    grubert2009_cfs_r036   -300 V   (2.48x breakdown)   was accelerating
    grubert2009_cfs_250V   -250 V   (2.07x breakdown)   ramp 300->250 by 5.0 us
    grubert2009_cfs_200V   -200 V   (1.66x breakdown)   ramp 300->200 by 4.5 us

`_overnight/STATUS.txt` has a line per arm every 10 minutes: t, n_e,max, Te,max,
Te,min, dt, non-converged count, ranks, aborts, clamp hits.

## What I killed, and why

  * `grubert2009_cfs_400V` -- ignited violently, dt fell to ~5e-13, throughput
    worked out at 23 h per 0.1 us and 167 days to its endTime. No path to an
    answer.
  * `grubert2009_cfs_long` (r=0, -300 V) -- same fate: ignited at 4.276 us,
    dt froze at 1.56e-12. Its remaining job was to be the r=0 control for
    `cfs_r036`, and the two had already tracked each other closely up to ~4 us,
    so little was lost. Its data is intact on disk.

Both were killed on your standing authorisation. Their case directories and
logs are untouched if you want to look.

## What to look for

  1. **Does any arm settle?** The signature is the growth rate turning over AND
     STAYING down. I twice called a turnover from a single interval and was
     wrong both times -- it dipped and re-accelerated. Want several consecutive
     intervals.
  2. **dt is the health indicator.** ~1e-10 healthy; ~1e-12 means it has frozen
     in an avalanche and will not finish.
  3. **A cathode fall.** NONE had formed as of the last check. Run
     `python3 ../tools/news.py <case>` for the evolution and
     `python3 ../tools/cathode_fall.py <case> <time> 100` after
     `reconstructPar -time <t>` + `postProcess -func writeCellCentres -time <t>`.
     The tool now WARNS if the zero crossing is really the anode-side reversal
     rather than a cathode fall -- it caught itself reporting one at 87% of the
     gap.
  4. **If an arm reaches a steady glow**, extract phi_c, p*d_c and j/p^2 and
     compare against Engel & Steenbeck (Grubert Fig. 5). Note his own LMEA sits
     2-20x BELOW that data, so a shortfall is not automatically our defect.

## The honest headline so far

CompleteFlux fixed a real numerical failure -- the case that aborted at 1.680 us
on standard/ROUNDF completed 2.0 us with Te,max 11.38 eV and zero clamp hits.
It has NOT yet produced the target physics: no cathode fall has formed in any
arm, so the discharge is still Townsend-like and the similarity comparison is
out of reach. Stability is not the same as reproducing the discharge.
