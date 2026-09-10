# Why this case cannot run in parallel yet -- two defects, measured 2026-09-07

Prompted by the user asking whether the refined uniform benchmark
(`validation/grubert2009_unif2000_300`, 10,000 cells) could be run in parallel to
speed it up. It cannot, and the reasons are in the framework, not the case.

Setup used: `decomposePar` with `method scotch`, `numberOfSubdomains 4`
(2520 cells per subdomain), Open MPI 4.1.6, `mpirun -np 4 soPlasmaFoam -parallel`.
`decomposePar` itself succeeded.

## DEFECT 1 -- the boundary-role library reads the PER-PROCESSOR path

`boundaryRoleLibrary::caseDeclaration()` builds

    const fileName p(runTime.path()/"configuration"/"boundaries");

`Time::path()` returns `<case>/processorN` in a parallel run, so every rank looks
for its own copy and all four die:

    Cannot read the case's boundary description:
        ".../grubert2009_unif2000_300_par/processor3/configuration/boundaries"

**FIX: use `runTime.globalPath()`**, which returns the case root in both serial
and parallel. One line, in
`src/models/electromagnetics/boundaryRoleLibrary/boundaryRoleLibrary.C` (~line
149). This is LAYER 1 of the G2 architecture, so the defect blocks parallel
execution for EVERY case that uses the boundary-role library -- i.e. all new
cases. Worth fixing on the next rebuild.

Workaround used to get past it for diagnosis only:
`for d in processor*; do ln -sfn ../configuration $d/configuration; done`
That is a HACK and must not be shipped -- it hides the defect.

## DEFECT 2 -- MPI message-size mismatch, and this one needs real work

Past defect 1, the run aborts with

    MPI_ERR_TRUNCATE: message truncated
    MPI_ERRORS_ARE_FATAL

A truncation error means a receive buffer was smaller than the matching send:
somewhere the ranks disagree about how much data to exchange. That is a genuine
parallel-correctness bug, not a configuration problem, and it is NOT diagnosed
here. Candidates worth checking first, in order of suspicion:

* the per-species / per-patch exchanges added with the wall-flux closures --
  anything sized from a LOCAL patch count rather than a reduced global one;
* `plasmaStepAudit` and `verifyChargeDensity`, both added 2026-09-06 and never
  exercised in parallel;
* the unit-potential / `C_self` surface integrals in `floatingElectrode` and the
  discharge-current diagnostic, which integrate over patches that may be split
  across ranks.

## Consequence, stated plainly

**Every result in this thread is from a SERIAL run, and parallel execution of
this solver is unverified.** Memory records a case where a fix applied to one of
two sibling models left electrostatics dead for three weeks
([[parallel-models-guard-divergence]]), so a parallel path that merely RUNS is
not evidence that it is right. When defect 2 is fixed, a parallel arm must be
verified against a serial reference on the same case before any parallel result
is used -- the natural check is the early trajectory, since a serial reference
already exists for this case.

## Practical effect on the benchmark

`grubert2009_unif2000_300` therefore runs SERIAL, at ~5.2 h to reach the current
reversal at t = 1.121 us. That is acceptable and it is what is running. The
parallel copy `grubert2009_unif2000_300_par` is left in place, not running, as
the reproducer for defect 2.
