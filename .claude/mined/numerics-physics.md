# numerics-physics

## Summary
This dimension covers what a numerical-analysis agent and a physics-validator agent must know about SoPlasma's solver: the equation set actually assembled (singleRegionPoisson + per-species driftDiffusion nEqn() + LMEA nEps eEqn(), coupled three ways), the two outer-solver strategies (segregated Picard+Aitken vs PETSc SNES/JFNK), the timestep-limiter stack and which limiter names itself, and the exact PETSc flag/dict-key vocabulary with measured iteration counts. It also carries the project's physics vocabulary (canonical field names, SI units and sign conventions, all confirmed against source), the sanity checks that have actually caught errors here (where the extremum is, quasineutral ppm cancellation, floor-vs-result, trivial-solution checks), and the closure's structural absences (no e-e collisions in the fluid pipeline, no EEDF shape/trapping, single-group LMEA, two-term tables keyed on one variable). The single most important meta-fact: several design docs carry `Status: DESIGN ... Nothing implemented`, and citing one as a capability is a recorded failure mode here. The second: this project's own numbers have repeatedly been invalidated by uncontrolled second variables and by knobs that never moved, so every claim below is stated with the measurement that backs it and the date it was taken.

## Facts

### Design docs carry an explicit Status line and several are DESIGN-ONLY. `docs/design/adaptive-dt-design.md:3` = 'Status: DESIGN + step-1 measurement, 2026-09-06. Nothing implemented.'; `docs/design/steady-and-stability-design.md:3` = 'Status: DESIGN for review, 2026-09-06. Nothing implemented.'; `docs/design/schur-semiimplicit-poisson-preconditioner.md:3` = '**Status: DERIVED, not yet implemented.**'; `docs/design/numerics-generator-plan.md:3` = 'Status: DESIGN, 2026-09-06. Nothing implemented.'; `docs/design/stationary-solver-plan.md:3` = 'Status: PARTLY SUPERSEDED, 2026-09-11'. `newton-outer-solver-design.md`, `newton-ignition-experiments.md`, `flux-schemes-theory-and-implementation.md`, `steady-mode-spec.md` and `theory/scharfetter_gummel.md` have NO Status line and are mixed implemented/aspirational — read the per-section verdict, not the title.
**Evidence:** grep -n 'Status' over /home/kkourtza/soplasma-scratch/docs/design/*.md and docs/theory/*.md

**Rule:** Before citing any design doc as a capability, run `head -5 <doc>` and read the Status line; if it says DESIGN/DERIVED, say 'designed, not implemented'. For the un-statused docs, quote the section heading's own verdict (headings encode it: 'BUILT AND PASSING', 'REFUTED', 'RETRACTED', 'SUPERSEDES').

**Cost:** Misciting a DESIGN doc as a capability is a recorded failure mode in this project.


### The actual equation set, confirmed by reading the assembly code, not assumed: Poisson is `fvm::laplacian(epsilon_, ePotential_) == -chargeDensity_` (singleRegionPoisson.C:187-191), or under `poissonScheme semiImplicit` the same operator with effective permittivity `epsilon + dt*electricalConductivity` and RHS `-chargeDensity - dt*diffusiveChargeSource` (singleRegionPoisson.C:214-254). Species: `virtual tmp<fvScalarMatrix> nEqn() const` (plasmaTransportModel.H:146), `driftDiffusion::nEqn()` = `fvm::ddt(n) + fvm::div(phi,n) - fvm::laplacian(D,n)` with NO reaction terms (driftDiffusion.C:328-400); chemistry sources are added by `plasmaTransport::mechanismSourceTerms()` (plasmaTransport.C:3394+). Energy: `localEnergyEnergyModel::eEqn()` (localEnergyEnergyModel.C:1312+).
**Evidence:** docs/design/newton-outer-solver-design.md:557-600 ('Real system's actual equation structure, confirmed by reading the code'); src/models/plasmaModels/plasmaTransport/plasmaTransportModels/driftDiffusion/driftDiffusion.C:328-400

**Rule:** When reasoning about the residual or the Jacobian, treat `nEqn()` as transport-only and remember chemistry enters separately via chemP_/chemL_ — a 'species equation' block in the Pmat that carries fvm::Sp(chemL) is still missing -dP/dn.

**Cost:** The missing -dP/dn term (the ionisation frequency) cost the preconditioner its knowledge of the avalanche; adding it bought +41% advance/step.


### The outer sequence under Picard is, per timestep inside `while (pimple.loop())` (soPlasmaFoam.C:463-522): `em->solve()` -> `transport.solve(...)` (all species in a flat loop, plasmaTransport.C:1266-1286) -> `energy->solveSpeciesEnergy()` -> `species.updateChargeDensity()` -> `transport.updateSurfaceCharge()`. Aitken relaxation (`plasmaOuterRelaxation`) is invoked from INSIDE the species and energy solves (plasmaTransport.C:1311-1314, plasmaEnergy.C:156-167), auto-triggering a JOINT relaxation once every enrolled field has contributed.
**Evidence:** docs/design/newton-outer-solver-design.md:557-620

**Rule:** Never describe the species loop as nested per-species Picard; it is one flat loop. And never add manual `relaxationFactors` — Aitken is already running and fixed factors fight it.

**Cost:** Adding relaxationFactors 0.3 to a working pseudo-transient took it from 9472/9472 converged to 0/10 converged (steady-mode-spec.md 1.2).


### Picard's failure at a negative-differential-resistance point is NOT a crash: the contraction factor `rho` crosses 1.0 and climbs while the clock keeps advancing with the coupling unconverged. Measured on grubert2009_pseudo t=0->4e-8 (40,000 steps, dt 1e-12 fixed, 953 s): rho crossed 1.0 at step ~18,878 and reached ~1321, omega stayed 0.8-1.0, PIMPLE kept its fixed 4 correctors, no FatalError. At ignition on grubert2009_ballast400 Picard instead DIES: dt collapses 4.6e-13 -> 7.2e-15 -> ~1e-15 with 31-42 correctors/step.
**Evidence:** docs/design/newton-outer-solver-design.md:1210-1224; docs/design/newton-ignition-experiments.md:60-84

**Rule:** Two distinct Picard failure modes exist — silent non-contraction (rho>1, clock advances) and dt collapse with corrector saturation. Diagnose by reading rho AND achieved dt AND correctors/step, never by 'did it crash'.

**Cost:** A run that 'completed' with rho ~ 1321 produced 40,000 steps of unconverged coupling that read as success.


### `rho [contraction]` IS A KNOWN-UNRELIABLE SIGNAL and governs nothing. Correlated over 12 cases and ~1.35 MILLION steps: grubert2009 reports median rho = 378.2 with ZERO discarded steps; grubert2009_iset reports rho = 22.05 with zero discards, on a run verified against the analytic I_set/C_gap slope to 1.33% and igniting within 4% of the predicted voltage.
**Evidence:** docs/design/adaptive-dt-design.md:46-83 (STEP 1 RESULT)

**Rule:** Do not call a run divergent on `rho` alone. Corroborate with achieved dt, correctors/step, discard count, and a field-level check. Use `omega [coupling margin]` (the actuated Aitken damping) in preference to rho.

**Cost:** rho>1 fires on demonstrably accurate runs; a governor keyed on it would throttle healthy cases.


### Three dt controllers run simultaneously and disagree by five orders of magnitude. Measured on grubert2009_iset, ~16,000 steps: the temporal-error PI wants 9.85e-09 (752x larger than actual) and binds on 18 steps; the ENERGY Courant `maxEnergyConvectiveCo`, HAND-SET to 15, wants 1.31e-11 and binds on 15,723 steps; the robustness signal wants 1.92e-12 and never binds. The number that decides is the one a user typed — the stated G1 violation.
**Evidence:** docs/design/adaptive-dt-design.md:5-24

**Rule:** When asked why dt is what it is, read `deltaT set by:` from the log (plasmaTimeControl.C:1996) rather than assuming the accuracy controller is in charge. The source states the health test: 'if one of them still names itself in `deltaT set by`, the controller is not in charge and that is the thing to investigate' (plasmaTimeControl.C:1289-1292).

**Cost:** Tuning the accuracy target while a hand-set Courant cap binds 97% of steps changes nothing.


### The complete set of dt limiter names, each bound through `bind()` in plasmaTimeControl.C so it self-identifies in `deltaT set by`: 'maxDeltaT' (1188), 'contraction ceiling' (1282), 'temporal error (PI)' (1322), 'dielectric relaxation' (1345), 'Co_conv' (1424), 'energy relaxation' (1509), 'Co_diff' (1530), 'Co_chem' (1591), 'voltage rise' (1619), 'coupling margin (backoff)' (1644), 'coupling margin (hold)' (1649), 'outer loop not converged' (1664), plus non-bind assignments 'growth cap (errorMaxGrow)'/'growth cap (1.2x)' (1691-1693), 'rejection memory' (1706) and 'minDeltaT floor' (1749).
**Evidence:** grep -n 'bind(' and 'dtLimiterName_ =' /home/kkourtza/soplasma-scratch/src/tools/controls/timeControl/plasmaTimeControl.C

**Rule:** Treat this list as the closed set of legal `deltaT set by` values. A dt you cannot attribute to one of these is a bug in the reporting, not physics — three limiters once set dt with a bare min() instead of bind(), so the log named the PREVIOUS limiter.

**Cost:** The log said `deltaT set by: Co_conv` while printing Co_conv = 0.224 against a cap of 1.5 — a confident wrong answer that is worse than no diagnostic.


### The dielectric relaxation limiter is `dt <= maxDielectricRelaxationRatio * eps0 / max|sigma|`, i.e. tau = epsilon/sigma, and it is OFF BY DEFAULT: `limitDielectricRelaxationRatio` defaults false, `printDielectricRelaxationRatio` defaults false, `maxDielectricRelaxationRatio` defaults 1.0.
**Evidence:** src/tools/controls/timeControl/plasmaTimeControl.C:36-38, 279-287, 1333-1346, 2105-2114

**Rule:** To measure dt/tau_dielectric on any run, set `printDielectricRelaxationRatio true` in the `plasmaTimeControl` sub-dict of `system/plasmaSimulationControls` — printing does not limit. Do not assume the ratio is bounded by default.

**Cost:** The AP claim is stated in units of this ratio; without printing it you cannot say what dt/tau a run achieved.


### Co_chem is NOT a stability limit. It bounds how far chemistry moves the STATE in one step — a SPLITTING-error constraint, because the field and the rate coefficients read from it are frozen while the chemistry integrates over dt (`implicitRate` is unconditionally stable and the ODE path substeps under its own error control). Its source is `transport.maxChemStateRate()`, the per-species fractional net rate of change, falling back to the legacy air-fitted `k_eff` only when no mechanism assembles P/L. Growth only, never mag(): a negative net rate is attachment-dominated decay and is benign.
**Evidence:** src/tools/controls/timeControl/plasmaTimeControl.C:1560-1594

**Rule:** Never justify Co_chem as 'stability'. If Co_chem binds, the fix is tighter coupling or a smaller splitting window, not a different ODE integrator.

**Cost:** Limiting on L rather than on state movement pinned Co_chem at its cap every step (measured 2026-08-20, plasmaTransport.H:1122).


### MEASURED TIMESCALE SEPARATIONS. Grubert 100 Pa argon, 1 cm gap: electron transit tau_e = 3.48e-08 s vs ion transit tau_i = 1.07e-05 s, ratio 307 (at 100 V). Avalanche e-folding 149 ns shortening to 33 ns while ion transit across the gap is 4-7 us — 30-45x slower, which is WHY no cathode fall has ever formed in this project. Ionisation e-folding gamma^-1 measured 0.400 ns in the runaway; 1/nu_i = 14.19 ns at 431 Td and 0.73 ns at 2450 Td from the argon ionisation table.
**Evidence:** memory density-floor-is-a-source; memory avalanche-outruns-ion-transit; docs/design/adaptive-dt-design.md:199-224; memory reducede-table-axis-is-si (table cross-check values)

**Rule:** Before claiming a discharge structure has formed, check that the run covered at least one ion transit time (4-7 us for this gap), not just one avalanche time. State both timescales in any stiffness argument.

**Cost:** Every Grubert run to date reached at most a few us; a cathode fall is a space-charge structure and cannot form on a 30-45x separation.


### The density floor `minNumberDensity` is a per-step CLAMP, not an initial condition, and is therefore a CONTINUOUS VOLUMETRIC SOURCE amplified by the ion/electron timescale ratio: ions per clamped electron = N*k_iz*tau_i, measured 691 at the Grubert case's 100 V, so a 1e11 floor bootstraps n_Arp ~ 6.9e13 (measured 3.6e13-1.6e14). At a 1e11 floor the floor's own ions gave E_sc/E_applied = 1.251 and drove a sub-breakdown gap into local breakdown at -104 V against a computed 120.8 V. `0/n_e` is written `internalField uniform 0` and the clamp supplies the seed, so the floor and the initial condition are conflated in the case files.
**Evidence:** memory density-floor-is-a-source; memory grubert-operating-current-is-1e-6-A

**Rule:** Choose the floor by CONVERGENCE (1e9 and 1e7 agreed to 2.4% while both differed from 1e11, so 1e9 was picked). Before reading any density result, check `n_e min` — if it IS the floor, max/min is a ratio to the floor, not a physical dynamic range, and max/min = 1.22 means an undeveloped state, not a result.

**Cost:** A whole-gap distributed weak cathode, ionising volumetrically instead of closing the cathode-fall/gamma loop; profiles sat at exactly 1.000e11 through the cathode region in every arm.


### THE AP CLAIM, MEASURED. The semi-implicit Poisson is INCONSISTENT inside a Newton residual: its derivation PREDICTS rho^{n+1} ~ rho^n - dt*div(sigma*E) and substitutes it into Gauss's law, but a Newton residual already carries rho^{n+1} as an unknown, so the charge relaxation is counted twice and F(u)=0 has no root. The spurious term's size is exactly dt*sigma/eps, the dielectric relaxation ratio. Across 3167 steps in two runs, NO step whose dielectric relaxation ratio exceeded 0.3278 ever converged (in one run 0 of 161 successes at or above the lowest failing ratio). Replacing it with the plain (explicit) Poisson residual: SNES failures 4.92% -> 0% (0 of 19), retries 10 -> 0, dt ~3.4e-12 -> 7.40e-11 (22x), max dielectric ratio 0.60 -> 19.52, later confirmed to 308 over 152 steps with 0 failures.
**Evidence:** docs/design/newton-ignition-experiments.md:692-764 (section 23); memory ap-proof-newton-removes-dielectric-constraint

**Rule:** Under `outerSolver newton`, set the Poisson scheme to `explicit`, never `semiImplicit`. Confirm from the run's OWN log line ('Model: singleRegionPoisson Poisson scheme: ...'), because these cases set it through a `$poissonScheme` variable.

**Cost:** 110-400x the timestep on Grubert: explicit arms reached further in 214-261 steps than semi-implicit arms did in 5,000-11,400 steps.


### The accuracy of the 22x larger explicit-Poisson step is VERIFIED, not assumed: explicit Poisson with dt free vs explicit Poisson capped at 2e-12 (417 steps, 0 failures) at a common time gave n_e max 1.572e17 vs 1.581e17 (0.57%), n_e min 6.661e15 vs 6.749e15 (1.3%), Emag max 9.657e04 vs 9.707e04 V/m (0.52%). The large difference from Picard is NOT a large-dt error — it is that the semi-implicit form modifies Poisson by O(dt*sigma/eps), which was 20-60% in the Picard runs.
**Evidence:** docs/design/newton-ignition-experiments.md, section 23 ('CONFIRMED over a longer run, and THE ACCURACY IS VERIFIED')

**Rule:** When a solver change buys dt, the required control is the SAME solver at capped dt, compared at a common TIME — not against the previous solver, which differs in the equation as well as the step.

**Cost:** Without it the 22x reads as speed bought out of the solution.


### THE JFNK COST MODEL, and the 1500x defect it hid. One Krylov iteration = one full nonlinear residual assembly. Profiled: a residual evaluation costs 3.07 ms (unpack 0.17 s, derived 2.27, transport 4.62, CHEM 41.11, energy 23.59, residual 20.34 over 30,020 calls). The 92 s per timestep came from 30,020 calls where ~20 were justified; 30020/3 ~ 10,000 = exactly the number of unknowns (5 fields x 2000 cells), the signature of finite-differencing a Jacobian COLUMN BY COLUMN. Cause: `SNESSetUseMatrixFree(snes, PETSC_TRUE, ...)` with `mf_operator = TRUE`, which forces PETSc to build a real Pmat via SNESComputeJacobianDefault. Fixed by `MatCreateSNESMF` + `SNESSetJacobian(snes, Jmf, Jmf, MatMFFDComputeJacobian, nullptr)`. Result: 30,020 -> 16 residual calls/step, 97.88 s -> 0.195 s/step.
**Evidence:** docs/design/newton-outer-solver-design.md:950-1055

**Rule:** Before optimising anything inside the residual, count residual evaluations per timestep and divide by the SNES iteration count. If the quotient is near the DOF count, PETSc is assembling a Jacobian column-by-column — that is the bug, not the residual cost.

**Cost:** Lagging chemistry (45% of a call) would have bought under 2x while the real defect was a factor of 1500.


### The 449k-cell cost arithmetic that explains the 1.15M wall: 50 SNES x 200 KSP = 10,000 evaluations x 4.16 s = 11.6 h/timestep against 10.06 h observed. Measured at 449k, dt 1e-11, streamer bed: Picard 14.6 s/step vs Newton 142 s/step (~10x). Warm-started Picard survives dt = 5e-11, so Newton must hold ~5e-10 merely to BREAK EVEN; nothing measured to date shows Newton taking a larger step at all. Newton converges at 40k, 211k and 449k but not at 1.15M — proven to be MESH SIZE by running 1.15M both warm (10 h CPU/rank, still step 1) and cold (3 h, stuck after handover), both at 99.9% CPU.
**Evidence:** memory newton-vs-picard-benchmark-state

**Rule:** Measure Newton-vs-Picard as SPEED PER NS OF SIMULATED TIME AT EQUAL ACCURACY, never per step — Newton costs ~4-5x per step, so a per-step comparison is rigged for Picard. State the break-even dt explicitly.

**Cost:** A cold-start 449k ladder was briefly read as 'JFNK does not pay'; it answers nothing because it never reaches the regime (peak n_e still at the 1.3e13 seed, dt/tau ~ 2e-7, Picard using 2-4 correctors against a cap of 20).


### A HARDCODED `-ksp_max_it 100` produced EVERY negative JFNK verdict. Measured on grubert2009_pseudo (2000 cells, dt=1e-10, 400 steps, same 398 SNES solves): kspMaxIt 100 -> 10.6% failure rate, 38 DIVERGED_LINEAR_SOLVE, KSP median/p90/max 4/14/89, 121,161 evals, 326.1 s. kspMaxIt 1000 -> 1.0% rate, 0 linear-solve failures, 4/14/115, 123,234 evals, 334.2 s. The 38 killed solves needed 104 to 115 iterations — 4 to 15 short, 0.92% of 4,235 solves. PETSc reports it as DIVERGED_LINEAR_SOLVE, which reads as failure but is a budget expiring.
**Evidence:** memory ksp-cap-was-the-bug; commits f17a489, aed5dbe; snesNewtonSolver.C:1530

**Rule:** Raise `kspMaxIt` to 1000 before diagnosing any DIVERGED_LINEAR_SOLVE. Never read a tail statistic off a CAPPED distribution — the 'bimodal preconditioner cliff' was the cap truncating the distribution; the true one is median 4, p90 14, max 115.

**Cost:** Invalidated: the 449k streamer ladder's negative Newton verdict, the pseudo-transient dt ceiling of 1e-11 (at dt=1e-10 the rate is now 1.0% against 21%, so ~10x on dt and a ~12 h run becomes ~1 h), and the needleDBD multi-region result.


### The physics-based PCSHELL is the STRONGER preconditioner, contrary to an earlier reading. Same bed: PCSHELL (`assembledPmat false`) gives KSP median/p90/max 3/6/22 against PCFIELDSPLIT's 4/14/115 (5x tighter tail) and 113,208 vs 121,161 residual evaluations (8% fewer), with 0 linear-solve failures. Its 4.3% headline rate is ALL line-search failures (17 of them, reason -6) — a different defect. Both `kspMaxIt` and `assembledPmat` are now reachable from the `newtonSolver` dict, defaulting to prior behaviour.
**Evidence:** memory ksp-cap-was-the-bug (table); commits f17a489, aed5dbe; docs/CAPABILITIES.md

**Rule:** Compare preconditioners by the REASON BREAKDOWN (-3 linear solve / -5 maxIt / -6 line search) and the KSP iteration distribution, never by a single failure-rate number.

**Cost:** The headline rate hid that the shell was the better preconditioner, and the wrong one was left as default.


### THE SCHUR COMPLEMENT WAS INERT for most of the effort. PCFIELDSPLIT on splits {phi, transport} forms S = A_tt - A_tphi*A_phiphi^-1*A_phit. A_phit (charge density, species->Poisson) was assembled; A_tphi (d(drift)/d(phi), Poisson->species) was ABSENT, i.e. zero, so S = A_tt exactly and the whole Pmat was doing the work of a block-triangular preconditioner. The derivative: phiE = -snGrad(ePotential)*magSf (singleRegionPoisson.C:39) and the drift term is div(Z*mu_f*phiE*n), so d/dphi[div(Z*mu*phiE*n)] = -laplacian(Z*mu*n, dphi). The energy row gets the same with Ze*mu_eps*nEps. The (phi, energy) block is CORRECTLY absent: nEps carries no charge, so d(F_phi)/d(nEps) is identically zero — that coupling is genuinely one-way.
**Evidence:** docs/design/newton-ignition-experiments.md:364-421 (section 10)

**Rule:** Before tuning any sub-preconditioner, enumerate which Jacobian blocks the Pmat actually contains and check each against the residual it is supposed to differentiate. A missing block cannot be preconditioned around.

**Cost:** hypre, bjacobi, selfp, 5x the Krylov budget, Eisenstat-Walker and three line searches all failed — every one tested against a matrix missing the physics.


### THE COUPLING SIGNS, derived from this project's own residual and NOT to be copied from a survey. residualCallback computes F_0 = lap(effEps, phi) - rhsSource with rhsSource = -chargeDensity - dt*(...), i.e. F_0 = lap + chargeDensity + dt*(...), and chargeDensity = sum_i Z_i*e*n_i. So dF_0/dn_i = +Z_i*e, a DIAGONAL block. (An external survey quoted -Z_i*e; that is a different sign convention and would be silently wrong here.) The ionisation derivative: F = ddt + div - lap - P + L*n, so dF/dn = L - dP/dn; for electron-impact ionisation P_e = k_iz*n_e*n_Ar so dP_e/dn_e = k_iz*n_Ar = THE IONISATION FREQUENCY, approximated as P_s/n_s (exact whenever production is first order in the species itself) and deliberately NEGATIVE on the diagonal.
**Evidence:** docs/design/newton-outer-solver-design.md:1592-1622 (design item 4); docs/design/newton-ignition-experiments.md:422-446 (section 12)

**Rule:** Derive every Jacobian sign from the residual as this codebase writes it. Also: every Pmat entry in block (row r, col c) must be multiplied by sX[c]/sF[r], because everything PETSc sees is scaled — getting that wrong yields a valid matrix of the WRONG system.

**Cost:** A wrong sign or missing scaling 'fails in the least obvious way possible' and is invisible except through preconditioner behaviour.


### EISENSTAT-WALKER (`adaptiveForcing`, `-snes_ksp_ew`) is the largest single option win, and it FAILED when first tested. Measured at the pre-ignition state, 90 steps: baseline mean dt 1.773e-10 (2.2% KSP failure) vs EW 3.190e-10 (+80%, 0.2% failure — a 13x drop). It failed in the first sweep because it ran against the BROKEN Pmat: an adaptive forcing term chooses the linear tolerance from how well the previous Newton step reduced ||F||, so it needs a Jacobian whose steps mean something. CAVEAT, measured at the deep-avalanche state (t=1.90026e-6, N=30): EW gives only +4.7% on dt while causing 5x more DIVERGED_MAX_IT and 50% more retries.
**Evidence:** docs/design/newton-ignition-experiments.md:480-509 (15), 604-635 (20)

**Rule:** RE-TEST OPTIONS AFTER FIXING THE MODEL THEY WERE TESTED AGAINST. And measure at BOTH an easy and a hard state before making anything a default — an option tuned on the easy phase is not automatically right for the hard one.

**Cost:** EW, hypre, maxIt and the line searches were all dismissed on evidence gathered with a broken Pmat, and one of them was the biggest win available.


### With EW on, `chemJacobian` adds NOTHING, measured at equal sample size: EW + chemJacobian on gives mean dt 3.653e-10 @16 / 3.325e-10 @30; EW + chemJacobian OFF gives 3.653e-10 / 3.316e-10 — identical. The +41% attributed to the ionisation derivative was measured WITHOUT EW; EW subsumes it. `chemJacobian` stays default TRUE for correctness reasons (a genuinely missing Jacobian term, nearly free, worth +41% whenever EW is off), not performance.
**Evidence:** docs/design/newton-ignition-experiments.md:559-586 (section 18)

**Rule:** ALWAYS COMPARE AT EQUAL N. `ew_nochem` appeared to beat `ew_v2` (3.416e-10 vs 2.955e-10) purely because its average covered 34 steps against 93. This unequal-sample trap occurred four times in one session.

**Cost:** A commit message overstated the ionisation derivative's standalone importance in the final configuration.


### The matrix-free differencing parameter must be built from the accuracy F is actually evaluable to, not machine epsilon. F integrates a per-cell ADAPTIVE STIFF CHEMISTRY ODE whose internal step sequence changes discontinuously with its input. Measured sweep of `-mat_mffd_err` at t=2e-9: ~1.5e-8 (PETSc default, PETSC_SQRT_MACHINE_EPSILON) DIVERGED_BREAKDOWN with true residual 2.4e9; 1e-6 DIVERGED_BREAKDOWN; 1e-4 CONVERGED in 3 iterations. After the meanE-ordering fix the valid measurement (20 steps) is: 1e-8/1e-7/1e-6 fail at step 1, 1e-5 gives 20/20 converged. `newtonSolver/mffdErr` therefore ships as a documented WORKAROUND, default 1e-5. F IS STILL NOT SMOOTH AT THE 1e-6 LEVEL.
**Evidence:** docs/design/newton-outer-solver-design.md:1007-1043, 1124-1157 (including the RETRACTED first sweep); snesNewtonSolver.C:1136

**Rule:** A case that diverges under Newton should raise `mffdErr` first. If `mffdErr` can eventually drop toward PETSc's default, that is the objective measure that F has become a clean function of u.

**Cost:** A whole sweep was invalidated because `mffdErr_` still defaulted to 1e-4 and setPetscOptions() inserted `-mat_mffd_err 1e-4` AFTER the swept PETSC_OPTIONS value — all three runs were the same configuration, which is why they agreed to the last digit.


### FGMRES, not GMRES, at both levels, and for a stated reason: PCFIELDSPLIT+Schur solves its blocks iteratively, so the preconditioner VARIES between Krylov iterations and plain GMRES loses the Arnoldi relation (Knoll & Keyes 3.5). `-ksp_type fgmres -ksp_gmres_restart 100` removed a DIVERGED_BREAKDOWN at iteration 30 (GMRES's default restart length) and turned the failure mode into an honest DIVERGED_ITS. The identical defect recurred one level down: the phi inner solve took 55,65,71,71,82,85 iterations on successive applications, so S was not a fixed linear operator, GMRES's Arnoldi recurrence divided by a vanishing norm and the next MatMult raised SIGFPE.
**Evidence:** docs/design/newton-ignition-experiments.md:122-130 (section 3), 861-943 (section 25c)

**Rule:** Any inner solve inside a Schur application must be a FIXED linear operator. Either make it exact (`-fieldsplit_phi_ksp_type preonly -fieldsplit_phi_pc_type lu`) or fix its iteration count and skip its convergence test (`-fieldsplit_transport_ksp_type richardson -ksp_max_it 5 -ksp_convergence_test skip`). Keep FGMRES outside as the safety net.

**Cost:** A SIGFPE inside MatMult_SeqAIJ, ten log lines after the only inner-solve breakdown in ~40,000 inner solves.


### THE `newtonSolver` DICT KEYS AND THEIR DEFAULTS, read from the constructor: rtol 1e-8, maxIt 50, kspMaxIt 100, assembledPmat true, mffdErr 1e-5, bounded false, petscOptions (string::null), rebalanceScales true, chemJacobian true, chemJacobianElectronOnly true, chemCrossJacobian false, jouleJacobian false, schurOnPhi false, perCellScaling false, perCellScaleFloor 1e-6, jouleJacobianRatioMax 10.0, adaptiveForcing true, scaleDiag false, extrapolateGuess false, extrapolatePotential true. Selection is `outerCoupling.outerSolver picard | newton` plus `outerCoupling.newtonSolver { type SNES; }` in `system/plasmaSimulationControls`.
**Evidence:** src/numerics/newtonSolverPETSc/snesNewtonSolver.C:1132-1157; docs/design/newton-outer-solver-design.md:632-660

**Rule:** Set numerics through this dict, not through `-snes_*` on the command line (PetscInitialize does not strip its own flags from argc/argv, so OpenFOAM's argList rejects them). `petscOptions` is the case-level escape hatch for raw PETSc flags.

**Cost:** Two 'Newton' arms once ran PICARD for 24,000 steps because the case had no `outerSolver` key for the sed to match.


### `chemCrossJacobian` HURTS and is the author's own refuted idea: approximating ALL of a species' production as first order in n_e gave mean dt 1.358e-10 against 1.904e-10 with it off — MINUS 29%. The approximation is right for electron-impact channels and wrong for ion-neutral and recombination ones. `bounded` (SNESVINEWTONRSLS) also FAILS: 360 steps, 0 fatals, 2383 SNES calls looks like success but advanced only 4e-11 in total with dt COLLAPSED to 5.0e-15 and 135 retries — 10^5x worse than unbounded Newton on the identical state. Raising `maxIt` 50->200 HURT at the easy state (-16%).
**Evidence:** docs/design/newton-ignition-experiments.md:510-521 (16), 339-363 (9), 587-603 (19)

**Rule:** READ THE ACHIEVED dt, NEVER THE STEP COUNT. A high step count with a collapsed dt is grinding, not progress. Default chemCrossJacobian false and bounded false; do not re-derive either.

**Cost:** Bounded Newton's 360 steps read as a success until the achieved time was checked.


### THE CLAMP-OUTSIDE-THE-EQUATIONS DEFECT. `clampNumberDensities()` raises densities after the solve, so the discrete system being solved is not F(u)=0 and no Newton method can converge to it. Measured at a cold start: every term exactly 0 at step 1 -> ||F|| = 0 -> SNES reports CONVERGED_FNORM_ABS with 0 iterations (fake convergence); at step 2, ||F|| equals ||ddt|| to 0.03% in every transported block with div/lap/chem 100-1000x smaller — the residual was ENTIRELY the time-derivative of the clamp's own 0 -> 1e11 jump. A second instance one level down: `updateDerived()` writes a CLAMP back into nEps_, so the state F is evaluated at is not the state PETSc perturbed.
**Evidence:** docs/design/newton-outer-solver-design.md:788-827 (DEFECT A), 1145-1157

**Rule:** Any state change made OUTSIDE the equations makes F=0 unreachable. Check for clamps, one-shot seeds and post-hoc corrections inside or around the residual before diagnosing a Newton convergence problem.

**Cost:** 18,400 residual evaluations bought one Newton iteration and a factor-of-2 reduction.


### THE 'IMPOSSIBLE MEASUREMENT' AND ITS REAL CAUSE. Under PCNONE the preconditioned and true residual norms MUST be identical and were not (iter 0: 178.16/178.16; iter 1: 35.96/348.57). Three successive attributions were made and two were RETRACTED: 'genuine nonlinear stiffness' (wrong), 'conditioning, dt is the lever' (wrong, and its dt sweep was non-monotonic at a second state). The actual cause: `localEnergyEnergyModel::correct()` carries a ONE-SHOT LFA seed `if (seedFromLFA_) { seedFromLFA_ = false; nEps_ == meanE0*n_e; }` which fired on residual call #1 and no other, so F was a DIFFERENT FUNCTION at the base point than at every perturbed point. Measured: nEps block scaled residual 172.46 on call #1 against 0.00086 on call #2, a factor of 2e5 and ~97% of the reported initial norm.
**Evidence:** docs/design/newton-outer-solver-design.md:1405-1448

**Rule:** When the preconditioned and true residuals disagree under PCNONE, the operator is not a fixed linear map — look for state mutation inside the residual (one-shot initialisations, lazy construction, history buffers), not for stiffness or conditioning.

**Cost:** Every diagnosis built on the residual was invalidated; and the same seed is a PRODUCTION bug — `freshStart = timeIndex == startTimeIndex` is evaluated in the CONSTRUCTOR so it is ALWAYS true, meaning every LMEA restart with `initialMeanEnergy` silently discards its stored energy state, on the PICARD path too.


### PER-BLOCK SCALING IS MANDATORY AND IS STILL IMPERFECT. Unscaled, the combined SNES norm was 100% species: Poisson block ||F|| ~ 1.2e-3 against species ~2e24, combined 3.163e+24 — with rtol 1e-8 the target sits nineteen orders ABOVE the entire Poisson residual, so SNES could not see Poisson at all. Fix: sX[block] = rms(that block's field), sF[transported] = sX/dt, sF[Poisson] = rms(lap(effEps, ePotential)), computed ONCE per solveOuterStep and applied at four sites (packing x, unpacking x, writing F, and in the PC). Even after that the blocks spanned 0.0077 to 111 (a 14,000x spread), fixed by `rebalanceScales` (default true) normalising sF by the MEASURED initial residual — which made every block 44.72 and did NOT improve convergence (KSP 24/49/100 vs 25/49/100).
**Evidence:** docs/design/newton-outer-solver-design.md:828-861 (DEFECT B); docs/design/newton-ignition-experiments.md:284-338 (7, 8)

**Rule:** State per-block |F_scaled| in any Newton diagnosis. A block contributing ~1e-4 of the norm is invisible to the Krylov method however wrong it is. But do not expect rebalancing alone to fix convergence — it removes a confound, it is not a cure.

**Cost:** Block imbalance explains why every PETSc knob failed, and yet fixing it changed nothing — the two facts must be reported together or the next reader repeats the sweep.


### VOCABULARY TABLE, every entry confirmed in source. ePotential [V] (electromagneticsModel.H:92); E [V/m] volVector, `-fvc::grad(ePotential)` or `fvc::reconstruct(phiE)` per `poisson.EScheme grad|reconstruct`, default reconstruct; Emag [V/m] = mag(E); phiE [V*m^2] surfaceScalar = `-fvc::snGrad(ePotential)*mesh.magSf()` (singleRegionPoisson.C:39); reducedE [V*m^2] = Emag/N_background, dimensionSet(1,4,-3,0,0,-1,0); chargeDensity [C/m^3] dimensionSet(0,-3,1,0,0,1,0) = sum_i Z_i*e*n_i; surfCharge [C/m^2] dimensionSet(0,-2,1,0,0,1,0), owned by the gas; epsilon [F/m] dimensionSet(-1,-3,4,0,0,2,0); n_<species> [1/m^3]; nEps_<e> [eV/m^3] carried as dimless/dimVolume, the TRANSPORTED energy variable = n_e*eps_bar; meanE [eV] carried as DIMLESS, bare (not per-species) name, derived as nEps/n_e; T_<e> [K] = (2/3)*eps_bar/k_B; mu_<species> [m^2/V/s] dimensionSet(-1,0,2,0,0,1,0); D_<species> [m^2/s]; S_iz [1/m^3/s] dimensionSet(0,-3,-1,...); alpha [1/m]; alphaDx [dimless]; Sph [photoionisation source]; dSdEps_lmea [1/s].
**Evidence:** src/models/electromagnetics/electromagneticsModel.H:80-130 and .C:100-152; localEnergyEnergyModel.C:46-120; plasmaTransport.C:125-170; plasmaMobilityModel.C:50; plasmaDiffusivityModel.C:50

**Rule:** Use these exact names when grepping fields or writing diagnostics; species fields are `n_`+name, energy is `nEps_`+name, mobility `mu_`+name, diffusivity `D_`+name, temperature `T_`+name, while `meanE` and `reducedE` are GLOBAL bare names because TabulatedProperty1D resolves `lookupVariable` from the registry BY NAME and builds table paths as quantity+'_vs_'+lookupVariable.

**Cost:** `meanE_e` would look for muN_vs_meanE_e and fail; a wrong prefix silently finds no field.


### reducedE's TABLE AXIS IS SI, V*m^2, NOT TOWNSEND — and the failure is silent and plausible-looking. Table headers say so: `k[EI_AR_ION_AR] [m^3/s] vs reducedE [V m^2]`. `ENmax 55000` Td becomes an axis running to 5.5e-17. Feed a value in Td (12, 431, 2450) and every one is >= 5.5e-17, so the interpolator clamps and returns the TABLE MAXIMUM — for argon ionisation 2.728e-13 m^3/s, a perfectly reasonable-looking rate. Conversion: E/N [V m^2] = E/N [Td] * 1e-21. The source itself warns: 'Suspect a units mismatch: reducedE is in V m^2, while E/N' (plasmaReactionRates.C:514).
**Evidence:** memory reducede-table-axis-is-si; src/models/plasmaModels/plasmaReactionRates/plasmaReactionRates/plasmaReactionRates.C:514

**Rule:** THE TELL is identical output for 12 Td and 2450 Td. If a rate does not change across two decades of field, the axis units are wrong, not the physics. Cross-check any new table reader against the two recorded values: 1/nu_i = 14.19 ns at 431 Td and 0.73 ns at 2450 Td (validation/grubert2009_spike/COMPARE.md).

**Cost:** Every rate coefficient in the run returns the table maximum while looking sensible.


### THREE CURRENTS, THREE MEANINGS, and they are not checks on each other. I_total = I_cond + I_disp is Sato's EXTERNAL-CIRCUIT current, a volume-weighted integral over the whole domain (postProcessing/dischargeCurrent/current.csv) — what a series ballast drops voltage across. I_cond is the CONDUCTION current alone, and an external circuit needs THIS one, not I_total: using I_total in V = V_source - R*I makes that relation the RC differential equation in disguise, and evaluating it explicitly amplifies by R*C/dt per step — measured 885 at dt=1e-10 on the Grubert case, which diverged to 1e17 V in five steps. I_collected (floating.csv) is the NET CHARGE FLUX onto ONE conductor, a conduction current only; in a DBD I_total is dominated by displacement while I_collected is zero.
**Evidence:** src/tools/diagnostics/plasmaDischargeCurrent/plasmaDischargeCurrent.H:226-252; src/models/electromagnetics/floatingElectrode/floatingElectrode.H:115-132

**Rule:** Name which current you mean in every statement. Never compare I_collected against I_total as a consistency check — they are different quantities with different units of meaning.

**Cost:** Divergence to 1e17 V in five steps.


### SIGN CONVENTION FOR DRIFT: `convectivePhi() = Z * fvc::interpolate(mu) * em().phiE()` with `phiE = -snGrad(ePotential)*magSf` and Z the species charge NUMBER. So for electrons (Z = -1) the carrier flux is +mu_e*snGrad(phi)*magSf, i.e. drift velocity = Z*mu*E is ANTIPARALLEL to E for electrons — electrons move toward higher potential, i.e. toward the ANODE. The current-source regulator is `dV = dt(I_set - I_cond)/(C + |g|dt)` and I_cond is NEGATIVE on a negative electrode (plasmaExternalCircuit.C:273, 506-522).
**Evidence:** src/models/plasmaModels/plasmaTransport/plasmaTransportModels/driftDiffusion/driftDiffusion.C:78-92; src/models/electromagnetics/electrostaticModels/singleRegionPoisson/singleRegionPoisson.C:39

**Rule:** First sanity check on ANY result: confirm the electron density front moves toward the anode and the positive-ion front toward the cathode. If they are swapped, suspect the charge number or the phiE sign, not the physics.

**Cost:** A sign error here would invert every space-charge structure while still 'running'.


### QUASINEUTRALITY IS A CATASTROPHIC-CANCELLATION AMPLIFIER, measured at ~1000x. In the Grubert bulk the plasma is quasi-neutral to 6 ppm (net/n_e = 5.96e-06), so rho = e(n_Arp + n_Ar2p - n_e) is a near-total cancellation — AND rho IS the Poisson source. A tolerance-IMMUNE mesh-asymmetry-induced Ey of 3.0e-07 of |E| was amplified through rho into a 25% n_e / 41% nEps_e / 23% chargeDensity lateral asymmetry, then grown by the ionisation feedback which e-folds the perturbation every 1.46 ns against the mean's 2.32 ns (1.59x faster). Under Newton the same effect shows as chargeDensity being the loosest field in a Picard-vs-Newton comparison (8.197e-05 max rel diff against 4.755e-06 for n_e).
**Evidence:** docs/design/grubert-lateral-asymmetry.md (link 4 of the chain, with the decisive testPoissonSymmetry experiment); docs/design/newton-outer-solver-design.md:1056-1076

**Rule:** Expect chargeDensity to be the loosest-agreeing field in any comparison and do not read that as an error. For lateral asymmetry: use NY = 1 for 1-D validation (no antisymmetric mode can exist) and normal refinement for 2-D/3-D (Ey ~ NY^-1.7). Do NOT tighten the Poisson tolerance — it improves phi and Ex but leaves Ey completely unchanged, and GAMG at 1e-16 costs 2000 iterations for nothing.

**Cost:** A 0.69% phi asymmetry becomes a 29% field asymmetry through differentiation (41x), then grows exponentially.


### MAGNITUDE SCALES, and the trap they set. A streamer channel/head is ~1.14e20 m^-3 (peak n_e reported 4.7e18 to 1.5e19 at ignition on Grubert, ~1e20 for a developed streamer head); the residual anode wall-flux layer under `includeDriftFlux false` is ~1e14 m^-3, SIX orders below the channel, and it is MODEL PHYSICS, not a bug. A global max over the domain has locked onto the wrong feature three separate times: a '112x physics difference' in peak n_e that was a boundary feature 9 mm from the streamer; a '23x stronger streamer' that was the anode surface; an 'anode layer is quasineutral, ratio 1.1' that divided an electron max at the anode by an ion max 2 mm away (locally the ratio is 4.9).
**Evidence:** memory wall-flux-two-models; memory never-read-a-number-without-its-control (instances 1-3); memory lfa-vs-lmea-2ns-results (trap 1, with the --head-min/--head-max windowing fix)

**Rule:** PRINT WHERE THE EXTREMUM IS, always, in coordinates. Window every head metric (--head-min/--head-max) and report the EXCLUDED maximum rather than silently masking it. A max over a region is not a local property on a mesh with features four decades apart nine millimetres apart.

**Cost:** Three would-be-published findings were artefacts of an unlocated global maximum.


### Te BOUNDS AND THE CLAMP. `meanEnergyMin` defaults to 1.5*kB_eV*Tgas (~0.039 eV at 300 K); `meanEnergyMax` must be set explicitly or FALLS BACK TO 100 eV while the Grubert argon tables run to 2644 eV — so a case that does not set it silently discards the top of its own tables. A reported Te,max of 1762.9 K/eV is the tell: 1762.9 = (2/3)*2644.46, i.e. the mean-energy clamp pinned. Healthy bulk values: the r=0 Grubert arm gives Te,max 11.58 eV / Te,min 2.66 eV; the LMEA streamer runs 7.76 -> 12.35 eV.
**Evidence:** localEnergyEnergyModel.C:249-257, 535, 557-559; localEnergyEnergyModel.H:263-270; memory grubert-anode-electron-reflection-runaway; memory lfa-vs-lmea-2ns-results

**Rule:** Check whether meanE is pinned at its clamp before interpreting any Te. `T_e = (2/3)*meanE/k_B`, so Te[K] = 7736*meanE[eV]; if Te,max equals (2/3)*meanEnergyMax you are reading the clamp, not the physics. Also note readings taken on NON-CONVERGED steps are not physics.

**Cost:** n_e ~3e18 / Te ~1762 eV were reported as results when they were clamp values on non-converged steps.


### FIELDS WITH A TRIVIAL SOLUTION TO CHECK AGAINST — the 'is this field sitting on its own trivial solution?' invariant. (a) S_iz sitting at ~0 was read as 'no ionisation yet' for most of a day; checked against the mechanism table at the SAME cell and time it was wrong by 250 ORDERS OF MAGNITUDE (~1.9e-234 against an expected 7.8e19 at meanE~23 eV, n_e~9.8e10). (b) n_e min sitting at exactly the floor means the field is undeveloped, not converged. (c) The vacuum/Laplace field is the control for any Poisson result; reducedE is left at ZERO where there is no background gas, deliberately and with a warning (singleRegionPoisson.C:81-108). (d) meanE at the cold floor at a wall is what breaks the LFA/LMEA mobility cancellation (wall/interior mu: LFA 1.000, LMEA 1.722).
**Evidence:** memory s_iz-frozen-under-option4-stiff-chemistry; memory density-floor-is-a-source; singleRegionPoisson.C:81-108; memory lfa-vs-lmea-2ns-results

**Rule:** For every diagnostic field, ask what its trivial value is (0, the floor, the clamp, the initial condition) and check it is NOT sitting there before drawing any conclusion. Cross-check S_iz against the mechanism table's own k at the same cell and time.

**Cost:** Any case using LMEA + `energySource chemistry` + `solver adaptive*` predating 2026-09-09 has broken S_iz/k_eff/alpha, which the AMR criterion, the photoionisation source and the Coulomb-heating term all read.


### SCHARFETTER-GUMMEL: exact with no source, first order with one. Measured on verification/fluxScheme1D, steady 1D DDR at Pe=1e4: with s=0 the SG error is 5.4e-16 to 9.7e-15 (machine precision); with s=1 it is 2.500e-6 to 3.125e-7, a UNIFORM offset (L2 == Linf to all digits) equal to the source flux dropped over half a cell — S*h^2/(8D) at low grid Peclet (O(h^2)), S*h/(2v) at high (O(h)), both confirmed to six digits. Measured order 1.00 over five refinements from P_grid 500 down to 15.6. The Complete Flux Scheme restores second order exactly there (order 2.00-2.02, 3390x more accurate at N=640).
**Evidence:** docs/design/flux-schemes-theory-and-implementation.md sections 2.2, 3.5

**Rule:** SG's first-order behaviour and its source-blindness are ONE defect seen from two angles. Do not report an asymptotic order without stating P_grid — h->0 drives P_grid->0 by construction, so a convergence table at domain Pe <= 100 never samples the regime where SG degrades.

**Cost:** Testing CFS at low Peclet showed nothing; and a uniform source makes div(Gamma^i)==0 so CFS collapses onto SG in the interior (measured 0.2% apart) — the uniform source chosen to ISOLATE the source term is the least favourable case for the scheme that fixes it.


### SG IS THE WRONG LEVER FOR JFNK CONDITIONING — refuted by measurement. Only the scheme varied: dt=1e-12 SG FAIL (as standard does); dt=5e-13 standard 20/20 converged; dt=5e-13 SG FAIL at step 1. SG's Bernoulli flux introduces exponential dependence on the local Peclet number, making F MORE nonlinear and harder to finite-difference. The same reasoning, plus COMSOL's own guide ('This makes it more numerically stable but INCREASES THE NONLINEARITY of the equation system'), is why log(n) is expected to make the current blocker worse, not better. SG and CFS are also UNBOUNDED and can undershoot a density negative; Pasolari & Kourtzanidis (arXiv 2607.05137) find SG 'excessively diffusive' on a positive streamer and recommend ROUNDF.
**Evidence:** docs/design/newton-outer-solver-design.md:1361-1386; docs/design/flux-schemes-theory-and-implementation.md section 7

**Rule:** SG remains the right choice for positivity and monotonicity in a drift-dominated cell (and is the remedy the wall-flux singularity error text itself names), but never propose it to improve JFNK conditioning. Neither SG nor CFS should become a default: neither verification bed has extrema, so boundedness is never exercised.

**Cost:** Re-running a refuted experiment; and adopting an unbounded scheme as default on smooth-solution evidence.


### `fluxScheme` IS READ FROM THE TRANSPORT MODEL, NOT FROM PATCH DICTIONARIES. Editing `fluxScheme` into the patch dictionaries of a restart time directory is silently stripped at startup — the BC reads `ddModel.fluxScheme()`. The real knob is per-species: `constant/plasmaSpeciesProperties: driftDiffusionCoeffs { fluxScheme $driftDiffusionFluxScheme; }` fed from `configuration/config: driftDiffusionFluxScheme standard|ScharfetterGummel|CompleteFlux`. The electron's setting propagates automatically to derived ions (`ionFluxScheme` defaults to the electron's value) and to the LMEA energy equation.
**Evidence:** docs/design/newton-ignition-experiments.md:674-691 (section 22); docs/design/flux-schemes-theory-and-implementation.md section 5

**Rule:** Confirm the scheme from the run's OWN species table at startup before trusting any A/B. Changing it changes the transport DISCRETISATION everywhere, not just at the wall, so it is a physics-affecting change whose accuracy needs checking against the existing 2 ns results.

**Cost:** A run looked like it was testing SG and was testing nothing.


### `zeroGradient` ON AN ELECTRODE IS NOT NO-FLUX. It zeroes only the diffusion term; the drift term still carries carriers through the patch. Under LFA a flat wall mobility makes the in- and out-fluxes cancel exactly; under LMEA mu = mu(meanE), meanE falls to the cold floor at the wall and the table returns a LARGER mu there (wall/interior mu: LFA 1.000, LMEA 1.722), breaking the cancellation. Verified at matched time t~5e-10: n_e 4.107e18 -> 1.0000e13, n_e/n_pos 19196 -> 1.0, chargeDensity -0.66 -> -3.9e-08 C/m^3 after replacing it with `electronDDWallFluxMixed`/`energyDDWallFluxMixed` (commit 3adbb00). It was never a sheath: 16300 electrons per positive ion, and ionisation makes PAIRS.
**Evidence:** docs/theory/scharfetter_gummel.md section 3.5; memory lfa-vs-lmea-2ns-results; memory electrode-bc-zerogradient-trap

**Rule:** Check any case pulled from an archive for `electronDDWallFluxMixed` before using it — do not assume. And when a near-electrode population appears, check the electron/positive-ion ratio: ionisation makes pairs, so a ratio of 19196 is a BC artefact and a ratio near 1 with both signs of chargeDensity is a real sheath.

**Cost:** Six of eight archived LMEA cases predate the fix and carry the artefact.


### THE WALL-FLUX SINGULARITY THAT ENDED THE NEWTON IGNITION RUN, and the two legitimate wall models behind it. `electronDDWallFluxMixed` imposes a total flux n*(uDrift_n + uEff) with `includeDriftFlux` selecting the model: false (the electron DEFAULT) gives uEff = u_th - uDrift_n, thermal only; true gives uEff = u_th + max(0,-uDrift_n), Hagelaar & Kroesen (2000). Both are correct. The valueFraction f = uEff/(D/delta + uEff) is SINGULAR when the drift into the wall exceeds the thermal speed by more than D/delta — the error text names the remedy: `fluxScheme ScharfetterGummel`, whose denominator D/delta*Bern(Pe) + uAbs has both terms non-negative for any reflection coefficient r and cannot invert. Measured margin on the streamer anode cell at 0.5 ns: u_th 1.0225e5, uDrift_n 1.0113e5, D/delta 4.44e4.
**Evidence:** memory wall-flux-two-models; docs/design/newton-ignition-experiments.md:636-673 (section 21)

**Rule:** When a run dies with 'wall-flux condition ... has become singular', that is a BOUNDARY-CONDITION MODEL limit, not a solver failure — switch the flux scheme or refine near the wall, and do not attribute it to the outer solver.

**Cost:** The Newton verdict run stopped at t=1.949675e-6 on this, and reading it as a Newton convergence failure would have been wrong.


### THE HONEST NEWTON-vs-PICARD VERDICT AT IGNITION: Newton did NOT beat Picard in time reached. Picard died at t = 1.971771e-06 (dt floor, 150/150 correctors, dt 1.7986e-14, a 72,000x collapse, n_e max 4.61e18). Newton stopped earlier at t = 1.949675e-06 but on a different cause (wall-flux BC singular) with dt 9.74e-14 (5.4x larger), a ~2,000x collapse (36x less severe), and n_e max 1.50e19 (3x higher). What this establishes is that the SOLVER is no longer what stops the run — a boundary-condition scheme is.
**Evidence:** docs/design/newton-ignition-experiments.md:636-673 (section 21)

**Rule:** State both the time reached AND the cause of stopping AND the dt collapse factor. 'Newton beats Picard' is not supported by this data and must not be claimed.

**Cost:** An overclaim that a reviewer or the next session would have to retract.


### A TRUE `ddtSchemes steadyState` MODE SEGFAULTS STRUCTURALLY, and the reason generalises: ddt IS the diagonal of a segregated species-transport equation. With steadyState the V/dt diagonal is gone and the implicit chemistry loss contributes nothing either (log: 'implicitRate, 2000 active cells, max(L*dt) = 0'), so a species with no loss term has NO diagonal at all — Gauss-Seidel divides by zero on the first sweep. simpleFoam can drop ddt because SIMPLE's pressure equation and under-relaxation supply diagonal dominance; a species equation has no such mechanism, and OpenFOAM's equation relaxation only scales a diagonal that must already exist. The working alternative TODAY, with no code changes: ddtSchemes backward (BDF2, 5.5x faster than Euler here) + adjustTimeStep false + fixed deltaT 1e-12 + currentSource + NO relaxationFactors + omit onNonConvergence — measured 9472 consecutive converged steps, 0 non-converged, 4 correctors/step against a 150 cap.
**Evidence:** docs/design/steady-mode-spec.md sections 1.1-1.4

**Rule:** Do not propose `ddtSchemes steadyState` for this solver. Propose the pseudo-transient recipe, and treat its deltaT as a RELAXATION PARAMETER, not physics. Note `onNonConvergence retryStep` is REJECTED when adjustTimeStep is off, since it retries by shortening deltaT.

**Cost:** A SIGFPE on iteration 1, twice measured.


### THE SEMI-IMPLICIT POISSON OPERATOR IS THE SCHUR COMPLEMENT, derived: eliminating the transport block gives S_f = A_ff - A_ft*A_tt^-1*A_tf, and with A_tt ~ I/dt, A_ft = q_i (diagonal) and A_tf = -div(Z*mu*n*grad(.)), this equals div((eps + dt*sigma)grad(.)) term for term, since sigma = q*Z*mu*n. The reaction terms cancel EXACTLY to all orders in dt because charge conservation makes q^T a LEFT NULL VECTOR of the reaction Jacobian (sum_i q_i*S_i(n) = 0 identically, so q^T J_S = 0). Only the TRANSPORT part of A_tt is approximate. FALSIFIABLE PREDICTION: the operator's quality as a Schur preconditioner should degrade with the transport CFL/diffusion numbers and be INDIFFERENT to k_eff.
**Evidence:** docs/design/schur-semiimplicit-poisson-preconditioner.md (Status: DERIVED, not yet implemented)

**Rule:** Resolve the apparent contradiction correctly: the semi-implicit Poisson is WRONG in the residual (double-counts rho^{n+1}) and RIGHT in the preconditioner (it IS the Schur complement). Cite it as a derivation, not a capability — it is not implemented.

**Cost:** An earlier version of the same note claimed the reaction rate should appear in the operator as div((eps + sigma/(1/dt + k_eff))grad); that was wrong and was corrected by the conservation argument.


### STRUCTURALLY ABSENT FROM THE MODEL — do not claim physics the closure cannot represent. (1) NO electron-electron (Coulomb) collisions in SoPlasma's table pipeline: SoEEDF implements them (src/CoulombTerms.C, Rockwood isotropic + Hagelaar 2016 anisotropic, `coulomb` option, validated) but the SoPlasma generator never invokes them and never passes local n_e/N — tables are 1D in E/N or mean energy only. This is an INTEGRATION GAP (a second table axis), not a research problem. (2) NO spatially non-local trapping: a cold electron confined in a potential well escaping via a rare Coulomb kick depends on the SHAPE of the potential over a finite extent and on the electron's history — NO local closure, however many local parameters it is given, can represent it by construction. (3) LMEA is SINGLE-GROUP, so Eliseev's H_secondary Coulomb-heating term (correctly implemented and validated against the paper's own numbers) cannot help: S_iz is computed from the SAME local meanE Joule heating sets, so it structurally cannot be large where the local field has collapsed — the one regime the term exists to correct. Measured: epsilon_eff*S_iz came out 7-8 orders below Psrc at every step checked. (4) Non-local kinetic theory is NOT a clean fit either: Kortshagen's argon criterion gives N0*R = 2.414e22 m^-3 cm for Grubert, 8x above the nonlocal threshold and 4x below the local one — squarely intermediate, where the full spatially-dependent kinetic PDE is needed.
**Evidence:** docs/design/electron-electron-collisions-gap.md; docs/design/nonlocal-kinetics-assessment.md; memory no-electron-electron-collisions-in-soplasma; memory coulomb-heating-term-structurally-mismatched-to-lmea

**Rule:** Never claim EEDF-shape effects, ionisation-degree-dependent coefficients, or negative-glow electron trapping from this model. Any case resembling a negative glow (Grubert 2009, the Carlsson/JC-PIC benchmark) may plateau at the wrong bulk Te/density regardless of run duration or numerics.

**Cost:** A whole session implementing and validating a Coulomb-heating term that has nothing to act on.


### THE VERIFICATION LADDER THAT DOES NOT EXIST YET, and why it matters: validation/ holds 60 cases, 58 of them Grubert variants, and NONE is a verification case with a known answer. The highest-value Tier-1 analytic tests, in the recommended order, are: D ambipolar diffusion (exact parabola n = (G0*l^2/8Da)[1-(2x/l)^2]+1; the natural fluxScheme test since `standard` and `ScharfetterGummel` must BOTH reproduce it); H Teunissen's semi-implicit artifact test (1D, 10 mm, dx=20 um, 10 kV, SOURCE ZERO, mu_e=0.03, De=0.1, n_e=n_p=1e20 in 4-6 mm — his result is spuriously LARGER field peaks at the ionised region's boundaries, which is our anode signature with no chemistry to blame); C Debye screening (exact V = V0*exp(-|x|/lambda_De)); G convergence order (Teunissen measures the semi-implicit field treatment at roughly FIRST order); K normal cathode fall similarity (V_n, (pd_c)_n, j_n/p^2 from Raizer ch. 8); L Carlsson/JC-PIC He benchmark (3.5 Torr, 0.62 cm gap, -211 V, gamma = 0.28, 2 eV emitted electrons, matched by three kinetic codes to ~2% on E0).
**Evidence:** docs/design/verification-map.md (from DeChant 2023 PSST 32 044006 and Teunissen 2020 PSST 29 015010)

**Rule:** When asked to validate a numerical change, prefer a Tier-1 bed with a known answer over a discharge case. Note MMS verifies the code solves the equations it CLAIMS to solve and cannot catch a wrong equation — the Hagelaar eq. (6.15) defect was a wrong equation faithfully implemented and MMS would have passed it.

**Cost:** Every defect found on 2026-09-07 was found by hand-instrumenting a diverging discharge; each would have been caught faster, or at all, by one of these.


### THE STEP-DISCARD INVARIANT, stated by the user as crucial: every piece of solution state must be IDEMPOTENT WITHIN A TIMESTEP — a step can be thrown away and re-run and the result must not depend on the discarded attempt. The mechanism: the time INDEX is held fixed across retries (soPlasmaFoam.C ~line 327), oldTime() rotates on the index, `storeOldTimes()` keys on `timeIndex_ != time().timeIndex()`, and `GeometricField::oldTime()` creates a level on demand and `storeOldTime()` recurses, so a third BDF2 level comes free. The authoritative list of what must be restored is `discardStep()`: species densities (plasmaTransport::discardStep, always), nEps_e (localEnergyEnergyModel::discardStep, LMEA only), T_gas (plasmaEnergy::discardStep, when solveGasEnergy_).
**Evidence:** memory step-discard-invariant; memory temporal-error-control-design

**Rule:** Any new field, cache or counter must be added to the discard path or be idempotent by construction. Enumerate from `discardStep()`, never hardcode. Under Newton specifically, peak-hold diagnostics (clampRaw_, advisoryLpeak_, chemStiffnessPeak_, chemPicardPeak_) are mutated once per matvec and are therefore MEANINGLESS.

**Cost:** A Phase D exclusion ('the clamp is not binding') came from one of those peak-hold diagnostics and is not sound evidence.


### TOLERANCE POLICY FOR ORDER STUDIES, measured with opposite answers for the two tolerance families. ODE `odeCoeffs absTol/relTol` = 1e-12/1e-10 is MANDATORY for LMEA order studies and is a THRESHOLD, not a gradient: 1e-4, 1e-6 and 1e-8 all give p = 0.748 to 5-6 s.f., and only 1e-10 moves it (1.004). NOT needed for LFA (the control measures 1.945 at the shipped 1e-4). Conversely linear-solver 1e-13 / outer-gate 1e-11 is measured WASTE: identical to 6-7 s.f. for +34% wall clock (487 vs 363 s). Baseline is linear 1e-10 / gate 1e-8. Counter-intuitive cost structure: the FINE dt arm is cheaper PER STEP (0.21 s/step at dt=1e-13 vs 0.73 s/step at 4e-13) because stiff-ODE work scales with step size, so shortening endTime is the biggest lever.
**Evidence:** memory order-study-tolerance-policy

**Rule:** For an LMEA order study set odeCoeffs 1e-12/1e-10 and leave the linear/outer tolerances at 1e-10/1e-8. Below the ODE threshold the integrator's own error swamps the O(dt) term under test, so every future fix would look like a failure.

**Cost:** A loosened study reports ~0.75 whatever the fix did.


### THE COLD-START HANDOVER AND ITS TWO BUGS. Newton is degenerate at a cold start: the bounded solver's own report read 'REDUCED system (2000 of 10000 unknowns)' — 8000 of 10000 unknowns pinned at their bound, leaving only the potential block free. The fix is a Picard warm-up that hands over the moment NO species sits on its minNumberDensity floor; measured, the handover fires at t = 3e-12, i.e. after only TWO Picard steps, and immediately on an established restart. The two bugs it shipped with: `anySpeciesOnFloor()` originally tested `gMin(n) <= floor` ('is ANY cell at the floor?'), which is true FOREVER in a real discharge because the quiescent far field always rests on the floor — on the ballast case Arp's min sits at exactly 1e11 for the whole run while its max reaches 5.7e18, so Newton could never be handed to on the one case it was wanted for; now `gMax(n) <= floor`. And the handover did not LATCH, flipping back to Picard silently whenever a species touched its floor.
**Evidence:** docs/design/newton-outer-solver-design.md:1640-1681; docs/design/newton-ignition-experiments.md:143-155

**Rule:** Verify a guard in BOTH directions — that it fires when it should AND stays silent when it should. A 'Phase D 40,000-step Newton success' was recorded while Newton was actually off for 21,174 of those steps.

**Cost:** Four separate 'measurements that turned out to be measuring nothing' in one day, now rule 42.


### NEWTON'S CURRENT RESTRICTIONS, verified: it supports only `singleRegionPoisson` (needleDBD refuses because its regionProperties declares a `dielectric` MESH REGION) and only `driftDiffusion` transport models (FatalError at construction otherwise). Surface charge ITSELF works under Newton — proved with `thinDielectricOnElectrode` (pmma, 100 um, Vb=0) on the coarse streamer bed, exercising the surfCharge field, the chargingSurface wall flux and the Robin potential condition. The photoionisation refusal guard was keyed on `photoionization_.valid()`, and that autoPtr is ALWAYS valid (a case with none holds the null object `noPhotoionization`, TypeName "none"), so the guard refused every case in existence until commit 4547e95.
**Evidence:** docs/CAPABILITIES.md section on surface charge under Newton; docs/design/newton-ignition-experiments.md:974-997 (section 26); memory newton-vs-picard-benchmark-state

**Rule:** Before offering Newton for a new case, check for multiRegionPoisson, non-driftDiffusion transport models, photoionisation and the legacy Townsend `!rates_` source. Dielectric cases under Newton are currently UNGUARDED — a correctness hole, not an enhancement.

**Cost:** Newton silently solving a different system than the one the case describes.


### THE LEGACY TOWNSEND FIT IS NOT MISCALIBRATED — a claim of the author's own that was committed twice and then corrected. `alpha = A*exp(-2.73e7/Emag)` expects RAW V/m and is correct for air at 1 atm: alpha = 19 cm^-1 at the 30 kV/cm breakdown field, 1.06e3 cm^-1 at 100 kV/cm, mu = 2.398*E^-0.26 -> 0.036 m^2/V/s at 1e7 V/m, eta = 3.4 cm^-1. The fit is gas- and pressure-SPECIFIC, not miscalibrated. The 2026-09-09 measurement that motivated gating it off (S_iz ~250 orders low) stands, but because it was taken in a 100 Pa ARGON glow — outside the fit's regime.
**Evidence:** docs/design/newton-ignition-experiments.md:974-997 (section 26); commit ef0a59b

**Rule:** Distinguish 'a correlation applied outside its regime' from 'a correlation with wrong units'. Check the gas and pressure the fit was built for before concluding a normalisation error.

**Cost:** A wrong comment shipped twice in plasmaTransport.C and would have misled the next reader.


## Traps
- A DESIGN doc cited as a capability. Detector: `head -5 <doc>` — adaptive-dt-design, steady-and-stability-design, numerics-generator-plan all say 'Nothing implemented'; schur-semiimplicit-poisson-preconditioner says 'DERIVED, not yet implemented'.
- A swept knob that never moved. Four incidents in one day, now rule 42: `-fieldsplit_0_pc_type hypre` (the splits are named `phi` and `transport`, there is no split 0); `maxDeltaT` in system/controlDict (plasmaTimeControl reads its own dict and silently overrides with default GREAT); `sed s/outerSolver.*/newton/` on a case with no outerSolver key (two 'Newton' arms ran Picard for 24,000 steps); `fluxScheme` edited into patch dictionaries (the BC reads ddModel.fluxScheme()). Detector: `-options_left`, `-ksp_view`, and the run's own start-up banner — confirm the value CHANGED before believing any comparison.
- Comparing arms at unequal sample size or from different start times. Four misreadings in one session, in both directions: `ew_nochem` appeared to beat `ew_v2` (3.416e-10 vs 2.955e-10) purely because it averaged 34 steps against 93. Detector: print N and the common start time beside every mean; truncate both arms to the same N.
- Reading the STEP COUNT instead of the ACHIEVED dt. Bounded Newton: 360 steps, 0 fatals, 2383 SNES calls — and advanced 4e-11 total at dt = 5e-15. Detector: always report advance-per-step and mean dt, never steps.
- Reading a tail statistic off a CAPPED distribution. `kspMaxIt 100` truncated the KSP iteration distribution and made it look bimodal ('a preconditioner cliff'); uncapped it is median 4, p90 14, max 115 — unremarkable. Detector: count how many solves land EXACTLY on the cap (38 of 4,235 did).
- `DIVERGED_LINEAR_SOLVE` read as a solver failure when it is a budget expiring. 39 solves needed 104-115 iterations and all converged once the cap was raised. Detector: raise kspMaxIt to 1000 and re-count; cost was +2.5% wall clock for a 10.5x drop in failures.
- A global max locking onto the wrong feature. Three published-grade errors: '112x physics difference' (a boundary feature 9 mm away), '23x stronger streamer' (the anode surface), 'quasineutral, ratio 1.1' (an electron max at the anode divided by an ion max 2 mm away; locally 4.9). Detector: print the extremum's COORDINATES and window every head metric, reporting the excluded maximum.
- A field sitting on its own trivial solution, read as physics. S_iz sat ~1.9e-234 for a day and was read as 'no ionisation yet'; against the mechanism table at the same cell and time it was 250 ORDERS low. Detector: for every diagnostic field, compute what its trivial value is (0, floor, clamp, initial) and check it is not sitting there.
- `n_e min` IS the density floor, so max/min is a ratio to the floor. max/min = 1.22 means the whole field is within 22% of the floor — an undeveloped state, not a result. Detector: compare n_e min against minNumberDensity before quoting any dynamic range.
- `floor hits` in the failure block counts TIMESTEP-floor hits (minDeltaT), NOT density clipping (plasmaTimeControl.C:1557). Misreading it produced a whole false hypothesis about CFS driving densities negative. Detector: read the code that increments any counter before inferring a mechanism from it.
- Te pinned at the mean-energy clamp, read as physics. Te,max 1762.9 = (2/3)*2644.46 in three separate runs. And `meanEnergyMax` falls back to 100 eV when not stated while the argon tables reach 2644 eV. Detector: check whether Te,max equals (2/3)*meanEnergyMax, and whether the step it came from converged.
- A clamp or one-shot initialisation inside or around the residual makes F(u)=0 unreachable. `clampNumberDensities()` after the solve; `updateDerived()` writing a clamp back into nEps_; a one-shot LFA seed firing on residual call #1 only (172.46 vs 0.00086, a factor of 2e5). Detector: under PCNONE the preconditioned and true residual norms must be identical — if they are not, the function changed between evaluations.
- A matrix-free Pmat makes PCFIELDSPLIT degrade SILENTLY to `PC type: none`; `-fieldsplit_*_pc_type hypre` on a shell densely probes it, ~10,000 residual evaluations per setup — the 1500x defect reintroduced. Detector: `-ksp_view` and check each split's PC type is what you asked for.
- An inner Schur solve with a varying iteration count makes S a non-fixed operator and breaks the outer Krylov method. The phi inner solve took 55,65,71,71,82,85 iterations on successive applications, then SIGFPE in MatMult_SeqAIJ. Detector: after the fix the phi solve reads '1 iteration x 692,813 applications'; any spread in that count is the bug.
- `PETSC_INFINITY` is PETSC_MAX_REAL/4 ~ 4.5e307, so xu - xl ~ 9e307 and further arithmetic overflows to inf, which OpenFOAM TRAPS under FOAM_SIGFPE. Use a finite 1e30. Detector: a SIGFPE immediately on entering a bounded solve.
- `zeroGradient` on an electrode is NOT no-flux — it zeroes only diffusion while drift carries carriers through. Under LMEA the wall mobility is 1.722x the interior value, breaking the cancellation LFA gets for free. Detector: check n_e/n_pos near the electrode; ionisation makes PAIRS, so 19196:1 is a BC artefact, not a sheath.
- `wmake` in src/numerics exits 0 having silently NOT built libplasmaNewtonSolverPETSc (its Make/files is one directory down). Detector: check the library mtime and grep the binary for a string you just added — verify the artefact, not the exit code.
- A guard verified in only one direction. The photoionisation refusal was keyed on `photoionization_.valid()`, always true, so it refused every case in existence; it had been 'verified' only by making it fire. Detector: test that a guard STAYS SILENT when it should, as well as that it fires.
- A restart that is not the restart you think. A run REWRITES the snapshot it restarts from including `uniform/time`'s deltaT, so a second restart from 'the same' snapshot starts at a different dt (the FIELDS stay byte-identical, so it is invisible unless the time state is checked). And `timePrecision 6` cannot address a directory OpenFOAM wrote with 7 digits; `startFrom latestTime` does NOT help, raising timePrecision does.
- Peak-hold diagnostics (clampRaw_, advisoryLpeak_, chemStiffnessPeak_, chemPicardPeak_) are mutated once per matvec and are MEANINGLESS in Newton mode. A Phase D exclusion ('the clamp is not binding') rested on one of them. Detector: check whether a diagnostic is written inside residualCallback before citing it under Newton.
- Manual `relaxationFactors` added on top of the adaptive Aitken relaxation that is already running. Detector: 9472/9472 converged became 0/10 converged. The generator must never emit them.

## Open questions
- Whether Newton pays at all remains UNANSWERED. The stated benchmark — a WARM-STARTED ladder from a developed streamer at 449k, measuring the largest dt each solver survives at equal accuracy with Picard-era limiters off — has not been run, and every prior ladder is invalidated by the `kspMaxIt 100` cap. Break-even requires Newton to hold ~5e-10 against warm Picard's 5e-11.
- The cause of the 1.15M-cell Newton wall is still unknown. The 'starved Schur complement' hypothesis was REFUTED by measurement (2 vs 8 inner cycles: -5% SNES, -4% Krylov, +14% wall clock).
- F is still not a smooth function of u at the 1e-6 level, so `mffdErr` ships as a documented workaround at 1e-5 rather than PETSc's default. Whether the remaining non-smoothness sources (the nEps clamp inside updateDerived, chemistry ODE step-sequence discontinuity) can be removed is the objective test, and it has not been done.
- `d(Psrc)/d(phi)` — Joule heating's response to the potential — is still missing from the energy row of the Pmat. It is direct and probably stronger than the retained drift term but is not laplacian-shaped, so it needs its own derivation. It remains a suspect for the LINE-SEARCH failures, but on the argument that an incomplete row gives a wrong direction, NOT on magnitude evidence (that claim was mis-assigned to a retry's assembly and corrected).
- Per-CELL scaling is untested: sX is still one scalar per FIELD while Arp spans 1e11-5.7e18 within its block. 42 usage sites, and a wrong refactor would silently corrupt every result recorded so far. A `perCellScaling` switch exists and defaults false.
- `maxIt` is state-dependent and unresolved: raising it 50->200 HURT at the easy pre-ignition state (-16%) while the deep-avalanche run fails with DIVERGED_MAX_IT at maxIt 50. A sweep AT the hard state has not been run.
- Pmat staleness is untested — it is assembled once per outer step and held for the whole SNES solve, and at ignition the state moves fast within a step.
- Chemistry accuracy under Newton is UNVERIFIED: the `adaptiveError` nOuterCorrectors > 1 guard is bypassed for `outerSolver newton` only.
- Whether the explicit-Poisson Newton trajectory is CORRECT and not merely self-consistent. Two explicit arms agree with each other to ~1% and the semi-implicit form differs by O(dt*sigma/eps) = 20-60%, but this has not been checked against the PUBLISHED benchmark, so 'better physics rather than merely different' is not yet established — and the AP proof is explicitly deferred until it is.
- The reaction-cancellation prediction for the Schur preconditioner is falsifiable and untested: iteration counts should degrade with the transport CFL/diffusion numbers and be INDIFFERENT to k_eff. A sweep raising the ionisation rate at fixed dt and mesh has not been run.
- Whether electron reflection r = 0.36 at the anode is the correct physics or the closure's r-handling is wrong. r is a PHYSICAL input; r = 0 is a diagnostic, not an established fix, and r = 0 only delayed the failure (1.125 -> 1.680 us).
- Grubert's own stated LMEA validity conditions (10a) nu_E,M >> dz*E and (10b) nu_E,M >> dt*E have NEVER been evaluated on this project's fields.
- Whether an LMEA analogue of Teunissen's parallel-diffusion artefact exists. He states the mechanism is LFA-specific; under LMEA ionisation keys on meanE rather than local E, so it should be damped by energy relaxation — but whether electrons carrying nEps into the high-field region reproduce it is open, and the f_eps diagnostic has not been computed.
- The wall-cell Poisson residual is nonzero even at r = 0 (6.4 against a required 0.029) — unexplained; a converged FV Poisson solve should satisfy it.
- Whether `chemSrcPrev_`/`chemOuterCount_` mutation during residual evaluation can make the residual history-dependent. Reviewed and cleared for the current code (the only consumer, `finalOuterIteration()`, is dead code, and a freeze test showed no difference to 8 digits), but it would become a Tier-0 defect the moment anything consumes `chemPicardChange_` in Newton mode.
