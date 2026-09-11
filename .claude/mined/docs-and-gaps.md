# docs-and-gaps

## Summary
Directive D1 ("documentation is part of done") obliges a sweep over a surface of 82 markdown files plus two large READMEs, 8 tutorial_info.md files, 5 smoke-bed dictionary sets and every tutorial case comment — concretely, 45 .md under /home/kkourtza/soplasma-scratch/docs (two of them zero-byte stubs), 37 .md under /home/kkourtza/Projects/SoEEDF/docs (5 in archive/), SoPlasma README.md (47204 B, 987 lines) and SoEEDF README.md (79657 B, 1158 lines) — both sizes confirmed. The one mechanised honesty check in the project, tools/checkConfigReference.py, is currently RED and covers almost nothing: it exits 1 today and reports "3 documented, 1 verified against source; 301 options read by src/", i.e. 294 undocumented keys and exactly ONE machine-verified key in the whole reference. Two whole-framework inventories exist and have already diverged: SoPlasma's docs/CAPABILITIES.md (committed 2026-09-11, 88 mentions of Newton/SNES/JFNK) versus SoEEDF's docs/framework-state.md (2026-09-05) and docs/capability-overview.md (2026-08-21), which contain ZERO mentions of the Newton/PETSc outer solver that is the project's largest recent piece of work — so SoEEDF's inventories are stale by an entire solver. docs/INDEX.md is the best PROGRESS.md model in either tree (it carries a per-document status column, a "single source of truth" owner table, and dated "linked here"/"moved" stamps), but it is itself incomplete: it never mentions CAPABILITIES.md, docs/design/, docs/reference/, rules-postmortems.md or CLAUDE.md, and SoPlasma has no index at all. The gap that matters most operationally is that every reusable procedure in this project — the rename sweep, the build gate, the pre-launch gate, the knob-moved proof, the snapshot protocol — exists only as prose inside a 397-line post-mortem file and ~128 memories, and is therefore re-derived (and re-failed) by hand each session.

## Facts

### The D1 documentation surface is 82 markdown files plus 2 READMEs: 45 .md under soplasma-scratch/docs (docs/design 22, docs/reference 6, docs/models 8, docs/getting_started 4, docs/theory 1, docs/simulationManuals 1, CAPABILITIES.md), and 37 .md under SoEEDF/docs (5 in archive/). Both README sizes confirmed: SoPlasma 47204 B / 987 lines, SoEEDF 79657 B / 1158 lines.
**Evidence:** find /home/kkourtza/soplasma-scratch/docs /home/kkourtza/Projects/SoEEDF/docs -name '*.md' | wc -l; wc -l -c on both README.md (987 47204 / 1158 79657)

**Rule:** Treat the D1 sweep set as: both README.md, docs/**/*.md in BOTH trees, SoEEDF/mechanisms/FORMAT.md, SoEEDF/{data,meshes,Literature}/README.md, soplasma-scratch/{PR-DESCRIPTION.md,CLAUDE.md,validation/TUNE_*.md,validation/RESOLUTION_TEST.md}, the 8 tutorials/**/tutorial_info.md, every tutorial and smoke/ dictionary INCLUDING comment text, and etc/boundaryRoles.

**Cost:** Rule 1: anything left behind is not merely stale, it is an INSTRUCTION a later session reads and follows.


### A naive recursive grep over the live SoPlasma tree does NOT complete a rename sweep: `grep -rn ... .` with --exclude-dir for build/ThirdParty/validation/polyMesh was killed at 120 s (exit 143). `git grep` scoped with the pathspec ':!validation' returns the same answer in 0.067 s — a >1800x difference.
**Evidence:** Background task b1712bn05 exit 143 after 120 s; `time git grep -nE 'maxVoltageRisePerStep' -- ':!validation'` = real 0m0.067s, 29 hits. Unscoped `git grep -- .` = real 1m37s (user 0.04 s — pure I/O stall on the deleted validation/ index entries under WSL2).

**Rule:** NEVER use `grep -r` for a key sweep in soplasma-scratch. Always `git grep -nE '<old>|<new>' -- ':!validation'`, run once per tree. Tracked file count is 4914, so git grep sees the whole real surface.

**Cost:** The sweep is abandoned or truncated mid-way, which is exactly how the smoke beds accumulated four rejected keys.


### tools/checkConfigReference.py is RED right now and its coverage is one key. Default mode exits 1 with: MISSING READER `uniformValue` (docs/reference/changeDictionary.md:33), STALE DEFAULT `defaultSEEC` (documented 0.001, source has '<derived; see the reader>'), 294 UNDOCUMENTED, and the footer `3 documented, 1 verified against source; 301 options read by src/`.
**Evidence:** cd /home/kkourtza/soplasma-scratch && python3 tools/checkConfigReference.py; REAL EXIT=1. Footer line verbatim. --strict EXIT=1, --case tutorials/plasma/soPlasmaFoam/needleDBD EXIT=1.

**Rule:** Do not present the checker as a passing gate. Either fix the two reported defects and adopt it as a pre-commit gate for docs/reference/, or record in PROGRESS.md that its current baseline is 'exit 1, 2 defects, 1/301 verified' so a later session can tell regression from baseline.

**Cost:** A gate that is already red teaches everyone to ignore its exit code, which is how it stops catching the ePotentialControls class of error it was written for.


### The reference docs barely use the convention the checker parses, so the check is near-vacuous. Only 7 level-3 headings in docs/reference/ carry a backticked key, and the two largest reference documents have ZERO: plasmaSimulationControls.md (10113 B) and fvSchemes-fvSolution.md (7195 B) have no conforming `### \`key\`` heading at all.
**Evidence:** grep -cE '^###[[:space:]]+.*`[A-Za-z_][A-Za-z0-9_]*`' per file: README.md 0, changeDictionary.md 2, configuration-config.md 3, controlDict.md 1, fvSchemes-fvSolution.md 0, plasmaSimulationControls.md 0.

**Rule:** When adding or editing any docs/reference/ entry, write it as `### \`keyName\`` followed before the next heading by either `**REQUIRED -- no default.**` or `**Default:** \`value\``. That single formatting rule is what converts prose into a checkable claim.

**Cost:** 301 options read by src/ and 1 verified; the reference can drift arbitrarily and the tool reports OK-shaped output.


### docs/reference/README.md links a document that does not exist: `regions-and-materials.md` is the listed reference for constant/regionProperties and constant/<region>/electricalProperties, and there is no such file.
**Evidence:** Link-resolution loop over docs/reference/README.md prints `BROKEN: regions-and-materials.md`; ls docs/reference/ shows only README, changeDictionary, configuration-config, controlDict, fvSchemes-fvSolution, plasmaSimulationControls.

**Rule:** Add a broken-link check to every doc sweep (one-line shell loop, see exact_commands). Two more reference docs are openly marked *pending*: plasmaSpeciesProperties and plasmaTransportProperties.

**Cost:** The reference README is the entry point organised 'by the file you edit'; a dead row sends the reader nowhere for the two dictionaries that carry permittivity and gamma.


### SoEEDF/docs/INDEX.md is the strongest PROGRESS.md model in the project — a status column per document, a 'Where each subject is documented (single source of truth)' owner table, dated stamps for additions and moves — but it is cross-tree blind. It contains ZERO mentions of CAPABILITIES.md, docs/design/, docs/reference/, rules-postmortems.md or CLAUDE.md, and it omits rules-postmortems.md, rules-consolidation-proposal.md and all 5 archive/*.md from its own listing.
**Evidence:** grep -n 'CAPABILITIES\|docs/design\|docs/reference\|CLAUDE.md\|rules-postmortems' docs/INDEX.md returns nothing. Orphan loop reports 7 docs present but unlisted. INDEX.md last commit 2026-09-10, i.e. the SAME DAY CAPABILITIES.md was created, and it still does not link it. Its 'On the SoPLASMA side' table names only theory/, models/, simulationManuals/, getting_started/ and tutorial_info.md.

**Rule:** Model PROGRESS.md on INDEX.md's three devices (status column, owner table, dated stamps) but make it span BOTH trees and include docs/design/, docs/reference/ and CAPABILITIES.md. Keep exactly one such file so the two do not diverge the way the inventories already have.

**Cost:** INDEX.md's own opening claim — 'Every document in the framework' and 'each fact has one home' — is already false, and it is the document new readers are told to start from.


### SoPlasma has NO index equivalent. docs/ contains only CAPABILITIES.md plus six subdirectories; the closest thing is CAPABILITIES.md, which explicitly says 'This file POINTS. It does not duplicate.' The only cross-reference from SoPlasma's README to an index is a GitHub URL to SoEEDF's INDEX.md at README.md:269.
**Evidence:** ls /home/kkourtza/soplasma-scratch/docs/ -> CAPABILITIES.md design getting_started models reference simulationManuals theory. grep -n 'CAPABILITIES\|INDEX.md' README.md -> only line 269 (a github.com/kkourtza/SoEEDF/blob/master/docs/INDEX.md link).

**Rule:** Make CAPABILITIES.md the authoritative live inventory (it is the newest and the only one that knows about Newton), and have PROGRESS.md and any SoEEDF index POINT at it rather than restate it.

**Cost:** 45 SoPlasma docs, 22 of them design notes, with no listing anywhere — setting-up-a-new-case.md was already found unfindable for its entire existence (INDEX.md notes it was 'Linked here 2026-09-04 -- it was unreferenced from any index').


### DUPLICATION AND STALENESS, measured: SoEEDF/docs/framework-state.md (2026-09-05) and docs/capability-overview.md (2026-08-21) both inventory SoPlasma and both contain ZERO occurrences of 'newton', 'snes' or 'jfnk'; soplasma-scratch/docs/CAPABILITIES.md (2026-09-11) contains 88. The Newton/PETSc outer solver — ~6900 lines per rule 38's post-mortem — is invisible in the SoEEDF-side inventories.
**Evidence:** grep -cin 'newton\|snes\|jfnk': framework-state.md 0, capability-overview.md 0, CAPABILITIES.md 88. git log -1 --format=%cs: framework-state.md 2026-09-05, capability-overview.md 2026-08-21, CAPABILITIES.md 2026-09-11.

**Rule:** Declare soplasma-scratch/docs/CAPABILITIES.md the single authoritative capability inventory. Demote framework-state.md to a dated audit artefact (retitle its H1 to keep '2026-09-04' visible) and capability-overview.md to a commercial narrative that cites CAPABILITIES.md rather than listing capabilities itself.

**Cost:** Rule B4/R0 exists because a session proposed already-implemented work; two inventories, one of them a week behind and missing the flagship solver, guarantee that failure repeats.


### framework-state.md §5 is a working stale-document REGISTER — a table of 11 documents with 'what it says' vs 'what is true' — and it is partly self-cleaning: three of its entries have since been fixed IN PLACE as it prescribes (AMR.md now carries '(CORRECTED 2026-09-04: this sentence named plasmaDielectricFoam, a solver that no longer...)' at line 13; ddSolidSurfaceFlux.md opens with '# SUPERSEDED 2026-09-04 — read this first'; validation-summary.md:419 keeps the old 'Not modelled at all yet' text inside a parenthetical correction).
**Evidence:** grep -n 'plasmaDielectricFoam' docs/simulationManuals/AMR.md:13; head -12 docs/models/.../ddSolidSurfaceFlux.md; grep -n 'Not modelled at all yet' validation-summary.md:419.

**Rule:** Adopt §5's exact shape as a standing PROGRESS.md section ('document | what it says | what is true'), and keep rule 19's in-place marking: mark SUPERSEDED/CORRECTED with the date, never delete. A finding that quietly disappears leaves no way to recognise its stale copies elsewhere.

**Cost:** Without the register, the only mechanism that has ever caught these is a whole-framework sweep 'from source rather than from the other documents', which cost a full session.


### CAPABILITIES.md §1 instructs the reader to 'check the `Status:` line' on docs/design/*.md, but 14 of 22 design documents have no occurrence of the word 'status' in their first 30 lines — including the two largest and most consequential, newton-outer-solver-design.md (97615 B) and newton-ignition-experiments.md (53030 B).
**Evidence:** Loop over docs/design/*.md: NO 'status' IN FIRST 30 LINES for avalanche-vs-ion-transit, electron-electron-collisions-gap, electron-energy-balance, external-circuit-plan, flux-schemes-theory-and-implementation, grubert-instability-analysis, grubert-lateral-asymmetry, mesh-convergence-plan, newton-ignition-experiments, newton-outer-solver-design, nonlocal-kinetics-assessment, parallel-blockers, steady-mode-spec, verification-map. The 8 that do carry one are explicit and dated, e.g. stationary-solver-plan.md '**Status: PARTLY SUPERSEDED, 2026-09-11...'.

**Rule:** Require a dated `Status:` line in the first 10 lines of every docs/design/*.md, and make the one-line loop in exact_commands a pre-commit check. Until then, treat an absent Status line as UNKNOWN, never as 'implemented'.

**Cost:** A 97 KB design note with no status is read as a description of what exists; framework-state.md §5 already caught four plan documents still marked unimplemented after the work had landed.


### SoEEDF's standalone (no-OpenFOAM) constraint is enforced by the build, not by a test: BOLTZMANN_NO_OPENFOAM is set globally at CMakeLists.txt:10 (`add_compile_options(-O2 -DNDEBUG -DBOLTZMANN_NO_OPENFOAM)`) and again per target on all 10 targets (lines 33,38,43,48,53,58,63,68,73,78), and it guards three blocks in src/BoltzmannCoupling.C (lines 8, 68, 111) and two in include/BoltzmannCoupling.H (lines 17, 68).
**Evidence:** grep -rn 'BOLTZMANN_NO_OPENFOAM' over SoEEDF src/ include/ CMakeLists.txt docs/ README.md; documented at docs/installation.md:83 and README.md:447.

**Rule:** The CMake build IS the standalone test — `cmake --build build` failing is the only detector. Rule 40 is scoped OUT of SoEEDF for exactly this reason; never introduce an OpenFOAM idiom, type or header into SoEEDF src/ or include/, and never remove the per-target compile_definitions lines (the global add_compile_options alone would not survive a target added by another path).

**Cost:** The one property that makes SoEEDF embeddable and reusable is lost silently, with no test to notice.


### SoEEDF is embedded into SoPlasma by a LINK-TIME path variable, not a symlink. There is no symlink anywhere under soplasma-scratch pointing at SoEEDF. The coupling is `BOLTZMANN_DIR ?= $(HOME)/Projects/SoEEDF` in src/models/plasmaModels/plasmaBoltzmann/Make/options:10, with `-L$(BOLTZMANN_DIR)/build -lSoEEDF` at line 21 and an -Wl,-rpath in the generated linux64GccDPInt32Opt/options:51.
**Evidence:** find soplasma-scratch -maxdepth 6 -type l -lname '*SoEEDF*' returns nothing; grep -rn 'SoEEDF' in plasmaBoltzmann/Make/*; Make/options:8 comment explains 'the loader cannot find libSoEEDF.so at run time -- and the fix'.

**Rule:** After ANY change to SoEEDF's ABI or headers, rebuild SoEEDF's build/ FIRST, then rebuild SoPlasma — the rpath means a SoPlasma binary silently keeps loading whatever libSoEEDF.so is at $HOME/Projects/SoEEDF/build. Confirm with `ls -la $HOME/Projects/SoEEDF/build/libSoEEDF.so` and check its timestamp against the SoPlasma link.

**Cost:** The pre-run-checklist failure mode in its purest form: results attributed to new code produced by an old .so.


### SoGlobal does not exist and is explicitly recorded as not existing, with a check date. CLAUDE.md:6 states '(NOT YET IMPLEMENTED -- planned.) ... There is no SoGlobal tree yet: checked 2026-09-04', and framework-state.md:40 repeats '**There is no SoGlobal tree.**' What exists instead is plasmaChemistry0D, a utility inside the SoPlasma tree, documented in solver-0d.md and tutorial-0d-reactor.md; its stated gap is the self-consistent field (E/N as an eigenvalue) which is not implemented.
**Evidence:** grep -rni 'soglobal' across both CLAUDE.md and both docs trees returns exactly those two lines; framework-state.md §1.3.

**Rule:** Never scaffold or reference a SoGlobal directory. Any 0-D global work goes into soplasma-scratch's plasmaChemistry0D, and the eigenvalue-E/N gap is the item to cite when asked what is missing.

**Cost:** A third tree would fork the mechanism, the chemistry backends and the Boltzmann solve that plasmaChemistry0D currently shares with the CFD solver.


### The smoke beds are the documented worst case of a missed sweep, and the .gitignore now carries the post-mortem inline: five beds (smoke/{native,native-dt2,native-stiff,native-stiff-dt2,cantera}) were ignored wholesale, so 'no key-rename sweep ever reached them, and by 2026-09-02 all five carried FOUR settings the solver had come to reject and none of them could start'. Only run output is ignored now (614 MB), and the ignore patterns deliberately avoid `[0-9]*` because that glob 'has already destroyed 0.orig twice here'.
**Evidence:** soplasma-scratch/.gitignore lines 24-45; memory smoke-beds-lag-refactors.md lists the four in discovery order (limitVoltageRiseRate+maxVoltageRiseRate, dangling appliedVoltage, missing electronEnergyModel, unsplit chemistrySolver). git grep confirms all five beds still carry dated MIGRATED/SPLIT comments at system/plasmaSimulationControls:37-38 and constant/plasmaTransportProperties:39.

**Rule:** Every key sweep must include `smoke/` by path, and must read COMMENT text as well as key lines — the five beds' comments are the historical record and must be preserved (rule 19), not rewritten.

**Cost:** The smoke tests — the thing that is supposed to catch a broken build — were themselves broken and invisible for weeks.


### The claim-dating tool exists but is unusable synchronously on this machine: tools/audit_claims.py timed out at 110 s both repo-wide and scoped to docs/ alone. The cause is git blame cost — a single `git blame --line-porcelain -- docs/INDEX.md` takes 3.58 s real with ~0.01 s CPU (WSL2 I/O bound), so ~45 claim-bearing docs is ~160 s minimum.
**Evidence:** `timeout 110 python3 tools/audit_claims.py . --undated` EXIT=124; `timeout 110 python3 tools/audit_claims.py docs --undated` real 1m50, EXIT=124; `time git blame --line-porcelain -- docs/INDEX.md` real 0m3.581s. Its docstring records the baseline: 154 claims in src/ and 166 in docs/, of which 9 and 12 carried a date.

**Rule:** Always run audit_claims.py with run_in_background (rule 10 threshold is 1-2 minutes and this exceeds it at every scope), redirecting to a file; never inline in a turn. `git log -1 --format=%cs -- <file>` costs 0.002 s and answers 'is this document stale?' when per-line dating is not needed.

**Cost:** The tool gets abandoned as 'hung', and undated measurements keep being read as live conclusions — the exact 2026-08-30 failure that produced it.


### The comparison-contract tool lives in the WRONG tree for the cases it serves: tools/compare_cases.py and tools/audit_claims.py are in /home/kkourtza/Projects/SoEEDF/tools/, while every OpenFOAM case they analyse lives in /home/kkourtza/soplasma-scratch. soplasma-scratch/tools/ has no compare_cases.py (only compare_at_voltage.py).
**Evidence:** ls SoEEDF/tools -> audit_claims.py, co_sweep_compare.py, compare_cases.py, ...; ls soplasma-scratch/tools | grep -i compare -> compare_at_voltage.py only. compare_cases.py docstring: 'Compare OpenFOAM plasma cases against the baselines their COMPARE.md declares.'

**Rule:** Invoke it cross-tree with the absolute path: `python3 /home/kkourtza/Projects/SoEEDF/tools/compare_cases.py <caseDir>`. Record that path in CLAUDE.md/PROGRESS.md, because 'the tool does not exist here' is the failure that led to guessing baselines from directory names.

**Cost:** The 2026-08-30 wrong-control failure: every measured number correct, conclusion backwards by a factor of 9000.


### The memory base is larger than CLAUDE.md claims and the local pointer is stale. The directory holds 129 files (128 memories + MEMORY.md), not '~95'; and /home/kkourtza/.claude/projects/-home-kkourtza-soplasma-scratch/memory/MEMORY.md names `session-state-2026-09-06.md` as the newest status when session-state-2026-09-07.md, -09-07b, -09-07c and -09-08.md all exist.
**Evidence:** ls .../-home-kkourtza-Projects-SoEEDF/memory/ | wc -l -> 129; the soplasma-scratch MEMORY.md pointer text vs the actual session-state-* filenames.

**Rule:** Fix the two counts when writing CLAUDE.md/PROGRESS.md, and state the newest session-state by GLOB not by name: `ls -1 .../memory/session-state-*.md | sort | tail -1`.

**Cost:** A session reads the pointer, opens a five-day-old status file, and reconstructs the wrong state — precisely the failure the pointer was written to prevent.


### RULE POST-MORTEMS not in the consolidated rule text, compressed (rule 25 thread-switching, 26 probes, 27 invariants, 28 say-first are the well-known four; these are the rest). R29: a `break` in `while (pimple.loop())` skipped updateChargeDensity() on 364619 of 364670 steps (151 calls, 0.04%) leaving the field at the VACUUM value, then pimpleControl's counter accumulated across timesteps so 45-50% of ALL timesteps ran ZERO correctors and reported time ran 1.96-2.00x ahead of the physics; ten cases and a day of conclusions re-audited. R30: a custom relativeChange criterion was written because residualControl 'stalled' — it stalled only when gated on nEps_e (floor cells collapse the normalisation); gated on ePotential it had already measured ZERO rejections. R31: `if (pf.size())` around gAverage in plasmaExternalCircuit DEADLOCKED under `simple` decomposition and reported MPI_ERR_TRUNCATE under `scotch` — one cause, two symptoms; a sibling early-return in ddWallFluxMixed::updateCoeffs() produced a plausible PHYSICS number (-1e+300) instead of a crash and cost two hours in the wall-flux closure. R32: Boeuf & Pitchford 1995 was characterised from a secondary review — the paper actually solves ONE ambipolar equation and calls a self-consistent Te equation 'beyond the scope'. R35: a missing x10^-6 exponent on 3.7 in eq.14 was INFERRED from the NRL formulary instead of asking the user who had the PDF. R37: fvMatrix::residual() called standalone on a decomposed field misreports by ~21 ORDERS OF MAGNITUDE; the fix required a codebase-wide grep of every .residual() call site, not just the new JFNK one. R38: two days with no commit left 69 modified + 910 untracked paths and the ENTIRE new Newton solver (~6900 lines) untracked, while `git add tutorials/plasma` pulled in 436 paths of run output alongside 10 real dictionary changes. R39: PETSc printed 'Run with -malloc_debug' and that line was ignored through TWO wrong guesses with rebuilds; -malloc_debug then named the cause in ONE run (PCApply_Shell write past end of array — SNESVINEWTONRSLS hands the PC a SHORTER vector). R41: grubert2009_ballast400_picard reached the project's first genuine ignition (n_e 2.08e17 from a 1e11 seed) at t=1.97e-6 with writeInterval 5e-6 against endTime 20e-6 — ZERO time directories on disk; and hours after writing the rule, a 3.1 GB log led to KILLING a run (157 steps, 4.09e-9 s, ~2 h wall clock destroyed) when `: > the.log` would have solved it. R42: FOUR inert knobs in one day — -fieldsplit_0_pc_type hypre (splits are named 'phi'/'transport', so there is no split '0'), maxDeltaT written to system/controlDict (plasmaTimeControl reads its OWN dict, default GREAT, all three arms ran dt=1.29e-9 byte-identical), `sed s/outerSolver.*/newton/` on a case with no outerSolver key (two 'Newton' arms ran Picard for 24000 steps), and a 40000-step 'Newton success' where Newton was off for 21174 of them. R43: exact LU was made the DEFAULT Poisson-block PC justified as '2000 rows, nnz 9190, fill 3.1' — sparse LU is ~O(n^1.5) in 2-D; the real requirement was a FIXED linear operator, which a frozen-setup fixed AMG cycle count also satisfies at O(n) (GMRES/CG cannot serve: their Krylov polynomial depends on b; Richardson and Chebyshev can). R44: per-cell Newton scaling was queued on a claimed 8-decade species spread; measured in the regime where Newton actually runs it is 1.7e4 (FOUR decades), and implementation continued anyway across 46 usage sites in the core of a working solver.
**Evidence:** /home/kkourtza/Projects/SoEEDF/docs/rules-postmortems.md, 397 lines, read in full; rules 25-44 keyed by ORIGINAL rule number (consolidated 44 -> 19 on 2026-09-10, mapping in docs/rules-consolidation-proposal.md).

**Rule:** Put the NUMBERS, not the rule text, into CLAUDE.md/skill descriptions — '45-50% of timesteps ran zero correctors', '21 orders of magnitude', 'four inert knobs in one day', 'zero time directories at the only ignition'. The rule text is already consolidated; the evidence is what makes an agent actually run the check.

**Cost:** Consolidation 44->19 moved the evidence OUT of the governing file; an agent that reads only CLAUDE.md gets the instruction without the reason and skips it under time pressure.


### modelling-guidelines.md (2026-08-21, 393 lines) and framework-architecture.md (2026-09-04, 295 lines) are the two SoEEDF docs that do NOT duplicate CAPABILITIES.md — they are decision-recipes and contracts respectively, not inventories. modelling-guidelines §8 gives recipes by application (streamer, radical yield, ns-pulsed, electronegative, low pressure), §9 a pre-flight checklist, §10 known limits; framework-architecture gives §4 boundaries/contracts, §6 deliberate non-goals, §7 extension points. Their staleness risk is different: framework-architecture §8 was the source of one entry in framework-state.md's stale register ('No gas heating. ... This is the next development.' vs gasTemperatureEnergyModel existing).
**Evidence:** Heading outlines of both files; framework-state.md §5 row 1 names framework-architecture.md §8 specifically.

**Rule:** Keep both authoritative for their own subject and do NOT fold them into CAPABILITIES.md — INDEX.md's owner table already assigns 'how to choose settings for your problem' to modelling-guidelines.md. Audit ONLY their §8/§10 'limitations' sections against CAPABILITIES.md on each sweep: those sections are the ones that duplicate an inventory and therefore rot.

**Cost:** A limitations section is read as a capability statement; validation-summary.md and framework-architecture.md have already contradicted each other on gas heating within the same 'deliberately non-overlapping' document set.


### SKILL 1 — key-rename / API-change doc sweep (D1). The single highest-value procedure: it is mandated by rule 1, has a documented worst case (the smoke beds), and is currently re-derived by hand every time.
**Evidence:** rules-postmortems.md Rule 1; memory api-changes-need-doc-sweep.md (6-point coverage list); memory smoke-beds-lag-refactors.md; measured git grep timings above.

**Rule:** Skill reads FIRST: rules-postmortems.md Rule 1, memory api-changes-need-doc-sweep.md, memory smoke-beds-lag-refactors.md, .gitignore lines 24-45. Runs `git grep -nE '<old>|<new>' -- ':!validation'` in BOTH trees, then checkConfigReference.py, then the broken-link loop. NEVER `grep -r` from a tree root; NEVER delete a historical mention (rule 19 — mark SUPERSEDED in place with the date); NEVER skip smoke/ or tutorial comment text. Tier: SONNET, medium effort — the search is mechanical, but each hit needs the one judgement 'historical record or live instruction?'

**Cost:** Four rejected keys across five smoke beds, none of which could start, unnoticed for weeks.


### SKILL 2 — config-reference honesty check.
**Evidence:** docs/reference/README.md 'Keeping it honest' section; tools/checkConfigReference.py docstring (three real instances found 2026-09-01/02).

**Rule:** Reads FIRST: docs/reference/README.md, the checkConfigReference.py docstring (its 'THE CONVENTION IT PARSES' and 'CASE VARIABLES ARE NOT SOLVER OPTIONS' sections). Runs the three invocations in exact_commands and reports MISSING READER / STALE DEFAULT / DANGLING against the recorded baseline. NEVER treat --case UNRESOLVED as a failure (5 scheme variables are legitimately nested in compound values); NEVER add a key to the reference without grepping for its reader; NEVER edit source to satisfy the doc. Tier: HAIKU for the run and the diff-against-baseline; escalate to SONNET only when a new MISSING READER appears, since deciding whether to delete the doc line or wire the key is a design call.

**Cost:** `ePotentialControls { nonCoupledResidualControl }` was documented in two files and has never had a reader; a user set a tolerance and an iteration cap that did nothing.


### SKILL 3 — build verification gate (three checks, not one).
**Evidence:** memory pre-run-checklist.md §2: build-all.sh printed BUILD-COMPLETE with zero 'error:' lines while four components had FAIL and libplasmaTools.so did not exist on disk; build-all.sh also refuses to build while a solver is running, exiting 0 with NO BUILD-COMPLETE.

**Rule:** Reads FIRST: memory pre-run-checklist.md, memory rebuild-with-allwmake-not-piecemeal.md, memory openfoam-abi-partial-rebuild.md, memory allwmake-can-silently-skip-a-changed-directory.md. Asserts all three: `grep -c FAIL` == 0, `grep -c BUILD-COMPLETE` == 1, and `ls -la` on the changed .so showing a NEW timestamp. NEVER judge a build by the exit code of a grep pipeline (a clean build with unmatched output reports exit 1); NEVER skip killing running solvers first. Tier: HAIKU, low effort — three deterministic assertions with fixed pass criteria.

**Cost:** Symbols found in libraries came from an earlier manual wmake, not the build just 'run'; multiple wasted cycles.


### SKILL 4 — pre-launch case gate (Rule 36 + the 7-step pre-run checklist).
**Evidence:** rules-postmortems.md Rule 36 ('STATE THE FULL SET OF INITIAL AND MAIN PARAMETERS FIRST ... every relaunch, not only the first'); memory pre-run-checklist.md steps 1-7, especially step 6 (tables REGENERATED at run time from boltzmann/ENmax, so a hand-extended table never reached a single run) and step 7 (`grep -c 'beyond its range'` returned 0 all night; the real string is 'beyond the tabulated range').

**Rule:** Reads FIRST: rules-postmortems.md Rule 36, memory pre-run-checklist.md, memory config-variable-may-be-dangling.md, memory state-parameters-before-every-simulation-launch.md. Emits a parameter block (circuit type+values, applied V or I and its time profile, simulationType, gas/pressure/temperature, endTime, timestep/Courant strategy) BEFORE the launch command, and proves every edited config variable is referenced with `grep -rF '$'"$k"` (note: `"\$$k"` expands $$ to the shell PID and falsely reports ~47 dangling per case). NEVER launch on a parameter read from memory instead of from the file the case will actually use. Tier: SONNET, medium — the checklist is mechanical but step 7 ('prove the check CAN fail') is genuine reasoning about the instrument.


### SKILL 5 — knob-actually-moved proof (Rule 42).
**Evidence:** rules-postmortems.md Rule 42: four inert knobs in one day, each reported as a result before it was checked.

**Rule:** Reads FIRST: rules-postmortems.md Rule 42. Before ANY comparison is believed, runs the matching cheap check from the run's own output: `grep -c SNES <log>` before a Newton claim, `-ksp_view` before a preconditioner claim, the ACHIEVED dt before a timestep claim, `-options_left` to catch PETSc options nothing consumed, and `cmp` on the two arms' output before reporting a difference. Treats BYTE-IDENTICAL output across arms as CONCLUSIVE proof the knob is unwired — go read the code, do not re-run. NEVER report an arm comparison without the achieved-value line beside the result. Tier: HAIKU, low — four fixed greps with fixed pass criteria.


### SKILL 6 — snapshot / stop / restart protocol (Rule 41).
**Evidence:** rules-postmortems.md Rule 41: ignition at t=1.97e-6 with first snapshot due at 5e-6 -> zero restart points; a restart REWRITES the snapshot it restarts from including uniform/time deltaT (observed first steps 7.71e-12 vs 1.11e-11, each exactly twice its own stored deltaT) while the FIELDS stay byte-identical; and killing a run over a 3.1 GB log destroyed 157 steps and ~2 h.

**Rule:** Reads FIRST: rules-postmortems.md Rule 41, memory restart-from-snapshots-and-preserve-them.md, memory kill-runs-before-launching.md. Enforces writeInterval <= endTime/10 AT CASE CONFIGURATION time; before any intentional stop runs the three-step sequence (set `stopAt writeNow;` then `touch system/controlDict` because the configuration/ #include is only reprocessed when controlDict's own mtime changes, then CONFIRM the new time directory name matches the time the log last reported, THEN stop); copies a reused restart point aside read-only first. NEVER kill a run to fix a log — `: > the.log` truncates in place. NEVER restart for a claim that depends on the trajectory's history (temporal order, dt sequences, how a state was reached). Tier: SONNET, medium — the sequence is mechanical, the 'can a restart change this conclusion?' proviso is the judgement and must be stated explicitly before launch.


### SKILL 7 — stale-claim and design-status audit.
**Evidence:** tools/audit_claims.py (154 undated claims in src/, 166 in docs/; measured 3.58 s per git blame here); framework-state.md §5 register; 14 of 22 design docs with no Status line.

**Rule:** Reads FIRST: tools/audit_claims.py docstring, framework-state.md §5, rules 18-20 via memory dated-measurements-and-supersession.md. Runs the design-Status loop (instant), the `git log -1 --format=%cs` per-document freshness table (0.002 s each), and audit_claims.py IN THE BACKGROUND only. NEVER delete a superseded claim — mark it SUPERSEDED BY <what> in place with the date. NEVER treat age as wrongness: the output is a manifest, not a verdict. Tier: SONNET, medium — deciding which of two contradictory documents is true requires reading the source, not the documents.


### SKILL 8 — commit staging hygiene (Rule 38).
**Evidence:** rules-postmortems.md Rule 38: 69 modified + 910 untracked with the whole ~6900-line Newton solver untracked; `git add tutorials/plasma` staged 436 paths of run output for 10 real dictionary changes.

**Rule:** Reads FIRST: rules-postmortems.md Rule 38. Reports uncommitted work PLAINLY and as its own statement after any new file/dir, after a logically complete unit, and before anything long-running. Stages EXPLICITLY by path and re-reads `git diff --cached --name-only`, unstaging what does not belong: src/, docs/, tools/, etc/, build-all.sh and tutorial DICTIONARIES are the work; validation/ and tutorial time directories are run output and stay out. NEVER `git add -A`; NEVER `git add <dir>` on a tree holding case output; NEVER commit unasked — this is a PROMPT. Tier: HAIKU, low — a fixed include/exclude list plus a diff read-back.


### SKILL 9 — comparison contract (COMPARE.md + compare_cases.py).
**Evidence:** tools/compare_cases.py docstring and memory comparison-contract-on-disk.md: the 2026-08-30 wrong-control report, backwards by a factor of 9000 (1.564e18 -> 1.69e14 was a success, reported as a failure).

**Rule:** Reads FIRST: the compare_cases.py docstring (it contains the COMPARE.md format verbatim — a fenced ```compare block with question/baseline/varies/matches/time/field/region/reference), memory comparison-contract-on-disk.md. Writes COMPARE.md into the case dir WHEN THE CASE IS CREATED, with ABSOLUTE baseline paths. Invokes the tool by absolute cross-tree path. NEVER reconstruct a comparison set from directory names; NEVER report a conclusion in the same message as a prior number that failed to reproduce — say 'I cannot reproduce X' and STOP. Tier: SONNET, medium — the tool is mechanical, but 'what makes this row a valid control' is the judgement the tool prints and cannot make.


### SKILL 10 — case cold start from the layer-1 generators.
**Evidence:** memory reread-inventory-after-context-rebuild.md: a cold start was hand-rolled with nohup, two `exit 127`s and a nonexistent utility name (`plasmaCreateFields`; the real one is `plasmaCreateSpeciesFields`), with maxDeltaT as a literal in system/plasmaSimulationControls, while tools/make_mesh_arm.sh already encoded the mesh side including three measured traps (no `set -e`; `export` on separate lines; `rm -rf 0 && mkdir -p 0`).

**Rule:** Reads FIRST: docs/CAPABILITIES.md (the tooling inventory), memory two-layer-case-architecture.md (G2), memory boundary-role-library-layer1.md, tools/make_mesh_arm.sh, tools/plasmaSetupRegions.sh. Drives the real chain: mesh via make_mesh_arm.sh -> plasmaSetupRegions.sh -> configuration/boundaries + etc/boundaryRoles -> foamPlasmaSetupBoundaries -> foamPlasmaCreateSpeciesFields. NEVER hand-write a layer-2 OpenFOAM dictionary (G2: generated, header says DO NOT EDIT, gitignored); NEVER invent a utility name — the nine utilities are listed in framework-state.md §1.4 and CAPABILITIES.md. Tier: SONNET, low-medium — it is an encoded sequence; the only judgement is the `advanced { }` passthrough when an option is unmodelled.


### SUBAGENT 1 — parallel-correctness reviewer (collectives and decomposition).
**Evidence:** rules-postmortems.md Rule 31: two spellings, two full debugging sessions, one DEADLOCK and one plausible-but-fake physics number (-1e+300); 'reading the collectives got nowhere three times'.

**Rule:** Must read FIRST: Rule 31 in rules-postmortems.md, memory collective-behind-a-local-guard.md, memory parallel-models-guard-divergence.md, memory fvmatrix-residual-broken-in-parallel-use-fvc-instead.md. Reviews every reduce/gSum/gMax/gMin/gAverage/returnReduce for a rank-dependent guard or early return upstream, and every standalone .residual() call. MUST NEVER recommend guarding a collective harder — the fix is to HOIST IT OUT (gAverage over an empty local field is well defined) or DELETE it (FatalErrorInFunction aborts the whole job anyway and names the owning rank). MUST NEVER approve on inspection alone: it must require per-rank Pout probes bracketing each construction and a test under TWO decompositions (simple AND scotch — they failed differently on the same bug and either alone misleads). Tier: OPUS, high effort — the failure produces a plausible physical diagnosis rather than a crash, which is exactly the class where cheap models confidently agree with the wrong hypothesis.


### SUBAGENT 2 — loop-contract and per-step invariant auditor.
**Evidence:** rules-postmortems.md Rule 29 (the worst defect of the project: 151 updateChargeDensity calls in 364670 steps, then 45-50% of timesteps with ZERO correctors and the clock 1.96-2.00x ahead of the physics) and Rule 27 (invariants before physics: steps (a)-(c) together cost under a minute against a day lost).

**Rule:** Must read FIRST: Rules 27 and 29, memory raw-break-skips-loop-contract.md, memory step-discard-invariant.md, and the plasmaStepAudit source. Reviews any change to loop control flow and demands per-step invariant COUNTS (a call counter, not inspection) plus registration of any new per-step update with plasmaStepAudit. MUST NEVER accept 'the case runs' or 'it looks converged' as evidence; MUST NEVER approve a break/return/goto out of a control object's loop — the fix is a flag plus DRAINING the loop. Tier: OPUS, high — it reasons about hidden object state across iterations, and the cost of a miss is every conclusion drawn on the affected runs.


### SUBAGENT 3 — wiring auditor ('declared but unreachable').
**Evidence:** memory declared-but-unreachable.md: three defects of one family in one session — a fully implemented plasmaEnergy::eEqn() invoked from nowhere; argList::addOption("lmeaDt") never fetched, so an 18-run dt scan came back BIT-IDENTICAL; and nEfloor_=1.0 guarding a quantity whose species floor is 1e13, unreachable by thirteen orders of magnitude while its diagnostic printed the OPPOSITE of the truth.

**Rule:** Must read FIRST: memory declared-but-unreachable.md, memory single-source-of-truth-defaults.md, memory verify-wiring-not-just-compilation.md. For every new function greps for a CALLER outside its own definition; for every new key greps for the READER (getOrDefault/lookup/args.found), not the declaration; for every new guard/floor/cap asks which OTHER limit is already active on that quantity and whether this one can ever bind. Also flags any key read with getOrDefault in more than one component (the half-on defect: a single relaxed field took 191 correctors against 91 for the joint pair). MUST NEVER accept 'it compiles' or a passing run as wiring evidence. Tier: SONNET, high effort — the greps are mechanical but the binding-threshold comparison and the two-owner detection need real reasoning; the volume of call sites is what needs the effort, not the depth.


### SUBAGENT 4 — documentation-sweep auditor (the D1 gate).
**Evidence:** Rule 1; the measured surface above; 45+37 docs, 14/22 design docs with no Status line, 1 broken reference link, 7 docs orphaned from INDEX.md, INDEX.md blind to three whole SoPlasma doc subtrees.

**Rule:** Must read FIRST: Rule 1, memory api-changes-need-doc-sweep.md, docs/CAPABILITIES.md, docs/INDEX.md. Runs AFTER a change is otherwise complete and reports what the change did NOT touch: the git grep hits in both trees, checkConfigReference output, broken links, orphaned docs, missing Status lines, and whether CAPABILITIES.md was updated in the same commit. MUST NEVER write the documentation itself (that is the implementer's job and the reason the sweep is separable); MUST NEVER pass a change that touched a dictionary key without a corresponding docs/reference/ edit. Tier: SONNET, medium — enumeration plus a per-hit historical-vs-live call.


### SUBAGENT 5 — numerics-scaling reviewer (Rule 43).
**Evidence:** Rule 43: exact LU made the default Poisson-block PC on a 2000-cell justification; the naming of the real reasoning error (fixedness is a property of the ITERATION, not of exactness); GMRES/CG cannot serve as a fixed linear operator because their Krylov polynomial depends on b, Richardson and Chebyshev can. User: 'if this degrades badly in large domains with millions of cells then whats the point????'

**Rule:** Must read FIRST: Rule 43, docs/design/newton-outer-solver-design.md, docs/design/schur-semiimplicit-poisson-preconditioner.md (Status: DERIVED, not yet implemented), memory ksp-cap-was-the-bug.md, memory newton-vs-picard-benchmark-state.md. Requires a 2k -> 20k -> 200k sweep reporting KRYLOV ITERATIONS PER NEWTON STEP at each size, not wall-clock at one size. MUST NEVER endorse a default validated at a single problem size; MUST NEVER endorse an automatic size-dependent algorithm switch (it makes the solver untestable — a benchmark below the switch says nothing above it); the unscalable-but-exact option stays available as explicit opt-in only. Tier: OPUS, high — this is where a wrong endorsement becomes a default that silently makes the target 2-D/3-D cases unaffordable.


### SUBAGENT 6 — physics and literature reviewer (Rules 32, 34, 35).
**Evidence:** Rule 32 (Boeuf & Pitchford 1995 mischaracterised from a secondary review, then a novelty assessment begun on WebSearch blurbs for three more unopened papers); Rule 34 (JC-PIC/Boeuf is the standing reference, in Literature/reference-codes-and-manuals/JC-PIC_Boeuf/, and one reading surfaced three results at once including the absence of any electron-electron collision term in SoPLASMA); Rule 35 (a missing x10^-6 inferred rather than asked).

**Rule:** Must read FIRST: Rules 32/34/35, Literature/reference-codes-and-manuals/JC-PIC_Boeuf/README.md table of contents, memory jcpic-boeuf-standing-reference.md, memory no-electron-electron-collisions-in-soplasma.md. Checks JC-PIC's TOC BEFORE deriving a formula, choosing a validation target, or asserting no benchmark exists. MUST NEVER build an assessment, novelty claim or recommendation on an abstract, a search summary, or another paper's citation — if the full text is inaccessible it says so and asks for the PDF; MUST NEVER guess an ambiguous in-paper detail (equation, coefficient, exponent) when the user holds the PDF — it asks. Tier: OPUS, high — the failure mode is a confident, well-written, wrong physics claim, which is the single most expensive thing a cheap model produces here.


### SUBAGENT 7 — evidence and claims auditor (Rules 11-20, 42, 44).
**Evidence:** memory comparison-contract-on-disk.md ('Reminders do not fire; artifacts do'); memory never-read-a-number-without-its-control.md; Rule 44 (an 8-decade premise measured at 1.7e4 — four decades — and the work continued anyway across 46 usage sites in the core of a working solver).

**Rule:** Must read FIRST: docs/working-with-claude.md §2.6, memory comparison-contract-on-disk.md, memory dated-measurements-and-supersession.md, memory never-read-a-number-without-its-control.md, Rules 42 and 44. For every reported result checks: is the control named and valid (varies/matches), is the achieved value printed beside the requested one, is the measurement DATED, and has the expected payoff changed since the work was agreed. MUST NEVER let a conclusion share a message with an unreproduced prior number; MUST NEVER contradict an established project conclusion from recall — it re-reads source and docs and reconciles explicitly. Tier: OPUS, medium-to-high — pure judgement, but on a short, well-specified checklist.


### SUBAGENT 8 — case-configuration reviewer (layer 1 / G2).
**Evidence:** memory two-layer-case-architecture.md (G2, commit 782c305, 2026-09-04: the kind vocabulary, why drivenElectrode and groundedElectrode are separate kinds, why regionInterface is never declared, why material is a VOLUME property); memory config-variable-may-be-dangling.md (appliedVoltage dangling in five beds while a hardcoded table drove the electrode; a 'low-field' arm ran at the original 62.1 Td); checkConfigReference --case currently reports simulationType DANGLING in needleDBD.

**Rule:** Must read FIRST: memory two-layer-case-architecture.md, memory boundary-role-library-layer1.md, etc/boundaryRoles, docs/reference/configuration-config.md, and runs `checkConfigReference.py --case <dir>`. Checks every configuration/config variable resolves to a key src/ actually reads, that no rejected key is present, and that boundary `kind` values come from the agreed vocabulary. MUST NEVER approve a hand-edited layer-2 dictionary (G2); MUST NEVER put a value used in exactly one dictionary into configuration/config (measured on needleDBD: 22 of 37 variables were used once; config went 37 -> 15). Tier: SONNET, medium — a fixed rule set plus a tool, with one judgement (does this value belong in config or in the dictionary that reads it).


## Traps
- An index that claims completeness but is not: docs/INDEX.md opens with 'Every document in the framework' and 'each fact has one home ... if you find the same explanation in two places, one of them is a bug'. It omits 7 of its own docs and all of SoPlasma's design/, reference/ and CAPABILITIES.md. Detector: the two one-line loops in exact_commands (orphan check + broken-link check), <0.1 s each.
- Two capability inventories drifting: CAPABILITIES.md (2026-09-11) knows about the Newton/SNES solver 88 times; framework-state.md (2026-09-05) and capability-overview.md (2026-08-21) mention it zero times. Detector: `grep -cin 'newton\|snes\|jfnk'` across all inventory docs — any inventory returning 0 while another returns 88 is stale by a whole subsystem.
- A reference doc that the checker silently does not check. checkConfigReference.py only sees `### `key`` headings; plasmaSimulationControls.md (10113 B) and fvSchemes-fvSolution.md (7195 B) have ZERO, so they contribute nothing to 'documented' and their claims are unverified. Detector: `grep -cE '^###[[:space:]]+.*`[A-Za-z_][A-Za-z0-9_]*`' docs/reference/*.md` — a large file with count 0 is invisible to the gate.
- configuration-config.md is DELIBERATELY SKIPPED by the default checker run ('note: skipping configuration-config.md -- it documents case VARIABLES, not option keys'). Reading 'OK' after that note does NOT mean config is checked. Detector: the note itself appears on stdout; the only real check is `--case <caseDir>`.
- A rename sweep that appears to finish because grep was killed. `grep -r` over soplasma-scratch exits 143 at 120 s with partial output and no error banner; a wrapper that pipes to head or `|| true` turns that into a clean-looking empty result. Detector: always use git grep, and check ${PIPESTATUS[0]} not the pipeline's status.
- A doc marked 'plan / not implemented' long after the work landed, and the inverse. framework-state.md §5 caught four such plan documents at once. Detector: the design-Status loop plus `git log -1 --format=%cs` per file; 14 of 22 SoPlasma design notes currently have no Status line at all, so the CAPABILITIES.md instruction 'check the Status: line' silently returns nothing.
- checkConfigReference's MISSING READER can be a FALSE POSITIVE from its own hand-maintained allowlists. `uniformValue` is flagged today because it is an OpenFOAM-owned changeDictionary key absent from OPENFOAM_OWNED; NOT_USER_FACING only suppresses UNDOCUMENTED, not MISSING READER. Detector: before acting on a MISSING READER, grep src/ for the key yourself and check OPENFOAM_OWNED/COMPUTED_KEY_READS/NOT_USER_FACING near the top of the script.
- audit_claims.py reads as hung rather than slow (>110 s at every scope, ~0 CPU, pure git-blame I/O on WSL2). Detector: `time git blame --line-porcelain -- docs/INDEX.md` = 3.58 s; multiply by the file count before deciding whether to foreground it. Run it in the background.
- Documentation for a tool that lives in the other tree. compare_cases.py and audit_claims.py are in SoEEDF/tools/ while the cases they serve are in soplasma-scratch. A session that looks in soplasma-scratch/tools/ concludes the tool does not exist and reconstructs a comparison set by guessing directory names — the exact 2026-08-30 failure. Detector: `ls /home/kkourtza/Projects/SoEEDF/tools/` before concluding any tool is missing.
- A zero-byte documentation file that looks like coverage in a listing: docs/models/transport/drift_diffusion/mobilityModels.md and diffusivityModels.md are both 0 bytes. Detector: `find docs -name '*.md' -size -1c`.
- A historical mention deleted instead of marked. Rule 19 requires SUPERSEDED text be marked in place; the three fixed entries (AMR.md:13 CORRECTED, ddSolidSurfaceFlux.md SUPERSEDED header, validation-summary.md:419 parenthetical) are the pattern. Detector: after a sweep, `git diff` must show ADDED supersession markers, not removed lines, in files whose content is historical record — including the dated MIGRATED/SPLIT comments in all five smoke beds.

## Open questions
- Is docs/rules-consolidation-proposal.md (9596 B, the 44 -> 19 mapping) the current mapping, or was it superseded by what actually landed in CLAUDE.md on 2026-09-10? It is orphaned from INDEX.md, so nothing states its status. The orchestrator has read CLAUDE.md and can settle whether the consolidated numbering matches the proposal.
- Does a CI or pre-commit hook exist anywhere that runs checkConfigReference.py? The tool's own docstring says 'in CI, if there is one'. I found no .github/, no hook, and only settings.local.json under both .claude/ directories — but I did not enumerate git hooks (.git/hooks) in either tree.
- What are the current true values behind the two checkConfigReference defects? Specifically: is `uniformValue` genuinely an OpenFOAM-owned key that belongs in OPENFOAM_OWNED (making the MISSING READER a tool bug), and is `defaultSEEC`'s real default the '<derived; see the reader>' path at electronDDWallFluxImplicitFvPatchScalarField.C:224 (making the doc's 0.001 simply wrong)? Both need a source read I did not do.
- Are the two *pending* reference documents (plasmaSpeciesProperties, plasmaTransportProperties) queued anywhere, and was regions-and-materials.md ever written and lost, or only ever planned? docs/reference/README.md links it as if it exists.
- How many of framework-state.md §5's 11 stale-document rows are still open? I verified 3 have been fixed in place (AMR.md, ddSolidSurfaceFlux.md, validation-summary.md gas heating) and 1 is still live (INDEX.md does not carry a phase-by-phase status for gas-heating-plan.md). The other 7 I did not check individually.
- Does SoEEDF have any test that would FAIL if an OpenFOAM header were introduced, beyond the CMake build itself? I found BOLTZMANN_NO_OPENFOAM only as compile definitions and #ifndef guards; there is no grep-based or CI guard, so the standalone constraint is enforced only by whoever runs cmake.
- Does plasmaStepAudit have a registration list I could enumerate to check that every per-step update is actually registered (Rule 29's enforcement)? I located the rule text but did not read the class.
- What is the actual current state of the deferred wall-loss revalidation (framework-state.md item #24)? It is described as blocking the citation of any streamer or anode-layer number in any document, which would make several documents unquotable — but that claim is itself a week old.
