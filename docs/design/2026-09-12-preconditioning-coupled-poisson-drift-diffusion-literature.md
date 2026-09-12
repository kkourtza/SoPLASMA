# Preconditioning the fully-coupled Poisson + drift-diffusion system: what the literature already says

**Status: LITERATURE NOTE, 2026-09-12 (Part A) + 2026-09-12 second pass (Part B).** No code, no
numbers handed to the solver. **PART B, at the end of this file, SUPERSEDES parts of Part A** —
see its supersession list before quoting §0.3, §1.3, §3.2 or §4.2. This note
answers five questions put by the parent agent and records, for each claim, (a) whether it is
ESTABLISHED in a source I read in full, (b) OUR OWN measurement, or (c) MY INFERENCE. Nothing
here is concluded from an abstract, a search snippet, or another paper's citation of a work
(rule A4). Where I could not obtain a full text I say so and ask for the PDF.

Superseding/related notes, NOT edited by this one:
`docs/design/schur-semiimplicit-poisson-preconditioner.md` (our derivation, 2026-09-10),
`docs/design/newton-outer-solver-design.md`, `docs/design/newton-ignition-experiments.md`,
memory `ap-proof-newton-removes-dielectric-constraint`.

---

## 0. Evidence inventory — what I actually read, and what I could not get

### 0.1 Read IN FULL (or the complete section carrying the claim, with its stated assumptions)

| source | local path | what I read |
|---|---|---|
| Hagelaar, HDR *Modelling methods for low-temperature plasmas* | `SoEEDF/Literature/plasma-fluid-closures-LMEA-LFA/hdr-hagelaar.pdf` | Ch. 8 pp. 48-50 in full, **eqs. (8.3)-(8.16) checked on the RENDERED pages 49-50**, plus bibliography entries |
| Villa, Barbieri, Gondola & Malgesini, *JCP* **242** (2013) 86-102 | `SoEEDF/Literature/Villa_1-s2.0-S0021999113001265-main.pdf` | §2-3 (scheme + AP proof), §5 (algorithm), §6.2-6.3 (cost); **p. 99 checked on the RENDERED page** |
| Knoll & Keyes, *JCP* **193** (2004) 357-397 | `SoEEDF/Literature/Knoll_1-s2.0-S0021999103004340-main.pdf` | §3.3-3.6 in full (multigrid-, physics-based, structure-based, matrix-free, nonlinear preconditioning) + bibliography |
| Duarte Gomez, Deak & Bisetti, *JCP* **480** (2023) 112007 | `SoEEDF/Literature/preconditioning-and-coupled-solvers/DuarteGomez_Bisetti_2023_*.pdf` | abstract, §1, §5 (preconditioner) in full, §7.1-7.2.4, §7.3.3, §8, §9. Accepted manuscript, 56 pp, DOE/Elsevier user licence |
| Teunissen, *PSST* **29** (2020) 015010 | `SoEEDF/Literature/flux-schemes-and-numerics/Teunissen_2020_*.pdf` | §3 in full; **article pp. 5-6 checked on RENDERED pages** (also cross-checked against arXiv:1911.01190) |
| Chacón, *JCP* **526** (2025) 113789, preprint arXiv:2407.07031 | `.../preconditioning-and-coupled-solvers/Chacon_2025_*.pdf` | abstract, §1, §3-3.2 in full (parabolization ↔ Schur factorisation), §Conclusions |
| Manteuffel, Ruge & Southworth, *SIAM J. Sci. Comput.* **40** (2018) A4105-A4130, DOI 10.1137/17M1144350, preprint arXiv:1708.06065v4 | `.../preconditioning-and-coupled-solvers/Manteuffel_*.pdf` | abstract, §1 in full, §4.2.1-4.2.2, §5 Conclusions |
| Arslanbekov & Kolobov, arXiv:2003.03812 (2020) | `.../preconditioning-and-coupled-solvers/Arslanbekov_Kolobov_2020_*.pdf` | abstract, §3.3 (linear solver) and §3.4 (nonlinear convergence and time stepping) in full |
| Zhu, Wu, Cao, Wei & Xu, arXiv:2607.18123 (2026) | `.../preconditioning-and-coupled-solvers/Zhu_2026_*.pdf` | abstract, §1 in full, §3 efficiency passage |
| Pasolari & Kourtzanidis, arXiv:2607.05137 (our own) | `SoEEDF/Literature/SoPlasma-development/*.pdf` | §2.1.1, §6.3.1 |
| hypre user + reference manual, v3.2.0 | fetched | §5.7.6 "AMG for systems of PDEs"; `HYPRE_BoomerAMGSetNumFunctions` / `SetDofFunc` / `SetNodal` / `SetFilterFunctions` reference entries |
| PETSc 3.25 manual pages | fetched | `PCGAMG`, `MatSetBlockSize` Notes sections |

JC-PIC standing reference (`Literature/reference-codes-and-manuals/JC-PIC_Boeuf/`) was consulted
as required: `grep -i` over `library_article_book.txt` and `manual.txt` for *dielectric
relaxation*, *implicit field*, *semi-implicit*, *preconditioner*, *Newton* returns **nothing on
this topic** — it is an explicit PIC-MCC code and its book does not treat implicit field solves or
preconditioning. It contributes nothing here, which is itself worth recording so nobody re-checks.

### 0.2 OBTAINED but NOT read in full — do not build on these until read

* Chen, Chacón, Leibs, Knoll & Taitano, *JCP* **258** (2014) 555, preprint arXiv:1309.6243,
  `.../preconditioning-and-coupled-solvers/Chen_Chacon_2014_*.pdf` — abstract + targeted grep only.
* Brannick, Falgout, Kahl, Schroder & Shen, arXiv:2501.16248 (2025), nodal AMG coarsening for PDE
  systems — abstract + §1 only. It targets H(curl)/H(div)/Stokes, not drift-diffusion.
* Lapenta & Ju, arXiv:0804.0150 (2008), predictor-corrector preconditioners (semi-implicit scheme
  as predictor, fully implicit as corrector) — abstract only. **This one looks directly on-point
  for us and should be read next.**

### 0.3 COULD NOT OBTAIN — I am asking for these PDFs

I cannot get the full text of the following. Each is a place where the specific thing we want to
claim as new could already exist, so **until I can read them I will not assess novelty against
them, and I report only that they exist**:

1. **Bank, Chan, Coughran & Smith, "The alternate-block-factorization procedure for systems of
   partial differential equations", *BIT Numer. Math.* **29** (1989) 938-954**, DOI
   10.1007/BF01932753. Cited by Knoll & Keyes [10] as a structure-based split preconditioner for
   *charge transport in semiconductors*. **This is the single most likely prior home for our
   Schur-on-φ identity** and I could not read it.
2. **Kerkhoven & Saad, "On acceleration methods for coupled nonlinear elliptic systems",
   *Numer. Math.* **60** (1991/92) 525-548**, DOI 10.1007/BF01385735. Knoll & Keyes §3.4.2 describe
   it as the drift-diffusion case where "the outer Newton-Krylov method is regarded as the
   accelerator to the inner fixed point method" — i.e. **Gummel-as-preconditioner**, which is
   exactly our Picard-vs-Newton question.
3. **Kerkhoven, "A proof of convergence of Gummel's algorithm for realistic device geometries",
   *SIAM J. Numer. Anal.* **23** (1986) 1121-1137**; **Kerkhoven, "On the effectiveness of Gummel's
   method", *SIAM J. Sci. Stat. Comput.* **9** (1988) 48-60**; and a "spectral analysis of the
   decoupling algorithm for semiconductor simulation" in *SIAM J. Numer. Anal.* (DOI 10.1137/0725073).
   These are the published convergence theory of the Gummel map. **Q4 cannot be answered properly
   without at least one of them.**
4. **Chacón & Knoll, "A 2D high-β Hall MHD implicit nonlinear solver", *JCP* **188** (2003) 573-592**,
   DOI 10.1016/S0021-9991(03)00193-1, and **Chacón, *Phys. Plasmas* **15** (2008) 056103**. These are
   the two papers Chacón's 2025 review cites for the parabolization↔Schur connection.
5. **Lin & Shadid, "Performance of a parallel algebraic multilevel preconditioner for stabilized
   finite element semiconductor device modeling", *JCP* **228** (2009) 6250-6267**, DOI
   10.1016/j.jcp.2009.05.024; **Lin, Shadid et al., *IJNME* (2012), DOI 10.1002/nme.4315**; and
   **Lin, Shadid & Hoekstra (2008), "Performance of various algebraic multigrid based preconditioners
   for the drift-diffusion equations for semiconductor devices"**. These are the fully-coupled
   Newton-Krylov-AMG drift-diffusion papers and are the direct answer to Q3.
6. **Meza & Tuminaro, "A multigrid preconditioner for the semiconductor equations",
   *SIAM J. Sci. Comput.* **17** (1996) 118-132.**
7. **Ventzek, Sommerer, Hoekstra & Kushner, *Appl. Phys. Lett.* **63** (1993) 605-607**, and
   **Boeuf, "A two-dimensional model of dc glow discharges", *J. Appl. Phys.* **63** (1988) 1342-1349** —
   the two works cited as the ORIGIN of the semi-implicit Poisson. I have read neither; see §1.3.

**Ask: can you pull 1, 2, 3, 5 through institutional access?** 1 and 5 decide the novelty question;
2 and 3 decide Q4.

---

## 1. Q1 — Is `S_f = div((eps + dt*sigma) grad .)` known, and under what name?

### 1.1 The OPERATOR is completely standard in plasma fluid modelling. It is the semi-implicit Poisson equation. Do not claim it.

**ESTABLISHED.** Hagelaar's HDR, Ch. 8 "Semi-implicit method", p. 49 (equations checked on the
rendered page, not on extracted text):

> (8.7)  `-div(eps0 grad Phi^{k+1}) = sum_j [ q_j n_j^k - dt div( |q_j| mu_j n_j^k grad Phi^{k+1} - q_j D grad n^k ) ]`
> (8.8)  `-div( eps0 (1 + chi_e) grad Phi^{k+1} ) = sum_j q_j (2 n_j^k - n_j^{k-1}) - div( eps0 chi_e grad Phi^k )`
> (8.9)  `chi_e = (dt/eps0) sum_j |q_j| mu_j n_j^k = dt / tau_d`

with, at (8.5),

> `tau_d = eps0 / sum_j |q_j| mu_j n_j  ~=  eps0 / (e mu_e n_e)`,

and the explicit statement that "the numerical discretisation must respect ... `dt < tau_d` for
collisional drift-diffusion models", which the semi-implicit form removes: "This scheme remains
stable even if `dt >> tau_d`."

`eps0 (1 + chi_e) = eps0 + dt*sigma`. **So the operator IS `div((eps + dt sigma) grad .)` and its
coefficient IS literally `dt/tau`, written as `chi_e` in the standard reference.** The framing
"the ratio `dt/tau` is the relevant parameter" is therefore not ours — it is the coefficient of
the published operator.

Hagelaar attributes the method to [Ven93] = Ventzek, Sommerer, Hoekstra & Kushner, *APL* **63**
(1993) 605-607, and notes the analogous devices for collisionless fluid models (Crispel, Degond &
Vignal, *JCP* **223** (2007) 208) and for PIC (Langdon, Cohen & Friedman, *JCP* **51** (1983) 107).

Villa et al. 2013, eq. (9), give the same operator in the form
`div[ (1 + dt^k sigma^k/eps0) E^{k+1} ] = q^k/eps0` with `sigma^k = e(n_e mu_e + n_p mu_p + n_n mu_n)`
(their eq. 10), and state in the text immediately after:

> "The term `sigma^k/eps0` is the inverse of the characteristic relaxation time ... This parameter,
> as we will see, plays a crucial role. When the time integration period `Dt^k` is much smaller than
> the relaxation time, Eq. (9) reduces to the explicit method. On the contrary, for large time steps,
> Eq. (9) reduces to the current continuity equation `div J^{k+1} = 0`."

They also note the linear-algebra consequence explicitly: "The conditioning number of the discrete
versions of (9) depends on the diffusion coefficient and, in turn, on the variations of the electric
conductivity, while in (7) [plain Poisson] the conditioning number is constant."

Hagelaar makes the same limit statement in words (p. 50): "for larger time steps `chi_e >> 1`, the
modified Poisson equation (8.8) becomes equivalent to the quasi-neutral current conservation
equation (8.16)", where his (8.16) is `div(sigma grad Phi_d) = 0`.

**Our own solver paper already uses this operator and already cites Villa for it**: Pasolari &
Kourtzanidis, arXiv:2607.05137, eq. (2), `div[(eps + dt sigma^k) grad phi^{k+1}] = -rho^k - dt sum_j q_j div(D_j grad n_j^k)`,
"following the approach proposed in [49]" = Villa et al. So there is no room at all to claim the
operator.

### 1.2 The IDENTIFICATION of a semi-implicit operator with a Schur complement is also standard, and it has a name: *parabolization* / *physics-based (approximate block factorization) preconditioning*.

**ESTABLISHED, but demonstrated in MHD and shallow water, not in plasma drift-diffusion.**

Knoll & Keyes 2004 §3.4.1 works the canonical example (1D shallow water with a stiff gravity wave).
They eliminate the momentum update from the continuity update in delta form and obtain, their eq. (41),

> `dh/Dt + d/dx( -Dt g h^n d(dh)/dx ) = res_h + d/dx(Dt res_uh)`

— a scalar parabolic operator `I/Dt - div(Dt * c^2 grad .)` in exactly the shape of ours, to be
solved approximately (they suggest multigrid) as the preconditioner for the JFNK solve of the
untouched hyperbolic system. They then say, at the end of §3.4.1:

> "In [43] a connection is made between the concept of the semi-implicit method as a preconditioner
> and the Schur complement of the Jacobian."

([43] = Chacón & Knoll 2003, which I could not obtain — §0.3 item 4.)

Chacón's 2025 Hall-MHD paper (arXiv:2407.07031, §3, read in full) restates and generalises the
mechanism, and is itself a readable primary source for it:

* §3 shows the 2×2 hyperbolic block system factorising into exactly the parabolic operator
  `I - Dt^2 d_xx`, and says "The connection between parabolization and the Schur factorization
  [3, 18] enables the generalization of these ideas to more complicated hyperbolic systems."
* §3.2 forms the exact Schur complement `P_Schur = D_v - L M^{-1} U`, declares it impractical
  because of `M^{-1}`, and replaces it by the **"small-bulk-flow approximation" `M^{-1} ~ Dt I`**,
  giving `P_SF = D_v - Dt L U`.
* The keyword list of that paper is literally "physics-based preconditioning, approximate block
  factorization".

**That is our derivation, term for term, in a different physics.** Our
`A_tt^{-1} ~ dt*I  =>  S_f = A_ff - dt*A_ft A_tf = div((eps + dt sigma) grad .)` is the same move
(`M^{-1} ~ Dt I`, then `D - Dt L U`), and Chacón's system is hyperbolic-stiff where ours is
relaxation-stiff. The general principle is 20+ years old and named.

### 1.3 What I could NOT find in any full text I read

I did **not** find, in anything I was able to read in full, the specific statement that
*the exact Schur complement of the coupled plasma / semiconductor drift-diffusion Jacobian on the
potential block is the semi-implicit Poisson operator, and that using it as the `SchurPre` matrix
in a fieldsplit preconditioner is what fixes convergence above `tau`.*

The two halves are separately standard (§1.1, §1.2); the composition, for THIS system, I did not
find. **My confidence that it is genuinely absent from the literature is LOW-to-MODERATE**, for
three concrete reasons, and I want that stated rather than softened:

* Bank, Chan, Coughran & Smith 1989 (alternate-block-factorization, *for semiconductors*) is
  precisely a block-factorisation preconditioner for this system and I could not read it.
* Lin, Shadid & Hoekstra (2008) is described by its own publisher record as covering "approximate
  block factorization and physics-based preconditioners" for drift-diffusion. I could not read it.
* Kerkhoven & Saad 1992 is the drift-diffusion split-as-preconditioner paper. I could not read it.

Any one of those three could contain it. **Get me those PDFs before anybody writes the word "novel"
in a paper.** What I *can* say now with confidence is the weaker and still useful thing: the
*operator* is not ours, the *principle* is not ours, and the *`dt/tau` ratio* is not ours (it is
Hagelaar's `chi_e`). If anything is ours it is narrow: the observation that the DEFAULT PETSc
`a11` Schur approximation is exactly the `chi_e -> 0` limit of the published operator, so that a
fieldsplit solver silently inherits the pre-1993 dielectric-relaxation time-step limit even though
the solver is "fully coupled". I have not seen that said anywhere, but see the caveat above.

### 1.4 A measured result that CUTS AGAINST our story and must not be buried

**ESTABLISHED (Duarte Gomez, Deak & Bisetti, *JCP* **480** (2023) 112007, §5.2 read in full).**
Their JFNK streamer preconditioner does the OPPOSITE of what we propose, and works.

Their approximate Jacobian, their eq. (85), deliberately **zeroes the entire `d(F_species)/d(phi)`
column** ("the linearization ... neglects the coupling between `dn, du, dv` and `dphi` because the
electric field ... is approximated by that associated with the potential `phi*`"). The resulting
preconditioner is block lower-triangular, and the potential is recovered LAST, their eq. (99), from
`L w_phi = (w_u - w_n - w_v) - r_phi`, where they state explicitly: "the resulting matrix `P_phi`
does not depend on the state vector `q`". `P_phi` is a **plain, constant Laplacian**, factorised
ONCE by MUMPS LU for the whole run.

In our language: their preconditioner's Schur complement on `phi` is exactly `A_ff` with the
conductivity term omitted — i.e. structurally the same thing as PETSc's default `a11`
approximation, which our diagnosis says must fail by the factor `dt/tau`.

And they report it works: 2-4 Newton iterations per step throughout, `GMRES_o/Newton` of 8.4-13.2
on the parallel-plate case; Table 3, BDF2/BDF3 at `h` = 10, 20, 30 ps on a 2.1M-point mesh
(8.4M unknowns), air at 300 K and 1 atm, background field 15 kV/cm, `u0` = 5e18 m^-3 seed. Their
abstract and §9 claim the method "overcomes traditional restrictions in the time step size due to
processes such as electron drift, electron diffusion, **and dielectric relaxation**".

**Their conditions, which matter:** finite differences on a structured/curvilinear mesh;
non-dimensionalised equations; 3 charged species and NO electron energy equation (LFA, ionisation
rate a function of reduced field via a fit); `P_n, P_u, P_v` preconditioned by HYPRE Euclid ILU(1),
`P_phi` by direct LU (MUMPS); 46% of total cost is the `P_phi` solve and ~85% of inner GMRES
iterations are the electron block. **They never report a measured `dt/tau_dielectric` ratio
anywhere** — the intro says the dielectric limit is "O(1e-11 - 1e-12 s)" and the runs use
10-30 ps, but they report `h/h_CFL-A` and `h/h_CFL-D` (Table 6) and not `h/tau`.

**MY INFERENCE, flagged as such and NOT to be quoted as theirs:** at 10-30 ps against a stated
dielectric scale of 10-12 ps, their `dt/tau` is order 1-3 — i.e. right around, and only modestly
above, where our Picard wall sits (0.86 to 1.7) and nowhere near the `dt/tau` = 10-300 regime our
Newton path is being asked to reach. **This is a hypothesis, not a reading**, and it is exactly the
kind of quantity the authors should be asked for rather than reconstructed. It is also the single
most important open question raised by this note: *does the plain-Laplacian Schur approximation
survive to `dt/tau` >> 1, or did Bisetti et al. simply never go there?* Their §7.2.4 remark that
"larger time steps result in a greater number of iterations as expected" is consistent with
degradation having already begun.

---

## 2. Q2 — What people actually use to precondition the fully-coupled system

Each entry: what it is, the regime it was DEMONSTRATED on, what it cost. All ESTABLISHED unless
marked.

### 2.1 Physics-based / approximate-block-factorization (Schur) preconditioning — Knoll, Keyes, Chacón lineage

* **What.** Take an existing segregated/semi-implicit solver, put it in delta form, and use it as
  the preconditioner for a JFNK solve of the untouched fully-implicit residual. Formally: eliminate
  the non-stiff block, approximate `A_tt^{-1} ~ dt I`, get a scalar parabolic operator on the
  remaining unknown, invert it approximately (multigrid). Knoll & Keyes §3.4.1; Chacón 2025 §3.
* **Demonstrated on.** 1D/2D shallow water with stiff gravity waves (Mousseau, Knoll & Reisner,
  *Mon. Weather Rev.* **130** (2002) 2611); non-equilibrium radiation diffusion (Mousseau, Knoll &
  Rider, *JCP* **160** (2000) 743); reduced resistive MHD (Chacón, Knoll & Finn, *JCP* **178** (2002) 15);
  Hall MHD (Chacón & Knoll 2003; Chacón 2025, up to 16384 MPI tasks in 2D); low-Mach reactor flows
  (Park et al., *JCP* **228** (2009) 9131). **Not, in anything I read, on plasma drift-diffusion +
  Poisson.**
* **Cost.** Knoll & Keyes make the key economic point directly: "Since the action of the true
  operator is maintained in the evaluation of the full nonlinear residual ... the inverse Jacobian
  used in the preconditioner can be further weakened without compromise to the solution in the
  interest of minimizing execution time. For instance, a few cycles of a multigrid method ... can be
  used to approximate the solution of Eq. (38)." I.e. the parabolic Schur solve is meant to be
  approximate and cheap, not converged. They also warn: "A JFNK wrapper can provide implicit balance
  to an inner operator-split solver ... **even if the operator split preconditioner does not allow
  one to use the large time step size that is otherwise achievable by the outer JFNK solver.**"

### 2.2 Structure-based operator splitting (transport block vs reaction block)

* **What.** Knoll & Keyes §3.4.2. With `J = aI + S + R` (`a = 1/dt`, `S` = convection-diffusion,
  `R` = pointwise reaction), apply `(aI + R)^{-1} a (aI + S)^{-1}`, which approximates `J` by
  `aI + S + R + a^{-1} S R` — first-order splitting error, small at small `dt`. The reaction block
  is pointwise-block-diagonal, the transport block is independent scalar solves per component.
* **Demonstrated on.** Radiation transport, semiconductor charge transport (their [10] = Bank et al.
  1989, [86] = Kerkhoven & Saad), subsurface two-phase flow, fluid flow.
* **Cost / warning.** "In cases in which there are strong couplings across different components
  stored at different mesh points, the type of phenomenon-based operator splitting described above
  is not expected to be successful." That is precisely the `n <-> phi` coupling, which is non-local.

### 2.3 The one paper that does exactly our problem: Duarte Gomez, Deak & Bisetti 2023

The most directly comparable published solver. Structure in §1.4 above. Concretely, per Krylov
iteration they solve five inner systems: `P_n`, `P_u`, `P_v` (one per species, each containing that
species' own drift + diffusion + self-reaction, ILU(1) via HYPRE Euclid), `Q` (pointwise 3×3 blocks
carrying the remaining inter-species kinetics, direct), and `P_phi` (constant Laplacian, MUMPS LU).
Cost breakdown, their Table 7: `P_n` 12%, `P_u` 8%, `P_v` 8%, `P_phi` **46%**, `Q` 2%, residual
evaluation 20%, GMRES orthogonalisation 3% — preconditioning is ~75% of total.

For scalability they replace the inner ILU/LU by GAMG / block-Jacobi / GASM for species and
**BoomerAMG for the potential** (their Table 6, 1.2M to 19.2M grid points, 56-896 cores, `h` = 10 ps,
homogeneous-Neumann species BCs): Newton/step stays 2.0-2.4 and `GMRES_o`/Newton grows only
2.1 -> 4.6 over a 16× problem-size increase, while `h/h_CFL-A` goes 2 -> 8 and `h/h_CFL-D` 1.1 -> 17.8.
Block Jacobi and GASM beat AMG on wall clock at these sizes (9.4 s vs 17.9 s per step at the
smallest size).

Their own honest limits: with physically realistic nonlinear electrode BCs (cathode sheath), the
pin-to-pin case degrades from `GMRES_o`/Newton ≈ 10-30 to **40-70** after ignition, and replacing
the anode BC by homogeneous Neumann gave "nearly identical streamer ignition and propagation with a
**sixfold decrease** in overall wall-clock time". Match that register when quoting them.

### 2.4 Full-Newton with an assembled Jacobian and an ILU-class preconditioner

* **What.** Arslanbekov & Kolobov (arXiv:2003.03812), §3.3: drift-diffusion (Scharfetter-Gummel)
  for electrons + ions, Poisson, AND an electron-energy equation with Joule heating — the same
  system as SoPlasma. Assembled Jacobian, allocated per time step, updated each Newton sub-iteration.
  Linear solve by a CGS-type Krylov method with high-order incomplete-decomposition preconditioning
  plus Kershaw diagonal modification, claimed complexity ~`O(N^{5/4})`. **No multigrid, no Schur, no
  fieldsplit.** They note TCAD/plasma linear systems are "highly ill-conditioned and therefore quite
  challenging for direct and preconditioned iterative solvers", and that nonsymmetric permutations
  with scaling reorderings gave the best results.
* **Demonstrated on.** 2D DC Ar cell, 400 mTorr, 2 cm gap, 200 V, gamma = 0.1, adaptive Cartesian
  mesh; 15000 steps, 10 fixed nonlinear sub-iterations, `dt` ramped to 200 ns over the first 5000.
* **Cost.** Not reported as wall clock per unit simulated time.

### 2.5 Dual-time / pseudo-transient with localized block-implicit relaxation, plus semi-implicit Poisson

* **What.** Zhu, Wu, Cao, Wei & Xu, arXiv:2607.18123 (2026), §1 and §2.4. BDF physical-time residual
  driven to zero by pseudo-time iteration; a **cell-local block-implicit (LU-SGS) relaxation** for
  the coupled species + electron-energy unknowns, so no global transport Jacobian is ever assembled;
  the potential advanced by the **semi-implicit Poisson** each inner iteration.
* **Their stated motivation, which is a fair criticism of our path.** "implicit Newton-Krylov methods
  robustly address the nonlinear coupling ... permitting large physical time steps. However, these
  methods require solving global nonlinear systems, computing intricate Jacobian (or Jacobian-vector)
  products, and designing highly specialized preconditioners. This inevitably leads to massive memory
  footprints and complicates the solver implementation."
* **Demonstrated on.** 1D argon RF CCP benchmark run on a 2D slab, plus a genuinely 2D case with a
  grounded transverse wall; pseudo-CFL 10000, `dt = T/100`, 100 pseudo-iterations/step; periodic
  steady state in ~1000 s wall clock. Raising inner iterations 100 -> 200 cost 1404 s -> 2679 s
  (1.9×) for a 0.82% change in cycle-averaged `n_e`.

### 2.6 Nonlinear / Krylov-as-preconditioner and inexact-Newton control

Knoll & Keyes §3.5-3.6: FGMRES with an inner Krylov preconditioner (their [156] = Saad 1993) is the
fully matrix-free option, and "In the JFNK context ... plain GMRES is surprisingly forgiving of
mildly inconsistent preconditioning." They also cover pseudo-transient continuation theory
(Kelley & Keyes, *SIAM J. Numer. Anal.* **35** (1998) 508) and nonlinear preconditioning (ASPIN).
**This is where our own measured Eisenstat-Walker result lives** — memory
`newton-rate-is-the-EW-forcing-term`, our measurement, not literature.

### 2.7 Vanka-style point-block smoothers — NOT found for this system

I found no full text demonstrating Vanka/Braess-Sarazin point-block relaxation on the coupled
Poisson + drift-diffusion plasma or semiconductor system. Where Vanka is demonstrated, it is
saddle-point physics: Stokes, resistive MHD (Adler, Benson, Cyr, MacLachlan & Tuminaro, *SIAM J.
Sci. Comput.* **38** (2016) B1, and the follow-on arXiv:2006.15700), poroelasticity, immersed boundary.
**I did not read those papers in full and so I report only that they exist and what they target.**
The structural reason Vanka fits those and may not fit ours is that they are saddle-point systems
with a constraint (incompressibility, solenoidality) and ours is not — that is MY INFERENCE and I
have no source for it.

---

## 3. Q3 — Multigrid on the COUPLED (phi, n) pair, and AMG on advection-dominated species

### 3.1 Classical AMG on advection-dominated operators: measured to FAIL, and the fix is known but scalar-only

**ESTABLISHED. Manteuffel, Ruge & Southworth, *SIAM J. Sci. Comput.* **40** (2018) A4105-A4130.**

Their §1 sets the framing in the same terms the parent agent used: with grid Reynolds number
`R_h = beta h / kappa`, "For `R_h > 1` the problem is advection-dominated ... which is well-discretized
by upwinding, resulting in a triangular or block-triangular matrix ... **AMG (and many other
iterative methods) are designed for elliptic, diffusion-like problems**".

Their measured verdict, Remark 3 of §4.2.1, quoted in full because the hedge matters:

> "With careful tuning of parameters, scalar, classical AMG methods do converge on most problems
> that are strongly diffusion dominated, `kappa >> 1/h`. However, such results were found to be very
> sensitive to parameter tuning, particularly for increased finite-element order and problem
> dimension, and typically did not outperform NSA or `lAIR`. **Classical AMG was unable to converge
> for all advection-dominated problems tested**, except for linear finite elements in two dimensions,
> with `kappa ~ 1/h`, where convergence factors were still worse than `lAIR` or NSA."

The fix they propose, `lAIR` (local approximate ideal restriction), is "robust from purely advective
to purely diffusive", 3-8× faster than nonsymmetric smoothed aggregation in work-per-digit on the
advection-dominated cases.

**Conditions, and they are restrictive for us.** The model problem is deliberately SCALAR
("a scalar model problem is chosen here intentionally to isolate the effects of nonsymmetry. Linear
systems with block structure resulting from a system of PDEs can be difficult for AMG, often
requiring individual attention that is outside the scope of this work"). The discretisations tested
are upwind DG and SUPG finite elements on structured and unstructured meshes — **not a finite-volume
Scharfetter-Gummel scheme**, and their strongest results lean on the triangular structure that upwind
DG produces. The implementation used was serial. And they note the opposite-direction failure: block
`lAIR` degrades in 3D diffusion-dominated cases.

**MY INFERENCE, labelled:** at our Courant numbers (20-55 per the memory
`newton-rate-is-the-EW-forcing-term`) the species block sits deep in the `R_h >> 1` regime where
this paper measures classical AMG to fail outright. That is a strong reason not to expect
`-fieldsplit_transport_pc_type gamg/hypre` to work on the species block, and a reason to look at
`lAIR` (available in hypre as `HYPRE_BoomerAMGSetRestriction`, which I have NOT verified in the
manual and should check before anyone tries it).

### 3.2 Is a coupled (phi, n) multigrid actually demonstrated? — I cannot answer this properly

**HONEST ANSWER: I could not obtain the papers that would settle it.** The three candidates
(Meza & Tuminaro 1996; Lin, Shadid & Hoekstra 2008; Lin & Shadid, *JCP* **228** (2009) 6250) are all
in §0.3. Their own publisher records describe fully-coupled Newton-Krylov with an algebraic
multilevel preconditioner applied to the full Jacobian of the drift-diffusion system — which, if it
reads as advertised, means **it is demonstrated and not aspirational**, at least in the
semiconductor (steady/quasi-steady, moderate-Courant) setting. **I will not conclude that from a
publisher blurb.** Get me those PDFs.

What I *can* report from full texts:

* `lAIR`'s own conclusions, 2018: "a future research direction of interest is extending `lAIR` to
  systems of PDEs. **Systems of PDEs remain something of an open question for AMG in general.**"
* Knoll & Keyes §3.3: multigrid as a JFNK *preconditioner* is more robust than multigrid as a
  *solver*, and a suboptimal multigrid (e.g. piecewise-constant prolongation, violating the
  `m_P + m_R > 2` order rule) that fails as a solver can still give a scalable preconditioner
  (their [100], [110]). **That is the right way to attempt this: never as a standalone solver.**
* Duarte Gomez et al. use BoomerAMG successfully — but only on the DECOUPLED, constant, scalar
  Poisson block `P_phi` (their Table 6), not on any coupled operator. That is a single-field
  elliptic AMG solve and should not be cited as coupled multigrid.
* hypre's own manual, §5.7.6, on the coupled case: the default unknown-based approach "can work well
  for weakly coupled PDEs, but **will be problematic for strongly coupled PDEs. For such problems, we
  recommend to use hypre's multigrid reduction (MGR) solver**." That is the library authors telling
  us systems-AMG is not the recommended tool for a strongly coupled system.

### 3.3 The dof-ordering question — the premise is half right

**ESTABLISHED (primary documentation, read directly):**

* **hypre BoomerAMG does NOT require interlaced ordering.** `HYPRE_BoomerAMGSetDofFunc(solver, dof_func)`
  "Sets the mapping that assigns the function to each variable, if using the systems version. **If no
  assignment is made** and the number of functions is `k > 1`, **the mapping generated is
  `(0,1,...,k-1,0,1,...,k-1,...)`**" — i.e. interlaced is only the DEFAULT, and an arbitrary
  (including field-major) ordering is supported by supplying `dof_func` explicitly. The one place
  the manual does hard-require it is `HYPRE_BoomerAMGSetFilterFunctions`, whose note reads: "This
  option assumes that variables are stored in an interleaved format."
* **PETSc's GAMG effectively DOES require it.** `PCGAMG` Notes: "To obtain good performance for
  PCGAMG for vector valued problems you must call `MatSetBlockSize()` to indicate the number of
  degrees of freedom per grid point", and `MatSetBlockSize` documents block storage
  (`MATBAIJ`/`MATSBAIJ` "ALWAYS have square block storage"), which means the `k` dofs of a grid
  point must be contiguous — that is interlaced ordering.

**So: the parent agent's statement "we would need to switch to interlaced for systems-AMG to work at
all" is TRUE for PETSc GAMG and for hypre's function-filtering option, and FALSE as a blanket
statement for hypre BoomerAMG.** If the route chosen is BoomerAMG systems/nodal AMG, the ordering
change may be avoidable by supplying `dof_func`; if the route is GAMG, it is not. I have seen no
paper comment on the ordering either way — this is documentation, not literature, and I have
labelled it so.

---

## 4. Q4 — Is the Gummel-vs-Newton crossover characterised?

### 4.1 What is established from full texts I read

* **The Gummel map IS the physics-based preconditioner archetype.** Knoll & Keyes §3.4 list
  "Gummel's method for the semiconductor drift-diffusion equations [70]" alongside SIMPLE and the
  semi-implicit method as the canonical examples of a legacy segregated algorithm that should be
  reused as a preconditioner rather than replaced. [70] = Gummel, *IEEE Trans. Electron Devices*
  **ED-11** (1964) 455. **So "our segregated outer loop IS the Gummel map" is correct and standard,
  and the recommended disposition of a Gummel map in a JFNK code is: keep it, as the preconditioner.**
* **The general mechanism of the wall is stated in the plasma literature.** Arslanbekov & Kolobov,
  §3.4: "non-linear convergence within each time step is largely controlled by the selected time step
  (`Dt`), **which defines the diagonal dominance of the Jacobian matrix**." And: "it is not possible to
  provide a general recipe for selecting the maximum allowed time step `Dt` a priori, contrary to the
  explicit solvers where the time step can be most often estimated accurately from the CFL condition
  and from the dielectric relaxation time controlling the particle-field coupling."
  **Read carefully, that is an explicit statement that no published a-priori criterion exists for the
  implicit case** — which is the closest thing to a direct answer to Q4 that I have from a full text.

### 4.2 What I cannot answer, and why

**A published criterion of the form "switch from Gummel to Newton when `dt/tau` exceeds X" — I did
not find one, and I could not read the three papers most likely to contain the steady-state
analogue** (Kerkhoven 1986, Kerkhoven 1988, Kerkhoven & Saad 1992; §0.3 items 2-3). The classic
qualitative wisdom the parent agent cites (Gummel at low injection / weak coupling, Newton at high
injection / strong coupling) I encountered only in search summaries and secondhand descriptions,
**which under A4 I may not build on and am not building on.** I am not going to dress a search
snippet up as a literature answer.

What I will note, as **MY INFERENCE and nothing more**, is that the two criteria are dimensionally
the same object. Our `dt/tau = dt*sigma/eps0` measures how far the `-dt div(sigma grad .)` term in
the Schur complement dominates `div(eps grad .)`; "high injection" in a device means large carrier
density, hence large `sigma`, hence the same ratio at the (effectively infinite) steady-state `dt`.
**That equivalence is a conjecture until Kerkhoven is read.** Our measurement stands on its own:

**OUR OWN MEASUREMENT (parent agent, this session, restated not re-derived):** on a developed
streamer, 81,640 cells, peak `n_e` 1.19e19, `tau` = 1.16e-10 s, limiters off, common endpoint —
Picard converges in 2-4 correctors at `dt/tau <= 0.86` and SIGFPEs on step 1 at `dt/tau >= 1.7`.
That is a sharply located crossover and I have found nothing published that predicts its location
for a transient plasma fluid model. **If Kerkhoven turns out not to contain a transient criterion,
locating this crossover at `dt/tau` with a derivation from the Schur complement is a publishable
result in its own right — arguably a cleaner one than the preconditioner.**

---

## 5. Q5 — Does beating the dielectric relaxation limit actually pay? The literature says: only sometimes, and it names the condition.

This is the question the parent agent said matters most, and it has the best-supported answer in
this note.

### 5.1 The published criterion

**ESTABLISHED. Teunissen, *PSST* **29** (2020) 015010, §3.6** (local file
`Literature/flux-schemes-and-numerics/Teunissen_2020_*.pdf`, article p. 5-6, checked on the rendered
page):

> "Avoiding the time step restriction due to the dielectric relaxation time `tau` is beneficial when
> `tau` becomes smaller than other time step constraints, in particular when `tau < tau_CFL`."

with the two situations he lists where that happens: high-density low-pressure discharges (high
conductivity), and a localized high-density/high-conductivity region such as near an electrode.

**That is the whole economics in one sentence: beating `tau` buys nothing unless `tau` is the
BINDING constraint.**

### 5.2 The same point, measured, in the paper the user already owns

**ESTABLISHED. Villa et al., *JCP* **242** (2013) 86-102, §6.3**, rendered page 99:

> "Till now we have described some test cases where some small spatial scales are used. The time
> step, in those cases, is quite small. In fact, the satisfaction of the Courant Friedrichs Lewy
> condition implies that the time step is linearly constrained by the size of the mesh. In other
> words, **in all the cases we have shown till now the difference of performances between the
> implicit and explicit time integration schemes is not really relevant.**"

Only in their §6.3 benchmark — a 1 m box, i.e. a coarse mesh where the CFL is loose — does the
implicit scheme show a decisive advantage (the explicit scheme "in a few steps develops a global
instability"). Their streamer runs sit at `dt` = 1e-9 to 1e-10 s against `tau ~ 1e-11 s`, i.e.
`dt/tau` ~ 10-100 (p. 99, checked on the rendered page).

### 5.3 The same point again, from the full-Newton plasma side

**ESTABLISHED. Arslanbekov & Kolobov, arXiv:2003.03812 §3.4:**

> "Since an implicit time step costs more than an explicit time advance, one can expect that the FNM
> scheme will be most efficient for steady-state or slowly-varying problems by allowing large time
> steps. At the same time, **transient problems evolving on the (fast) electron-time scale can be
> more efficiently solved by explicit techniques.**"

Their own runs bear it out: during the avalanche/breakdown phase they are forced to `dt` = 0.05-0.1 ns
"because of rapidly evolving plasma dynamics", and only once the dynamics move to the ion time scale
do they ramp to a few ns (RF) or 200-500 ns (DC).

### 5.4 The accuracy half of the argument, which is usually left out

**ESTABLISHED. Teunissen 2020 §3.4, figure 3** (checked on the rendered page). 1D test, zero source,
`mu_e` = 0.03 m^2/(V s), `D_e` = 0.1 m^2/s, 10 mm domain, `dx` = 20 um, 10 kV applied,
`n_e = n_p` = 1e20 m^-3 over 4-6 mm, giving `tau` ≈ 18.4 ps:

> "The semi-implicit method shows roughly first-order convergence. Errors are significantly larger
> than with the current-limited approach, **also for time steps larger than the dielectric relaxation
> time `tau`**."

From his figure 3, the semi-implicit RMSE in `E` runs ~1e-3 to 3e-2 MV/m across `dt` = 1-100 ps,
against ~1e-8 to 1e-3 MV/m for his current-limited scheme — several orders worse at equal `dt`, on
both sides of `tau`. **So the semi-implicit device buys stability past `tau` at a real and measured
accuracy price, and it is only first-order.** (Our own paper reaches a compatible conclusion from a
different direction, arXiv:2607.05137 §6.3.1: "the semi-implicit Poisson coupling alone cannot
overcome the need for outer PIMPLE iterations", with the same corrector sensitivity as the explicit
case, at a case where the dielectric relaxation ratio was only ~0.3 and convective Courant ~0.23.)

### 5.5 The one measured result that points the other way

**ESTABLISHED. Duarte Gomez et al. 2023 §7.2.4, figure 8 and surrounding text:**

> "The measurement of the wall-clock time per step shows that smaller time steps ease the task for
> the iterative linear solver, but such a positive improvement does not scale linearly, as shown by
> the overall wall-clock time required to advance the solution to 4 ns. **Smaller time steps imply a
> more costly simulation overall**, although as expected, it is easier to offset such cost through
> the reuse of factorizations when compared to larger time steps."

Their optimum: `h × f` ≈ 0.6 at `h` = 0.3 (dimensionless), i.e. ~15 preconditioner factorisations
per nanosecond during streamer propagation. **Conditions: `h` restricted to 10-30 ps; BDF2/BDF3;
2.1M points; the accuracy ceiling on `h` was set by their own temporal convergence study (BDF2/BDF3
recover their design order over 10-30 ps against a 3 ps BDF3 reference).** They are trading inside a
band where larger `dt` is still accurate; the relevant lever in their study is factorisation REUSE,
which is a cost our matrix-free path does not have in the same form.

### 5.6 Verdict on our measurement

**OUR OWN MEASUREMENT (memory `newton-rate-is-the-EW-forcing-term`):** at a common physical window,
`dt` = 1e-9 vs 2e-10 on the same bed, cost per ns of simulated time is identical to ~1%.

**MY INFERENCE, clearly labelled:** that is exactly what Teunissen's criterion predicts for a bed
that is ALREADY past its advective CFL. The same memory records that bed running at `Co_conv(e)`
20-37 and `Co_conv(energy)` 30-55 at `dt` = 1e-9. If `tau_CFL` is already violated by an order of
magnitude, then `tau` is not the binding constraint and, per §5.1, beating it "buys nothing" —
which is precisely the ~1% that was measured. **This is a hypothesis with a cheap test attached:
repeat the `dt` sweep on a bed where `tau < tau_CFL` (a high-conductivity, coarse-mesh, low-pressure
configuration — exactly Teunissen's two listed cases) and see whether the payoff appears.** If it
does, the preconditioner work is motivated; if it does not, it is not. That test is worth more than
any further preconditioner development and it should be run first.

Note also the confound to control for (A2): the `dt` = 1e-9 arm converged only 8.3% of its steps and
the `dt` = 2e-10 arm 100%. Comparing cost-per-ns across arms with a 12× difference in step-success
rate is comparing two different algorithms, not two time steps.

---

## 6. Consolidated citation list

Full texts read (§0.1) are marked ✅; obtained-not-read (§0.2) ◐; not obtained (§0.3) ✗.

* ✅ G. J. M. Hagelaar, *Modelling methods for low-temperature plasmas*, HDR thesis, Université
  Paul Sabatier, Toulouse. Ch. 8, eqs. (8.5)-(8.16), pp. 48-50.
  `/home/kkourtza/Projects/SoEEDF/Literature/plasma-fluid-closures-LMEA-LFA/hdr-hagelaar.pdf`
* ✅ A. Villa, L. Barbieri, M. Gondola, R. Malgesini, "An asymptotic preserving scheme for the
  streamer simulation", *J. Comput. Phys.* **242** (2013) 86-102. Eqs. (9)-(11), §6.2-6.3.
  `/home/kkourtza/Projects/SoEEDF/Literature/Villa_1-s2.0-S0021999113001265-main.pdf`
* ✅ D. A. Knoll, D. E. Keyes, "Jacobian-free Newton-Krylov methods: a survey of approaches and
  applications", *J. Comput. Phys.* **193** (2004) 357-397. §3.3-3.6, eqs. (32), (38)-(42).
  `/home/kkourtza/Projects/SoEEDF/Literature/Knoll_1-s2.0-S0021999103004340-main.pdf`
* ✅ A. Duarte Gomez, N. Deak, F. Bisetti, "Jacobian-free Newton-Krylov method for the simulation of
  non-thermal plasma discharges with high-order time integration and physics-based preconditioning",
  *J. Comput. Phys.* **480** (2023) 112007, DOI 10.1016/j.jcp.2023.112007. §5, eqs. (84)-(99),
  Tables 3, 6, 7.
  `.../Literature/preconditioning-and-coupled-solvers/DuarteGomez_Bisetti_2023_JFNK-streamer-preconditioner_JCP480-112007_acceptedMS.pdf`
* ✅ J. Teunissen, "Improvements for drift-diffusion plasma fluid models with explicit time
  integration", *Plasma Sources Sci. Technol.* **29** (2020) 015010. §3.4-3.6, figure 3.
  `/home/kkourtza/Projects/SoEEDF/Literature/flux-schemes-and-numerics/Teunissen_2020_Plasma_Sources_Sci._Technol._29_015010.pdf`
* ✅ L. Chacón, "A scalable multidimensional fully implicit solver for Hall magnetohydrodynamics",
  *J. Comput. Phys.* **526** (2025) 113789; preprint arXiv:2407.07031. §3-3.2, eqs. (12)-(18).
  `.../Literature/preconditioning-and-coupled-solvers/Chacon_2025_HallMHD-implicit-Schur-preconditioner_arXiv-2407.07031.pdf`
* ✅ T. A. Manteuffel, J. Ruge, B. S. Southworth, "Nonsymmetric algebraic multigrid based on local
  approximate ideal restriction (lAIR)", *SIAM J. Sci. Comput.* **40** (2018) A4105-A4130,
  DOI 10.1137/17M1144350; preprint arXiv:1708.06065v4. §1, §4.2, Remark 3, §5.
  `.../Literature/preconditioning-and-coupled-solvers/Manteuffel_Ruge_Southworth_2018_lAIR-nonsymmetric-AMG_arXiv-1708.06065.pdf`
* ✅ R. Arslanbekov, V. Kolobov, "Implicit and coupled multi-fluid solver for collisional
  low-temperature plasma", arXiv:2003.03812 (2020). §3.3-3.4. *(A related paper by the same authors
  exists as PSST **30** (2021) 045013 under a different title; I have NOT verified that they are the
  same work and am citing only the preprint I read.)*
  `.../Literature/preconditioning-and-coupled-solvers/Arslanbekov_Kolobov_2020_full-Newton-multifluid-LTP_arXiv-2003.03812.pdf`
* ✅ Y. Zhu, H. Wu, J. Cao, Y. Wei, K. Xu, "A low-storage implicit dual-time finite-volume framework
  for radio-frequency capacitively coupled plasma fluid simulations", arXiv:2607.18123 (2026). §1, §3.
  `.../Literature/preconditioning-and-coupled-solvers/Zhu_2026_implicit-dual-time-CCP-fluid_arXiv-2607.18123.pdf`
* ✅ I. Pasolari, K. Kourtzanidis, "SoPlasmaFoam: ...", arXiv:2607.05137. Eqs. (2)-(3), §6.3.1.
* ✅ hypre Documentation, Release 3.2.0, §5.7.6 and the BoomerAMG reference entries; PETSc 3.25
  manual pages `PCGAMG`, `MatSetBlockSize`.
* ◐ G. Chen, L. Chacón, C. Leibs, D. Knoll, W. Taitano, "Fluid preconditioning for Newton-Krylov-based,
  fully implicit, electrostatic particle-in-cell simulations", *J. Comput. Phys.* **258** (2014) 555;
  preprint arXiv:1309.6243. **Obtained, not yet read.**
* ◐ G. Lapenta, S. Ju, "Predictor-corrector preconditioners for Newton-Krylov solvers in fluid
  problems", arXiv:0804.0150 (2008). **Obtained, not yet read — read this next.**
* ◐ J. Brannick, R. Falgout, K. Kahl, J. Schroder, T. Shen, "Nodal AMG coarsening and interpolation
  for PDE systems", arXiv:2501.16248 (2025). **Obtained, abstract + §1 only.**
* ✗ R. Bank, T. Chan, W. Coughran, R. Smith, "The alternate-block-factorization procedure for systems
  of partial differential equations", *BIT Numer. Math.* **29** (1989) 938-954, DOI 10.1007/BF01932753.
* ✗ T. Kerkhoven, Y. Saad, "On acceleration methods for coupled nonlinear elliptic systems",
  *Numer. Math.* **60** (1991) 525-548, DOI 10.1007/BF01385735.
* ✗ T. Kerkhoven, *SIAM J. Numer. Anal.* **23** (1986) 1121-1137; *SIAM J. Sci. Stat. Comput.* **9**
  (1988) 48-60; *SIAM J. Numer. Anal.* DOI 10.1137/0725073.
* ✗ L. Chacón, D. A. Knoll, "A 2D high-β Hall MHD implicit nonlinear solver", *J. Comput. Phys.*
  **188** (2003) 573-592, DOI 10.1016/S0021-9991(03)00193-1; L. Chacón, *Phys. Plasmas* **15** (2008) 056103.
* ✗ P. T. Lin, J. N. Shadid, "Performance of a parallel algebraic multilevel preconditioner for
  stabilized finite element semiconductor device modeling", *J. Comput. Phys.* **228** (2009) 6250-6267,
  DOI 10.1016/j.jcp.2009.05.024; P. T. Lin et al., *Int. J. Numer. Meth. Engng* (2012),
  DOI 10.1002/nme.4315; P. T. Lin, J. N. Shadid, R. J. Hoekstra (2008).
* ✗ J. C. Meza, R. S. Tuminaro, "A multigrid preconditioner for the semiconductor equations",
  *SIAM J. Sci. Comput.* **17** (1996) 118-132.
* ✗ P. L. G. Ventzek, T. J. Sommerer, R. J. Hoekstra, M. J. Kushner, *Appl. Phys. Lett.* **63** (1993)
  605-607; J.-P. Boeuf, *J. Appl. Phys.* **63** (1988) 1342-1349.
* ✗ H. K. Gummel, *IEEE Trans. Electron Devices* **ED-11** (1964) 455-465. *(Cited by Knoll & Keyes
  [70]; not read here.)*

Also cited from within full texts read, but not themselves read: Crispel, Degond & Vignal,
*JCP* **223** (2007) 208; Langdon, Cohen & Friedman, *JCP* **51** (1983) 107; Scharfetter & Gummel,
*IEEE Trans. Electron Devices* **ED-16** (1969) 64; Adler, Benson, Cyr, MacLachlan & Tuminaro,
*SIAM J. Sci. Comput.* **38** (2016) B1.

---

## 7. What this does NOT say

* It does not say our Schur derivation is wrong. §1.1-1.2 say each of its two ingredients is
  published; the derivation in `docs/design/schur-semiimplicit-poisson-preconditioner.md` is
  untouched by this note and its falsifiable prediction (quality degrades with the transport CFL,
  is indifferent to `k_eff`) is still the right test to run.
* It does not establish that our proposed `SchurPre` will work. Duarte Gomez et al. (§1.4) reach
  10-30 ps on 8.4M unknowns with the OPPOSITE choice, and until somebody reports their `dt/tau`
  that contrast is unexplained.
* It does not assess novelty. Three unread papers (§1.3) could each contain the claim.
* It does not say implicit is or is not worth it for SoPlasma. It says the literature's condition
  for "worth it" is `tau < tau_CFL` (§5.1) and that our benchmark bed does not appear to satisfy it
  — a hypothesis with a stated test, not a conclusion.

---

# PART B — SECOND PASS, 2026-09-12. The five paywalled PDFs, read in full.

**Status: LITERATURE NOTE (appended), 2026-09-12.** No code, no numbers handed to the solver.
The user obtained five PDFs specifically to close §0.3 items 1, 2 and 5. All five are now read
**in full**, and this part supersedes the parts of Part A that rested on not having them.
Every equation quoted below was checked on the **rendered page**, not on `pdftotext` output —
Bank et al. is an OCR'd scan whose text layer flattens sub/superscripts (it renders `(z⁰)ᵀ` as
`(z°) r`, and eq. (30) as `A = [-qr  T + eS]`), so nothing from its text layer is quoted here.

**Supersession notices, per D2 (originals left in place, unedited):**

* §0.3 item 1 (Bank et al.) — **SUPERSEDED BY §B.2.** Obtained and read in full.
* §0.3 item 2 (Kerkhoven & Saad) — **SUPERSEDED BY §B.3.** Obtained and read in full.
* §0.3 item 5 (Lin et al.) — **SUPERSEDED IN PART BY §B.4.** The *JCP* 228 (2009) paper and the
  2008 WCCM8 item are read; the 2012 *IJNME* paper is still not obtained.
* §0.2 (Lapenta & Ju, "read this next") — **SUPERSEDED BY §B.6.** Read in full.
* §1.3 ("my confidence that it is genuinely absent is LOW-to-MODERATE") — **SUPERSEDED BY §B.2.5.**
  Two of the three named risks are now discharged.
* §3.2 ("I cannot answer this properly") — **SUPERSEDED BY §B.4.** It can now be answered.
* §4.2 — **SUPERSEDED IN PART BY §B.5.** Two quantitative *steady-state* criteria now exist in
  hand; a transient one still does not.

Also newly present in the folder and **NOT read, not built on**:
`Boeuf_1988_2D-model-dc-glow-discharge_JAP63-1342.pdf` (§0.3 item 7) and
`Wang_Fan_1995_existence-uniqueness-steady-semiconductor-equations_AMS15-180.pdf`.

---

## B.1 Evidence inventory for this pass

| source | pages | what I read | rendered-page checks |
|---|---|---|---|
| Bank, Chan, Coughran & Smith, *BIT* **29** (1989) 938-954 | 17 | **all of it**, §1-§5 + refs | journal pp. 940, 941, 942, 943, 944, 945, 946, 947, 948 |
| Kerkhoven & Saad, *Numer. Math.* **60** (1992) 525-548 | 24 | **all of it**, §1-§7 + App. A, B, C | journal p. 544 (Thm A.2) |
| Lin, Shadid, Sala, Tuminaro, Hennigan & Hoekstra, *JCP* **228** (2009) 6250-6267 | 18 | **all of it**, §1-§7 | journal p. 6263 (the §6.3 elimination) |
| Lin, Shadid & Hoekstra, WCCM8/ECCOMAS 2008 | 2 | **all of it** | — |
| Lapenta & Ju, arXiv:0804.0150 (2008) | 13 | **all of it**, §1-§5 | — |

**First correction, and it matters.** Item 4 of the brief, `Lin_Shadid_Hoekstra_2008_...pdf`, is a
**two-page conference extended ABSTRACT** (its PDF Title is literally `Microsoft Word -
paul-abstract.doc`). It contains **no results at all**: "Results *will be* presented...",
"Additional results *will be* presented...". It therefore **cannot** settle Q2 and Part A's §1.3
should not have listed it as a novelty risk on the strength of its publisher description. What it
*does* contribute is one dated negative datum, quoted in full:

> "A block Jacobi, block Gauss–Seidel, and block SOR preconditioner has already been implemented.
> **More advanced preconditioners such as use of the Gummel iteration [4] as a preconditioner as
> well as use of operator-splitting techniques are presently under consideration.**"

**ESTABLISHED:** as of mid-2008, the Sandia group had *not* built Gummel-as-preconditioner or an
operator-split/approximate-block-factorization preconditioner for drift-diffusion; they were
considering it. (Their reference [1] for the motivation is Knoll & Keyes 2004.)

---

## B.2 Q1 — BANK et al. 1989. VERDICT: **not prior art for our composition; a genuine ancestor of one half of it.**

### B.2.1 What they actually factor — and it is not a Schur complement

**ESTABLISHED.** The system is eqs. (3)-(5), p. 939, **steady and elliptic**, primitive variables:

> (3) `L1(u,n,p) = -∇²u + n - p - N(x) = 0`
> (4) `L2(u,n,p) = ∇·(n∇u - ∇n) = 0`
> (5) `L3(u,n,p) = -∇·(p∇u + ∇p) = 0`

with a quasi-Fermi variant (8)-(10) via `n = e^{u-v}`, `p = e^{w-u}`. **There is no time derivative
anywhere in the paper.** `grep -i` over the full text returns zero hits for *Schur*, *dielectric*,
*semi-implicit*, *time step*; "transient" occurs only inside two reference titles.

The method, defined at p. 942-943 and checked on those rendered pages. Write the Newton correction
system in PDE-blocked form (17), `A_ij ∈ R^{ν×ν}`. Then

> (18) `D` is the block matrix with `D_ij = diag(A_ij)`
> (19) `(A D^{-1})(D x) = b`
> "`D^{-1}` is the ABF postconditioner."
> (22) `D̃ = P D P^T` — the same `D`, re-blocked by GRID POINT into `ν` dense `m×m` blocks

So **ABF is a RIGHT preconditioner (they call it a postconditioner) by the inverse of the matrix of
pointwise DIAGONALS of the Jacobian's blocks**, applied so that a block Gauss-Seidel iteration can
then be run on the postconditioned system. Its cost is "the inversion of `ν` matrices of order `m`"
(p. 941). **It never inverts, nor approximates the inverse of, any block `A_tt`.** It is a local
`m×m` change of variables at each grid point — "a temporary local change of variables" (p. 941),
justified by "if we assume that the coupling between the PDEs is largely localized by grid point"
(p. 942). Their own placement of it in the literature (§5, p. 953) is: similar to element-by-element
preconditioners, and "**a special instance of τ-transforming smoothers**" (Wittum, *Numer. Math.*
**54** (1989) 543). They never mention Schur complements, block factorization in the
`D - L A^{-1} U` sense, or any physics-based/semi-implicit operator.

For the one-carrier drift-diffusion pair they give, eq. (27), p. 944:

> `A = [ -Δ   I ;  -M   C ]`

with `-Δ` the discrete Laplacian (Poisson row), `I = ∂L1/∂n`, `-M = ∂L2/∂u` the linearised **drift**
operator `∇·(n̄∇·)` (`M` spd because `n̄ > 0`), `C = ∂L2/∂n` the convection-diffusion operator.
Their postconditioned matrix, eq. (28)-(29):

> (28) `A D^{-1} = [ (-Δ diag(C) + diag(M))δ ,  (-diag(Δ) + Δ)δ ;  (-M diag(C) + C diag(M))δ ,  (-C diag(Δ) + M)δ ]`
> (29) `δ = (-diag(Δ) diag(C) + diag(M))^{-1}`

### B.2.2 The family resemblance, stated precisely — MY INFERENCE, with the algebra

The `(1,1)` block of (28) is **the exact Schur complement on the potential block with `A_21` and
`A_22` replaced by their diagonals, right-scaled**. With `S_exact = A_11 - A_12 A_22^{-1} A_21 =
-Δ + C^{-1} M` and `S_diag = A_11 - A_12 diag(A_22)^{-1} diag(A_21) = -Δ + diag(C)^{-1} diag(M)`:

```
   -Δ diag(C) + diag(M)  ==  S_diag · diag(C)          and    δ^{-1} == diag(S_diag) · diag(C)
```

both identities exact (verified symbolically/numerically to 1.4e-14 on random `5×5` blocks;
script at `/tmp/claude-1000/.../scratchpad/check_abf.py`). The same is **not** true of `S_exact`
(residual 2.3, i.e. not an identity). **This is MY INFERENCE from their equations, labelled as
such. Bank et al. never write it and never use the word Schur.**

So the *elimination* underlying our derivation has a 1989 ancestor for this exact PDE system, in
diagonal-approximated, postconditioning form. That is a real citation obligation and it must go in
any paper we write.

### B.2.3 What plays the role of our `dt·σ` — nothing, and the operator is a different one

**ESTABLISHED.** Their model problem, eq. (30), p. 945 (rendered):

> `A = [ T   I ;  -ηT   T + εS ]`,  `T = h^{-2}[-1  2  -1]`,  `S = h^{-1}[-1  1  0]`
> "`η` corresponds to a carrier density, say `n`, while `-ε` corresponds to the electric field, `E = -∇u`."

and the postconditioned matrix, eq. (35):

> `A D^{-1} = σ^{-1} [ (1 + εh/2)T + ηI ,  I - h²T/2 ;  εη(S - hT/2) ,  (1 + ηh²/2)T + εS ]`,
> `σ = (2/h²)(1 + εh/2 + ηh²/2)`

**The operator ABF produces on the potential row is `(1 + εh/2)T + ηI` — a Laplacian plus a
ZEROTH-ORDER (mass/screening) term whose coefficient is the carrier density `η`.** Ours is
`∇·((ε₀ + dt σ)∇·)` — a **modified diffusion coefficient inside the divergence**. These are
structurally different operators, and **MY INFERENCE** of why is exact and checkable: their `A_tt`
is the steady transport operator `C` (second order), so `A_ft A_tt^{-1} A_tf → diag(C)^{-1}diag(M)`
is a *number* (`~ η`), whereas ours is `A_tt = I/dt + transport`, so `A_tt^{-1} ≈ dt·I` and
`A_ft A_tt^{-1} A_tf ≈ -(dt/ε₀)∇·(σ∇·)` **retains the differential operator**. The entire difference
comes from the presence of `∂/∂t`. Their steady analogue of the dielectric-relaxation term is
**Debye screening**; ours is **dielectric relaxation**. `dt·σ` has no counterpart in their paper
because there is no `dt`.

### B.2.4 Conditions of their results — state them, do not generalise them

**ESTABLISHED**, §3-§4. The spectral analysis (§3) is an *ad hoc* 1D Fourier analysis on an
equispaced mesh that "ignores boundary conditions and other important issues"; they say so
themselves. Fig. 1 is drawn at a **fixed frequency** (`c = cos(kπh) = 0.9`) and `εh = 2`, plotted
against `ηh²`. They state "for realistic semiconductor simulations" `ηh²` can be `10³` or more
while `εh = O(1)` (p. 951). The numerical experiments (§4) use **sparse DIRECT solvers for every
linear system** — no multigrid, no Krylov method is used anywhere in their results. Cases: a 1D
one-carrier resistive bar; a 2D two-carrier resistive slab, scaled doping `10⁴`-`10⁸`, coupled-ABF
6 nonlinear iterations vs plug-in ~50; and Table 1, a small 2D MOS transistor in saturation at
**1163 and 2765 vertices** — coupled-ABF 41 nonlinear / 123 linear at 1163 and 26 / 78 at 2765,
against plug-in 42 / 217 and 41 / 200. Their own limits, quoted so we match their register: on a
forward-biased pn junction with **low** doping, where "diffusion effects are more prominent",
"the coupled-ABF approach has more difficulty... It is sometimes necessary to do 12 or 15 block
Gauss-Seidel iterations for each Newton iteration. Further experiments on pn junctions and bipolar
transistors are needed." And: "like other preconditioners, **ABF is not a panacea**" (p. 950, after
the breakdown example (53)-(57), where an infinite eigenvalue appears at `α = √(1+μ)`).

### B.2.5 VERDICT ON Q1

**Our composition is GENUINELY DISTINCT from Bank et al., and it is NOT a special case of theirs.**
Four independent reasons, each sufficient:

1. **No time discretisation exists in their paper.** Nothing can play the role of `dt·σ`, and the
   `dt/τ` framing has no referent in a steady elliptic system.
2. **ABF is not a Schur-complement preconditioner.** It inverts only pointwise `m×m` diagonals; it
   never forms, nor approximates, `A_ft A_tt^{-1} A_tf` as an operator.
3. **The operator produced is different in kind** — a zeroth-order screening term `ηI` versus our
   modified diffusion coefficient `∇·((ε₀+dt σ)∇·)` (§B.2.3).
4. **They never identify it with any pre-existing physics-based/semi-implicit algorithm.** Their own
   analogies are element-by-element preconditioners and τ-transforming smoothers.

**But the novelty claim must be narrowed, and this is the real finding of this pass.** The
*elimination itself*, applied to the semiconductor drift-diffusion Jacobian to produce a better
operator on the potential row, is 1989 prior art in diagonal-approximated form (§B.2.2). What
remains ours to claim is the **transient** composition: that with `A_tt = I/dt + transport` the
same elimination yields exactly the **published semi-implicit Poisson operator**, that its
coefficient is exactly Hagelaar's `χ_e = dt/τ_d`, and that PETSc's default `a11` Schur
approximation is precisely its `χ_e → 0` limit. Phrase it that way, cite Bank et al. as the
steady-state ancestor, and it survives.

---

## B.3 Q1 (continued) — KERKHOVEN & SAAD 1992. VERDICT: **not prior art; it is the Gummel-as-NONLINEAR-preconditioner paper.**

**ESTABLISHED.** The system, eqs. (1.1)-(1.3), is again **steady**, in quasi-Fermi variables. `T` is
the Gummel map: "a nonlinear block Gauß-Seidel iteration on the above system... `v` is updated by
solving (1.1), then `w` by solving (1.2), and finally `u` from the third equation" (p. 526), named
as Gummel's method again in Appendix A. `NLGMR` is an inexact Newton method applied to
`u - T(u) = 0`, Jacobian-free via the difference quotient (3.1).

The paper contains **no Schur complement, no block factorization, and no linear preconditioner at
all** (`grep`: zero hits for "Schur"). Their central claim is the opposite of preconditioning the
linear system, and they say so, p. 527:

> "the solution of Newton's equations `[I - T_u] du = -[u - T]` by GMRES **does not require
> preconditioning**... this version of Newton's method for the solution of the nonlinear elliptic
> system is **already preconditioned in that an implicit form of preconditioning is embedded in the
> operator `T`**. Our techniques can be viewed alternatively as a preconditioned nonlinear Krylov
> subspace method **where the preconditioning consists of one step of the nonlinear fixed point
> algorithm**."

The theory: `T_c` is compactly differentiable, so by Theorem 2.1 (from Chatelin) `σ(T_u)` accumulates
only at 0, hence `σ(I - T_u)` clusters at 1; their Theorem 5.1 then proves that a residual-minimising
Krylov method converges **superlinearly** under the clustering rate `M(ε) ≈ Kε^{-r}` (5.3), with the
bound (5.4). **That is the real content of the paper and it is a statement about the Gummel map's
spectrum, not about an operator on the potential block.**

**Conditions, and they are narrow.** 2D N-MOSFET, geometry of their Fig. 1; **20 × 27 mesh**;
Scharfetter-Gummel discretisation with the Poisson nonlinear terms "lumped" on the diagonal;
constant mobility 820 cm²/(V·s), `n_i` = 1.4e10 cm⁻³, `ε_Si` = 11.7, 300 K, background doping
3e15 cm⁻³, device 3 μm × 2.8 μm, oxide 250 Å; biases backgate 0, source 0.5 V, gate 6.5 V, drain
6.5 V; zero generation/recombination; IBM 4381. Results over 100 iterations: unaccelerated factor
`1e-1`; Chebyshev `ρ_ex` = 0.922 (theory 0.9036 from the fitted ellipse `d` = 0.475, `c` = 0.475,
`a` = 0.505); stationary 2nd-order essentially identical; **NLGMR `ρ_ex` = 0.838**, overall "a factor
larger than 8" in computational time. `m` adapted between 2 and 25.

**Relevance to us, and it is real but different from what Part A guessed.** Kerkhoven & Saad is the
published ancestor of **nonlinear preconditioning by a Gummel sweep** — which is Lapenta & Ju's
construction (§B.6) and *not* our `PCSHELL`. It is not a competitor to the Schur claim.

---

## B.4 Q2 — is a COUPLED `(φ, n)` multigrid demonstrated, or aspirational? **DEMONSTRATED — under conditions that are not ours.**

Source: Lin, Shadid, Sala, Tuminaro, Hennigan & Hoekstra, *JCP* **228** (2009) 6250-6267, read in full.

### B.4.1 Do they coarsen the coupled system? YES.

**ESTABLISHED**, §5.2. They aggregate on the graph of the **nodal blocks of the full Jacobian**:

> "A nodal block refers to the submatrix which couples all degrees of freedom defined at the same
> grid node. **In our case, the nodal block dimension is m = 3 corresponding to the electrostatic
> potential, electron concentration, and hole concentration unknowns.**"

`P_ℓ` is built from (12) as **piecewise-constant interpolation over each aggregate for each of the
three solution components** (nonsmoothed aggregation), columns normalised by a QR of each `I_ℓ^s B_ℓ`;
`R_ℓ = P_ℓ^T`; coarse operators by Galerkin projection `A_{ℓ+1} = R_ℓ A_ℓ P_ℓ`. Aggregates come from
METIS/ParMETIS on the nodal-block graph. So this is genuine **systems/nodal AMG on the fully-coupled
`(ψ, n, p)` Jacobian** — not block-separate preconditioning. It answers Part A §3.2 directly.

### B.4.2 What smoother? **Not Gauss-Seidel, not Vanka — one-level additive Schwarz with ILU.**

**ESTABLISHED**, §5.2 and §6: `S_ℓ(A_ℓ,u,b): repeat m times  u ← u + M̃_AS(b - A_ℓ u)` with `M̃_AS`
being their eq. (11) additive Schwarz with the local inverses replaced by `ILU(k)`. Production
settings throughout §6: **ILU(2), one level of overlap**, on fine and medium levels; **KLU direct
solve** on the coarse level; three levels; aggressive coarsening (50-175 nodes per aggregate, optimum
80-125). Their own explanation, and it is the sentence to quote:

> "**The drift-diffusion equations require a heavyweight smoother such as ILU.** Without aggressive
> coarsening, the next coarser level will be large, and ILU will also be expensive on this level...
> **We are pursuing physics-based preconditioning methods which will allow us to use less expensive
> smoothers than ILU** (for examples of physics-based preconditioners, see [10])."

([10] = Knoll & Keyes 2004.) So as of 2009 the physics-based/Schur route was explicitly **future
work** at Sandia too. That is now two dated statements (2008 abstract, 2009 paper) that the group
best placed to have done it had not done it.

### B.4.3 The CONDITIONS — mesh, regime, transient or steady, and `dt`

**ESTABLISHED.** Discretisation: **stabilised finite element** — SUPG-type plus a nonisotropic
discontinuity-capturing term, eqs. (7)-(9), on **uniform quadrilateral** meshes (they note the code
handles unstructured meshes but the test cases do not use them). **Not** a Scharfetter-Gummel finite
volume scheme: "A FV Scharfetter-Gummel technique is currently under development, and a future study
will evaluate the performance of the multilevel preconditioner on this discretization technique."
They also disclose that in the scaling studies "there is sufficient resolution so that the
**discontinuity capturing terms above are not employed**".

Three cases:

1. **Steady-state 2D NPN BJT**, 2 × 1.5 μm, silicon, **0.3 V bias**, max donor doping `1e19`, max
   acceptor `1e16`; initial guess = the nonlinear-Poisson solution. Weak scaling 110K → **112M
   unknowns**, 4 → 4096 cores of Red Storm. **All runs required seven Newton steps.**
2. **Transient 2D diode**, 1 × 0.5 μm, **sinusoidal potential on one contact, amplitude 0.5, period
   1.0** (scaled), first-order **backward Euler**, **fixed** `Δt` = 0.01, 0.05, 0.1 → **50, 10 and 5
   time steps** to half a period. Initial condition = the steady state at zero bias. 637K / 2.54M /
   10.2M unknowns.
3. **Pseudo-1D NP diode**, 1 × 0.125 μm, symmetric Gaussian doping, **zero bias**, doping swept
   `1e16`-`1e20` (and `1e16`-`1e22` for the BJT).

**This is near-equilibrium / quasi-steady device operation throughout.** Nothing in the paper is an
avalanche, a streamer, or a fast transient. And, verified by `grep` over the full text: the words
**"dielectric", "relaxation time", "Courant" and "CFL" do not appear**. "Debye" appears exactly once,
defining `λ` in `λ² = εV₀/(q x₀² C₀)`, "the minimal Debye length of the device".

**Consequently their `dt/τ` is NOT reportable from this paper.** `Δt` is given only in the scaled
unit `t₀ = x₀²/D₀` with `D₀ = max(D_n, D_p)`, and **neither `D₀` nor the mobilities are stated
anywhere in the paper.** I am not going to reconstruct it from another paper's silicon mobility —
that is exactly the inference rule A4 forbids. **If we need their `dt/τ`, the authors should be
asked, or the companion paper [28] (Lin, Shadid et al. on variational multiscale vs
Scharfetter-Gummel) obtained.**

### B.4.4 The performance, honestly — it is NOT optimal and they say so

**ESTABLISHED.** Weak scaling, Table 5 (steady BJT, W(1,1), agg85, ILU(2)/overlap 1):

| cores | fine unknowns | 1-level DD ILU iter/Newton | 3-level iter/Newton | 3-level lin-solve time/Newton |
|---|---|---|---|---|
| 4 | 110K | 68 | 36 | 3.8 s |
| 64 | 1.75M | 287 | 111 | 9.5 s |
| 1024 | 27.9M | 1145 | 219 | 25 s |
| 4096 | 112M | 2264 | **295** | 53 s |

> "While the three-level method is **clearly not scaling optimally, both in iteration count and CPU
> time**, these results are encouraging" — and about 20× faster than one-level at 112M. The one-level
> method shows the theoretical `N^{1/2}` growth.

Transient, Table 7, 3-level agg60 V(1,1): iterations/Newton **84-154**, Newton steps/`Δt` 2.9 (at
`Δt`=0.01) to 4.8 (at `Δt`=0.1), ~10× faster than one-level at 10.2M. Their statement:

> "While for a fixed time step an h independent convergence rate is not achieved, only a moderate
> increase is evident. **For a given resolution the iteration count is relatively constant as a
> function of time step size.**"

**Conditions on that last sentence, which is the one we would most want to quote:** it covers a **10×
`Δt` range (0.01 to 0.1)**, three meshes, one device, near equilibrium, with the Newton count itself
rising 2.9 → 4.8 over that range. It is not a statement about large `dt`.

Doping sensitivity: 1D diode (Table 8) iterations *decrease* with doping, 198 → 113 at 25.2M over
`1e16`→`1e20`, "only a factor of about two variation... over four orders of magnitude"; 2D BJT
(Table 9, 27.9M) is non-monotonic, peaking at **356** iterations at `1e19` doping and falling to 136
at `1e22`, with Newton steps rising as "the problem is becoming more nonlinear".

### B.4.5 Does it survive Manteuffel, Ruge & Southworth Remark 3? **The paper does not answer the question, and Remark 3 does not directly apply.**

**MY INFERENCE, labelled, with the reasons checkable:**

* Remark 3 of Manteuffel et al. 2018 is about **classical (Ruge-Stüben C/F) AMG on SCALAR** model
  problems. Lin et al. use **nonsmoothed aggregation on a 3-dof nodal system with a heavyweight
  ILU-Schwarz smoother**. Different AMG family, different problem class — Remark 3 is not a
  counter-example to their result and must not be quoted as one.
* **Lin et al. never place their runs on the advection-dominance axis.** No grid Reynolds number,
  no Courant number, no CFL, no `Pe_h` is reported anywhere. Their SUPG stabilisation explicitly
  "improves the conditioning of, and therefore the iterative solution of, the Jacobian matrices"
  (§3) — i.e. added numerical diffusion is part of why AMG works here, and its magnitude is not
  reported.
* They concede the nonsymmetric gap in their own words (§5.2): `R_ℓ = P_ℓ^T` "is almost always done
  for symmetric problems though **for nonsymmetric systems alternatives may be warranted** (see
  [41])... We intend to explore other possibilities in a future paper based on the AMG algorithm
  described in [42]."

**So Q2's honest answer is: coupled `(φ, n, p)` AMG is DEMONSTRATED at 10⁸ unknowns, non-optimally,
on steady and slow-transient silicon devices with a stabilised-FE discretisation and an ILU smoother
— and it is UNTESTED in the advection-dominated regime where our species block sits.** The two
papers together (Manteuffel's scalar failure, Lin's system success with heavy stabilisation and ILU)
do not contradict each other; they simply do not overlap. Anyone claiming either result transfers to
SoPlasma's Scharfetter-Gummel FV species block at `Co` 20-55 would be quoting both above their tier.

### B.4.6 One thing in Lin et al. that IS close to our derivation, and must be cited

**ESTABLISHED**, §6.3, the displayed equation on journal p. 6263 (checked on the rendered page;
note the extracted text corrupts `λ` into `k`). To explain doping sensitivity they linearise the
electron equation about `(ψ̄, n̄, p̄)`, 1D, constant `(λ, μ_n, D_n)`, and **substitute `∇²ψ` from the
potential equation into the electron equation**:

> `∂n/∂t + ū_n ∂n/∂x - D_n ∂²n/∂x² + ( ∂G/∂n|_{n̄,p̄} - (μ_n/λ²)[p̄ - n̄ + C] ) n = 0`
> with `ū_n = μ_n ∂ψ̄/∂x = -μ_n Ē_x = -(μ_n/λ²)∫[p̄ - n̄ + C]dx`
> "In this formulation the linearized electron equation is in a standard transient
> convection–diffusion-**reaction** form."

**This is the same elimination we perform, in the OPPOSITE direction** — they remove `ψ` from the
transport row; we remove `n` from the potential row. The reaction coefficient `(μ_n/λ²)[p̄-n̄+C]` is,
in dimensional terms, `μ_n ρ_net/ε` = the divergence of the drift velocity, i.e. an inverse
dielectric-relaxation rate built on the **net space charge** (not on `n_e`; it vanishes where the
plasma is quasineutral). **MY INFERENCE** that it is the dielectric-relaxation rate in disguise —
they never call it that, the words do not occur in the paper, and the `[p̄-n̄+C]` factor genuinely
distinguishes it from `σ/ε`.

They use it **only as a qualitative conditioning argument**, never to build a preconditioner:
"However the nonlinearity and the coupling of the original drift-diffusion system would make a
complete analytical understanding of this system difficult to obtain. For this reason a numerical
study will be presented." **No `dt` appears in the coefficient.** It is prior art for *the idea of
substituting Poisson into the linearised transport equation*, and it is not prior art for the Schur
operator on the potential block.

---

## B.5 Q3 — Gummel↔Newton crossover. **Two published QUANTITATIVE criteria exist. Neither is transient, and neither is `dt/τ`.**

### B.5.1 The qualitative statement, and the standard architecture

**ESTABLISHED**, Bank et al. p. 940 (rendered):

> "If the PDEs in (2) are weakly coupled, the convergence can be quite rapid. **When the PDEs are
> strongly coupled, the convergence of the outer plug-in iteration can be quite slow, or even
> diverge.** For the semiconductor problem in quasi-Fermi variables, (8)-(10), **the plug-in method
> converges well when the iterates are far from the solution [14]; device simulation programs often
> use a few initial plug-in iterations to improve an initial guess before switching to a coupled
> approach.**"

([14] = Kerkhoven, *SIAM J. Sci. Stat. Comput.* **9** (1988) 48-60, still not obtained.) Kerkhoven &
Saad say the same from the other side, p. 527: "Iteration with the nonlinear mapping `T`... converges
rapidly far away from the solution `u*`, but slows down once `u*` is approached... NLGMR
acceleration... is likely to be needed close to the solution only."

**This is direct published support for the Picard-warm-up → Newton-handover architecture we already
run**, and it is a stronger statement than "Newton needs a good initial guess": in this system the
segregated map is *better* far from the solution and *worse* near it, which is the reverse of the
usual intuition. It is also why `anySpeciesOnFloor()`-style handover gating is the right shape of
mechanism even though our particular gate is unsatisfactory (PROGRESS Task 2).

### B.5.2 Criterion 1 — Bank et al. eq. (50): diagonal dominance, parameter `ηh²`

**ESTABLISHED**, p. 947 (rendered). For block Gauss-Seidel applied directly to the **primitive-variable**
Newton correction matrix `A` of (30) — which they state at p. 941 gives "the same linear equations as
those solved in the plug-in method if Newton's method were used to deal with the scalar PDEs", so it
is the Gummel rate for the primitive formulation:

> (50) `ρ_PV = | h²η / ( (1-c)(2+hε) + hεis ) |`
> "which is **monotonically increasing in `η`**. In order for `ρ_PV < 1`, it is necessary that
> `ε = δη` for `δ` sufficiently large or, in other words, **`A` must be diagonally dominant**. These
> arguments suggest that **Newton-Gauss-Seidel applied to (3)-(5) can diverge**, which is in
> agreement with the results of [21]."

Contrast, eq. (52), the **quasi-Fermi** formulation: `ρ_QF < 1` always, and `ρ_ABF < ρ_QF`. Their Fig. 1
plots all three against `ηh²` at fixed `c = 0.9`, `εh = 2`, and the dashed `ρ_PV` curve crosses 1
near `ηh² ≈ 1`.

**This is a real quantitative divergence criterion for the primitive-variable segregated map — which
is the formulation our Picard loop uses — and its parameter is carrier density × mesh² with the field
`ε` in the denominator. There is no `dt` in it.**

### B.5.3 Criterion 2 — Kerkhoven & Saad Theorem A.2: a bound that never exceeds 1

**ESTABLISHED**, journal p. 544 (checked on the rendered page; the text layer mangles it). For the 1D,
constant-doping, steady, quasi-Fermi system (A.1) on `[0,L]` with `u(0)=0`, `u(L)=V_B`:

> **Theorem A.2.** `ρ  ≤  [ (N+P) / ((π/L)² + N + P) ] · [ 1 / √(1 + (2π/V_B)²) ]`

with `N`, `P` the scaled electron and hole densities at the solution (the sentence immediately above
the theorem says the bound is "in terms of the densities `n` and `p` at the solution and the device
length"; the paper switches to capitals inside the theorem without redefining them — I read them as
the same quantities, and flag the inconsistency rather than hide it). They add: "By Kitchen's Theorem
([34]), **Gummel's method converges locally for arbitrary variation of the 'bias potentials'**."

**Read carefully, this bound is always `< 1`. It predicts SLOWDOWN — `ρ → 1` as `(N+P)` grows against
`(π/L)²` and as `V_B` grows — and it never predicts divergence.** That is consistent with Bank's
`ρ_QF < 1`: in quasi-Fermi variables Gummel converges, slowly; in primitive variables it can diverge.
Our Picard loop is in primitive variables and it **diverges** (SIGFPE), so **Bank's (50), not
Kerkhoven's A.2, is the right published analogue.**

**MY INFERENCE, labelled, and easy to check or refute:** in the scaling that makes (1.3) and (3)
have unit coefficients, lengths are in intrinsic Debye lengths, so `(N+P)/((π/L)² + N+P) =
1/(1 + π²λ_D²/L²)` and Bank's `ηh² = (Δx/λ_D)²`. Both published criteria are then **spatial**
ratios — device length, or mesh spacing, measured in Debye lengths. Our `dt/τ_d` is the **temporal**
ratio, and the two are related exactly by the Einstein relation, since `τ_d = ε₀/(e μ n) = λ_D²/D`:

```
    dt/τ_d  =  (dt·D/Δx²) · (Δx/λ_D)²  =  Fo_diff × (Bank's ηh²)
```

i.e. they are the same object up to the diffusive Fourier number. **Neither paper states this, and
neither paper states the scaling I used** — Bank et al. never say their lengths are in Debye lengths
(Lin et al. do define `λ` as the minimal Debye length, in a different scaling). **This is a
conjecture with a stated derivation, not a reading, and it should be confirmed against Selberherr
or Bank-Rose-Fichtner before it goes in a paper.**

### B.5.4 VERDICT ON Q3

* **Is the crossover characterised anywhere in these five texts? Qualitatively yes; quantitatively,
  only in STEADY spatial parameters (`ηh²`, `(L/λ_D)²`, bias), never in `dt` or `dt/τ`.**
* **Is `dt/τ` the stated parameter anywhere? NO.** None of the five papers contains the phrase
  "dielectric relaxation"; Bank and Kerkhoven contain no time step at all; Lin et al. report `Δt` but
  never a `τ`; Lapenta & Ju report `Δt` for a diffusion and a cavity problem with no `τ` in them.
* **The bridge that nobody has built.** Bank et al.: the segregated map diverges when the Jacobian
  loses diagonal dominance (eq. 50). Arslanbekov & Kolobov (Part A §4.1, read in full): "non-linear
  convergence within each time step is largely controlled by the selected time step (`Δt`), which
  defines the diagonal dominance of the Jacobian matrix", and "**it is not possible to provide a
  general recipe for selecting the maximum allowed time step `Δt` a priori**". **Both halves are
  published; the quantitative join — that for the transient Poisson + drift-diffusion system the
  diagonal dominance is lost at `dt ≈ τ_d`, derivable from the Schur complement — is not.**

**OUR OWN MEASUREMENT (PROGRESS §1, 2026-09-12, restated not re-derived):** `$HOME/streamer-warm`,
81,640 cells, warm start `t` = 1e-9 from a developed streamer, peak `n_e` = 1.19e19,
`τ = ε₀/(e μ_e n_e)` = 1.16e-10 s, limiters OFF, common endpoint `t` = 2e-9 — Picard reaches the
endpoint at `dt/τ` = 0.09, 0.43, 0.86 and **SIGFPEs on step 1** at `dt/τ` = 1.7 and 4.3. The
named control is the `dt` = 1e-11 arm reaching the endpoint, which proves the restart is sound.
On the 449k bed the absolute `dt` ceiling moves (that bed is 2.35× finer, so the drift Courant number
is 2.35× higher at equal `dt`) — **only the ratio is claimed to transfer, and that has not yet been
verified on the 449k bed.**

**So: locating this crossover at `dt/τ ≈ 1`, deriving it from the Schur complement, and connecting
it to Bank's steady `ηh²` criterion through `τ_d = λ_D²/D`, remains unclaimed in the literature I
have read in full — and on the evidence of this pass it is a cleaner and better-defended
contribution than the preconditioner.** Two caveats before anyone writes it: it rests on one mesh
and one gas, and Kerkhoven 1988 (both papers) is still unread.

---

## B.6 LAPENTA & JU 2008 — what they retrofit, and whether it is our `assembledPmat false` PCSHELL

### B.6.1 What they do

**ESTABLISHED**, §3. They require a semi-implicit scheme that is **LINEAR in the new state**, eq. (6):

> `A x¹ + f_SI(x⁰) = 0`, "where `A` is a linear operator (matrix) and the function `f_SI` depends
> only on the old state `x⁰`".

They then note that the **classic** way to use it is eq. (7), `A δx = r_k`, "the matrix `A` of the
semi-implicit scheme becomes the preconditioner matrix `P` for the Jacobian matrix `J`" — and that
this is what they are **NOT** doing. Their construction is eq. (8):

> (P) `A x¹ + f_SI(x*) = 0`   (C) `r = f(x⁰, x¹)`

**The Newton unknown is changed from `x¹` to `x*`, a fictitious modified OLD state.** GMRES is run
**unpreconditioned** on the composed residual `F(x*) = f(x⁰, -A^{-1} f_SI(x*))`. They prove
first-order equivalence to right preconditioning, eqs. (9)-(12), ending at `J A^{-1} A δx = r_k`:
"To first order in the Taylor series expansion, the new approach is identical to applying the
traditional preconditioners directly to the Jacobian equation. However, to higher order this might
be a better approach as it reduces the distance between the initial guess (`x⁰`) and the solution
for `x*`." The selling point is software engineering, in their words: "representing perhaps a greater
advance in software engineering than in computational science", because the existing code "operates
on full fields, with their boundary conditions, **not on variations**".

### B.6.2 Does it correspond to our `assembledPmat false` PCSHELL? **NO — it is the construction they explicitly contrast themselves against.**

**ESTABLISHED + MY INFERENCE, separated:**

* **ESTABLISHED:** our PCSHELL applies one segregated sweep as `P^{-1}` acting on a Krylov **variation**
  `δx` inside the linear solve. That is *verbatim* their eq. (7), the "classic implementation" they
  describe and then set aside. Their own method is at a different level: a **nonlinear** right
  preconditioning / nonlinear elimination, with no preconditioner inside GMRES at all.
* **ESTABLISHED:** their predictor is a **single linear solve** with a fixed matrix `A` (a tridiagonal
  solve in §4.1; a linearised vorticity update plus one stream-function elliptic solve in §4.2). It is
  **not** a Picard iteration driven to convergence.
* **MY INFERENCE:** their eq. (12) is the standard statement that the preconditioned operator is
  `J A^{-1}`. If the *linearised* fixed-point iteration for the same splitting diverges — i.e.
  `ρ(I - A^{-1}J) > 1` — then `A^{-1}J` has eigenvalues far from 1 and **GMRES preconditioned by `A`
  must be expected to converge slowly or stall**, which is the shape of what we measured. So our
  structural explanation is *consistent with* standard splitting theory, but **Lapenta & Ju neither
  state nor test it**: the words "diverge", "Poisson", "plasma", "drift" and "dielectric" do not occur
  anywhere in their paper (`grep`-verified), and they never report a failure of the predictor.
* **A caution against our own story, and it should be recorded:** a divergent *nonlinear* Picard
  iteration and a poor *linear* preconditioner `A` are not the same object. One application of
  `A^{-1}` to a Krylov vector can be perfectly well defined at a `dt` where iterating
  `x ← -A^{-1}f_SI(x)` blows up. **"The PCSHELL IS a Picard sweep and Picard SIGFPEs, so it inherits
  the divergence" is therefore a plausible mechanism, not an established one** — the SIGFPE is a
  property of the nonlinear iterate trajectory (negative densities, clamps), while `DIVERGED_ITS` at
  1000 with a frozen SNES residual is a property of the spectrum of `A^{-1}J`. The discriminating
  instrument would be to measure the eigenvalue spread of `A^{-1}J` (or simply the GMRES residual
  history, which will be flat for a stalling preconditioner and slowly decreasing for a merely bad
  one) rather than to argue from the Picard failure. That is a `diagnostician`/`numerical-analyst`
  job, and it is cheap.

### B.6.3 Their conditions, which are far from ours

**ESTABLISHED.** Two benchmarks, **neither a plasma, neither Poisson-coupled drift-diffusion**:

1. **1D nonlinear diffusion**, `∂φ/∂t = ∂/∂x(D(φ)∂φ/∂x)`, `D = α₀ + α₁φ`, `α₀ = α₁/10`, `α₁ = 1`,
   `L = 4`, Crank-Nicolson corrector, predictor lagging `D` at the **old time level** `φ⁰`, `Δt` = 0.1
   fixed to `t` = 1.0, 100-800 cells, `η_r = η_a = 1e-5`, Eisenstat-Walker, Kelley's NK code.
   Preconditioned: Newton 2.82-3.00, **Krylov 3.18-4.67, essentially mesh-independent**; unpreconditioned
   Krylov 15.7 → 109.4 over the same meshes; 7× CPU saving at 800 cells.
2. **2D driven cavity**, vorticity-streamfunction, `Re` = 1000, 10-60 cells, `Δt` = 0.01-0.1, backward
   Euler. Preconditioned: "GMRES is never actually called", 1 preconditioner application per Newton
   iteration, Newton 2.00 → 5.31; unpreconditioned Krylov 79 → 2732 and **at `N`=60, `Δt`=0.1 the
   unpreconditioned run did not converge at all** while the preconditioned one did.

One detail of theirs that bears directly on our design: in the cavity case the predictor's lagged
velocity is taken at the **previous Newton iterate**, not the old time level — "We remark that using
the old velocity rather than the previous guess from the Newton iteration results in **much poorer
performances**." If we ever build their construction, the semi-implicit `σ` must be evaluated at the
current Newton iterate, not at `t^n`.

Note also their own measured trend, which is the opposite of a "large `dt` is free" story: Newton
iterations rise with `Δt` at fixed mesh (2.00 → 5.31 at `N`=60 over `Δt` 0.01 → 0.1).

**VERDICT:** Lapenta & Ju is a real, readable primary source for *retrofitting a semi-implicit code
as a JFNK preconditioner*, it is a legitimate citation for the PCSHELL family, and its construction
is an **alternative we have not built** (nonlinear preconditioning on `x*`) rather than a description
of what we have. **It neither predicts nor contradicts our PCSHELL failure above `τ`** — it contains
no stiff-source, no elliptic-constraint, and no relaxation-stiff problem at all.

---

## B.7 What these five texts do NOT settle, and the exact papers that would

1. **A transient Gummel↔Newton criterion.** The two remaining Kerkhoven papers are Bank's [14] and
   [15] and are the quantitative steady theory: **Kerkhoven, "On the effectiveness of Gummel's
   method", *SIAM J. Sci. Stat. Comput.* **9** (1988) 48-60**, and **Kerkhoven, "A spectral analysis
   of the decoupling algorithm for semiconductor simulation", *SIAM J. Numer. Anal.* **25** (1988)
   1299-1312** (note: Part A §0.3 gave this second one as *SIAM J. Numer. Anal.* **23** (1986)
   1121-1137 / DOI 10.1137/0725073 — **Bank's reference list gives volume 25, pp. 1299-1312, 1988**,
   and the 1986 item is a separate paper on realistic device geometries. Both citations should be
   checked before use). If neither contains a `dt`-parameterised criterion, §B.5.4 stands.
2. **Whether Sandia later built the physics-based/block-factorization preconditioner they said they
   were pursuing.** The follow-on is **Lin, Shadid et al., *Int. J. Numer. Meth. Engng* (2012), DOI
   10.1002/nme.4315**, plus their companion discretisation paper, their reference [28]
   (variational-multiscale FE vs Scharfetter-Gummel FV for drift-diffusion). **These are now the
   single largest remaining novelty risk for Q1** and I am asking for them.
3. **The primary source for "semi-implicit method ≡ Schur complement of the Jacobian".** Still
   **Chacón & Knoll, *JCP* **188** (2003) 573-592**. Knoll & Keyes attribute it there and I have only
   Chacón's 2025 restatement.
4. **Lin et al.'s `dt/τ_dielectric`.** Not reportable from *JCP* 228 — `D₀` and the mobilities are not
   in the paper. Ask the authors, or get [28].
5. **Meza & Tuminaro, *SIAM J. Sci. Comput.* **17** (1996) 118-132**, still not obtained.

## B.8 Bottom line for the three decisions

* **Q1 — novelty.** The Schur-on-`φ` *composition for the transient system* survives Bank et al. and
  Kerkhoven & Saad: neither has a time step, neither forms a Schur complement, and Bank's operator is
  a screening term where ours is a modified diffusivity. **But the claim must be narrowed and must
  cite Bank et al. 1989 as the steady-state ancestor of the same elimination** (§B.2.2) **and Lin et
  al. 2009 §6.3 as the same substitution in the opposite direction** (§B.4.6). Do not write "novel"
  until item B.7.2 is read.
* **Q2 — coupled AMG.** Demonstrated, not aspirational, at 10⁸ unknowns — with nonsmoothed
  aggregation, an ILU(2)/Schwarz smoother, aggressive coarsening, a stabilised-FE discretisation, and
  near-equilibrium device physics; non-optimal scaling (36 → 295 iterations over a 1000× size range);
  and **no reported advection-dominance measure anywhere**, so it transfers to our species block only
  as a hypothesis.
* **Q3 — the crossover.** Two published *steady, spatial* criteria (Bank eq. 50; Kerkhoven & Saad
  Thm A.2), **no published transient criterion, and `dt/τ` appears nowhere.** On this pass's evidence
  the `dt/τ ≈ 1` crossover, derived from the Schur complement, is the stronger of our two candidate
  contributions.
