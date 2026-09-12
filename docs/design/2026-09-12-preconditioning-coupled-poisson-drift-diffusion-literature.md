# Preconditioning the fully-coupled Poisson + drift-diffusion system: what the literature already says

**Status: LITERATURE NOTE, 2026-09-12.** No code, no numbers handed to the solver. This note
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
