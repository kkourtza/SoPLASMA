# The semi-implicit Poisson operator IS the Schur complement

**Status: DERIVED, not yet implemented.** Design note, 2026-09-10. The user's
reaction on being shown the idea: "i like the semi-implicit poisson as a
preconditioner idea!"

## The claim

The semi-implicit Poisson operator

    div( (eps + dt*sigma) grad phi )

is not merely *analogous* to a Schur complement of the coupled Newton system.
It **is** that Schur complement, to leading order in `dt`, when the transport
unknowns are the ones eliminated. This is an identity to be checked, not a
metaphor.

## Derivation

Order the coupled Newton system as `[phi; t]`, with `t` the transport block
(all species densities plus, under LMEA, the electron energy):

    J = [ A_ff   A_ft ]        f = phi
        [ A_tf   A_tt ]        t = transport

with, from the residuals this solver actually assembles:

* `A_ff = d(F_phi)/d(phi) = div(eps grad .)` -- the Poisson operator.
* `A_ft = d(F_phi)/d(n_i) = q_i` -- the Gauss residual carries `rho = sum q_i n_i`,
  so this block is multiplication by the charge per particle. **Diagonal.**
* `A_tf = d(F_n)/d(phi) = -div(Z_i mu_i n_i grad .)` -- the drift term's field
  derivative. This is the block whose absence made the Schur complement inert
  (experiment log section 10) and which is now assembled as
  `-fvm::laplacian(Z*mu*n, dphi)`.
* `A_tt = d(F_n)/d(n)` -- advection-diffusion-reaction, with `ddt` contributing
  `1/dt` on the diagonal.

**Two Schur complements are available**, depending on which block is eliminated:

    eliminate phi :  S_t = A_tt - A_tf * A_ff^-1 * A_ft     (acts on transport)
    eliminate t   :  S_f = A_ff - A_ft * A_tt^-1 * A_tf     (acts on phi)

PETSc currently forms the FIRST (splits are registered `phi` then `transport`,
and its Schur eliminates split 0). The second is the interesting one.

Take `A_tt` in the `dt`-dominant regime, where the time derivative dominates
transport and reaction:

    A_tt ~ I/dt      =>      A_tt^-1 ~ dt*I

Substitute into `S_f`, for a single charged species, using the blocks above:

    A_ft * A_tt^-1 * A_tf  ~  q * dt * ( -div(Z mu n grad .) )
                           =  -dt * div( q Z mu n grad . )
                           =  -dt * div( sigma grad . )

because `sigma = q Z mu n` is exactly the electrical conductivity (summed over
species for the multi-species case -- the sum falls out of the species sum in
`rho`). Hence

    S_f  ~  div(eps grad .) + dt*div(sigma grad .)
         =  div( (eps + dt*sigma) grad . )                    <-- QED

which is the semi-implicit Poisson operator, term for term.

## Why this resolves the apparent contradiction in section 23

Section 23 established, and measured, that the semi-implicit Poisson is
**inconsistent inside a Newton residual**: across 3167 steps no step whose
dielectric relaxation ratio exceeded 0.3278 ever converged, and removing it
took SNES failures 4.92% -> 0% with dt 3.4e-12 -> 7.4e-11.

Both facts are now explicable at once, and they are not in tension:

* **As a RESIDUAL it is wrong** because the semi-implicit substitution PREDICTS
  `rho^{n+1} ~ rho^n - dt*div(sigma E)`, and a Newton residual already carries
  `rho^{n+1}` as an unknown. The prediction and the unknown double-count the
  same charge relaxation, so `F(u) = 0` has no root. It is a SEGREGATED device:
  its purpose is to avoid needing `rho^{n+1}`, which a segregated solver does
  not have.
* **As a PRECONDITIONER it is exactly right**, because approximating the Schur
  complement is precisely what a fieldsplit Schur preconditioner is for, and
  the operator IS that Schur complement to leading order in `dt`.

So the semi-implicit Poisson was never a wrong idea -- it was in the wrong
place. Wrong operator to put in the residual; right operator to put in the
preconditioner.

**This sharpens the deferred AP proof rather than competing with it.** The AP
property of the semi-implicit scheme and the AP property of the fully coupled
Newton solve become the same fact seen from two sides: the coupled system
*contains* the semi-implicit operator as its Schur complement, so a solver that
inverts the coupled system accurately is doing what the semi-implicit
reformulation does explicitly -- without paying its consistency error. See
memory `ap-proof-newton-removes-dielectric-constraint`.

## Why it should also be FAST

The current arrangement puts the expensive object in the wrong place, and this
is measured, not supposed (experiment log 25c, `-ksp_view` plus counts on
`grubert_long`):

| | current: Schur on transport | proposed: Schur on phi |
|---|---|---|
| Schur block | `S_t`, (nSpecies+1) x nCells | `S_f`, 1 x nCells |
| its preconditioner | ILU(0), natural ordering, on `A_tt` | **assembled `div((eps+dt sigma)grad)`, AMG-ready** |
| inner (eliminated) solve | `A_ff` -- Poisson, exact LU, 2000 rows, cheap | `A_tt` -- advection-diffusion-reaction |
| measured cost | 167 solves/step, median 230 its, 2.6% hit the 10,000 cap | -- |

`S_f` is a scalar elliptic operator on one field. That is the ideal case for
algebraic multigrid, and unlike `S_t` it can be ASSEMBLED explicitly rather
than only applied -- which is the whole reason `S_t` currently runs on ILU(0)
of `A11` (PETSc's `a11` default) and grinds.

## Implementation plan

1. **Reverse the split registration order** in `snesBridge.C`, so `transport`
   is split 0 and `phi` is split 1. PETSc's Schur eliminates split 0, so this
   is what selects `S_f` over `S_t`. One-line change; verify with `-ksp_view`
   that the line reads `S = A11 - A10 inv(A00) A01` with A11 now 2000 rows
   (rule 42 -- confirm the knob moved, do not assume).

2. **Assemble the semi-implicit operator as the Schur PC matrix.** Everything
   needed exists:
   * `transport.electricalConductivity()` -- `plasmaTransport.H:863`.
   * `em.epsilon()`.
   * `blockMatrixCOO` already converts an `fvScalarMatrix` to COO triplets;
     instantiate it with `nFields = 1` for a phi-only matrix.
   * the operator is `fvm::laplacian(em.epsilon() + dt*sigma, dePotential)`,
     with `dePotential`'s homogeneous BCs -- the same field the existing
     `-laplacian(...)` blocks are built on.

   Then `PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, Sp)`. There
   is no `SchurPre` call in the tree today, so this is new but small.

3. **Solve it with AMG**: `-fieldsplit_phi_pc_type hypre` (or `gamg`). This is
   the payoff -- a scalar elliptic operator with a physics-based approximation
   is exactly AMG's home ground.

4. **THE TRAP TO AVOID, and it is the same one that caused today's SIGFPE.**
   `A_tt^-1` now becomes the INNER solve, applied on every Schur application.
   It is bigger than the phi block, so exact LU is not affordable and it must
   be iterative -- which reintroduces precisely the defect of section 25c: an
   inner solve taking a varying number of iterations makes the Schur operator
   NOT a fixed linear operator, and the outer Krylov method breaks down.

   Mitigation, and it must be in from the start rather than discovered again:
   make the inner solve a FIXED linear operator by giving it a fixed iteration
   count and skipping its convergence test --

       -fieldsplit_transport_ksp_type richardson
       -fieldsplit_transport_ksp_max_it 5
       -fieldsplit_transport_ksp_convergence_test skip

   -- so every application performs identical arithmetic. FGMRES on the outside
   remains the safety net.

5. **Measure against Tier 0** (`-fieldsplit_transport_ksp_rtol 1e-2`,
   `max_it 200`, `restart 100`), which is the cheap no-code alternative and may
   capture much of the gain on its own. Tier 3 only earns its complexity if it
   beats Tier 0 measurably.

## The reaction terms cancel EXACTLY -- charge conservation does it

**An earlier version of this note claimed the opposite and was wrong.** It said
that at an ignition front `A_tt ~ 1/dt + dP/dn`, so the reaction rate should
appear in the operator as `div((eps + sigma/(1/dt + k_eff)) grad .)`, and
offered that as a "testable refinement". The user asked the right question --
in the classical semi-implicit models the reaction source terms cancel exactly,
so why would they not here? They do. Here is why.

Any charge-conserving mechanism satisfies, IDENTICALLY as a function of `n`,

    sum_i q_i S_i(n) = 0

-- electron-impact ionisation makes an electron AND an ion, attachment moves
charge from one carrier to another, recombination removes a pair. Differentiate
that identity with respect to each `n_j`:

    sum_i q_i dS_i/dn_j = 0   for every j        <=>   q^T J_S = 0

so **`q^T` is a LEFT NULL VECTOR of the reaction Jacobian.** The Schur
complement needs `q^T A_tt^-1`, and with `A_tt = I/dt - J_S`,

    q^T (I/dt - J_S)^-1 = dt q^T (I - dt J_S)^-1
                        = dt q^T ( I + dt J_S + dt^2 J_S^2 + ... )
                        = dt q^T

because every term past the first carries a factor `q^T J_S = 0`. **The
reaction Jacobian drops out exactly, to all orders in dt.** No `k_eff`
correction exists to be made, however violent the ionisation.

This is the same cancellation the classical semi-implicit Poisson relies on --
`d(rho)/dt = -div(J) + sum_i q_i S_i` and the source sum vanishes -- and it
survives into the Schur complement untouched. It also STRENGTHENS the identity
of this note: `div((eps + dt*sigma) grad .)` is not merely the dt-dominant
limit, it is exact in the chemistry.

## What is actually approximate, then

Only the TRANSPORT part of `A_tt`. Writing `A_tt = I/dt + T - J_S` with `T` the
advection-diffusion operator, `q^T T != 0` -- transport MOVES charge between
cells rather than creating it, so it has no reason to annihilate `q^T`. The
error in the operator is therefore governed by how much of `A_tt` is transport
rather than `ddt`, i.e. by the transport CFL/diffusion numbers, and NOT by the
reaction rate at all.

**This is the sharp, falsifiable prediction to test.** The quality of
`div((eps + dt*sigma) grad .)` as a Schur preconditioner should degrade with
the transport numbers and be INDIFFERENT to `k_eff`. A sweep that raises the
ionisation rate at fixed dt and fixed mesh should leave the iteration count
alone; one that raises dt (or refines the mesh) at fixed chemistry should not.
If that holds, it is direct evidence for the identity above rather than for a
plausible-looking fit -- and one more caution against "refinements" derived
from the shape of an operator instead of from its conservation properties.
