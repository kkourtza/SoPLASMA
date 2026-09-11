/*---------------------------------------------------------------------------*\
Class
    Foam::andersonRelaxation

Description
    Anderson acceleration of a Picard (fixed-point) outer loop, after Anderson,
    J. ACM 12 (1965) 547-560, in the difference form of Walker & Ni, SIAM J.
    Numer. Anal. 49 (2011) 1715-1735.

    WHY THIS EXISTS, AND WHY AITKEN IS NOT ENOUGH. The n_e / nEps_e coupling
    enters a period-2 limit cycle above a critical deltaT. aitkenRelaxation
    damps it with a single SCALAR factor, which is the best a scalar can do:
    measured on the synthetic two-field problem of testAitken CASE 7, at loop
    gain 1.43, joint Aitken needs 91 correctors, per-field Aitken 191 with
    omega pinned on its floor, and relaxing both fields separately does not
    converge at all. Damping buys stability by making every mode slower,
    including the ones that were already fine.

    Anderson instead builds the next iterate from a least-squares combination
    of the last m residuals, so it can CANCEL the offending mode rather than
    merely shrink it. On a linear problem with m at least the number of active
    modes it is equivalent to GMRES and terminates in m steps; the two-field
    cycle has ONE active mode, so even m = 1 should see it.

    THE METHOD. With residual r_k = G(x_k) - x_k and the last m differences

        dR = [ r_{k-m+1}-r_{k-m} ... r_k - r_{k-1} ]
        dX = [ x_{k-m+1}-x_{k-m} ... x_k - x_{k-1} ]

    solve the small least-squares problem

        gamma = argmin || r_k - dR gamma ||_2

    and take

        x_{k+1} = x_k + beta r_k - (dX + beta dR) gamma.

    With no history this is exactly x_{k+1} = x_k + beta r_k, i.e. plain
    under-relaxation, which is what makes the first corrector of a step safe.

    THE SAFEGUARDS, AND WHY EACH IS THERE.

    * TIKHONOV REGULARISATION. dR becomes rank-deficient as the iteration
      converges (successive residuals become parallel), and the normal
      equations are then singular. lambda is scaled to the trace, so it is
      dimensionless in the problem's own units rather than an absolute number
      that is wrong for n_e ~ 1e19.

    * STEP LIMITER. An Anderson step is not a descent step and can be
      arbitrarily long on a non-smooth iteration -- and this solver HAS a
      non-smooth iteration: the density clamp is a projection, and
      testAitken CASE 5 already records that extrapolation across it is not
      valid. If the proposed step exceeds maxStepRatio times the plain
      relaxed step, the history is dropped and the plain step taken. A
      method that can silently take a 100x step has no place near a
      production solver.

    * RESTART EVERY TIME STEP. reset() clears the history, because the fixed
      point moves between steps and a difference spanning two of them
      describes neither.

Note
    Every inner product is a collective reduction over VALUES, so the m x m
    system and hence gamma are identical on every rank by construction, exactly
    as in aitkenRelaxation. Cell indices are never communicated.

SourceFiles
    andersonRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef andersonRelaxation_H
#define andersonRelaxation_H

#include "scalarField.H"
#include "word.H"
#include "DynamicList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class andersonRelaxation Declaration
\*---------------------------------------------------------------------------*/

class andersonRelaxation
{
    // Private data

        //- Field (or concatenated field group) this instance accelerates
        word name_;

        //- History depth. m = 0 degenerates to plain relaxation by beta.
        label m_;

        //- Mixing factor on the residual. beta = 1 is pure Anderson; below 1
        //  it under-relaxes the part of the step Anderson does not cancel.
        scalar beta_;

        //- Tikhonov parameter, RELATIVE to the trace of the normal matrix.
        scalar lambda_;

        //- Largest tolerated ratio of the Anderson step to the plain step
        //  beta*r. Above it the history is dropped and the plain step taken.
        scalar maxStepRatio_;

        //- Difference columns, newest last, at most m_ of each
        DynamicList<scalarField> dX_;
        DynamicList<scalarField> dR_;

        //- Previous iterate and residual, for forming the next difference
        scalarField xPrev_;
        scalarField rPrev_;
        bool havePrev_;

    // Diagnostics, per step

        //- Correctors taken since reset()
        label corr_;

        //- How many of them fell back to the plain step, and why
        label nFallbackStep_;
        label nFallbackSolve_;

        //- Largest history depth actually used
        label mUsedMax_;


    // Private Member Functions

        //- Solve the small dense system A gamma = b by Gauss elimination with
        //  partial pivoting. Returns false if it is singular to working
        //  precision, in which case the caller falls back.
        static bool solveDense
        (
            List<scalarField>& A,
            scalarField& b,
            scalarField& gamma
        );

        //- Drop all history. Used by reset() and by the step limiter.
        void clearHistory();


public:

    // Constructors

        andersonRelaxation
        (
            const word& name,
            const label m = 3,
            const scalar beta = 1.0,
            const scalar lambda = 1e-10,
            const scalar maxStepRatio = 10.0
        );


    // Member Functions

        //- Apply the accelerated update to x IN PLACE, given the residual
        //  r = G(x) - x evaluated at the current x.
        void correct(scalarField& x, const scalarField& r);

        //- Clear the history. MUST be called at the start of every time step:
        //  the fixed point moves, and a difference spanning two steps
        //  describes neither.
        void reset();

        //- Diagnostics for the step just taken
        label corr() const { return corr_; }
        label nFallbackStep() const { return nFallbackStep_; }
        label nFallbackSolve() const { return nFallbackSolve_; }
        label mUsedMax() const { return mUsedMax_; }
        const word& name() const { return name_; }
};


} // End namespace Foam

#endif

// ************************************************************************* //
