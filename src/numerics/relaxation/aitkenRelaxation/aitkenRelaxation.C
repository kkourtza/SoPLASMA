#include "aitkenRelaxation.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::aitkenRelaxation::aitkenRelaxation
(
    const word& name,
    const scalar omega0,
    const scalar omegaMin,
    const scalar omegaMax,
    const scalar descentLimit
)
:
    name_(name),
    rPrev_(),
    omegaPrev_(omega0),
    corr_(0),
    omegaMin_(omegaMin),
    omegaMax_(omegaMax),
    // CLAMPED into the bounds. The first corrector of every step returns
    // omega0_ directly, so an unclamped seed silently bypasses the ceiling --
    // MEASURED: with relaxOmegaMax = 0.8 the reported per-step maximum was
    // still 1.0, because the seed never went through the bounding branch.
    omega0_(max(omegaMin, min(omegaMax, omega0))),
    descentLimit_(descentLimit),
    omegaMinSeen_(omegaMax),
    omegaMaxSeen_(0)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::aitkenRelaxation::reset()
{
    rPrev_.clear();
    omegaPrev_ = omega0_;
    corr_ = 0;
    omegaMinSeen_ = omegaMax_;
    omegaMaxSeen_ = 0;
}


Foam::scalar Foam::aitkenRelaxation::omega(const scalarField& r)
{
    ++corr_;

    // First corrector of the step: no history, so there is no Aitken estimate
    // to make. Deliberately NOT 1 by default -- the first corrector after a
    // time-step change is where the coupling is furthest from its fixed point.
    if (rPrev_.empty() || rPrev_.size() != r.size())
    {
        rPrev_ = r;
        omegaPrev_ = omega0_;
        omegaMinSeen_ = min(omegaMinSeen_, omegaPrev_);
        omegaMaxSeen_ = max(omegaMaxSeen_, omegaPrev_);
        return omegaPrev_;
    }

    // dr = r_k - r_{k-1}, the change in the residual between correctors.
    // Collective sums over VALUES; identical on every rank by construction.
    const scalarField dr(r - rPrev_);

    scalar den = sumSqr(dr);
    scalar num = sumProd(rPrev_, dr);

    reduce(den, sumOp<scalar>());
    reduce(num, sumOp<scalar>());

    scalar w = omegaPrev_;

    // A vanishing denominator means the residual stopped changing between
    // correctors -- the iteration has converged (or stagnated) and Aitken has
    // nothing to estimate from. Keeping the previous factor is the safe
    // answer; dividing would be a zero-over-zero.
    if (den > VSMALL)
    {
        w = -omegaPrev_*num/den;
    }

    // RATE-LIMIT THE DESCENT before bounding. One outlier estimate must not
    // be able to slam omega onto its floor; a sustained stiffening still walks
    // it down over successive correctors. Applied before the hard bounds so
    // the floor remains a genuine last resort rather than the usual outcome.
    if (descentLimit_ > 0 && w < descentLimit_*omegaPrev_)
    {
        w = descentLimit_*omegaPrev_;
    }

    // Bound it. An unbounded Aitken factor can go negative or huge on a
    // non-smooth iteration, and either would be worse than no relaxation.
    if (!(w > omegaMin_))     // also catches NaN
    {
        w = omegaMin_;
    }
    else if (w > omegaMax_)
    {
        w = omegaMax_;
    }

    rPrev_ = r;
    omegaPrev_ = w;
    omegaMinSeen_ = min(omegaMinSeen_, w);
    omegaMaxSeen_ = max(omegaMaxSeen_, w);

    return w;
}

// ************************************************************************* //
