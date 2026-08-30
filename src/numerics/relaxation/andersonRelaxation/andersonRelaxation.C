#include "andersonRelaxation.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::andersonRelaxation::andersonRelaxation
(
    const word& name,
    const label m,
    const scalar beta,
    const scalar lambda,
    const scalar maxStepRatio
)
:
    name_(name),
    m_(max(m, label(0))),
    beta_(beta),
    lambda_(lambda),
    maxStepRatio_(maxStepRatio),
    dX_(),
    dR_(),
    xPrev_(),
    rPrev_(),
    havePrev_(false),
    corr_(0),
    nFallbackStep_(0),
    nFallbackSolve_(0),
    mUsedMax_(0)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::andersonRelaxation::clearHistory()
{
    dX_.clear();
    dR_.clear();
    havePrev_ = false;
}


bool Foam::andersonRelaxation::solveDense
(
    List<scalarField>& A,
    scalarField& b,
    scalarField& gamma
)
{
    const label n = b.size();
    gamma.setSize(n, 0.0);

    if (n == 0) return true;

    // Gauss elimination with PARTIAL PIVOTING. Without pivoting the normal
    // matrix of a nearly rank-deficient dR -- which is what convergence looks
    // like -- produces a tiny leading pivot and a gamma that is numerical
    // noise amplified by 1/eps.
    for (label k = 0; k < n; ++k)
    {
        label piv = k;
        scalar best = mag(A[k][k]);
        for (label i = k + 1; i < n; ++i)
        {
            if (mag(A[i][k]) > best) { best = mag(A[i][k]); piv = i; }
        }

        // Singular to working precision. Reported as a failure rather than
        // patched: the caller has a safe fallback and should take it.
        if (best < SMALL) return false;

        if (piv != k)
        {
            Swap(A[piv], A[k]);
            Swap(b[piv], b[k]);
        }

        for (label i = k + 1; i < n; ++i)
        {
            const scalar f = A[i][k]/A[k][k];
            if (f == 0) continue;
            for (label j = k; j < n; ++j) { A[i][j] -= f*A[k][j]; }
            b[i] -= f*b[k];
        }
    }

    for (label i = n - 1; i >= 0; --i)
    {
        scalar s = b[i];
        for (label j = i + 1; j < n; ++j) { s -= A[i][j]*gamma[j]; }
        gamma[i] = s/A[i][i];
    }

    return true;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::andersonRelaxation::reset()
{
    clearHistory();
    corr_ = 0;
    nFallbackStep_ = 0;
    nFallbackSolve_ = 0;
    mUsedMax_ = 0;
}


void Foam::andersonRelaxation::correct(scalarField& x, const scalarField& r)
{
    ++corr_;

    // A size change means the field was reallocated under us; the stored
    // differences describe a different vector space.
    if (havePrev_ && (xPrev_.size() != x.size() || rPrev_.size() != r.size()))
    {
        clearHistory();
    }

    // The plain relaxed step, which is both the m = 0 answer and the fallback
    // whenever the accelerated one cannot be trusted.
    const scalarField plain(beta_*r);

    if (m_ > 0 && havePrev_)
    {
        dX_.append(scalarField(x - xPrev_));
        dR_.append(scalarField(r - rPrev_));

        while (dX_.size() > m_)
        {
            // Oldest first out. Anderson's window is a moving one; keeping
            // stale differences describes an iteration matrix that no longer
            // applies.
            for (label i = 1; i < dX_.size(); ++i)
            {
                dX_[i - 1] = dX_[i];
                dR_[i - 1] = dR_[i];
            }
            dX_.setSize(dX_.size() - 1);
            dR_.setSize(dR_.size() - 1);
        }
    }

    xPrev_ = x;
    rPrev_ = r;
    havePrev_ = true;

    const label mUsed = dR_.size();
    mUsedMax_ = max(mUsedMax_, mUsed);

    if (mUsed == 0)
    {
        x += plain;
        return;
    }

    // Normal equations of  min || r - dR gamma ||.  All inner products are
    // collective sums over VALUES, so A and b -- and therefore gamma -- are
    // bit-identical on every rank.
    List<scalarField> A(mUsed, scalarField(mUsed, 0.0));
    scalarField b(mUsed, 0.0);

    for (label i = 0; i < mUsed; ++i)
    {
        for (label j = i; j < mUsed; ++j)
        {
            scalar aij = sumProd(dR_[i], dR_[j]);
            reduce(aij, sumOp<scalar>());
            A[i][j] = aij;
            A[j][i] = aij;
        }
        scalar bi = sumProd(dR_[i], r);
        reduce(bi, sumOp<scalar>());
        b[i] = bi;
    }

    // TIKHONOV, scaled to the trace. dR goes rank-deficient exactly as the
    // iteration converges, so an ABSOLUTE lambda would be far too small for a
    // density field at 1e19 and far too large for a normalised one.
    scalar trace = 0;
    for (label i = 0; i < mUsed; ++i) { trace += A[i][i]; }

    if (trace <= VSMALL)
    {
        // The residual stopped changing between correctors: converged or
        // stagnated, and there is nothing to extrapolate from.
        x += plain;
        return;
    }

    const scalar reg = lambda_*trace/scalar(mUsed);
    for (label i = 0; i < mUsed; ++i) { A[i][i] += reg; }

    scalarField gamma;
    if (!solveDense(A, b, gamma))
    {
        ++nFallbackSolve_;
        clearHistory();
        x += plain;
        return;
    }

    // x_{k+1} = x_k + beta r - (dX + beta dR) gamma
    scalarField step(plain);
    for (label i = 0; i < mUsed; ++i)
    {
        step -= gamma[i]*(dX_[i] + beta_*dR_[i]);
    }

    // STEP LIMITER. An Anderson step is not a descent step, and this solver's
    // iteration is NOT smooth -- the density clamp is a projection, and
    // testAitken CASE 5 records that extrapolating across it is invalid.
    // Compared in the 2-norm, on the same collective reduction as everything
    // else.
    scalar nStep = sumSqr(step);
    scalar nPlain = sumSqr(plain);
    reduce(nStep, sumOp<scalar>());
    reduce(nPlain, sumOp<scalar>());

    if (nStep > sqr(maxStepRatio_)*nPlain + VSMALL)
    {
        ++nFallbackStep_;
        clearHistory();
        x += plain;
        return;
    }

    x += step;
}


// ************************************************************************* //
