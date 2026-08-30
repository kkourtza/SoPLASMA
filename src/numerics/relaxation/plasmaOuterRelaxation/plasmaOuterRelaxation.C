#include "plasmaOuterRelaxation.H"
#include "fvMesh.H"
#include "Time.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"

namespace Foam
{
    defineTypeNameAndDebug(plasmaOuterRelaxation, 0);
}

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::plasmaOuterRelaxation::plasmaOuterRelaxation
(
    const fvMesh& mesh,
    const dictionary& oc,
    const bool defaultActive
)
:
    regIOobject
    (
        IOobject
        (
            "plasmaOuterRelaxation",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    active_(oc.getOrDefault<Switch>("adaptiveRelaxation", defaultActive)),
    fields_(),
    prev_(),
    pending_(),
    contributed_(),
    nContributed_(0),
    // Read ONCE, here. Everything else -- the estimator below, the start-up
    // banner, any future reporter -- uses these resolved values.
    omegaStart_(oc.getOrDefault<scalar>("relaxOmegaStart", 0.8)),
    omegaFloor_(oc.getOrDefault<scalar>("relaxOmegaMin", 0.05)),
    descentLimit_(oc.getOrDefault<scalar>("relaxOmegaDescentLimit", 0.0)),

    aitken_
    (
        "joint",

        // THE SINGLE USER-FACING KNOB: relaxOmegaStart, the factor each time
        // step BEGINS with. From the second corrector onwards Aitken computes
        // omega itself, so this is the one number a user ever needs to set.
        //
        // Default 0.8, not 1.0. MEASURED 2026-08-21 on the LMEA streamer case,
        // 15 steps each: seeding at 0.8 reached t = 5.99e-11 against 5.21e-11
        // from an unrelaxed (1.0) seed -- about 15%, and level with a
        // hand-tuned FIXED factor of 0.8. An unrelaxed first corrector injects
        // a disturbance the rest of the step then has to fight.
        //
        // 0.8 IS A STARTING POINT, NOT A UNIVERSAL CONSTANT, and it is
        // CASE-SPECIFIC: the right value depends on the coupling strength,
        // which varies with mechanism, mesh and deltaT. **Lower it** (0.5-0.7)
        // if the run reports frequent `min(omega)` floor hits or
        // `deltaT GOVERNED` -- both mean the step starts too loose and
        // correctors are spent recovering.
        //
        // CAVEAT: the figure comes from SINGLE runs with a 10-15% spread and
        // no repeats. Indicative, not established.
        omegaStart_,

        // Floor: bounds how hard Aitken may damp before the deltaT governor
        // takes over instead. Advanced; rarely needs changing.
        omegaFloor_,

        // Ceiling FIXED at 1.0 and deliberately NOT exposed. Capping it was
        // measured to cost ~10% and buy nothing, and a second knob that mostly
        // duplicates the first is a trap.
        1.0,

        // Descent rate limit. OFF by default -- measured HARMFUL: it kept
        // omega off its floor, which also suppressed the deltaT governor's
        // trigger, and the run grew deltaT into a SIGFPE at step 11 where the
        // unlimited run survived 15 steps.
        descentLimit_
    ),
    // DEFAULT STAYS `aitken` even though Anderson measured 30x better on the
    // unit bed (testAitken CASE 10: 91 correctors -> 3). That result is on
    // synthetic problems; it has not yet been shown on a case with real
    // chemistry, a Poisson solve and spatially varying coefficients. A default
    // is a promise, and this one has not been earned yet.
    scheme_(oc.getOrDefault<word>("outerScheme", "aitken")),
    // DEPTH 4, not the unit bed's 2. The synthetic two-field cycle has ONE
    // oscillatory pair, so a 2-dimensional Krylov space spans it exactly and
    // deeper windows only span more of the density clamp's active-set changes.
    // A real residual carries many more modes. MEASURED on the fast LMEA bed,
    // one variable, against joint Aitken at 49 steps / 477 correctors:
    //     m = 2   51 steps, 537 correctors   -- WORSE than Aitken
    //     m = 4   39 steps, 411 correctors   -- -20% steps, -24% wall clock
    //     m = 8   49 steps, 532 correctors   -- level with Aitken
    // Shipping 2 as the default would hand anyone who enables Anderson the one
    // setting measured worse than not enabling it.
    //
    // The optimum is NON-MONOTONE and problem-dependent -- the same pattern
    // appears on the clamped unit bed (m = 2,3,4,5 -> 8, 11, 27, 8 correctors)
    // -- and 4 rests on ONE case. Sweep it before trusting it on a new problem.
    andersonDepth_(oc.getOrDefault<label>("andersonDepth", 4)),
    andersonBeta_(oc.getOrDefault<scalar>("andersonBeta", 1.0)),
    // Sticky window. 5 is a starting point, not a measured optimum.
    andersonStickySteps_(oc.getOrDefault<label>("andersonStickySteps", 5)),
    lastEscalatedIndex_(-1000000),
    anderson_("joint", andersonDepth_, andersonBeta_),
    timeIndex_(-1),
    omega_(1.0),
    omegaMinStep_(1.0),
    omegaMaxStep_(0.0),
    nRelaxed_(0),
    rNormPrev_(0),
    contractionMaxStep_(0),
    escalatedTimeIndex_(-1)
{
    if (scheme_ != "aitken" && scheme_ != "anderson" && scheme_ != "adaptive")
    {
        FatalErrorInFunction
            << "outerCoupling/outerScheme = `" << scheme_
            << "` is not recognised." << nl
            << "    Use aitken | anderson | adaptive." << exit(FatalError);
    }

    if (!active_) return;

    if (scheme_ == "adaptive")
    {
        WarningInFunction
            << "outerScheme `adaptive` MEASURED WORSE than both alternatives"
               " on the only case tried." << nl << nl
            << "    Aitken normally, escalating to Anderson (depth "
            << andersonDepth_ << ") on a RETRIED step at the" << nl
            << "    same deltaT. The mechanism works -- it reached the target"
               " in the FEWEST steps" << nl
            << "    of any scheme -- but each escalation is preceded by a"
               " complete FAILED Aitken" << nl
            << "    attempt at maxCorrectors, and that detection cost exceeds"
               " the saving." << nl << nl
            << "    Measured, fast LMEA bed, total correctors INCLUDING"
               " discarded attempts:" << nl
            << "        aitken     49 steps,  597 correctors, 76.0 s" << nl
            << "        anderson   39 steps,  431 correctors, 57.7 s" << nl
            << "        adaptive   31 steps,  624 correctors, 92.9 s" << nl << nl
            << "    The premise -- that Aitken is cheaper on the easy steps and"
               " Anderson should be" << nl
            << "    spent only where needed -- is FALSE on that bed: Aitken is"
               " dearer overall"  << nl
            << "    (597 vs 431) because it fails four times as often. Prefer"
               " `anderson`." << nl
            << "    Kept because a cheaper trigger (escalate mid-step, or stay"
               " escalated for" << nl
            << "    several steps) might change this, and because"
               " `retryStep` is required either way."
            << endl;
    }

    if (scheme_ == "anderson")
    {
        if (andersonDepth_ < 2)
        {
            FatalErrorInFunction
                << "outerCoupling/andersonDepth = " << andersonDepth_
                << " cannot converge this coupling." << nl
                << "    The n_e/nEps_e oscillation is a COMPLEX CONJUGATE"
                   " PAIR, and a real Krylov" << nl
                << "    method cannot represent half of one. Measured"
                   " (testAitken CASE 10): m = 1" << nl
                << "    does not converge at all where m = 2 takes 3"
                   " correctors. Use 2 or more." << exit(FatalError);
        }

        Info<< "plasmaOuterRelaxation: ANDERSON acceleration, depth "
            << andersonDepth_ << ", beta " << andersonBeta_ << nl
            << "    The deltaT coupling-margin governor is DISABLED under"
               " Anderson, deliberately:" << nl
            << "    Anderson does not damp, so it has no omega, and it"
               " converged at every" << nl
            << "    stiffness swept on both the smooth and the clamped unit"
               " bed. It fails as a" << nl
            << "    CLIFF, and `onNonConvergence retryStep` is the mechanism"
               " for that -- make" << nl
            << "    sure it is on. `rho [contraction]` remains the reported"
               " diagnostic." << endl;
    }
}


// * * * * * * * * * * * * * * * *  Selectors  * * * * * * * * * * * * * * * //

Foam::plasmaOuterRelaxation& Foam::plasmaOuterRelaxation::New
(
    const fvMesh& mesh,
    const dictionary& oc,
    const bool defaultActive
)
{
    plasmaOuterRelaxation* p = lookup(mesh);
    if (p) return *p;

    plasmaOuterRelaxation* fresh =
        new plasmaOuterRelaxation(mesh, oc, defaultActive);
    fresh->store();                     // the registry takes ownership
    return *fresh;
}


Foam::plasmaOuterRelaxation* Foam::plasmaOuterRelaxation::lookup
(
    const fvMesh& mesh
)
{
    return mesh.getObjectPtr<plasmaOuterRelaxation>("plasmaOuterRelaxation");
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::plasmaOuterRelaxation::find(const word& name) const
{
    forAll(fields_, i)
    {
        if (fields_[i].name() == name) return i;
    }
    return -1;
}


void Foam::plasmaOuterRelaxation::enrol(volScalarField& f)
{
    if (!active_) return;
    if (find(f.name()) >= 0) return;        // already enrolled

    const label i = fields_.size();
    fields_.resize(i + 1);
    fields_.set(i, &f);

    prev_.resize(i + 1);
    pending_.resize(i + 1);
    contributed_.resize(i + 1, false);

    // Seed the previous iterate with the INITIAL CONDITION. enrol() is called
    // before the time loop, so this is the only point at which the pre-solve
    // value is available without a before-solve hook.
    prev_[i] = f.primitiveField();
    pending_[i] = f.primitiveField();

    Info<< "plasmaOuterRelaxation: enrolled `" << f.name()
        << "` (" << fields_.size() << " field(s) in the joint residual)"
        << endl;
}


void Foam::plasmaOuterRelaxation::newStepCheck()
{
    if (timeIndex_ == mesh_.time().timeIndex()) return;

    timeIndex_ = mesh_.time().timeIndex();

    // A residual may never be carried across a time step: it would relate
    // iterates of two different problems. prev_ is NOT re-seeded -- it already
    // holds the accepted solution of the previous step, which is exactly where
    // this step's Picard iteration starts.
    aitken_.reset();
    anderson_.reset();
    omegaMinStep_ = 1.0;
    omegaMaxStep_ = 0.0;
    nRelaxed_ = 0;
    contributed_ = false;
    nContributed_ = 0;

    // A residual norm may no more be carried across a step than a residual.
    rNormPrev_ = 0;
    contractionMaxStep_ = 0;
}


void Foam::plasmaOuterRelaxation::contribute(volScalarField& f)
{
    if (!active_) return;

    newStepCheck();

    const label i = find(f.name());
    if (i < 0) return;                      // not taking part

    if (contributed_[i])
    {
        FatalErrorInFunction
            << "field `" << f.name() << "` contributed twice within one outer"
            << " corrector." << nl
            << "    The joint Aitken factor is assembled from exactly one"
            << " residual per enrolled field" << nl
            << "    per corrector; a second contribution means the call is in"
            << " the wrong place." << nl
            << exit(FatalError);
    }

    pending_[i] = f.primitiveField();
    contributed_[i] = true;
    ++nContributed_;

    if (nContributed_ == fields_.size())
    {
        applyJoint();
    }
}


Foam::word Foam::plasmaOuterRelaxation::effectiveScheme() const
{
    if (scheme_ != "adaptive") return scheme_;

    const label ti = mesh_.time().timeIndex();

    // Escalated for THIS step, or still inside the sticky window opened by a
    // recent one. The window is what makes the detection cost amortise: without
    // it every hard step pays a full failed Aitken attempt to be recognised.
    if (escalatedTimeIndex_ == ti) return word("anderson");

    if (andersonStickySteps_ > 0
     && ti - lastEscalatedIndex_ <= andersonStickySteps_)
    {
        return word("anderson");
    }

    return word("aitken");
}


bool Foam::plasmaOuterRelaxation::escalateForRetry()
{
    if (!active_ || scheme_ != "adaptive") return false;

    const label ti = mesh_.time().timeIndex();

    // NOTHING TO ESCALATE TO. If Anderson is ALREADY driving this step -- this
    // step was escalated, or a recent one was and the sticky window still
    // covers it -- then the failure happened UNDER Anderson, and re-running at
    // the same deltaT with the same solver cannot do better. Let the retry
    // ladder shorten the step, which is the mechanism that always works.
    //
    // MEASURED before this guard existed: with andersonStickySteps 1000 the
    // first escalation was at step 19 and the window then covered the whole
    // run, yet steps 24, 28, 29 and 30 "escalated" again -- four guaranteed
    // wasted retries at maxCorrectors each. It also renewed the window every
    // time, which is why windows of 5, 20 and 1000 produced BIT-IDENTICAL runs
    // (32 steps, 521 correctors, 76.2 s) and looked like a knob that did
    // nothing.
    if (effectiveScheme() == "anderson") return false;

    escalatedTimeIndex_ = ti;
    lastEscalatedIndex_ = ti;

    Info<< "  outer loop did not converge -- ESCALATING to Anderson (depth "
        << andersonDepth_ << ") and re-running" << nl
        << "    at the SAME deltaT. If it fails again the step is shortened as"
           " usual." << nl;

    if (andersonStickySteps_ > 0)
    {
        Info<< "    Staying on Anderson for the next " << andersonStickySteps_
            << " step(s): stiffness is temporally correlated, and" << nl
            << "    re-detecting it costs a full failed attempt each time."
            << endl;
    }
    else
    {
        Info<< endl;
    }

    return true;
}


void Foam::plasmaOuterRelaxation::applyJoint()
{
    // Concatenate the residuals into ONE vector. The unstable mode lives in
    // the coupling between the fields, so it is invisible to any single
    // field's own residual -- see the class comment for the measurement.
    label nTot = 0;
    forAll(fields_, i) { nTot += pending_[i].size(); }

    scalarField rj(nTot);
    label off = 0;
    forAll(fields_, i)
    {
        const scalarField& p = pending_[i];
        const scalarField& q = prev_[i];
        forAll(p, c) { rj[off + c] = p[c] - q[c]; }
        off += p.size();
    }

    // CONTRACTION, measured before anything else touches omega.
    //
    // rho = ||r_k|| / ||r_{k-1}||: below 1 the outer loop is contracting, at
    // or above 1 it is not. Unlike omega this is defined for ANY relaxation
    // scheme -- including Anderson, which has no omega -- so it is the signal
    // a time-step governor can keep using across a change of scheme.
    //
    // Collective, over VALUES, exactly as the Aitken reductions are, so every
    // rank sees the same number.
    scalar rNorm = sumSqr(rj);
    reduce(rNorm, sumOp<scalar>());
    rNorm = Foam::sqrt(rNorm);

    if (rNormPrev_ > VSMALL)
    {
        // The WORST ratio over the step, because one corrector that fails to
        // contract is what a governor needs to see; an average over a step
        // that recovered would hide it.
        contractionMaxStep_ = max(contractionMaxStep_, rNorm/rNormPrev_);
    }
    rNormPrev_ = rNorm;

    // NO RELAXATION ON THE FINAL OUTER CORRECTOR.
    //
    // The accepted iterate of a time step must be what the equations actually
    // SOLVED, not a blend of the solve and the previous iterate. Relaxation is
    // a device for getting the Picard iteration to its fixed point; once the
    // step is being accepted it must not bias the answer.
    //
    // The fixed-factor path has always done this -- a case writes
    // `n_eFinal 1.0` and fvMatrix::relax() picks it up through
    // psi_.select(isFinalIteration()). This path did not, and the cost was
    // MEASURED on the order study: with the outer loop exiting in 2-3
    // correctors, the residual (1 - omega) bias does not vanish at the right
    // rate as deltaT falls and the observed temporal order collapsed to
    // p = 0.53 -- for the SPECIES as well as the energy, on a bed where the
    // same arm without LMEA measures p = 1.95.
    //
    // At a genuinely converged fixed point r -> 0 and this changes nothing;
    // it matters exactly when the loop stops early, which is the normal case.
    // THE SAME RULE BINDS ANDERSON. Its update is also a blend -- of the
    // history rather than of two iterates -- so accepting it as the step's
    // answer biases the result exactly as a relaxed iterate does. On the final
    // corrector both schemes take the solve unchanged.
    const bool finalIter = mesh_.data().isFinalIteration();

    if (effectiveScheme() == "anderson" && !finalIter)
    {
        // Anderson works on the WHOLE joint vector, so gather prev into one
        // contiguous x, let it update in place, and scatter back. rj is
        // already the joint residual r = pending - prev evaluated at x.
        scalarField x(nTot);
        label o = 0;
        forAll(fields_, i)
        {
            const scalarField& q = prev_[i];
            forAll(q, c) { x[o + c] = q[c]; }
            o += q.size();
        }

        anderson_.correct(x, rj);

        o = 0;
        forAll(fields_, i)
        {
            volScalarField& f = fields_[i];
            scalarField& fi = f.primitiveFieldRef();
            forAll(fi, c) { fi[c] = x[o + c]; }
            o += fi.size();

            f.correctBoundaryConditions();
            prev_[i] = f.primitiveField();
        }

        // omega has no meaning here; report 1 so the margin reads "nothing
        // damped" rather than a stale number from a scheme that is not running.
        omega_ = 1.0;
        omegaMinStep_ = min(omegaMinStep_, omega_);
        omegaMaxStep_ = max(omegaMaxStep_, omega_);
        ++nRelaxed_;

        contributed_ = false;
        nContributed_ = 0;
        return;
    }

    omega_ = finalIter ? 1.0 : aitken_.omega(rj);

    omegaMinStep_ = min(omegaMinStep_, omega_);
    omegaMaxStep_ = max(omegaMaxStep_, omega_);
    ++nRelaxed_;

    // Apply the SAME factor to every field: x <- prev + omega (pending - prev)
    forAll(fields_, i)
    {
        volScalarField& f = fields_[i];
        scalarField& fi = f.primitiveFieldRef();
        const scalarField& p = pending_[i];
        const scalarField& q = prev_[i];

        forAll(fi, c)
        {
            fi[c] = q[c] + omega_*(p[c] - q[c]);
        }

        f.correctBoundaryConditions();

        // The relaxed value is the accepted iterate this corrector hands on.
        prev_[i] = f.primitiveField();
    }

    contributed_ = false;
    nContributed_ = 0;
}


void Foam::plasmaOuterRelaxation::discardStep()
{
    if (!active_) return;

    // The attempt is being thrown away, so its iterate history is meaningless.
    // Re-seed from the fields as they stand AFTER the restore, and force the
    // next contribute() to treat this as a fresh step.
    aitken_.reset();

    // THE STEP-DISCARD INVARIANT applies to Anderson's history exactly as to
    // Aitken's: a difference taken across a discarded attempt describes a step
    // that no longer exists.
    anderson_.reset();

    forAll(fields_, i)
    {
        prev_[i] = fields_[i].primitiveField();
        pending_[i] = prev_[i];
    }
    contributed_ = false;
    nContributed_ = 0;
    omegaMinStep_ = 1.0;
    omegaMaxStep_ = 0.0;
    nRelaxed_ = 0;
    timeIndex_ = -1;

    // THE STEP-DISCARD INVARIANT. These are per-pass accumulators like every
    // other diagnostic here: without this the discarded attempt's contraction
    // would be reported alongside, and eventually govern, the kept one.
    rNormPrev_ = 0;
    contractionMaxStep_ = 0;
}

// ************************************************************************* //
