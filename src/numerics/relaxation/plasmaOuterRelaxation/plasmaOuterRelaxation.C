#include "plasmaOuterRelaxation.H"
#include "fvMesh.H"
#include "Time.H"

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
    timeIndex_(-1),
    omega_(1.0),
    omegaMinStep_(1.0),
    omegaMaxStep_(0.0),
    nRelaxed_(0)
{}


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
    omegaMinStep_ = 1.0;
    omegaMaxStep_ = 0.0;
    nRelaxed_ = 0;
    contributed_ = false;
    nContributed_ = 0;
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

    omega_ = aitken_.omega(rj);
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
}

// ************************************************************************* //
