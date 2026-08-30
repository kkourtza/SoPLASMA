/*---------------------------------------------------------------------------*\
Class
    Foam::plasmaOuterRelaxation

Description
    JOINT Aitken under-relaxation of the fields that take part in the PIMPLE
    outer (Picard) loop, held on the mesh registry so the models that solve
    those fields do not have to know about each other.

    WHY JOINT, AND NOT ONE FACTOR PER FIELD. The instability this exists to
    damp is a period-2 cycle BETWEEN the electron density and the electron
    energy: the loop gain is a property of the PAIR, and neither field's own
    residual contains it. MEASURED in testAitken on the two-field model
    problem (composite gain -1.43):

        only n_e relaxed .................. 191 correctors, omega floors at 0.05
        both relaxed, INDEPENDENT factors . never converges, omega floors
        JOINT factor ...................... 91 correctors, omega 0.371, no floor

    Independent per-field factors are not merely weaker -- they are WORSE than
    relaxing a single field, because each estimate is chasing a map the other
    field keeps moving. Aitken for a vector fixed point takes the whole
    residual vector (Kuettler & Wall, Comput. Mech. 43 (2008) 61-72), which is
    what this class assembles.

    OWNERSHIP OF THE PREVIOUS ITERATE. The coordinator keeps its own copy of
    each enrolled field's last ACCEPTED iterate, rather than relying on
    GeometricField::storePrevIter(). That makes it independent of the order in
    which the models solve, needs no before-solve hook, and cannot be
    invalidated by another part of the solver storing its own prevIter.

    USAGE
      - enrol(f) once, at model construction, BEFORE the time loop: this also
        seeds the previous iterate with the initial condition.
      - contribute(f) once per outer corrector, immediately after that field's
        solve. When every enrolled field has contributed, the joint factor is
        computed and applied to all of them.
      - discardStep() when a step is thrown away and retried, so a stale
        iterate from the abandoned attempt is not carried into the retry.

SourceFiles
    plasmaOuterRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaOuterRelaxation_H
#define plasmaOuterRelaxation_H

#include "aitkenRelaxation.H"
#include "andersonRelaxation.H"
#include "volFields.H"
#include "regIOobject.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class plasmaOuterRelaxation Declaration
\*---------------------------------------------------------------------------*/

class plasmaOuterRelaxation
:
    public regIOobject
{
    // Private data

        const fvMesh& mesh_;

        //- Off unless the case asks for it. Every validated result in this
        //  project predates this class.
        bool active_;

        //- Fields taking part, in enrolment order
        UPtrList<volScalarField> fields_;

        //- Last ACCEPTED iterate of each field
        List<scalarField> prev_;

        //- This corrector's solved (pre-relaxation) value of each field
        List<scalarField> pending_;

        //- Which fields have contributed in the current corrector
        boolList contributed_;
        label nContributed_;

        //- The RESOLVED settings. Declared BEFORE aitken_ so the estimator is
        //  constructed FROM them: the dictionary keys are then read exactly
        //  ONCE in the whole program, and anything that wants to report the
        //  configuration asks for these rather than re-reading with a default
        //  of its own -- which is how a half-configured feature happens.
        scalar omegaStart_;
        scalar omegaFloor_;
        scalar descentLimit_;

        //- The joint estimator
        aitkenRelaxation aitken_;

        //- `aitken` (default), `anderson`, or `adaptive`.
        //
        //  ADAPTIVE runs Aitken normally and ESCALATES TO ANDERSON on a
        //  retried step, at the SAME deltaT, before the retry ladder is
        //  allowed to halve.
        //
        //  Why that trigger and not a predictive one: testAitken CASE 16
        //  measured NO LEAD in any available signal -- omega, rho and the
        //  corrector count all look the same on the last converging step and
        //  the first failing one -- so a step cannot be known to be hard until
        //  it fails. Under `retryStep` a failed step is already discarded and
        //  re-run, so the FAILURE ITSELF is the trigger and no prediction is
        //  needed.
        //
        //  Why it should pay: on the fast LMEA bed Anderson m = 4 rejected 3
        //  steps where Aitken rejected 12, at the same tolerance. Anderson's
        //  measured value is rescuing steps Aitken cannot converge, not
        //  running the easy ones faster -- correctors per step were 9.73 for
        //  Aitken against ~10.5 for every Anderson arm. Adaptive spends
        //  Anderson only where it demonstrably helps.
        //
        //  Switching at a STEP BOUNDARY on a discarded step is sound; switching
        //  mid-step would not be, because Anderson's differences would describe
        //  an Aitken-relaxed map rather than the one it is extrapolating.
        //  discardStep() clears the history in any case.
        word scheme_;

        //- Anderson history depth. 4 by default -- measured, not guessed, and
        //  NOT the unit bed's answer; see the constructor for the numbers.
        //  Legacy note on why 2 is the floor rather than the default:
        //  the n_e/nEps_e oscillation is a COMPLEX CONJUGATE PAIR, so a real
        //  Krylov method needs a 2-dimensional subspace and m = 1 cannot
        //  converge it at all (testAitken CASE 10). Deeper windows span more
        //  active-set changes of the density clamp, over which the linear
        //  model does not hold, and measured WORSE: on the clamped bed m = 2
        //  took 8 correctors, m = 4 took 27 with six solve fallbacks.
        label andersonDepth_;

        //- Mixing factor on the residual; 1 is pure Anderson.
        scalar andersonBeta_;

        andersonRelaxation anderson_;

        //- Time index the history was last reset at
        label timeIndex_;

        //- Factor used by the last corrector, and the smallest this step
        scalar omega_;
        scalar omegaMinStep_;
        scalar omegaMaxStep_;

        //- Correctors relaxed since the last reset
        label nRelaxed_;

        //- 2-norm of the joint residual at the previous corrector, and the
        //  WORST contraction ratio ||r_k||/||r_{k-1}|| seen this step.
        //
        //  WHY, ALONGSIDE omega. omega answers "how hard did the relaxation
        //  have to damp", which is meaningful only for a method that damps.
        //  Anderson acceleration does not damp -- it cancels the offending
        //  mode -- so it has no omega at all, and a time-step governor keyed
        //  on omega would silently stop governing the moment it were enabled.
        //
        //  The contraction ratio asks the question the governor actually
        //  wants: IS THE OUTER LOOP CONTRACTING, AND HOW FAST. It is defined
        //  for any method, including none, so it can carry the governor across
        //  a change of relaxation scheme.
        //
        //  It is NOT a drop-in replacement for omega and is currently reported
        //  only. rho is measured on the RELAXED iteration, so a small rho with
        //  a small omega means "contracting only because it is being damped
        //  hard" -- the pair says more than either alone, and which of them
        //  should drive the governor is a question for measurement, not for
        //  assertion.
        scalar rNormPrev_;
        scalar contractionMaxStep_;

        //- How many FURTHER steps stay on Anderson after an escalation.
        //
        //  The escalation trigger is a failed Aitken step, which costs a full
        //  maxCorrectors attempt to detect. MEASURED on the fast LMEA bed: 12
        //  escalations, 11 of which converged at the SAME deltaT that Aitken
        //  had just failed at -- a 92% success rate -- but 260 correctors were
        //  burned detecting them, making the scheme a net loss on wall clock.
        //
        //  Stiffness is temporally correlated, so a step that needed Anderson
        //  is a good predictor that the next few will too. Staying escalated
        //  amortises one detection over several steps instead of re-paying it.
        //  0 restores the pay-per-step behaviour.
        label andersonStickySteps_;

        //- Last time index at which Anderson was in use, for the sticky window.
        label lastEscalatedIndex_;

        //- Time index this step was escalated to Anderson at, or -1.
        //
        //  Keyed on the INDEX rather than a bool that has to be reset:
        //  prepareRetry() re-runs the same step with runTime_.setTime(..., ti),
        //  which KEEPS the index, so the escalation survives exactly the retry
        //  it was made for and lapses by itself on the next step. A bool would
        //  need clearing in newStepCheck() -- which discardStep() forces to
        //  re-fire by setting timeIndex_ = -1, so it would clear the escalation
        //  immediately and the feature would never engage.
        label escalatedTimeIndex_;




    // Private Member Functions

        //- Reset the history when a new time step starts
        void newStepCheck();

        //- Compute the joint factor and apply it to every enrolled field
        void applyJoint();

        //- Index of a field by name, or -1
        label find(const word& name) const;

        //- No copying
        plasmaOuterRelaxation(const plasmaOuterRelaxation&) = delete;
        void operator=(const plasmaOuterRelaxation&) = delete;


public:

    TypeName("plasmaOuterRelaxation");

    // Constructors

        //- `defaultActive` is the default for `adaptiveRelaxation` when the
        //  case does not set it. It is TRUE for LMEA (an electron-energy
        //  equation is being solved, so the n_e / nEps_e Picard pair exists
        //  and can enter the period-2 cycle this damps) and FALSE otherwise
        //  -- with no energy equation there is no coupled pair, so there is
        //  nothing for a joint factor to do.
        plasmaOuterRelaxation
        (
            const fvMesh& mesh,
            const dictionary& oc,
            const bool defaultActive = false
        );


    //- Destructor
    virtual ~plasmaOuterRelaxation() = default;


    // Selectors

        //- The instance for this mesh, creating it on first use. The
        //  dictionary is read only on creation.
        static plasmaOuterRelaxation& New
        (
            const fvMesh& mesh,
            const dictionary& oc,
            const bool defaultActive = false
        );

        //- The existing instance, or nullptr if none has been created
        static plasmaOuterRelaxation* lookup(const fvMesh& mesh);


    // Member Functions

        bool active() const { return active_; }

        //- Take part. Call once, before the time loop.
        void enrol(volScalarField& f);

        //- Offer this corrector's solved value. Call once per corrector,
        //  immediately after the field's solve.
        void contribute(volScalarField& f);

        //- Abandon the current attempt's history (step rejected / retried)
        void discardStep();

        //- Smallest joint factor this step. 1 means nothing needed damping,
        //  and is what is returned when the class is inactive, so a case that
        //  does not use it imposes no time-step limit.
        //- The time-step controller's coupling stability margin.
        //
        //  Under ANDERSON this deliberately returns 1, i.e. no limit. omega
        //  answers "how hard did the relaxation have to damp", which is
        //  meaningless for a scheme that does not damp -- Anderson cancels the
        //  mode instead. MEASURED, testAitken CASE 17: Anderson m = 2
        //  converged at EVERY stiffness swept, on both the linear and the
        //  clamped bed, so there is no boundary for a pre-emptive governor to
        //  pre-empt. It fails as a cliff, and `retryStep` is what handles a
        //  cliff. Inventing a substitute governor would throttle deltaT for no
        //  reason; see the start-up note.
        scalar margin() const
        {
            if (!active_ || scheme_ == "anderson") return 1.0;
            return omegaMinStep_;
        }

        //- Which scheme the case asked for
        const word& scheme() const { return scheme_; }

        //- The scheme actually driving THIS corrector. Under `adaptive` that
        //  is Aitken until the step has been escalated, Anderson after.
        word effectiveScheme() const;

        //- Escalate this step to Anderson for its retry. Returns true if the
        //  escalation happened, which is the caller's signal to re-run at the
        //  SAME deltaT instead of halving. False when the scheme is not
        //  adaptive, or when this step has already been escalated -- a second
        //  failure means Anderson could not converge it either, and then the
        //  retry ladder should shorten the step as usual.
        bool escalateForRetry();

        //- Worst joint-residual contraction ratio this step, ||r_k||/||r_k-1||.
        //  Below 1 the outer loop contracted; at or above 1 it did not.
        //  Returns 0 when nothing was measured (fewer than two correctors, or
        //  inactive), which the reporter distinguishes from a real 0.
        scalar contraction() const
        {
            return active_ ? contractionMaxStep_ : 0.0;
        }

        scalar omega() const { return omega_; }

        // The RESOLVED configuration, for reporting. Ask, do not re-read.
        scalar omegaStart() const { return omegaStart_; }
        scalar omegaFloor() const { return omegaFloor_; }
        scalar descentLimit() const { return descentLimit_; }

        //- Largest joint factor this step (see margin() for the smallest)
        scalar omegaMaxStep() const { return active_ ? omegaMaxStep_ : 1.0; }
        label nRelaxed() const { return nRelaxed_; }
        label nFields() const { return fields_.size(); }

        //- regIOobject interface: nothing is written
        virtual bool writeData(Ostream&) const { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
