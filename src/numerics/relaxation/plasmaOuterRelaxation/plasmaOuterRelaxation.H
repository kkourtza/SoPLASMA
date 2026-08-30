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
        scalar margin() const { return active_ ? omegaMinStep_ : 1.0; }

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
