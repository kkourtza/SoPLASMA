/*---------------------------------------------------------------------------*\
Class
    Foam::aitkenRelaxation

Description
    Aitken \f$\Delta^2\f$ dynamic under-relaxation for a Picard (fixed-point)
    outer loop, after Irons & Tuck (1969) and as used for partitioned coupling
    by Kuettler & Wall, Comput. Mech. 43 (2008) 61-72.

    WHY THIS EXISTS. The species and electron-energy equations of this solver
    couple through a Picard iteration inside the PIMPLE outer loop. Above a
    critical time step that iteration has a gain > 1 and enters a period-2
    limit cycle: measured on the LMEA streamer benchmark at deltaT = 6.19 ps,
    n_e's residual alternates ~0.03 / ~1.0 in antiphase with nEps_e and grows
    until the Poisson solve is handed a diverged charge density.

    A FIXED under-relaxation factor damps it, but the usable value is both
    case-specific and TIME-VARYING: measured on that same benchmark, 0.9 was
    stable at 6.19 ps and unstable at 8.03 ps IN THE SAME RUN. Worse, a fixed
    factor only moves the stability boundary -- the adaptive time-step
    controller then grows deltaT until it finds the new one. Hence a factor
    that is derived from the iteration itself, per corrector.

    THE METHOD. With residual \f$r_k = \psi_k - \psi_{k-1}\f$ between successive
    outer iterates,

    \f[
        \omega_k = -\omega_{k-1}
        \frac{r_{k-1} \cdot (r_k - r_{k-1})}{|r_k - r_{k-1}|^2}
    \f]

    and the accepted iterate is \f$\psi_{k-1} + \omega_k r_k\f$. For a
    period-2 oscillation successive residuals nearly reverse, so the
    denominator is large and \f$\omega\f$ falls; for a cleanly converging
    iteration successive residuals are nearly parallel and \f$\omega\f$ returns
    towards 1, so no accuracy is paid for damping that is not needed.

    OMEGA IS ALSO THE MARGIN. omega ~ 1 means the coupling is benign at this
    deltaT; a small omega means the iteration is only stable because it is
    being damped hard, i.e. deltaT is near the coupling's stability limit. The
    minimum over a step is therefore a usable signal for the time-step
    controller, which otherwise has no knowledge of this mode at all.

Note
    All reductions are collective and take VALUES, never cell indices, so the
    factor is identical on every rank by construction.

SourceFiles
    aitkenRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef aitkenRelaxation_H
#define aitkenRelaxation_H

#include "scalarField.H"
#include "word.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class aitkenRelaxation Declaration
\*---------------------------------------------------------------------------*/

class aitkenRelaxation
{
    // Private data

        //- Field this instance relaxes, for diagnostics
        word name_;

        //- Residual of the previous outer iteration
        scalarField rPrev_;

        //- Factor used by the previous outer iteration
        scalar omegaPrev_;

        //- Outer iteration counter within the current time step
        label corr_;

        //- Bounds. omega is never allowed to reach 0 (the iterate would
        //  freeze) nor to exceed omegaMax (over-relaxation, which is what we
        //  are trying to prevent).
        scalar omegaMin_;
        scalar omegaMax_;

        //- Factor for the first corrector, where there is no history
        scalar omega0_;

        //- Most the factor may FALL in one corrector, as a fraction of the
        //  previous value. A single bad estimate can otherwise slam omega onto
        //  its floor, and MEASURED on the LMEA case it did: the trace read
        //  ... 0.561 0.184 0.05 0.05 0.539 ... where omega = 0.05 implies a
        //  loop gain of ~19, while a CONSTANT factor of 0.8 runs the same case
        //  stably (gain < 0.25). Those cannot both be true, so the floor hits
        //  were noise -- correctors where the residual is not a clean linear
        //  mode (after the stiff chemistry integration, or across a clamp).
        //
        //  Limiting the DESCENT only: a sustained stiffening still walks omega
        //  down over several correctors, but one outlier cannot. Recovery
        //  upwards is deliberately unlimited -- there is no risk in relaxing
        //  less, and prompt recovery is what keeps the accuracy cost low.
        scalar descentLimit_;

        //- Smallest factor seen since the last reset -- the stability margin
        scalar omegaMinSeen_;

        //- Largest factor seen since the last reset. Reported alongside the
        //  minimum because a single number cannot tell "omega sat at 0.5 all
        //  step" apart from "omega swung between 1.0 and 0.18", and those say
        //  different things about the coupling.
        scalar omegaMaxSeen_;


public:

    // Constructors

        aitkenRelaxation
        (
            const word& name = "undefined",
            const scalar omega0 = 1.0,
            const scalar omegaMin = 0.05,
            const scalar omegaMax = 1.0,
            const scalar descentLimit = 0.0
        );


    // Member Functions

        //- Start a new time step: forget the previous step's history.
        //  MUST be called once per time step. Carrying a residual across a
        //  step would relate iterates of two different problems.
        void reset();

        //- Factor for this corrector, from the residual r = psi - psiPrevIter.
        //  Collective. Call exactly once per corrector per field.
        scalar omega(const scalarField& r);

        //- Smallest factor since reset(); 1 means nothing needed damping.
        scalar margin() const { return omegaMinSeen_; }

        scalar omegaMax() const { return omegaMaxSeen_; }

        //- Correctors taken since reset()
        label nCorrectors() const { return corr_; }

        const word& name() const { return name_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
