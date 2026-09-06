/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.
    Copyright (C) 2026
    GNU General Public License v3 or later.

Class
    Foam::plasmaStepAudit

Description
    PER-STEP CALL AUDIT: asserts that every update which MUST run once per
    timestep actually ran.

    WHY THIS EXISTS. On 2026-09-06 `species.updateChargeDensity()` sat behind
    `if (pimple.finalIter())`, which was correct only because
    `pimpleControl::loop()` guarantees a final iteration when residualControl
    converges. A `relativeChange` convergence measure added the same day exits
    the loop with a raw `break` at the TOP of an iteration, bypassing that
    guarantee -- so the update ran on 151 of 364670 timesteps (0.04%).

    Poisson's source was therefore frozen. The field stayed at the VACUUM
    value (E/N measured 1104 Td against V/L = 1106 Td, 0.2% apart) while the
    ion density grew three decades past the point where space charge must
    dominate, so nothing screened and nothing arrested ionisation. A day went
    into wall-flux boundary forms, Scharfetter-Gummel, ion mobility tables,
    the mean-energy clamp, ballast RC timescales and load-line algebra --
    every hypothesis plausible, every one wrong.

    THE GENERAL DEFECT, which is why this is a registry and not one check: a
    per-step update reached through a GUARD stops being reached when any exit
    path changes, and nothing in the solution announces it. The symptom is a
    stale derived field, which looks exactly like physics.

    This turns that class of failure into one line, at the step it happens:

        *** PER-STEP UPDATE DID NOT RUN: updateChargeDensity ***

    See CLAUDE.md rule 27 and doc/case-monitoring-plan.md.

    Usage: declare what must run once at set-up, record each call, and close
    the step.

        plasmaStepAudit::expect("updateChargeDensity");
        ...
        plasmaStepAudit::record("updateChargeDensity");   // in the callee
        ...
        plasmaStepAudit::endStep(runTime.timeIndex());    // once per step

    Cost is a hash lookup per call, so it is left ON always. A diagnostic that
    has to be enabled is a diagnostic nobody has enabled.

    HEADER-ONLY, deliberately. `plasmaTools` already depends on
    `plasmaModels` (plasmaTimeControl takes a plasmaTransport), so a linked
    class here could not be called from plasmaSpecies without a circular
    dependency. C++17 inline statics remove the question entirely: any module
    that can see this header can use it, and nothing links.

\*---------------------------------------------------------------------------*/

#ifndef plasmaStepAudit_H
#define plasmaStepAudit_H

#include "HashTable.H"
#include "word.H"
#include "label.H"
#include "error.H"
#include "IOstreams.H"

namespace Foam
{

class plasmaStepAudit
{
    //- Calls recorded since the last endStep().
    inline static HashTable<label> thisStep_;

    //- Total calls over the run, for the closing summary.
    inline static HashTable<label> total_;

    //- Names that MUST be recorded at least once per step.
    inline static HashTable<label> required_;

    //- Steps closed so far (audited steps, not solver timeIndex).
    inline static label nSteps_ = 0;

    //- Off until the first expect(); a solver that registers nothing pays
    //  nothing and is not nagged.
    inline static bool active_ = false;

public:

    //- Declare that `name` must run at least once per timestep.
    inline static void expect(const word& name)
    {
        required_.set(name, 1);
        thisStep_.set(name, 0);
        total_.set(name, 0);
        active_ = true;
    }

    //- Record one call. Cheap enough to sit in the callee unconditionally.
    inline static void record(const word& name)
    {
        if (!active_) return;

        auto it = thisStep_.find(name);
        if (it.good()) { ++(*it); } else { thisStep_.set(name, 1); }

        auto tt = total_.find(name);
        if (tt.good()) { ++(*tt); } else { total_.set(name, 1); }
    }

    //- Close the step: every required name must have been recorded. A name
    //  with zero calls is FATAL -- it means a derived quantity is stale, and
    //  a stale derived quantity is indistinguishable from a physics result.
    inline static void endStep(const label timeIndex)
    {
        if (!active_) return;

        forAllConstIters(required_, iter)
        {
            const word& name = iter.key();

            if (thisStep_.lookup(name, 0) == 0)
            {
                // FATAL, not a warning. A warning is worse than nothing here:
                // the run continues, the derived field goes stale, and the
                // result looks like physics. Measured 2026-09-06 -- 364670
                // steps with a frozen Poisson source.
                FatalErrorInFunction
                    << "PER-STEP UPDATE DID NOT RUN: `" << name << "`"
                    << nl << nl
                    << "    It is declared as required once per timestep and"
                    << " was recorded ZERO times" << nl
                    << "    in step " << timeIndex << ". Whatever it derives"
                    << " is now STALE, and a stale" << nl
                    << "    derived field is indistinguishable from a physics"
                    << " result -- which is" << nl
                    << "    exactly how this class of bug survives." << nl << nl
                    << "    LOOK FOR AN EXIT PATH THAT SKIPS IT. The known"
                    << " cause is a call guarded by" << nl
                    << "    `pimple.finalIter()` in a corrector loop left with"
                    << " a raw `break`:" << nl
                    << "    pimpleControl::loop() promises a final iteration"
                    << " when residualControl" << nl
                    << "    converges, and a `break` voids that promise."
                    << nl << nl
                    << "    See CLAUDE.md rule 27 and"
                    << " doc/case-monitoring-plan.md." << nl
                    << exit(FatalError);
            }
        }

        forAllIters(thisStep_, iter) { *iter = 0; }
        ++nSteps_;
    }

    //- End-of-run summary: calls per step for each audited update.
    inline static void report()
    {
        if (!active_ || nSteps_ == 0) return;

        Info<< "plasmaStepAudit: " << nSteps_ << " steps audited" << nl;

        forAllConstIters(total_, iter)
        {
            const scalar perStep = scalar(iter())/scalar(nSteps_);

            Info<< "    " << iter.key() << ": " << iter() << " calls, "
                << perStep << " per step";

            // Print it rather than only asserting on it: 0.00041 was the
            // measured value when this was broken.
            if (perStep < 1.0) Info<< "  <-- BELOW ONE PER STEP";
            Info<< nl;
        }
        Info<< endl;
    }

    //- Calls recorded for `name` in the current (open) step.
    inline static label countThisStep(const word& name)
    {
        return thisStep_.lookup(name, 0);
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
