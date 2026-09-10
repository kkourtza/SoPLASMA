/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::plasmaNewtonSolver

Description
    Runtime-selectable interface for a genuine Newton-type outer solver,
    replacing the segregated Picard sweep (Poisson -> species -> electron
    energy) for cases where that fixed-point iteration cannot converge --
    see docs/design/newton-outer-solver-design.md for the full diagnosis (a
    negative-differential-resistance region of the discharge's own
    current-voltage characteristic) and the JFNK/PETSc-SNES proofs of
    concept this implementation is built on.

    DELIBERATELY PETSc-FREE. This header (and the base library it belongs
    to) never sees a PETSc header, so `soPlasmaFoam` links against it
    UNCONDITIONALLY with no new dependency -- exactly like every other
    runtime-selected SoPLASMA model. The actual PETSc-based implementation
    lives in a SEPARATE, optionally-built library
    (`libplasmaNewtonSolverPETSc.so`, only built when this project's own
    PETSc environment is available), following the SAME pattern OpenFOAM's
    own `petsc4Foam` uses for its (optional, dynamically-loaded) PETSc
    linear-solver plugin: a case that wants `outerSolver newton` adds that
    library to its `controlDict`'s `libs` entry; one that does not is
    completely unaffected, and selecting `newton` without the library
    loaded fails with a clear, actionable error rather than a link error
    or a crash.

    Exposed as a `system/plasmaSimulationControls` numerics switch
    (`outerCoupling.outerSolver picard | newton`), not a generator/
    semantic-layer concept -- see the design doc's exposure decision.

SourceFiles
    plasmaNewtonSolver.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaNewtonSolver_H
#define plasmaNewtonSolver_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

// Forward declarations -- by reference/pointer only in the interface below,
// so the real headers (and everything THEY pull in) are never needed here.
class fvMesh;
class electromagneticsModel;
class plasmaSpecies;
class plasmaTransport;
class plasmaEnergy;

/*---------------------------------------------------------------------------*\
                     Class plasmaNewtonSolver Declaration
\*---------------------------------------------------------------------------*/

class plasmaNewtonSolver
{
public:

    //- Runtime type information
    TypeName("plasmaNewtonSolver");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,
        plasmaNewtonSolver,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (mesh, dict)
    );

    //- Selector. FatalErrors, naming the missing library, if `dict`'s
    //  `type` entry names an implementation that never registered itself
    //  (i.e. the library that provides it was not loaded via `libs`).
    static autoPtr<plasmaNewtonSolver> New
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~plasmaNewtonSolver() = default;

    // Member Functions

        //- Replace ONE Picard sweep (the existing, unchanged
        //  Poisson -> species -> electron energy sequence
        //  soPlasmaFoam.C's `picard` branch runs every `pimple.loop()`
        //  iteration) with one Newton (outer) iteration. Implementations
        //  reach into these SAME live, already-constructed objects --
        //  nothing is reconstructed or duplicated.
        //
        //  `energy` is nullptr under LFA (no transported electron energy
        //  field to include in the coupled system).
        virtual void solveOuterStep
        (
            electromagneticsModel& em,
            plasmaSpecies& species,
            plasmaTransport& transport,
            plasmaEnergy* energy
        ) = 0;

        //- Did the LAST solveOuterStep() converge?
        //
        //  REPORTED RATHER THAN FATAL, so a non-converged Newton step can be
        //  DISCARDED AND RETRIED AT A SHORTER dt through exactly the machinery
        //  Picard already uses (`outerCoupling/onNonConvergence retryStep`),
        //  instead of killing the run.
        //
        //  The asymmetry this removes was real and it was the binding one:
        //  Picard could back off from a step it could not solve and Newton
        //  could not, so on any handover Newton inherited whatever dt Picard's
        //  easy phase had wound up to and died on contact with it. Measured
        //  2026-09-10 on grubert2009_ballast400: Newton is handed dt = 1.29e-9
        //  and fails immediately, while the SAME state solves cleanly at
        //  dt <= 1e-11 (154 and 256 steps, zero failures). The step was simply
        //  too large, and a solver that cannot say so has no way to find that
        //  out.
        //
        //  Default true, so an implementation that never fails need not
        //  override it.
        virtual bool lastSolveConverged() const
        {
            return true;
        }
};

} // End namespace Foam

#endif

// ************************************************************************* //
