/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::snesNewtonSolver

Description
    Concrete `plasmaNewtonSolver` implementation using PETSc's matrix-free
    SNES -- see doc/newton-outer-solver-design.md for the full design and
    the proofs of concept this is built on.

    NEVER includes any PETSc header directly -- petscmath.h and OpenFOAM's
    own math operator overloads are ambiguous when both are visible in one
    translation unit (a real build obstacle found and fixed in the proof of
    concept). All PETSc calls live in `snesBridge.C`, reached through the
    plain-types interface in `snesBridge.H`, exactly as the proof of concept
    established.

    Physics-based (Knoll & Keyes) preconditioned matrix-free JFNK: one
    block Gauss-Seidel sweep (Poisson, then species, then energy -- the
    SAME order and operators the real Picard sequence uses) as the
    preconditioner, each equation solving a linearised CORRECTION against
    a synthetic GMRES right-hand side. Bare/unpreconditioned JFNK was
    tried first (matching how the single-equation proof of concept was
    validated) and found, 2026-09-09, to reliably hit
    SNES_DIVERGED_LINEAR_SOLVE on the real system's first genuine Newton
    iteration -- expected, matching the design doc's diagnosis that this
    coupled system is inherently stiff (fast electron transport against
    Poisson) independent of the later negative-differential-resistance
    region.

    Field block, one instance per rank, flat-concatenated in this fixed
    order: [ePotential][n_0]...[n_{nSpecies-1}][nEps_e if LMEA]. The outer
    residual for EVERY block is computed via EXPLICIT fvc:: evaluation,
    never fvm::+.residual() -- confirmed broken in parallel for that
    specific usage, independent of PETSc (see
    [[fvmatrix-residual-broken-in-parallel-use-fvc-instead]]).

SourceFiles
    snesNewtonSolver.C

\*---------------------------------------------------------------------------*/

#ifndef snesNewtonSolver_H
#define snesNewtonSolver_H

#include "plasmaNewtonSolver.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class snesNewtonSolver Declaration
\*---------------------------------------------------------------------------*/

class snesNewtonSolver
:
    public plasmaNewtonSolver
{
    const fvMesh& mesh_;

    //- `-snes_rtol`/`-snes_max_it`-equivalent, read from the dict rather
    //  than expecting PETSc command-line flags (which argList would reject
    //  as unrecognised -- a real obstacle found in the proof of concept).
    scalar rtol_;
    label maxIt_;

public:

    TypeName("SNES");

    snesNewtonSolver(const fvMesh& mesh, const dictionary& dict);

    virtual ~snesNewtonSolver() = default;

    virtual void solveOuterStep
    (
        electromagneticsModel& em,
        plasmaSpecies& species,
        plasmaTransport& transport,
        plasmaEnergy* energy
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
