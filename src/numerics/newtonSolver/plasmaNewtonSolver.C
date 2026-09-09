/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

\*---------------------------------------------------------------------------*/

#include "plasmaNewtonSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(plasmaNewtonSolver, 0);
    defineRunTimeSelectionTable(plasmaNewtonSolver, dictionary);
}

// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::plasmaNewtonSolver> Foam::plasmaNewtonSolver::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word solverType(dict.get<word>("type"));

    auto* ctorPtr = dictionaryConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown/unregistered plasmaNewtonSolver type '" << solverType
            << "'.\n"
            << "    If this names a PETSc-based implementation (e.g. 'SNES'),"
            << " the library that provides it\n"
            << "    was most likely not loaded -- add it to this case's"
            << " controlDict, e.g.:\n"
            << "        libs (\"libplasmaNewtonSolverPETSc.so\");\n"
            << "    Registered types: "
            << (dictionaryConstructorTablePtr_
                  ? dictionaryConstructorTablePtr_->sortedToc()
                  : wordList())
            << nl << exit(FatalIOError);
    }

    return autoPtr<plasmaNewtonSolver>(ctorPtr(mesh, dict));
}

// ************************************************************************* //
