/*---------------------------------------------------------------------------*\
  File: plasmaDischargeCurrent.C
  Part of: SoPLASMA
  Copyright (C) 2026
  License: GNU General Public License v3 or later
\*---------------------------------------------------------------------------*/

#include "plasmaDischargeCurrent.H"
#include "electromagneticsModel.H"
#include "plasmaConstants.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasmaDischargeCurrent::plasmaDischargeCurrent
(
    const fvMesh& mesh,
    const dictionary& dict,
    const electromagneticsModel& em
)
:
    mesh_(mesh)
{
    if (!dict.found("dischargeCurrent")) return;

    const dictionary& cd = dict.subDict("dischargeCurrent");

    enabled_ = cd.getOrDefault<Switch>("enabled", false);
    if (!enabled_) return;

    drivenPatch_     = cd.get<word>("drivenPatch");
    groundedPatches_ = cd.getOrDefault<wordList>("groundedPatches", wordList());
    perSpecies_      = cd.getOrDefault<Switch>("perSpecies", false);
    writeInterval_   = cd.getOrDefault<label>("writeInterval", 1);
    printInterval_   = cd.getOrDefault<label>("printInterval", 0);

    if (mesh_.boundaryMesh().findPatchID(drivenPatch_) < 0)
    {
        FatalErrorInFunction
            << "dischargeCurrent/drivenPatch `" << drivenPatch_
            << "` is not a patch of this mesh." << nl
            << "    Available: " << mesh_.boundaryMesh().names() << nl
            << exit(FatalError);
    }

    // Every grounded patch must exist too. A misspelt name would otherwise
    // silently become a zeroGradient boundary on the weighting field, which
    // changes C_g without any error -- the failure would surface only as a
    // wrong current, long after the fact.
    for (const word& p : groundedPatches_)
    {
        if (mesh_.boundaryMesh().findPatchID(p) < 0)
        {
            FatalErrorInFunction
                << "dischargeCurrent/groundedPatches names `" << p
                << "`, which is not a patch of this mesh." << nl
                << "    Available: " << mesh_.boundaryMesh().names() << nl
                << exit(FatalError);
        }
    }

    if (groundedPatches_.empty())
    {
        FatalErrorInFunction
            << "dischargeCurrent needs at least one grounded patch." << nl
            << "    Without a reference the weighting-field problem is"
            << " singular: psi is determined only up to a constant, so"
            << " e_hat and C_g are meaningless." << nl
            << exit(FatalError);
    }

    computeWeightingField(em);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::plasmaDischargeCurrent::appliedVoltage
(
    const electromagneticsModel& em
) const
{
    const label patchi = mesh_.boundaryMesh().findPatchID(drivenPatch_);
    const fvPatchScalarField& pf = em.ePotential().boundaryField()[patchi];

    // Area-weighted, so a non-uniform electrode potential still gives the
    // single number the circuit sees.
    const scalar a = gSum(mesh_.boundary()[patchi].magSf());

    if (a <= VSMALL) return 0;

    return gSum(mesh_.boundary()[patchi].magSf()*pf)/a;
}


// ************************************************************************* //
