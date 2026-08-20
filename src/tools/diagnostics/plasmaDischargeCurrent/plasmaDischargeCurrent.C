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
#include "multiRegionPoisson.H"

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

    // Patch validation spans ALL regions, not just the gas.
    //
    // In a real DBD the electrodes sit on different meshes -- the driven one
    // against the gas, the grounded one behind the dielectric slab -- so
    // checking either name against the gas mesh alone would reject a
    // perfectly valid setup. Found on the plate2D two-region case, where
    // `left` is a gas patch and `right` belongs to the dielectric.
    wordList allPatches;
    {
        DynamicList<word> names;

        for (const word& w : mesh_.boundaryMesh().names())
        {
            names.append(w);
        }

        if (isA<multiRegionPoisson>(em))
        {
            const multiRegionPoisson& mrp =
                refCast<const multiRegionPoisson>(em);

            for (label i = 0; i < mrp.nDielectrics(); ++i)
            {
                for (const word& w : mrp.dielectric(i).mesh().boundaryMesh().names())
                {
                    names.append(w);
                }
            }
        }

        allPatches.transfer(names);
    }

    if (!allPatches.found(drivenPatch_))
    {
        FatalErrorInFunction
            << "dischargeCurrent/drivenPatch `" << drivenPatch_
            << "` is not a patch of any region." << nl
            << "    Available: " << allPatches << nl
            << exit(FatalError);
    }

    // A misspelt grounded patch would otherwise silently keep its cloned
    // boundary condition, changing C_g with no error at all -- the failure
    // would surface only as a wrong current, long after the fact.
    for (const word& p : groundedPatches_)
    {
        if (!allPatches.found(p))
        {
            FatalErrorInFunction
                << "dischargeCurrent/groundedPatches names `" << p
                << "`, which is not a patch of any region." << nl
                << "    Available: " << allPatches << nl
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
