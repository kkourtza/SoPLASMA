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
#include "wedgePolyPatch.H"

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

Foam::scalar Foam::plasmaDischargeCurrent::wedgeRevolutionFactor() const
{
    // On a wedge mesh, mesh.V() holds only the theta-slice volumes, so every
    // extensive quantity -- charge, capacitance, current -- must be scaled to
    // the full 2*pi body of revolution.
    //
    // The wedge half-angle comes from the patch normal against the wedge
    // axis-plane: OpenFOAM's wedgePolyPatch exposes cosAngle() as the cosine
    // of the angle between the patch normal and the CENTRE plane normal, so
    // the full opening angle is 2*acos(cosAngle).
    //
    // Taken from the geometry rather than from a user entry, because a wedge
    // angle typed into a dictionary is exactly the kind of thing that silently
    // disagrees with the mesh.
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchi)
    {
        if (isA<wedgePolyPatch>(pbm[patchi]))
        {
            const wedgePolyPatch& wp =
                refCast<const wedgePolyPatch>(pbm[patchi]);

            const scalar halfAngle =
                Foam::acos(min(max(wp.cosAngle(), scalar(-1)), scalar(1)));

            const scalar theta = 2.0*halfAngle;

            if (theta > SMALL)
            {
                return constant::mathematical::twoPi/theta;
            }
        }
    }

    // No wedge patch: a 3-D or planar-2-D mesh, where the volumes are already
    // what they claim to be.
    return 1.0;
}


void Foam::plasmaDischargeCurrent::computeWeightingField
(
    const electromagneticsModel& em
)
{
    // UNIT applied voltage, deliberately. Solving for psi = 1 on the driven
    // electrode makes the field time-independent even under a ramped or AC
    // drive, which is what makes computing it once legitimate. V_a re-enters
    // only through dV_a/dt in the displacement term.
    wordList patchTypes(mesh_.boundary().size(), zeroGradientFvPatchScalarField::typeName);

    forAll(mesh_.boundary(), patchi)
    {
        const word& pname = mesh_.boundary()[patchi].name();
        const word& ptype = mesh_.boundary()[patchi].type();

        // Constraint patches (wedge, empty, symmetry, processor) must keep
        // their own type or the solve is ill-posed on an axisymmetric mesh.
        if
        (
            ptype == "wedge" || ptype == "empty" || ptype == "symmetry"
         || ptype == "symmetryPlane" || ptype == "processor"
         || ptype == "cyclic" || ptype == "processorCyclic"
        )
        {
            patchTypes[patchi] = ptype;
        }
        else if (pname == drivenPatch_ || groundedPatches_.found(pname))
        {
            patchTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField psiHat
    (
        IOobject
        (
            "psiHat",
            mesh_.time().constant(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0),
        patchTypes
    );

    forAll(mesh_.boundary(), patchi)
    {
        const word& pname = mesh_.boundary()[patchi].name();

        if (pname == drivenPatch_)
        {
            psiHat.boundaryFieldRef()[patchi] == 1.0;
        }
        else if (groundedPatches_.found(pname))
        {
            psiHat.boundaryFieldRef()[patchi] == 0.0;
        }
    }

    // MULTI-REGION: refuse rather than mislead.
    //
    // Dielectrics here are SEPARATE MESHES (dielectricRegion owns its own
    // regionMesh_, ePotential_ and uniform epsilon_), coupled to the gas at
    // the interfaces. So psi_hat must be solved as a COUPLED multi-region
    // problem -- continuous psi with continuous eps dpsi/dn across every
    // interface -- exactly like the real Poisson solve, and C_g must integrate
    // over the gas mesh AND every dielectric mesh with its own eps.
    //
    // Solving the gas region alone would produce a perfectly plausible field
    // and a quietly wrong C_g, with nothing to indicate it. Not implemented
    // yet, so it is refused.
    if (isA<multiRegionPoisson>(em))
    {
        const multiRegionPoisson& mrp = refCast<const multiRegionPoisson>(em);

        if (mrp.nDielectrics() > 0)
        {
            FatalErrorInFunction
                << "dischargeCurrent is not yet implemented for multi-region"
                << " (dielectric) cases." << nl
                << "    This case carries " << mrp.nDielectrics()
                << " dielectric region(s)." << nl
                << "    Sato's weighting field must be solved as a COUPLED"
                << " multi-region problem, and the gap capacitance integrated"
                << " over the dielectrics too; solving the gas region alone"
                << " gives a plausible but WRONG current." << nl
                << "    Disable dischargeCurrent, or use a single-region case."
                << nl
                << exit(FatalError);
        }
    }

    // eps-weighted, NOT the vacuum Laplace field. multiRegionPoisson exists,
    // so dielectrics are real here, and the vacuum field would misweight every
    // contribution in a dielectric's shadow.
    //
    // NOTE: electromagneticsModel currently exposes permittivity as a single
    // dimensionedScalar for the gas region. That is exact for a bare gap --
    // the streamer case -- and for a uniform dielectric. A multi-region case
    // with contrasting permittivities needs eps as a FIELD here; see
    // docs/discharge-current-sato-plan.md Sect. 6, which flags the dielectric
    // treatment as unverified against a DBD source.
    const dimensionedScalar& eps = em.epsilon();

    for (label i = 0; i <= 2; ++i)     // non-orthogonal correctors
    {
        fvScalarMatrix psiEqn
        (
            fvm::laplacian(eps, psiHat)
        );

        psiEqn.solve();
    }

    eHat_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "eHat",
                mesh_.time().constant(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -fvc::grad(psiHat)
        )
    );

    // AXISYMMETRY. On a wedge mesh, mesh_.V() holds the volumes of the WEDGE
    // only -- a theta/2pi slice of the real body of revolution. Every volume
    // integral therefore has to be scaled to the full revolution or the
    // current comes out low by that same factor (a 5 degree wedge understates
    // it by 72x). Charge, capacitance and current are all extensive, so this
    // is not optional.
    revolutionFactor_ = wedgeRevolutionFactor();

    // C_g = INT eps |e_hat|^2 dV. Falls out of the same field and is a
    // constant of the geometry. On a parallel plate this must equal
    // eps0 A / d, which is the cheapest decisive test of the whole solve --
    // no plasma involved.
    Cg_ = revolutionFactor_*eps.value()*gSum
    (
        mesh_.V()*magSqr(eHat_().primitiveField())
    );

    Info<< "plasmaDischargeCurrent: Sato weighting field solved" << nl
        << "    driven patch      " << drivenPatch_ << nl
        << "    grounded patches  " << groundedPatches_ << nl
        << "    gap capacitance   " << Cg_ << " F" << nl
        << "    max |e_hat|       " << gMax(mag(eHat_().primitiveField()))
        << " 1/m" << nl
        << "    revolution factor " << revolutionFactor_
        << (revolutionFactor_ > 1.0 ? "  (wedge -> full 2pi)" : "")
        << endl;
}


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
