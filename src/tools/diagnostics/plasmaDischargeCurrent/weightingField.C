/*---------------------------------------------------------------------------*\
  File: weightingField.C
  Part of: SoPLASMA -- plasmaDischargeCurrent

  Sato's geometric weighting field, solved MONOLITHICALLY across the gas and
  every dielectric region.

  Copyright (C) 2026
  License: GNU General Public License v3 or later
\*---------------------------------------------------------------------------*/

#include "plasmaDischargeCurrent.H"
#include "electromagneticsModel.H"
#include "multiRegionPoisson.H"
#include "fixedValueFvPatchFields.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * Static Helpers  * * * * * * * * * * * * * * * //

namespace Foam
{

//- A psi_hat field on one region, CLONED from that region's ePotential.
//
//  Cloning rather than building a patch-type list by hand is the whole trick,
//  and getting it wrong is silent.
//
//  The monolithic assembly couples regions through `mappedPatchBase` patches
//  carrying `useImplicit true` (multiRegionPoisson asserts this and refuses
//  mixed coupling). Those boundary conditions cannot be reconstructed from a
//  type NAME alone -- they need their dictionary entries (sample mode, the
//  neighbour region and patch). Constructing psi_hat from a hand-built
//  wordList of types would therefore produce default-constructed interface
//  BCs, i.e. a DECOUPLED solve that still runs and still returns a
//  plausible-looking field.
//
//  So: copy the real field, then overwrite only what must differ -- the
//  electrode patches become fixedValue 1 or 0, and the interior starts at
//  zero. Everything else, interfaces and constraint patches included, is the
//  same object the Poisson solve uses.
static autoPtr<volScalarField> clonePsiHat
(
    const volScalarField& ePot,
    const word& name,
    const word& drivenPatch,
    const wordList& groundedPatches
)
{
    auto psi = autoPtr<volScalarField>::New
    (
        IOobject
        (
            name,
            ePot.mesh().time().constant(),
            ePot.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ePot                            // clone: BCs come across intact
    );

    psi->primitiveFieldRef() = 0.0;
    psi->dimensions().reset(dimless);

    volScalarField::Boundary& bf = psi->boundaryFieldRef();

    forAll(bf, patchi)
    {
        const word& pname = bf[patchi].patch().name();

        const bool driven   = (pname == drivenPatch);
        const bool grounded = groundedPatches.found(pname);

        if (!driven && !grounded)
        {
            // Not an electrode: leave the cloned BC exactly as it is. On a
            // dielectric region this is how the interface coupling survives.
            bf[patchi] == 0.0;
            continue;
        }

        // An electrode. Replace the BC outright rather than assigning into
        // it: the real ePotential electrode may be a time-varying condition
        // (a ramp, a table, a coded BC), which would re-evaluate on
        // updateCoeffs and silently overwrite the unit value this field
        // depends on.
        bf.set
        (
            patchi,
            new fixedValueFvPatchScalarField
            (
                bf[patchi].patch(),
                psi->internalField()
            )
        );

        bf[patchi] == (driven ? 1.0 : 0.0);
    }

    return psi;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::plasmaDischargeCurrent::wedgeRevolutionFactor
(
    const fvMesh& m
)
{
    // On a wedge mesh, m.V() holds only the theta-slice volumes, so every
    // extensive quantity -- charge, capacitance, current -- must be scaled to
    // the full 2*pi body of revolution. Taken from the geometry rather than a
    // dictionary entry, which is exactly the kind of thing that silently
    // disagrees with the mesh.
    const polyBoundaryMesh& pbm = m.boundaryMesh();

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

    return 1.0;      // 3-D or planar 2-D: volumes are already what they claim
}


void Foam::plasmaDischargeCurrent::computeWeightingField
(
    const electromagneticsModel& em
)
{
    // UNIT applied voltage, deliberately: psi_hat = 1 on the driven electrode
    // makes the field time-independent even under a ramped or AC drive, which
    // is what makes computing it once legitimate. V_a re-enters only through
    // dV_a/dt in the displacement term.
    psiHatGas_ = clonePsiHat
    (
        em.ePotential(), "psiHat", drivenPatch_, groundedPatches_
    );

    const bool multi = isA<multiRegionPoisson>(em);

    const multiRegionPoisson* mrp =
        multi ? &refCast<const multiRegionPoisson>(em) : nullptr;

    const label nDiel = mrp ? mrp->nDielectrics() : 0;

    psiHatDiel_.setSize(nDiel);

    for (label i = 0; i < nDiel; ++i)
    {
        psiHatDiel_.set
        (
            i,
            clonePsiHat
            (
                mrp->dielectric(i).ePotential(),
                "psiHat",
                drivenPatch_,
                groundedPatches_
            )
        );
    }

    const dimensionedScalar& epsGas = em.epsilon();

    // MONOLITHIC across regions when the solver is coupled -- which is the
    // default here, and the solver paper's headline contribution (Pasolari &
    // Kourtzanidis, arXiv:2607.05137): a monolithic multi-region Poisson
    // assembly over arbitrarily many curved conforming interfaces.
    //
    // Assembling psi_hat the same way the real Poisson solve is assembled
    // gives continuous psi with continuous eps dpsi/dn across every interface
    // BY CONSTRUCTION, rather than by any interface handling written here.
    const bool coupled = (nDiel > 0) && mrp->coupled();

    for (label nonOrth = 0; nonOrth <= em.nNonOrthCorr(); ++nonOrth)
    {
        if (coupled)
        {
            fvMatrix<scalar> assembly
            (
                psiHatGas_(),
                dimensionSet(0, 0, 1, 0, 0, 1, 0)
            );

            fvScalarMatrix gasEqn(fvm::laplacian(epsGas, psiHatGas_()));
            assembly.addFvMatrix(gasEqn);

            for (label i = 0; i < nDiel; ++i)
            {
                fvScalarMatrix dielEqn
                (
                    fvm::laplacian(mrp->epsilon(i), psiHatDiel_[i])
                );
                assembly.addFvMatrix(dielEqn);
            }

            assembly.solve();

            psiHatGas_().correctBoundaryConditions();
            for (label i = 0; i < nDiel; ++i)
            {
                psiHatDiel_[i].correctBoundaryConditions();
            }
        }
        else
        {
            // Single region, or a multi-region case the user configured as
            // segregated. Segregated psi_hat would need its own outer loop to
            // converge the interfaces; rather than half-implement that, it is
            // refused below when dielectrics are present.
            fvScalarMatrix psiEqn(fvm::laplacian(epsGas, psiHatGas_()));
            psiEqn.solve();
            psiHatGas_().correctBoundaryConditions();
        }
    }

    if (nDiel > 0 && !coupled)
    {
        FatalErrorInFunction
            << "dischargeCurrent requires the MONOLITHIC (coupled) Poisson"
            << " path when dielectric regions are present." << nl
            << "    This case has " << nDiel << " dielectric region(s) with"
            << " region coupling disabled, so Sato's weighting field cannot"
            << " be made continuous across the interfaces here." << nl
            << "    Enable the coupled solve, or disable dischargeCurrent."
            << nl
            << exit(FatalError);
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
            -fvc::grad(psiHatGas_())
        )
    );

    revolutionFactor_ = wedgeRevolutionFactor(mesh_);

    // C_g = INT eps |e_hat|^2 dV over the WHOLE domain, dielectrics included.
    //
    // Restricting the integral to the gas would put a bounding surface on the
    // dielectric face, where psi' is not zero -- and that surface term is
    // precisely what Morrow & Sato's derivation drops. Keeping the dielectric
    // interior keeps the bounding surface at electrodes and infinity.
    Cg_ = revolutionFactor_*epsGas.value()
        * gSum(mesh_.V()*magSqr(eHat_().primitiveField()));

    for (label i = 0; i < nDiel; ++i)
    {
        const fvMesh& rm = mrp->dielectric(i).mesh();
        const volVectorField eR(-fvc::grad(psiHatDiel_[i]));

        // Each region carries its own revolution factor: nothing guarantees a
        // dielectric mesh was built with the same wedge angle as the gas.
        Cg_ += wedgeRevolutionFactor(rm)*mrp->epsilon(i).value()
             * gSum(rm.V()*magSqr(eR.primitiveField()));
    }

    Info<< "plasmaDischargeCurrent: Sato weighting field solved" << nl
        << "    driven patch      " << drivenPatch_ << nl
        << "    grounded patches  " << groundedPatches_ << nl
        << "    regions           " << (nDiel + 1)
        << (coupled ? "  (monolithic)" : "  (single region)") << nl
        << "    revolution factor " << revolutionFactor_
        << (revolutionFactor_ > 1.0 ? "  (wedge -> full 2pi)" : "") << nl
        << "    gap capacitance   " << Cg_ << " F" << nl
        << "    max |e_hat|       " << gMax(mag(eHat_().primitiveField()))
        << " 1/m" << endl;
}


// ************************************************************************* //
