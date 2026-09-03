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
#include "mappedPatchBase.H"

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
// clonePsiHat MOVED to Foam::unitPotentialField (2026-09-03), where the
// floating electrode shares it. See unitPotentialField.C for the interface
// retargeting and why it is load-bearing.

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// wedgeRevolutionFactor MOVED to unitPotentialField as a static member.

void Foam::plasmaDischargeCurrent::computeWeightingField
(
    const electromagneticsModel& em
)
{
    // UNIT applied voltage, deliberately: psi_hat = 1 on the driven electrode
    // makes the field time-independent even under a ramped or AC drive, which
    // is what makes computing it once legitimate. V_a re-enters only through
    // dV_a/dt in the displacement term.
    //
    // The solve itself lives in unitPotentialField, shared with the floating
    // electrode. Everything Sato-specific stays here: which patch is the
    // measuring electrode, e_hat, and the ENERGY form of the capacitance.
    psiHat_.reset
    (
        new unitPotentialField
        (
            mesh_,
            "psiHat",
            drivenPatch_,
            groundedPatches_,
            string("Or disable dischargeCurrent.")
        )
    );

    psiHat_->solve(em);

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
            -fvc::grad(psiHat_->gas())
        )
    );

    revolutionFactor_ = unitPotentialField::wedgeRevolutionFactor(mesh_);

    // C_g = INT eps |e_hat|^2 dV over the WHOLE domain, dielectrics included.
    //
    // Restricting the integral to the gas would put a bounding surface on the
    // dielectric face, where psi' is not zero -- and that surface term is
    // precisely what Morrow & Sato's derivation drops. Keeping the dielectric
    // interior keeps the bounding surface at electrodes and infinity.
    Cg_ = psiHat_->energyCapacitance(em);

    Info<< "plasmaDischargeCurrent: Sato weighting field solved" << nl
        << "    driven patch      " << drivenPatch_ << nl
        << "    grounded patches  " << groundedPatches_ << nl
        << "    regions           " << (psiHat_->nDielectrics() + 1)
        << (psiHat_->coupled() ? "  (monolithic)" : "  (single region)") << nl
        << "    revolution factor " << revolutionFactor_
        << (revolutionFactor_ > 1.0 ? "  (wedge -> full 2pi)" : "") << nl
        << "    gap capacitance   " << Cg_ << " F" << nl
        << "    max |e_hat|       " << gMax(mag(eHat_().primitiveField()))
        << " 1/m" << endl;
}


// ************************************************************************* //


void Foam::plasmaDischargeCurrent::writeFields() const
{
    if (psiHat_ && psiHat_->solved())
    {
        const volScalarField& pg = psiHat_->gas();
        pg.write();

        Info<< "  psiHat (gas):        min " << gMin(pg.primitiveField())
            << "  max " << gMax(pg.primitiveField())
            << "  volume " << gSum(pg.mesh().V()) << " m^3" << endl;
    }

    if (eHat_)
    {
        eHat_->write();
    }

    if (psiHat_ && psiHat_->solved())
    {
        const PtrList<volScalarField>& pd = psiHat_->dielectrics();

        forAll(pd, i)
        {
            pd[i].write();

            Info<< "  psiHat (dielectric " << i << "): min "
                << gMin(pd[i].primitiveField())
                << "  max " << gMax(pd[i].primitiveField())
                << "  volume " << gSum(pd[i].mesh().V()) << " m^3" << endl;
        }
    }
}


// ************************************************************************* //
