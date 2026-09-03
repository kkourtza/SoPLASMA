/*---------------------------------------------------------------------------*\
  File: floatingElectrode.C
  Part of: SoPLASMA
\*---------------------------------------------------------------------------*/

#include "floatingElectrode.H"
#include "electromagneticsModel.H"
#include "multiRegionPoisson.H"
#include "floatingPotentialInterface.H"
#include "plasmaConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floatingElectrode::floatingElectrode
(
    const fvMesh& gasMesh,
    const electromagneticsModel& em
)
:
    mesh_(gasMesh),
    patchName_(word::null),
    patchi_(-1),
    Q0_(0),
    Q_(0),
    Vf_(0),
    Cself_(0),
    psi_(nullptr),
    cacheValid_(false),
    builtAtMeshEvent_(-1),
    reported_(false)
{
    wordList driven, grounded, floating;
    em.classifyElectrodePatches(driven, grounded, floating);

    if (floating.empty())
    {
        FatalErrorInFunction
            << "floatingElectrode constructed with no floating conductor."
            << nl << exit(FatalError);
    }

    if (floating.size() > 1)
    {
        FatalErrorInFunction
            << "More than one floating conductor: " << floating << nl << nl
            << "    TWO OR MORE FLOATING CONDUCTORS COUPLE TO EACH OTHER, so a"
               " single self-capacitance" << nl
            << "    is not enough: the constraint becomes an N x N MUTUAL"
               " CAPACITANCE MATRIX, requiring" << nl
            << "    N unit-potential solves and a linear solve of that matrix"
               " each step." << nl << nl
            << "    That is not implemented. One floating conductor is"
               " supported." << nl
            << exit(FatalError);
    }

    patchName_ = floating[0];

    patchi_ = mesh_.boundaryMesh().findPatchID(patchName_);

    if (patchi_ < 0)
    {
        // The classification scans every region, so a floating patch on a
        // DIELECTRIC mesh is found and named -- but the correction below adds
        // dV_f*psi using the gas mesh's patch, and unitPotentialField's
        // surface integral is taken on the gas mesh. Rather than return a
        // wrong self-capacitance, say so.
        FatalErrorInFunction
            << "The floating conductor `" << patchName_ << "` is not a patch"
            << " of the gas region `" << mesh_.name() << "`." << nl << nl
            << "    A floating electrode buried in a DIELECTRIC region is not"
               " supported: its charge" << nl
            << "    would have to be integrated on that region's mesh, and"
               " the self-capacitance with it." << nl
            << "    Only a conductor in contact with the gas is handled." << nl
            << exit(FatalError);
    }

    // Q0 comes from the boundary condition, which is where the user states it.
    const fvPatchScalarField& pf = em.ePotential().boundaryField()[patchi_];

    const floatingPotentialInterface* fep =
        dynamic_cast<const floatingPotentialInterface*>(&pf);

    if (!fep)
    {
        FatalErrorInFunction
            << "Patch `" << patchName_ << "` was classified as a floating"
            << " conductor but its condition does not implement"
            << " floatingPotentialInterface." << nl
            << "    Its type is `" << pf.type() << "`." << nl
            << exit(FatalError);
    }

    Q0_ = fep->initialCharge();
    Q_  = Q0_;

    // A restart resumes at the potential it had reached.
    Vf_ = fep->floatingPotential();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::floatingElectrode::build
(
    const electromagneticsModel& em,
    const volScalarField* effEpsGas
)
{
    wordList driven, grounded, floating;
    em.classifyElectrodePatches(driven, grounded, floating);

    // psi = 1 on the floating conductor, 0 on every electrode whose potential
    // IS known. The other electrodes must be at zero, not at their actual
    // potentials: psi is the homogeneous solution, and their contribution is
    // already carried by the solved field this corrects.
    wordList zeros(driven.size() + grounded.size());
    label k = 0;
    for (const word& w : driven)   zeros[k++] = w;
    for (const word& w : grounded) zeros[k++] = w;

    if (zeros.empty())
    {
        FatalErrorInFunction
            << "A floating conductor needs at least one electrode with a KNOWN"
            << " potential." << nl
            << "    Without one, the problem is singular: the potential is"
               " determined only up to a" << nl
            << "    constant, and so is the floating potential." << nl
            << exit(FatalError);
    }

    psi_.reset
    (
        new unitPotentialField
        (
            mesh_,
            "psiFloat",
            patchName_,
            zeros,
            string("A floating electrode requires it.")
        )
    );

    if (effEpsGas)
    {
        psi_->solve(em, *effEpsGas);
    }
    else
    {
        psi_->solve(em);
    }

    Cself_ = psi_->surfaceCapacitance(em, patchName_);

    if (Cself_ <= VSMALL)
    {
        FatalErrorInFunction
            << "The self-capacitance of floating conductor `" << patchName_
            << "` came out as " << Cself_ << " F, which cannot be." << nl
            << "    A conductor raised to 1 V with every other electrode"
               " grounded must carry POSITIVE" << nl
            << "    charge. A non-positive value means the unit-potential"
               " solve or the surface" << nl
            << "    integral is wrong, not that the geometry is unusual." << nl
            << exit(FatalError);
    }

    cacheValid_ = true;
    builtAtMeshEvent_ = mesh_.topoChanging() ? -1 : mesh_.time().timeIndex();
}


Foam::scalar Foam::floatingElectrode::measureCharge
(
    const electromagneticsModel& em
) const
{
    // Q = + INT eps snGrad(V) dA, the SAME sign convention as
    // unitPotentialField::surfaceCapacitance -- see the derivation there. The
    // two must agree, because dV_f = (Q_target - Q)/C_self divides one by the
    // other: a sign slip in either would produce a correction that drives the
    // charge AWAY from its target, and the run would diverge rather than
    // give a slightly wrong answer.
    const fvPatchScalarField& pf = em.ePotential().boundaryField()[patchi_];

    const scalarField& magSf = pf.patch().magSf();

    return unitPotentialField::wedgeRevolutionFactor(mesh_)
         * em.epsilon().value()
         * gSum(pf.snGrad()*magSf);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::floatingElectrode::present(const electromagneticsModel& em)
{
    wordList driven, grounded, floating;
    em.classifyElectrodePatches(driven, grounded, floating);
    return !floating.empty();
}


void Foam::floatingElectrode::correct
(
    electromagneticsModel& em,
    const volScalarField* effEpsGas
)
{
    // CACHE VALIDITY, decided rather than assumed.
    //
    // psi solves div(op grad psi) = 0, so it is a constant of the OPERATOR and
    // the MESH. It may be reused only when neither has moved:
    //
    //   * semiImplicit -- the operator is eps + dt*sigma. Both dt (adaptive)
    //     and sigma (the plasma) change every step, so REBUILD EVERY STEP.
    //     This is the default scheme, so this is the normal path.
    //   * explicit -- the operator is eps. Reusable while the mesh is static.
    //   * a changing mesh invalidates it either way.
    const bool operatorMoves = (effEpsGas != nullptr);

    const bool meshMoved =
        mesh_.changing()
     || mesh_.topoChanging()
     || builtAtMeshEvent_ < 0;

    if (!cacheValid_ || operatorMoves || meshMoved)
    {
        build(em, effEpsGas);
    }

    // THE PLASMA CHARGE LEDGER IS NOT IMPLEMENTED YET (2026-09-03).
    //
    // Q(t) = Q0 + INT I_plasma dt' requires the net charge flux onto the
    // conductor from the species wall fluxes. addCharge() is the entry point
    // and NOTHING CALLS IT, so Q is held at Q0 for the whole run.
    //
    // With no plasma that is exactly right and the constraint is exact --
    // validated against the analytic closed form to 4.5e-09 of the voltage
    // scale. WITH a plasma it is WRONG, and wrong quietly: the electrode would
    // hold its initial charge while the discharge deposited charge on it, and
    // the resulting floating potential would look entirely plausible.
    //
    // So the plasma case is REFUSED rather than approximated. A floating
    // electrode charges negative in a plasma -- that is the whole physical
    // point of one -- and a run that cannot represent it must not pretend to.
    const bool haveSpaceCharge =
        gMax(mag(em.chargeDensity().primitiveField())) > 0;

    if (effEpsGas || haveSpaceCharge)
    {
        FatalErrorInFunction
            << "A floating electrode in a PLASMA run is not supported yet."
            << nl << nl
            << "    Detected: "
            << (effEpsGas ? "a semi-implicit Poisson operator (so a"
                            " conductivity was supplied)" : "")
            << (effEpsGas && haveSpaceCharge ? " and " : "")
            << (haveSpaceCharge ? "non-zero space charge" : "")
            << "." << nl << nl
            << "    The constraint itself is exact and validated, but the"
               " CHARGE LEDGER is not wired:" << nl
            << "        Q(t) = Q0 + INT I_plasma dt'" << nl
            << "    needs the net charge flux onto the conductor from the"
               " species wall fluxes, and" << nl
            << "    nothing supplies it. Q would stay at Q0 for the whole run"
               " while the discharge" << nl
            << "    deposited charge on the electrode -- giving a plausible"
               " but wrong floating" << nl
            << "    potential rather than an error." << nl << nl
            << "    A floating electrode in a plasma charges NEGATIVE, to"
               " roughly -(kTe/2e) ln(2 pi me/mi)," << nl
            << "    and that behaviour is exactly what the missing term"
               " produces." << nl << nl
            << "    Electrostatics-only runs are fully supported. See"
               " docs/models/poisson_equation/floating-electrode.md" << nl
            << exit(FatalError);
    }

    // The field has just been solved with the patch at Vf_. Whatever charge
    // that implies is generally NOT the charge the electrode actually holds.
    const scalar Qmeasured = measureCharge(em);

    const scalar dVf = (Q_ - Qmeasured)/Cself_;

    Vf_ += dVf;

    // EXACT, not iterative. L(psi) = 0, so V + dVf*psi is still a solution of
    // the same equation, and its value on the floating patch is exactly
    // Vf_old + dVf. Applying it to EVERY region matters: psi is continuous
    // across the interfaces by construction of the monolithic solve, and
    // correcting only the gas would leave a jump at the dielectric face.
    volScalarField& ePot = em.ePotentialRef();

    ePot.primitiveFieldRef() += dVf*psi_->gas().primitiveField();

    if (isA<multiRegionPoisson>(em))
    {
        multiRegionPoisson& mrp = refCast<multiRegionPoisson>(em);

        const PtrList<volScalarField>& psiD = psi_->dielectrics();

        for (label i = 0; i < mrp.nDielectrics(); ++i)
        {
            mrp.dielectric(i).ePotential().primitiveFieldRef()
                += dVf*psiD[i].primitiveField();
        }
    }

    // Publish the new equipotential to the boundary condition, then let the
    // boundary conditions re-evaluate so the corrected field is consistent.
    {
        fvPatchScalarField& pf = ePot.boundaryFieldRef()[patchi_];

        dynamic_cast<floatingPotentialInterface&>(pf)
            .setFloatingPotential(Vf_);
    }

    ePot.correctBoundaryConditions();

    if (isA<multiRegionPoisson>(em))
    {
        multiRegionPoisson& mrp = refCast<multiRegionPoisson>(em);

        for (label i = 0; i < mrp.nDielectrics(); ++i)
        {
            mrp.dielectric(i).ePotential().correctBoundaryConditions();
        }
    }

    if (!reported_)
    {
        reported_ = true;

        Info<< "floatingElectrode: `" << patchName_ << "`" << nl
            << "    self-capacitance  " << Cself_ << " F" << nl
            << "    initial charge    " << Q0_ << " C" << nl
            << "    floating potential " << Vf_ << " V" << nl
            << "    psi rebuilt       "
            << (operatorMoves ? "every step (semiImplicit operator)"
                              : "once (explicit operator, static mesh)")
            << endl;
    }
}


// ************************************************************************* //
