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
    reported_(false),
    Icollected_(0),
    IperSpecies_(),
    speciesNames_(),
    chargeFed_(false),
    Qbase_(0),
    QbaseIndex_(-1),
    closure_(0),
    csv_(nullptr),
    headerWritten_(false),
    pendingIndex_(-1),
    pendingTime_(0),
    pendingVf_(0),
    pendingQ_(0),
    pendingI_(0),
    pendingClosure_(0)
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


Foam::floatingElectrode::~floatingElectrode()
{
    if (Pstream::master())
    {
        flushPending();
    }
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

    // THE LEDGER MUST ACTUALLY BE FED.
    //
    // Q(t) = Q0 + INT I_plasma dt' is supplied from OUTSIDE this class, by
    // plasmaTransport::updateSurfaceCharge. If space charge exists but nothing
    // has ever called addCollectedCurrent(), the electrode is silently holding Q0
    // while a discharge deposits charge on it -- which produces an entirely
    // plausible and entirely wrong floating potential rather than an error.
    //
    // Checked at RUNTIME rather than trusted, because "the call site exists"
    // and "the call site is reached" are different claims, and the second one
    // is the one that matters.
    if (gMax(mag(em.chargeDensity().primitiveField())) > 0 && !chargeFed_)
    {
        if (mesh_.time().timeIndex() > mesh_.time().startTimeIndex() + 1)
        {
            FatalErrorInFunction
                << "The floating electrode's charge ledger is NOT BEING FED."
                << nl << nl
                << "    There is space charge in the domain, but nothing has"
                   " called addCollectedCurrent()," << nl
                << "    so Q is still Q0 = " << Q0_ << " C after "
                << mesh_.time().timeIndex() << " steps." << nl << nl
                << "    Q(t) = Q0 + INT I_plasma dt' is supplied by"
                   " plasmaTransport::updateSurfaceCharge." << nl
                << "    A solver that does not call it cannot represent a"
                   " floating electrode in a plasma:" << nl
                << "    the conductor would keep its initial charge while the"
                   " discharge charged it up." << nl
                << exit(FatalError);
        }
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

    // CLOSURE OF GAUSS'S LAW, verified rather than assumed.
    //
    // After the correction the conductor must hold EXACTLY the charge the
    // ledger says it holds. Re-measuring is cheap and it is the one check that
    // catches the subtle failure mode of this whole design: if psi were built
    // from a DIFFERENT operator than the solve used -- the true eps under a
    // semi-implicit scheme, say -- then V + dV_f*psi is NOT a solution, the
    // charge does not land on target, and everything else still looks
    // plausible.
    //
    // Normalised by C_self*1V, i.e. expressed as the potential error it
    // corresponds to, so the number is readable against the voltage scale
    // rather than being an absolute charge nobody can judge.
    {
        const scalar Qafter = measureCharge(em);

        closure_ = mag(Qafter - Q_)/max(Cself_, VSMALL);

        // Loose enough not to trip on the linear solver's own tolerance,
        // tight enough that a wrong operator (which gives an O(1) relative
        // error) cannot pass.
        const scalar Vscale = max(mag(Vf_), scalar(1));

        if (closure_ > 1e-3*Vscale)
        {
            FatalErrorInFunction
                << "Gauss's law does not close on floating conductor `"
                << patchName_ << "`." << nl << nl
                << "    charge target   " << Q_ << " C" << nl
                << "    charge measured " << Qafter << " C" << nl
                << "    discrepancy     " << closure_
                << " V equivalent (Q error / C_self)" << nl << nl
                << "    The superposition V += dV_f*psi is exact ONLY if psi"
                   " solves the SAME operator as" << nl
                << "    the field it corrects. Under `scheme semiImplicit`"
                   " that operator is" << nl
                << "    eps + dt*sigma, NOT eps. A discrepancy of order the"
                   " potential itself means psi" << nl
                << "    was built from the wrong one." << nl
                << exit(FatalError);
        }
    }

    write();

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




void Foam::floatingElectrode::addCollectedCurrent
(
    const wordList& names,
    const scalarList& I,
    const scalar dt
)
{
    speciesNames_ = names;
    IperSpecies_  = I;

    scalar Itot = 0;
    forAll(I, i) { Itot += I[i]; }

    const label ti = mesh_.time().timeIndex();

    // A NEW step: the charge standing now is the accepted baseline.
    if (ti != QbaseIndex_)
    {
        QbaseIndex_ = ti;
        Qbase_      = Q_;
    }

    // BASELINED, not accumulated -- see Qbase_. A retried step recomputes.
    Icollected_ = Itot;
    Q_          = Qbase_ + Itot*dt;
    chargeFed_  = true;

    // KEEP THE CSV ROW SELF-CONSISTENT.
    //
    // correct() runs at the START of a step, using the charge standing at the
    // END of the previous one, and buffers the row. This -- the ledger update
    // -- runs at the END of the same step. Without this the row for step k
    // carried step k's TIME beside step k-1's CHARGE: a one-step lead of the
    // label over the data, which is exactly the trap that made an earlier
    // reconciliation of Q against INT I dt miss by 11.5% and look like a
    // quadrature error.
    if (pendingIndex_ == ti)
    {
        pendingQ_    = Q_;
        pendingI_    = Icollected_;
        pendingPerSpecies_ = IperSpecies_;
        pendingTime_ = mesh_.time().timeOutputValue();
    }
}


void Foam::floatingElectrode::discardStep()
{
    // Undo this step's charge, and forget the row buffered for it: both belong
    // to an attempt that is being thrown away.
    Q_          = Qbase_;
    Icollected_ = 0;
    pendingIndex_ = -1;
}


void Foam::floatingElectrode::flushPending()
{
    if (pendingIndex_ < 0 || !csv_) return;

    if (!headerWritten_)
    {
        headerWritten_ = true;

            csv_() << "# floating conductor `" << patchName_ << "`" << nl
                   << "# C_self = " << Cself_ << " F, Q0 = " << Q0_ << " C" << nl
                   << "#" << nl
                   << "# I_collected IS NOT THE DISCHARGE CURRENT." << nl
                   << "#   I_collected  net charge flux onto THIS ONE conductor,"
                      " summed over species [A]." << nl
                   << "#                A CONDUCTION current only. It is what"
                      " charges this electrode," << nl
                   << "#                and it is what Q(t) integrates." << nl
                   << "#   Sato's I_total, in"
                      " postProcessing/dischargeCurrent/current.csv," << nl
                   << "#                is the current in the EXTERNAL CIRCUIT:"
                      " weighted over the whole" << nl
                   << "#                domain and INCLUDING the displacement"
                      " term. In a barrier discharge" << nl
                   << "#                it is dominated by displacement while"
                      " I_collected can be zero." << nl
                   << "#                The two are different quantities and"
                      " neither checks the other." << nl
                   << "#" << nl
                   << "# one row per TIMESTEP, converged values" << nl
                   << "# Q is integrated as Q += I*dt, a FIRST-ORDER rectangle"
                      " rule, so reconciling this" << nl
                   << "# file against a trapezoid of I will differ at O(dt)."
                   << nl
                   << "# closure is |Q_measured - Q_target|/C_self, in volts"
                   << nl
                   << "time,V_f,Q,I_collected,closure";

            // PER-SPECIES BREAKDOWN, so "do the ions actually contribute?" is a
            // question the output answers rather than one the reader has to infer
            // from a total that moved.
            forAll(speciesNames_, i)
            {
                csv_() << ",I_" << speciesNames_[i];
            }

            csv_() << endl;
    }

    csv_() << pendingTime_ << ','
           << pendingVf_ << ','
           << pendingQ_ << ','
           << pendingI_ << ','
           << pendingClosure_;

    forAll(pendingPerSpecies_, i)
    {
        csv_() << ',' << pendingPerSpecies_[i];
    }

    csv_() << endl;
}


void Foam::floatingElectrode::write()
{
    if (!Pstream::master()) return;

    if (!csv_)
    {
        const fileName dir
        (
            mesh_.time().globalPath()/"postProcessing"/"floatingElectrode"
        );

        mkDir(dir);

        csv_.reset(new OFstream(dir/"floating.csv"));

        // The floating potential is a RESULT, and it is the single most
        // informative number about such an electrode -- written every step so
        // the settling behaviour can be READ rather than assumed. A floating
        // electrode in a plasma should charge NEGATIVE and settle near
        //     V_f - V_plasma ~ -(kTe/2e) ln(2 pi me/mi)
        // which is a few times -kTe/e. One that charges POSITIVE, or that
        // keeps drifting without settling, means the sign of I_collected is wrong
        // or the electron flux to the wall is unresolved.
    }

    const label ti = mesh_.time().timeIndex();

    if (ti != pendingIndex_)
    {
        // The step has advanced, so the buffered row is now final.
        flushPending();
        pendingIndex_ = ti;
    }

    pendingTime_    = mesh_.time().timeOutputValue();
    pendingVf_      = Vf_;
    pendingQ_       = Q_;
    pendingI_       = Icollected_;
    pendingClosure_ = closure_;
    pendingPerSpecies_ = IperSpecies_;
}


// ************************************************************************* //
