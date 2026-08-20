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


void Foam::plasmaDischargeCurrent::update
(
    const plasmaTransport& transport,
    const plasmaSpecies& species,
    const electromagneticsModel& em
)
{
    if (!enabled_) return;

    const scalar t  = mesh_.time().value();
    const scalar Va = appliedVoltage(em);

    // CONDUCTION. Morrow & Sato Eq. (24)'s bracket is exactly SUM_s q_s
    // Gamma_s, with q_s the charge NUMBER: the paper writes it out as
    // (N_p W_p - N_e W_e - N_n W_n - D_p grad N_p + D_e grad N_e
    //  + D_n grad N_n), which is that sum once the signs are collected.
    //
    // particleFlux_ is already Gamma_s . Sf: convectiveFlux (fvm::div(phi,n)
    // .flux(), the drift) plus diffusiveFlux (-fvm::laplacian(D,n).flux(),
    // i.e. -D grad n . Sf). So the diffusive sign the paper asks for is
    // carried by construction and must NOT be re-applied here.
    const scalar qe = constant::plasma::eCharge.value();

    const List<scalar>& qs = species.speciesChargeNumbers();

    scalar Icond = 0;

    perSpeciesCurrent_.setSize(qs.size(), 0.0);

    forAll(transport.particleFlux(), s)
    {
        perSpeciesCurrent_[s] = 0.0;

        if (s >= qs.size() || mag(qs[s]) < SMALL) continue;   // neutrals

        // IMMOBILE species carry no flux field at all -- the entry in the
        // PtrList is never set, and dereferencing it aborts the run (found
        // exactly that way: "Cannot dereference nullptr at index 6 in range
        // [0,13)", the ions on the streamer case being `immobile`).
        //
        // Skipping them is not merely defensive, it is the physics: by
        // Shockley-Ramo only charge IN MOTION induces current in the external
        // circuit, so a stationary population contributes exactly zero.
        if (!transport.particleFlux().set(s)) continue;

        // Face flux -> cell-centred Gamma. Established idiom in
        // plasmaTransport.C.
        const volVectorField G(fvc::reconstruct(transport.particleFlux(s)));

        const scalar Is =
            qs[s]*gSum(mesh_.V()*(G.primitiveField() & eHat_().primitiveField()));

        perSpeciesCurrent_[s] = qe*revolutionFactor_*Is;
        Icond += Is;
    }

    Icond *= qe*revolutionFactor_;

    // DISPLACEMENT. In Sato's form this does NOT need dE/dt anywhere: the
    // space-charge part cancels out of the derivation and what survives is
    // the gap capacitance times the rate of change of the APPLIED voltage.
    // That is the whole numerical advantage over the electrode surface
    // integral, which must difference conduction against displacement where
    // the two nearly cancel.
    scalar Idisp = 0;

    if (seeded_ && t > prevTime_)
    {
        Idisp = Cg_*(Va - prevVoltage_)/(t - prevTime_);
    }

    prevVoltage_ = Va;
    prevTime_ = t;
    seeded_ = true;

    const scalar Itot = Icond + Idisp;

    // OUTPUT. The cadence is independent of the field write interval and may
    // be far shorter: this is one scalar per step, while a nanosecond current
    // pulse is destroyed by sampling at the field cadence.
    const label ti = mesh_.time().timeIndex();

    if (writeInterval_ > 0 && (ti % writeInterval_) == 0)
    {
        if (!os_ && Pstream::master())
        {
            // globalPath, NOT path: in a parallel run time().path() resolves
            // to the RANK's processorN directory, so the file landed in
            // processor0/postProcessing where nobody would look for it.
            const fileName dir
            (
                mesh_.time().globalPath()/"postProcessing"/"dischargeCurrent"
            );
            mkDir(dir);

            os_.reset(new OFstream(dir/"current.csv"));

            os_() << "# Sato's equation (Morrow & Sato, J Phys D 32 (1999) L20)"
                  << nl
                  << "# C_g = " << Cg_ << " F, revolution factor "
                  << revolutionFactor_ << nl
                  << "time,V_applied,I_total,I_cond,I_disp";

            if (perSpecies_)
            {
                forAll(qs, s)
                {
                    if (mag(qs[s]) < SMALL) continue;
                    os_() << "," << "I_" << species.speciesNames()[s];
                }
            }
            os_() << endl;
        }

        if (Pstream::master())
        {
            os_() << t << ',' << Va << ',' << Itot << ',' << Icond << ','
                  << Idisp;

            if (perSpecies_)
            {
                forAll(qs, s)
                {
                    if (mag(qs[s]) < SMALL) continue;
                    os_() << ',' << perSpeciesCurrent_[s];
                }
            }
            os_() << endl;
        }
    }

    if (printInterval_ > 0 && (ti % printInterval_) == 0)
    {
        Info<< "  discharge current:        " << Itot << " A"
            << "   (conduction " << Icond << ", displacement " << Idisp
            << ", V = " << Va << ")" << endl;
    }
}


// ************************************************************************* //
