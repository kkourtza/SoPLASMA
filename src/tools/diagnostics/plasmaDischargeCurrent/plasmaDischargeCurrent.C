/*---------------------------------------------------------------------------*\
  File: plasmaDischargeCurrent.C
  Part of: SoPLASMA
  Copyright (C) 2026
  License: GNU General Public License v3 or later
\*---------------------------------------------------------------------------*/

#include "mappedPatchBase.H"
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
    mesh_(mesh),
    Itot_(0.0),
    Icond_(0.0),
    measured_(false)
{
    // ON BY DEFAULT, and the sub-dictionary is OPTIONAL.
    //
    // Sato's discharge current is the primary measurable of almost every
    // discharge simulation -- the one number an experiment can be compared
    // against -- so a case should not have to ask for it. It costs ONE extra
    // Poisson solve at start-up and two surface integrals per write.
    //
    // `dischargeCurrent/enabled false` remains the explicit opt-out.
    dictionary cd;
    if (dict.found("dischargeCurrent"))
    {
        cd = dict.subDict("dischargeCurrent");
    }

    enabled_ = cd.getOrDefault<Switch>("enabled", true);
    if (!enabled_) return;

    // Both patch entries are now OPTIONAL and DERIVED when absent: the
    // potential's own boundary conditions already say which patch is driven
    // and which are grounded, and restating it is a second source of truth
    // that can disagree. An explicit entry still wins.
    drivenPatch_     = cd.getOrDefault<word>("drivenPatch", word::null);
    groundedPatches_ = cd.getOrDefault<wordList>("groundedPatches", wordList());
    perSpecies_      = cd.getOrDefault<Switch>("perSpecies", false);
    writeInterval_   = cd.getOrDefault<label>("writeInterval", 1);
    crossCheck_      = cd.getOrDefault<Switch>("crossCheck", false);
    printInterval_   = cd.getOrDefault<label>("printInterval", 0);

    if (drivenPatch_.empty() || groundedPatches_.empty())
    {
        deriveElectrodePatches(em);
    }

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

    // A REGION INTERFACE CANNOT BE A GROUNDED PATCH.
    //
    // computeWeightingField() solves a Laplace problem with the grounded
    // patches forced to fixedValue, REPLACING whatever condition is there. On a
    // coupled interface that destroys the coupling, and the neighbour side's
    // coupledElectricPotential then aborts in refCast with
    //   "Attempt to cast type fixedValue to type coupledElectricPotential"
    // -- from inside THIS constructor, with an error that reads like a Poisson
    // problem and says nothing about groundedPatches. Measured 2026-08-30 on a
    // needle-DBD case; only the stack trace identified it.
    //
    // Checked on the PATCH TYPE rather than the field's boundary type, because
    // the interface is a mappedWall whatever condition sits on it, and the
    // check must work before any field has been constructed.
    for (const word& pName : groundedPatches_)
    {
        const label pi = mesh_.boundaryMesh().findPatchID(pName);

        if (pi >= 0 && isA<mappedPatchBase>(mesh_.boundaryMesh()[pi]))
        {
            FatalErrorInFunction
                << "dischargeCurrent/groundedPatches names `" << pName
                << "`, which is a REGION INTERFACE (mappedWall)." << nl << nl
                << "    The weighting field is solved with the grounded patches"
                   " forced to fixedValue," << nl
                << "    which replaces the interface coupling and makes the"
                   " neighbour side abort." << nl << nl
                << "    Name a real electrode. If the true ground sits behind"
                   " a dielectric -- a DBD -- then" << nl
                << "    it is in ANOTHER REGION, and that is fine: the"
                   " weighting field is solved" << nl
                << "    MONOLITHICALLY across every region, so a ground on a"
                   " dielectric mesh is" << nl
                << "    reached correctly. Name that patch."
                << nl << nl
                << "    (SUPERSEDED 2026-09-03: this message used to say the"
                   " single-region weighting" << nl
                << "    solve could not reach another region and advised"
                   " disabling the diagnostic." << nl
                << "    That has been false since the monolithic multi-region"
                   " assembly landed, which" << nl
                << "    is validated to 1.1e-15 against analytic on a"
                   " two-region series stack.)"
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

const Foam::volScalarField*
Foam::plasmaDischargeCurrent::findPotentialForPatch
(
    const electromagneticsModel& em,
    const word& patchName,
    label& patchi
) const
{
    patchi = em.ePotential().mesh().boundaryMesh().findPatchID(patchName);

    if (patchi >= 0)
    {
        return &em.ePotential();
    }

    if (isA<multiRegionPoisson>(em))
    {
        const multiRegionPoisson& mrp = refCast<const multiRegionPoisson>(em);

        for (label i = 0; i < mrp.nDielectrics(); ++i)
        {
            const volScalarField& ePot = mrp.dielectric(i).ePotential();

            patchi = ePot.mesh().boundaryMesh().findPatchID(patchName);

            if (patchi >= 0)
            {
                return &ePot;
            }
        }
    }

    patchi = -1;
    return nullptr;
}


void Foam::plasmaDischargeCurrent::deriveElectrodePatches
(
    const electromagneticsModel& em
)
{
    // THE CLASSIFICATION RULE HAS ONE OWNER:
    // electromagneticsModel::classifyElectrodePatches(). It is shared with the
    // floating electrode, which needs the floating list and needs these two as
    // the patches its unit-potential problem holds at zero.
    //
    // Moved there 2026-09-03 when floatingElectrodePotential arrived: that
    // condition IS a fixedValue by inheritance, so a copy of the rule living
    // here would have classified a floating conductor as a DRIVEN electrode
    // and measured the discharge current at it.
    wordList drivenCandidates;
    wordList grounded;
    wordList floating;

    em.classifyElectrodePatches(drivenCandidates, grounded, floating);

    // --- the driven electrode -----------------------------------------------
    if (drivenPatch_.empty())
    {
        if (drivenCandidates.size() == 1)
        {
            drivenPatch_ = drivenCandidates[0];
        }
        else
        {
            FatalErrorInFunction
                << "dischargeCurrent is ON BY DEFAULT but the DRIVEN electrode"
                << " could not be derived." << nl << nl
                << "    Driven candidates found: " << drivenCandidates << nl
                << "    Grounded patches found:  " << grounded << nl
                << "    Floating conductors:     " << floating
                << "   (never a current-measuring electrode)" << nl << nl
                << (
                       drivenCandidates.empty()
                     ? "    NONE were found. A driven electrode is a Dirichlet"
                       " `ePotential` condition that is\n"
                       "    either time-varying (uniformFixedValue with a"
                       " table/sine/ramp) or a non-zero\n"
                       "    constant. If this case genuinely has no driven"
                       " electrode, there is no\n"
                       "    discharge current to measure.\n"
                     : "    MORE THAN ONE was found, so which one the current"
                       " is measured at is a\n"
                       "    physical choice, not something to guess.\n"
                   )
                << nl
                << "    Set it explicitly:" << nl
                << "        dischargeCurrent { drivenPatch <name>; }" << nl
                << "    or turn the diagnostic off:" << nl
                << "        dischargeCurrent { enabled false; }" << nl
                << exit(FatalError);
        }
    }

    // --- the grounded electrodes --------------------------------------------
    if (groundedPatches_.empty())
    {
        // The driven patch is never also a ground, even if a derivation quirk
        // put it in both lists.
        DynamicList<word> g;
        for (const word& w : grounded)
        {
            if (w != drivenPatch_) g.append(w);
        }


        if (g.empty())
        {
            FatalErrorInFunction
                << "dischargeCurrent is ON BY DEFAULT but no GROUNDED"
                << " electrode could be derived." << nl << nl
                << "    Driven patch: " << drivenPatch_ << nl
                << "    A ground is a `fixedValue` `ePotential` condition equal"
                << " to zero." << nl << nl
                << "    Without a reference the weighting-field problem is"
                   " SINGULAR: psi is determined" << nl
                << "    only up to a constant, so e_hat and C_g are"
                   " meaningless." << nl << nl
                << "    Set it explicitly:" << nl
                << "        dischargeCurrent { groundedPatches ( <name> ); }"
                << nl
                << "    or turn the diagnostic off:" << nl
                << "        dischargeCurrent { enabled false; }" << nl
                << exit(FatalError);
        }

        groundedPatches_.transfer(g);
    }

    Info<< "plasmaDischargeCurrent: electrode patches DERIVED from the"
        << " ePotential boundary conditions" << nl
        << "    driven   " << drivenPatch_
        << "   (imposed, time-varying or non-zero)" << nl
        << "    grounded " << groundedPatches_ << "   (imposed, zero)" << nl
        << "    Override either with dischargeCurrent/{drivenPatch,"
        << "groundedPatches}." << endl;
}


Foam::scalar Foam::plasmaDischargeCurrent::appliedVoltage
(
    const electromagneticsModel& em
) const
{
    // Searched across ALL regions. This used to be a gas-mesh lookup whose
    // -1 was fed straight into boundaryField()[patchi] -- so a drive behind a
    // barrier, which is an ordinary DBD, was undefined behaviour rather than
    // an error.
    label patchi = -1;
    const volScalarField* ePotPtr =
        findPotentialForPatch(em, drivenPatch_, patchi);

    if (!ePotPtr)
    {
        FatalErrorInFunction
            << "dischargeCurrent/drivenPatch `" << drivenPatch_
            << "` was not found on any region while reading the applied"
            << " voltage." << nl
            << exit(FatalError);
    }

    const fvPatchScalarField& pf = ePotPtr->boundaryField()[patchi];

    // Area-weighted, so a non-uniform electrode potential still gives the
    // single number the circuit sees.
    const scalarField& magSf = pf.patch().magSf();

    const scalar a = gSum(magSf);

    if (a <= VSMALL) return 0;

    return gSum(magSf*pf)/a;
}


// ************************************************************************* //


void Foam::plasmaDischargeCurrent::surfaceIntegralCurrent
(
    const plasmaTransport& transport,
    const plasmaSpecies& species,
    const electromagneticsModel& em,
    scalar& Icond,
    scalar& Idisp
)
{
    Icond = 0;
    Idisp = 0;

    const label patchi = mesh_.boundaryMesh().findPatchID(drivenPatch_);

    if (patchi < 0) return;      // electrode lives on another region's mesh

    const scalar qe = constant::plasma::eCharge.value();
    const List<scalar>& qs = species.speciesChargeNumbers();

    // CONDUCTION at the electrode. particleFlux_ on the patch is already
    // Gamma_s . Sf with the OUTWARD normal, so this is the particle current
    // leaving the gas into the electrode.
    forAll(transport.particleFlux(), s)
    {
        if (s >= qs.size() || mag(qs[s]) < SMALL) continue;
        if (!transport.particleFlux().set(s)) continue;      // immobile

        const scalarField& pf =
            transport.particleFlux()[s].boundaryField()[patchi];

        Icond += qs[s]*gSum(pf);
    }

    Icond *= qe*revolutionFactor_;

    // DISPLACEMENT at the electrode: eps dE/dt . Sf.
    //
    // This is the term Sato's volume form exists to avoid. It differentiates
    // the field in time AT THE BOUNDARY, where the condition is imposed, and
    // near the electrode it very nearly cancels the conduction term -- so the
    // difference of two large numbers is what produces the small answer. Any
    // disagreement between the two routes is expected to live here.
    const scalar eps = em.epsilon().value();

    const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
    const vectorField Ep(em.E().boundaryField()[patchi]);

    scalarField eFlux(eps*(Ep & Sfp));

    const scalar dt = mesh_.time().deltaTValue();

    if (eFluxSeeded_ && prevEFlux_.size() == eFlux.size() && dt > 0)
    {
        Idisp = revolutionFactor_*gSum(eFlux - prevEFlux_)/dt;
    }

    prevEFlux_ = eFlux;
    eFluxSeeded_ = true;
}


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

    // Kept for plasmaExternalCircuit, which drops the ballast voltage across
    // THIS current rather than deriving its own.
    Itot_ = Itot;
    Icond_ = Icond;
    measured_ = true;

    // CROSS-CHECK. The identity is
    //     INT_electrode J_tot . dS = - INT_V e_hat . J_tot dV
    // with the outward normal, so the surface route carries the OPPOSITE sign
    // to Sato's volume route. It is negated here on that basis -- from the
    // derivation, not by fitting it to the answer -- so that agreement means
    // agreement.
    scalar IsurfC = 0, IsurfD = 0;

    if (crossCheck_)
    {
        surfaceIntegralCurrent(transport, species, em, IsurfC, IsurfD);
        IsurfC = -IsurfC;
        IsurfD = -IsurfD;
    }

    const scalar Isurf = IsurfC + IsurfD;

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
                  << "#" << nl
                  << "# I_total IS THE EXTERNAL-CIRCUIT CURRENT: the current an"
                     " ammeter in series with" << nl
                  << "# the electrode would read. It is weighted over the WHOLE"
                     " domain by Sato's" << nl
                  << "# equation and it INCLUDES the displacement term"
                     " (I_total = I_cond + I_disp)." << nl
                  << "# In a barrier discharge it is dominated by I_disp." << nl
                  << "#" << nl
                  << "# DO NOT CONFUSE IT with I_collected in" << nl
                  << "# postProcessing/floatingElectrode/floating.csv, which is"
                     " the net charge flux" << nl
                  << "# onto ONE conductor -- a conduction current only, local"
                     " to that patch, and the" << nl
                  << "# thing that charges a floating electrode. Neither is a"
                     " check on the other." << nl
                  << "#" << nl
                  << "time,V_applied,I_total,I_cond,I_disp";

            if (crossCheck_)
            {
                os_() << ",I_surface,I_surface_cond,I_surface_disp";
            }

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

            if (crossCheck_)
            {
                os_() << ',' << Isurf << ',' << IsurfC << ',' << IsurfD;
            }

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

        if (crossCheck_)
        {
            const scalar rel =
                mag(Itot - Isurf)/max(max(mag(Itot), mag(Isurf)), VSMALL);

            Info<< "    surface-integral check: " << Isurf << " A"
                << "   relative difference " << rel << endl;
        }
    }
}


// ************************************************************************* //
