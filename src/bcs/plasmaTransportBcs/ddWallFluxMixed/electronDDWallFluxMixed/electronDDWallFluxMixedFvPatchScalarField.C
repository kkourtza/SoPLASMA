/*---------------------------------------------------------------------------*\
  File: electronDDWallFluxMixedFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::electronDDWallFluxMixedFvPatchScalarField.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "mappedPatchBase.H"
#include "IOdictionary.H"

#include "plasmaTransport.H"
#include "electronDDWallFluxMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    electronDDWallFluxMixedFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

dimensionedScalar
electronDDWallFluxMixedFvPatchScalarField::defaultTValue() const
{
    // 1 eV. Electrons arriving at a wall are not thermalised with the gas --
    // using T_gas (300 K, 0.026 eV) understates the thermal velocity by a
    // factor of ~6 -- and an LFA case has no T_e field to follow because the
    // default `gasTemperature` energy model registers none. 1-2 eV is the
    // usual quoted range; 1 eV is the conservative end.
    //
    // Reported, not silent: a wall temperature nobody chose should not go
    // unannounced.
    const scalar T1eV =
        constant::plasma::eCharge.value()
       /constant::plasma::kappaBoltzmann.value();

    Info<< "    " << patch().name() << "/" << this->internalField().name()
        << ": no `T` or `TeV` given, using the default wall electron"
           " temperature 1 eV (" << T1eV << " K)." << nl
        << "      Set `TeV <eV>` to change it, or `T <fieldName>` to follow a"
           " field." << endl;

    return dimensionedScalar("T", dimTemperature, T1eV);
}


tmp<scalarField>
electronDDWallFluxMixedFvPatchScalarField::normalisingDensity() const
{
    // For n_e the normalising density IS this field.
    return tmp<scalarField>::New(static_cast<const scalarField&>(*this));
}


tmp<scalarField>
electronDDWallFluxMixedFvPatchScalarField::wallLossSpeed
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal,
    scalarField& gRatio
) const
{
    // A -- the thermal base, eq. (6.6) = vT/sqrt(pi) = 2 x eq. (6.3).
    tmp<scalarField> tW = calcThermalVelocity(m, T);
    scalarField& W = tW.ref();

    const scalarField A(W);

    gRatio.setSize(W.size());
    gRatio = Zero;

    // Dd -- the drift addition of eq. (6.8), OUTSIDE the max by construction.
    scalarField Dd(W.size(), Zero);
    if (includeDriftFlux_)
    {
        Dd = max(scalar(0), uDriftNormal);
    }

    // Gc/n -- the creation term of eq. (6.2), normalised by the density that
    // eq. (6.6) actually divides by.
    scalarField GcOverN(W.size(), Zero);
    if (emittedValid_ && emitted_.size() == W.size())
    {
        const tmp<scalarField> tn = normalisingDensity();
        const scalarField& n = tn();

        forAll(GcOverN, faceI)
        {
            // Guarded on the density: with n at or below zero there is no
            // meaningful Gamma_w/n, and the loss speed is left as the thermal
            // one rather than being handed an infinity.
            if (n[faceI] > VSMALL)
            {
                GcOverN[faceI] = emitted_[faceI]/n[faceI];
            }
        }
    }

    const scalar r = electronReflection_;

    // THE CLAMP INSIDE hagelaarClosure() IS PHYSICS, NOT A NUMERICAL GUARD --
    // it is eq. (6.8)'s max, whose floor is Dd because the drift is added
    // outside it. Hagelaar: it fires "for electrons at the cathode due to
    // secondary emission by ion impact, and it is then indeed appropriate to
    // set w_w = 0 because the secondary emission coefficients given in the
    // literature have generally been deduced NEGLECTING THERMAL ELECTRON LOSS
    // TO THE CATHODE."
    forAll(W, faceI)
    {
        hagelaarClosure
        (
            A[faceI], Dd[faceI], GcOverN[faceI], r, W[faceI], gRatio[faceI]
        );
    }

    // n is the lagged wall value, so Gc/n is one evaluation behind -- the
    // relation is implicit in n and there is nothing else to use.

    return tW;
}


tmp<scalarField>
electronDDWallFluxMixedFvPatchScalarField::calcAbsorptionVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal
) const
{
    // The absorption speed IS the total loss speed W.
    scalarField gRatio;
    return wallLossSpeed(m, T, uDriftNormal, gRatio);
}


tmp<scalarField>
electronDDWallFluxMixedFvPatchScalarField::calcEffectiveWallVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal
) const
{
    // The mixed condition imposes n_p*(uDrift_n + uEff), so uEff is defined by
    // the bookkeeping identity uEff = W - uDrift_n. ONE expression feeds both
    // routines: they were provably equal before (on BOTH includeDriftFlux
    // branches) and keeping them separate is how a fix lands on one sibling
    // and not the other.
    scalarField gRatio;
    tmp<scalarField> tVel = wallLossSpeed(m, T, uDriftNormal, gRatio);

    tVel.ref() -= uDriftNormal;

    return tVel;
}


namespace
{

//- Gamma from the MATERIAL of the region behind an interface patch, or -1.
//
//  Secondary emission is a surface property, but for a MESHED barrier the
//  emitting surface IS that region, so gamma belongs to its material
//  declaration in constant/<region>/electricalProperties -- not to a boundary
//  condition on the far side of the interface, where it would be the same
//  number written twice.
//
//  dielectricRegion registers `gammaSEE` on its own mesh registry alongside
//  `epsilonR`. An interface patch is a mappedPatchBase, so the gas side can
//  reach the neighbour's registry through it and read what the material says.
//
//  Returns -1 when there is nothing to read: not an interface patch, or the
//  neighbour's material does not state a gamma. The caller then keeps its own
//  value rather than being handed a fabricated one.
Foam::scalar neighbourMaterialSEEC(const Foam::fvPatch& p)
{
    using namespace Foam;

    const auto* mppPtr = isA<mappedPatchBase>(p.patch());
    if (!mppPtr || !mppPtr->sameWorld())
    {
        return -1;
    }

    const polyMesh& nbrMesh = mppPtr->sampleMesh();

    if (!nbrMesh.foundObject<IOdictionary>("dielectricProperties"))
    {
        return -1;
    }

    const IOdictionary& props =
        nbrMesh.lookupObject<IOdictionary>("dielectricProperties");

    return props.getOrDefault<scalar>("gammaSEE", -1);
}

} // End anonymous namespace

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
electronDDWallFluxMixedFvPatchScalarField::
electronDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxMixedFvPatchScalarField(p, iF),
    enableSurfaceCharging_(false),
    includeDriftFlux_(false),
    enableSEE_(false),
    seeReport_(false),
    seeInertReported_(false),
    defaultSEEC_(0.0),
    speciesSEEC_(dictionary::null),
    emissionDict_(dictionary::null),
    electronReflection_(0.0),
    material_(word::null),
    emitted_(p.size(), Zero),
    emittedValid_(false),
    emission_(),
    emissionReported_(false),
    seec_(0),
    mapped_(false)
{}

// Dictionary Constructor
electronDDWallFluxMixedFvPatchScalarField::
electronDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ddWallFluxMixedFvPatchScalarField(p, iF, dict),
    enableSurfaceCharging_
            (dict.lookupOrDefault<bool>("enableSurfaceCharging", false)),
    includeDriftFlux_(dict.lookupOrDefault<bool>("includeDriftFlux", false)),
    enableSEE_(dict.lookupOrDefault<bool>("enableSEE", false)),
    seeReport_(dict.lookupOrDefault<bool>("seeReport", false)),
    seeInertReported_(false),
    // 0.001, not 0.05.
    //
    // gamma -- the probability that an ion striking a surface releases an
    // electron -- spans roughly 1e-3 for contaminated oxides and dielectric
    // barriers up to 1e-1 for clean metals in vacuum. 0.05 is a clean-metal
    // figure, and it is the wrong end of that range to default to: the surfaces
    // this solver is aimed at are air-exposed electrodes and dielectric
    // barriers, both contaminated.
    //
    // It matters because gamma is not a detail. In a barrier discharge it is
    // what makes the discharge SELF-SUSTAINING rather than a single avalanche,
    // so a default fifty times too high manufactures sustainment the case never
    // asked for -- and the run looks entirely healthy.
    //
    // Default changed 2026-09-02. Any case that relied on the old default
    // should state its own value; a case that already states one is unaffected.
    // PRECEDENCE: an explicit `defaultSEEC` wins; otherwise the neighbouring
    // region's MATERIAL; otherwise 0.001 (contaminated oxide / barrier).
    //
    // So a meshed barrier needs gamma stated ONCE, with its material, and the
    // interface picks it up -- while a hand-written case can still override per
    // patch, and a patch with no material behind it keeps the default.
    defaultSEEC_
    (
        dict.found("defaultSEEC")
      ? dict.get<scalar>("defaultSEEC")
      : (neighbourMaterialSEEC(p) > 0 ? neighbourMaterialSEEC(p) : 0.001)
    ),
    speciesSEEC_(dict.subOrEmptyDict("speciesSEEC")),
    emissionDict_(dict.subOrEmptyDict("emission")),
    electronReflection_(dict.getOrDefault<scalar>("electronReflection", 0.0)),
    material_(dict.getOrDefault<word>("material", word::null)),
    emitted_(p.size(), Zero),
    emittedValid_(false),
    emission_(),
    emissionReported_(false),
    seec_(0), 
    mapped_(false)
{}

// Mapping Constructor
electronDDWallFluxMixedFvPatchScalarField::
electronDDWallFluxMixedFvPatchScalarField
(
    const electronDDWallFluxMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ddWallFluxMixedFvPatchScalarField(ptf, p, iF, mapper),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_),
    includeDriftFlux_(ptf.includeDriftFlux_),
    enableSEE_(ptf.enableSEE_),
    seeReport_(ptf.seeReport_),
    seeInertReported_(false),
    defaultSEEC_(ptf.defaultSEEC_),
    speciesSEEC_(ptf.speciesSEEC_),
    emissionDict_(ptf.emissionDict_),
    electronReflection_(ptf.electronReflection_),
    material_(ptf.material_),
    emitted_(ptf.emitted_),
    emittedValid_(false),
    emission_(),
    emissionReported_(false),
    seec_(ptf.seec_),
    mapped_(ptf.mapped_)
{}

// Copy Constructor
electronDDWallFluxMixedFvPatchScalarField::
electronDDWallFluxMixedFvPatchScalarField
(
    const electronDDWallFluxMixedFvPatchScalarField& ptf
)
:
    ddWallFluxMixedFvPatchScalarField(ptf),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_),
    includeDriftFlux_(ptf.includeDriftFlux_),
    enableSEE_(ptf.enableSEE_),
    seeReport_(ptf.seeReport_),
    seeInertReported_(false),
    defaultSEEC_(ptf.defaultSEEC_),
    speciesSEEC_(ptf.speciesSEEC_),
    emissionDict_(ptf.emissionDict_),
    electronReflection_(ptf.electronReflection_),
    material_(ptf.material_),
    emitted_(ptf.emitted_),
    emittedValid_(false),
    emission_(),
    emissionReported_(false),
    seec_(ptf.seec_),
    mapped_(ptf.mapped_)
{}

// Copy Constructor (with new internal field)
electronDDWallFluxMixedFvPatchScalarField::
electronDDWallFluxMixedFvPatchScalarField
(
    const electronDDWallFluxMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxMixedFvPatchScalarField(ptf, iF),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_),
    includeDriftFlux_(ptf.includeDriftFlux_),
    enableSEE_(ptf.enableSEE_),
    seeReport_(ptf.seeReport_),
    seeInertReported_(false),
    defaultSEEC_(ptf.defaultSEEC_),
    speciesSEEC_(ptf.speciesSEEC_),
    emissionDict_(ptf.emissionDict_),
    electronReflection_(ptf.electronReflection_),
    material_(ptf.material_),
    emitted_(ptf.emitted_),
    emittedValid_(false),
    emission_(),
    emissionReported_(false),
    seec_(ptf.seec_),
    mapped_(ptf.mapped_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void electronDDWallFluxMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Initialize the SEE entries
    if (enableSEE_ && !mapped_)
    {
        if (!db().foundObject<plasmaTransport>("plasmaTransport"))
        {
            ddWallFluxMixedFvPatchScalarField::updateCoeffs();
            return;
        }

        const plasmaSpecies& speciesDB =
            db().lookupObject<plasmaTransport>("plasmaTransport").species();

        seec_.setSize(speciesDB.nSpecies(), defaultSEEC_);

        forAll(speciesDB.speciesNames(), specI)
        {
            const word& name = speciesDB.speciesNames()[specI];
            if (speciesSEEC_.found(name))
            {
                seec_[specI] =
                    speciesSEEC_.get<scalar>(name);
            }
        }

        mapped_ = true;
    }

    // ---- THE EMITTED FLUX, COMPUTED FIRST ---------------------------------
    //
    // Eq. (6.8) needs Gamma_w INSIDE the loss speed, and the loss speed is
    // built by the base class call below -- so it has to exist before that
    // call, not after. Stored in emitted_ and consumed by
    // calcEffectiveWallVelocity/calcAbsorptionVelocity.
    emitted_.setSize(this->patch().size());
    emitted_ = Zero;
    emittedValid_ = false;

    if (db().foundObject<plasmaTransport>("plasmaTransport"))
    {
        const plasmaTransport& tr0 =
            db().lookupObject<plasmaTransport>("plasmaTransport");

        const plasmaSpecies& sdb0 = tr0.species();

        const scalarField& magSf0 = this->patch().magSf();

        // ION-INDUCED, the built-in path.
        if (enableSEE_ && mapped_)
        {
            for (const label specI : sdb0.positiveIonSpeciesIDs())
            {
                const word fn("particleFlux_" + sdb0.speciesName(specI));

                if
                (
                    !this->patch().boundaryMesh().mesh()
                        .foundObject<surfaceScalarField>(fn)
                ) continue;

                const fvsPatchScalarField& phiI =
                    this->patch().lookupPatchField
                    <surfaceScalarField, scalar>(fn);

                emitted_ += seec_[specI]*max(scalar(0), phiI/magSf0);
            }
        }

        // THE EMISSION MODELS.
        if (!emissionDict_.empty())
        {
            if (emission_.empty())
            {
                emission_.setSize(emissionDict_.size());

                label k = 0;
                for (const word& type : emissionDict_.toc())
                {
                    // The patch's material is INHERITED: a model is handed
                    // only its own sub-dictionary, so a `material` declared
                    // beside `emission` must be merged in.
                    dictionary md(emissionDict_.subDict(type));

                    if (!material_.empty() && !md.found("material"))
                    {
                        md.add("material", material_);
                    }

                    emission_.set
                    (
                        k++, emissionModel::New(type, this->patch(), md)
                    );
                }
            }

            forAll(emission_, i)
            {
                emitted_ += emission_[i].emittedFlux();
            }

            if (!emissionReported_)
            {
                emissionReported_ = true;

                Info<< "    electronDDWallFluxMixed on patch `"
                    << this->patch().name() << "`: emission mechanisms" << nl;

                forAll(emission_, i)
                {
                    emission_[i].report(Info);
                }

                Info<< "            emitted flux at the FIRST evaluation, max"
                       " over the patch: " << gMax(emitted_) << " 1/m2/s"
                    << nl;

                if (gMax(emitted_) <= 0)
                {
                    Info<< "            (zero is EXPECTED here for a"
                           " flux-driven mechanism: the species fields" << nl
                        << "            are not solved yet. It is not evidence"
                           " that the mechanism is inert.)" << nl;
                }

                Info<< endl;
            }
        }

        emittedValid_ = true;
    }

    // Call the standard updateCoeffs from the base class
    ddWallFluxMixedFvPatchScalarField::updateCoeffs();

    // ---- THE SEPARATE -Gamma_w OF EQ. (6.1) --------------------------------
    //
    // This condition imposes n*(uDrift_n + uEff) = n*w_w. Eq. (6.1) is
    //     Gamma = n*w_w - Gamma_w
    // so the -Gamma_w must be applied here, in addition to the -Gamma_w/n
    // already inside w_w (see emitted_). Together, unclamped and at r = 0,
    // they give Gamma = n*thermal + n*drift - 2*Gamma_w -- eq. (6.6) closed
    // with (6.1).
    //
    // emitted_ was computed BEFORE the base call and carries BOTH the built-in
    // ion-SEE path and every emission model, so the two uses cannot disagree
    // about what Gamma_w is.
    if (emittedValid_ && gMax(emitted_) > 0)
    {
        const plasmaTransport& trE =
            db().lookupObject<plasmaTransport>("plasmaTransport");

        const plasmaSpecies& sdbE = trE.species();

        const plasmaTransportModel& bmE =
            trE.model(sdbE.speciesID(resolveSpeciesName()));

        const driftDiffusion& ddmE = refCast<const driftDiffusion>(bmE);

        const tmp<scalarField> tDe(this->patchDiffusivity(ddmE));
        const scalarField& De = tDe();

        const scalar DeFloor = SMALL*max(gMax(De), SMALL);

        forAll(emitted_, faceI)
        {
            if (De[faceI] > DeFloor)
            {
                this->refGrad()[faceI] += emitted_[faceI]/De[faceI];
            }
        }
    }

    if (!enableSEE_) return;

    // Shared with the base, so a field that is a PROPERTY of a species
    // rather than its density -- the LMEA energy density -- can say which
    // species it belongs to instead of having it guessed from its own name.
    const word speciesName = resolveSpeciesName();

    // Modify refValue or refGrad for secondary electron emission
    const fvPatch& p = patch();

    const scalarField& magSf = p.magSf();
    scalarField totalSEE(p.size(), 0.0);

    // Registry lookup for the plasmaTransport object
    const plasmaTransport& transport = 
                        db().lookupObject<plasmaTransport>("plasmaTransport");

    const plasmaSpecies& speciesDB = transport.species();

    // SEE NEEDS AN ION FLUX, AND AN IMMOBILE ION HAS NONE.
    //
    // `particleFlux_<species>` is registered only for mobileSpeciesIDs, and
    // `immobile::nEqn()` is `fvm::ddt(n)` alone -- no div, no laplacian -- so
    // an immobile ion has no transport equation, no wall flux, and no field
    // here to find. Every lookup below then fails, totalSEE stays 0, and
    // `enableSEE true` does NOTHING. Measured on the shipped streamer case,
    // whose ions are all `transportModel immobile`: 0 emitted at BOTH
    // electrodes. Announced once, because a setting that silently does nothing
    // is worse than one that is rejected.
    label nFlux = 0;
    for (const label specI : speciesDB.positiveIonSpeciesIDs())
    {
        if (p.boundaryMesh().mesh().foundObject<surfaceScalarField>
            ("particleFlux_" + speciesDB.speciesName(specI)))
        {
            ++nFlux;
        }
    }

    if (nFlux == 0 && !seeInertReported_)
    {
        seeInertReported_ = true;

        WarningInFunction
            << "`enableSEE true` on patch " << p.name()
            << ", but NO positive ion has a wall flux," << nl
            << "    so no secondary electrons can be emitted and the setting"
               " has no effect." << nl
            << "    Cause: the positive ions are `transportModel immobile`,"
               " which solves ddt only" << nl
            << "    (no transport, hence no wall flux). Secondary emission"
               " requires the ions to be" << nl
            << "    transported -- `driftDiffusion` -- which is a physics and"
               " cost decision, not a" << nl
            << "    switch. Set `enableSEE false` to say so deliberately."
            << endl;
    }

    for (const label specI : speciesDB.positiveIonSpeciesIDs())
    {
        const word fluxName =
            "particleFlux_" + speciesDB.speciesName(specI);

        if (!p.boundaryMesh().mesh().foundObject<surfaceScalarField>(fluxName))
            continue;

        const fvsPatchScalarField& phiI =
            p.lookupPatchField<surfaceScalarField, scalar>(fluxName);

        totalSEE += seec_[specI] * max(0.0, phiI / magSf);
    }

    // THE SAME diffusivity the wall flux above was built from -- through the
    // virtual, so the energy condition divides by D_eps and not by D_e. It
    // used to look up "D_<species>" directly, which handed the energy field
    // the electron's coefficient.
    const plasmaTransportModel& baseModelD =
        transport.model(speciesDB.speciesID(speciesName));
    const driftDiffusion& ddModelD = refCast<const driftDiffusion>(baseModelD);
    const tmp<scalarField> tDf(this->patchDiffusivity(ddModelD));
    const scalarField& Df = tDf();

    if (seeReport_)
    {
        // BOTH SIGNS of the raw flux are printed, because max(0,phi) is what
        // drives emission: a patch whose phi is negative everywhere MUST report
        // zero emission. Emission is expected at the CATHODE (positive ions
        // accelerated into it) and NOT at the anode (they are repelled), so a
        // non-zero anode total localises the defect instead of inferring it.
        Info<< "    SEE[" << p.name() << "] emitted "
            << gSum(totalSEE*magSf) << " 1/s (patch sum)" << endl;

        for (const label specI : speciesDB.positiveIonSpeciesIDs())
        {
            const word fluxName =
                "particleFlux_" + speciesDB.speciesName(specI);

            if (!p.boundaryMesh().mesh().foundObject<surfaceScalarField>
                 (fluxName))
            {
                continue;
            }

            const fvsPatchScalarField& phiR =
                p.lookupPatchField<surfaceScalarField, scalar>(fluxName);

            const scalarField fn(phiR/magSf);
            Info<< "      " << speciesDB.speciesName(specI)
                << ": phi/|Sf| min " << gMin(fn) << " max " << gMax(fn)
                << " ; positive part " << gSum(max(fn, scalar(0))*magSf)
                << "  (positive = INTO the wall, so it emits)" << endl;
        }
    }

    // Guarded on the diffusivity's OWN scale. D is zero until the transport
    // models are first corrected, and `+ VSMALL` (1e-300) turns that into an
    // overflow rather than a guard -- which is how this family used to die in
    // the first Poisson solve.
    const scalar Dfloor = SMALL*max(gMax(Df), SMALL);

    // NO refGradient CONTRIBUTION HERE ANY MORE. The ion-SEE flux is part of
    // emitted_, which was computed before the base call and applied above --
    // once. Adding it again here would double it, and the totalSEE computed in
    // this block now serves ONLY the `seeReport` diagnostics above.
    //
    // Kept as a guarded no-op rather than deleted so the diagnostic block, the
    // inert-SEE warning and the Dfloor reasoning above stay attached to the
    // code they describe.
    (void)Dfloor;
}

void electronDDWallFluxMixedFvPatchScalarField::write(Ostream& os) const
{
    ddWallFluxMixedFvPatchScalarField::write(os);   

    // ROUND-TRIP INVARIANT: write ALL of what read() accepts.
    os.writeEntry("electronReflection", electronReflection_);

    if (!material_.empty())
    {
        os.writeEntry("material", material_);
    }
    // ROUND-TRIP INVARIANT: write ALL of what read() accepts.
    //
    // `emission` was missing here, so the solver's own rewrite of a field
    // silently DROPPED every mechanism the case had configured -- which is
    // exactly how the first equivalence test came to compare a run WITH
    // emission against a run whose emission had been erased, and report
    // agreement. Same class of bug as thinDielectricPotential's, 2026-09-03.
    if (!emissionDict_.empty())
    {
        os.beginBlock("emission");
        emissionDict_.write(os, false);
        os.endBlock();
    }

    os.writeEntry("enableSurfaceCharging", enableSurfaceCharging_);
    os.writeEntry("includeDriftFlux", includeDriftFlux_);
    os.writeEntry("enableSEE", enableSEE_);
    os.writeEntry("defaultSEEC", defaultSEEC_);
    if (!speciesSEEC_.empty())
    {
        os.writeEntry("speciesSEEC", speciesSEEC_);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
