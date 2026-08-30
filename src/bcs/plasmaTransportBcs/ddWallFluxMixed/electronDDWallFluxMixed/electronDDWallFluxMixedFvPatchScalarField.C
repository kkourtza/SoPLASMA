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
electronDDWallFluxMixedFvPatchScalarField::calcAbsorptionVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal
) const
{
    tmp<scalarField> tVel = calcThermalVelocity(m, T);
    scalarField& uAbs = tVel.ref();

    // If drift flux is enabled, add the directed motion component
    if (includeDriftFlux_)
    {
        uAbs += max(0.0, uDriftNormal);
    }
    
    return tVel;
}

tmp<scalarField> 
electronDDWallFluxMixedFvPatchScalarField::calcEffectiveWallVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal
) const
{
    tmp<scalarField> tVel = calcThermalVelocity(m, T);
    scalarField& uWall = tVel.ref();

    // If drift flux is enabled, add the directed motion component
    if (includeDriftFlux_)
    {
        uWall += max(0.0, -uDriftNormal);
    }
    else
    {
        uWall -= uDriftNormal;
    }

    return tVel;
}

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
    defaultSEEC_(dict.lookupOrDefault<scalar>("defaultSEEC", 0.05)),
    speciesSEEC_(dict.subOrEmptyDict("speciesSEEC")),
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

    // Call the standard updateCoeffs from the base class
    ddWallFluxMixedFvPatchScalarField::updateCoeffs();

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

    forAll(Df, faceI)
    {
        if (Df[faceI] > Dfloor)
        {
            this->refGrad()[faceI] += totalSEE[faceI]/Df[faceI];
        }
    }
}

void electronDDWallFluxMixedFvPatchScalarField::write(Ostream& os) const
{
    ddWallFluxMixedFvPatchScalarField::write(os);   

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
