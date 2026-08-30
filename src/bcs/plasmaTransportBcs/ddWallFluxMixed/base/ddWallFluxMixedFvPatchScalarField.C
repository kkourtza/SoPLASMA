/*---------------------------------------------------------------------------*\
  File: ddWallFluxMixedFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::ddWallFluxMixedFvPatchScalarField.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvPatchFieldMapper.H"

#include "ddWallFluxMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(ddWallFluxMixedFvPatchScalarField, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

//- Thermal velocity for a single constant temperature (dimensionedScalar)
dimensionedScalar ddWallFluxMixedFvPatchScalarField::calcThermalVelocity
(
    const dimensionedScalar& m,
    const dimensionedScalar& T
) const
{
    return 0.25 * sqrt 
    (
        (8.0 * constant::plasma::kappaBoltzmann * T)
        /
        (constant::mathematical::pi * m)
    );
}

//- Thermal velocity for a field temperature (scalarField)
tmp<scalarField> ddWallFluxMixedFvPatchScalarField::calcThermalVelocity
(
    const dimensionedScalar& m,
    const scalarField& T
) const
{
    return 0.25 * sqrt 
    (
        (8.0 * constant::plasma::kappaBoltzmann.value() * T)
        /
        (constant::mathematical::pi * m.value())
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
ddWallFluxMixedFvPatchScalarField::ddWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("none"),
    speciesNameOverride_(word::null),
    TValue_("T", dimTemperature, 300.0)
{
    this->refValue()      = 0.0;
    this->refGrad()       = 0.0;
    this->valueFraction() = 0.0;
}

// Dictionary Constructor
ddWallFluxMixedFvPatchScalarField::ddWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("none"),
    speciesNameOverride_(dict.lookupOrDefault<word>("species", word::null)),
    TValue_("T", dimTemperature, 300.0)
{
    // THE WALL TEMPERATURE. Three ways in, and exactly one applies:
    //
    //   T     <fieldName>   -- follow a solved field (T_e under LMEA)
    //   TeV   <electronVolts> -- a FIXED value, in the unit electron
    //                            temperatures are actually quoted in
    //   (neither)           -- the condition's own default; 1 eV for electrons
    //
    // `TeV` exists because requiring Kelvin here made the user do the
    // conversion, and the two conventions in circulation (eps = 3/2 kT vs the
    // temperature equivalent of eps) differ by 1.5. `TeV 1` means kT = 1 eV,
    // which is what "Te = 1 eV" means everywhere in the literature.
    const bool haveT = dict.found("T");
    const bool haveTeV = dict.found("TeV");

    if (haveT && haveTeV)
    {
        FatalIOErrorInFunction(dict)
            << "both `T` and `TeV` are set. They are the same control." << nl
            << "    Use `T <fieldName>` to follow a field, or `TeV <eV>` for a"
               " fixed value." << exit(FatalIOError);
    }

    if (haveTeV)
    {
        const scalar TeV = dict.get<scalar>("TeV");

        if (TeV <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "TeV must be positive; got " << TeV << "." << nl
                << "    It is a wall temperature in eV (kT), not an energy"
                   " offset." << exit(FatalIOError);
        }

        // kT = TeV [eV]  =>  T = TeV * e / k_B  (1 eV = 11604.5 K)
        TValue_.value() =
            TeV*constant::plasma::eCharge.value()
               /constant::plasma::kappaBoltzmann.value();
        TName_ = "none";
    }
    else if (haveT)
    {
        const entry& e = dict.lookupEntry("T", keyType::LITERAL);
        ITstream& is = e.stream();

        if (is.peek().isWord())
        {
            is >> TName_;
        }
        else if (is.peek().isNumber())
        {
            is >> TValue_.value();
            TName_ = "none";
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Entry 'T' must be a word (field name) or a scalar in"
                   " KELVIN." << nl
                << "    For a fixed electron temperature use `TeV <eV>`"
                   " instead." << exit(FatalIOError);
        }
    }
    else
    {
        // No field, no explicit value: the condition's own default rather
        // than a hard error, so a case that has no electron-temperature field
        // at all (LFA runs the `gasTemperature` energy model and registers
        // none) still starts.
        TValue_ = this->defaultTValue();
        TName_ = "none";
    }

    if (this->readMixedEntries(dict))
    {
        // Full restart or values provided in dictionary
    }
    else
    {
        this->refValue()      = 0.0;
        this->refGrad()       = 0.0;
        this->valueFraction() = 0.0;
    }

    // INITIALISE THE PATCH VALUE. Without this the boundary field is the
    // uninitialised memory handed out by `mixedFvPatchScalarField(p, iF)`:
    // that constructor allocates the patch Field and does NOT set it, and
    // nothing here read the `value` entry, so a case supplying
    // `value $internalField;` was silently ignored.
    //
    // The consequence is not subtle. The first assembly that touches this
    // field -- chargeDensity, before any species is solved -- picks up that
    // garbage, and the Poisson solve diverges to nan (measured: ePotential
    // initial residual 1, final nan, 2000 iterations) or trips SIGFPE. It
    // looks like a physics blow-up and is not one.
    //
    // The IMPLICIT sibling (ddWallFluxImplicit) never had this: it passes the
    // dictionary down as `fvPatchScalarField(p, iF, dict)`, which reads
    // `value`. This is the same contract, restored for the mixed family.
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        // No `value` given: the internal field is the only defensible
        // starting point, and is what OpenFOAM's own mixed BC falls back to.
        fvPatchScalarField::operator=(this->patchInternalField());
    }
}

// Mapping Constructor
ddWallFluxMixedFvPatchScalarField::ddWallFluxMixedFvPatchScalarField
(
    const ddWallFluxMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    speciesNameOverride_(ptf.speciesNameOverride_),
    TValue_(ptf.TValue_)
{}

// Copy Constructor (from another patch field)
ddWallFluxMixedFvPatchScalarField::ddWallFluxMixedFvPatchScalarField
(
    const ddWallFluxMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    speciesNameOverride_(ptf.speciesNameOverride_),
    TValue_(ptf.TValue_)
{}

// Copy Constructor (from patch field and new internal field)
ddWallFluxMixedFvPatchScalarField::ddWallFluxMixedFvPatchScalarField
(
    const ddWallFluxMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    speciesNameOverride_(ptf.speciesNameOverride_),
    TValue_(ptf.TValue_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word
Foam::ddWallFluxMixedFvPatchScalarField::resolveSpeciesName() const
{
    if (!speciesNameOverride_.empty())
    {
        return speciesNameOverride_;
    }

    word name = this->internalField().name();
    if (name.startsWith("n_"))
    {
        name.erase(0, 2);
    }
    return name;
}


void ddWallFluxMixedFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const fvPatch& p = patch();
    if (p.size() == 0)
    {
        mixedFvPatchScalarField::updateCoeffs();
        return;
    }

    // Access the normal vector and delta coeffs
    const scalarField& delta = p.deltaCoeffs();

    // Determine species name (e.g., n_e -> e), or take the override.
    const word speciesName = resolveSpeciesName();

    // Lookup Transport Registry and Species Data
    if (!db().foundObject<plasmaTransport>("plasmaTransport"))
    {
        FatalErrorInFunction
            << "plasmaTransport not found in registry." << nl
            << exit(FatalError);
    }
    const plasmaTransport& transport =
                          db().lookupObject<plasmaTransport>("plasmaTransport");

    const plasmaSpecies& speciesDB = transport.species();

    const label speciesID = speciesDB.speciesID(speciesName);
    const dimensionedScalar& m = speciesDB.speciesMass(speciesID);
    const scalar Z = speciesDB.speciesChargeNumber(speciesID);

    // Access the drift-diffusion model
    const plasmaTransportModel& baseModel = transport.model(speciesID);

    if (!isA<driftDiffusion>(baseModel))
    {
        FatalErrorInFunction
            << "Species '" << speciesName << "' must use the driftDiffusion "
            << "transport model for this boundary condition." << nl
            << "Current model: " << baseModel.type() << nl
            << exit(FatalError);
    }

    const driftDiffusion& ddModel = refCast<const driftDiffusion>(baseModel);

    // Access the patch mobility, diffusivity and electric field.
    // Through the virtuals, so the energy condition can substitute its own.
    const tmp<scalarField> tmuf(this->patchMobility(ddModel));
    const tmp<scalarField> tDf(this->patchDiffusivity(ddModel));
    const scalarField& muf = tmuf();
    const scalarField& Df = tDf();
    const surfaceScalarField& phiE =
            p.boundaryMesh().mesh().lookupObject<surfaceScalarField>("phiE");

    const scalarField& phiEp = phiE.boundaryField()[p.index()];
    const scalarField uDrift_n(Z * muf * (phiEp / p.magSf()));

    // Physics Calculations
    tmp<scalarField> tUEff;
    tmp<scalarField> tUAbs;

    if (TName_ == "none")
    {
        const scalarField Tconst(p.size(), TValue_.value());
        
        tUEff = this->calcEffectiveWallVelocity(m, Tconst, uDrift_n);
        tUAbs = this->calcAbsorptionVelocity(m, Tconst, uDrift_n);
    }
    else
    {
        const auto* TPtr = db().findObject<volScalarField>(TName_);

        if (!TPtr)
        {
            FatalErrorInFunction
                << "Temperature field '" << TName_ << "' not found in registry." << nl
                << "Either set T to a constant scalar value (e.g., T 300;)," << nl
                << "or ensure the volScalarField exists in the mesh registry." << nl
                << exit(FatalError);
        }

        const scalarField& TField = TPtr->boundaryField()[p.index()];

        tUEff = this->calcEffectiveWallVelocity(m, TField, uDrift_n);
        tUAbs = this->calcAbsorptionVelocity(m, TField, uDrift_n);
    }

    const scalarField& uEff = tUEff();
    const scalarField& uAbs = tUAbs();

    // Set mixed b.c. parameters
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->operator==(this->patchInternalField());
    scalarField& f = this->valueFraction();

    // Set the D/Δ 
    const scalarField D_delta(Df * delta);
    const word scheme = ddModel.fluxScheme();

    if (scheme == "ScharfetterGummel")
    {
        const auto Bern = [](scalar x) -> scalar
        {
            const scalar ax = mag(x);
            if (ax < 1e-4) return 1.0 - 0.5*x + (x*x)/12.0;
            if (x > 100.0)  return 0.0;
            if (x < -100.0) return -x;
            return x / (Foam::exp(x) - 1.0);
        };

        forAll(p, faceI)
        {
            const scalar Pe = uDrift_n[faceI] / (D_delta[faceI] + VSMALL);
            const scalar num = uEff[faceI];
            const scalar den = D_delta[faceI] * Bern(Pe) + uAbs[faceI];
            f[faceI] = num / (den + VSMALL);
        }
    }
    else // Standard schemes (e.g. upwind(div) + central(laplacian))
    {
        // WELL-POSEDNESS OF THE STANDARD BRANCH.
        //
        // The condition imposes a TOTAL wall flux n_p*(uDrift_n + uEff) and
        // sets n_p = (1 - f)*n_c with f = uEff/(D/delta + uEff). Under
        // `includeDriftFlux false` -- the electron default -- uEff is
        // u_th - uDrift_n, chosen so the drift the convective term already
        // carries cancels and the total flux is the THERMAL flux alone,
        // Gamma = (1/4) n v_th. That subtraction is deliberate, not a defect.
        //
        // Its consequence is that uEff goes NEGATIVE once the drift into the
        // wall exceeds the thermal speed. That much is still self-consistent:
        // holding the total flux at (1/4) n v_th while the field drives more
        // than that into the wall requires diffusion to carry the surplus back
        // out, which needs n_p > n_c. The wall value legitimately rises above
        // the interior.
        //
        // What is NOT self-consistent is the denominator reaching zero, at
        //     uDrift_n = u_th + D/delta,
        // where f diverges and changes sign: the boundary value becomes
        // unbounded and the discretisation has inverted. Nothing detected
        // that -- `+ SMALL` turns the singularity into an overflow rather
        // than a guard, exactly as `+ VSMALL` did on the diffusivity.
        //
        // MEASURED margin on the positive-streamer anode cell at 0.5 ns:
        // u_th 1.0225e5, uDrift_n 1.0113e5, D/delta 4.44e4 m/s, so the
        // denominator is 4.55e4 -- positive, but only because D/delta is
        // carrying it. A stronger field or a coarser near-wall cell closes it.
        //
        // The ScharfetterGummel branch above cannot reach this: its
        // denominator is D/delta*Bern(Pe) + uAbs, and both terms are
        // non-negative with uAbs >= u_th > 0.
        //
        // Reported before it bites, fatal once it has. Auto-correcting an
        // inverted boundary condition would be inventing a wall flux nobody
        // asked for; the two real remedies are named in the message.
        scalar minDenFrac = GREAT;

        forAll(p, faceI)
        {
            const scalar den = D_delta[faceI] + uEff[faceI];

            minDenFrac = min(minDenFrac, den/(D_delta[faceI] + VSMALL));

            f[faceI] = uEff[faceI] / (den + SMALL);
        }

        reduce(minDenFrac, minOp<scalar>());

        if (minDenFrac <= 0)
        {
            FatalErrorInFunction
                << "wall-flux condition on patch " << p.name()
                << " for field " << this->internalField().name()
                << " has become singular." << nl << nl
                << "    The drift into the wall now exceeds the thermal speed"
                   " by more than the" << nl
                << "    near-wall diffusive velocity D/delta, so"
                   " D/delta + uEff <= 0 and the boundary" << nl
                << "    value it implies is unbounded. Worst face:"
                   " (D/delta + uEff)/(D/delta) = "
                << minDenFrac << "." << nl << nl
                << "    Two remedies, both physical:" << nl
                << "      * `includeDriftFlux true` on this patch -- the"
                   " Hagelaar & Kroesen (2000) wall" << nl
                << "        flux Gamma = (1/4) n v_th + mu n E, whose uEff is"
                   " never below u_th, so this" << nl
                << "        failure mode does not exist. The present default,"
                   " `false`, imposes the" << nl
                << "        thermal flux alone and is what allows uEff to go"
                   " negative." << nl
                << "      * refine the near-wall cell: D/delta grows as the"
                   " cell shrinks, which is the" << nl
                << "        term holding the denominator open."
                << exit(FatalError);
        }

        if (minDenFrac < 0.1 && !nearSingularReported_)
        {
            nearSingularReported_ = true;

            WarningInFunction
                << "wall-flux condition on patch " << p.name()
                << " for field " << this->internalField().name()
                << " is approaching its singularity:" << nl
                << "    (D/delta + uEff)/(D/delta) = " << minDenFrac
                << " at the worst face (fatal at 0)." << nl
                << "    The drift into the wall is close to the thermal"
                   " speed. See `includeDriftFlux`" << nl
                << "    and near-wall refinement. Reported once."
                << endl;
        }
    }

    mixedFvPatchField<scalar>::updateCoeffs();
}

dimensionedScalar ddWallFluxMixedFvPatchScalarField::defaultTValue() const
{
    // The base has no business guessing a temperature for an arbitrary
    // species: an ion or neutral wall is at the GAS temperature, an electron
    // wall is not, and silently picking one would be wrong half the time.
    FatalErrorInFunction
        << "no wall temperature given for patch " << patch().name()
        << " on field " << this->internalField().name() << "." << nl
        << "    Set `T <fieldName>` to follow a field, or `TeV <eV>` for a"
           " fixed value." << nl
        << "    (Only the electron conditions carry a default.)" << nl
        << exit(FatalError);

    return dimensionedScalar("T", dimTemperature, 0.0);
}


tmp<scalarField> ddWallFluxMixedFvPatchScalarField::patchMobility
(
    const driftDiffusion& ddModel
) const
{
    return tmp<scalarField>::New(ddModel.mobility().muPatch(patch().index()));
}


tmp<scalarField> ddWallFluxMixedFvPatchScalarField::patchDiffusivity
(
    const driftDiffusion& ddModel
) const
{
    return tmp<scalarField>::New
    (
        ddModel.diffusivity().DPatch(patch().index())
    );
}


void ddWallFluxMixedFvPatchScalarField::write(Ostream& os) const
{
    // Write standard Mixed BC entries (valueFraction, refValue, etc.)
    mixedFvPatchScalarField::write(os);

    // Write our custom entries so the simulation can be restarted
    if (TName_ != "none")
    {
        os.writeEntry("T", TName_);
    }
    else
    {
        os.writeEntry("T", TValue_.value());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
