/*---------------------------------------------------------------------------*\
  File: ddWallFluxImplicitFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::ddWallFluxImplicitFvPatchScalarField.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvPatchFieldMapper.H"

#include "ddWallFluxImplicitFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(ddWallFluxImplicitFvPatchScalarField, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

//- THE CONSISTENT WALL LOSS SPEED, Hagelaar HDR chapter 6.
//
//  This returns eq. (6.6)'s thermal term, vT/sqrt(pi) with vT = sqrt(2eT/m),
//  written here as 0.5*sqrt(8 k T/(pi m)) so the mean-speed form is visible.
//
//  IT IS TWICE WHAT THIS FUNCTION USED TO RETURN, and the factor is not a
//  correction bolted on -- it is the difference between two derivations:
//
//    eq. (6.3)  w_w = (1/2) vT/sqrt(pi) = (1/4) n <v> / n
//               The half-space integral of a MAXWELLIAN. Hagelaar calls this
//               "the simplest approach" and notes it "does not account for the
//               effects of the electric field and particle density gradient
//               and gives a bad description in case of significant directed
//               motion". This is what the code computed before 2026-09-04.
//
//    eq. (6.6)  w_w = max( vT/sqrt(pi) - Gamma_w/n , 0 )
//               The SHIFTED-Maxwellian treatment, integrated over the
//               half-space and closed self-consistently with eq. (6.1). In
//               Hagelaar's words: "without reflection or wall creation the
//               effective loss speed is TWICE AS LARGE as (6.3), but it is
//               reduced as Gamma_w increases."
//
//  WHY TWICE: the one-way flux from a FULL Maxwellian of density n is
//  (1/4) n <v>. At an absorbing wall the outgoing half of velocity space is
//  empty, so the density the fluid solves for is only the inward half -- half
//  of the full-Maxwellian density that produces that flux. Expressed in terms
//  of the ACTUAL wall density, the flux is therefore (1/2) n <v>.
//
//  Hagelaar's eq. (6.8) -- thermal-and-creation clamped at zero, plus a
//  separately clamped drift term -- "has consistent limits for all particle
//  species and is RECOMMENDED INSTEAD OF EQUATION (6.4)", and (6.4) is the
//  form this family implemented. Hence the change, and hence that it applies
//  to EVERY species and to the energy condition, not only to electrons: the
//  derivation is purely kinematic, with no charge or mass in it.
//
//  STILL TO DO, and deliberately not done in this step: the -Gamma_w/n term
//  and its max(...,0) clamp, which move secondary emission INSIDE the loss
//  speed. Emission is currently still added separately to refGradient. See the
//  note in electronDDWallFluxMixed.
//
//  Source: Hagelaar, G. J. M., HDR thesis, chapter 6, eqs. (6.1)-(6.8);
//  Literature/hdr-hagelaar.pdf. Footnote 30 there notes that many of the
//  chapter's equations "are not standard".
//- Thermal velocity for a single constant temperature (dimensionedScalar)
dimensionedScalar ddWallFluxImplicitFvPatchScalarField::calcThermalVelocity
(
    const dimensionedScalar& m,
    const dimensionedScalar& T
) const
{
    return 0.5 * sqrt 
    (
        (8.0 * constant::plasma::kappaBoltzmann * T)
        /
        (constant::mathematical::pi * m)
    );
}

//- Thermal velocity for a field temperature (scalarField)
tmp<scalarField> ddWallFluxImplicitFvPatchScalarField::calcThermalVelocity
(
    const dimensionedScalar& m,
    const scalarField& T
) const
{
    return 0.5 * sqrt 
    (
        (8.0 * constant::plasma::kappaBoltzmann.value() * T)
        /
        (constant::mathematical::pi * m.value())
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchScalarField(p, iF),
    TName_("none"),
    TValue_("T", dimTemperature, 300.0)
{}

// Dictionary Constructor
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchScalarField(p, iF, dict),
    TName_("none"),
    TValue_("T", dimTemperature, 300.0)
{
    // Record WHAT THE CASE STATED before anything defaults. See suppliedKeys_
    // in plasmaWallBC for why a defaulted value must not be written back.
    suppliedKeys_ = dict.toc();

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
            << "Entry 'T' must be either a word (field name) or a scalar value."
            << exit(FatalIOError);
    }
}

// Mapping Constructor
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const ddWallFluxImplicitFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    TValue_(ptf.TValue_)
{
    // A cloned/mapped field inherits WHAT THE CASE STATED, so a defaulted
    // value still is not written back. See plasmaWallBC::suppliedKeys_.
    suppliedKeys_ = ptf.suppliedKeys_;
}

// Copy Constructor (from another patch field)
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const ddWallFluxImplicitFvPatchScalarField& ptf
)
:
    fvPatchScalarField(ptf),
    TName_(ptf.TName_),
    TValue_(ptf.TValue_)
{
    // A cloned/mapped field inherits WHAT THE CASE STATED, so a defaulted
    // value still is not written back. See plasmaWallBC::suppliedKeys_.
    suppliedKeys_ = ptf.suppliedKeys_;
}

// Copy Constructor (from patch field and new internal field)
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const ddWallFluxImplicitFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    TValue_(ptf.TValue_)
{
    // A cloned/mapped field inherits WHAT THE CASE STATED, so a defaulted
    // value still is not written back. See plasmaWallBC::suppliedKeys_.
    suppliedKeys_ = ptf.suppliedKeys_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ddWallFluxImplicitFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Set n_face = n_cell (consistent with valueInternalCoeffs = 1 for drift)
    // The actual flux is controlled by the coefficient methods, not n_face.
    // This ensures post-processing and fvc operators see a sensible face value.
    fvPatchScalarField::operator==(this->patchInternalField());

    fvPatchField<scalar>::evaluate();
}

void ddWallFluxImplicitFvPatchScalarField::updateCoeffs()
{
    // Coefficients are computed on-the-fly in the coefficient methods.
    // Nothing to precompute here.
    if (updated()) return;
    fvPatchScalarField::updateCoeffs();
}

tmp<Field<scalar>>
ddWallFluxImplicitFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    const fvPatch& p = patch();

    // Access the normal vector and delta coeffs
    const vectorField nf = p.nf();

    // Determine species name (e.g., n_e -> e)
    const word fieldName = this->internalField().name();
    word speciesName = fieldName;
    if (speciesName.startsWith("n_"))
    {
        speciesName.erase(0, 2);
    }

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

    // Access the patch mobility, diffusivity and electric field
    const driftDiffusion& ddModel = refCast<const driftDiffusion>(baseModel);

    const scalarField& muf = ddModel.mobility().muPatch(p.index());
    const surfaceScalarField& phiE =
            p.boundaryMesh().mesh().lookupObject<surfaceScalarField>("phiE");

    const scalarField& phiEp = phiE.boundaryField()[p.index()];
    const scalarField Ef(phiEp / p.magSf());

    // Physics Calculations
    word scheme = ddModel.fluxScheme();
    const scalarField uDrift_n(Z * muf * Ef);

    return tmp<Field<scalar>>(new scalarField(pos(uDrift_n)));
}

tmp<Field<scalar>>
ddWallFluxImplicitFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // No explicit RHS contribution from drift
    return tmp<Field<scalar>>
    (
        new scalarField(this->size(), Zero)
    );
}

tmp<Field<scalar>>
ddWallFluxImplicitFvPatchScalarField::gradientInternalCoeffs() const
{
    const fvPatch& p = patch();

    // Determine species name (e.g., n_e -> e)
    const word fieldName = this->internalField().name();
    word speciesName = fieldName;
    if (speciesName.startsWith("n_"))
    {
        speciesName.erase(0, 2);
    }

    // Lookup Temperature (T)
    tmp<scalarField> tT;
    if (TName_ == "none")
    {
        tT = tmp<scalarField>::New(p.size(), TValue_.value());
    }
    else if (db().foundObject<volScalarField>(TName_))
    {
        tT = p.lookupPatchField<volScalarField, scalar>(TName_);
    }
    else
    {
        FatalErrorInFunction
            << "Temperature field '" << TName_ << "' not found in registry." 
            << nl << "Either set TName to 'none' and provide a constant TValue,"
            << " or ensure the field exists in the mesh registry." << nl
            << exit(FatalError);
    }
    const scalarField& T = tT();

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

    // Access the patch diffusivity
    const scalarField& Df  = ddModel.diffusivity().DPatch(p.index());

    // Physics Calculations
    tmp<scalarField> tUth = this->calcWallThermalVelocity(m, T);
    const scalarField& uth4 = tUth();

    // Implicit thermal flux: gradCoeff = -u_th/4 / D
    // (negative so fvm::laplacian adds positive diagonal contribution)
    return tmp<Field<scalar>>
    (
        new scalarField(-uth4 / (Df + VSMALL))
    );
}

tmp<Field<scalar>>
ddWallFluxImplicitFvPatchScalarField::gradientBoundaryCoeffs() const
{
    // No explicit RHS contribution from diffusion
    return tmp<Field<scalar>>
    (
        new scalarField(this->size(), Zero)
    );
}

void ddWallFluxImplicitFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    // Write our custom entries so the simulation can be restarted
    if (TName_ != "none")
    {
        os.writeEntry("T", TName_);
    }
    else
    {
        os.writeEntry("T", TValue_.value());
    }

    fvPatchScalarField::writeValueEntry(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
