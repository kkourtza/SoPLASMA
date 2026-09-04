/*---------------------------------------------------------------------------*\
  File: energyDDWallFluxMixedFvPatchScalarField.C
  Part of: SoPLASMA
  Copyright (C) 2026
  License: GNU General Public License v3 or later
\*---------------------------------------------------------------------------*/

#include "energyDDWallFluxMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

namespace Foam
{

defineTypeNameAndDebug(energyDDWallFluxMixedFvPatchScalarField, 0);

addToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    energyDDWallFluxMixedFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

word energyDDWallFluxMixedFvPatchScalarField::resolveSpeciesName() const
{
    // An explicit `species` entry wins: a case may name its electron something
    // other than the registry's electronSpeciesID.
    const word base(ddWallFluxMixedFvPatchScalarField::resolveSpeciesName());
    if (base != this->internalField().name())
    {
        return base;                    // the user set `species`
    }

    // Otherwise the electron, by construction: this condition is only ever
    // applied to the electron-energy field. The base class would return the
    // FIELD name here -- "nEps_e" does not start with "n_" -- which is then
    // looked up as a species and fails.
    const auto& transport =
        db().lookupObject<plasmaTransport>(plasmaTransport::typeName);

    const plasmaSpecies& speciesDB = transport.species();
    return speciesDB.speciesName(speciesDB.electronSpeciesID());
}


tmp<scalarField>
energyDDWallFluxMixedFvPatchScalarField::normalisingDensity() const
{
    // "n_" + the electron species name -- DERIVED, so a case that names its
    // electron something unusual still resolves. See the header for the
    // measured defect this override fixes.
    const word nName("n_" + this->resolveSpeciesName());

    if (!this->patch().boundaryMesh().mesh().foundObject<volScalarField>(nName))
    {
        FatalErrorInFunction
            << "Electron density field `" << nName << "` not found on patch `"
            << this->patch().name() << "`." << nl
            << "The electron-energy wall condition needs n_e -- not n_eps --"
            << " to form Gamma_w/n_e in Hagelaar eqs (6.6) and (6.15)." << nl
            << exit(FatalError);
    }

    return tmp<scalarField>::New
    (
        this->patch().lookupPatchField<volScalarField, scalar>(nName)
    );
}


tmp<scalarField>
energyDDWallFluxMixedFvPatchScalarField::energyWeight
(
    const scalarField& gRatio
) const
{
    // eq. (6.15): eps_w/eps = 5/3 - (2/3)*(Gamma_w/n_e)/A.
    //
    // The 2/3 is the T_e -> eps conversion, eps = (3/2) T_e, which is the
    // SAME Maxwellian assumption that fluxEnergyFactor_'s 5/3 default carries.
    // A user who overrides the factor for a non-Maxwellian EEDF is overriding
    // only the zero-creation value; the shape of the correction is (6.15)'s.
    tmp<scalarField> tW = tmp<scalarField>::New(gRatio.size(), Zero);
    scalarField& w = tW.ref();

    forAll(w, faceI)
    {
        w[faceI] = hagelaarEnergyWeight(fluxEnergyFactor_, gRatio[faceI]);
    }

    return tW;
}


tmp<scalarField>
energyDDWallFluxMixedFvPatchScalarField::calcAbsorptionVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal
) const
{
    // The electron condition's TOTAL loss speed W = (1-r) w_w, weighted by
    // eq. (6.15). The energy leaves on the SAME electrons, so reusing W keeps
    // the two conditions consistent by construction: any change to the
    // thermal base, the drift-flux option, reflection, emission or the surface
    // treatment is inherited rather than duplicated here.
    scalarField gRatio;
    tmp<scalarField> tVel =
        electronDDWallFluxMixedFvPatchScalarField::wallLossSpeed
        (
            m, T, uDriftNormal, gRatio
        );

    tVel.ref() *= energyWeight(gRatio)();

    return tVel;
}


tmp<scalarField>
energyDDWallFluxMixedFvPatchScalarField::calcEffectiveWallVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal
) const
{
    // uEff = W_eps - uDrift_n, the same bookkeeping identity as the electron
    // condition. uDriftNormal here is the ENERGY's drift velocity, which is
    // what the mixed condition's n_p*(uDrift_n + uEff) requires.
    scalarField gRatio;
    tmp<scalarField> tVel =
        electronDDWallFluxMixedFvPatchScalarField::wallLossSpeed
        (
            m, T, uDriftNormal, gRatio
        );

    tVel.ref() *= energyWeight(gRatio)();
    tVel.ref() -= uDriftNormal;

    return tVel;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

energyDDWallFluxMixedFvPatchScalarField::
energyDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    electronDDWallFluxMixedFvPatchScalarField(p, iF),
    fluxEnergyFactor_(5.0/3.0),
    secondaryElectronEnergy_(0.0)
{}


energyDDWallFluxMixedFvPatchScalarField::
energyDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    electronDDWallFluxMixedFvPatchScalarField(p, iF, dict),
    // 4/3, not 5/2 and not 1 -- see the header. 5/2 is the INFLOW value for a
    // Maxwellian reservoir and would over-drain at every wall; 1 is the
    // population mean and would under-drain, because faster electrons reach
    // the wall more often.
    fluxEnergyFactor_(dict.getOrDefault<scalar>("fluxEnergyFactor", 5.0/3.0)),
    // USER-DEFINED, with a default rather than a derivation, because there is
    // no rigid calculation available: the birth energy of a secondary electron
    // depends on the surface material and work function and on the incident
    // ion energy, none of which this solver models. 2 eV is a conventional
    // value for metals and dielectrics in the 1-5 eV range reported in the
    // literature. Anything claiming to compute it from first principles here
    // would be inventing the number, not deriving it.
    secondaryElectronEnergy_
    (
        dict.lookupOrDefault<scalar>("secondaryElectronEnergy", 2.0)
    )
{}


energyDDWallFluxMixedFvPatchScalarField::
energyDDWallFluxMixedFvPatchScalarField
(
    const energyDDWallFluxMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    electronDDWallFluxMixedFvPatchScalarField(ptf, p, iF, mapper),
    fluxEnergyFactor_(ptf.fluxEnergyFactor_),
    secondaryElectronEnergy_(ptf.secondaryElectronEnergy_)
{}


energyDDWallFluxMixedFvPatchScalarField::
energyDDWallFluxMixedFvPatchScalarField
(
    const energyDDWallFluxMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    electronDDWallFluxMixedFvPatchScalarField(ptf, iF),
    fluxEnergyFactor_(ptf.fluxEnergyFactor_),
    secondaryElectronEnergy_(ptf.secondaryElectronEnergy_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void energyDDWallFluxMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // The secondary electrons that enter here are the SAME ones the electron
    // condition already counted; each carries secondaryElectronEnergy_, so the
    // energy inflow is simply that inflow times that energy.
    //
    // WAS A DELTA, AND THE DELTA WAS WRONG. The previous version captured
    // refGrad BEFORE calling the parent and treated it as a baseline:
    //
    //     grad0 = refGrad;  parent();  refGrad = grad0 + s*(refGrad - grad0);
    //
    // That assumes the parent ACCUMULATES onto refGrad. It does not -- the
    // base `ddWallFluxMixed::updateCoeffs` RESETS `refGrad() = 0.0` before
    // adding the emission term. So grad0 was never a baseline; it was the
    // value left over from the PREVIOUS corrector, which the parent had
    // already discarded. The recursion that produced,
    //
    //     x_{n+1} = s*SEE_{n+1} + (1 - s)*x_n,
    //
    // does not settle on s*SEE at all: at the shipped s = 2 it alternates
    // sign every corrector and never converges. Nothing downstream can damp a
    // boundary condition that flips sign on each pass.
    //
    // After the parent returns, refGrad IS the secondary-electron particle
    // inflow (the base zeroed it first), so the energy version is one
    // multiplication. The ion loop, the per-species coefficients and the
    // diffusivity lookup are still inherited, not duplicated.
    // The parent has now built the wall flux from the ENERGY's own mu_eps and
    // D_eps (see patchMobility/patchDiffusivity below), and its refGrad is the
    // secondary-electron inflow expressed against D_eps. Those are the same
    // electrons the electron condition emitted; each carries
    // secondaryElectronEnergy_, so the energy inflow is that inflow times that
    // energy -- one multiplication, no diffusivity ratio anywhere.
    electronDDWallFluxMixedFvPatchScalarField::updateCoeffs();

    this->refGrad() *= secondaryElectronEnergy_;
}


tmp<scalarField> energyDDWallFluxMixedFvPatchScalarField::patchMobility
(
    const driftDiffusion& ddModel
) const
{
    return energyCoefficient("muEps", ddModel, true);
}


tmp<scalarField> energyDDWallFluxMixedFvPatchScalarField::patchDiffusivity
(
    const driftDiffusion& ddModel
) const
{
    return energyCoefficient("DEps", ddModel, false);
}


//- The energy's own coefficient, published by localEnergyEnergyModel from the
//  same expression its laplacian/div use. NOT derived from the electron
//  coefficient by a factor: the 5/3 is the Maxwellian limit, while these are
//  tabulated from the EEDF (the model applies a factor of 1 in that case).
tmp<scalarField> energyDDWallFluxMixedFvPatchScalarField::energyCoefficient
(
    const word& fieldName,
    const driftDiffusion& ddModel,
    const bool isMobility
) const
{
    if (!this->db().foundObject<volScalarField>(fieldName))
    {
        FatalErrorInFunction
            << "energyDDWallFluxMixed needs `" << fieldName << "`, which"
               " localEnergyEnergyModel registers." << nl
            << "    It was not found: this condition is only meaningful with"
               " an LMEA electron-energy" << nl
            << "    equation enabled." << nl
            << exit(FatalError);
    }

    return tmp<scalarField>::New
    (
        this->db().lookupObject<volScalarField>(fieldName)
            .boundaryField()[this->patch().index()]
    );
}


void energyDDWallFluxMixedFvPatchScalarField::write(Ostream& os) const
{
    electronDDWallFluxMixedFvPatchScalarField::write(os);

    os.writeEntry("fluxEnergyFactor", fluxEnergyFactor_);
    os.writeEntry("secondaryElectronEnergy", secondaryElectronEnergy_);
}


} // End namespace Foam

// ************************************************************************* //
