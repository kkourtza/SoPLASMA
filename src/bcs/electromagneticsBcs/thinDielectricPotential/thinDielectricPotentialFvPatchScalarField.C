/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "thinDielectricPotentialFvPatchScalarField.H"
#include "electromagneticsModel.H"
#include "SoPLASMAConstants.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar thinDielectricPotentialFvPatchScalarField::layerCapacitance() const
{
    // FREE-STANDING: no conductor behind the layer, so there is nothing to be a
    // capacitor to and C is exactly zero. Signalled by a negative thickness,
    // which is what the constructor stores when `thickness` is absent -- as
    // opposed to a defaulted number, which would invent a capacitance.
    if (thickness_ <= 0)
    {
        return 0;
    }

    return constant::plasma::epsilon0.value()*epsilonR_/thickness_;
}


tmp<scalarField> thinDielectricPotentialFvPatchScalarField::surfaceCharge() const
{
    auto tsigma = tmp<scalarField>::New(patch().size(), Zero);

    if (surfChargeName_ == "none")
    {
        return tsigma;
    }

    if (!db().foundObject<volScalarField>(surfChargeName_))
    {
        // A WARNING, not a failure: the field is created by the transport
        // model, which may not have registered it yet on the first construction
        // of this patch field. It is looked up again every updateCoeffs.
        return tsigma;
    }

    tsigma.ref() =
        db().lookupObject<volScalarField>(surfChargeName_)
            .boundaryField()[patch().index()];

    return tsigma;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thinDielectricPotentialFvPatchScalarField::
thinDielectricPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    epsilonR_(1.0),
    thickness_(-1.0),
    backingPotential_(0.0),
    surfChargeName_("none")
{
    this->refValue() = Zero;
    this->refGrad()  = Zero;
    this->valueFraction() = 0.0;
}


thinDielectricPotentialFvPatchScalarField::
thinDielectricPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    epsilonR_(dict.getOrDefault<scalar>("epsilonR", 1.0)),
    thickness_(dict.getOrDefault<scalar>("thickness", -1.0)),
    backingPotential_(dict.getOrDefault<scalar>("backingPotential", 0.0)),
    surfChargeName_(dict.getOrDefault<word>("surfCharge", "surfCharge"))
{
    // `thickness` AND `backingPotential` ARE A PAIR, meaningless apart.
    //
    // A thickness with no conductor behind it is a capacitance to nothing, and
    // a backing potential with no thickness has no layer to drop across. So
    // one without the other is rejected rather than defaulted -- the two case
    // `kind`s (thinDielectricSurface, thinDielectricOnElectrode) exist for
    // exactly this reason.
    const bool haveT = dict.found("thickness");
    const bool haveV = dict.found("backingPotential");

    if (haveT != haveV)
    {
        FatalIOErrorInFunction(dict)
            << "`thickness` and `backingPotential` must be given TOGETHER or"
            << " not at all, on patch " << p.name() << "." << nl << nl
            << "    Given: " << (haveT ? "thickness" : "backingPotential")
            << " only." << nl << nl
            << "    WITH BOTH -- a barrier bonded to an electrode. The layer is"
            << " a capacitor to that electrode:" << nl
            << "        epsilon dV/dn = sigma - C (V - Vb),  C = eps0*epsilonR/d"
            << nl << nl
            << "    WITH NEITHER -- a free-standing insulating sheet, open gas"
            << " or vacuum behind it. There is no" << nl
            << "    conductor to be a capacitor to, so C = 0 and the condition"
            << " is pure Neumann:" << nl
            << "        epsilon dV/dn = sigma" << nl << nl
            << "    A thickness alone would be a capacitance to nothing; a"
            << " backing potential alone has no layer to" << nl
            << "    drop across. Neither is a thing to default." << nl
            << exit(FatalIOError);
    }

    if (haveT && thickness_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "`thickness` must be positive, not " << thickness_
            << ", on patch " << p.name() << "." << nl
            << "    Omit it entirely for the free-standing case." << nl
            << exit(FatalIOError);
    }

    if (epsilonR_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "`epsilonR` must be positive, not " << epsilonR_
            << ", on patch " << p.name() << "." << nl
            << exit(FatalIOError);
    }

    // `epsilonR` DOES NOTHING ON A FREE-STANDING LAYER, so refuse it there.
    //
    // The layer's permittivity enters ONLY through C = eps0*epsilonR/d. With no
    // conductor behind the sheet, C = 0 and the condition is
    //     eps_g dV/dn = sigma
    // in which the only permittivity is the GAS one -- read from the region,
    // not from here. Accepting `epsilonR` in that branch would let a user set
    // 4.6, read the case back, and believe the barrier material mattered.
    if (!haveT && dict.found("epsilonR"))
    {
        FatalIOErrorInFunction(dict)
            << "`epsilonR` has NO EFFECT on a free-standing layer, on patch "
            << p.name() << "." << nl << nl
            << "    Given without `thickness`, so this patch is free-standing:"
            << " open gas or vacuum behind" << nl
            << "    the sheet, no conductor to be a capacitor to. The condition"
            << " is then pure Neumann," << nl
            << "        eps_g dV/dn = sigma" << nl
            << "    whose only permittivity is the GAS permittivity, taken from"
            << " the region. The layer's" << nl
            << "    own epsilonR appears only in C = eps0*epsilonR/d, and there"
            << " is no d here." << nl << nl
            << "    REMOVE `epsilonR` -- a free-standing patch needs only"
            << " `surfCharge`." << nl
            << "    Or, if the sheet IS bonded to an electrode, add BOTH"
            << " `thickness` and `backingPotential`," << nl
            << "    which is what makes epsilonR mean something." << nl
            << exit(FatalIOError);
    }

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (!this->readMixedEntries(dict))
    {
        this->refValue() = *this;
        this->refGrad()  = Zero;
        this->valueFraction() = 0.0;
    }
}


thinDielectricPotentialFvPatchScalarField::
thinDielectricPotentialFvPatchScalarField
(
    const thinDielectricPotentialFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    epsilonR_(ptf.epsilonR_),
    thickness_(ptf.thickness_),
    backingPotential_(ptf.backingPotential_),
    surfChargeName_(ptf.surfChargeName_)
{}


thinDielectricPotentialFvPatchScalarField::
thinDielectricPotentialFvPatchScalarField
(
    const thinDielectricPotentialFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    epsilonR_(ptf.epsilonR_),
    thickness_(ptf.thickness_),
    backingPotential_(ptf.backingPotential_),
    surfChargeName_(ptf.surfChargeName_)
{}


thinDielectricPotentialFvPatchScalarField::
thinDielectricPotentialFvPatchScalarField
(
    const thinDielectricPotentialFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    epsilonR_(ptf.epsilonR_),
    thickness_(ptf.thickness_),
    backingPotential_(ptf.backingPotential_),
    surfChargeName_(ptf.surfChargeName_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thinDielectricPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // The permittivity of the region this patch lives in -- the GAS side. Read
    // the same way coupledElectricPotential reads it, so the two cannot
    // disagree: the region's own registered dielectricProperties, else the
    // electromagneticsModel's value.
    scalar epsilonRgas = 1.0;
    if (db().foundObject<IOdictionary>("dielectricProperties"))
    {
        epsilonRgas =
            db().lookupObject<IOdictionary>("dielectricProperties")
                .getOrDefault<scalar>("epsilonR", 1.0);
    }
    else if
    (
        db().time().foundObject<electromagneticsModel>
        (
            "electromagneticsProperties"
        )
    )
    {
        epsilonRgas =
            db().time().lookupObject<electromagneticsModel>
            (
                "electromagneticsProperties"
            ).epsilonR();
    }

    const scalar epsG = constant::plasma::epsilon0.value()*epsilonRgas;
    const scalarField sigma(surfaceCharge());
    const scalar C = layerCapacitance();

    if (C <= 0)
    {
        // FREE-STANDING: pure Neumann, epsG dV/dn = sigma.
        //
        // Implemented as refGradient rather than through refValue, and this is
        // not a style choice: refValue = Vb + sigma/C DIVERGES as C -> 0. The
        // product f*refValue stays finite -- it tends to sigma/(epsG*Delta) --
        // so the limit is correct but unreachable through that expression.
        valueFraction() = 0.0;
        refGrad() = sigma/epsG;
        refValue() = *this;                 // unused when f = 0
    }
    else
    {
        // BACKED: the Robin condition epsG dV/dn = sigma - C (V - Vb).
        //
        //     f        = C / (C + epsG*Delta)
        //     refValue = Vb + sigma/C
        //
        // Safe here because C > 0 is guaranteed by the branch.
        const scalarField& delta = patch().deltaCoeffs();

        valueFraction() = C/(C + epsG*delta);
        refValue() = backingPotential_ + sigma/C;
        refGrad() = Zero;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void thinDielectricPotentialFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    // ROUND-TRIP INVARIANT: write only what read() will ACCEPT.
    //
    // All three keys are backed-only -- `epsilonR` is rejected on a
    // free-standing patch because it does nothing there -- so writing
    // `epsilonR` unconditionally made a written field UNREADABLE: a restart
    // from `startFrom latestTime` hit that guard and aborted. Caught 2026-09-03
    // on needleDBD, where the generated field carried `epsilonR 1` on a
    // free-standing patch that had never been given one.
    if (thickness_ > 0)
    {
        os.writeEntry("epsilonR", epsilonR_);
        os.writeEntry("thickness", thickness_);
        os.writeEntry("backingPotential", backingPotential_);
    }

    os.writeEntry("surfCharge", surfChargeName_);

    fvPatchField<scalar>::writeValueEntry(os);
    this->refValue().writeEntry("refValue", os);
    this->refGrad().writeEntry("refGradient", os);
    this->valueFraction().writeEntry("valueFraction", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    thinDielectricPotentialFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
