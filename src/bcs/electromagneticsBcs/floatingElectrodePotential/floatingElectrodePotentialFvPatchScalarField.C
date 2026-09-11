/*---------------------------------------------------------------------------*\
  File: floatingElectrodePotentialFvPatchScalarField.C
  Part of: SoPLASMA
\*---------------------------------------------------------------------------*/

#include "floatingElectrodePotentialFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floatingElectrodePotentialFvPatchScalarField::
floatingElectrodePotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    initialCharge_(0),
    floatingPotential_(0)
{}


Foam::floatingElectrodePotentialFvPatchScalarField::
floatingElectrodePotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    initialCharge_(dict.getOrDefault<scalar>("initialCharge", 0)),
    floatingPotential_(dict.getOrDefault<scalar>("floatingPotential", 0))
{
    // A POTENTIAL IS NOT AN INPUT HERE.
    //
    // The whole point of this condition is that V is the unknown. Accepting a
    // potential would let a case ask for a contradiction -- a conductor whose
    // potential is both prescribed and solved for -- and the run would look
    // fine while quietly ignoring one of the two.
    //
    // `value` is exempt: OpenFOAM requires it as the field's stored patch
    // value, and it is overwritten on the first solve. `floatingPotential` is
    // exempt because it is what THIS class writes out, so a restart must be
    // able to read its own output back (the round-trip invariant).
    for (const word& banned : {"potential", "uniformValue", "V", "voltage"})
    {
        if (dict.found(banned))
        {
            FatalIOErrorInFunction(dict)
                << "`" << banned << "` is not an input to"
                << " floatingElectrodePotential, on patch " << p.name() << "."
                << nl << nl
                << "    A floating electrode is defined by its CHARGE, not its"
                   " potential: the potential" << nl
                << "    is the UNKNOWN this condition solves for. Prescribing"
                   " one would make this a" << nl
                << "    driven electrode -- use `uniformFixedValue` for that."
                << nl << nl
                << "    The only genuinely free parameter is `initialCharge`"
                   " [C], and 0 is the honest" << nl
                << "    default: an electrode that has never been connected to"
                   " anything starts neutral." << nl
                << exit(FatalIOError);
        }
    }

    // The stored value must be the equipotential from the outset, so the very
    // first Poisson solve sees a consistent conductor rather than whatever
    // non-uniform field happened to be in the file.
    fvPatchScalarField::operator==(floatingPotential_);
}


Foam::floatingElectrodePotentialFvPatchScalarField::
floatingElectrodePotentialFvPatchScalarField
(
    const floatingElectrodePotentialFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    initialCharge_(ptf.initialCharge_),
    floatingPotential_(ptf.floatingPotential_)
{}


Foam::floatingElectrodePotentialFvPatchScalarField::
floatingElectrodePotentialFvPatchScalarField
(
    const floatingElectrodePotentialFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    initialCharge_(ptf.initialCharge_),
    floatingPotential_(ptf.floatingPotential_)
{}


Foam::floatingElectrodePotentialFvPatchScalarField::
floatingElectrodePotentialFvPatchScalarField
(
    const floatingElectrodePotentialFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    initialCharge_(ptf.initialCharge_),
    floatingPotential_(ptf.floatingPotential_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::floatingElectrodePotentialFvPatchScalarField::setFloatingPotential
(
    const scalar Vf
)
{
    floatingPotential_ = Vf;

    // UNIFORM by construction. The equipotential condition is not "nearly
    // constant over the patch" -- it is exactly one scalar, and assigning it
    // face-by-face from a field would let discretisation error accumulate into
    // a fictitious surface potential gradient on a metal.
    fvPatchScalarField::operator==(Vf);
}


void Foam::floatingElectrodePotentialFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    os.writeEntry("initialCharge", initialCharge_);

    // Written so a restart resumes at the potential it reached, and so the
    // settling behaviour can be READ. The floating potential is the single
    // most informative number about such an electrode.
    os.writeEntry("floatingPotential", floatingPotential_);

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        floatingElectrodePotentialFvPatchScalarField
    );
}

// ************************************************************************* //
