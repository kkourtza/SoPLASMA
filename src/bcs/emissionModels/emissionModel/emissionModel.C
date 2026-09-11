/*---------------------------------------------------------------------------*\
  File: emissionModel.C
  Part of: SoPLASMA
\*---------------------------------------------------------------------------*/

#include "emissionModel.H"
#include "volFields.H"
#include "materialLibrary.H"

// * * * * * * * * * * * * * * Runtime Information * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(emissionModel, 0);
    defineRunTimeSelectionTable(emissionModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emissionModel::emissionModel
(
    const fvPatch& p,
    const dictionary& dict
)
:
    patch_(p),
    dict_(dict)
{}


Foam::autoPtr<Foam::emissionModel> Foam::emissionModel::New
(
    const word& modelType,
    const fvPatch& p,
    const dictionary& dict
)
{
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown emission model `" << modelType << "` on patch `"
            << p.name() << "`." << nl << nl
            << "    Available: " << dictionaryConstructorTablePtr_->sortedToc()
            << nl << exit(FatalIOError);
    }

    return autoPtr<emissionModel>(ctorPtr(p, dict));
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::emissionModel::patchField() const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    if (!mesh.foundObject<volScalarField>("Emag"))
    {
        FatalErrorInFunction
            << "A field-driven emission model on patch `" << patch_.name()
            << "` needs `Emag`, which does not exist." << nl << nl
            << "    Not defaulted to zero: a field-driven mechanism with no"
               " field cannot be evaluated," << nl
            << "    and returning zero would switch it off invisibly." << nl
            << exit(FatalError);
    }

    return tmp<scalarField>::New
    (
        mesh.lookupObject<volScalarField>("Emag")
            .boundaryField()[patch_.index()]
    );
}


Foam::tmp<Foam::scalarField> Foam::emissionModel::patchTemperature() const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    if (!mesh.foundObject<volScalarField>("T_gas"))
    {
        FatalErrorInFunction
            << "A temperature-driven emission model on patch `"
            << patch_.name() << "` needs `T_gas`, which does not exist." << nl
            << nl
            << "    T_gas exists as a FIELD only when the gas energy equation"
               " is solved. Set" << nl
            << "        backgroundGas { energy { solve true; ... } }" << nl
            << "    in constant/plasmaSpeciesProperties, and set THAT field's"
               " boundary condition on" << nl
            << "    this patch to the wall temperature you want -- its boundary"
               " value IS the wall" << nl
            << "    temperature this model reads." << nl
            << exit(FatalError);
    }

    return tmp<scalarField>::New
    (
        mesh.lookupObject<volScalarField>("T_gas")
            .boundaryField()[patch_.index()]
    );
}


Foam::scalar Foam::emissionModel::workFunction() const
{
    // PRECEDENCE, as everywhere else: an explicit number WINS and is an
    // override; otherwise the named material supplies it.
    if (dict_.found("workFunction"))
    {
        return dict_.get<scalar>("workFunction");
    }

    if (dict_.found("material"))
    {
        const word mat(dict_.get<word>("material"));

        return materialLibrary::get
        (
            mat,
            "phi",
            "emission model on patch `" + patch_.name() + "`",
            "workFunctions"
        );
    }

    FatalIOErrorInFunction(dict_)
        << "This emission model on patch `" << patch_.name()
        << "` needs a work function." << nl << nl
        << "    Give it either way:" << nl
        << "        material        copper;      // from"
           " $SoPLASMA_ETC/materials/workFunctions" << nl
        << "        workFunction    4.65;        // [eV], an explicit override"
        << nl << nl
        << "    It is not defaulted. Field and thermionic emission depend"
           " EXPONENTIALLY on the work" << nl
        << "    function, so a default would be a fabricated current, not a"
           " reasonable starting point." << nl
        << exit(FatalIOError);

    return 0;
}


// ************************************************************************* //
