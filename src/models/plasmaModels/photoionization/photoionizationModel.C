/*---------------------------------------------------------------------------*\
  File: photoionizationModel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::photoionizationModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "photoionizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(photoionizationModel, 0);
defineRunTimeSelectionTable(photoionizationModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from mesh
photoionizationModel::photoionizationModel(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "photoionizationProperties",
            mesh.time().constant(),
            mesh,
            // READ_IF_PRESENT, not MUST_READ. A case that does not model
            // photoionization should not need a file to say so -- see New().
            // A model that DOES need coefficients still fails clearly on its
            // own missing `<type>Coeffs` sub-dictionary.
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    Sph_
    (
        IOobject
        (
            "Sph",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    solveInterval_(getOrDefault<label>("solveInterval", 1))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<photoionizationModel> photoionizationModel::New
(
  const fvMesh& mesh
)
{
    IOdictionary tmpDict
    (
        IOobject
        (
            "photoionizationProperties",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    // DEFAULTS TO `none`, and the file is therefore OPTIONAL.
    //
    // This was MUST_READ with a bare get<word>(), so every case had to carry a
    // constant/<region>/photoionizationProperties whose entire content was
    // `photoionizationModel none;` -- and in a multi-region case, one per
    // region. Omitting it stopped the run with
    //     cannot find file ".../constant/gas/photoionizationProperties"
    // which says nothing about photoionization being the thing at issue.
    //
    // A file required only in order to say "off" is a file that should not be
    // required: `none` is the honest default, because photoionization is an
    // ADDITIONAL source and leaving it out changes nothing else. Contrast
    // `electronEnergyModel`, which is deliberately required with no default --
    // there, LFA and LMEA are two different physics and neither can be assumed.
    // Here there is one sensible default and it is the absence of a model.
    //
    // Said out loud rather than silently, because a photoionization model that
    // is off without the case saying so is worth one line at start-up.
    const word modelName
    (
        tmpDict.getOrDefault<word>("photoionizationModel", "none")
    );

    if (!tmpDict.found("photoionizationModel"))
    {
        Info<< "photoionizationModel: none (default -- no"
            << " constant/" << mesh.name() << "/photoionizationProperties)."
            << nl
            << "    Photoionization is an ADDITIONAL electron source, so"
            << " omitting it changes nothing else. State" << nl
            << "    `photoionizationModel nTermHelmholtz;` (or"
            << " threeGroupEddington | threeGroupSP3) to enable it." << endl;
    }

    // Look up the constructor in the table
    auto* ctorPtr = dictionaryConstructorTable(modelName);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(tmpDict)
            << "Unknown photoionizationModel type '" << modelName << "'\n"
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return autoPtr<photoionizationModel>(ctorPtr(mesh));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
