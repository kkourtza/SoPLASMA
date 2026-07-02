/*---------------------------------------------------------------------------*\
  File: threeGroupEddington.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::threeGroupEddington.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fixedValueFvPatchFields.H"
#include "mixedFvPatchFields.H"

#include "SoPLASMAConstants.H"
#include "threeGroupEddington.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(threeGroupEddington, 0);
addToRunTimeSelectionTable
(
    photoionizationModel,
    threeGroupEddington,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

threeGroupEddington::threeGroupEddington
(
    const fvMesh& mesh
)
:
    photoionizationModel(mesh),
    coefficientUnits_(word::null),
    boundaryConditionType_(word::null),
    p_("p", dimensionSet(1, 0, -2, 0, 0, 0, 0), 0.0),
    sourceMode_(word::null),
    sourceName_(word::null),
    quenchingFactor_(1.0),
    xiNuRatio_(1.0),
    C_(1.0),
    A_(),
    lambda_(),
    c_("c", dimensionSet(0, 1, -1, 0, 0, 0, 0), constant::plasma::cLight.value()),
    PsiStar_j_(),
    Sph_j_()
{
    const word coeffsName(type() + "Coeffs");

    if (!found(coeffsName))
    {
        FatalIOErrorInFunction(*this)
            << "Missing required dictionary '" << coeffsName << "' in "
            << objectPath() << nl << exit(FatalIOError);
    }

    const dictionary& coeffs(subDict(coeffsName));

    // Boundary condition type
    boundaryConditionType_ =
        coeffs.getOrDefault<word>("boundaryConditionType", "marshak");

    if
    (
        boundaryConditionType_ != "marshak"
     && boundaryConditionType_ != "zero"
    )
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown boundaryConditionType '" << boundaryConditionType_
            << "'." << nl
            << "Valid options are: (marshak | zero)" << nl
            << exit(FatalIOError);
    }

    // Unit conversion factors to SI
    coefficientUnits_ = coeffs.getOrDefault<word>("coefficientUnits", "SI");

    scalar pToSI = 1.0;
    scalar lambdaToSI = 1.0;
    scalar AToSI = 1.0;

    if (coefficientUnits_ == "TorrCm")
    {
        const scalar torrToPa = 101325.0/760.0;
        const scalar cmToM    = 0.01;

        pToSI      = torrToPa;
        lambdaToSI = 1.0 / (cmToM * torrToPa);
        AToSI      = 1.0 / (cmToM * torrToPa);
    }
    else if (coefficientUnits_ != "SI")
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown coefficientUnits '" << coefficientUnits_ << "'." << nl
            << "Valid options are: (SI | TorrCm)" << nl
            << exit(FatalIOError);
    }

    // Gas pressure
    p_.value() = coeffs.get<scalar>("p") * pToSI;

    // Source mode and emission prefactor C
    sourceMode_ = coeffs.get<word>("sourceMode");

    if (sourceMode_ == "emissionRate")
    {
        // source = I already includes (pq/(p+pq))*xi*(nu_u/nu_i)*Siz
        sourceName_      = coeffs.get<word>("I");
        quenchingFactor_ = 1.0;
        xiNuRatio_       = 1.0;
        C_               = 1.0;
    }
    else if (sourceMode_ == "ionizationRate")
    {
        // source = Siz; prefactor C applied internally
        sourceName_      = coeffs.get<word>("Siz");
        quenchingFactor_ = coeffs.get<scalar>("quenchingFactor");
        xiNuRatio_       = coeffs.get<scalar>("xiNuRatio");
        C_               = quenchingFactor_ * xiNuRatio_;
    }
    else
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown sourceMode '" << sourceMode_ << "'." << nl
            << "Valid options: (emissionRate | ionizationRate)" << nl
            << exit(FatalIOError);
    }

    if (!mesh_.foundObject<volScalarField>(sourceName_))
    {
        FatalIOErrorInFunction(coeffs)
            << "Source field '" << sourceName_
            << "' not found in the mesh object registry." << nl
            << exit(FatalIOError);
    }

    // Fitting parameters (exactly nGroups_ pairs required)
    const List<Tuple2<scalar, scalar>> fitting
    (
        coeffs.lookup("fittingParameters")
    );

    if (fitting.size() != nGroups_)
    {
        FatalIOErrorInFunction(coeffs)
            << "The three-group Eddington model requires exactly "
            << nGroups_ << " (lambda  A) pairs in 'fittingParameters', "
            << "but " << fitting.size() << " were given." << nl
            << exit(FatalIOError);
    }

    A_.setSize(nGroups_);
    lambda_.setSize(nGroups_);

    forAll(fitting, j)
    {
        lambda_[j] = fitting[j].first()  * lambdaToSI;
        A_[j]      = fitting[j].second() * AToSI;
    }

    // Patch types for Psi*_j fields
    // Mixed (Robin) for Marshak; fixedValue (zero) otherwise.
    // Constraint patches (empty, symmetry, cyclic) are always preserved.
    const polyBoundaryMesh& bmesh = mesh_.boundaryMesh();

    const word nonConstraintPatchType =
    (
        boundaryConditionType_ == "marshak"
    ? mixedFvPatchScalarField::typeName
    : fixedValueFvPatchScalarField::typeName
    );

    wordList patchTypes(bmesh.size(), nonConstraintPatchType);

    forAll(bmesh, patchi)
    {
        if (polyPatch::constraintType(bmesh[patchi].type()))
        {
            patchTypes[patchi] = bmesh[patchi].type();
        }
    }

    // Psi*_j and Sph_j fields
    PsiStar_j_.resize(nGroups_);

    forAll(PsiStar_j_, j)
    {
        PsiStar_j_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "PsiStar_" + Foam::name(j + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero", dimensionSet(0, -3, 0, 0, 0, 0, 0), 0.0
                ),
                patchTypes
            )
        );
    }

    Sph_j_.resize(nGroups_);

    forAll(Sph_j_, j)
    {
        Sph_j_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "Sph_" + Foam::name(j + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", Sph_.dimensions(), 0.0)
            )
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void threeGroupEddington::correct()
{
    const label timeIndex = mesh_.time().timeIndex();

    if (!mesh_.changing() && (timeIndex % solveInterval_ != 0))
    {
        return;
    }

    Info<< "Solving photoionization (three-group Eddington, "
        << boundaryConditionType_ << " BC)..." << endl;

    const volScalarField& source =
        mesh_.lookupObject<volScalarField>(sourceName_);

    const polyBoundaryMesh& bmesh = mesh_.boundaryMesh();

    // Solve each Eddington group equation independently
    forAll(A_, j)
    {
        // 3*(lambda_j*p)^2  [m^-2]
        const dimensionedScalar mu2
        (
            "mu2",
            dimensionSet(0, -2, 0, 0, 0, 0, 0),
            3.0*sqr(lambda_[j] * p_.value())
        );

        // 3*lambda_j*p/c * C  [m^-2 s]
        // C = 1 in emissionRate mode; C = pq/(p+pq)*xi*(nu_u/nu_i) otherwise
        const dimensionedScalar kappa
        (
            "kappa",
            dimensionSet(0, -2, 1, 0, 0, 0, 0),
            3.0*lambda_[j]*p_.value()/c_.value() * C_
        );

        // A_j*p*c  [s^-1]
        const dimensionedScalar Ajpc
        (
            "Ajpc",
            dimensionSet(0, 0, -1, 0, 0, 0, 0),
            A_[j]*p_.value()*c_.value()
        );

        // Apply Marshak BC: grad(Psi*_j) & n = -(3/2)*lambda_j*p * Psi*_j
        // Implemented as OpenFOAM mixed BC with refValue=0, refGrad=0,
        // valueFraction = k/(k + deltaCoeffs), k = (3/2)*lambda_j*p.
        // Updated every correct() call; zero BC patches need no update.
        if (boundaryConditionType_ == "marshak")
        {
            forAll(bmesh, patchi)
            {
                if (!polyPatch::constraintType(bmesh[patchi].type()))
                {
                    mixedFvPatchScalarField& pf =
                        refCast<mixedFvPatchScalarField>
                        (
                            PsiStar_j_[j].boundaryFieldRef()[patchi]
                        );

                    const scalar k = 1.5*lambda_[j]*p_.value();

                    const scalarField& deltaCoeffs =
                        pf.patch().deltaCoeffs();

                    pf.refValue()      = scalarField(pf.size(), 0.0);
                    pf.refGrad()       = scalarField(pf.size(), 0.0);
                    pf.valueFraction() = k/(k + deltaCoeffs);
                }
            }
        }

        fvScalarMatrix PsiStarjEqn
        (
            fvm::laplacian(PsiStar_j_[j])
          - fvm::Sp(mu2, PsiStar_j_[j])
         ==
          - kappa * source
        );

        PsiStarjEqn.solve(mesh_.solver(PsiStar_j_[j].name()));

        // Recover partial photoionization source: S_j = A_j*p*c * Psi*_j
        Sph_j_[j] = Ajpc*PsiStar_j_[j];
    }

    // Accumulate partial sources into the total photoionization source
    Sph_ == dimensionedScalar(Sph_.dimensions(), Zero);

    forAll(Sph_j_, j)
    {
        Sph_ += Sph_j_[j];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
