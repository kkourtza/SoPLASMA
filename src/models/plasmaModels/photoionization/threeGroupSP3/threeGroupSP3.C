/*---------------------------------------------------------------------------*\
  File: threeGroupSP3.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::threeGroupSP3.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fixedValueFvPatchFields.H"
#include "mixedFvPatchFields.H"

#include "SoPLASMAConstants.H"
#include "threeGroupSP3.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(threeGroupSP3, 0);
addToRunTimeSelectionTable
(
    photoionizationModel,
    threeGroupSP3,
    dictionary
);

// * * * * * * * * * * * * * * SP3 Closure Constants  * * * * * * * * * * * //

// kappa^2_{1,2} = 3/7 -+ (2/7)*sqrt(6/5)
const scalar threeGroupSP3::kappa1Sqr_ =
    3.0/7.0 - (2.0/7.0)*Foam::sqrt(6.0/5.0);

const scalar threeGroupSP3::kappa2Sqr_ =
    3.0/7.0 + (2.0/7.0)*Foam::sqrt(6.0/5.0);

// gamma_n = (5/7)*[1 + (-1)^n * 3*sqrt(6/5)]
const scalar threeGroupSP3::gamma1_ =
    (5.0/7.0)*(1.0 - 3.0*Foam::sqrt(6.0/5.0));

const scalar threeGroupSP3::gamma2_ =
    (5.0/7.0)*(1.0 + 3.0*Foam::sqrt(6.0/5.0));

// * * * * * * * * * * * * * * Larsen BC Constants  * * * * * * * * * * * * //

// alpha_{1,2} = (5/96)*(34 +- 11*sqrt(6/5))
const scalar threeGroupSP3::alpha1_ =
    (5.0/96.0)*(34.0 + 11.0*Foam::sqrt(6.0/5.0));

const scalar threeGroupSP3::alpha2_ =
    (5.0/96.0)*(34.0 - 11.0*Foam::sqrt(6.0/5.0));

// beta_{1,2} = (5/96)*(2 -+ sqrt(6/5))
const scalar threeGroupSP3::beta1_ =
    (5.0/96.0)*(2.0 - Foam::sqrt(6.0/5.0));

const scalar threeGroupSP3::beta2_ =
    (5.0/96.0)*(2.0 + Foam::sqrt(6.0/5.0));

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

threeGroupSP3::threeGroupSP3
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
    c_
    (
        "c",
        dimensionSet(0, 1, -1, 0, 0, 0, 0),
        constant::plasma::cLight.value()
    ),
    larsenMaxIters_(1),
    larsenTol_(1e-6),
    phi1_j_(),
    phi2_j_(),
    Sph_j_()
{
    const word coeffsName(type() + "Coeffs");

    if (!found(coeffsName))
    {
        FatalIOErrorInFunction(*this)
            << "Missing required dictionary '" << coeffsName << "' in "
            << objectPath() << nl
            << exit(FatalIOError);
    }

    const dictionary& coeffs(subDict(coeffsName));

    // Boundary condition type
    boundaryConditionType_ =
        coeffs.getOrDefault<word>("boundaryConditionType", "larsen");

    if
    (
        boundaryConditionType_ != "larsen"
     && boundaryConditionType_ != "zero"
    )
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown boundaryConditionType '"
            << boundaryConditionType_ << "'." << nl
            << "Valid options: (larsen | zero)" << nl
            << exit(FatalIOError);
    }

    if (boundaryConditionType_ == "larsen")
    {
        larsenMaxIters_ =
            coeffs.getOrDefault<label>("larsenMaxIters", 1);
        larsenTol_ =
            coeffs.getOrDefault<scalar>("larsenTol", 1e-6);
    }

    // Unit conversion factors to SI
    coefficientUnits_ =
        coeffs.getOrDefault<word>("coefficientUnits", "SI");

    scalar pToSI      = 1.0;
    scalar lambdaToSI = 1.0;
    scalar AToSI      = 1.0;

    if (coefficientUnits_ == "TorrCm")
    {
        const scalar torrToPa = 101325.0/760.0;
        const scalar cmToM    = 0.01;

        pToSI      = torrToPa;
        lambdaToSI = 1.0/(cmToM*torrToPa);
        AToSI      = 1.0/(cmToM*torrToPa);
    }
    else if (coefficientUnits_ != "SI")
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown coefficientUnits '" << coefficientUnits_ << "'." << nl
            << "Valid options: (SI | TorrCm)" << nl
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
            << "It must be registered before photoionizationModel::New()."
            << nl << exit(FatalIOError);
    }

    // Fitting parameters (exactly nGroups_ pairs required)
    const List<Tuple2<scalar,scalar>> fitting
    (
        coeffs.lookup("fittingParameters")
    );

    if (fitting.size() != nGroups_)
    {
        FatalIOErrorInFunction(coeffs)
            << "Three-group SP3 requires exactly " << nGroups_
            << " (lambda A) pairs, but " << fitting.size()
            << " were given." << nl
            << exit(FatalIOError);
    }

    A_.setSize(nGroups_);
    lambda_.setSize(nGroups_);

    forAll(fitting, j)
    {
        lambda_[j] = fitting[j].first()  * lambdaToSI;
        A_[j]      = fitting[j].second() * AToSI;
    }

    // Patch types for phi1_j and phi2_j fields
    // Mixed (Robin) for Larsen; fixedValue (zero) otherwise.
    // Constraint patches (empty, symmetry, cyclic) are always preserved.
    const polyBoundaryMesh& bmesh = mesh_.boundaryMesh();

    const word nonConstraintType =
    (
        boundaryConditionType_ == "larsen"
      ? mixedFvPatchScalarField::typeName
      : fixedValueFvPatchScalarField::typeName
    );

    wordList patchTypes(bmesh.size(), nonConstraintType);

    forAll(bmesh, patchi)
    {
        if (polyPatch::constraintType(bmesh[patchi].type()))
        {
            patchTypes[patchi] = bmesh[patchi].type();
        }
    }

    // phi1_j, phi2_j and Sph_j fields
    phi1_j_.resize(nGroups_);
    phi2_j_.resize(nGroups_);

    forAll(phi1_j_, j)
    {
        phi1_j_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "phi1_" + Foam::name(j + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero",
                    dimensionSet(0, -3, 0, 0, 0, 0, 0),
                    0.0
                ),
                patchTypes
            )
        );

        phi2_j_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "phi2_" + Foam::name(j + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero",
                    dimensionSet(0, -3, 0, 0, 0, 0, 0),
                    0.0
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

void threeGroupSP3::correct()
{
    const label timeIndex = mesh_.time().timeIndex();

    if (!mesh_.changing() && (timeIndex % solveInterval_ != 0))
    {
        return;
    }

    Info<< "Solving photoionization (three-group SP3, "
        << boundaryConditionType_ << " BC)..." << endl;

    const volScalarField& source =
        mesh_.lookupObject<volScalarField>(sourceName_);

    const polyBoundaryMesh& bmesh = mesh_.boundaryMesh();

    forAll(A_, j)
    {
        // (lambda_j*p)^2/kappa1^2  [m^-2]
        const dimensionedScalar mu1Sqr
        (
            "mu1Sqr",
            dimensionSet(0, -2, 0, 0, 0, 0, 0),
            sqr(lambda_[j]*p_.value())/kappa1Sqr_
        );

        // (lambda_j*p)^2/kappa2^2  [m^-2]
        const dimensionedScalar mu2Sqr
        (
            "mu2Sqr",
            dimensionSet(0, -2, 0, 0, 0, 0, 0),
            sqr(lambda_[j]*p_.value())/kappa2Sqr_
        );

        // (lambda_j*p)/(kappa1^2*c) * C  [m^-2 s]
        // C = 1 in emissionRate mode; C = pq/(p+pq)*xi*(nu_u/nu_i) otherwise
        const dimensionedScalar kappaSrc1
        (
            "kappaSrc1",
            dimensionSet(0, -2, 1, 0, 0, 0, 0),
            (lambda_[j]*p_.value())/(kappa1Sqr_*c_.value())*C_
        );

        // (lambda_j*p)/(kappa2^2*c) * C  [m^-2 s]
        const dimensionedScalar kappaSrc2
        (
            "kappaSrc2",
            dimensionSet(0, -2, 1, 0, 0, 0, 0),
            (lambda_[j]*p_.value())/(kappa2Sqr_*c_.value())*C_
        );

        // A_j*p*c  [s^-1]
        const dimensionedScalar Ajpc
        (
            "Ajpc",
            dimensionSet(0, 0, -1, 0, 0, 0, 0),
            A_[j]*p_.value()*c_.value()
        );

        // Apply diagonal part of Larsen BC once per group per correct() call.
        // valueFraction = k/(k + deltaCoeffs), refValue = refGrad = 0.
        // Zero BC patches need no update.
        if (boundaryConditionType_ == "larsen")
        {
            const scalar k1 = lambda_[j]*p_.value()*alpha1_;
            const scalar k2 = lambda_[j]*p_.value()*alpha2_;

            forAll(bmesh, patchi)
            {
                if (!polyPatch::constraintType(bmesh[patchi].type()))
                {
                    mixedFvPatchScalarField& pf1 =
                        refCast<mixedFvPatchScalarField>
                        (
                            phi1_j_[j].boundaryFieldRef()[patchi]
                        );

                    mixedFvPatchScalarField& pf2 =
                        refCast<mixedFvPatchScalarField>
                        (
                            phi2_j_[j].boundaryFieldRef()[patchi]
                        );

                    const scalarField& deltaCoeffs = pf1.patch().deltaCoeffs();

                    pf1.refValue() = 0.0;
                    pf1.refGrad()  = 0.0;
                    pf1.valueFraction() = k1/(k1 + deltaCoeffs);

                    pf2.refValue() = 0.0;
                    pf2.refGrad()  = 0.0;
                    pf2.valueFraction() = k2/(k2 + deltaCoeffs);
                }
            }
        }

        // Picard iteration on Larsen beta cross-coupling.
        // Iteration 0: diagonal only (refGrad = 0).
        // Iterations 1+: beta cross-term added to refGrad from previous solve.
        List<scalarField> phi1BoundaryPrev(bmesh.size());
        List<scalarField> phi2BoundaryPrev(bmesh.size());

        for (label iter = 0; iter < larsenMaxIters_; ++iter)
        {
            if (boundaryConditionType_ == "larsen" && larsenMaxIters_ > 1)
            {
                Info<< "    SP3 group " << j + 1 << ", Larsen iter "
                    << iter + 1 << "/" << larsenMaxIters_ << endl;
            }

            // Update beta cross-coupling from previous iteration's solution
            if (boundaryConditionType_ == "larsen" && iter > 0)
            {
                forAll(bmesh, patchi)
                {
                    if (!polyPatch::constraintType(bmesh[patchi].type()))
                    {
                        mixedFvPatchScalarField& pf1 =
                            refCast<mixedFvPatchScalarField>
                            (
                                phi1_j_[j].boundaryFieldRef()[patchi]
                            );

                        mixedFvPatchScalarField& pf2 =
                            refCast<mixedFvPatchScalarField>
                            (
                                phi2_j_[j].boundaryFieldRef()[patchi]
                            );

                        pf1.refGrad() =
                            -lambda_[j]*p_.value()*beta2_
                           *phi2_j_[j].boundaryField()[patchi];

                        pf2.refGrad() =
                            -lambda_[j]*p_.value()*beta1_
                           *phi1_j_[j].boundaryField()[patchi];
                    }
                }
            }

            // Snapshot boundary values before solve for convergence check
            if (boundaryConditionType_ == "larsen" && larsenMaxIters_ > 1)
            {
                forAll(bmesh, patchi)
                {
                    if (!polyPatch::constraintType(bmesh[patchi].type()))
                    {
                        phi1BoundaryPrev[patchi] =
                            phi1_j_[j].boundaryField()[patchi]
                           .patchInternalField();
                        phi2BoundaryPrev[patchi] =
                            phi2_j_[j].boundaryField()[patchi]
                           .patchInternalField();
                    }
                }
            }

            fvScalarMatrix phi1jEqn
            (
                fvm::laplacian(phi1_j_[j])
              - fvm::Sp(mu1Sqr, phi1_j_[j])
             ==
              - kappaSrc1*source
            );

            phi1jEqn.solve(mesh_.solver(phi1_j_[j].name()));

            fvScalarMatrix phi2jEqn
            (
                fvm::laplacian(phi2_j_[j])
              - fvm::Sp(mu2Sqr, phi2_j_[j])
             ==
              - kappaSrc2*source
            );

            phi2jEqn.solve(mesh_.solver(phi2_j_[j].name()));

            // Check max relative change in boundary values across all patches.
            // Normalised against the global field maximum.
            if (boundaryConditionType_ == "larsen" && iter > 0)
            {
                const scalar globalRefPhi1 = 
                            gMax(mag(phi1_j_[j].primitiveField())) + SMALL;
                const scalar globalRefPhi2 = 
                            gMax(mag(phi2_j_[j].primitiveField())) + SMALL;

                scalar relDelta1 = 0.0;
                scalar relDelta2 = 0.0;

                forAll(bmesh, patchi)
                {
                    if (!polyPatch::constraintType(bmesh[patchi].type()))
                    {
                        const scalarField phi1New =
                            phi1_j_[j].boundaryField()[patchi]
                           .patchInternalField();
                        const scalarField phi2New =
                            phi2_j_[j].boundaryField()[patchi]
                           .patchInternalField();

                        relDelta1 = max
                        (
                            relDelta1,
                            max(mag(phi1New - phi1BoundaryPrev[patchi]))
                           /globalRefPhi1
                        );
                        relDelta2 = max
                        (
                            relDelta2,
                            max(mag(phi2New - phi2BoundaryPrev[patchi]))
                           /globalRefPhi2
                        );
                    }
                }

                scalar globalDelta1 = returnReduce(relDelta1, maxOp<scalar>());
                scalar globalDelta2 = returnReduce(relDelta2, maxOp<scalar>());

                if (globalDelta1 < larsenTol_ && globalDelta2 < larsenTol_)
                {
                    break;
                }
            }
        }

        // SP3 closure: Psi_j = (gamma2*phi1_j - gamma1*phi2_j)/(gamma2-gamma1)
        // Recover partial source: S_j = A_j*p*c * Psi_j
        Sph_j_[j] =
            Ajpc
           *(gamma2_*phi1_j_[j] - gamma1_*phi2_j_[j])
           /(gamma2_ - gamma1_);
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
