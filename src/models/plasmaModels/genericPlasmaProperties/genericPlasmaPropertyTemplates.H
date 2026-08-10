/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::ConstantProperty
    Foam::TabulatedProperty1D
    Foam::TabulatedProperty2D
    Foam::PowerLawProperty
    Foam::Function1Property

Description:
    Generic templated implementations for plasma physical properties.

    This header provides a library of mathematical logic that is independent of
    the specific physical quantity being modeled (e.g., Mobility, Diffusivity,
    Thermal Conductivity etc.).

    By using template specialization, these models are "stamped" onto physical 
    base classes (like plasmaMobilityModel). This ensures code reuse, 
    mathematical consistency across different physics, and centralized 
    maintenance.

\*---------------------------------------------------------------------------*/

#ifndef genericPlasmaPropertyTemplates_H
#define genericPlasmaPropertyTemplates_H

#include "fvMesh.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "Function1.H"
#include "interpolationTable.H"
#include "interpolation2DTable.H"
#include "OSspecific.H"

namespace Foam
{

// Used by all non-constant templates

namespace detail
{

template<class Mesh>
tmp<volScalarField> getScalarField
(
    const Mesh& mesh,
    const word& name,
    const word& propertyName
)
{
    if (mesh.template foundObject<volScalarField>(name))
    {
        return mesh.template lookupObject<volScalarField>(name);
    }

    if (mesh.template foundObject<volVectorField>(name))
    {
        return mag(mesh.template lookupObject<volVectorField>(name));
    }

    const wordList availableScalars =
        mesh.template names<volScalarField>();
    const wordList availableVectors =
        mesh.template names<volVectorField>();

    FatalErrorInFunction
        << "plasmaPropertyEvaluator '" << propertyName
        << "': lookupVariable '" << name
        << "' not found in mesh registry." << nl
        << "Available volScalarFields: " << availableScalars << nl
        << "Available volVectorFields: " << availableVectors << nl
        << exit(FatalError);

    return tmp<volScalarField>(nullptr);
}

} // End namespace detail

/*---------------------------------------------------------------------------*\
                      Class ConstantProperty Declaration
\*---------------------------------------------------------------------------*/

template<class Base>
class ConstantProperty
:
    public Base
{
    scalar val_;

public:
    ConstantProperty
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const word& propertyName,
        const dimensionSet& dims
    )
    :
        Base(modelType, dict, mesh, propertyName, dims),
        val_(0.0)
    {
        if (!dict.readIfPresent("value", val_))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': constant model requires entry 'value'." << nl
                << exit(FatalIOError);
        }
    }

    //- Destructor
    virtual ~ConstantProperty() = default;

    //- Member Functions

    bool isUniform() const override { return true; }
    bool isConstant() const override { return true; }

    dimensionedScalar evaluateScalar() const override
    {
        return dimensionedScalar(this->propertyName(), this->dims_, val_);
    }

    tmp<volScalarField> evaluate() const override
    {
        return tmp<volScalarField>::New
        (
            IOobject
            (
                this->propertyName(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar(this->propertyName(), this->dims_, val_)
        );
    }

    void correct(volScalarField& field) const override
    {
        field.primitiveFieldRef() = val_;

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                field.boundaryFieldRef()[patchi] = val_;
            }
        }

        field.correctBoundaryConditions();
    }

    void evaluatePatch(scalarField& pField, const label) const override
    {
        pField = val_;
    }
};

/*---------------------------------------------------------------------------*\
                      Class TabulatedProperty1D Declaration
\*---------------------------------------------------------------------------*/

template<class Base>
class TabulatedProperty1D
:
    public Base
{
    //- The interpolation table (handles file reading and interpolation)
    interpolationTable<scalar> table_;

    //- The name of the field used as the lookup coordinate (e.g., "E" or "T_e")
    word lookupVariable_;

public:

    //- Constructor
    TabulatedProperty1D
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const word& propertyName,
        const dimensionSet& dims
    )
    :
        Base(modelType, dict, mesh, propertyName, dims),
        table_(dict),
        lookupVariable_(dict.get<word>("lookupVariable"))
    {
        if (!dict.found("file"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': tabulated1D requires entry 'file'." << nl
                << exit(FatalIOError);
        }

        const fileName tableFile =
            dict.get<fileName>("file").expand();

        if (!isFile(tableFile))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': table file not found:" << nl
                << "    " << tableFile << nl
                << exit(FatalIOError);
        }
    }

    //- Destructor
    virtual ~TabulatedProperty1D() = default;

    //- Member Functions

    bool isUniform() const override { return false; }
    bool isConstant() const override { return false; }

    dimensionedScalar evaluateScalar() const override
    {
        FatalErrorInFunction
            << "evaluateScalar() called on tabulated1D evaluator '"
            << this->propertyName() << "' which is not uniform." << nl
            << exit(FatalError);
        return dimensionedScalar(this->propertyName(), this->dims_, 0.0);
    }

    tmp<volScalarField> evaluate() const override
    {
        auto tField = tmp<volScalarField>::New
        (
            IOobject
            (
                this->propertyName(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar(this->propertyName(), this->dims_, 0.0)
        );
        correct(tField.ref());
        return tField;
    }

    void correct(volScalarField& field) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const volScalarField& var = tVar();

        field.primitiveFieldRef() =
            table_.interpolateValues(var.primitiveField());

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                scalarField& pField = field.boundaryFieldRef()[patchi];
                const scalarField& pVar = var.boundaryField()[patchi];   // reuse, no re-fetch
                forAll(pField, i)
                {
                    pField[i] = table_.interpolateValue(pVar[i]);
                }
            }
        }

        field.correctBoundaryConditions();
    }

    void evaluatePatch
    (
        scalarField& pField,
        const label patchi
    ) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const scalarField& pVar = tVar().boundaryField()[patchi];

        forAll(pField, i)
        {
            pField[i] = table_.interpolateValue(pVar[i]);
        }
    }
};

/*---------------------------------------------------------------------------*\
                      Class TabulatedProperty2D Declaration
\*---------------------------------------------------------------------------*/

template<class Base>
class TabulatedProperty2D
:
    public Base
{
    //- The 2D interpolation table
    interpolation2DTable<scalar> table_;

    //- The names of the two fields used for lookup (e.g., {"E", "T_e"})
    wordList lookupVariables_;

public:

    //- Constructor
    TabulatedProperty2D
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const word& propertyName,
        const dimensionSet& dims
    )
    :
        Base(modelType, dict, mesh, propertyName, dims),
        table_(dict),
        lookupVariables_(dict.get<wordList>("lookupVariables"))
    {
        if (lookupVariables_.size() != 2)
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': tabulated2D requires exactly 2 lookupVariables." << nl
                << exit(FatalIOError);
        }

        if (!dict.found("file"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': tabulated2D requires entry 'file'." << nl
                << exit(FatalIOError);
        }

        const fileName tableFile =
            dict.get<fileName>("file").expand();

        if (!isFile(tableFile))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': table file not found:" << nl
                << "    " << tableFile << nl
                << exit(FatalIOError);
        }
    }

    //- Destructor
    virtual ~TabulatedProperty2D() = default;

    //- Member Functions

    bool isUniform() const override { return false; }
    bool isConstant() const override { return false; }

    dimensionedScalar evaluateScalar() const override
    {
        FatalErrorInFunction
            << "evaluateScalar() called on tabulated2D evaluator '"
            << this->propertyName() << "' which is not uniform." << nl
            << exit(FatalError);
        return dimensionedScalar(this->propertyName(), this->dims_, 0.0);
    }

    tmp<volScalarField> evaluate() const override
    {
        auto tField = tmp<volScalarField>::New
        (
            IOobject
            (
                this->propertyName(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar(this->propertyName(), this->dims_, 0.0)
        );
        correct(tField.ref());
        return tField;
    }

    void correct(volScalarField& field) const override
    {
        tmp<volScalarField> tX = detail::getScalarField
        (
            this->mesh_, lookupVariables_[0], this->propertyName()
        );
        tmp<volScalarField> tY = detail::getScalarField
        (
            this->mesh_, lookupVariables_[1], this->propertyName()
        );

        const volScalarField& X = tX();
        const volScalarField& Y = tY();

        forAll(field.primitiveField(), cellI)
        {
            field.primitiveFieldRef()[cellI] = table_(X[cellI], Y[cellI]);
        }

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                scalarField& pField = field.boundaryFieldRef()[patchi];
                const scalarField& pX = X.boundaryField()[patchi];       // reuse
                const scalarField& pY = Y.boundaryField()[patchi];       // reuse
                forAll(pField, i)
                {
                    pField[i] = table_(pX[i], pY[i]);
                }
            }
        }

        field.correctBoundaryConditions();
    }

    void evaluatePatch
    (
        scalarField& pField,
        const label patchi
    ) const override
    {
        tmp<volScalarField> tX = detail::getScalarField
        (
            this->mesh_, lookupVariables_[0], this->propertyName()
        );
        tmp<volScalarField> tY = detail::getScalarField
        (
            this->mesh_, lookupVariables_[1], this->propertyName()
        );

        const scalarField& pX = tX().boundaryField()[patchi];
        const scalarField& pY = tY().boundaryField()[patchi];

        forAll(pField, i)
        {
            pField[i] = table_(pX[i], pY[i]);
        }
    }
};

/*---------------------------------------------------------------------------*\
                      Class PowerLawProperty Declaration
\*---------------------------------------------------------------------------*/

template<class Base>
class PowerLawProperty
:
    public Base
{
    //- Amplitude (A)
    scalar amplitude_;

    //- Exponent (b)
    scalar exponent_;

    //- Lookup field name (e.g., "E")
    word lookupVariable_;

public:

    //- Constructor
    PowerLawProperty
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const word& propertyName,
        const dimensionSet& dims
    )
    :
        Base(modelType, dict, mesh, propertyName, dims),
        amplitude_(0.0),
        exponent_(0.0),
        lookupVariable_("none")
    {
        if (!dict.found("amplitude"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': powerLaw requires entry 'amplitude'." << nl
                << exit(FatalIOError);
        }
        amplitude_ = dict.get<scalar>("amplitude");

        if (!dict.found("exponent"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': powerLaw requires entry 'exponent'." << nl
                << exit(FatalIOError);
        }
        exponent_ = dict.get<scalar>("exponent");

        if (!dict.found("lookupVariable"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': powerLaw requires entry 'lookupVariable'." << nl
                << "Provide the name of the field to use as the independent "
                << "variable (e.g. 'Emag', 'reducedE', 'Te')." << nl
                << exit(FatalIOError);
        }
        lookupVariable_ = dict.get<word>("lookupVariable");
    }

    //- Destructor
    virtual ~PowerLawProperty() = default;

    //- Member Functions

    bool isUniform() const override { return false; }
    bool isConstant() const override { return false; }

    dimensionedScalar evaluateScalar() const override
    {
        FatalErrorInFunction
            << "evaluateScalar() called on powerLaw evaluator '"
            << this->propertyName() << "' which is not uniform." << nl
            << exit(FatalError);
        return dimensionedScalar(this->propertyName(), this->dims_, 0.0);
    }

    tmp<volScalarField> evaluate() const override
    {
        auto tField = tmp<volScalarField>::New
        (
            IOobject
            (
                this->propertyName(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar(this->propertyName(), this->dims_, 0.0)
        );
        correct(tField.ref());
        return tField;
    }

    void correct(volScalarField& field) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const volScalarField& var = tVar();

        field.primitiveFieldRef() =
            amplitude_ * pow(max(var.primitiveField(), SMALL), exponent_);

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                field.boundaryFieldRef()[patchi] =
                    amplitude_ * pow(max(var.boundaryField()[patchi], SMALL), exponent_);
            }
        }

        field.correctBoundaryConditions();
    }

    void evaluatePatch
    (
        scalarField& pField,
        const label patchi
    ) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const scalarField& pVar = tVar().boundaryField()[patchi];

        pField = amplitude_ * pow(max(pVar, SMALL), exponent_);
    }
};

/*---------------------------------------------------------------------------*\
                      Class Function1Property Declaration
\*---------------------------------------------------------------------------*/

template<class Base>
class Function1Property
:
    public Base
{
    autoPtr<Function1<scalar>> functionPtr_;

    word lookupVariable_;

public:

    //- Constructor
    Function1Property
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const word& propertyName,
        const dimensionSet& dims
    )
    :
        Base(modelType, dict, mesh, propertyName, dims),
        functionPtr_(nullptr),
        lookupVariable_("none")
    {
        if (!dict.found("lookupVariable"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': function1 requires entry 'lookupVariable'." << nl
                << "Provide the name of the field to use as the independent "
                << "variable (e.g. 'Emag', 'reducedE', 'Te')." << nl
                << exit(FatalIOError);
        }
        lookupVariable_ = dict.get<word>("lookupVariable");

        if (!dict.found("function"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': function1 requires entry 'function'." << nl
                << "Provide a valid OpenFOAM Function1<scalar> entry." << nl
                << exit(FatalIOError);
        }
        functionPtr_ = Function1<scalar>::New("function", dict, &mesh);
    }

    //- Destructor
    virtual ~Function1Property() = default;

    //- Member Functions

    bool isUniform() const override { return false; }
    bool isConstant() const override { return false; }

    dimensionedScalar evaluateScalar() const override
    {
        FatalErrorInFunction
            << "evaluateScalar() called on function1 evaluator '"
            << this->propertyName() << "' which is not uniform." << nl
            << exit(FatalError);
        return dimensionedScalar(this->propertyName(), this->dims_, 0.0);
    }

    tmp<volScalarField> evaluate() const override
    {
        auto tField = tmp<volScalarField>::New
        (
            IOobject
            (
                this->propertyName(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar(this->propertyName(), this->dims_, 0.0)
        );
        correct(tField.ref());
        return tField;
    }

    void correct(volScalarField& field) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const volScalarField& var = tVar();

        field.primitiveFieldRef() =
            functionPtr_->value(var.primitiveField());

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                field.boundaryFieldRef()[patchi] =
                    functionPtr_->value(var.boundaryField()[patchi]);
            }
        }

        field.correctBoundaryConditions();
    }

    void evaluatePatch
    (
        scalarField& pField,
        const label patchi
    ) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const scalarField& pVar = tVar().boundaryField()[patchi];

        pField = functionPtr_->value(pVar);
    }
};

/*---------------------------------------------------------------------------*\
                  Class MechanismProperty Declaration
\*---------------------------------------------------------------------------*/

//- Electron transport straight from the Boltzmann solver's own tables.
//
//  The EEDF solve already produces the electron mobility and diffusivity that
//  are consistent with the cross sections the reaction rates came from. Using
//  them here means a case cannot describe its electrons with one set of physics
//  and its chemistry with another -- which is what happens with a hand-fitted
//  `powerLaw`, and is easy to miss because both look reasonable in isolation.
//
//  This is OPTIONAL on purpose. A user may legitimately want a published fit,
//  an experimental correlation or a coded expression, for comparison against a
//  reference calculation or because the mechanism's mixture is not the one being
//  simulated. Any other evaluator type remains available per species and per
//  property; nothing forces the mechanism's values on a case.
//
//  IMPLICATION OF CHOOSING THIS. The tables are REDUCED quantities -- mu*N and
//  D*N -- because that is what is similarity-invariant. They are divided by the
//  local gas density here, so the result tracks density automatically. But they
//  were computed at ONE mixture and ONE gas temperature, and the EEDF is a
//  functional of the whole mixture: as excited-state fractions build up during a
//  discharge, the true mobility drifts away from the tabulated one. Under the
//  local-field approximation this is the same assumption already made by every
//  rate coefficient in the run, so it adds no new error -- but it is not exact,
//  and it is why `eedfRefresh` exists.
//
//  Dictionary:
//      mobility
//      {
//          type            fromMechanism;
//          quantity        muN;         // muN | DLN | DTN
//          tableDir        "constant/plasmaTables";
//          lookupVariable  reducedE;    // must match the key the tables use
//      }
template<class Base>
class MechanismProperty
:
    public Base
{
    //- Reduced quantity from the sweep, e.g. mu*N [1/(V m s)]
    interpolationTable<scalar> table_;

    word lookupVariable_;

    //- Name of the registered field holding the gas number density, used to
    //  turn the reduced quantity into the physical one
    word densityField_;

    //- Fallback when no density field is registered: the uniform background
    //  density, which is the usual case for a fixed-composition run
    scalar uniformN_;

public:

    MechanismProperty
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const word& propertyName,
        const dimensionSet& dims
    )
    :
        Base(modelType, dict, mesh, propertyName, dims),
        table_(),
        lookupVariable_(dict.getOrDefault<word>("lookupVariable", "reducedE")),
        densityField_(dict.getOrDefault<word>("densityField", "N_gas")),
        uniformN_(dict.getOrDefault<scalar>("gasNumberDensity", 0.0))
    {
        const word quantity = dict.getOrDefault<word>("quantity", "muN");
        const fileName dir =
            dict.getOrDefault<fileName>("tableDir", "constant/plasmaTables");

        // The table is produced by the in-process Boltzmann sweep, which runs
        // when plasmaReactionRates is constructed. If it is missing, the cause
        // is almost always that the case has no `chemistry` dictionary, so say
        // that rather than just naming a missing file.
        const fileName path = dir/(quantity + "_vs_" + lookupVariable_);
        if (!isFile(path))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': fromMechanism needs " << path << ", which does not"
                << " exist." << nl
                << "    These tables are written by the Boltzmann sweep, which"
                << " runs only when the case has a `chemistry` dictionary in"
                << " plasmaTransportProperties." << nl
                << "    Add one, or choose an evaluator that does not depend on"
                << " the mechanism (powerLaw, tabulated1D, function1)." << nl
                << exit(FatalIOError);
        }

        dictionary td;
        td.add("file", path);
        td.add("outOfBounds", dict.getOrDefault<word>("outOfBounds", "clamp"));
        table_ = interpolationTable<scalar>(td);
    }

    virtual ~MechanismProperty() = default;

    bool isUniform() const override { return false; }
    bool isConstant() const override { return false; }

    dimensionedScalar evaluateScalar() const override
    {
        FatalErrorInFunction
            << "evaluateScalar() called on fromMechanism evaluator '"
            << this->propertyName() << "' which is not uniform." << nl
            << exit(FatalError);
        return dimensionedScalar(this->dims_, 0.0);
    }

    //- Local gas number density, for reducing mu*N to mu
    scalar densityAt(const label celli) const
    {
        if (this->mesh_.template foundObject<volScalarField>(densityField_))
        {
            return this->mesh_.template lookupObject<volScalarField>
                (densityField_).primitiveField()[celli];
        }
        return uniformN_;
    }

    tmp<volScalarField> evaluate() const override
    {
        auto tField = tmp<volScalarField>::New
        (
            IOobject
            (
                this->propertyName(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar(this->propertyName(), this->dims_, 0.0)
        );
        correct(tField.ref());
        return tField;
    }

    void correct(volScalarField& field) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const volScalarField& var = tVar();

        if (uniformN_ <= 0
         && !this->mesh_.template foundObject<volScalarField>(densityField_))
        {
            FatalErrorInFunction
                << "fromMechanism evaluator '" << this->propertyName()
                << "' has no gas density: field '" << densityField_
                << "' is not registered and `gasNumberDensity` was not given."
                << nl
                << "    The tables hold REDUCED quantities (mu*N, D*N), so a"
                << " density is required to recover the physical value."
                << nl << exit(FatalError);
        }

        scalarField& f = field.primitiveFieldRef();
        const scalarField& x = var.primitiveField();
        forAll(f, celli)
        {
            f[celli] = table_.interpolateValue(x[celli])/densityAt(celli);
        }

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                scalarField& pF = field.boundaryFieldRef()[patchi];
                const scalarField& pX = var.boundaryField()[patchi];
                const scalar N = (uniformN_ > 0) ? uniformN_ : densityAt(0);
                forAll(pF, i)
                {
                    pF[i] = table_.interpolateValue(pX[i])/N;
                }
            }
        }

        field.correctBoundaryConditions();
    }

    void evaluatePatch(scalarField& pField, const label patchi) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const scalarField& pX = tVar().boundaryField()[patchi];
        const scalar N = (uniformN_ > 0) ? uniformN_ : densityAt(0);
        forAll(pField, i)
        {
            pField[i] = table_.interpolateValue(pX[i])/N;
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
