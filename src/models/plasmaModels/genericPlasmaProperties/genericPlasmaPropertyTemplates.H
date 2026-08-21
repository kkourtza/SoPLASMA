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
#include <fstream>
#include <string>

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

//- Bounded 1-D table lookup: interpolate inside the range, HOLD below it,
//  POWER-LAW extrapolate above it.
//
//  Written because OpenFOAM's interpolationTable does NOT extrapolate at all.
//  Its CLAMP and WARN both `return Zero` beyond the range and REPEAT wraps
//  modulo the span -- so a loss coefficient silently becomes ZERO past the
//  table, removing the only sink from the electron energy equation. That was
//  measured as a runaway: mean energy 5.95 -> 9.2 -> 26 -> 60 -> clamp, with
//  J/L reaching 1.4e4.
//
//  THE METHOD, chosen from the data rather than assumed. Fitting log(y)
//  against log(x) over the top decade of the air tables gives R^2 = 0.999 for
//  PelasticN, 0.999 for DLN, 0.996 for muEpsN -- these are power laws. The
//  rate coefficients fit far worse (nuEffN 0.90, k_ION 0.95) because they are
//  THRESHOLD processes that saturate once eps >> eps_threshold.
//
//  Using the LOCAL slope over the last few points handles both automatically:
//
//      table              top-decade p     last-5 p
//      muN                    -0.380        -0.364     smooth, unchanged
//      PelasticN               1.347         1.232     smooth, unchanged
//      PinelasticN             2.010         1.183     saturating
//      alphaN                  3.069         0.867     saturating
//      nuEffN                  4.452         1.621     saturating
//
//  So the fit self-limits exactly where the physics does, with no per-variable
//  special-casing. Extrapolating alphaN at its top-decade slope of 3.07 rather
//  than its local 0.87 would over-ionise by orders of magnitude.
//
//  Anchored at the last table point, y = y_last (x/x_last)^p, so the result is
//  CONTINUOUS at the boundary.
//
//  BELOW the range the value is HELD at the first entry, not extrapolated: the
//  tables start at zero field where the physics is already the thermal limit,
//  and a power law through the origin is meaningless there.
class boundedTable1D
{
    scalarField x_;
    scalarField y_;
    scalar slopeHi_ = 0;
    bool logSafe_ = false;
    mutable bool warned_ = false;
    word name_;

public:

    void init(const List<Tuple2<scalar, scalar>>& tbl, const word& nm)
    {
        name_ = nm;
        const label n = tbl.size();
        x_.setSize(n); y_.setSize(n);
        for (label i = 0; i < n; ++i)
        {
            x_[i] = tbl[i].first();
            y_[i] = tbl[i].second();
        }

        // Local log-log slope over the last few points. Needs strictly
        // positive x and y; anything else falls back to holding the last
        // value, which is bounded and cannot invert a sign.
        const label nFit = min(label(5), n);
        logSafe_ = (n >= 2);
        for (label i = n - nFit; i < n; ++i)
        {
            if (i < 0 || x_[i] <= 0 || y_[i] <= 0) { logSafe_ = false; break; }
        }

        if (logSafe_)
        {
            scalar sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (label i = n - nFit; i < n; ++i)
            {
                const scalar lx = log(x_[i]), ly = log(y_[i]);
                sx += lx; sy += ly; sxx += lx*lx; sxy += lx*ly;
            }
            const scalar den = nFit*sxx - sx*sx;
            slopeHi_ = (mag(den) > SMALL) ? (nFit*sxy - sx*sy)/den : 0.0;

            // BOUND the slope. Nothing here is mechanism-specific -- the fit
            // is taken from each table's own data -- but a mechanism we have
            // not seen could have a process whose threshold sits just below
            // the sweep maximum, leaving the last few points on a near
            // vertical rise. Extrapolating THAT slope would explode.
            //
            // For air every local slope is |p| < 1.7, so this never binds on
            // known data; it exists so an unfamiliar mechanism degrades to a
            // steep-but-finite power law instead of a divergence.
            const scalar pMax = 6.0;

            if (mag(slopeHi_) > pMax)
            {
                WarningInFunction
                    << "table `" << nm << "` has a local log-log slope of "
                    << slopeHi_ << " at its upper end, which is too steep to"
                    << " extrapolate safely; limiting to "
                    << sign(slopeHi_)*pMax << "." << nl
                    << "    This usually means a process threshold sits just"
                    << " below the top of the sweep. Extend the E/N range so"
                    << " the region is tabulated rather than extrapolated."
                    << endl;

                slopeHi_ = sign(slopeHi_)*pMax;
            }
        }
    }

    scalar lookup(const scalar x) const
    {
        const label n = x_.size();
        if (n == 0) return 0;
        if (n == 1) return y_[0];

        // BELOW: hold the first value. No extrapolation, by design.
        if (x <= x_[0]) return y_[0];

        if (x >= x_[n-1])
        {
            if (!warned_)
            {
                warned_ = true;
                WarningInFunction
                    << "table `" << name_ << "` evaluated at " << x
                    << ", beyond its range [" << x_[0] << ", " << x_[n-1]
                    << "]." << nl
                    << "    Extrapolating as a power law with the local"
                    << " log-log slope " << slopeHi_ << "." << nl
                    << "    Extend the EEDF sweep if this region matters:"
                    << " extrapolated coefficients are not supported by the"
                    << " Boltzmann solution." << endl;
            }

            if (!logSafe_ || x_[n-1] <= 0) return y_[n-1];

            return y_[n-1]*pow(x/x_[n-1], slopeHi_);
        }

        // Inside: linear interpolation, matching the table's own behaviour.
        label lo = 0, hi = n - 1;
        while (hi - lo > 1)
        {
            const label mid = (lo + hi)/2;
            if (x < x_[mid]) hi = mid; else lo = mid;
        }
        // LOG-LOG interpolation, falling back to linear where a log is
        // undefined (a zero or negative endpoint, which the rate tables have
        // at their bottom end).
        //
        // MEASURED on these tables by dropping every second point and
        // predicting it from its neighbours:
        //
        //     table              span      linear    log-log
        //     muN               5.2e1      0.053%     0.037%
        //     PelasticN         1.5e5      0.035%     0.025%
        //     PinelasticN       5.1e6      0.166%     0.125%
        //     k_EI_N2_ION_N2    2.6e9      1.907%     1.103%
        //
        // At 52 points/decade linear is already good to ~0.1% for the smooth
        // quantities -- the method barely matters there. It matters for the
        // STEEP rate coefficients, which span nine decades and where linear
        // interpolation of a convex function overshoots systematically.
        //
        // It also makes interpolation and EXTRAPOLATION the same method, so
        // there is no change of behaviour at the last table point.
        if (x_[lo] > 0 && x_[hi] > 0 && y_[lo] > 0 && y_[hi] > 0)
        {
            const scalar lx0 = log(x_[lo]), lx1 = log(x_[hi]);
            const scalar dlx = lx1 - lx0;

            if (mag(dlx) > SMALL)
            {
                const scalar p =
                    (log(y_[hi]) - log(y_[lo]))/dlx;

                return y_[lo]*pow(x/x_[lo], p);
            }
        }

        const scalar w = (x - x_[lo])/max(x_[hi] - x_[lo], VSMALL);
        return y_[lo] + w*(y_[hi] - y_[lo]);
    }
};


template<class Base>
class TabulatedProperty1D
:
    public Base
{
    //- The interpolation table (handles file reading and interpolation)
    interpolationTable<scalar> table_;

    //- Bounded lookup: hold below, power-law above. See boundedTable1D.
    boundedTable1D bounded_;

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

        // Bounded lookup replaces interpolationTable's out-of-range
        // behaviour, which returns ZERO rather than extrapolating.
        // See boundedTable1D.
        bounded_.init(table_, propertyName);
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

        // Per cell, through the bounded lookup: interpolate inside the
        // table, hold below it, power-law above it.
        {
            scalarField& f = field.primitiveFieldRef();
            const scalarField& v = var.primitiveField();
            forAll(f, celli)
            {
                f[celli] = bounded_.lookup(v[celli]);
            }
        }

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                scalarField& pField = field.boundaryFieldRef()[patchi];
                const scalarField& pVar = var.boundaryField()[patchi];   // reuse, no re-fetch
                forAll(pField, i)
                {
                    pField[i] = bounded_.lookup(pVar[i]);
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
            pField[i] = bounded_.lookup(pVar[i]);
        }
    }

    //- Scalar lookup for the chemistry-ODE energy coupling. Returns the RAW
    //  table value -- see plasmaPropertyEvaluator::atValue() on why no
    //  post-processing is applied here.
    virtual bool hasScalarLookup() const override { return true; }

    virtual scalar atValue(const scalar x) const override
    {
        return bounded_.lookup(x);
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
//  The gas number density is NOT given here. It is a property of the gas, not
//  of a species' mobility, and it is already declared once in
//  `backgroundGas` -- as `pressure` and `T`, closed by the ideal gas law. It
//  reaches this evaluator through the registered `totalNeutralDensity` field,
//  so every species and every property necessarily uses the same value, and a
//  case cannot end up with two densities that disagree. Reading it from the
//  field also means a spatially varying density (a heated channel, a shock)
//  works with no further change.
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

    //- Bounded lookup: hold below, power-law above. See boundedTable1D.
    boundedTable1D bounded_;

    word lookupVariable_;

    //- Registered field holding the gas number density, used to turn the
    //  reduced quantity into the physical one. Overridable by name only --
    //  never by value, so that the density has exactly one definition.
    word densityField_;

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
        densityField_
        (
            dict.getOrDefault<word>("densityField", "totalNeutralDensity")
        )
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

        // Bounded lookup replaces interpolationTable's out-of-range
        // behaviour, which returns ZERO rather than extrapolating -- fatal
        // for a loss coefficient. See boundedTable1D.
        bounded_.init(table_, this->propertyName());

        // Echo the table's provenance line.
        //
        // These tables are generated -- by the EEDF sweep, or by ionmob from an
        // LXCat download -- and the generator makes choices that change
        // results: which measured dataset was used when several disagree, which
        // background gas it was measured in, what E/N range it covers, and
        // whether its physics checks passed. Recorded in the file header, that
        // information reached only whoever ran the generator. Printing it here
        // puts it in the run log, where the person running the SIMULATION sees
        // it, which is the only place it can actually be reviewed.
        {
            std::ifstream f(path);
            std::string first;
            if (f && std::getline(f, first) && first.compare(0, 2, "//") == 0)
            {
                Info<< "    " << propertyName << ": "
                    << first.substr(2).c_str() << endl;
            }
        }
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

    //- The gas number density field. Fatal if absent: guessing a density would
    //  scale every mobility in the run by an unknown factor while looking
    //  perfectly healthy.
    const volScalarField& gasDensity() const
    {
        if (!this->mesh_.template foundObject<volScalarField>(densityField_))
        {
            FatalErrorInFunction
                << "fromMechanism evaluator '" << this->propertyName()
                << "' cannot find the gas number density field '"
                << densityField_ << "'." << nl
                << "    The tables hold REDUCED quantities (mu*N, D*N), so a"
                << " density is required to recover the physical value." << nl
                << "    It comes from `backgroundGas` in"
                << " plasmaSpeciesProperties -- give `pressure` and `T`, or"
                << " `numberDensity` directly." << nl
                << exit(FatalError);
        }
        return this->mesh_.template lookupObject<volScalarField>(densityField_);
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

        const volScalarField& N = gasDensity();

        scalarField& f = field.primitiveFieldRef();
        const scalarField& x = var.primitiveField();
        const scalarField& n = N.primitiveField();
        forAll(f, celli)
        {
            f[celli] = bounded_.lookup(x[celli])/n[celli];
        }

        // Boundaries take the density from the SAME field's boundary values
        // rather than a cell value standing in for all of them, so a patch
        // sitting in gas at a different density is handled correctly.
        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                scalarField& pF = field.boundaryFieldRef()[patchi];
                const scalarField& pX = var.boundaryField()[patchi];
                const scalarField& pN = N.boundaryField()[patchi];
                forAll(pF, i)
                {
                    pF[i] = bounded_.lookup(pX[i])/pN[i];
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
        const scalarField& pN = gasDensity().boundaryField()[patchi];
        forAll(pField, i)
        {
            pField[i] = bounded_.lookup(pX[i])/pN[i];
        }
    }

    //- Scalar lookup for the chemistry-ODE energy coupling. Returns the RAW
    //  table value -- see plasmaPropertyEvaluator::atValue() on why no
    //  post-processing is applied here.
    virtual bool hasScalarLookup() const override { return true; }

    virtual scalar atValue(const scalar x) const override
    {
        return bounded_.lookup(x);
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
