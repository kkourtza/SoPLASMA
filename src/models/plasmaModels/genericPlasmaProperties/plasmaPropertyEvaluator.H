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
    Foam::plasmaPropertyEvaluator

Description:
    Base class for mathematical property evaluation.

    The plasmaPropertyEvaluator class provides a unified interface for 
    mathematical models (interpolation, analytical functions, etc.) that are 
    decoupled from specific physical identities. 

SourceFiles:
    plasmaPropertyEvaluator.C
\*---------------------------------------------------------------------------*/

#ifndef plasmaPropertyEvaluator_H
#define plasmaPropertyEvaluator_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class plasmaPropertyEvaluator Declaration
\*---------------------------------------------------------------------------*/

class plasmaPropertyEvaluator
{
protected:

    // Protected Data

        //- The specific math type (e.g., "tabulated1D")
        const word modelType_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- The name of the evaluated field (for prints)
        const word propertyName_;

        //- The dimensionSet of the evaluated field
        const dimensionSet dims_;

public:

    //- Runtime type information
    TypeName("plasmaPropertyEvaluator");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,
        plasmaPropertyEvaluator,
        dictionary,
        (
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const word& propertyName,
            const dimensionSet& dims
        ),
        (modelType, dict, mesh, propertyName, dims)
    );

    // Constructors

        //- Construct from components
        plasmaPropertyEvaluator
        (
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const word& propertyName,
            const dimensionSet& dims
        );

    //- Selector
    static autoPtr<plasmaPropertyEvaluator> New
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const word& propertyName,
        const dimensionSet& dims
    );

    //- Destructor
    virtual ~plasmaPropertyEvaluator() = default;

    // Member Functions

        //- Returns true if the property is spatially uniform
        virtual bool isUniform() const = 0;

        //- Returns true if the property is constant in time
        virtual bool isConstant() const = 0;

        // Generate and return a temporary field property
        virtual tmp<volScalarField> evaluate() const = 0;

        //- Return the scalar value for uniform field properties
        virtual dimensionedScalar evaluateScalar() const = 0;

        //- Update an existing field in place
        virtual void correct(volScalarField& field) const = 0;

        //- Whether this property can be evaluated at a SCALAR lookup value,
        //  outside the field machinery.
        //
        //  Needed because the electron energy is integrated inside the
        //  chemistry ODE (docs/lmea-option4-plan.md), where the coefficients
        //  must be re-evaluated at an eps that is CHANGING during the substep.
        //  A field-at-a-time evaluator cannot serve that: there is no field,
        //  only the integrator's current state.
        //
        //  Only the tabulated 1-D forms support it. Everything else returns
        //  false, and asking anyway is a programming error rather than a
        //  recoverable condition.
        virtual bool hasScalarLookup() const { return false; }

        //- The RAW tabulated value at `x`, in the units the TABLE carries.
        //
        //  NOT post-processed: `fromMechanism` divides by the gas density when
        //  it fills a field, and this does NOT. The caller therefore receives
        //  mu*N for a `muN` table, and is responsible for the division. That is
        //  deliberate -- the division needs a gas density, which changes with
        //  the gas temperature, and capturing one in a lambda is how a stale
        //  value gets baked in.
        virtual scalar atValue(const scalar x) const
        {
            FatalErrorInFunction
                << "property '" << propertyName_ << "' (type " << modelType_
                << ") has no scalar lookup; check hasScalarLookup() first."
                << exit(FatalError);
            return 0;
        }

        //- Evaluate for a patch field
        virtual void evaluatePatch
        (
            scalarField& pField, 
            const label patchi
        ) const = 0;



        //- Interface compatibility: Keyword to look for in 'constant' type
        static word coeffKeyword() { return "value"; }

        //- Interface compatibility: Returns the physical model category
        static word propertyType() { return "evaluator"; }

        //- Access to model type string
        word modelType() const { return modelType_; }

        //- Access the name of the field that is evaluated (for prints)
        const word& propertyName() const { return propertyName_; }

        //- Access the dimensionSet of the evaluated field
        const dimensionSet& dims() const { return dims_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
