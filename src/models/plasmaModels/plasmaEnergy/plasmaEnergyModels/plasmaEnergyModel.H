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
    Foam::plasmaEnergyModel

Description
    Abstract base class for per-species plasma energy behaviour.

    This class defines the minimal interface required for describing how the 
    energy in a single plasma species is transported within the simulation.

    The plasmaEnergyModel does not implement any energy physics itself.
    Instead, it serves as a common framework from which concrete energy
    models are derived.  Each derived class implements its own formulation
    for updating the energy state of a species during the solution loop.

    The run-time selection mechanism is used to construct the appropriate
    model for each species based on the 'energyModel' entry within the
    species dictionary.

SourceFiles
    plasmaEnergyModelI.H
    plasmaEnergyModel.C
    
\*---------------------------------------------------------------------------*/

#ifndef plasmaEnergyModel_H
#define plasmaEnergyModel_H

#include "DynamicList.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

#include "plasmaSpecies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class plasmaEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class plasmaEnergyModel
{
protected:

    //- Model name (e.g. isothermal, localField, localEnergy)
    const word modelName_;

    //- Mesh reference
    const fvMesh& mesh_;

    //- Reference to species database
    const plasmaSpecies& species_;

    //- Dictionary for this species
    const dictionary& dict_;

    //- Index of the species this model applies to
    const label specieIndex_;

    //- Reference to the electric field
    const volVectorField& E_;

public:

    //- Runtime type information
    TypeName("plasmaEnergyModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,                        // storage class
        plasmaEnergyModel,           // base class
        dictionary,                     // lookup keyword
        (                               // constructor signature
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex,
            const volVectorField& E
        ),
        (modelName, dict, mesh, species, specieIndex, E)
    );

    // Constructors

        //- Construct base model
        plasmaEnergyModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex,
            const volVectorField& E
        );

    //- Selector
    static autoPtr<plasmaEnergyModel> New
    (
        const word& modelName,
        const dictionary& dict,
        const fvMesh& mesh,
        const plasmaSpecies& species,
        const label specieIndex,
        const volVectorField& E
    );

    //- Destructor
    virtual ~plasmaEnergyModel() {}

    // Member Functions

        //- Correct model
        virtual void correct() = 0;

        //- The electron-energy density field, when this model carries one.
        //  Returns nullptr otherwise (e.g. the local-field model, which has
        //  no transported energy). Used by plasmaOuterRelaxation to enrol the
        //  energy in the JOINT outer-loop residual: it is one half of the
        //  period-2 cycle relaxation exists to damp, and relaxing the density
        //  alone was measured to floor the factor without converging.
        virtual volScalarField* nEpsPtr() { return nullptr; }

        //- Throw away the step just solved and return to oldTime().
        //
        //  A model that SOLVES a transported field must restore it here, or a
        //  retried step re-solves from a mixed state: plasmaTransport rolls
        //  the species back to t^n while this model's field is still at the
        //  discarded t^{n+1}. Under LMEA that mismatch is formed into
        //  meanE = nEps_e/n_e, and it SIGFPE'd on the first discard.
        //
        //  Default no-op: a model with no transported state of its own
        //  (isothermal, gasTemperature, localField) has nothing to undo.
        virtual void discardStep() {}

        //- Append every field this model TRANSPORTS to `fields`.
        //
        //  Default no-op, for the same reason discardStep() above is one: a
        //  model with no transported state of its own (isothermal,
        //  gasTemperature, localField) has nothing to contribute.
        //
        //  WHY THIS MIRRORS discardStep(). Whatever must be RESTORED when a
        //  step is thrown away is exactly the transported state, so the two
        //  lists must never diverge. Keeping them adjacent and identically
        //  shaped is deliberate: a model that adds a transported field and
        //  forgets one of the two will be obvious. retryStep shipped with
        //  nEps_e missing from discardStep() and it took a SIGFPE at step 6 to
        //  find, because every validation bed was LFA.
        virtual void appendTransportedFields
        (
            DynamicList<const volScalarField*>& fields,
            DynamicList<word>& names
        ) const
        {}

        //- Return the energy equation (Matrix)
        virtual tmp<fvScalarMatrix> eEqn() const = 0;

    // Access
    
        //- Returns true if the temperature is uniform (constant in space)
        virtual bool isUniform() const = 0;

        //- Returns true if the temperature is constant in time
        virtual bool isConstant() const = 0;

        //- Return a tmp T field
        //- It can be a reference to a real field or a created uniform field
        virtual tmp<volScalarField> T() const = 0;

        //- Return the scalar value for isothermal models.
        virtual const dimensionedScalar& Tvalue() const;


        inline const word& modelName() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaEnergyModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

