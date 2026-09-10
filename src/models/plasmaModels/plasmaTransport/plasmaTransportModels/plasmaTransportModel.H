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
    Foam::plasmaTransportModel

Description
    Abstract base class for per-species plasma transport behaviour.

    This class defines the minimal interface required for describing how a
    single plasma species is transported within the simulation.

    The plasmaTransportModel does not implement any transport physics itself.
    Instead, it serves as a common framework from which concrete transport
    models are derived. Each derived class implements its own formulation
    for updating the transport state of a species during the solution loop.

    The run-time selection mechanism is used to construct the appropriate
    model for each species based on the 'transportModel' entry within the
    species dictionary.

SourceFiles
    plasmaTransportModelI.H
    plasmaTransportModel.C
    
\*---------------------------------------------------------------------------*/

#ifndef plasmaTransportModel_H
#define plasmaTransportModel_H

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
                     Class plasmaTransportModel Declaration
\*---------------------------------------------------------------------------*/

class plasmaTransportModel
{
private:

    //- Disallow default bitwise copy construct
    plasmaTransportModel(const plasmaTransportModel&) = delete;
    void operator=(const plasmaTransportModel&) = delete;

protected:

    //- Model name (e.g. driftDiffusion, fullMomentum, immobile, etc.)
    const word modelName_;

    //- Mesh reference
    const fvMesh& mesh_;

    //- Reference to species database
    plasmaSpecies& species_;

    //- Dictionary for this species
    const dictionary& dict_;

    //- Index of the species
    const label specieIndex_;
    
public:

    //- Runtime type information
    TypeName("plasmaTransportModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,                        // storage class
        plasmaTransportModel,           // base class
        dictionary,                     // lookup keyword
        (                               // constructor signature
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            plasmaSpecies& species,
            const label specieIndex
        ),
        (modelName, dict, mesh, species, specieIndex)
    );

    // Constructors

        //- Construct base model
        plasmaTransportModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            plasmaSpecies& species,
            const label specieIndex
        );

    //- Selector
    static autoPtr<plasmaTransportModel> New
    (
        const word& modelName,
        const dictionary& dict,
        const fvMesh& mesh,
        plasmaSpecies& species,
        const label specieIndex
    );

    //- Destructor
    virtual ~plasmaTransportModel() = default;

    // Member Functions

        //- Update transport coefficients
        virtual void correct() = 0;

        //- Return the transport equation (Matrix)
        virtual tmp<fvScalarMatrix> nEqn() const = 0;

        //- Hand this species' NET volumetric source to the transport model.
        //
        //  VIRTUAL WITH A NO-OP DEFAULT, so plasmaTransport can call it on
        //  every species without knowing which model or which flux scheme is
        //  underneath -- the same reason updateFluxes() is virtual here rather
        //  than reached through a driftDiffusion cast.
        //
        //  Only the CompleteFlux scheme reads it: CFS derives its face flux
        //  from the local two-point problem WITH the source retained, where
        //  Scharfetter-Gummel sets it to zero. Every other model and scheme
        //  ignores this.
        virtual void setNetSource(const volScalarField&) const
        {}

        //- Update the fluxes in place (1/s)
        virtual void updateFluxes
        (
            const fvScalarMatrix& nEqnMatrix,
            surfaceScalarField& convectiveFlux,
            surfaceScalarField& diffusiveFlux,
            surfaceScalarField& particleFlux
        ) const = 0;

        //- Return electrical conductivity contribution [S/m]
        virtual tmp<volScalarField> electricalConductivity() const = 0;

        //- Return diffusive charge source contribution [C/m3/s]
        virtual tmp<volScalarField> diffusiveChargeSource() const = 0;

        //- Particle mobility [driftDiffusion] / diffusivity, for models that
        //  have one. Default FatalErrors, naming the model, rather than
        //  silently returning a meaningless field -- for a model with no
        //  drift/diffusion coefficient (immobile) there is nothing correct
        //  to return.
        //
        //  Added so an explicit (fvc::) residual for the JFNK outer solver
        //  (docs/design/newton-outer-solver-design.md) can reconstruct THIS species'
        //  convective/diffusive flux at an ARBITRARY trial state, exactly
        //  matching nEqn()'s own formula, without calling nEqn()+.residual()
        //  itself (confirmed broken in parallel for that specific usage --
        //  see [[fvmatrix-residual-broken-in-parallel-use-fvc-instead]]).
        virtual tmp<volScalarField> mu() const
        {
            FatalErrorInFunction
                << "Model '" << modelName_ << "' has no particle mobility"
                << " (mu()) -- the Newton outer solver does not yet support"
                << " this transport model." << nl << exit(FatalError);
            return nullptr;
        }
        virtual tmp<volScalarField> D() const
        {
            FatalErrorInFunction
                << "Model '" << modelName_ << "' has no diffusivity"
                << " (D()) -- the Newton outer solver does not yet support"
                << " this transport model." << nl << exit(FatalError);
            return nullptr;
        }

    // Access

        inline const word& modelName() const;
        inline const fvMesh& mesh() const;
        inline const plasmaSpecies& species() const;
        inline label specieIndex() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaTransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

