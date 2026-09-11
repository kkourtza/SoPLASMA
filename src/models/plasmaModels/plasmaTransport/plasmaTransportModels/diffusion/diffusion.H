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
    Foam::diffusion

Description
    Transport by diffusion alone, for species that do not drift:

        dn/dt = div(D grad(n)) + Sources

    UNCHARGED SPECIES ONLY, and the constructor enforces it -- for the simple
    reason that this model has NO DRIFT TERM, and a charged species drifts.
    There is nothing subtler to it than that.

    In particular this is not a statement about flux schemes. `driftDiffusion`
    offers both treatments, and which one is used is the case's choice:

      fluxScheme standard             drift and diffusion are SEPARATE terms,
                                      fvm::div(phi, n) and fvm::laplacian(D, n),
                                      so the convective scheme is whatever
                                      fvSchemes names for div(phi_s,n_s) -- the
                                      ROUND family included -- and the diffusion
                                      scheme is chosen independently.

      fluxScheme ScharfetterGummel    the two are combined into a single
                                      exponentially weighted face flux, which is
                                      what stays stable when drift dominates.
                                      They cannot be separated in this mode, and
                                      splitSGFlux() only recovers the pieces for
                                      diagnostics.

    So drift and diffusion are separable when the case asks for it. What a
    charged species cannot do is omit the drift term altogether, which is what
    choosing this model would mean.

    Between this and `immobile` the question is physical, not numerical:
    compare the diffusion length sqrt(D*tau) over the species' CHEMICAL
    lifetime with the scale of the gradient it sits in. At 1 bar, D ~ 2e-5
    m^2/s, so a state living 1 ns spreads 0.14 um -- far below a streamer
    head, and `immobile` is then correct rather than merely cheap. A radical
    living 1 ms spreads 140 um and must diffuse. `plasmaTransport` reports both
    numbers per species at start-up so the choice can be checked rather than
    assumed.

    The model is activated by specifying:

        transportModel diffusion;

    in the corresponding species dictionary, with a `diffusivity` sub-dictionary
    exactly as `driftDiffusion` takes.

SourceFiles
    diffusion.C

\*---------------------------------------------------------------------------*/

#ifndef diffusion_H
#define diffusion_H

#include "plasmaTransportModel.H"
#include "plasmaDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class diffusion Declaration
\*---------------------------------------------------------------------------*/

class diffusion
:
    public plasmaTransportModel
{
private:

    // Private data

        //- Diffusivity model
        autoPtr<plasmaDiffusivityModel> diffusivityModel_;

    //- Disallow default bitwise copy construct
    diffusion(const diffusion&) = delete;
    void operator=(const diffusion&) = delete;

public:

    //- Runtime type information
    TypeName("diffusion");

    // Constructors

        //- Construct from components
        diffusion
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            plasmaSpecies& species,
            const label specieIndex
        );

    //- Destructor
    virtual ~diffusion() = default;

    // Member Functions

        //- Update the diffusivity
        virtual void correct() override;

        //- Return ddt(n) - laplacian(D, n)
        virtual tmp<fvScalarMatrix> nEqn() const override;

        //- Diffusive flux only; the convective flux is identically zero
        virtual void updateFluxes
        (
            const fvScalarMatrix& nEqnMatrix,
            surfaceScalarField& convectiveFlux,
            surfaceScalarField& diffusiveFlux,
            surfaceScalarField& particleFlux
        ) const override;

        //- Zero: the species is uncharged, enforced at construction
        virtual tmp<volScalarField> electricalConductivity() const override;

        //- Zero, for the same reason
        virtual tmp<volScalarField> diffusiveChargeSource() const override;

        //- Access to the diffusivity, for the start-up advisory report
        const plasmaDiffusivityModel& diffusivity() const
        {
            return *diffusivityModel_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
