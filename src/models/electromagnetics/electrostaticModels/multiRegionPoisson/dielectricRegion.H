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
    Foam::dielectricRegion
    
Description
    A container class for dielectric region fields and properties.

    It is designed to be used in multi-region Poisson solvers to eliminate the
    need of multiple synchronized pointer lists.

    Holds per-region storage for ePotential, E, phiE, and permittivity. Each
    region uses a single uniform permittivity (no spatial variation), passed in
    by the owning model, which reads it from
    constant/<region>/electricalProperties.

\*---------------------------------------------------------------------------*/

#ifndef dielectricRegion_H
#define dielectricRegion_H

#include "fvCFD.H"

#include "plasmaConstants.H"
#include "materialLibrary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dielectricRegion Declaration
\*---------------------------------------------------------------------------*/

class dielectricRegion
{

protected:
   
    //- Reference to the mesh for this specific dielectric region
    const fvMesh& regionMesh_;

    //- Electric potential field [V]
    volScalarField ePotential_;

    //- Electric field [V/m]
    volVectorField E_;

    //- Flux of the electric field [V * m2]
    surfaceScalarField phiE_;

    //- Surface charge density (owned by dielectric) [C/m2]
    volScalarField surfCharge_;

    //- Dielectric permittivity [F/m]
    dimensionedScalar epsilon_;

    //- Relative permittivity of the medium
    scalar epsilonR_;

    //- Numerical scheme for computing E
    word EScheme_;

public:

    // Constructors

        //- Construct from mesh, relative permittivity and E scheme.
        //
        //  Both are passed in rather than read here: epsilonR because its file
        //  location is the owning model's business (and has a legacy fallback),
        //  EScheme because it is a GLOBAL numerical setting -- forming E by a
        //  different scheme on either side of an interface would make the
        //  coupled flux inconsistent, so it cannot be per-region.
        dielectricRegion
        (
            const fvMesh& regionMesh,
            const scalar epsilonR,
            const word& EScheme
        )
        :
            regionMesh_(regionMesh),
            ePotential_
            (
                IOobject
                (
                    "ePotential",
                    regionMesh_.time().timeName(),
                    regionMesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                regionMesh_
            ),
            E_
            (
                IOobject
                (
                    "E",
                    regionMesh_.time().timeName(),
                    regionMesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                regionMesh_,
                dimensionedVector
                (
                    "zero",
                    dimensionSet(1, 1, -3, 0, 0, -1, 0),
                    vector::zero
                )
            ),
            phiE_
            (
                IOobject
                (
                    "phiE",
                    regionMesh_.time().timeName(),
                    regionMesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                -fvc::snGrad(ePotential_) * regionMesh_.magSf()
            ),
            surfCharge_
            (
                IOobject
                (
                    "surfCharge",
                    regionMesh_.time().timeName(),
                    regionMesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                regionMesh_,
                dimensionedScalar
                (
                    "zero",
                    dimensionSet(0, -2, 1, 0, 0, 1, 0),
                    0.0
                )
            ),
            epsilon_
            (
                "epsilon",
                dimensionSet(-1, -3, 4, 0, 0, 2, 0),
                0.0
            ),
            epsilonR_(epsilonR),
            EScheme_(EScheme)
        {
            epsilon_ = dimensionedScalar("epsilon", epsilonR_ * constant::plasma::epsilon0);

            // EScheme is validated once, globally, by
            // electromagneticsModel::readPoissonNumerics().

            // Register a lightweight in-memory dictionary in this region's mesh
            // registry. This is NOT a file: it is how the coupled-potential BC
            // finds the permittivity of the region it lives on, since a patch
            // field can reach its own registry but not the owning model.
            IOdictionary* epsilonDict = new IOdictionary
            (
                IOobject
                (
                    "dielectricProperties",
                    regionMesh_.time().constant(),
                    regionMesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                )
            );
            epsilonDict->add("epsilonR", epsilonR_);

            // GAMMA TRAVELS WITH THE MATERIAL, not with the patch.
            //
            // Secondary emission IS a surface property, unlike permittivity --
            // but for a MESHED barrier the emitting surface IS this region, so
            // its gamma belongs to the region's material declaration and not to
            // a boundary condition on the far side of the interface.
            //
            // Registered here so the gas-side wall-flux condition can read it
            // through the mapped patch, which is what stops gamma being stated
            // twice: once as the material and once as `defaultSEEC`.
            //
            // Read from the SAME file as epsilonR. Absent means "not stated",
            // and the wall condition keeps its own default rather than being
            // handed a fabricated value.
            {
                IOdictionary props
                (
                    IOobject
                    (
                        "electricalProperties",
                        regionMesh_.time().constant(),
                        regionMesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    )
                );

                if (props.found("gammaSEE"))
                {
                    // An explicit value WINS: it is an override.
                    epsilonDict->add
                    (
                        "gammaSEE", props.get<scalar>("gammaSEE")
                    );
                }
                else if (props.found("material"))
                {
                    // NAMED MATERIAL (G2). Gamma travels with the material,
                    // and the material is stated ONCE -- the same entry that
                    // supplied epsilonR. Resolved through the SAME library as
                    // epsilonR, so the two cannot come from different places.
                    const word mat(props.get<word>("material"));

                    epsilonDict->add
                    (
                        "gammaSEE",
                        materialLibrary::get
                        (
                            mat,
                            "gammaSEE",
                            "constant/" + regionMesh_.dbDir()
                          + "/electricalProperties (region `"
                          + regionMesh_.name() + "`)"
                        )
                    );
                }
            }

            regionMesh_.thisDb().store(epsilonDict);
        }

    //- Destructor
    virtual ~dielectricRegion() = default;

    // Member functions

        //- Update the electric field
        void updateE()
        {
            phiE_ = -fvc::snGrad(ePotential_) * regionMesh_.magSf();

            if (EScheme_ == "reconstruct")
            {
                E_ = fvc::reconstruct(phiE_);
            }
            else // "grad"
            {
                E_ = -fvc::grad(ePotential_);
            }
            E_.correctBoundaryConditions();
        }

        //- Return the LHS of the Poisson equation
        tmp<fvScalarMatrix> PoissonMatrix() const
        {
            return fvm::laplacian(epsilon_, ePotential_);
        }

        tmp<fvScalarMatrix> PoissonEquation() const
        {
            // In a dielectric, the Poisson equation reduces to Laplace
            // (no free charges in the bulk; surface charges enter via BCs)
            return PoissonMatrix();
        }

    // Access

        inline const fvMesh& mesh() const { return regionMesh_; }
        inline const volScalarField& ePotential() const { return ePotential_; }
        inline volScalarField& ePotential() { return ePotential_; }
        inline const volVectorField& E() const { return E_; }
        inline const surfaceScalarField& phiE() const { return phiE_; }
        inline const volScalarField& surfCharge() const { return surfCharge_; }
        inline volScalarField& surfCharge() { return surfCharge_; }
        inline const dimensionedScalar& epsilon() const { return epsilon_; }
        inline scalar epsilonR() const { return epsilonR_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

