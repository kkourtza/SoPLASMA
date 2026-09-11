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
    Foam::threeGroupEddington
    
Description

    Three-group Eddington (P1) photoionization model for plasma simulations.

    Derived from the first-order spherical harmonic approximation of the
    radiative transfer equation. For each spectral group j, the xi-weighted
    isotropic photon distribution Psi*_j satisfies:

        laplacian(Psi*_j) - 3*(lambda_j*p)^2 * Psi*_j
            = -(3*lambda_j*p/c) * C * source

    where:
        Psi*_j   : xi-weighted isotropic photon distribution, group j [m^-3]
        lambda_j : inverse absorption length, group j               [m^-1 Pa^-1]
        A_j      : amplitude coefficient, group j                   [m^-1 Pa^-1]
        p        : gas pressure (see note)                          [Pa]
        c        : speed of light                                   [m/s]
        C        : emission prefactor (see sourceMode)              [-]
        source   : registered source field (see sourceMode)         [m^-3 s^-1]

    The photoionization source contribution of each group is:
        S_j = A_j * p * c * Psi*_j

    and the total photoionization source is:
        S_ph = sum_{j=1}^{3} S_j

    The coefficients A_j and lambda_j correspond to the three-group fit of
    g(R)/p directly by exponentials:
        g(R)/p = sum_{j=1}^{3} A_j * exp(-lambda_j * p * R)

    Note: A_j has units [m^-1 Pa^-1] here, unlike the Helmholtz model where
    A_j has units [m^-2 Pa^-2], reflecting the different fitting form.

    Source modes ('sourceMode' keyword):
      - emissionRate   : source = I = (pq/(p+pq))*xi*(nu_u/nu_i)*Siz
                         C = 1. A_j used directly from the literature.
      - ionizationRate : source = Siz (electron impact ionization rate).
                         C = quenchingFactor * xiNuRatio applied internally.
                         A_j used directly from the literature.

    Boundary condition types ('boundaryConditionType' keyword):
      - marshak : Robin condition from radiative transfer theory:
                  grad(Psi*_j) & n = -(3/2)*lambda_j*p * Psi*_j
      - zero    : Zero Dirichlet on all Psi*_j

    Note on pressure:
        p must be consistent with the fitting parameters. For the Bourdon
        2007 air parametrisation, p is the O2 partial pressure
        (p_O2 = 150 Torr at atmospheric pressure).

    Unit conventions ('coefficientUnits' keyword):
      - SI     : p [Pa],   lambda [m^-1 Pa^-1],   A [m^-1 Pa^-1]
      - TorrCm : p [Torr], lambda [cm^-1 Torr^-1], A [cm^-1 Torr^-1]
    All quantities are converted to and stored in SI internally.

Usage
    Example of the model in constant/photoionizationProperties:
    \verbatim
    photoionizationModel    threeGroupEddington;
    solveInterval           10;

    threeGroupEddingtonCoeffs
    {
        coefficientUnits        TorrCm;
        boundaryConditionType   marshak;
        p                       150;            // O2 partial pressure [Torr]

        // Option A: emission rate (I already includes quenching prefactor)
        sourceMode              emissionRate;
        I                       Iph;

        // Option B: ionization rate (prefactor applied internally)
        // sourceMode           ionizationRate;
        // Siz                  Siz;
        // quenchingFactor      0.03797;
        // xiNuRatio            0.06;

        fittingParameters               // Bourdon 2007, Table 3
        (
            (0.0447   0.0067)
            (0.1121   0.0346)
            (0.5994   0.3059)
        );
    }
    \endverbatim

Reference
    Segur P, Bourdon A, Marode E, Bessieres D, and Paillol J H (2006)
    The use of an improved Eddington approximation to facilitate the
    calculation of photoionization in streamer discharges.
    Plasma Sources Science and Technology, 15(4), 648-660.

    Bourdon A, Pasko V P, Liu N Y, Celestin S, Segur P, and Marode E (2007)
    Efficient models for photoionization produced by non-thermal gas
    discharges in air based on radiative transfer and the Helmholtz
    equations. Plasma Sources Science and Technology, 16(3), 656-678.

SourceFiles
    threeGroupEddington.C

\*---------------------------------------------------------------------------*/

#ifndef threeGroupEddington_H
#define threeGroupEddington_H

#include "PtrList.H"
#include "scalarList.H"
#include "Tuple2.H"
#include "volFields.H"

#include "photoionizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class threeGroupEddington Declaration
\*---------------------------------------------------------------------------*/

class threeGroupEddington
:
    public photoionizationModel
{
private:

    //- Disallow copy
    threeGroupEddington
    (
        const threeGroupEddington&
    ) = delete;

    void operator=(const threeGroupEddington&) = delete;

    // Private Member Data

        //- Unit convention ("SI" or "TorrCm")
        word coefficientUnits_;

        //- Boundary condition type ("marshak" or "zero")
        word boundaryConditionType_;

        //- Number of spectral groups (fixed = 3)
        static const label nGroups_ = 3;

       //- Gas pressure consistent with the fitting parameters [Pa]
        dimensionedScalar p_;

        //- Source field mode ("emissionRate" or "ionizationRate")
        word sourceMode_;

        //- Name of the source field in the object registry
        word sourceName_;

        //- Quenching factor pq/(p+pq) [-]; ionizationRate mode only
        scalar quenchingFactor_;

        //- Photoionization efficiency xi*nu_u/nu_i[-]; ionizationRate mode only
        scalar xiNuRatio_;

        //- Combined prefactor C = quenchingFactor_ * xiNuRatio_ [-]
        //  C = 1 in emissionRate mode
        scalar C_;

        //- Amplitude coefficients A_j [m^-1 Pa^-1]
        scalarList A_;

        //- Inverse absorption lengths lambda_j [m^-1 Pa^-1]
        scalarList lambda_;

        //- Speed of light [m/s]
        dimensionedScalar c_;

        //- Photon distribution components Psi*_j, one per group [m^-3]
        PtrList<volScalarField> PsiStar_j_;

        //- Partial photoionization source fields S_j, one per group [m^-3 s^-1]
        PtrList<volScalarField> Sph_j_;

public:

    //- Runtime type information
    TypeName("threeGroupEddington");

    // Constructors

        //- Construct from mesh
        threeGroupEddington(const fvMesh& mesh);

    //- Destructor
    virtual ~threeGroupEddington() = default;


    // Member Functions

        //- Solve 3 Eddington equations & sum into total photoionization source
        virtual void correct() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

