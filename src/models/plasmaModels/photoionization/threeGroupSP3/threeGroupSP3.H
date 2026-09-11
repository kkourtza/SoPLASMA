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
    Foam::threeGroupSP3

Description

    Three-group SP3 (third-order simplified spherical harmonics) photoionization
    model for plasma simulations.

    Improves on the Eddington model when absorption is weak and source gradients
    are large. For each spectral group j, two auxiliary fields phi1_j and phi2_j
    satisfy independent Helmholtz PDEs:

        laplacian(phi1_j) - (lambda_j*p)^2/kappa1^2 * phi1_j
            = -(lambda_j*p)/(kappa1^2*c) * C * source

        laplacian(phi2_j) - (lambda_j*p)^2/kappa2^2 * phi2_j
            = -(lambda_j*p)/(kappa2^2*c) * C * source

    where:
        phi1_j, phi2_j : auxiliary photon distribution components, group j[m^-3]
        lambda_j       : inverse absorption length, group j         [m^-1 Pa^-1]
        A_j            : amplitude coefficient, group j             [m^-1 Pa^-1]
        p              : gas pressure (see note)                    [Pa]
        c              : speed of light                             [m/s]
        C              : emission prefactor (see sourceMode)        [-]
        source         : registered source field (see sourceMode)   [m^-3 s^-1]
        kappa1^2       = 3/7 - (2/7)*sqrt(6/5)
        kappa2^2       = 3/7 + (2/7)*sqrt(6/5)

    The photoionization source contribution of each group is recovered via
    the SP3 closure:
        Psi_j = (gamma2*phi1_j - gamma1*phi2_j) / (gamma2 - gamma1)
        S_j   = A_j * p * c * Psi_j

    where:
        gamma_n = (5/7) * [1 + (-1)^n * 3*sqrt(6/5)]

    and the total source is:
        S_ph = sum_{j=1}^{3} S_j

    The SP3 model is mathematically equivalent to a six-exponential Helmholtz
    model and generally provides better agreement with the Zheleznyak integral
    than either the three-exponential Helmholtz or three-group Eddington models,
    particularly near domain boundaries and at large p_O2*R values.

    The coefficients A_j and lambda_j are the same as for the Eddington model:
        g(R)/p = sum_{j=1}^{3} A_j * exp(-lambda_j * p * R)

    Source modes ('sourceMode' keyword):
      - emissionRate   : source = I = (pq/(p+pq))*xi*(nu_u/nu_i)*Siz
                         C = 1. A_j used directly from the literature.
      - ionizationRate : source = Siz (electron impact ionization rate).
                         C = quenchingFactor * xiNuRatio applied internally.
                         A_j used directly from the literature.

    Boundary condition types ('boundaryConditionType' keyword):
      - larsen : Coupled Robin BC from radiative transfer theory (Liu 2007):
                 grad(phi1_j) & n = -lambda_j*p*(alpha1*phi1_j + beta2*phi2_j)
                 grad(phi2_j) & n = -lambda_j*p*(alpha2*phi2_j + beta1*phi1_j)
                 Solved via Picard iteration on the beta cross-coupling.
      - zero   : Zero Dirichlet on both phi1_j and phi2_j

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
    photoionizationModel    threeGroupSP3;
    solveInterval           10;

    threeGroupSP3Coeffs
    {
        coefficientUnits        TorrCm;
        boundaryConditionType   larsen;
        larsenMaxIters          10;
        larsenTol               1e-6;
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

    Liu N Y, Celestin S, Bourdon A, Pasko V P, Segur P, and Marode E (2007)
    Application of photoionization models based on radiative transfer and
    the Helmholtz equations to studies of streamers in weak electric
    fields. Appl. Phys. Lett., 91 (21), 211501.

    Larsen E W, Thommes G, Klar A, Seaid M, and Gotz T (2002)
    Simplified PN approximations to the equations of radiative heat
    transfer and applications. J. Comput. Phys., 183 (2), 652-675.

SourceFiles
    threeGroupSP3.C

\*---------------------------------------------------------------------------*/

#ifndef threeGroupSP3_H
#define threeGroupSP3_H

#include "PtrList.H"
#include "scalarList.H"
#include "Tuple2.H"
#include "volFields.H"

#include "photoionizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class threeGroupSP3 Declaration
\*---------------------------------------------------------------------------*/

class threeGroupSP3
:
    public photoionizationModel
{
private:

    //- Disallow copy
    threeGroupSP3(const threeGroupSP3&) = delete;

    void operator=(const threeGroupSP3&) = delete;

    // Private Member Data

        //- Unit convention ("SI" or "TorrCm")
        word coefficientUnits_;

        //- Boundary condition type ("larsen" or "zero")
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

        //- SP3 closure constants kappa1^2 = 3/7 - (2/7)*sqrt(6/5)
        static const scalar kappa1Sqr_;

        //- SP3 closure constants kappa2^2 = 3/7 + (2/7)*sqrt(6/5)
        static const scalar kappa2Sqr_;

        //- SP3 closure weight gamma1 = (5/7)*[1 - 3*sqrt(6/5)]
        static const scalar gamma1_;

        //- SP3 closure weight gamma2 = (5/7)*[1 + 3*sqrt(6/5)]
        static const scalar gamma2_;

        //- Larsen BC coefficient alpha1 = (5/96)*(34 + 11*sqrt(6/5))
        static const scalar alpha1_;

        //- Larsen BC coefficient alpha2 = (5/96)*(34 - 11*sqrt(6/5))
        static const scalar alpha2_;

        //- Larsen BC coefficient beta1 = (5/96)*(2 - sqrt(6/5))
        static const scalar beta1_;

        //- Larsen BC coefficient beta2 = (5/96)*(2 + sqrt(6/5))
        static const scalar beta2_;

        //- Maximum Picard iterations on Larsen beta cross-coupling
        label larsenMaxIters_;

        //- Relative tolerance for Larsen Picard convergence
        scalar larsenTol_;

        //- Auxiliary photon distribution phi1_j, one per group [m^-3]
        PtrList<volScalarField> phi1_j_;

        //- Auxiliary photon distribution phi2_j, one per group [m^-3]
        PtrList<volScalarField> phi2_j_;

        //- Partial photoionization source fields S_j, one per group [m^-3 s^-1]
        PtrList<volScalarField> Sph_j_;

public:

    //- Runtime type information
    TypeName("threeGroupSP3");

    // Constructors

        //- Construct from mesh
        threeGroupSP3(const fvMesh& mesh);

    //- Destructor
    virtual ~threeGroupSP3() = default;


    // Member Functions

        //- Solve 6 SP3 sub-equations (2 per group), apply SP3 closure,
        //  and sum into total photoionization source
        virtual void correct() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

