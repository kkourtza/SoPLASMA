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
    Foam::nTermHelmholtz
    
Description

    N-term Helmholtz photoionization model for plasma simulations in air.

    Approximates the non-local Zheleznyak photoionization integral by
    decomposing the absorption kernel into N exponential terms, yielding
    N independent screened-Poisson (Helmholtz) PDEs
    
        laplacian(Sph_j) - (lambda_j*p)^2 * Sph_j = -A_j * p^2 * C * source

    where:
        Sph_j    : photoionization source component, term j  [m^-3 s^-1]
        lambda_j : inverse absorption length, term j         [m^-1 Pa^-1]
        A_j      : amplitude coefficient, term j             [m^-2 Pa^-2]
        p        : gas pressure (see note)                   [Pa]
        C        : emission prefactor (see sourceMode)       [-]
        source   : registered source field (see sourceMode)  [m^-3 s^-1]

    The total photoionization source is:
        Sph = sum_{j=1}^{N} Sph_j

    The coefficients A_j and lambda_j are obtained by fitting:
        g(R)/p / (p*R) = sum_{j=1}^{N} A_j * exp(-lambda_j * p * R)
    to the Zheleznyak absorption function.

    Source modes ('sourceMode' keyword):
      - emissionRate   : source = I = (pq/(p+pq))*xi*(nu_u/nu_i)*Siz
                         C = 1. A_j used directly from the literature.
      - ionizationRate : source = Siz (electron impact ionization rate).
                         C = quenchingFactor * xiNuRatio applied internally.
                         A_j used directly from the literature.

    Note on pressure:
        p must be consistent with the fitting parameters. For the Bourdon
        2007 air parametrisation, p is the O2 partial pressure
        (p_O2 = 150 Torr at atmospheric pressure).

    Unit conventions ('coefficientUnits' keyword):
      - SI     : p [Pa],   lambda [m^-1 Pa^-1],   A [m^-2 Pa^-2]
      - TorrCm : p [Torr], lambda [cm^-1 Torr^-1], A [cm^-2 Torr^-2]
    All quantities are converted to and stored in SI internally.

Usage
    Example of the model in constant/photoionizationProperties:
    \verbatim
    photoionizationModel    nTermHelmholtz;
    solveInterval           10;

    nTermHelmholtzCoeffs
    {
        coefficientUnits    TorrCm;
        p                   150;            // O2 partial pressure [Torr]

        // Option A: emission rate (I already includes quenching prefactor)
        sourceMode          emissionRate;
        I                   Iph;

        // Option B: ionization rate (prefactor applied internally)
        // sourceMode       ionizationRate;
        // Siz              Siz;
        // quenchingFactor  0.03797;        // pq/(p+pq) = 30/(760+30)
        // xiNuRatio        0.06;           // xi * nu_u/nu_i

        nTerms              3;
        fittingParameters               // Bourdon 2007, Table 2
        (
            (0.0553   1.986e-4)
            (0.1460   5.100e-3)
            (0.8900   4.886e-1)
        );
    }
    \endverbatim

Reference
    Zheleznyak M B, Mnatsakanyan A Kh, and Sizykh S V (1982)
    Photoionization of nitrogen and oxygen mixtures by radiation
    from a gas discharge. High Temperature Science, 20(3), 357-362.

    Bourdon A, Pasko V P, Liu N Y, Celestin S, Segur P, and Marode E (2007)
    Efficient models for photoionization produced by non-thermal gas
    discharges in air based on radiative transfer and the Helmholtz
    equations. Plasma Sources Science and Technology, 16(3), 656-678.    

SourceFiles
    nTermHelmholtz.C

\*---------------------------------------------------------------------------*/

#ifndef nTermHelmholtz_H
#define nTermHelmholtz_H

#include "PtrList.H"
#include "scalarList.H"
#include "Tuple2.H"
#include "volFields.H"

#include "photoionizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class nTermHelmholtz Declaration
\*---------------------------------------------------------------------------*/

class nTermHelmholtz
:
    public photoionizationModel
{
private:

    //- Disallow copy and assignment
    nTermHelmholtz
    (
        const nTermHelmholtz&
    ) = delete;

    void operator=(const nTermHelmholtz&) = delete;

    // Private Member Data

        //- Unit convention ("SI" or "TorrCm")
        word coefficientUnits_;

        //- Number of Helmholtz terms
        label nTerms_;

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

        //- Amplitude coefficients A_j [m^-2 Pa^-2]
        scalarList A_;

        //- Inverse absorption lengths lambda_j [m^-1 Pa^-1]
        scalarList lambda_;

        //- Partial photoionization source fields, one per Helmholtz term
        PtrList<volScalarField> Sph_j_;

public:

    //- Runtime type information
    TypeName("nTermHelmholtz");

    // Constructors

        //- Construct from mesh
        nTermHelmholtz(const fvMesh& mesh);

    //- Destructor
    virtual ~nTermHelmholtz() = default;

    // Member Functions

        //- Solve N Helmholtz equations & sum into total photoionization source
        virtual void correct() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

