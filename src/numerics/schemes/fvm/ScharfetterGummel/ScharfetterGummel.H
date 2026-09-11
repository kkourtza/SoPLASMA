/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025
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
    Foam::fvm::ScharfetterGummel

Description

    Scharfetter-Gummel finite volume matrix operator for drift-diffusion.

    This function builds a coupled fvMatrix for the transport equation:
        div(J) = Source

    Where the total advective-diffusive flux J is defined as:
        J = phi * rho - D * grad(rho)

    This scheme solves the 1D convection-diffusion equation analytically 
    between cell centers (assuming constant E and D locally) to determine 
    the exponential profile of the scalar field. This ensures positivity 
    of the species density even in regions of high electric fields (sheaths) 
    where standard central or upwind schemes may fail or produce oscillations.

    The matrix coefficients are derived from the flux expression:
        J_f = D/h * [ B(-Pe) * rho_P  -  B(Pe) * rho_N ]

    Where:
      - Pe is the local Peclet number: phi / (D/h)
      - B(x) is the Bernoulli function: x / (exp(x) - 1)
      - rho_P, rho_N are the cell-center values of owner and neighbour.

    For a detailed theoretical derivation, please refer to:
    docs/theory/scharfetter_gummel.md

References
    Scharfetter, D. L., & Gummel, H. K. (1969).
    "Large-Signal Analysis of a Silicon Read Diode Oscillator."
    IEEE Transactions on Electron Devices, 16(1), 64–77.
    
\*---------------------------------------------------------------------------*/

#ifndef ScharfetterGummel_H
#define ScharfetterGummel_H

#include "fvMatrix.H"
#include "fvm.H"
#include "fvc.H"
#include "correctedSnGrad.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Function ScharfetterGummel
\*---------------------------------------------------------------------------*/

namespace fvm
{

    template<class Type>
    tmp<fvMatrix<Type>> ScharfetterGummel
    (
        const GeometricField<Type, fvPatchField, volMesh>& rho,
        const surfaceScalarField& phi, 
        const volScalarField& D
    )
    {
        const fvMesh& mesh = rho.mesh();

        tmp<fvMatrix<Type>> tfvm
        (
            new fvMatrix<Type>
            (
                rho,
                phi.dimensions() * rho.dimensions()
            )
        );
        fvMatrix<Type>& fvm = tfvm.ref();
    
        surfaceScalarField Df = fvc::interpolate(D);

        auto Bern = [](scalar x) -> scalar
        {
            const scalar ax = mag(x);
            if (ax < 1e-4)
            {
                return 1.0 - 0.5*x + (x*x)/12.0 - pow4(x)/720.0;;
            }
            if (x > 200.0)  return 0.0;
            if (x < -200.0) return -x;
            return x / (Foam::exp(x) - 1.0);
        };

        // NON-ORTHOGONAL DELTA COEFFICIENTS for the implicit part.
        //
        // Using mesh.deltaCoeffs() here leaves the scheme with NO
        // non-orthogonal treatment at all. MEASURED 2026-09-08 on
        // verification/fluxScheme2Dnonortho, an 11.3 deg sheared mesh: SG and
        // CFS sat on a FIXED error floor (~1.0e-2) that refinement did not
        // touch -- order 0.00 over four refinements -- while `standard` with a
        // corrected laplacian converged on the same mesh with the same
        // harness. They also "converged" in 2 non-orthogonal correctors
        // against standard's 15-16, because their matrix had no correction
        // term for the correctors to act on.
        //
        // THE SPLIT THAT MAKES THE FIX WELL-DEFINED. With the Bernoulli
        // identity B(-z) = B(z) + z the SG face flux factorises exactly:
        //
        //   Gamma = coeffP n_P - coeffN n_N
        //         = phi n_P - Df B(Pe) magSf snGrad(n)
        //
        // so SG's diffusive part is an ORDINARY diffusion flux with effective
        // diffusivity D*B(Pe). It therefore takes the standard non-orthogonal
        // treatment: nonOrthDeltaCoeffs in the implicit part, plus the
        // explicit snGrad correction scaled by Df*B(Pe).
        surfaceScalarField diffCond =
            Df * mesh.magSf() * mesh.nonOrthDeltaCoeffs();
        const scalarField& diffCondInternal = diffCond.primitiveField();
        const scalarField& convFlux = phi.primitiveField();

        const scalarField Pe = convFlux / (diffCondInternal + VSMALL);

        scalarField coeffP(diffCond.size());
        scalarField coeffN(diffCond.size());

        forAll(Pe, facei)
        {
            coeffP[facei] = diffCondInternal[facei] * Bern(-Pe[facei]);
            coeffN[facei] = diffCondInternal[facei] * Bern(Pe[facei]);
        }

        fvm.lower() = -coeffP; 
        fvm.upper() = -coeffN;
        fvm.negSumDiag();

// Boundaries 
forAll(rho.boundaryField(), patchi)
{
    const fvPatchField<Type>& ptf = rho.boundaryField()[patchi];
    const scalarField& patchPhi   = phi.boundaryField()[patchi];

    scalarField patchDiffCond =
          Df.boundaryField()[patchi]
        * mesh.magSf().boundaryField()[patchi]
        * mesh.nonOrthDeltaCoeffs().boundaryField()[patchi];

    scalarField patchPe = patchPhi / (patchDiffCond + VSMALL);

    scalarField pCoeffP(ptf.size()), pCoeffB(ptf.size());
    forAll(patchPe, facei)
    {
        pCoeffP[facei] = patchDiffCond[facei] * Bern(-patchPe[facei]);
        pCoeffB[facei] = patchDiffCond[facei] * Bern( patchPe[facei]);
    }

    if (ptf.coupled())
    {
        // Processor/cyclic: fvMatrix multiplies boundaryCoeffs by rho_N automatically
        fvm.internalCoeffs()[patchi] = pCoeffP;
        fvm.boundaryCoeffs()[patchi] = pCoeffB;
    }
    else
    {
        // Physical boundary: get alpha, beta from BC
        tmp<Field<Type>> tAlpha = ptf.valueInternalCoeffs(patchPhi);
        tmp<Field<Type>> tBeta  = ptf.valueBoundaryCoeffs(patchPhi);

        fvm.internalCoeffs()[patchi] = pCoeffP - pCoeffB * tAlpha();
        fvm.boundaryCoeffs()[patchi] =           pCoeffB * tBeta();
    }
}

        // ---- EXPLICIT NON-ORTHOGONAL CORRECTION -------------------------
        //
        // The diffusive flux above is -Df*B(Pe)*magSf*snGrad(n). Replacing
        // snGrad by its corrected form adds
        //
        //     dGamma = -Df*B(Pe)*magSf*snGradCorrection(n)
        //
        // and since the matrix represents div(Gamma), that explicit piece goes
        // to the source as  source -= V*div(dGamma).
        //
        // APPLIED UNCONDITIONALLY, which is safe and self-disabling: on an
        // orthogonal mesh the non-orthogonal correction vectors are ZERO, so
        // the term vanishes identically and the result is unchanged. The only
        // cost there is one gradient evaluation -- the same cost OpenFOAM's own
        // corrected laplacian already pays. There is deliberately NO user knob
        // and no corrector count to set: the correction is either needed by the
        // mesh or it is identically zero.
        {
            surfaceScalarField DeffB
            (
                IOobject
                (
                    "DeffB", mesh.time().timeName(), mesh,
                    IOobject::NO_READ, IOobject::NO_WRITE
                ),
                Df
            );

            scalarField& DeffBi = DeffB.primitiveFieldRef();
            forAll(DeffBi, facei)
            {
                DeffBi[facei] *= Bern(Pe[facei]);
            }

            forAll(DeffB.boundaryField(), patchi)
            {
                const scalarField& pPhi = phi.boundaryField()[patchi];
                const scalarField pdc
                (
                    Df.boundaryField()[patchi]
                  * mesh.magSf().boundaryField()[patchi]
                  * mesh.nonOrthDeltaCoeffs().boundaryField()[patchi]
                );
                scalarField& pD = DeffB.boundaryFieldRef()[patchi];
                forAll(pD, facei)
                {
                    pD[facei] *= Bern(pPhi[facei]/(pdc[facei] + VSMALL));
                }
            }

            const surfaceScalarField corrFlux
            (
              - DeffB*mesh.magSf()
               *fv::correctedSnGrad<Type>(mesh).correction(rho)()
            );

            fvm.source() -=
                mesh.V()*fvc::div(corrFlux)().primitiveField();
        }

        return tfvm;
    }

    // Overload for diffusivity as a dimensionedScalar
    template<class Type>
    tmp<fvMatrix<Type>> ScharfetterGummel
    (
        const GeometricField<Type, fvPatchField, volMesh>& rho,
        const surfaceScalarField& phi,
        const dimensionedScalar& D    
    )
    {
        volScalarField DField
        (
            IOobject
            (
                "D_tmp", 
                rho.mesh().time().timeName(), 
                rho.mesh(), 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ),
            rho.mesh(),
            D
        );

        // Call the core implementation
        return ScharfetterGummel(rho, phi, DField);
    }

    // Overload for Electric Potential and Mobility (mu) instead of flux
    template<class Type>
    tmp<fvMatrix<Type>> ScharfetterGummel
    (
        const GeometricField<Type, fvPatchField, volMesh>& rho,
        const volScalarField& ePotential,
        const volScalarField& mu,
        const volScalarField& D,
        const scalar chargeNumber
    )
    {
        const fvMesh& mesh = rho.mesh();
        
        surfaceScalarField phi = chargeNumber * fvc::interpolate(mu)
                                    * (-fvc::snGrad(ePotential))
                                    * mesh.magSf();   
  
        // Call the core implementation
        return ScharfetterGummel(rho, phi, D);
    }

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
