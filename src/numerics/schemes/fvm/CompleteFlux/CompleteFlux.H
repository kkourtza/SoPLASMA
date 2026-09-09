/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.
    Copyright (C) 2026
    GNU General Public License v3 or later.

Function
    Foam::fvm::CompleteFlux

Description
    The COMPLETE FLUX SCHEME (CFS) for the drift-diffusion-reaction equation,
    after ten Thije Boonkkamp & Anthonissen (2011) J. Sci. Comput. 46:47-70 and
    Liu, van Dijk, ten Thije Boonkkamp, Mihailova & van der Mullen (2014),
    Plasma Sources Sci. Technol. 23:015023 (Literature/).

    Scharfetter-Gummel derives its face flux from the local two-point problem

        d/dx ( v n - D dn/dx ) = s,   x_j < x < x_j+1

    WITH s SET TO ZERO. CFS keeps s. The flux therefore splits as Liu eq. (18):

        Gamma_j+1/2 = Gamma^h_j+1/2 + Gamma^i_j+1/2                     (18a)
        Gamma^h     = alpha n_j - beta n_j+1                            (18b)
        Gamma^i     = h ( gamma s_j + delta s_j+1 )                     (18c)

    with, from eqs (19a)-(19d), (22), (23),

        alpha = (D/h) B(-Pe),      beta  = (D/h) B(Pe)
        gamma = max(1/2 - W(Pe), 0),   delta = min(1/2 - W(Pe), 0)
        B(z)  = z/(exp(z) - 1)
        W(z)  = (exp(z) - 1 - z)/(z(exp(z) - 1)) = (1 - B(z))/z

    Gamma^h IS EXACTLY THE SG FLUX, so this operator is assembled as the
    shipped fvm::ScharfetterGummel matrix plus the divergence of Gamma^i. And
    since Gamma^i depends on s ONLY -- not on n -- it is an EXPLICIT
    contribution, which is also how Liu treats it ("the source terms in these
    DDR equations are treated explicitly").

    WHY THIS MATTERS, measured in verification/fluxScheme1D on 2026-09-07:
    with s = 0, SG is exact to machine precision (5e-16). With s != 0 its error
    is a UNIFORM offset equal to the source flux it drops over half a cell --
    S h^2/(8D) at low grid Peclet and S h/(2v) at high. That linear-in-h
    scaling is why SG degrades to exactly first order for |Pe_grid| >> 1
    (measured 1.00 over five refinements at domain Pe = 1e4). CFS exists to
    remove precisely that term.

    NOTE ON gamma/delta: only ONE of them is non-zero at a time -- they are the
    positive and negative parts of the same quantity (1/2 - W). The scheme
    therefore takes the source from the upwind side when 1/2 - W > 0 and from
    the downwind side otherwise.

SourceFiles
    CompleteFlux.H (header only, templated)

\*---------------------------------------------------------------------------*/

#ifndef CompleteFlux_H
#define CompleteFlux_H

#include "ScharfetterGummel.H"
#include "fvMatrix.H"
#include "fvm.H"
#include "fvc.H"
#include "surfaceFields.H"
#include "volFields.H"

namespace Foam
{

namespace fvc
{

//- The INHOMOGENEOUS part of the complete flux, Liu eq. (18c), returned as an
//  OpenFOAM face flux (i.e. already multiplied by |Sf|).
inline tmp<surfaceScalarField> completeFluxInhomogeneous
(
    const surfaceScalarField& phi,   // convective face flux, v*Sf
    const volScalarField& D,         // diffusivity
    const volScalarField& s          // NET volumetric source (production-loss)
)
{
    const fvMesh& mesh = D.mesh();

    // W(z) = (1 - B(z))/z, with the series used near z = 0 where that form
    // is 0/0. B(z) = z/(exp(z)-1); W(0) = 1/2, so (1/2 - W) -> 0 there.
    auto Wfun = [](const scalar z) -> scalar
    {
        const scalar az = mag(z);

        if (az < 1e-4)
        {
            // W = 1/2 - z/12 + z^3/720 - ...
            return 0.5 - z/12.0 + z*z*z/720.0;
        }
        if (z >  200.0) return 1.0/z;          // B -> 0
        if (z < -200.0) return (1.0 + z)/z;    // B -> -z

        const scalar B = z/(Foam::exp(z) - 1.0);
        return (1.0 - B)/z;
    };

    // V(z) -- the SECOND-MOMENT weight, needed for a source that VARIES across
    // the interval. Liu eq. (18c) is a two-point formula in s that is EXACT
    // only for a CONSTANT source; with s linear across [x_j, x_j+1] the exact
    // inhomogeneous flux is
    //
    //   Gamma^i = h [ s_L (1/2 - W(Pe)) + (s_R - s_L)(1/8 - V(Pe)/2) ]
    //   V(z)    = [2 - exp(-z)(2 + 2z + z^2)] / [z^2 (1 - exp(-z))]
    //
    // derived from Gamma^i = S(h/2) - [lambda/(1-e^-Pe)] Int_0^h S(x)e^{-lambda x}dx
    // with S(x) = Int_0^x s. At large Pe, V -> 0, so the exact flux tends to
    // h[s_L/2 + (s_R-s_L)/8] while (18c) gives only h s_L/2 -- it DROPS
    // h(s_R-s_L)/8 ~ h^2 s', which is an O(h) error in the divergence and caps
    // the scheme at FIRST order exactly where SG already is. MEASURED
    // 2026-09-07 before this term was added: CFS tracked SG at p = 0.99-1.00
    // for Pe_grid = 31..500, gaining only ~15% in absolute error.
    //
    // Checked at Pe -> 0 against pure diffusion (-D n'' = s), where the exact
    // inhomogeneous flux is -h(s_R - s_L)/24; this expression gives
    // 1/8 - V(0)/2 = 1/8 - 1/6 = -1/24. Agrees.
    auto Vfun = [](const scalar z) -> scalar
    {
        const scalar az = mag(z);

        if (az < 1e-2)
        {
            // V = 1/3 - z/12 + z^2/360 - ...  (the direct form loses all
            // significance here: the numerator is O(z^3) against a leading 2)
            return 1.0/3.0 - z/12.0 + z*z/360.0;
        }
        if (z >  200.0) return 2.0/(z*z);   // exp(-z) -> 0
        if (z < -200.0) return 1.0;         // exp(-z) dominates both parts

        const scalar e = Foam::exp(-z);
        return (2.0 - e*(2.0 + 2.0*z + z*z))/(z*z*(1.0 - e));
    };

    tmp<surfaceScalarField> tGi
    (
        new surfaceScalarField
        (
            IOobject
            (
                "GammaI(" + s.name() + ")",
                mesh.time().timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(phi.dimensions()*s.dimensions()*dimTime, Zero)
        )
    );
    surfaceScalarField& Gi = tGi.ref();

    const surfaceScalarField Df(fvc::interpolate(D));
    const surfaceScalarField diffCond(Df*mesh.magSf()*mesh.deltaCoeffs());

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    const scalarField& sI  = s.primitiveField();
    const scalarField& phiI = phi.primitiveField();
    const scalarField& dcI = diffCond.primitiveField();
    const scalarField& dCoI = mesh.deltaCoeffs().primitiveField();
    const scalarField& magSfI = mesh.magSf().primitiveField();

    scalarField& GiI = Gi.primitiveFieldRef();

    forAll(GiI, f)
    {
        const scalar Pe = phiI[f]/(dcI[f] + VSMALL);
        const scalar h  = 1.0/(dCoI[f] + VSMALL);

        const scalar sL = sI[own[f]];
        const scalar sR = sI[nei[f]];

        // EXACT for a source varying linearly across the face -- Liu's (18c)
        // first term, plus the gradient term it omits. See Vfun above.
        GiI[f] = h*( sL*(0.5 - Wfun(Pe))
                   + (sR - sL)*(0.125 - 0.5*Vfun(Pe)) )*magSfI[f];
    }

    // Boundary faces: no neighbour cell, so the downwind source is the patch
    // value and h is the cell-to-face distance.
    forAll(Gi.boundaryField(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        if (p.coupled() || Gi.boundaryField()[patchi].size() == 0) continue;

        const labelUList& fc = p.faceCells();
        const scalarField& phiB = phi.boundaryField()[patchi];
        const scalarField& sB   = s.boundaryField()[patchi];
        const scalarField dcB
        (
            Df.boundaryField()[patchi]*p.magSf()*p.deltaCoeffs()
        );

        scalarField& GiB = Gi.boundaryFieldRef()[patchi];

        forAll(GiB, f)
        {
            const scalar Pe = phiB[f]/(dcB[f] + VSMALL);
            const scalar h  = 1.0/(p.deltaCoeffs()[f] + VSMALL);

            const scalar sC = sI[fc[f]];   // at x = 0 of the local interval
            const scalar sF = sB[f];       // at x = H, the boundary face

            // THE BOUNDARY COEFFICIENTS ARE NOT THE INTERIOR ONES.
            //
            // For an INTERIOR face the local two-point problem spans the two
            // cell centres and the face sits at its MIDPOINT, so the flux is
            // evaluated at x = h/2 and S(h/2) = s h/2 supplies the 1/2 in
            // (1/2 - W). For a BOUNDARY face the interval runs from the cell
            // centre to the face and the flux is evaluated at its END, x = H,
            // where S(H) = s H -- giving 1, not 1/2:
            //
            //   Gamma^i_b = H [ s_C (1 - W(P)) + (s_F - s_C)(1/2 - V(P)/2) ]
            //
            // Using the interior coefficients here halves the boundary term.
            // MEASURED 2026-09-07 on the uniform-source case at Pe = 1e4,
            // where Gamma^i is exact and SG's error is analytically S h/(2v):
            // CFS removed EXACTLY half of it (ratio 0.536, 0.518, 0.509,
            // 0.505 -> 1/2 as h -> 0). In a drift-dominated problem the
            // inflow boundary sets the error, so this term is not a detail.
            GiB[f] = h*( sC*(1.0 - Wfun(Pe))
                       + (sF - sC)*(0.5 - 0.5*Vfun(Pe)) )*p.magSf()[f];
        }
    }

    return tGi;
}

} // End namespace fvc


namespace fvm
{

//- CFS matrix: the SG (homogeneous) matrix plus div of the explicit
//  inhomogeneous flux. Equation convention matches driftDiffusion::nEqn(),
//  where the returned matrix is ADDED and represents div - laplacian.
inline tmp<fvScalarMatrix> CompleteFlux
(
    const volScalarField& n,
    const surfaceScalarField& phi,
    const volScalarField& D,
    const volScalarField& s
)
{
    tmp<fvScalarMatrix> tm(fvm::ScharfetterGummel(n, phi, D));

    // += a volScalarField moves it to the LHS (source() -= V*field), which is
    // what `ddt(n) + div(Gamma^h) + div(Gamma^i) = s` requires.
    tm.ref() += fvc::div(fvc::completeFluxInhomogeneous(phi, D, s)());

    return tm;
}

} // End namespace fvm

} // End namespace Foam

#endif
