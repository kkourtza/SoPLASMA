/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.
    Copyright (C) 2026
    GNU General Public License v3 or later.

Application
    testFluxScheme

Description
    Verification bed for the INTERIOR flux schemes -- `standard`
    (fvm::div + fvm::laplacian) and `ScharfetterGummel` -- against an EXACT
    solution of the steady 1D drift-diffusion equation WITH A SOURCE:

        d/dx ( v n - D dn/dx ) = S,     n(0) = a,  n(L) = b

    whose solution is, with the Peclet number Pe = v L / D,

        n(x) = C1 + C2 exp(Pe x/L) + (S/v) x
        C2   = (b - a - S L / v) / (exp(Pe) - 1),    C1 = a - C2.

    WHY A NON-ZERO SOURCE.  This is the whole point of the bed, not a
    decoration.  Scharfetter-Gummel derives its face flux from the same local
    two-point problem WITH s SET TO ZERO; the Complete Flux Scheme (Liu et al
    2014, PSST 23 015023; ten Thije Boonkkamp & Anthonissen 2011) keeps s and
    gains an inhomogeneous, source-carrying part. A source-free test cannot
    tell the two apart, so it would be worthless as the acceptance test for
    CFS. With S != 0 the exact solution above separates them.

    It exercises the SHIPPED operators on a real mesh. A reimplementation of
    the scheme formulas here would verify nothing.

    Peclet is swept by the CALLER (Allrun) through -v and -D; the mesh is
    swept by regenerating blockMesh. Reports L2 and Linf against the exact
    solution, so the caller can fit a convergence order.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ScharfetterGummel.H"
#include "CompleteFlux.H"
#include "emptyFvPatch.H"
#include <cmath>
#include <cstdio>

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addOption("v", "scalar", "drift velocity [m/s] (non-zero)");
    argList::addOption("D", "scalar", "diffusivity [m2/s]");
    argList::addOption("S", "scalar", "source amplitude [1/s]");
    argList::addOption("mms", "scalar",
        "1D manufactured-solution amplitude A; 0 (default) = uniform source");
    argList::addOption("mms2d", "scalar",
        "2D manufactured-solution amplitude; imposes the exact value on EVERY"
        " patch, so the bed is valid on ANY mesh (sheared, graded, skewed)");
    argList::addOption("a", "scalar", "n at the left boundary");
    argList::addOption("b", "scalar", "n at the right boundary");
    argList::addOption("scheme", "word",
        "standard | ScharfetterGummel | CompleteFlux");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const scalar v  = args.getOrDefault<scalar>("v", 1.0);
    const scalar Dv = args.getOrDefault<scalar>("D", 1.0);
    const scalar S  = args.getOrDefault<scalar>("S", 1.0);

    // MANUFACTURED SOLUTION, and it exists because a UNIFORM source cannot
    // test CFS at all: Gamma^i = h(gamma s_P + delta s_N) is then IDENTICAL on
    // every interior face, so div(Gamma^i) == 0 and CFS collapses onto SG in
    // the interior. Measured 2026-09-07: the two came out 0.2% apart, and that
    // residue is boundary-only.
    //
    // With `-mms A` the solution is chosen and the source DERIVED from it:
    //     n(x) = a + (b-a) x + A sin(pi x)        (n(0)=a, n(1)=b exactly)
    //     s(x) = v n' - D n''
    //          = v[(b-a) + A pi cos(pi x)] + D A pi^2 sin(pi x)
    // so s VARIES in space, div(Gamma^i) != 0, and both the drift and the
    // diffusion terms are exercised at comparable magnitude -- DeChant 2023's
    // stated requirement for a useful manufactured solution.
    const scalar A = args.getOrDefault<scalar>("mms", 0.0);
    const bool mms = (mag(A) > VSMALL);

    // 2D MANUFACTURED SOLUTION -- the only form valid on a NON-ORTHOGONAL mesh.
    //
    //   n(x,y) = 1 + x + A2 sin(pi x) sin(pi y)
    //   s      = v dn/dx - D lap(n)
    //          = v[1 + A2 pi cos(pi x) sin(pi y)]
    //          + 2 D A2 pi^2 sin(pi x) sin(pi y)
    //
    // and the EXACT VALUE IS IMPOSED ON EVERY BOUNDARY FACE, evaluated at that
    // face's own centre. That is what makes it mesh-agnostic. The earlier
    // sheared bed failed precisely here: it kept the 1D solution n(x) while
    // shearing turned the left/right patches into SLANTED planes, so "n = a at
    // the left boundary" stopped being a condition at constant x. The boundary
    // data and the geometry then disagreed by O(shear), and the measured error
    // floor tracked the shear ~1:1 (0.01 -> 8.0e-3, 0.10 -> 1.07e-1) for BOTH
    // schemes at EVERY Peclet, with zero convergence. That was the test, not
    // the schemes.
    const scalar A2 = args.getOrDefault<scalar>("mms2d", 0.0);
    const bool mms2 = (mag(A2) > VSMALL);
    const scalar a  = args.getOrDefault<scalar>("a", 1.0);
    const scalar b  = args.getOrDefault<scalar>("b", 1.0);
    const word scheme
    (
        args.getOrDefault<word>("scheme", "standard")
    );

    if (mag(v) < VSMALL)
    {
        FatalErrorInFunction
            << "v must be non-zero: the exact solution divides by it." << nl
            << exit(FatalError);
    }

    const scalar Lx = 1.0;   // blockMeshDict fixes the domain at [0,1]

    volScalarField n
    (
        IOobject("n", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    );

    volScalarField D
    (
        IOobject("D", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("D", dimensionSet(0,2,-1,0,0,0,0), Dv)
    );

    volScalarField Ssrc
    (
        IOobject("S", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("S", dimensionSet(0,0,-1,0,0,0,0), S)
    );

    if (mms2)
    {
        const scalar pi = constant::mathematical::pi;

        auto nex = [&](const scalar x, const scalar y)
        {
            return 1.0 + x + A2*std::sin(pi*x)*std::sin(pi*y);
        };
        auto sex = [&](const scalar x, const scalar y)
        {
            return v*(1.0 + A2*pi*std::cos(pi*x)*std::sin(pi*y))
                 + 2.0*Dv*A2*pi*pi*std::sin(pi*x)*std::sin(pi*y);
        };

        const volVectorField& CC = mesh.C();
        scalarField& si = Ssrc.primitiveFieldRef();
        forAll(si, c) si[c] = sex(CC[c].x(), CC[c].y());

        forAll(Ssrc.boundaryFieldRef(), pI)
        {
            const fvPatch& p = mesh.boundary()[pI];
            scalarField& sb = Ssrc.boundaryFieldRef()[pI];
            forAll(sb, f) sb[f] = sex(p.Cf()[f].x(), p.Cf()[f].y());
        }

        // IMPOSE THE EXACT SOLUTION ON EVERY PATCH, at each face's own centre.
        forAll(n.boundaryFieldRef(), pI)
        {
            const fvPatch& p = mesh.boundary()[pI];
            if (p.coupled() || isA<emptyFvPatch>(p)) continue;
            scalarField nb(p.size());
            forAll(nb, f) nb[f] = nex(p.Cf()[f].x(), p.Cf()[f].y());
            n.boundaryFieldRef()[pI] == nb;
        }
    }
    else if (mms)
    {
        const scalar pi = constant::mathematical::pi;
        const volScalarField Cxs(mesh.C().component(0));
        const scalarField& xs = Cxs.primitiveField();
        scalarField& si = Ssrc.primitiveFieldRef();

        forAll(si, c)
        {
            si[c] = v*((b - a) + A*pi*std::cos(pi*xs[c]))
                  + Dv*A*pi*pi*std::sin(pi*xs[c]);
        }
        Ssrc.correctBoundaryConditions();

        forAll(Ssrc.boundaryFieldRef(), pi_)
        {
            const fvPatch& p = mesh.boundary()[pi_];
            scalarField& sb = Ssrc.boundaryFieldRef()[pi_];
            const vectorField& Cf = p.Cf();
            forAll(sb, f)
            {
                sb[f] = v*((b - a) + A*pi*std::cos(pi*Cf[f].x()))
                      + Dv*A*pi*pi*std::sin(pi*Cf[f].x());
            }
        }
    }

    // Uniform drift in +x: phi = v * Sf_x
    surfaceScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        v*(mesh.Sf().component(0))
       *dimensionedScalar("u", dimensionSet(0,1,-1,0,0,0,0), 1.0)
       /dimensionedScalar("one", dimless, 1.0)
    );

    // NON-ORTHOGONAL CORRECTORS -- the equation is RE-ASSEMBLED each pass.
    //
    // OpenFOAM's `corrected` laplacian puts the non-orthogonal part in the
    // matrix SOURCE at ASSEMBLY time, computed from the current field. Simply
    // calling solve() again re-solves the SAME matrix and changes nothing --
    // measured: identical digits. The equation has to be rebuilt so the
    // correction is recomputed against the updated n.
    //
    // MEASURED 2026-09-08 on the sheared bed with a single assembly: BOTH SG
    // and CFS stalled at ~1.00e-2 with ZERO convergence at shear = 0.2, while
    // shear = 0 converged at order 2 -- an orthogonal mesh has no correction to
    // make, which is what localised the fault to the harness rather than the
    // schemes.
    int nCorr = 0;

    for (int nonOrth = 0; nonOrth < 50; ++nonOrth)
    {
        ++nCorr;

        tmp<fvScalarMatrix> tEqn;

        if (scheme == "ScharfetterGummel")
        {
            tEqn = (fvm::ScharfetterGummel(n, phi, D) == Ssrc);
        }
        else if (scheme == "CompleteFlux")
        {
            tEqn = (fvm::CompleteFlux(n, phi, D, Ssrc) == Ssrc);
        }
        else if (scheme == "standard")
        {
            tEqn = (fvm::div(phi, n, "div(phi,n)")
                  - fvm::laplacian(D, n, "laplacian(D,n)") == Ssrc);
        }
        else
        {
            FatalErrorInFunction
                << "unknown scheme '" << scheme << "'." << nl
                << "Valid: (standard | ScharfetterGummel | CompleteFlux)" << nl
                << exit(FatalError);
        }

        const scalar res = tEqn.ref().solve().initialResidual();

        if (res < 1e-13) break;
    }

    // ---- exact solution -------------------------------------------------
    const scalar Pe = v*Lx/Dv;

    // NUMERICALLY STABLE FORM.  The textbook expression multiplies
    // C2 = K/(exp(Pe)-1) by exp(Pe x), and at Pe = 100 that is 1e-44 times
    // 1e43 -- representable, but it throws away most of the mantissa. Written
    // as exp(Pe(x-1))/(1-exp(-Pe)) both factors are O(1) for every Pe > 0.
    const scalar K = b - a - S*Lx/v;

    // HOLD THE tmp.  `mesh.C().component(0)()` returns a reference INTO a
    // temporary that dies at the end of the statement; binding a
    // const scalarField& to it dangles, and the resulting garbage x produced
    // L2 = 1e53 and SIGFPEs in the exact solution -- with the linear solve
    // converging to 7e-16 the whole time. Caught 2026-09-07 by this bed's own
    // L2rel column reading exactly 1.0.
    const volScalarField Cx(mesh.C().component(0));
    const scalarField& xc = Cx.primitiveField();
    const volScalarField Cy(mesh.C().component(1));
    const scalarField& ycc = Cy.primitiveField();

    scalar l2 = 0.0, linf = 0.0, nrm = 0.0;

    forAll(xc, c)
    {
        const scalar xi = xc[c]/Lx;
        const scalar ex =
            mms2
          ? 1.0 + xc[c]
              + A2*std::sin(constant::mathematical::pi*xc[c])
                  *std::sin(constant::mathematical::pi*ycc[c])
          : mms
          ? a + (b - a)*xi
              + A*std::sin(constant::mathematical::pi*xi)
          : a + K*(std::exp(Pe*(xi - 1.0)) - std::exp(-Pe))
                 /(1.0 - std::exp(-Pe))
              + (S/v)*xc[c];
        const scalar e  = n[c] - ex;
        l2   += e*e;
        nrm  += ex*ex;
        linf  = max(linf, mag(e));
    }

    l2  = std::sqrt(l2/max(scalar(xc.size()), scalar(1)));
    nrm = std::sqrt(nrm/max(scalar(xc.size()), scalar(1)));

    const scalar h        = Lx/scalar(mesh.nCells());
    const scalar PeGrid   = v*h/Dv;

    std::printf
    (
        "RESULT scheme=%-18s N=%5d h=%.6e Pe=%.6e PeGrid=%.6e "
        "L2=%.8e Linf=%.8e L2rel=%.8e nCorr=%d\n",
        scheme.c_str(), mesh.nCells(), h, Pe, PeGrid, l2, linf,
        l2/max(nrm, VSMALL), nCorr
    );

    Info<< "End\n" << endl;
    return 0;
}
