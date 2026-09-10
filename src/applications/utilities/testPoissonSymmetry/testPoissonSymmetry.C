/*---------------------------------------------------------------------------*\
    testPoissonSymmetry

    ISOLATES THE POISSON SOLVE from the plasma coupling, to answer one
    question: given a charge density that is EXACTLY mirror-symmetric, does
    the Poisson solve return a symmetric potential and a symmetric E?

    The whole difficulty in diagnosing the grubert asymmetry (2026-09-10) is
    that rho -> E -> species -> rho is a LOOP: every field in the coupled
    state is asymmetric, so no measurement of the running case can say which
    link introduces it. Stated by the user: "even with measuring E you cant
    really understand whats wrong cause the space charge asymmetry also
    reflects to E - theres a two way coupling".

    This breaks the loop. rho is symmetrised HERE, in memory, to bit equality;
    one Poisson solve follows, with the case's own fvSchemes and fvSolution;
    and the mirror asymmetry of phi and of E = -grad(phi) is reported. Any
    asymmetry in the output was manufactured by the solve.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // ---- mirror map about the y midplane
    const vectorField& C = mesh.C();
    scalar ymin = GREAT, ymax = -GREAT;
    forAll(C, c) { ymin = min(ymin, C[c].y()); ymax = max(ymax, C[c].y()); }
    const scalar W = ymin + ymax;

    labelList mirror(mesh.nCells(), -1);
    forAll(C, i)
    {
        scalar best = GREAT; label bj = -1;
        forAll(C, j)
        {
            const scalar d = magSqr(vector(C[i].x()-C[j].x(),
                                           (W-C[i].y())-C[j].y(),
                                           C[i].z()-C[j].z()));
            if (d < best) { best = d; bj = j; }
        }
        mirror[i] = bj;
    }
    Info<< "mirror map built for " << mesh.nCells() << " cells, midplane y = "
        << 0.5*W << endl;

    volScalarField ePotential
    (
        IOobject("ePotential", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    );

    volScalarField rho
    (
        IOobject("chargeDensity", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );

    // ---- SYMMETRISE rho to bit equality
    {
        scalarField& r = rho.primitiveFieldRef();
        scalarField r0(r);
        forAll(r, c) { r[c] = 0.5*(r0[c] + r0[mirror[c]]); }
        rho.correctBoundaryConditions();

        scalar worst = 0;
        forAll(r, c) { worst = max(worst, mag(r[c] - r[mirror[c]])); }
        Info<< "rho asymmetry AFTER symmetrising (must be 0): " << worst << endl;
    }

    const dimensionedScalar epsilon0("eps0", dimensionSet(-1,-3,4,0,0,2,0), 8.8541878128e-12);

    Info<< "solving laplacian(eps0, phi) == -rho, case's own schemes" << endl;
    for (int it = 0; it < 3; ++it)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(epsilon0, ePotential, "laplacian(epsilon,ePotential)")
         ==
           -rho
        );
        pEqn.solve();
    }

    const volScalarField& phi = ePotential;
    const volVectorField E("E", -fvc::grad(ePotential));

    scalar aphi = 0, aEx = 0, aEy = 0;
    scalar sphi = 0, sE = 0;
    forAll(phi, c)
    {
        sphi = max(sphi, mag(phi[c]));
        sE   = max(sE,   mag(E[c]));
    }
    forAll(phi, c)
    {
        const label m = mirror[c];
        aphi = max(aphi, mag(phi[c] - phi[m]));
        aEx  = max(aEx,  mag(E[c].x() - E[m].x()));
        aEy  = max(aEy,  mag(E[c].y() + E[m].y()));   // Ey is ODD
    }

    Info<< nl << "======== RESULT: symmetric rho in, what comes out ========" << nl
        << "  |phi| scale                 " << sphi << nl
        << "  |E|   scale                 " << sE << nl
        << "  phi mirror asymmetry (rel)  " << aphi/max(sphi,SMALL) << nl
        << "  Ex  mirror asymmetry (rel)  " << aEx/max(sE,SMALL) << nl
        << "  Ey  anti-symmetry violation " << aEy/max(sE,SMALL) << nl
        << nl
        << "  ~1e-16 : the Poisson solve is symmetric; the asymmetry is NOT here." << nl
        << "  >>1e-16: the Poisson solve itself breaks symmetry." << nl
        << endl;

    return 0;
}
