/*---------------------------------------------------------------------------*\
  File: testDischargeCurrent.C
  Part of: SoPLASMA

Application
    testDischargeCurrent

Description
    Validates Sato's geometric weighting field against an ANALYTIC gap
    capacitance. No plasma is involved, which is the point: it exercises the
    weighting-field solve -- including the MONOLITHIC multi-region assembly --
    on a geometry whose answer is known in closed form.

    For a parallel plate of area A made of slabs (d_i, eps_r,i) in series,

        C = eps0 A / SUM_i ( d_i / eps_r,i )

    and Sato's C_g must reproduce it. The identity is exact, not approximate:
    with D = eps_i E_i continuous across the interfaces and SUM_i E_i d_i = 1
    for the unit-voltage weighting field,

        INT eps |e_hat|^2 dV = A * D * SUM_i E_i d_i = A * D = C

    So a discrepancy is a defect in the solve, the interface coupling or the
    volume scaling -- not modelling error. That is what makes it a decisive
    test rather than a plausibility check.

    Usage:
        testDischargeCurrent -case <caseDir> [-region <gasRegion>] \
            -expected <C_analytic>

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "regionProperties.H"
#include "dynamicFvMesh.H"
#include "electromagneticsModel.H"
#include "plasmaDischargeCurrent.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "expected",
        "F",
        "analytic gap capacitance to compare against [F]"
    );
    argList::addOption
    (
        "tolerance",
        "rel",
        "relative tolerance on the comparison (default 1e-3)"
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"

    // The regions are built exactly as the solvers build them, so the model
    // sees the same topology -- and so the MONOLITHIC assembly is exercised
    // rather than a single-region shortcut.
    UPtrList<fvMesh> dielectricMeshes(dielectricRegions.size());
    forAll(dielectricRegions, i)
    {
        dielectricMeshes.set(i, &dielectricRegions[i]);
    }

    fvMesh& mesh = gasMesh();

    Info<< nl << "testDischargeCurrent: Sato weighting field vs analytic C"
        << nl << "==========================================================="
        << nl << endl;

    autoPtr<electromagneticsModel> em =
        electromagneticsModel::New(gasMesh(), dielectricMeshes);

    IOdictionary plasmaDict
    (
        IOobject
        (
            "plasmaSimulationControls",
            runTime.system(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    plasmaDischargeCurrent current(mesh, plasmaDict, em());

    if (!current.enabled())
    {
        FatalErrorInFunction
            << "dischargeCurrent is not enabled in"
            << " system/plasmaSimulationControls." << nl
            << "    Add a `dischargeCurrent` sub-dictionary with"
            << " `enabled true;` and the electrode patch names." << nl
            << exit(FatalError);
    }

    current.writeFields();

    const scalar Cg = current.Cg();

    Info<< nl << "  C_g (Sato)     " << Cg << " F" << endl;

    int nFail = 0;

    if (args.found("expected"))
    {
        const scalar Cref = args.get<scalar>("expected");
        const scalar tol  = args.getOrDefault<scalar>("tolerance", 1e-3);
        const scalar rel  = mag(Cg - Cref)/max(mag(Cref), VSMALL);

        Info<< "  C analytic     " << Cref << " F" << nl
            << "  relative error " << rel
            << "   (tolerance " << tol << ")" << nl << endl;

        if (rel <= tol)
        {
            Info<< "  [ok  ] weighting field reproduces the analytic"
                << " capacitance" << endl;
        }
        else
        {
            Info<< "  [FAIL] C_g does not match the analytic value." << nl
                << "         Check, in order: that every electrode patch is"
                << " named in dischargeCurrent; that the interface patches"
                << " carry useImplicit true so the assembly is monolithic;"
                << " and the wedge revolution factor on an axisymmetric mesh."
                << endl;
            ++nFail;
        }
    }
    else
    {
        Info<< "  (no -expected given: reporting only)" << endl;
    }

    Info<< nl << (nFail ? "FAILED" : "PASSED") << nl << endl;

    return nFail ? 1 : 0;
}


// ************************************************************************* //
