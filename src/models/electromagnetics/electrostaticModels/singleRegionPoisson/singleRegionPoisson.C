/*---------------------------------------------------------------------------*\
  File: singleRegionPoisson.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::singleRegionPoisson.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "singleRegionPoisson.H"
#include "plasmaSimulationProfiler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(singleRegionPoisson, 0);
addToRunTimeSelectionTable
(
    electromagneticsModel,
    singleRegionPoisson,
    dictionary
);

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void singleRegionPoisson::updateDerivedFields()
{
    plasmaSimulationProfiler::start("Electromagnetics", "Calc phiE");
    phiE_ = -fvc::snGrad(ePotential_) * mesh_.magSf();
    plasmaSimulationProfiler::stop("Electromagnetics", "Calc phiE");

    plasmaSimulationProfiler::start("Electromagnetics", "Calc E");
    if (EScheme_ == "reconstruct")
    {
        E_ = fvc::reconstruct(phiE_);
    }
    else // "grad"
    {
        E_ = -fvc::grad(ePotential_);
    }
    plasmaSimulationProfiler::stop("Electromagnetics", "Calc E");

    plasmaSimulationProfiler::start("Electromagnetics", "Correct E boundary conditions");
    E_.correctBoundaryConditions();
    plasmaSimulationProfiler::stop("Electromagnetics", "Correct E boundary conditions");


    plasmaSimulationProfiler::start("Electromagnetics", "Calc Emag");
    Emag_ = mag(E_);
    Emag_.correctBoundaryConditions();
    plasmaSimulationProfiler::stop("Electromagnetics", "Calc Emag");


    plasmaSimulationProfiler::start("Electromagnetics", "Calc reducedE");

    // REGRESSION, measured 2026-09-03. Before fe827ae (2026-08-11, "One owner
    // for the gas density") backgroundDensityUniform_ defaulted to 2.5e25 via
    // `coeffs.getOrDefault<scalar>("backgroundDensity", 2.5e25)`. That commit
    // removed the default and made the key fatal -- correctly, since a silent
    // 2.5e25 in a case that never stated one is invented physics -- but left
    // THIS division unguarded. The only caller of setBackgroundDensity() is
    // plasmaSpecies, which does not exist in singleRegionElectrostaticFoam, so
    // every pure-electrostatics case has died with SIGFPE in Emag_/0 since.
    // Verified 2026-09-03 on the SHIPPED, unmodified tutorial
    // tutorials/electrostatics/singleRegionElectrostaticFoam/plate2D_timeVaryingBC:
    // the Poisson solve converged (residual 7.98e-10) and the crash was here.
    //
    // multiRegionPoisson::updateDerivedFields() was guarded at the time; this
    // is the same guard, so the two models cannot diverge again.
    //
    // Leaving reducedE at ZERO is the honest answer: |E|/N has no meaning where
    // there is no gas, and a fabricated N would propagate into every rate
    // coefficient and mobility. Said once, because a field that silently stays
    // zero is its own trap.
    if (backgroundDensityFieldPtr_)
    {
        reducedE_ = Emag_ / *backgroundDensityFieldPtr_;
        reducedE_.correctBoundaryConditions();
    }
    else if (backgroundDensityUniform_.value() > SMALL)
    {
        reducedE_ = Emag_ / backgroundDensityUniform_;
        reducedE_.correctBoundaryConditions();
    }
    else if (!noGasReported_)
    {
        noGasReported_ = true;

        Info<< "singleRegionPoisson: no background gas density, so `reducedE`"
               " stays ZERO." << nl
            << "    |E|/N is undefined without a gas. This is expected for a"
               " pure electrostatics" << nl
            << "    case; if this run has a plasma, its species model has not"
               " published a density" << nl
            << "    and every rate coefficient keyed on reducedE would be"
               " wrong." << endl;
    }
    plasmaSimulationProfiler::stop("Electromagnetics", "Calc reducedE");

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

singleRegionPoisson::singleRegionPoisson
(
    const fvMesh& mesh,
    const UPtrList<fvMesh>& dielectricMeshes
)
:
    electromagneticsModel(mesh, dielectricMeshes),
    EScheme_("reconstruct"),
    PoissonScheme_("explicit"),
    nNonOrthCorr_(0),
    backgroundDensityUniform_
    (
        "backgroundDensity",
        dimensionSet(0, -3, 0, 0, 0, 0, 0),
        0.0
    ),
    backgroundDensityFieldPtr_(nullptr)
{
    if (dielectricMeshes.size() > 0)
    {
        FatalErrorInFunction
            << "Dielectric regions were detected (" << dielectricMeshes.size()
            << "), but you are using the 'singleRegionPoisson' model." << nl
            << "This model ignores dielectrics and solves only in the gas."
            << nl << "Switch to 'multiRegionPoisson' model in your "
            << "properties file or remove the dielectric regions." << nl
            << exit(FatalError);
    }

    // The `<type>Coeffs` sub-dictionary of constant/electromagneticsProperties
    // is DEPRECATED and therefore optional. The numerics now come from
    // system/plasmaSimulationControls (`poisson`) and the permittivity from
    // constant/<region>/electricalProperties; this dictionary is consulted only
    // as a fallback so an unmigrated case keeps running, and doing so prints
    // where each setting has moved to.
    const dictionary& coeffs(subOrEmptyDict(type() + "Coeffs"));

    epsilonR_ = readEpsilonR(mesh_, coeffs);
    epsilon_ = dimensionedScalar
                            ("epsilon", epsilonR_ * constant::plasma::epsilon0);

    const poissonNumerics num(readPoissonNumerics(mesh_, coeffs));

    EScheme_ = num.EScheme;
    PoissonScheme_ = num.scheme;
    nNonOrthCorr_ = num.nNonOrthogonalCorrectors;

    // NOT read from this dictionary. The gas density has one owner --
    // `backgroundGas` in plasmaSpeciesProperties, where it is either stated or
    // closed from pressure and temperature -- and plasmaSpecies pushes it in
    // through setBackgroundDensity(). A second copy here silently evaluated
    // reducedE, and therefore every rate coefficient and mobility in the run,
    // at a different gas state than the case declared.
    if (coeffs.found("backgroundDensity"))
    {
        FatalIOErrorInFunction(coeffs)
            << "`backgroundDensity` in " << coeffs.name() << " is no longer"
            << " read: it was a second, independent definition of the gas"
            << " density and it did not agree with the one in"
            << " plasmaSpeciesProperties." << nl
            << "    Remove it. The density comes from `backgroundGas` there,"
            << " as `numberDensity` or as `pressure` + `T`." << nl
            << exit(FatalIOError);
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

//- Explicit Poisson branch
void singleRegionPoisson::solve()
{
    for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        plasmaSimulationProfiler::start("Electromagnetics", "Build ePotentialEqn");
        fvScalarMatrix ePotentialEqn
        (
            fvm::laplacian(epsilon_, ePotential_)
         == -chargeDensity_
        );
        plasmaSimulationProfiler::stop("Electromagnetics", "Build ePotentialEqn");

        if (nonOrth < nNonOrthCorr_)
        {
            ePotentialEqn.relax();
        }
        plasmaSimulationProfiler::start("Electromagnetics", "Solve ePotentialEqn");
        ePotentialEqn.solve();
        plasmaSimulationProfiler::stop("Electromagnetics", "Solve ePotentialEqn");
    }

    // plasmaSimulationProfiler::start("emupdateDerivedFields");
    updateDerivedFields();
    // plasmaSimulationProfiler::stop("emupdateDerivedFields");
}

//- Semi-implicit Poisson branch
void singleRegionPoisson::solve
(
    const volScalarField& electricalConductivity,
    const volScalarField& diffusiveChargeSource
)
{
    if (PoissonScheme_ == "explicit")
    {
        this->solve();
        return;
    }

    const dimensionedScalar deltaT = mesh_.time().deltaT();

    const volScalarField effEps(epsilon_ + deltaT * electricalConductivity);
    const volScalarField rhsSource(-chargeDensity_ - deltaT * diffusiveChargeSource);

    for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        fvScalarMatrix ePotentialEqn
        (
            fvm::laplacian
            (
                effEps, ePotential_
            )
         == rhsSource
        );

        if (nonOrth < nNonOrthCorr_)
        {
            ePotentialEqn.relax();
        }

        ePotentialEqn.solve();
    }

    updateDerivedFields();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
