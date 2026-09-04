/*---------------------------------------------------------------------------*\
  File: multiRegionPoisson.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::multiRegionPoisson.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "mappedPatchBase.H"

#include "multiRegionPoisson.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(multiRegionPoisson, 0);
addToRunTimeSelectionTable
(
    electromagneticsModel,
    multiRegionPoisson,
    dictionary
);

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void multiRegionPoisson::updateDerivedFields()
{
    phiE_ = -fvc::snGrad(ePotential_) * mesh_.magSf();

    if (EScheme_ == "reconstruct")
    {
        E_ = fvc::reconstruct(phiE_);
    }
    else // "grad"
    {
        E_ = -fvc::grad(ePotential_);
    }
    E_.correctBoundaryConditions();

    Emag_ = mag(E_);
    Emag_.correctBoundaryConditions();

    // reducedE = |E|/N NEEDS A GAS. Without one the divisor is zero and this
    // is a guaranteed SIGFPE, which is exactly what the shipped
    // multiRegionElectrostaticFoam tutorials do: they have no species model, so
    // nothing supplies a density.
    //
    // REGRESSION, dated. Before fe827ae (2026-08-11, "One owner for the gas
    // density") backgroundDensityUniform_ defaulted to 2.5e25 via
    // `coeffs.getOrDefault<scalar>("backgroundDensity", 2.5e25)`. That commit
    // removed the default and deprecated the key -- correctly, since a silent
    // 2.5e25 in a case that never stated one is invented physics -- but left
    // this division unguarded, so every gas-free multi-region case has crashed
    // since. Nothing caught it because no multi-region case is exercised.
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

        Info<< "multiRegionPoisson: no background gas density, so `reducedE`"
               " stays ZERO." << nl
            << "    |E|/N is undefined without a gas. This is expected for a"
               " pure electrostatics" << nl
            << "    case; if this run has a plasma, its species model has not"
               " published a density" << nl
            << "    and every rate coefficient keyed on reducedE would be"
               " wrong." << endl;
    }

    for (dielectricRegion& reg : dielectrics_)
    {
        reg.updateE();
    }
}

void multiRegionPoisson::solveCoupled()
{
    Info<< "Solving for ePotential in coupled regions "
        << "(monolithically)" << endl;

    // ePotential_.correctBoundaryConditions();
    // for (dielectricRegion& reg : dielectrics_)
    // {
    //     reg.ePotential().correctBoundaryConditions();
    // }

    for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        fvScalarMatrix gasEqn
        (
            fvm::laplacian(epsilon_, ePotential_)
         == -chargeDensity_
        );

        fvMatrixAssemblyPtr_->addFvMatrix(gasEqn);

        for (dielectricRegion& reg : dielectrics_)
        {
            fvMatrixAssemblyPtr_->addFvMatrix(reg.PoissonMatrix().ref());
        }

        fvMatrixAssemblyPtr_->solve();

        ePotential_.correctBoundaryConditions();
        for (dielectricRegion& reg : dielectrics_)
        {
            reg.ePotential().correctBoundaryConditions();
        }

        fvMatrixAssemblyPtr_->clear();
    }
}

void multiRegionPoisson::solveCoupled
(
    const volScalarField& electricalConductivity,
    const volScalarField& diffusiveChargeSource
)
{
    Info<< "Solving for ePotential in coupled regions "
        << "(monolithically)" << endl;

    const dimensionedScalar deltaT = mesh_.time().deltaT();

    // ePotential_.correctBoundaryConditions();
    // for (dielectricRegion& reg : dielectrics_)
    // {
    //     reg.ePotential().correctBoundaryConditions();
    // }

    const volScalarField effEps(epsilon_ + deltaT * electricalConductivity);
    const volScalarField rhsSource(-chargeDensity_ - deltaT * diffusiveChargeSource);


    for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        fvScalarMatrix gasEqn
        (
            fvm::laplacian(effEps, ePotential_) == rhsSource
        );

        fvMatrixAssemblyPtr_->addFvMatrix(gasEqn);

        for (dielectricRegion& reg : dielectrics_)
        {
            fvMatrixAssemblyPtr_->addFvMatrix(reg.PoissonMatrix().ref());
        }

        fvMatrixAssemblyPtr_->solve();

        ePotential_.correctBoundaryConditions();
        for (dielectricRegion& reg : dielectrics_)
        {
            reg.ePotential().correctBoundaryConditions();
        }

        fvMatrixAssemblyPtr_->clear();
    }
}

void multiRegionPoisson::solveSegregated()
{
    Info<< "Solving for ePotential in non-coupled regions "
        << "(segregated)" << endl;

    for (int iter = 1; iter <= maxNonCoupledIterations_; ++iter)
    {
        if (maxNonCoupledIterations_ > 1)
        {
            Info<< "  -ePotential iteration (" << iter << "/"
                << maxNonCoupledIterations_ << ")" << endl;
        }

        bool allOK = true;

        {
            scalar resGas = 0.0;

            for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
            {
                fvScalarMatrix gasEqn
                (
                    fvm::laplacian(epsilon_, ePotential_)
                 == -chargeDensity_
                );

                if (nonOrth < nNonOrthCorr_)
                    gasEqn.relax();

                Info<< "    -Solving for ePotential (gas region: "
                    << mesh_.name() << ")" << endl;

                auto sp = gasEqn.solve();

                if (nonOrth == 0)
                    resGas = sp.initialResidual();

                ePotential_.correctBoundaryConditions();
            }

            if (resGas >= nonCoupledTolerance_) allOK = false;
        }

        for (dielectricRegion& reg : dielectrics_)
        {
            scalar resDiel = 0.0;

            reg.ePotential().correctBoundaryConditions();

            for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
            {
                fvScalarMatrix dielEqn
                (
                    fvm::laplacian(reg.epsilon(), reg.ePotential())
                );

                if (nonOrth < nNonOrthCorr_)
                    dielEqn.relax();

                Info<< "    -Solving for ePotential (dielectric: "
                    << reg.mesh().name() << ")" << endl;

                auto sp = dielEqn.solve();

                if (nonOrth == 0)
                    resDiel = sp.initialResidual();

                reg.ePotential().correctBoundaryConditions();
            }

            if (resDiel >= nonCoupledTolerance_) allOK = false;
        }

        ePotential_.correctBoundaryConditions();

        if (allOK)
        {
            if (maxNonCoupledIterations_ > 1)
            {
                Info<< ">>> ePotential converged in " << iter
                    << " iterations." << endl;
            }
            break;
        }
        else if (iter == maxNonCoupledIterations_)
        {
            Info<< ">>> WARNING: ePotential did NOT converge after "
                << iter << " iterations." << endl;
        }
    }
}

void multiRegionPoisson::solveSegregated
(
    const volScalarField& electricalConductivity,
    const volScalarField& diffusiveChargeSource
)
{
    Info<< "Solving for ePotential in non-coupled regions "
        << "(segregated)" << endl;

    const dimensionedScalar deltaT = mesh_.time().deltaT();

    const volScalarField effEps(epsilon_ + deltaT * electricalConductivity);
    const volScalarField rhsSource(-chargeDensity_ - deltaT * diffusiveChargeSource);


    for (int iter = 1; iter <= maxNonCoupledIterations_; ++iter)
    {
        if (maxNonCoupledIterations_ > 1)
        {
            Info<< "  -ePotential iteration (" << iter << "/"
                << maxNonCoupledIterations_ << ")" << endl;
        }

        bool allOK = true;

        for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
        {
            fvScalarMatrix gasEqn
            (
                fvm::laplacian(effEps, ePotential_) == rhsSource
            );

            if (nonOrth < nNonOrthCorr_)
            {
                gasEqn.relax();
            }

            Info<< "    -Solving for ePotential (gas region: "
                << mesh_.name() << ")" << endl;

            scalar resGas = gasEqn.solve().initialResidual();

            ePotential_.correctBoundaryConditions();

            if (nonOrth == nNonOrthCorr_)
            {
                if (resGas >= nonCoupledTolerance_) allOK = false;
            }

            for (dielectricRegion& reg : dielectrics_)
            {
                fvScalarMatrix dielEqn
                (
                    fvm::laplacian(reg.epsilon(), reg.ePotential())
                );

                if (nonOrth < nNonOrthCorr_)
                {
                    dielEqn.relax();
                }

                Info<< "    -Solving for ePotential (dielectric: "
                    << reg.mesh().name() << ")" << endl;

                scalar resDiel = dielEqn.solve().initialResidual();

                reg.ePotential().correctBoundaryConditions();

                if (nonOrth == nNonOrthCorr_)
                {
                    if (resDiel >= nonCoupledTolerance_) allOK = false;
                }
            }
        }

        if (allOK)
        {
            if (maxNonCoupledIterations_ > 1)
            {
                Info<< ">>> ePotential converged in " << iter
                    << " iterations." << endl;
            }
            break;
        }
        else if (iter == maxNonCoupledIterations_)
        {
            Info<< ">>> WARNING: ePotential did NOT converge after "
                << iter << " iterations." << endl;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
multiRegionPoisson::multiRegionPoisson
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
    backgroundDensityFieldPtr_(nullptr),
    dielectrics_(dielectricMeshes.size()),
    coupled_(false),
    fvMatrixAssemblyPtr_(nullptr),
    maxNonCoupledIterations_(100),
    nonCoupledTolerance_(1e-6)
{
    if (dielectricMeshes.size() == 0)
    {
        FatalErrorInFunction
            << "No dielectric regions were detected, but you are using the "
            << "'multiRegionPoisson' model." << nl
            << "This model requires at least one dielectric mesh." << nl
            << "Switch to 'singleRegionPoisson' model in your properties "
            << "file if you only have a gas region." << nl
            << exit(FatalError);
    }

    // The `<type>Coeffs` sub-dictionary of constant/electromagneticsProperties
    // is DEPRECATED and therefore optional. The numerics now come from
    // system/plasmaSimulationControls (`poisson`), uniform over every region
    // the potential is solved on, and each region's permittivity from its own
    // constant/<region>/electricalProperties. This dictionary is consulted only
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
    maxNonCoupledIterations_ = num.maxNonCoupledIterations;
    nonCoupledTolerance_ = num.nonCoupledTolerance;

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

    // Build dielectric regions.
    //
    // Each one now reads its own permittivity from
    // constant/<region>/electricalProperties, so ADDING A REGION MEANS ADDING A
    // DIRECTORY -- nothing global is edited. Previously every region needed a
    // hand-written sub-dictionary in this file, which was a per-region property
    // stored in a global file and was fatal if forgotten.
    forAll(dielectricMeshes_, i)
    {
        dielectrics_.set
        (
            i,
            new dielectricRegion
            (
                dielectricMeshes_[i],
                readEpsilonR(dielectricMeshes_[i], coeffs),
                EScheme_
            )
        );
    }

    // Detect implicit coupling
    bool anyImplicit = false;

    forAll(ePotential_.boundaryField(), patchI)
    {
        if (ePotential_.boundaryField()[patchI].useImplicit())
        {
            anyImplicit = true;
            break;
        }
    }

    if (!anyImplicit)
    {
        for (const auto& reg : dielectrics_)
        {
            forAll(reg.ePotential().boundaryField(), patchI)
            {
                if (reg.ePotential().boundaryField()[patchI].useImplicit())
                {
                    anyImplicit = true;
                    break;
                }
            }
            if (anyImplicit) break;
        }
    }

    if (anyImplicit)
    {
        coupled_ = true;
        Info<< "    Implicit coupling detected." << endl;

        // THE LINEAR SOLVER MUST BE ASSEMBLY-AWARE, and the failure otherwise
        // does not name the cause.
        //
        // Monolithic coupling assembles every region into one matrix whose mesh
        // is an `lduPrimitiveMeshAssembly`, not an `fvMesh`. GAMG's DEFAULT
        // agglomerator (`faceAreaPair`) does `refCast<const fvMesh>` on it and
        // the run aborts with
        //
        //     Attempt to cast type lduPrimitiveMeshAssembly to type fvMesh
        //
        // which reads like an internal error and is a solver-configuration
        // error. OpenFOAM ships `assembledFaceAreaPair` for exactly this case
        // (TypeName in src/finiteVolume/lduPrimitiveMeshAssembly).
        //
        // Checked here because monolithic is now the DEFAULT for a region
        // interface, so a case that never mentioned coupling can meet this --
        // and the abort would give it no way to connect the two.
        {
            const dictionary& sol =
                mesh_.solution().subOrEmptyDict("solvers");
            const dictionary eq = sol.subOrEmptyDict("ePotential");
            const word lin(eq.getOrDefault<word>("solver", word::null));

            if (lin == "GAMG")
            {
                // DERIVED, NOT REQUIRED (G1). The correct agglomerator is
                // fully determined by the coupling mode, so it is not a
                // choice the user can usefully make: monolithic coupling
                // needs `assembledFaceAreaPair` and nothing else works. If
                // the case did not state one, it is set here and reported.
                //
                // An EXPLICIT wrong value is still fatal. Silently overriding
                // what a user actually wrote would hide their mistake instead
                // of naming it, and the two cases are distinguishable because
                // the dictionary records which keys it was given.
                if (!eq.found("agglomerator"))
                {
                    // MUTATE THE CACHED SOLVER DICT, not the file-level one.
                    // `solution` keeps `solvers_` as a COPY taken in read(),
                    // and `fvMatrix::solve()` goes through solverDict(), so a
                    // set() on subDict("solvers") is silently ineffective --
                    // measured 2026-09-04: the derivation reported success and
                    // the solve still aborted on `faceAreaPair`.
                    //
                    // Caveat, stated rather than discovered later: a runtime
                    // re-read of fvSolution reloads `solvers_` from disk and
                    // would drop this. The entry is derived every construction,
                    // so it returns on the next rebuild of the model.
                    const_cast<dictionary&>
                    (
                        mesh_.solution().solverDict("ePotential")
                    ).set("agglomerator", word("assembledFaceAreaPair"));

                    Info<< "    ePotential: `agglomerator"
                        << " assembledFaceAreaPair` DERIVED -- monolithic"
                        << " coupling assembles every" << nl
                        << "      region into one lduPrimitiveMeshAssembly,"
                        << " which GAMG's default `faceAreaPair`" << nl
                        << "      cannot agglomerate. Not a case setting;"
                        << " state it only to override." << endl;
                }

                // NOTE `eq` is a BY-VALUE copy of the sub-dictionary, so it
                // does NOT see the set() above -- validating through it after
                // deriving compared against the stale copy and aborted on the
                // value it had just corrected. Measured 2026-09-04. Hence
                // derive-OR-validate, never both.
                const word agg
                (
                    eq.getOrDefault<word>("agglomerator", "faceAreaPair")
                );

                if (eq.found("agglomerator") && agg != "assembledFaceAreaPair")
                {
                    FatalErrorInFunction
                        << "ePotential STATES `solver GAMG` with"
                        << " `agglomerator " << agg
                        << "`, but the regions are coupled" << nl
                        << "    MONOLITHICALLY, so the matrix is an"
                        << " lduPrimitiveMeshAssembly rather than an fvMesh."
                        << nl << nl
                        << "    Either DELETE the `agglomerator` entry --"
                        << " it is derived when absent -- or set it to:" << nl
                        << "        agglomerator    assembledFaceAreaPair;"
                        << nl << nl
                        << "    Without it the solve aborts with \"Attempt to"
                        << " cast type lduPrimitiveMeshAssembly to type" << nl
                        << "    fvMesh\", which names neither the solver nor"
                        << " the coupling." << nl << nl
                        << "    Alternatives that need no agglomerator:"
                        << " PBiCGStab/DILU (the operator is NOT symmetric"
                        << nl
                        << "    under `poissonScheme semiImplicit`), or"
                        << " PCG/DIC under `explicit`." << nl << nl
                        << "    Or ask for segregated coupling explicitly with"
                        << " `useImplicit false` on the interface" << nl
                        << "    patches -- note that lags the interface and was"
                        << " measured to plateau the outer residual" << nl
                        << "    near 1e-7 on a needle geometry." << nl
                        << exit(FatalError);
                }
            }
        }

        // Validate all mapped patches are consistent
        forAll(ePotential_.boundaryField(), patchI)
        {
            const fvPatch& p = ePotential_.boundaryField()[patchI].patch();

            if (isA<mappedPatchBase>(p)
             && !ePotential_.boundaryField()[patchI].useImplicit())
            {
                FatalErrorInFunction
                    << "Mixed coupling detected on patch '" << p.name()
                    << "'." << nl << "In Monolithic mode, ALL interface "
                    << "patches must set 'useImplicit true'." << nl
                    << exit(FatalError);
            }
        }

        for (const auto& reg : dielectrics_)
        {
            forAll(reg.ePotential().boundaryField(), patchI)
            {
                const fvPatch& p =
                    reg.ePotential().boundaryField()[patchI].patch();

                if (isA<mappedPatchBase>(p)
                 && !reg.ePotential().boundaryField()[patchI].useImplicit())
                {
                    FatalErrorInFunction
                        << "Mixed coupling detected in region '"
                        << reg.mesh().name() << "' on patch '" << p.name()
                        << "'." << nl << "All mapped patches must be "
                        << "consistent (all implicit or all explicit)."
                        << exit(FatalError);
                }
            }
        }

        Info<< "Regions are COUPLED; "
            << "ePotential will be solved MONOLITHICALLY." << nl
            << "Assembling ePotential fvMatrixAssembly" << nl << endl;

    fvMatrixAssemblyPtr_.reset
    (
        new fvMatrix<scalar>
        (
            ePotential_,
            dimensionSet(0, 0, 1, 0, 0, 1, 0)
        )
    );
    
    }
    else
    {
        coupled_ = false;
        Info<< "Regions are NOT coupled; ePotential will be solved in "
            << "SEGREGATED way" << nl << endl;
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void multiRegionPoisson::solve()
{
    if (coupled_)
    {
        solveCoupled();
    }
    else
    {
        solveSegregated();
    }

    // Explicit scheme: the operator is eps, so psi needs no effective
    // permittivity. Before updateDerivedFields, which reads the potential.
    correctFloatingElectrode(nullptr);

    updateDerivedFields();
}


void multiRegionPoisson::solve
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

    if (coupled_)
    {
        solveCoupled(electricalConductivity, diffusiveChargeSource);
    }
    else
    {
        solveSegregated(electricalConductivity, diffusiveChargeSource);
    }

    // semiImplicit: rebuilt from the SAME operator the solve used. The
    // expression must match solveCoupled/solveSegregated exactly -- gas gets
    // eps + dt*sigma, dielectrics keep their constant eps (an insulator has no
    // conductivity), which unitPotentialField's overload does.
    {
        const volScalarField effEps
        (
            epsilon_ + mesh_.time().deltaT()*electricalConductivity
        );

        correctFloatingElectrode(&effEps);
    }

    updateDerivedFields();
}

const dielectricRegion& multiRegionPoisson::dielectric
(
    const word& name
) const
{
    forAll(dielectrics_, i)
    {
        if (dielectrics_[i].mesh().name() == name)
        {
            return dielectrics_[i];
        }
    }

    wordList availableNames(dielectrics_.size());
    forAll(dielectrics_, i)
    {
        availableNames[i] = dielectrics_[i].mesh().name();
    }

    FatalErrorInFunction
        << "Dielectric region '" << name << "' not found." << nl
        << "Available: " << availableNames << nl
        << exit(FatalError);

    return dielectrics_[0]; // unreachable
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
