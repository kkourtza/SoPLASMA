/*---------------------------------------------------------------------------*\
  File: unitPotentialField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::unitPotentialField.

    Extracted 2026-09-03 from plasmaDischargeCurrent::computeWeightingField,
    UNCHANGED in behaviour. The extraction baseline, which this must reproduce
    exactly, is recorded in BASELINE-unitPotential-refactor.md at the repository
    root: C_g = 1.4757e-12 F at a relative error of 1.0948e-15 against the
    analytic series-stack capacitance, max|e_hat| = 1.66667 1/m.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
\*---------------------------------------------------------------------------*/

#include "unitPotentialField.H"

#include "electromagneticsModel.H"
#include "multiRegionPoisson.H"
#include "fixedValueFvPatchFields.H"
#include "mappedPatchBase.H"
#include "wedgePolyPatch.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unitPotentialField::unitPotentialField
(
    const fvMesh& gasMesh,
    const word& name,
    const word& unitPatch,
    const wordList& zeroPatches,
    const string& hint
)
:
    mesh_(gasMesh),
    name_(name),
    unitPatch_(unitPatch),
    zeroPatches_(zeroPatches),
    gas_(nullptr),
    dielectrics_(),
    coupled_(false),
    nDiel_(0),
    solved_(false),
    hint_(hint)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::volScalarField>
Foam::unitPotentialField::clone
(
    const volScalarField& ePot,
    const word& regionName
) const
{
    auto psi = autoPtr<volScalarField>::New
    (
        IOobject
        (
            name_,
            ePot.mesh().time().constant(),
            ePot.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ePot                            // clone: BCs come across intact
    );

    psi->primitiveFieldRef() = 0.0;
    psi->dimensions().reset(dimless);

    volScalarField::Boundary& bf = psi->boundaryFieldRef();

    forAll(bf, patchi)
    {
        const word& pname = bf[patchi].patch().name();

        const bool unit     = (pname == unitPatch_);
        const bool grounded = zeroPatches_.found(pname);

        if (!unit && !grounded)
        {
            // Not an electrode. Interface conditions need RETARGETING, not
            // just copying, and this is the subtle failure the plate2D test
            // exposed.
            //
            // `coupledElectricPotential` resolves its neighbour by NAME:
            //     phiNbrName_(dict.getOrDefault<word>("phiNbr", "ePotential"))
            // so a cloned interface BC on psi couples it to the REAL POTENTIAL
            // in the neighbour region instead of to psi. It also carries a
            // surface-charge source, which does not belong in a
            // unit-potential problem at all.
            //
            // MEASURED before the fix: the interior field was exactly right
            // (|e_hat| = 1.66667 everywhere) while the single cell layer
            // against the interface reached 20, inflating C_g by 3.76x. A
            // wrong answer with a perfect-looking interior.
            //
            // So the BC is REBUILT from its own dictionary with the field name
            // retargeted and the surface charge disabled ("none" is the
            // documented default that zeroes the term). Rebuilding through the
            // dictionary keeps the type, `useImplicit` and all sampling
            // information, none of which can be reconstructed by hand.
            if (isA<mappedPatchBase>(bf[patchi].patch().patch()))
            {
                OStringStream os;
                bf[patchi].write(os);

                IStringStream is(os.str());
                dictionary bcDict(is);

                bcDict.set("phiNbr", name_);
                bcDict.set("surfCharge", word("none"));
                bcDict.set("surfChargeNbr", word("none"));

                bf.set
                (
                    patchi,
                    fvPatchScalarField::New
                    (
                        bf[patchi].patch(),
                        psi->internalField(),
                        bcDict
                    )
                );
            }

            bf[patchi] == 0.0;
            continue;
        }

        // An electrode. Replace the BC outright rather than assigning into it:
        // the real ePotential electrode may be a time-varying condition (a
        // ramp, a table, a coded BC), which would re-evaluate on updateCoeffs
        // and silently overwrite the unit value this field depends on.
        bf.set
        (
            patchi,
            new fixedValueFvPatchScalarField
            (
                bf[patchi].patch(),
                psi->internalField()
            )
        );

        bf[patchi] == (unit ? 1.0 : 0.0);
    }

    // regionName is carried for diagnostics only; the field's own mesh already
    // determines where it lives.
    (void)regionName;

    return psi;
}


Foam::dictionary Foam::unitPotentialField::solverDict() const
{
    // NOT borrowed from ePotential, which was tried and failed: that entry is
    // tuned for TIME-STEPPED continuation, where each step starts from the
    // previous solution and a few GaussSeidel sweeps suffice. psi is a
    // ONE-SHOT cold solve from a zero field, and the tutorial's
    // smoothSolver/GaussSeidel stalled at residual 7e-3 after its 2000-
    // iteration cap -- giving a field with max|e_hat| 11.5 against a physical
    // 1.67, and a C_g wrong by 91%.
    //
    // A symmetric Laplacian wants a Krylov method. Defaults are set here, in
    // code, so a user gets a correct field without configuring something they
    // did not ask for; an explicit block named after the field in fvSolution
    // overrides them.
    dictionary defaults;
    defaults.add("solver", word("PCG"));
    defaults.add("preconditioner", word("DIC"));
    defaults.add("tolerance", 1e-12);
    defaults.add("relTol", 0.0);
    defaults.add("maxIter", 5000);

    const dictionary& solvers = mesh_.solution().subDict("solvers");

    return solvers.found(name_) ? solvers.subDict(name_) : defaults;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::unitPotentialField::wedgeRevolutionFactor(const fvMesh& m)
{
    // On a wedge mesh, m.V() holds only the theta-slice volumes, so every
    // extensive quantity -- charge, capacitance, current -- must be scaled to
    // the full 2*pi body of revolution. Taken from the geometry rather than a
    // dictionary entry, which is exactly the kind of thing that silently
    // disagrees with the mesh.
    const polyBoundaryMesh& pbm = m.boundaryMesh();

    forAll(pbm, patchi)
    {
        if (isA<wedgePolyPatch>(pbm[patchi]))
        {
            const wedgePolyPatch& wp =
                refCast<const wedgePolyPatch>(pbm[patchi]);

            const scalar halfAngle =
                Foam::acos(min(max(wp.cosAngle(), scalar(-1)), scalar(1)));

            const scalar theta = 2.0*halfAngle;

            if (theta > SMALL)
            {
                return constant::mathematical::twoPi/theta;
            }
        }
    }

    return 1.0;      // 3-D or planar 2-D: volumes are already what they claim
}


const Foam::volScalarField& Foam::unitPotentialField::gas() const
{
    if (!gas_)
    {
        FatalErrorInFunction
            << "unitPotentialField `" << name_ << "` has not been solved."
            << nl << "    Call solve() before asking for the field." << nl
            << exit(FatalError);
    }

    return gas_();
}


void Foam::unitPotentialField::solve(const electromagneticsModel& em)
{
    solveImpl(em, nullptr);
}


void Foam::unitPotentialField::solve
(
    const electromagneticsModel& em,
    const volScalarField& effEpsGas
)
{
    solveImpl(em, &effEpsGas);
}


void Foam::unitPotentialField::solveImpl
(
    const electromagneticsModel& em,
    const volScalarField* effEpsGas
)
{
    // UNIT potential, deliberately: psi = 1 on the chosen patch makes the field
    // time-independent even under a ramped or AC drive, which is what makes
    // computing it once legitimate. The applied voltage re-enters only where
    // the caller uses it.
    gas_ = clone(em.ePotential(), mesh_.name());

    const bool multi = isA<multiRegionPoisson>(em);

    const multiRegionPoisson* mrp =
        multi ? &refCast<const multiRegionPoisson>(em) : nullptr;

    nDiel_ = mrp ? mrp->nDielectrics() : 0;

    dielectrics_.clear();
    dielectrics_.setSize(nDiel_);

    for (label i = 0; i < nDiel_; ++i)
    {
        dielectrics_.set
        (
            i,
            clone(mrp->dielectric(i).ePotential(), mrp->dielectric(i).mesh().name())
        );
    }

    const dimensionedScalar& epsGas = em.epsilon();

    // DISCRETISATION scheme, borrowed from the potential for the same reason:
    // psi IS the potential's Laplacian operator, so it wants the same scheme,
    // and a case cannot be expected to declare `laplacian(epsilon,psi)` for
    // something it did not ask for. Cases must define the potential's entry
    // anyway, for the Poisson solve.
    //
    // If a case sets `default` in laplacianSchemes this name resolves to it, so
    // the borrow is harmless there too.
    const word psiScheme("laplacian(" + epsGas.name() + ",ePotential)");

    const dictionary psiSolverDict(solverDict());

    // MONOLITHIC across regions when the solver is coupled -- which is the
    // default, and the solver paper's headline contribution (Pasolari &
    // Kourtzanidis, arXiv:2607.05137): a monolithic multi-region Poisson
    // assembly over arbitrarily many curved conforming interfaces.
    //
    // Assembling psi the same way the real Poisson solve is assembled gives
    // continuous psi with continuous eps dpsi/dn across every interface BY
    // CONSTRUCTION, rather than by any interface handling written here.
    coupled_ = (nDiel_ > 0) && mrp->coupled();

    for (label nonOrth = 0; nonOrth <= em.nNonOrthCorr(); ++nonOrth)
    {
        if (coupled_)
        {
            fvScalarMatrix gasEqn
            (
                effEpsGas
              ? fvm::laplacian(*effEpsGas, gas_(), psiScheme)
              : fvm::laplacian(epsGas, gas_(), psiScheme)
            );

            // Dimensions taken FROM the equation, not hardcoded. The real
            // Poisson assembly uses the charge-sourced potential's dimensions;
            // psi is dimensionless, so laplacian(eps, psi) differs and
            // addFvMatrix rejects the mismatch outright -- which is how this
            // was found.
            fvMatrix<scalar> assembly(gas_(), gasEqn.dimensions());

            assembly.addFvMatrix(gasEqn);

            for (label i = 0; i < nDiel_; ++i)
            {
                const word dielScheme
                (
                    "laplacian(" + mrp->epsilon(i).name() + ",ePotential)"
                );

                fvScalarMatrix dielEqn
                (
                    fvm::laplacian(mrp->epsilon(i), dielectrics_[i], dielScheme)
                );
                assembly.addFvMatrix(dielEqn);
            }

            assembly.solve(psiSolverDict);

            gas_().correctBoundaryConditions();
            for (label i = 0; i < nDiel_; ++i)
            {
                dielectrics_[i].correctBoundaryConditions();
            }
        }
        else
        {
            // Single region, or a multi-region case the user configured as
            // segregated. Segregated psi would need its own outer loop to
            // converge the interfaces; rather than half-implement that, it is
            // refused below when dielectrics are present.
            fvScalarMatrix psiEqn
            (
                effEpsGas
              ? fvm::laplacian(*effEpsGas, gas_(), psiScheme)
              : fvm::laplacian(epsGas, gas_(), psiScheme)
            );
            psiEqn.solve(psiSolverDict);
            gas_().correctBoundaryConditions();
        }
    }

    if (nDiel_ > 0 && !coupled_)
    {
        FatalErrorInFunction
            << "The unit-potential field `" << name_ << "` requires the"
            << " MONOLITHIC (coupled) Poisson path when dielectric regions"
            << " are present." << nl
            << "    This case has " << nDiel_ << " dielectric region(s) with"
            << " region coupling disabled, so psi cannot be made continuous"
            << " across the interfaces here." << nl
            << "    Enable the coupled solve."
            << (hint_.empty() ? "" : " " + hint_) << nl
            << exit(FatalError);
    }

    solved_ = true;
}


Foam::scalar Foam::unitPotentialField::energyCapacitance
(
    const electromagneticsModel& em
) const
{
    if (!solved_)
    {
        FatalErrorInFunction
            << "unitPotentialField `" << name_ << "` has not been solved."
            << nl << exit(FatalError);
    }

    // C = INT eps |grad psi|^2 dV over the WHOLE domain, dielectrics included.
    //
    // Restricting the integral to the gas would put a bounding surface on the
    // dielectric face, where psi' is not zero -- and that surface term is
    // precisely what Morrow & Sato's derivation drops. Keeping the dielectric
    // interior keeps the bounding surface at electrodes and infinity.
    const volVectorField eHat(-fvc::grad(gas_()));

    scalar C = wedgeRevolutionFactor(mesh_)*em.epsilon().value()
             * gSum(mesh_.V()*magSqr(eHat.primitiveField()));

    const multiRegionPoisson* mrp =
        isA<multiRegionPoisson>(em)
      ? &refCast<const multiRegionPoisson>(em)
      : nullptr;

    for (label i = 0; i < nDiel_; ++i)
    {
        const fvMesh& rm = mrp->dielectric(i).mesh();
        const volVectorField eR(-fvc::grad(dielectrics_[i]));

        // Each region carries its own revolution factor: nothing guarantees a
        // dielectric mesh was built with the same wedge angle as the gas.
        C += wedgeRevolutionFactor(rm)*mrp->epsilon(i).value()
           * gSum(rm.V()*magSqr(eR.primitiveField()));
    }

    return C;
}


Foam::scalar Foam::unitPotentialField::surfaceCapacitance
(
    const electromagneticsModel& em,
    const word& patchName
) const
{
    if (!solved_)
    {
        FatalErrorInFunction
            << "unitPotentialField `" << name_ << "` has not been solved."
            << nl << exit(FatalError);
    }

    const label patchi = mesh_.boundaryMesh().findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch `" << patchName << "` not found on mesh `"
            << mesh_.name() << "`." << nl
            << "    Available: " << mesh_.boundaryMesh().names() << nl
            << exit(FatalError);
    }

    // C_self = + INT eps snGrad(psi) dA over the patch.
    //
    // THE SIGN, derived rather than guessed -- the first version of this had it
    // backwards and gave exactly -C, which the two-electrode test caught on its
    // first run (2026-09-03).
    //
    // fvPatch::snGrad() is n_patch . grad(psi) where n_patch points OUT OF THE
    // FLUID, i.e. INTO the metal. A Gauss surface ENCLOSING the electrode has
    // its own outward normal pointing the other way, into the fluid:
    //     n_enclosing = -n_patch
    // so with D = eps E = -eps grad(psi),
    //     Q = INT D . n_enclosing dA
    //       = INT (-eps grad psi) . (-n_patch) dA
    //       = + INT eps snGrad(psi) dA
    // Hence NO leading minus. The physical check that fixes it independently of
    // the algebra: psi = 1 on this patch and 0 elsewhere must put POSITIVE
    // charge on it, so C_self > 0.
    //
    // NOTE the permittivity here is the TRUE eps, not any effective operator
    // permittivity: the physical charge on a conductor is the integral of the
    // real D = eps E, whatever operator was used to obtain the field.
    const scalarField& magSf = mesh_.magSf().boundaryField()[patchi];

    const scalarField snG(gas_().boundaryField()[patchi].snGrad());

    return wedgeRevolutionFactor(mesh_)*em.epsilon().value()
         * gSum(snG*magSf);
}


// ************************************************************************* //
