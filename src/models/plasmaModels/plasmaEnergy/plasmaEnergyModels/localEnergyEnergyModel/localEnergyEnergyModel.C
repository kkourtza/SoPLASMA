/*---------------------------------------------------------------------------*\
  File: localEnergyEnergyModel.C
  Part of: SoPLASMA
  Copyright (C) 2026
  License: GNU General Public License v3 or later
\*---------------------------------------------------------------------------*/

#include "localEnergyEnergyModel.H"
#include "plasmaSpecies.H"
#include "plasmaConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvc.H"

namespace Foam
{

defineTypeNameAndDebug(localEnergyEnergyModel, 0);

addToRunTimeSelectionTable
(
    plasmaEnergyModel,
    localEnergyEnergyModel,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

localEnergyEnergyModel::localEnergyEnergyModel
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex,
    const volVectorField& E
)
:
    plasmaEnergyModel(modelName, dict, mesh, species, specieIndex, E),

    nEps_
    (
        IOobject
        (
            "nEps_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimVolume, 0.0)
    ),

    // Registered under a name the tabulated evaluators can find. This is the
    // whole integration mechanism: TabulatedProperty1D resolves its
    // `lookupVariable` from the registry BY NAME, so a case moves from LFA to
    // LMEA by pointing lookupVariable at this field and at the *_vs_meanE
    // tables the Boltzmann solver already writes.
    meanE_
    (
        IOobject
        (
            "meanE_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),

    T_
    (
        IOobject
        (
            "T_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTemperature, 300.0)
    ),

    muEf_
    (
        IOobject("muE_lmea", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0, -1, 0), 0.0)
    ),
    DEf_
    (
        IOobject("DE_lmea", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
    ),
    PlossN_
    (
        IOobject("PlossN_lmea", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0, 0, 0), 0.0)
    )
{
    const dictionary& c = dict.subOrEmptyDict("localEnergyCoeffs");

    // See nEfloor_: eps_bar = n_eps/n_e is 0/0 in the quiescent far field.
    nEfloor_  = c.getOrDefault<scalar>("electronDensityFloor", 1.0);
    meanEmin_ = c.getOrDefault<scalar>("meanEnergyMin", 1e-3);
    meanEmax_ = c.getOrDefault<scalar>("meanEnergyMax", 100.0);

    if (!c.found("mobility") || !c.found("diffusivity"))
    {
        FatalIOErrorInFunction(dict)
            << "energyModel `localEnergy` needs `mobility` and `diffusivity`"
            << " in localEnergyCoeffs, tabulated against MEAN ENERGY." << nl
            << "    These must be the *_vs_meanE tables, not the *_vs_reducedE"
            << " ones: evaluating transport at E/N is precisely the LFA"
            << " closure that LMEA replaces." << nl
            << exit(FatalIOError);
    }

    muE_ = plasmaPropertyEvaluator::New
    (
        c.subDict("mobility").get<word>("type"), c.subDict("mobility"), mesh,
        "ElectronMobility", dimensionSet(0, 2, -1, 0, 0, -1, 0)
    );

    DE_ = plasmaPropertyEvaluator::New
    (
        c.subDict("diffusivity").get<word>("type"), c.subDict("diffusivity"),
        mesh, "ElectronDiffusivity", dimensionSet(0, 2, -1, 0, 0, 0, 0)
    );

    // Power losses. `tabulated` (the EEDF's own Pelastic/Pinelastic) is the
    // default because its error direction is the safe one: a subset mechanism
    // UNDER-counts inelastic loss, T_e rises, and ionisation is exponential in
    // T_e -- so the `mechanism` route runs away where the table merely
    // mis-states the composition.
    if (c.found("powerLoss"))
    {
        const dictionary& pl = c.subDict("powerLoss");

        if (pl.found("elastic"))
        {
            Pelastic_ = plasmaPropertyEvaluator::New
            (
                pl.subDict("elastic").get<word>("type"),
                pl.subDict("elastic"), mesh,
                "ElasticPowerLoss", dimensionSet(0, 3, -1, 0, 0, 0, 0)
            );
        }
        if (pl.found("inelastic"))
        {
            Pinelastic_ = plasmaPropertyEvaluator::New
            (
                pl.subDict("inelastic").get<word>("type"),
                pl.subDict("inelastic"), mesh,
                "InelasticPowerLoss", dimensionSet(0, 3, -1, 0, 0, 0, 0)
            );
        }
    }

    Info<< "  energyModel localEnergy (LMEA) for `"
        << species.speciesNames()[specieIndex] << "`:" << nl
        << "    transported n_eps, mean energy published as `"
        << meanE_.name() << "`" << nl
        << "    n_e floor " << nEfloor_ << " 1/m^3, mean energy clamped to ["
        << meanEmin_ << ", " << meanEmax_ << "] eV" << endl;

    updateDerived();
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void localEnergyEnergyModel::updateDerived()
{
    const volScalarField& ne = species_.numberDensity(specieIndex_);

    scalarField& me = meanE_.primitiveFieldRef();
    scalarField& Tf = T_.primitiveFieldRef();

    const scalarField& nEpsI = nEps_.primitiveField();
    const scalarField& neI = ne.primitiveField();

    // eV -> K: T = (2/3) eps / k_B, with eps in eV.
    const scalar eVtoK =
        2.0/3.0*constant::plasma::eCharge.value()
      / constant::plasma::kappaBoltzmann.value();

    forAll(me, c)
    {
        // Floor, not a VSMALL guard: dividing by a near-zero density and then
        // clamping is how the cross-term sink produced a SIGFPE. Compare
        // first, divide second.
        const scalar n = max(neI[c], nEfloor_);

        me[c] = min(max(nEpsI[c]/n, meanEmin_), meanEmax_);
        Tf[c] = me[c]*eVtoK;
    }

    meanE_.correctBoundaryConditions();
    T_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void localEnergyEnergyModel::correct()
{
    updateDerived();

    muE_->correct(muEf_);
    DE_->correct(DEf_);

    PlossN_ == dimensionedScalar("zero", PlossN_.dimensions(), 0.0);

    if (Pelastic_)
    {
        volScalarField Pe(PlossN_);
        Pelastic_->correct(Pe);
        PlossN_ += Pe;
    }
    if (Pinelastic_)
    {
        volScalarField Pi(PlossN_);
        Pinelastic_->correct(Pi);
        PlossN_ += Pi;
    }
}


tmp<fvScalarMatrix> localEnergyEnergyModel::eEqn() const
{
    const volScalarField& ne = species_.numberDensity(specieIndex_);

    // Energy flux moments of the two-term expansion: the energy drifts and
    // diffuses at 5/3 of the particle rates (Hagelaar & Kroesen Eqs. 10-13).
    const dimensionedScalar fiveThirds("5/3", dimless, 5.0/3.0);

    // Electrons drift AGAINST E, so the energy convective velocity is
    // -mu_e E. Sign carried here once, explicitly.
    const surfaceScalarField phiEps
    (
        "phiEps",
        fiveThirds*fvc::flux(-muEf_*E_)
    );

    const volScalarField DEps("DEps", fiveThirds*DEf_);

    // JOULE HEATING, explicit in phase A: -e Gamma_e . E, which for electrons
    // drifting against the field is a positive power input. Written from the
    // drift flux rather than from the full Gamma_e, because the diffusive part
    // carries E.grad(n_e) -- the sign-indefinite term that makes Hagelaar's
    // bracket not unconditionally damping.
    const volScalarField jouleHeating
    (
        "jouleHeating",
        muEf_*ne*magSqr(E_)
    );

    // LOSS, implicit. P_loss = n_e N (Pelastic + Pinelastic)(eps_bar) is a
    // genuine sink, so it goes on the diagonal as a rate: L = P_loss/n_eps.
    // Only the unambiguous sinks are made implicit here -- see the header on
    // why the full Hagelaar bracket is not adopted in phase A.
    // Background gas density: a dimensionedScalar here, since the background
    // is a uniform reservoir rather than a transported field.
    const volScalarField Ploss(PlossN_*ne*species_.backgroundDensity());

    // Rate form, guarded: below the floor the cell holds no electrons worth
    // damping and an unbounded L is exactly the anti-damping shape that broke
    // the Rosenbrock controller.
    volScalarField Lrate
    (
        IOobject("Lrate_lmea", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    );

    {
        scalarField& L = Lrate.primitiveFieldRef();
        const scalarField& P = Ploss.primitiveField();
        const scalarField& nE = nEps_.primitiveField();

        forAll(L, c)
        {
            L[c] = (nE[c] > VSMALL) ? max(P[c]/nE[c], 0.0) : 0.0;
        }
    }

    tmp<fvScalarMatrix> tEqn
    (
        new fvScalarMatrix
        (
            fvm::ddt(nEps_)
          + fvm::div(phiEps, nEps_)
          - fvm::laplacian(DEps, nEps_)
          + fvm::Sp(Lrate, nEps_)
         ==
            jouleHeating
        )
    );

    return tEqn;
}


tmp<volScalarField> localEnergyEnergyModel::T() const
{
    return T_;
}


const dimensionedScalar& localEnergyEnergyModel::Tvalue() const
{
    return plasmaEnergyModel::Tvalue();
}


} // End namespace Foam

// ************************************************************************* //
