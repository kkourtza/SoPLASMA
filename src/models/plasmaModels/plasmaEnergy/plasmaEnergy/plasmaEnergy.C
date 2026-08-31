/*---------------------------------------------------------------------------*\
  File: plasmaEnergy.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaEnergy.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaOuterRelaxation.H"
#include "plasmaEnergy.H"
#include "localEnergyEnergyModel.H"
#include "plasmaEnergyModel.H"

#include "plasmaSimulationProfiler.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaEnergy, 0);

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * *  //

Foam::scalar Foam::plasmaEnergy::maxEnergyRate() const
{
    // The ACTUAL energy-transport rate [1/s], from the energy flux itself.
    //
    // This replaces a 5/3 factor applied to the SPECIES rate, which was wrong
    // twice over: it assumed the Maxwellian mu_eps/mu_e = 5/3 (measured
    // 1.24-1.57 on this air set) and it used a mobility evaluated at the
    // reduced field rather than at the mean energy. Taking the rate from the
    // flux the equation actually convects with removes both.
    scalar r = 0;

    forAll(energyModels_, i)
    {
        if (isA<localEnergyEnergyModel>(energyModels_[i]))
        {
            r = max
            (
                r,
                refCast<const localEnergyEnergyModel>
                (
                    energyModels_[i]
                ).maxEnergyRate()
            );
        }
    }

    return r;
}


const Foam::localEnergyEnergyModel* Foam::plasmaEnergy::lmeaModel() const
{
    forAll(energyModels_, i)
    {
        if (isA<localEnergyEnergyModel>(energyModels_[i]))
        {
            return &refCast<const localEnergyEnergyModel>(energyModels_[i]);
        }
    }
    return nullptr;
}


Foam::scalar Foam::plasmaEnergy::maxEnergyRelaxationRate() const
{
    // Mirrors maxEnergyRate() -- the largest over every LMEA model present,
    // so a case with more than one transported energy is limited by the
    // stiffest of them rather than by whichever happens to be last.
    scalar r = 0;

    forAll(energyModels_, i)
    {
        if (isA<localEnergyEnergyModel>(energyModels_[i]))
        {
            r = max
            (
                r,
                refCast<const localEnergyEnergyModel>
                (
                    energyModels_[i]
                ).maxEnergyRelaxationRate()
            );
        }
    }

    return r;
}


void Foam::plasmaEnergy::solveSpeciesEnergy()
{
    forAll(energyModels_, i)
    {
        // Every model must be corrected -- that is what refreshes the derived
        // fields and the tabulated coefficients -- but only the models that
        // TRANSPORT energy return a matrix. The LFA family returns nullptr,
        // which is not a failure but a statement that their temperature is an
        // algebraic function of the local state.
        energyModels_[i].correct();

        tmp<fvScalarMatrix> tEqn = energyModels_[i].eEqn();

        if (tEqn.valid())
        {
            // LINEAR-SOLVER SETTINGS, borrowed from the species this energy
            // belongs to. Cases carry a regex entry `"n_.*"` covering every
            // transported species, which `nEps_e` does NOT match -- so without
            // this the run aborts with `Entry 'nEps_e' not found in
            // solvers`, and every case would have to declare a solver for an
            // equation it enabled with one keyword.
            //
            // Borrowing is right on the merits too: this is a transported
            // scalar on the same mesh with the same operator structure as the
            // species density, so it wants the same settings. An explicit
            // `nEps_<specie>` entry still wins, because subDict() resolves
            // exact names before regexes.
            const word sName = species_.speciesNames()[i];

            const dictionary& solvers =
                mesh_.solution().subDict("solvers");

            const word key =
                solvers.found("nEps_" + sName, keyType::REGEX)
              ? word("nEps_" + sName)
              : word("n_" + sName);

            tEqn.ref().relax();
            tEqn.ref().solve(solvers.subDict(key));

            // JOINT outer-loop relaxation: offer this corrector's energy
            // density. The coordinator applies ONE factor once every enrolled
            // field has contributed -- the electron density contributes from
            // plasmaTransport. This half was missing in the first
            // implementation, and relaxing the density alone was measured to
            // drive the factor to its floor without converging.
            //
            // BEFORE correct(): eps_bar and T_e are DERIVED from nEps_, so
            // they must follow the relaxed value, not the unrelaxed one.
            {
                plasmaOuterRelaxation* r =
                    plasmaOuterRelaxation::lookup(mesh_);
                if (r && r->active())
                {
                    volScalarField* ne = energyModels_[i].nEpsPtr();
                    if (ne) { r->contribute(*ne); }
                }
            }

            // Re-derive AFTER the solve: eps_bar and T_e follow from the new
            // energy density, and everything downstream reads them.
            energyModels_[i].correct();
        }
    }
}


bool plasmaEnergy::required(const plasmaSpecies& species)
{
    // Gas heating asked for?
    const dictionary& bg = species.backgroundDict();
    if (bg.subOrEmptyDict("energy").getOrDefault<bool>("solve", false))
    {
        return true;
    }

    // Or any species carrying its own energy equation? `gasTemperature` is the
    // default and needs none -- it is the historical "every species sits at
    // the gas temperature" behaviour. Anything else (localEnergy for LMEA)
    // does.
    forAll(species.speciesNames(), i)
    {
        if
        (
            species.speciesDict(i)
                .getOrDefault<word>("energyModel", "gasTemperature")
         != "gasTemperature"
        )
        {
            return true;
        }
    }

    return false;
}


void plasmaEnergy::constructModels()
{
    // Read the backgroundGas properties.
    //
    // subOrEmptyDict and getOrDefault, NOT subDict/get: since gas heating and
    // the per-species energy models were separated, a case may enable LMEA
    // (speciesProperties/<sp>/energyModel localEnergy) without wanting the gas
    // temperature solved at all, and such a case has no `energy` dictionary to
    // read. Demanding one would reinstate exactly the coupling that separation
    // removed.
    const dictionary& bgGasDict = species_.backgroundDict();
    const dictionary& bgGasEnergyDict = bgGasDict.subOrEmptyDict("energy");

    solveGasEnergy_ = bgGasEnergyDict.getOrDefault<bool>("solve", false);
    if (solveGasEnergy_)
    {
        isGasTempField_ = true;
        // READ_IF_PRESENT, not MUST_READ. A case that wants its own initial
        // condition and wall boundary conditions supplies 0/T_gas; one that
        // does not still runs, which is what makes gas heating switchable with
        // a single dictionary entry.
        //
        // The fallback is uniform at the dictionary temperature with
        // zeroGradient everywhere, i.e. ADIABATIC WALLS. That is a physical
        // assumption, not a neutral default, so it is announced rather than
        // assumed silently -- for a nanosecond pulse it is right (no heat
        // reaches a wall in 100 ns) and for a steady discharge it is not.
        const scalar T0 = bgGasEnergyDict.getOrDefault<scalar>("T", 300.0);
        const bool haveFile = IOobject
        (
            "T_gas", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE
        ).typeHeaderOk<volScalarField>(true);

        TgasFieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "T_gas",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("T_gas", dimTemperature, T0),
                "zeroGradient"
            )
        );

        if (!haveFile)
        {
            Info<< "plasmaEnergy: no 0/T_gas found; starting uniform at "
                << T0 << " K with zeroGradient (ADIABATIC) boundaries." << nl
                << "    Supply 0/T_gas to set your own initial and wall"
                << " conditions." << endl;
        }

        TgasValue_.value() = T0;

        eVibPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "e_vib",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimensionSet(1, -1, -2, 0, 0), Zero),
                "zeroGradient"
            )
        );

        kappaGas_   = bgGasEnergyDict.getOrDefault<scalar>("kappa", 0.026);
        tauVTfixed_ = bgGasEnergyDict.getOrDefault<scalar>("tauVT", -1);
        pGasAtm_    =
            bgGasEnergyDict.getOrDefault<scalar>("pressure", 101325.0)/101325.0;
    }
    else
    {
        if (bgGasEnergyDict.found("T"))
        {
            isGasTempField_ = false;
            TgasValue_.value() = bgGasEnergyDict.get<scalar>("T");
            TgasFieldPtr_.reset(nullptr);
        }
        else 
        {
            isGasTempField_ = true;
            TgasFieldPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "T_gas", 
                        mesh_.time().timeName(), 
                        mesh_, 
                        IOobject::MUST_READ, 
                        IOobject::AUTO_WRITE
                    ),
                    mesh_
                )
            );
            TgasValue_ = dimensionedScalar("Tgas_dummy", dimTemperature, 300.0);
        }
    }

    // Loop over species and create an energy model for each one
    for (label i = 0; i < species_.nSpecies(); ++i)
    {
        const word& sName = species_.speciesNames()[i];
        const dictionary& sDict = species_.speciesDict(sName);

        // DEFAULT rather than demand. Requiring the entry made constructing
        // this class break every case in the repository, because none of them
        // declare it -- and the historical behaviour it would be replacing is
        // exactly "every species sits at the gas temperature", which is what
        // gasTemperature means. So that is the default, and a case only writes
        // the entry when it wants something else.
        //
        // A required entry whose only sensible value is the historical one is
        // not a safety feature, it is a migration tax.
        const word modelName =
            sDict.getOrDefault<word>("energyModel", "gasTemperature");

        // Optional for the same reason energyModel is: the default model
        // needs no coefficients, so demanding the sub-dictionary would be a
        // second migration tax behind the first.
        const dictionary& modelDict =
            sDict.subOrEmptyDict("energyModelCoeffs");

        // Construct the model using the runtime selection system
        energyModels_.set
        (
            i,
            plasmaEnergyModel::New
            (
                modelName,
                modelDict, 
                mesh_, 
                species_, 
                i, 
                E_
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaEnergy::plasmaEnergy
(
    plasmaSpecies& species,
    const fvMesh& mesh,
    const volVectorField& E
)
:
    regIOobject
    (
        IOobject
        (
            "plasmaEnergy",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    solveGasEnergy_(false),
    isGasTempField_(false),
    TgasValue_("Tgas", dimTemperature, 300.0),
    mesh_(mesh),
    species_(species),
    E_(E),
    energyModels_(species.nSpecies())
{
    constructModels();

    // ---- JOINT outer-loop relaxation -------------------------------------
    //
    // plasmaEnergy is constructed BEFORE plasmaTransport, so it creates the
    // coordinator; plasmaTransport's New() then finds the same instance.
    // Enrolment must happen here, before the time loop, because enrol() seeds
    // the previous iterate from the field's INITIAL value -- the only point at
    // which a pre-solve value is available without a before-solve hook.
    {
        // REGION-SCOPED FIRST, THEN THE TOP-LEVEL FILE.
        //
        // Reading with `mesh_` as the database resolves to
        // system/<region>/plasmaSimulationControls, and with READ_IF_PRESENT
        // there is NO fallback. In a MULTI-REGION case that file does not
        // exist -- only fvSchemes and fvSolution are copied per region -- so
        // the whole outerCoupling block was silently EMPTY here while
        // plasmaTimeControl (which reads with `runTime` as the database) read
        // it correctly from the top-level file.
        //
        // The result was a PARTIALLY honoured block: target/tolerance/
        // maxCorrectors worked, `outerScheme` did not, and the "is not
        // recognised" validation in plasmaOuterRelaxation was unreachable.
        // VERIFIED 2026-08-31 on the needle-DBD bed: `outerScheme
        // BOGUS_SCHEME_NAME` was ACCEPTED and the run proceeded on Aitken.
        // A user asking for Anderson on any multi-region case got Aitken with
        // no warning.
        //
        // Region-scoped is still tried first so a genuine per-region override
        // keeps working; the top-level file is the fallback, which is where
        // users actually write these settings.
        IOdictionary regionControls
        (
            IOobject
            (
                "plasmaSimulationControls",
                mesh_.time().system(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );

        IOdictionary globalControls
        (
            IOobject
            (
                "plasmaSimulationControls",
                mesh_.time().system(),
                mesh_.time(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );

        const dictionary& controls =
            regionControls.found("outerCoupling")
          ? static_cast<const dictionary&>(regionControls)
          : static_cast<const dictionary&>(globalControls);
        const dictionary& oc = controls.subOrEmptyDict("outerCoupling");

        // LMEA is the case the joint relaxation exists for, so it is the
        // default there: an electron-energy equation means the n_e / nEps_e
        // Picard pair exists and can enter the period-2 cycle that killed this
        // benchmark. A case can still switch it off explicitly.
        bool hasLMEA = false;
        forAll(energyModels_, i)
        {
            if (energyModels_.set(i) && energyModels_[i].nEpsPtr())
            {
                hasLMEA = true;
                break;
            }
        }

        plasmaOuterRelaxation& r =
            plasmaOuterRelaxation::New(mesh_, oc, hasLMEA);

        if (r.active())
        {
            forAll(energyModels_, i)
            {
                if (energyModels_.set(i))
                {
                    volScalarField* ne = energyModels_[i].nEpsPtr();
                    if (ne) { r.enrol(*ne); }
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void Foam::plasmaEnergy::solveGasEnergy
(
    const volScalarField& Qgas,
    const volScalarField& rhoCv,
    const volScalarField& Pvib,
    const volScalarField& tauVT,
    const scalar dt
)
{
    if (!solveGasEnergy_ || !TgasFieldPtr_) return;

    volScalarField& T = *TgasFieldPtr_;
    volScalarField& ev = *eVibPtr_;

    // VIBRATIONAL RESERVOIR, integrated first and explicitly.
    //
    // Explicit is adequate because the reservoir is filled by the pulse and
    // drained on tau_VT, and the timestep that resolves the pulse is already
    // far below tau_VT -- microseconds against nanoseconds. If that ever stops
    // being true the drain should go implicit; the check is dt vs tau_VT and
    // it is reported below.
    scalarField& evI = ev.primitiveFieldRef();
    scalarField Qvt(mesh_.nCells(), Zero);

    // Baseline the reservoir once per timestep, then integrate FROM that
    // baseline rather than from wherever the field currently sits. Identical
    // arithmetic on the first call of a step; repeatable on any later call,
    // which is what makes discarding and re-running a step safe. See the
    // header for why this field is the only one that needed it.
    const label ti = mesh_.time().timeIndex();
    if (ti != eVibTimeIndex_ || eVibStart_.size() != evI.size())
    {
        eVibTimeIndex_ = ti;
        eVibStart_ = evI;
    }

    forAll(evI, celli)
    {
        const scalar tau = max(tauVT.primitiveField()[celli], SMALL);
        Qvt[celli] = eVibStart_[celli]/tau;
        evI[celli] =
            eVibStart_[celli]
          + (Pvib.primitiveField()[celli] - Qvt[celli])*dt;
        evI[celli] = max(evI[celli], scalar(0));
    }
    ev.correctBoundaryConditions();

    volScalarField Qtot
    (
        IOobject
        (
            "Qtot", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        Qgas
    );
    Qtot.primitiveFieldRef() += Qvt;

    const dimensionedScalar kappa
    (
        "kappa", dimensionSet(1, 1, -3, -1, 0), kappaGas_
    );

    // THE EQUATION, as an fvScalarMatrix rather than a per-cell update, so it
    // inherits OpenFOAM's ddt schemes, boundary conditions and linear solvers
    // -- and so conduction is present the moment a gradient exists. On a
    // uniform field the laplacian contributes nothing, which is what allows
    // this to be checked against the 0-D reactor exactly.
    fvScalarMatrix TEqn
    (
        rhoCv*fvm::ddt(T)
     ==
        Qtot
    );

    // CONDUCTION ONLY WHEN kappa > 0.
    //
    // The start-up guard offers `kappa 0` as a way to drop conduction, and
    // until 2026-08-21 the equation assembled fvm::laplacian(kappa, T)
    // unconditionally -- so a case that followed that advice got past the
    // friendly guard and then died on a raw
    //
    //     Entry 'laplacian(kappa,T_gas)' not found in ... laplacianSchemes
    //
    // from deep inside OpenFOAM. Assembling a laplacian with a zero
    // coefficient also costs a scheme lookup and a matrix contribution for a
    // term that is identically zero. Now the advice and the code agree.
    if (kappaGas_ > 0)
    {
        TEqn -= fvm::laplacian(kappa, T);
    }

    TEqn.relax();
    TEqn.solve();

    T.max(dimensionedScalar("Tmin", dimTemperature, 100.0));
}


void plasmaEnergy::correct()
{
    // forAll(energyModels_, i)
    // {
    //     energyModels_[i].correct();
    // }
}

tmp<volScalarField> plasmaEnergy::Tgas() const
{
    if (isGasTempField_)
    {
        return *TgasFieldPtr_;
    }
    else
    {
        return tmp<volScalarField>::New
        (
            IOobject
            (
                "Tgas_tmp", 
                mesh_.time().timeName(), 
                mesh_
            ),
            mesh_,
            TgasValue_
        );
    }
}

const dimensionedScalar& plasmaEnergy::TgasValue() const
{
    if (isGasTempField_)
    {
        FatalErrorInFunction
            << "Requested TgasValue() scalar accessor, but background "
            << "temperature is a spatial field." << nl
            << "Use Tgas() instead." << abort(FatalError);
    }

    return TgasValue_;
}

bool plasmaEnergy::writeData(Ostream& os) const
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

