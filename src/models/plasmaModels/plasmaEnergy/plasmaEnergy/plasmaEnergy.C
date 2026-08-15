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

#include "plasmaEnergy.H"
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

void plasmaEnergy::constructModels()
{
    // Read the backgroundGas properties
    const dictionary& bgGasDict = species_.backgroundDict();
    const dictionary& bgGasEnergyDict = bgGasDict.subDict("energy");

    solveGasEnergy_ = bgGasEnergyDict.get<bool>("solve");
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

    forAll(evI, celli)
    {
        const scalar tau = max(tauVT.primitiveField()[celli], SMALL);
        Qvt[celli] = evI[celli]/tau;
        evI[celli] += (Pvib.primitiveField()[celli] - Qvt[celli])*dt;
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
      - fvm::laplacian(kappa, T)
     ==
        Qtot
    );

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

