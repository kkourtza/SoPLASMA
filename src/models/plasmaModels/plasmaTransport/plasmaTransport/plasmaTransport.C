/*---------------------------------------------------------------------------*\
  File: plasmaTransport.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaTransport.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaTransport.H"
#include "plasmaReactionRates.H"
#include "plasmaChemistry.H"
#include <cmath>
#include "IFstream.H"
#include "plasmaTransportModel.H"
#include "plasmaWallBC.H"
#include "photoionizationModel.H"

// Remove these headers later
#include "interpolationTable.H"
#include "fvm.H"
#include "fvc.H"
#include "vibRelax.H"
#include "janafMixture.H"
#include "plasmaEnergy.H"
#include "plasmaSimulationProfiler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaTransport, 0);

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * *  //

void plasmaTransport::constructTransportModels()
{
    Info<< "Constructing plasma transport models" << endl;

    // Loop over species and create its transport model
    for (label i = 0; i < species_.nSpecies(); ++i)
    {
        const word& sName = species_.speciesName(i);
        const dictionary& sDict = species_.speciesDict(i);

        if (!sDict.found("transportModel"))
        {
            FatalIOErrorInFunction(sDict)
                << "Species '" << sName
                << "' is missing required entry 'transportModel' in "
                << species_.dictName() << nl << exit(FatalIOError);
        }

        const word modelName(sDict.get<word>("transportModel"));
        const word coeffsName(modelName + "Coeffs");

        const dictionary emptyDict;
        const dictionary& modelDict = sDict.found(coeffsName) 
                                    ? sDict.subDict(coeffsName)
                                    : emptyDict;

        // Construct the model using the RTS
        transportModels_.set
        (
            i,
            plasmaTransportModel::New
            (
                modelName,
                modelDict,
                mesh_,
                species_,
                i
            )
        );

        Info << "    " << sName << ": transport model '" << modelName 
             << "' successfully constructed." << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaTransport::plasmaTransport
(
    const fvMesh& mesh,
    plasmaSpecies& species
)
:
    regIOobject
    (
        IOobject
        (
            "plasmaTransport",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    species_(species),
    transportModels_(species.nSpecies()),
    convectiveFlux_(species.nSpecies()),
    diffusiveFlux_(species.nSpecies()),
    particleFlux_(species.nSpecies()),
    k_eff_
    (
        IOobject
        (
            "k_eff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ),
    alpha_
    (
        IOobject
        (
            "alpha",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE       // switch to AUTO_WRITE if you want it in output for paraview
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0)
    ),
    alphaDx_
    (
        IOobject
        (
            "alphaDx",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),                                          // ← closes alphaDx_ entirely
    S_iz_                                       // ← new sibling entry
    (
        IOobject
        (
            "S_iz",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    )
{
    constructTransportModels();

    particleFlux_.setSize(species.nSpecies());
    convectiveFlux_.setSize(species.nSpecies());
    diffusiveFlux_.setSize(species.nSpecies());

    const dimensionedScalar zeroFlux
    (
        "zero",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        0.0
    );

    for (const label i : species_.mobileSpeciesIDs())
    {
        const word& sName = species_.speciesName(i);

        particleFlux_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "particleFlux_" + sName, 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                zeroFlux
            )
        );

        convectiveFlux_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "convectiveFlux_" + sName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                zeroFlux
            )
        );

        diffusiveFlux_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "diffusiveFlux_" + sName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                zeroFlux
            )
        );
    }

            photoionization_ = photoionizationModel::New(mesh_);
        Info << "Photoionization model loaded: "
             << photoionization_->type() << endl;

    // Mechanism-driven chemistry, if the case asks for it. Reads its own
    // dictionary so no constructor signature changes and existing cases are
    // untouched.
    {
        IOdictionary transportDict
        (
            IOobject
            (
                "plasmaTransportProperties",
                mesh_.time().constant(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );
        readChemistry(transportDict);
    }
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

// Defined here, where janafMixture is complete: an autoPtr to a
// forward-declared type can only be destroyed where the type is known.
plasmaTransport::~plasmaTransport()
{}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void plasmaTransport::correctTransportModels()
{

    // Update transport coefficients for all models
    for (label i = 0; i < species_.nSpecies(); ++i)
    {
        transportModels_[i].correct();
    }
    
    Info << "Transport models corrected." << endl;
}

// This is for the positive streamer case
// * * * * * * * * * * * * * * * G2: gas energy * * * * * * * * * * * * * * //

void Foam::plasmaTransport::solveGasEnergy(const scalar dt)
{
    // `chem_` is deliberately NOT required. The prompt channels -- elastic,
    // rotational and the gas share of the inelastic defect -- come from
    // `rates_`, and only the heavy-reaction heat release needs the chemistry
    // object. Requiring it here made `chemistrySolver none` a SILENT no-op:
    // the run printed "gas heating ON, T solved" at start-up and then left
    // T_gas at exactly 300 K for the whole simulation, with no warning.
    if (!gasHeating_ || !rates_) return;

    // The EQUATION belongs to plasmaEnergy, which owns T_gas. What belongs
    // here is the SOURCE: only this class knows the reaction enthalpies, the
    // EEDF power channels and the mechanism thermodynamics. Splitting it that
    // way keeps the temperature in one place and the physics that heats it in
    // another, which is the division the plasmaEnergy interface already
    // assumed.
    if (!mesh_.foundObject<plasmaEnergy>("plasmaEnergy")) return;
    auto& energy =
        const_cast<plasmaEnergy&>(mesh_.lookupObject<plasmaEnergy>("plasmaEnergy"));
    if (!energy.solvesGasEnergy()) return;

    // Bind the borrowed view once. Done here rather than in the constructor
    // because plasmaEnergy is registered by the solver, which may construct it
    // after this class.
    if (!TgasField_)
    {
        TgasField_ = &energy.TgasField();

        // CHECK ONCE, AND FAIL WITH THE EXACT TEXT TO PASTE. Without this the
        // case dies mid-run on "Entry 'laplacian(kappa,T_gas)' not found",
        // which is accurate and tells the user nothing about what to write.
        //
        // The names are DERIVED FROM THE FIELD rather than written out. An
        // earlier version hard-coded `T`, plasmaEnergy calls its field
        // `T_gas`, and the checks therefore passed while the run died anyway.
        const word& Tn = TgasField_->name();

        const dictionary& ls =
            mesh_.schemes().subOrEmptyDict("laplacianSchemes");
        if (!ls.found("default") && !ls.found("laplacian(kappa," + Tn + ")"))
        {
            FatalErrorInFunction
                << "gas heating is on, but system/fvSchemes has no laplacian"
                << " scheme for the conduction term." << nl
                << "    Add to laplacianSchemes:" << nl << nl
                << "        laplacian(kappa," << Tn
                << ")  Gauss linear corrected;" << nl << nl
                << "    (or set kappa 0 in the energy dictionary to drop"
                << " conduction, which is exact for a 0-D comparison.)"
                << exit(FatalError);
        }

        const dictionary& sd = mesh_.solution().subOrEmptyDict("solvers");
        const word both = "\"(" + Tn + "|" + Tn + "Final)\"";
        if ((!sd.found(Tn) || !sd.found(Tn + "Final")) && !sd.found(both))
        {
            FatalErrorInFunction
                << "gas heating is on, but system/fvSolution has no linear"
                << " solver for " << Tn << "." << nl
                << "    Add to solvers:" << nl << nl
                << "        " << both << nl
                << "        {" << nl
                << "            solver          PCG;" << nl
                << "            preconditioner  DIC;" << nl
                << "            tolerance       1e-12;" << nl
                << "            relTol          0;" << nl
                << "        }" << nl << nl
                << "    A SYMMETRIC solver: the equation is ddt plus laplacian"
                << " with no convection, so DILU and the other asymmetric"
                << " preconditioners are rejected." << nl
                << "    Both " << Tn << " and " << Tn << "Final are needed:"
                << " it is solved on the final outer iteration."
                << exit(FatalError);
        }
    }

    const label nCells = mesh_.nCells();
    const label eIdx = species_.electronSpeciesID();
    const scalarField& ne =
        species_.numberDensities()[eIdx].primitiveField();

    static const scalar EVJ = 1.602176634e-19;

    // Power channels from plasmaReactionRates, re-interpolated every step
    // alongside k_j -- so a dynamic-EEDF rebuild reaches them too.
    const scalarField& Pel = rates_->PelasticN().primitiveField();
    const scalarField& Pgs = rates_->PgasN().primitiveField();
    const scalarField& Pvb = rates_->PvibN().primitiveField();

    volScalarField Qgas
    (
        IOobject
        (
            "Q_gas", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        mesh_, dimensionedScalar(dimensionSet(1, -1, -3, 0, 0), Zero)
    );
    volScalarField Pvib
    (
        IOobject
        (
            "P_vib", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        mesh_, dimensionedScalar(dimensionSet(1, -1, -3, 0, 0), Zero)
    );
    volScalarField rhoCv
    (
        IOobject
        (
            "rhoCv", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        mesh_, dimensionedScalar(dimensionSet(1, -1, -2, -1, 0), Zero)
    );
    volScalarField tauVT
    (
        IOobject
        (
            "tau_VT", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        mesh_, dimensionedScalar(dimTime, GREAT)
    );

    const volScalarField& T = energy.Tgas()();
    scalarField n(species_.nSpecies(), Zero);

    // THE BACKGROUND GAS CARRIES THE HEAT CAPACITY. In the CFD the bulk
    // neutrals are a RESERVOIR, not transported fields, so rho and c_v summed
    // over the transported species alone counted ~1e13 m^-3 of ions where
    // 2.4e25 m^-3 of N2 and O2 sit: rho*c_v came out at 3.5e-10 instead of
    // ~835 J/m^3/K, and half a nanosecond of a 62 Td field raised the gas by
    // 541 K instead of 5e-7 K. The identical code is correct in the 0-D
    // reactor, where every species IS in the state vector -- which is why the
    // two must be compared rather than assumed equivalent.
    const label iO  = species_.speciesNames().find("O");
    const label iN2 = species_.speciesNames().find("N2");
    const label iO2 = species_.speciesNames().find("O2");

    // The same two species as they appear in the BACKGROUND composition, which
    // is where they actually live whenever the case transports only charges.
    const label bgIN2 = bgNames_.find("N2");
    const label bgIO2 = bgNames_.find("O2");

    for (label celli = 0; celli < nCells; ++celli)
    {
        const scalar nEl = max(ne[celli], scalar(0));

        // LIVE heavy density: dissociation raises the particle count, and the
        // tabulated powers are per unit gas density.
        //
        // THE BACKGROUND RESERVOIR BELONGS HERE TOO -- the same trap that made
        // rho*c_v 1e9 too small, in the other factor of the same product. The
        // sweep tabulates power per electron PER UNIT GAS DENSITY, so this must
        // be the density of the whole gas. With `include charged` the only
        // transported heavies are ions, ~7e19 m^-3 against 2.4e25 m^-3 of N2
        // and O2 sitting in the reservoir, so the source came out ~3e5 times
        // too small: a streamer head at 310 Td and n_e = 6e19 heated by 3e-6 K
        // where the tables predict 0.5 K.
        scalar nHeavy = 0.0;
        forAll(n, sp)
        {
            n[sp] = species_.numberDensities()[sp].primitiveField()[celli];
            if (sp != eIdx) nHeavy += max(n[sp], scalar(0));
        }

        // Whatever the reservoir still holds once the transported neutrals are
        // subtracted, so a dissociated N2 is not counted both as a background
        // molecule and as a transported atom. Identical construction to the
        // rho*c_v term below, and for the identical reason.
        scalar nBgLeft = species_.backgroundDensity().value();
        for (const label sp : species_.neutralSpeciesIDs())
        {
            nBgLeft -= max(n[sp], scalar(0));
        }
        nBgLeft = max(nBgLeft, scalar(0));
        nHeavy += nBgLeft;

        // Prompt heat: elastic and rotational, plus the gas share of the
        // inelastic defect. The heavy reactions add fast gas heating on top,
        // from their own enthalpies -- the Popov two-step mechanism, which is
        // where most of the heat in air comes from.
        Qgas.primitiveFieldRef()[celli] =
            (Pel[celli] + Pgs[celli])*nEl*nHeavy*EVJ
          + (chem_ ? chem_->heavyHeatRelease(n, T[celli])*EVJ : 0.0);

        Pvib.primitiveFieldRef()[celli] = Pvb[celli]*nEl*nHeavy*EVJ;

        // Mole fractions for the V-T rate, over the FULL gas again. With
        // `include charged` neither N2 nor O2 is transported, so taking these
        // from the transported set alone handed vibRelax x_N2 = x_O2 = 0 --
        // a V-T time for a gas made of nothing. The reservoir supplies them.
        const scalar tot = max(nHeavy, scalar(1));
        const scalar nO  = (iO  >= 0) ? max(n[iO],  scalar(0)) : 0.0;
        const scalar nN2 = ((iN2 >= 0) ? max(n[iN2], scalar(0)) : 0.0)
                         + (bgIN2 >= 0 ? bgX_[bgIN2]*nBgLeft : 0.0);
        const scalar nO2 = ((iO2 >= 0) ? max(n[iO2], scalar(0)) : 0.0)
                         + (bgIO2 >= 0 ? bgX_[bgIO2]*nBgLeft : 0.0);
        const scalar xN2 = nN2/tot;
        const scalar xO2 = nO2/tot;

        tauVT.primitiveFieldRef()[celli] = (tauVTfixed_ > 0)
            ? tauVTfixed_
            : vibRelax::tauVT_N2(T[celli], pGasAtm_, xN2, xO2, nO);

        // Heat capacity per unit VOLUME from the mechanism's own NASA7 data,
        // over the FULL gas: transported species plus the background reservoir
        // left after they are subtracted. `nBgLeft` is the one computed for
        // nHeavy above -- the two must be the same gas, or the source and the
        // heat capacity describe different mixtures.
        scalar rc = thermo_.valid()
            ? thermo_().rho(n)*thermo_().cv(n, T[celli])
            : 0.0;

        if (bgThermo_.valid() && bgThermo_().valid())
        {
            scalarField nBg(bgX_*nBgLeft);
            rc += bgThermo_().rho(nBg)*bgThermo_().cv(nBg, T[celli]);
        }

        rhoCv.primitiveFieldRef()[celli] =
            (rc > VSMALL) ? rc : 1.16*1005.0/1.4;
    }

    {
        // One-off provenance line: which thermo the heat capacity came from,
        // and what it evaluates to. rho*c_v being wrong is silent in every
        // other diagnostic -- it shows up only as a temperature that is
        // implausible if you happen to estimate the right answer first.
        static bool reported = false;
        if (!reported)
        {
            reported = true;
            Info<< "  gas energy: rho*c_v = "
                << gMin(rhoCv.primitiveField()) << " .. "
                << gMax(rhoCv.primitiveField()) << " J/m^3/K"
                << " (transported thermo "
                << (thermo_.valid() && thermo_().valid() ? "on" : "off")
                << ", background thermo "
                << (bgThermo_.valid() && bgThermo_().valid() ? "on" : "off")
                << ")" << endl;
        }
    }

    Qgas.correctBoundaryConditions();
    Pvib.correctBoundaryConditions();
    rhoCv.correctBoundaryConditions();
    tauVT.correctBoundaryConditions();

    TgasMax_ = max(TgasMax_, gMax(T.primitiveField()));
    QgasMax_ = max(QgasMax_, gMax(Qgas.primitiveField()));

    energy.solveGasEnergy(Qgas, rhoCv, Pvib, tauVT, dt);

    eVibMax_ = max(eVibMax_, gMax(energy.eVib().primitiveField()));
}


void plasmaTransport::solve(const bool finalIter)
{
    // ── 0. Chemistry, first half-step ─────────────────────────────────────
    //
    // Strang splitting: C(dt/2) -> T(dt) -> C(dt/2). Second order in time,
    // where applying the chemistry once either side of the transport is only
    // first order. Rate coefficients are interpolated at the field as it
    // stands, which is what "the field does not move while the chemistry
    // integrates" means in practice.
    // Chemistry as a SOURCE, not as a field update.
    //
    // The first implementation advanced the species fields directly, Strang
    // split around the transport solve. That is wrong in OpenFOAM and the
    // measurement showed it: fvm::ddt(n) differences against n.oldTime(),
    // which is stored before the chemistry runs and does not follow it, so
    // the transport solve saw the chemical change as a transport change and
    // partially undid it. Audited directly -- oldTime kept the pre-chemistry
    // value at every step -- and it capped the coupled scheme at first order
    // while the unsplit path reached p = 1.94.
    //
    // Instead the ODE is integrated from the START-OF-STEP state over dt and
    // the mean rate (n_chem - n^n)/dt is handed to the transport equations as
    // an explicit source, exactly as reactingFoam does with RR. ddt and
    // oldTime then stay consistent by construction.
    //
    // Recomputed on EVERY outer iteration, from n^n rather than from the
    // partially updated field, so it is idempotent rather than cumulative --
    // which is what made the previous version scale with nOuterCorrectors --
    // and the chemistry-transport coupling becomes implicit as the outer loop
    // converges.
    if (chem_)
    {
        plasmaSimulationProfiler::start("Plasma Transport", "chemistry ODE");

        if (mesh_.time().timeIndex() != chemTimeIndex_)
        {
            chemTimeIndex_ = mesh_.time().timeIndex();
            chemN0_.setSize(species_.nSpecies());
            chemExt_.setSize(species_.nSpecies());
            forAll(chemN0_, s)
            {
                chemN0_[s] = species_.numberDensities()[s].primitiveField();
                chemExt_[s].setSize(mesh_.nCells(), Zero);
                chemExt_[s] = Zero;
            }
        }

        rates_->correct();
        if (chemistrySolver_ == csODE)
        {
            computeChemistrySources(mesh_.time().deltaTValue());
        }
        plasmaSimulationProfiler::stop("Plasma Transport", "chemistry ODE");
    }

    // ── 1. Update transport coefficients ──────────────────────────────────────
    plasmaSimulationProfiler::start("Plasma Transport", "correctTransportModels");
    correctTransportModels();
    plasmaSimulationProfiler::stop("Plasma Transport", "correctTransportModels");

    // ── 2. Species and fields ─────────────────────────────────────────────────
    plasmaSimulationProfiler::start("Plasma Transport", "chemistry");

    const label eIdx = species_.electronSpeciesID();

    // First positive ion, by charge rather than by name. The legacy Townsend
    // path needs somewhere to put its ionisation source and the tutorial called
    // that species `pIon`; a real mechanism has no such species. -1 when the
    // case carries no positive ion at all, which only the mechanism path can
    // cope with -- the legacy branch checks before using it.
    label iIdx = -1;
    forAll(species_.speciesChargeNumbers(), i)
    {
        if (species_.speciesChargeNumbers()[i] > 0) { iIdx = i; break; }
    }

    volScalarField& ne = species_.numberDensity(eIdx);

    // Calculate Electric field magnitude
    const volScalarField& Emag = species_.em().Emag();

    // alphaEff field (one allocation, reused)
    volScalarField alphaEff
    (
        IOobject("alphaEff", mesh_.time().timeName(), mesh_,
                IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0)
    );

    // constants as plain scalars — no dimensioned temporaries
    const scalar E_const = 2.73e7;
    const scalar E_pow   = 4.3666e26;

    scalarField& a = alpha_.primitiveFieldRef();
    const scalarField& E = Emag.primitiveField();

    forAll(a, c)
    {
        const scalar Ec  = max(E[c], 1.0);              // safeEmag inline
        const scalar inv = 1.0/Ec;
        a[c] = (1.1944e6 + E_pow*inv*inv*inv)            // pow(.,3) -> mult
             * Foam::exp(-E_const*inv);
        // NOTE: no "- eta" here; α stored plain so AMR can use it
    }
    alpha_.correctBoundaryConditions();
    // Update alpha*Dx for diagnostics / ParaView
    {
        const labelListList& cellPts = mesh().cellPoints();
        const pointField& pts = mesh().points();
        const Vector<label> gd = mesh().geometricD();
        const vector mask
        (
            gd.x() == 1 ? 1.0 : 0.0,
            gd.y() == 1 ? 1.0 : 0.0,
            gd.z() == 1 ? 1.0 : 0.0
        );

        scalarField& ad = alphaDx_.primitiveFieldRef();
        forAll(ad, c)
        {
            const labelList& cp = cellPts[c];
            scalar dMax = 0.0;
            forAll(cp, i)
            {
                for (label j = i + 1; j < cp.size(); ++j)
                {
                    const vector d = cmptMultiply(pts[cp[i]] - pts[cp[j]], mask);
                    dMax = max(dMax, magSqr(d));
                }
            }
            ad[c] = a[c] * Foam::sqrt(dMax);
        }
        alphaDx_.correctBoundaryConditions();
    }

    plasmaSimulationProfiler::stop("Plasma Transport", "chemistry");







    plasmaSimulationProfiler::start("Plasma Transport", "buildEquations");
    // ── 4. Build equations for ALL species ────────────────────────────────
    List<autoPtr<fvScalarMatrix>> eqns(species_.nSpecies());

    for (label i = 0; i < species_.nSpecies(); ++i)
    {
        eqns[i].reset(transportModels_[i].nEqn().ptr());
    }
    plasmaSimulationProfiler::stop("Plasma Transport", "buildEquations");
    

//    const volVectorField driftFlux = fvc::reconstruct(convectiveFlux_[eIdx]);

//    volScalarField explicitSource_ = alphaEff * mag(driftFlux);
    
//    k_eff_ = explicitSource_ / species_.numberDensity(eIdx);
//    k_eff_.correctBoundaryConditions();

    volScalarField explicitSource
    (
        IOobject
        (
            "explicitSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    );

    {
        const scalar mu_ref = 2.398;
        const scalar mu_exp = -0.26;
        const scalar eta    = 340.75;   // moved here

scalarField& src        = explicitSource.primitiveFieldRef();
scalarField& keff       = k_eff_.primitiveFieldRef();
scalarField& siz        = S_iz_.primitiveFieldRef();
const scalarField& aRaw = alpha_.primitiveField();
const scalarField& E    = Emag.primitiveField();
const scalarField& neI  = ne.primitiveField();

forAll(src, c)
{
    const scalar Ec  = max(E[c], 1.0);
    const scalar mu  = mu_ref * Foam::pow(Ec, mu_exp);
    const scalar vd  = mu * Ec;                          // drift speed
    const scalar S_i = aRaw[c] * vd * neI[c];            // α·μ·E·ne  (ionization)
    const scalar S_a = eta     * vd * neI[c];            // η·μ·E·ne  (attachment)

    siz[c]  = S_i;                                       // for photoionization
    src[c]  = S_i - S_a;                                 // net, for continuity
    keff[c] = (aRaw[c] - eta) * vd;                      // net rate, unchanged
}
S_iz_.correctBoundaryConditions();
    }
    k_eff_.correctBoundaryConditions();

    plasmaSimulationProfiler::start("Plasma Transport", "solveEquations");

    // Mechanism-driven sources replace the block below when a mechanism is
    // configured. Returns false with no `chemistry` dictionary, in which case
    // the legacy two-species Townsend fits apply exactly as before.
    if (!mechanismSourceTerms(eqns, ne, Emag))
    {
        // Legacy path: the effective-alpha fit produces one generic positive
        // ion, so it needs one to exist.
        if (iIdx < 0)
        {
            FatalErrorInFunction
                << "The legacy Townsend-fit source model requires a positive"
                << " ion species, and this case carries none." << nl
                << "    Add one to activeSpecies, or configure a `chemistry`"
                << " dictionary to use mechanism-driven sources." << nl
                << exit(FatalError);
        }

        // Solve Continuity Equations
        *eqns[eIdx] -= explicitSource;
        *eqns[iIdx] -= explicitSource;
    }
    
    // Photoionization source (if a photoionizationModel is loaded)
        photoionization_->correct();

        const volScalarField& Sph = photoionization_->Sph();
        *eqns[eIdx] -= Sph;

        // Photoionization of the gas produces a specific ion, not a generic
        // one: in air it is O2 that is ionized by the VUV emitted from excited
        // N2 (Zheleznyak et al. 1982), so the electron it releases is paired
        // with O2+. Falls back to the first positive ion when that species is
        // not carried.
        {
            label phIdx = iIdx;
            forAll(species_.speciesNames(), i)
            {
                if (species_.speciesNames()[i] == photoIonSpecies_)
                {
                    phIdx = i;
                    break;
                }
            }

            if (phIdx >= 0) { *eqns[phIdx] -= Sph; }
        }

    // Every species is solved, not just the electron and one ion. Assembling a
    // source term for a species whose equation is never solved would leave it
    // frozen at its initial value while looking entirely healthy.
    forAll(eqns, i)
    {
        eqns[i]->solve();
    }

    forAll(eqns, i)
    {
        species_.numberDensity(i).correctBoundaryConditions();
    }
    

    // plasmaSimulationProfiler::start("Clamp number densities");
    species_.clampNumberDensities();
    // plasmaSimulationProfiler::stop("Clamp number densities");

    // Transport rate implied by the solution just obtained, for the next outer
    // iteration's chemistry integration:
    //
    //     (n - n^n)/dt = T + RR      =>      T = (n - n^n)/dt - RR
    //
    // Iterating this to convergence is what makes the split second order: at
    // the fixed point the ODE integrates along the cell's true trajectory
    // rather than a closed-cell one, so its mean rate is the mean of the true
    // rate and not merely its value at the start of the step.
    // Report the Picard convergence once the timestep's outer loop is done.
    // Reported whether it passes or fails: a check whose only evidence of
    // success is silence cannot be told from one that never ran.
    if (chem_ && finalIter
     && (chemistrySolver_ == csImplicitRate || chemistrySolver_ == csAdaptive))
    {
        chemPicardPeak_ = max(chemPicardPeak_, chemPicardChange_);

        const label ti = mesh_.time().timeIndex();
        const bool firstReport = !chemPicardReported_;
        const bool duePicard =
            firstReport
         || (chemReportInterval_ > 0
          && ti - chemPicardLastReport_ >= chemReportInterval_);

        if (duePicard)
        {
            chemPicardReported_ = true;
            chemPicardLastReport_ = ti;

            // The PEAK since the last report, not this step's value: one
            // unresolved step is a defect in the solution whether or not the
            // step that happens to be sampled shows it.
            const scalar change = chemPicardPeak_;
            chemPicardPeak_ = 0;

            Info<< "plasmaChemistry: Picard change over the last outer"
                << " iteration: " << change << " (" << chemOuterCount_
                << " iterations, peak since last report)" << endl;

            // A configuration warning, so it is worth saying once and not on
            // every report -- nothing about it changes during the run.
            if (chemOuterCount_ < 2 && firstReport)
            {
                WarningInFunction
                    << "implicitRate with nOuterCorrectors = 1." << nl
                    << "    The source is then evaluated once, explicitly, and"
                    << " the scheme is FIRST order -- the second-order"
                    << " behaviour comes from the outer loop making it"
                    << " implicit. Raise nOuterCorrectors, or accept first"
                    << " order knowingly." << endl;
            }
            else if (change > 1e-2 && chemistrySolver_ == csAdaptive)
            {
                // In adaptive mode the stiff cells are already being
                // integrated, so a residual change here means the LIMIT is set
                // too high for this case rather than that the mode is wrong.
                WarningInFunction
                    << "the chemistry source is still changing by "
                    << change << " on the last outer iteration, with "
                    << chemNstiff_ << " cell(s) already integrated." << nl
                    << "    Lower `chemStiffnessLimit` (currently "
                    << chemStiffLimit_ << ") so more cells take the integrated"
                    << " path, or reduce deltaT." << endl;
            }
            else if (change > 1e-2)
            {
                WarningInFunction
                    << "the chemistry source is still changing by "
                    << change << " on the last outer iteration, so the"
                    << " Picard iteration has NOT converged." << nl
                    << "    max(L*dt) = " << chemStiffness_
                    << ": the chemistry is stiffer than this timestep, and"
                    << " linearising around the current state is not resolving"
                    << " it." << nl
                    << "    Either reduce deltaT, raise nOuterCorrectors, or"
                    << " switch to `chemistrySolver ode`, which integrates"
                    << " through the stiffness at the cost of first order."
                    << endl;
            }
        }
    }

    // csODE only: implicitRate never fills chemRR_, and indexing it here
    // segfaulted on the first timestep.
    if (chem_ && chemistrySolver_ == csODE && !chemN0_.empty() && !chemRR_.empty())
    {
        const scalar rdt = 1.0/mesh_.time().deltaTValue();
        forAll(chemExt_, s)
        {
            const scalarField& nNow = species_.numberDensities()[s].primitiveField();
            forAll(chemExt_[s], c)
            {
                chemExt_[s][c] = (nNow[c] - chemN0_[s][c])*rdt - chemRR_[s][c];
            }
        }
    }



    // ── 5. Update fluxes for mobile species ───────────────────────────────
    // for (const label i : species_.mobileSpeciesIDs())
    // {
    //     transportModels_[i].updateFluxes
    //     (
    //         *eqns[i],
    //         convectiveFlux_[i],
    //         diffusiveFlux_[i],
    //         particleFlux_[i]
    //     );
    // }
    // ── 6. Gas energy ─────────────────────────────────────────────────────
    //
    // On the FINAL outer iteration only. The energy equation is driven by the
    // converged chemistry and field of this timestep, and advancing it once
    // per outer iteration would apply the same heating several times -- the
    // same mistake the chemistry source made before it was reformulated to
    // integrate from the start-of-step state.
    if (finalIter && gasHeating_)
    {
        plasmaSimulationProfiler::start("Plasma Transport", "gas energy");
        solveGasEnergy(mesh_.time().deltaTValue());
        plasmaSimulationProfiler::stop("Plasma Transport", "gas energy");
    }

    plasmaSimulationProfiler::stop("Plasma Transport", "solveEquations");
}

// // Solve for SDBD
// void plasmaTransport::solve()
// {
//     // ── 1. Update transport coefficients ──────────────────────────────────────
//     correctTransportModels();

//     // ── 2. Species and fields ─────────────────────────────────────────────────
//     const label eIdx = species_.electronSpeciesID();
//     const label pIdx = species_.speciesID("pIon");
//     const label nIdx = species_.speciesID("nIon");

//     volScalarField& ne = species_.numberDensity(eIdx);
//     volScalarField& ni = species_.numberDensity(pIdx);
//     volScalarField& nn = species_.numberDensity(nIdx);

//     // ── 3. E/N and rate coefficients ──────────────────────────────────────────
//     const dimensionedScalar N_gas("N_gas", dimless/pow(dimLength,3), 2.4463e25);
//     const scalar N_val = N_gas.value();

//     const volScalarField safeEmag = max
//     (
//         species_.em().Emag(),
//         dimensionedScalar("minE", species_.em().Emag().dimensions(), 1.0)
//     );

//     static interpolationTable<scalar> tableAlpha, tableKatt, tableKei;
//     static bool tablesLoaded = false;
//     if (!tablesLoaded)
//     {
//         tableAlpha = interpolationTable<scalar>
//             (mesh_.time().constant()/"totalIonizationReducedTownsendCoeffs");
//         tableKatt  = interpolationTable<scalar>
//             (mesh_.time().constant()/"totalAttachmentRate");
//         tableKei   = interpolationTable<scalar>
//             (mesh_.time().constant()/"totalIonElectronRecombinationRate");
//         tablesLoaded = true;
//     }

//     volScalarField alpha
//     (
//         IOobject("alpha", mesh_.time().timeName(), mesh_,
//                  IOobject::NO_READ, IOobject::AUTO_WRITE),
//         mesh_,
//         dimensionedScalar("zero", dimensionSet(0,-1,0,0,0,0,0), 0.0)
//     );
//     volScalarField k_att
//     (
//         IOobject("k_att", mesh_.time().timeName(), mesh_,
//                  IOobject::NO_READ, IOobject::AUTO_WRITE),
//         mesh_,
//         dimensionedScalar("zero", dimensionSet(0,3,-1,0,0,0,0), 0.0)
//     );
//     volScalarField k_ei(k_att);
//     const dimensionedScalar k_ii("k_ii", dimensionSet(0,3,-1,0,0,0,0), 1.7e-12);

//     forAll(ne, cellI)
//     {
//         const scalar enKey = safeEmag[cellI] / N_val;

//         auto clamp = [](scalar x, const interpolationTable<scalar>& t) -> scalar
//         {
//             return max(t.first().first(), min(t.last().first(), x));
//         };

//         alpha[cellI]  = tableAlpha(clamp(enKey, tableAlpha)) * N_val;
//         k_att[cellI]  = tableKatt(clamp(enKey, tableKatt));
//         k_ei[cellI]   = tableKei(clamp(enKey, tableKei));
//     }

//     // ── 4. Build equations for all mobile species ─────────────────────────────
//     // Immobile species have no equation — mobileSpeciesIDs() skips them
//     List<autoPtr<fvScalarMatrix>> eqns(species_.nSpecies());
//     for (const label i : species_.mobileSpeciesIDs())
//         eqns[i].reset(transportModels_[i].nEqn().ptr());

//     // ── 5. Fill old-n fluxes — needed by chemistry below ─────────────────────
//     // updateFluxes is virtual on plasmaTransportModel — no driftDiffusion cast
//     // Called BEFORE solve: fvMatrix::flux() uses current (old) n
//     for (const label i : species_.mobileSpeciesIDs())
//     {
//         transportModels_[i].updateFluxes
//         (
//             *eqns[i],
//             convectiveFlux_[i],
//             diffusiveFlux_[i],
//             particleFlux_[i]
//         );
//     }
//     // ── 6. Townsend ionization from old-n electron fluxes ────────────────────
//     // Use plasmaTransport arrays directly — no driftDiffusion cast

//     const volVectorField driftVec    = fvc::reconstruct(convectiveFlux_[eIdx]);
//     const volVectorField particleVec = fvc::reconstruct(particleFlux_[eIdx]);

//     const dimensionedScalar smallFlux
//         ("small", driftVec.dimensions(), 1e-6);
//     const dimensionedScalar zeroFlux
//         ("zero",  driftVec.dimensions(), 0.0);

//     const volVectorField driftDir    = driftVec / (mag(driftVec) + smallFlux);

//     const volScalarField ionizationFlux = min
//     (
//         max(particleVec & driftDir, zeroFlux),
//         mag(driftVec)
//     );

//     const volScalarField S_iz = alpha * ionizationFlux;

//     // ── 7. Add chemistry to equations (reactive species only) ─────────────────
//     *eqns[eIdx] -= S_iz;
//     *eqns[eIdx] += (k_att * N_gas + k_ei * ni) * ne;

//     *eqns[pIdx] -= S_iz;
//     *eqns[pIdx] += (k_ei * ne + k_ii * nn) * ni;

//     *eqns[nIdx] -= k_att * N_gas * ne;
//     *eqns[nIdx] += k_ii * ni * nn;

//     // ── 8. Solve all mobile species ───────────────────────────────────────────
//     for (const label i : species_.mobileSpeciesIDs())
//         eqns[i]->solve();

//     // ── 9. Clamp densities ────────────────────────────────────────────────────
//     species_.clampNumberDensities();

//     for (const label i : species_.mobileSpeciesIDs())
//     {
//         transportModels_[i].updateFluxes
//         (
//             *eqns[i],
//             convectiveFlux_[i],
//             diffusiveFlux_[i],
//             particleFlux_[i]
//         );
//     }
// }

void plasmaTransport::updateSurfaceCharge()
{
    const scalar dt = mesh_.time().deltaTValue();

    volScalarField& sigma = species_.em().surfCharge();

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];

        if (p.coupled()) continue;

        scalarField& sigmaPatch = sigma.boundaryFieldRef()[patchi];

        // Always reset from previous time step value
        sigmaPatch = sigma.oldTime().boundaryField()[patchi];

        const scalarField& magSf = mesh_.magSf().boundaryField()[patchi];

        for (const label i : species_.mobileSpeciesIDs())
        {
            const fvPatchField<scalar>& pField =
                species_.numberDensity(i).boundaryField()[patchi];

            // Use plasmaWallBC interface 
            const plasmaWallBC* pBC =
                dynamic_cast<const plasmaWallBC*>(&pField);

            if (!pBC || !pBC->enableSurfaceCharging()) continue;

            if (!particleFlux_.set(i)) continue;

            sigmaPatch += species_.speciesCharge(i).value()
                       * particleFlux_[i].boundaryField()[patchi]
                       * dt / magSf;
        }

        sigma.correctBoundaryConditions();
    }

    Info << "Surface charge updated." << endl;
}

tmp<volScalarField> plasmaTransport::electricalConductivity() const
{
    tmp<volScalarField> tSigma
    (
        new volScalarField
        (
            IOobject
            (
                "electricalConductivity",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero", 
                dimensionSet(-1, -3, 3, 0, 0, 2, 0), 
                0.0
            )
        )
    );

    // Get reference
    volScalarField& sigma = tSigma.ref();

    for (const label i : species_.mobileSpeciesIDs())
    {
        sigma.ref() += transportModels_[i].electricalConductivity();
    }

    return tSigma;
}

tmp<volScalarField> plasmaTransport::diffusiveChargeSource() const
{
    tmp<volScalarField> tRhoDiff
    (
        new volScalarField
        (
            IOobject
            (
                "diffusiveChargeSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero", 
                dimensionSet(0, -3, 0, 0, 0, 1, 0), 
                0.0
            )
        )
    );

    // Get reference
    volScalarField& rhoDiff = tRhoDiff.ref();

    for (const label i : species_.mobileSpeciesIDs())
    {
        rhoDiff.ref() += transportModels_[i].diffusiveChargeSource();
    }

    return tRhoDiff;
}

bool plasmaTransport::writeData(Ostream& os) const
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * *  Mechanism-driven chemistry  * * * * * * * * * * * * //

void Foam::plasmaTransport::readChemistry(const dictionary& dict)
{
    if (!dict.found("chemistry"))
    {
        Info<< "plasmaTransport: no `chemistry` dictionary; using the legacy"
            << " hard-coded Townsend fits" << endl;
        return;
    }

    const dictionary& cd = dict.subDict("chemistry");

    const word model = cd.getOrDefault<word>("sourceModel", "reactions");
    if (model == "reactions")      sourceModel_ = smReactions;
    else if (model == "townsend")  sourceModel_ = smTownsend;
    else
    {
        FatalErrorInFunction
            << "Unknown sourceModel '" << model << "'" << nl
            << "Valid: reactions | townsend" << nl << exit(FatalError);
    }

    allowChargeNonConservation_ =
        cd.getOrDefault<Switch>("allowChargeNonConservation", false);

    photoIonSpecies_ = cd.getOrDefault<word>("photoIonSpecies", "O2p");

    // Operator-split chemistry.
    //
    //   none  (default) electron-impact sources are added to the transport
    //         equations, as before. Existing cases are unaffected.
    //   ode   the WHOLE mechanism -- electron-impact and heavy -- is
    //         integrated per cell by a stiff ODE solver, Strang-split around
    //         the transport solve. This is the only mode in which ions have a
    //         loss channel, because recombination lives in the heavy set.
    //
    // The two are mutually exclusive by construction: with `ode`,
    // mechanismSourceTerms() adds nothing, or the chemistry would be applied
    // twice.
    // DEFAULT: implicitRate. Second order, exact charge conservation,
    // guaranteed positivity, and far cheaper than the stiff substep -- one
    // O(n_reactions) rate evaluation per cell against a Jacobian and an LU
    // factorisation. `ode` remains for mechanisms whose fastest timescale is
    // far below dt, where linearising around the current state is not enough;
    // the run reports which situation it is in (see the convergence check
    // below) rather than leaving it to judgement.
    // DEFAULT: adaptive -- per cell, the accurate path where it is valid and
    // the robust one where it is not, so neither accuracy nor robustness is
    // traded for the other across a domain where the stiffness varies by
    // orders of magnitude.
    // ---- WHICH reactions (physics) and HOW they are integrated (numerics) --
    //
    // Two keys, because they are two decisions. The old single key
    // `chemistrySolver` is rejected outright rather than mapped onto either
    // axis: its value `none` meant "electron-impact reactions, integrated
    // explicitly", so a case written against it is making BOTH choices at
    // once and only the case author knows which half was intended.
    if (cd.found("chemistrySolver"))
    {
        FatalErrorInFunction
            << "`chemistrySolver` has been split into two keys, because it was"
            << " deciding two" << nl
            << "    different things: which reactions exist (physics) and how"
            << " they are integrated" << nl
            << "    (numerics)." << nl << nl
            << "    chemistry" << nl
            << "    {" << nl
            << "        reactions   all;       // all | electronImpact | none"
            << nl
            << "        solver      adaptive;  // adaptive | implicitRate |"
            << " ode | explicitSource" << nl
            << "    }" << nl << nl
            << "    Translating the old values:" << nl
            << "        chemistrySolver none      ->  reactions electronImpact;"
            << "  solver explicitSource;" << nl
            << "        chemistrySolver adaptive  ->  reactions all;"
            << "              solver adaptive;" << nl
            << "        chemistrySolver ode       ->  reactions all;"
            << "              solver ode;" << nl
            << "        (no old equivalent)       ->  reactions none;"
            << "             // no chemistry at all" << nl << nl
            << "    `none` in particular did NOT mean no chemistry: it solved"
            << " electron-impact" << nl
            << "    reactions and skipped the heavy set, so ions never decayed."
            << nl
            << exit(FatalError);
    }

    advisoryD_ = cd.getOrDefault<scalar>("advisoryDiffusivity", 2.0e-5);

    const word rxs = cd.getOrDefault<word>("reactions", "all");
    if      (rxs == "all")            chemReactions_ = rxAll;
    else if (rxs == "electronImpact") chemReactions_ = rxElectronImpact;
    else if (rxs == "none")           chemReactions_ = rxNone;
    else
    {
        FatalErrorInFunction
            << "Unknown chemistry/reactions '" << rxs << "'" << nl
            << "Valid: all | electronImpact | none" << nl
            << exit(FatalError);
    }

    const word cs = cd.getOrDefault<word>("solver", "adaptive");
    if      (cs == "adaptive")       chemistrySolver_ = csAdaptive;
    else if (cs == "implicitRate")   chemistrySolver_ = csImplicitRate;
    else if (cs == "ode")            chemistrySolver_ = csODE;
    else if (cs == "explicitSource") chemistrySolver_ = csExplicitSource;
    else
    {
        FatalErrorInFunction
            << "Unknown chemistry/solver '" << cs << "'" << nl
            << "Valid: adaptive | implicitRate | ode | explicitSource" << nl
            << exit(FatalError);
    }

    // The legacy explicit path loops over the TABULATED reactions only, which
    // are the electron-impact ones. It has no heavy branch at all, so pairing
    // it with `all` would silently drop recombination -- the same class of
    // silent omission this split exists to end.
    if (chemistrySolver_ == csExplicitSource && chemReactions_ == rxAll)
    {
        FatalErrorInFunction
            << "chemistry/solver `explicitSource` cannot integrate the heavy"
            << " reactions." << nl << nl
            << "    It applies tabulated electron-impact rates directly to the"
            << " transport equations" << nl
            << "    and has no heavy-reaction path, so `reactions all` would"
            << " silently drop" << nl
            << "    recombination, ion-ion neutralisation and detachment."
            << nl << nl
            << "    Use `reactions electronImpact` to keep the legacy"
            << " behaviour, or pick a solver" << nl
            << "    that handles the full set: adaptive | implicitRate | ode."
            << nl
            << exit(FatalError);
    }

    // `implicitRate` evaluates the rate at the current iterate, which only
    // becomes the t^{n+1} rate the scheme asks for once Picard has converged.
    // With one outer iteration there is no loop to converge: on the validated
    // streamer it overshot the saturated electron density before 1.5 ns and
    // then diverged. `adaptive` inherits the same linearisation.
    if
    (
        (chemistrySolver_ == csImplicitRate || chemistrySolver_ == csAdaptive)
     && chemReactions_ != rxNone
     && mesh_.solution().subOrEmptyDict("PIMPLE")
            .getOrDefault<label>("nOuterCorrectors", 1) < 2
    )
    {
        FatalErrorInFunction
            << "chemistry/solver `" << cs << "` needs nOuterCorrectors > 1."
            << nl << nl
            << "    It linearises the chemistry about the current iterate and"
            << " reaches the order it" << nl
            << "    claims only once the outer loop converges. With one outer"
            << " iteration there is no" << nl
            << "    loop to converge, and the source is evaluated at the"
            << " old state." << nl << nl
            << "    Either raise nOuterCorrectors in system/fvSolution, or"
            << " choose a solver that" << nl
            << "    does not linearise: `ode` integrates through the step,"
            << " `explicitSource` is the" << nl
            << "    legacy explicit source (with `reactions electronImpact`)."
            << nl
            << exit(FatalError);
    }

    // WHAT ORDER THIS COMBINATION CAN ACTUALLY REACH.
    //
    // Converging the outer loop is NECESSARY for second order but not
    // SUFFICIENT: `ode` splits the chemistry off and takes a mean rate over
    // the step, which is a Lie composition and first order by construction, so
    // it stays first order however well the loop converges. Measured on the
    // 40x40 benchmark (200 ps, dt 4:2:1, Richardson on n_e max and n_e sum):
    //
    //   lagged     + ode        p = 0.99
    //   converged  + ode        p = 0.99   <- converging did NOT help
    //   converged  + adaptive   p = 1.94
    //
    // Order is achieved, not selected, so the run says which it will get
    // rather than leaving the user to infer it from two separate keys.
    {
        IOdictionary controls
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
        const word target = controls.subOrEmptyDict("outerCoupling")
            .getOrDefault<word>("target", "converged");

        const bool secondOrderCapable =
            (chemistrySolver_ != csODE) && (chemReactions_ != rxNone);

        if (target != "converged")
        {
            Info<< "plasmaTransport: temporal order 1 -- outerCoupling target"
                << " is `" << target << "`, so the" << nl
                << "    Poisson-species coupling is lagged. `ddtSchemes"
                << " backward` does not repair that." << endl;
        }
        else if (!secondOrderCapable)
        {
            Info<< "plasmaTransport: temporal order 1 -- the outer loop is"
                << " converged, but chemistry/solver" << nl
                << "    `" << cs << "` splits the chemistry off with a mean"
                << " rate over the step, which is" << nl
                << "    first order by construction. Use `implicitRate` or"
                << " `adaptive` for second order." << endl;
        }
        else
        {
            Info<< "plasmaTransport: temporal order 2 -- outer loop converged"
                << " and chemistry/solver `" << cs << "`" << nl
                << "    is second-order capable (measured p = 1.94)." << endl;
        }
    }

    // The chemistry object is what integrates a reaction set. It is not needed
    // when there are no reactions, nor for the legacy explicit path, which
    // applies tabulated rates straight into the transport equations.
    if (chemReactions_ != rxNone && chemistrySolver_ != csExplicitSource)
    {
        // Skipping quiescent cells is not an optimisation to add later: a
        // stiff ODE in every cell of a 1.1M-cell mesh, where all the chemistry
        // happens in a streamer head a few hundred cells across, is the
        // difference between minutes and hours per timestep.
        // Defaults to ZERO: every cell is integrated unless the case says
        // otherwise. Correctness first, performance opt-in.
        //
        // A non-zero default is a trap. Set at 1e14 it silently excluded every
        // cell of a case whose background is 1e13, and the run then did no
        // chemistry at all while reporting nothing wrong -- the species simply
        // did not evolve. A threshold has to be chosen against the case's own
        // densities, which only the user knows.
        chemActivityThreshold_ =
            cd.getOrDefault<scalar>("chemActivityThreshold", 0.0);

        // L*dt at which `adaptive` stops linearising a cell and integrates it.
        //
        // L_s is the LOSS COEFFICIENT of species s -- dn_s/dt = P_s - L_s n_s
        // -- so 1/L_s is its loss timescale and L*dt is the number of loss
        // timescales crossed in one step:
        //
        //   L*dt << 1   the species barely decays in a step; linearise
        //   L*dt ~  1   one loss timescale per step; either path works
        //   L*dt >> 1   consumed many times over within a step; integrate,
        //               because linearising about a state the cell is leaving
        //               is extrapolation
        //
        // 1.0 is where linearising stops being obviously safe -- NOT a
        // stability boundary. The implicit sink guarantees positivity at any
        // L*dt; what degrades above 1 is the accuracy of linearising a
        // nonlinear system about a state that is being left behind.
        //
        // Raising it buys speed by linearising more cells, and the cost is
        // silent: it does not crash, it produces a less accurate source in
        // exactly the cells that are reacting. The Picard-change report is
        // what makes that visible, so raise it only with that number in view.
        chemStiffLimit_ =
            cd.getOrDefault<scalar>("chemStiffnessLimit", 1.0);

        // See plasmaTransport.H. Measured rather than assumed: in the 0-D
        // sweep 0.5 let a single accepted step leave an afterglow 190x off,
        // while 0.1 and below reproduced the stiff solver exactly.
        chemChangeLimit_ =
            cd.getOrDefault<scalar>("chemChangeLimit", 0.1);

        // How often the chemistry diagnostics are printed, in timesteps.
        //
        // NOT once at start-up, which is what this used to do and which is the
        // worst possible sampling point: at t=0 the densities are at their
        // floor and the chemistry is at its least stiff. A discharge becomes
        // stiff when it ignites, which is after that single sample is taken.
        //
        // Between reports the values are peak-held, so a stiff step is not
        // missed by falling between two samples. 0 restores start-up-only.
        chemReportInterval_ =
            cd.getOrDefault<label>("chemReportInterval", 200);

        dictionary ccfg(cd);
        ccfg.add("electronName", species_.speciesNames()[species_.electronSpeciesID()], true);
        ccfg.add("backgroundDensity", species_.backgroundDensity().value(), true);

        // The reaction-set choice is applied where the mechanism is READ, so
        // the ODE, the Jacobian, the heat release and every diagnostic all see
        // the same set. Filtering later would leave each consumer to remember
        // the same exclusion independently.
        ccfg.add
        (
            "includeHeavyReactions",
            bool(chemReactions_ == rxAll),
            true
        );

        scalarField q(species_.nSpecies());
        forAll(q, i) q[i] = species_.speciesChargeNumbers()[i];

        chem_.reset
        (
            new plasmaChemistry
            (
                cd.get<fileName>("mechanism"),
                species_.speciesNames(),
                q,
                ccfg
            )
        );

        // How many outer iterations a timestep takes, so the trailing
        // half-step lands on the last one.
        nOuterCorrectors_ = mesh_.solutionDict().subOrEmptyDict("PIMPLE")
            .getOrDefault<label>("nOuterCorrectors", 1);

        Info<< "plasmaTransport: chemistrySolver " << cs
            << ", activity threshold " << chemActivityThreshold_ << " 1/m3, "
            << nOuterCorrectors_ << " outer iteration(s) per step" << nl
            << "    The WHOLE mechanism is evaluated here -- electron-impact"
            << " and heavy -- so the legacy in-equation sources are not"
            << " applied; using both would count the chemistry twice." << nl;

        if (chemistrySolver_ == csImplicitRate)
        {
            Info<< "    implicitRate: instantaneous rate at the current"
                << " iterate, loss implicit. Second order once the outer loop"
                << " converges; the run reports whether it does." << endl;
        }
        else
        {
            Info<< "    ode: stiff substep, mean rate over the step. First"
                << " order, but it integrates THROUGH stiffness -- use it when"
                << " a chemical timescale is far below deltaT." << endl;
        }
    }

    rates_.reset
    (
        new plasmaReactionRates
        (
            mesh_,
            cd.get<fileName>("mechanism"),
            cd.get<fileName>("tableDir"),
            cd.getOrDefault<word>("lookupVariable", "reducedE"),
            cd.getOrDefault<word>("tableKey", "reducedE"),
            cd
        )
    );

    // ---- G2: gas energy ---------------------------------------------------
    tableDir_ = cd.get<fileName>("tableDir");
    {
        const dictionary& ed = species_.backgroundDict().subOrEmptyDict("energy");
        TgasConst_ = ed.getOrDefault<scalar>("T", 300.0);
        gasHeating_ = ed.getOrDefault<bool>("solve", false);   // existing key, see plasmaSpecies
        kappaGas_ = ed.getOrDefault<scalar>("kappa", 0.026);
        tauVTfixed_ = ed.getOrDefault<scalar>("tauVT", -1);
        pGasAtm_ = ed.getOrDefault<scalar>("pressure", 101325.0)/101325.0;

        if (gasHeating_)
        {
            // The FIELDS live in plasmaEnergy, which owns the temperature.
            // Only the mechanism thermodynamics are built here, because only
            // this class reads the mechanism.
            IFstream mis(cd.get<fileName>("mechanism"));
            dictionary mdict(mis);
            thermo_.reset
            (
                new janafMixture
                (
                    mdict,
                    species_.speciesNames(),
                    species_.speciesNames()[species_.electronSpeciesID()]
                )
            );
            // Background gas: whatever the mechanism's reference composition
            // names. These species are NOT transported and have no index in
            // plasmaSpecies, which is exactly why they need their own thermo
            // object -- and why omitting them made rho*c_v 1e9 too small.
            {
                IFstream cis(cd.get<fileName>("mechanism"));
                dictionary cdict(cis);
                const dictionary comp = cdict.subOrEmptyDict("composition");

                DynamicList<word> nm;
                DynamicList<scalar> xs;
                scalar xSum = 0;
                for (const entry& e : comp)
                {
                    nm.append(e.keyword());
                    xs.append(readScalar(e.stream()));
                    xSum += xs.last();
                }
                if (xSum > VSMALL)
                {
                    forAll(xs, i) xs[i] /= xSum;
                    bgNames_ = nm;
                    bgX_ = scalarField(xs);
                    bgThermo_.reset
                    (
                        new janafMixture(cdict, bgNames_, word::null)
                    );
                }
            }

            if (!thermo_().valid())
            {
                WarningInFunction
                    << "gasHeating is on but the mechanism carries no"
                    << " `speciesThermo` block, so the heat capacity falls"
                    << " back to a constant for AIR." << nl
                    << "    Recompile the mechanism with mechc." << endl;
            }

            Info<< "plasmaTransport: gas heating ON, T solved"
                << " (kappa = " << kappaGas_ << " W/m/K, thermo "
                << (thermo_().valid() ? "from mechanism" : "FALLBACK")
                << ")" << endl;

            // Say so when the heavy channel is absent. Without a chemistry
            // object there is no Popov two-step fast heating -- most of the
            // heat in air -- so the temperature is a LOWER BOUND. The prompt
            // channels still run, which is why this is a note and not an
            // error, but a temperature that omits the dominant term in air
            // should say so where the user reads it.
            if
            (
                chemReactions_ != rxAll
            )
            {
                Info<< "    NOTE chemistry/reactions is `" << rxs << "`, so"
                    << " heavy-reaction heat release is NOT included." << nl
                    << "    Only the prompt electron channels (elastic,"
                    << " rotational, inelastic gas share) heat the gas," << nl
                    << "    so T_gas is a LOWER BOUND: in air the two-step"
                    << " quenching of N2 electronic" << nl
                    << "    states by O2 is typically the larger term."
                    << endl;
            }
        }
    }

    {
        IFstream is(cd.get<fileName>("mechanism"));
        dictionary mech(is);
        if (mech.found("composition"))
        {
            mechComposition_ = mech.subDict("composition");
        }
    }

    Info<< "plasmaTransport: sourceModel " << model
        << ", " << rates_->size() << " electron-impact reactions" << endl;

    // The caveat below is TRUE ONLY OF THE LEGACY PATH. It used to print
    // unconditionally, which told a user running `adaptive` -- a mode that
    // does evaluate the heavy set -- that their ions had no loss channel.
    if (chemReactions_ == rxElectronImpact)
    {
        Info<< "    NOTE this path evaluates ELECTRON-IMPACT reactions only."
            << " Heavy chemistry (recombination," << nl
            << "    quenching, detachment) is NOT applied, so ion and"
            << " excited-state densities have NO" << nl
            << "    loss channel. Adequate for a short pulse; not for long"
            << " transients." << endl;
    }
    else if (chemReactions_ == rxNone)
    {
        Info<< "    NOTE chemistry/reactions is `none`: NO chemistry source"
            << " of any kind is applied." << nl
            << "    Drift-diffusion and Poisson only. Nothing ionises."
            << endl;
    }
}




bool Foam::plasmaTransport::finalOuterIteration()
{
    ++chemOuterCount_;
    return chemOuterCount_ >= nOuterCorrectors_;
}


void Foam::plasmaTransport::computeChemistrySources(const scalar dt)
{
    if (!chem_ || dt <= 0) return;

    chemRR_.setSize(species_.nSpecies());
    forAll(chemRR_, s)
    {
        chemRR_[s].setSize(mesh_.nCells(), Zero);
        chemRR_[s] = Zero;
    }


    const label nSp = species_.nSpecies();
    const label eIdx = species_.electronSpeciesID();
    const scalarField& neI = species_.numberDensity(eIdx).primitiveField();

    // Rate coefficients are interpolated per cell from the tables, then held
    // fixed for the substep: under operator splitting the field does not move
    // while the chemistry integrates.
    const label nTab = chem_->nTabulated();
    scalarField kTab(max(nTab, 1), Zero);
    scalarField n(nSp, Zero);
    scalarField ext(nSp, Zero);

    // Per CELL now, not a dictionary constant: this is the thermal-ionisation
    // feedback. See TgasCell().

    label nActive = 0;
    forAll(neI, celli)
    {
        if (neI[celli] < chemActivityThreshold_) continue;
        ++nActive;

        // FROM THE START-OF-STEP STATE, not from the partially updated
        // field: that is what makes recomputation idempotent.
        for (label s = 0; s < nSp; ++s)
        {
            n[s] = chemN0_[s][celli];
        }
        for (label j = 0; j < nTab; ++j)
        {
            kTab[j] = rates_->k(j).primitiveField()[celli];
        }

        // Charge conservation of the RHS, on the cell that is doing the most
        // chemistry. Switching to the split path bypassed the source-term
        // diagnostic, and dropping a conservation check when changing how
        // conservation is achieved is precisely the wrong moment to do it.
        if (!chemCellsReported_)
        {
            chemResidualMax_ =
                max(chemResidualMax_, chem_->chargeResidual(n, kTab, TgasCell(celli)));
        }

        // Transport rate for this cell, from the previous outer iteration.
        const bool haveExt = !chemExt_.empty() && !chemExt_[0].empty();
        if (haveExt)
        {
            for (label s = 0; s < nSp; ++s) ext[s] = chemExt_[s][celli];
        }

        chem_->integrate(n, kTab, TgasCell(celli), dt, haveExt ? &ext : nullptr);

        // The integrated change contains BOTH processes, so the transport part
        // is removed to leave the chemical source alone. Adding the whole
        // change would count transport twice -- once here and once in the
        // equation that is about to be solved.
        for (label s = 0; s < nSp; ++s)
        {
            chemRR_[s][celli] =
                (n[s] - chemN0_[s][celli])/dt - (haveExt ? ext[s] : 0.0);
        }

        // Project the source onto exact charge conservation.
        //
        // The reactions balance charge identically -- the right-hand-side
        // residual is 1e-16 -- but the INTEGRATED result does not, because the
        // ODE solver's own arithmetic does not preserve the linear invariant
        // exactly: measured 3.5e-06 with seulex and 5.0e-08 with rodas23. In a
        // plasma solver that residual is not cosmetic: net charge drives the
        // Poisson equation, so a per-step leak accumulates into a space-charge
        // field that nothing physical produced.
        //
        // The correction is distributed over the charged species in proportion
        // to their own contribution, so it cannot single one out, and it is
        // far below the ODE tolerance it is correcting. It CANNOT mask an
        // unbalanced mechanism: that would show in the right-hand-side
        // residual, which is checked separately and independently.
        {
            scalar netQ = 0, traffic = 0;
            for (label s = 0; s < nSp; ++s)
            {
                const scalar q = species_.speciesChargeNumbers()[s];
                netQ    += q*chemRR_[s][celli];
                traffic += mag(q*chemRR_[s][celli]);
            }

            if (traffic > VSMALL)
            {
                const scalar rel = mag(netQ)/traffic;
                chemSourceChargeMax_ = max(chemSourceChargeMax_, rel);

                // A large residual is structural, not roundoff. Correcting it
                // silently would turn a mechanism error into a plausible
                // answer, so it is reported and left uncorrected.
                if (rel < 1e-3)
                {
                    for (label s = 0; s < nSp; ++s)
                    {
                        const scalar q = species_.speciesChargeNumbers()[s];
                        if (mag(q) < SMALL) continue;
                        chemRR_[s][celli] -=
                            (netQ/traffic)*mag(chemRR_[s][celli])*sign(q);
                    }
                }
                else
                {
                    ++chemUnprojectedCells_;
                }

                // Measured AFTER the projection, so the claim that it is
                // exact is checked rather than asserted.
                if (!chemCellsReported_)
                {
                    scalar netAfter = 0, trafficAfter = 0;
                    for (label s = 0; s < nSp; ++s)
                    {
                        const scalar q = species_.speciesChargeNumbers()[s];
                        netAfter     += q*chemRR_[s][celli];
                        trafficAfter += mag(q*chemRR_[s][celli]);
                    }
                    if (trafficAfter > VSMALL)
                    {
                        chemSourceChargeAfter_ = max(chemSourceChargeAfter_,
                                                     mag(netAfter)/trafficAfter);
                    }
                }
            }
        }

        // ---- consistency checks on the SOURCE, not on the right-hand side --
        //
        // The RHS residual (chargeResidual) says the reactions balance. It
        // says nothing about whether the INTEGRATED source does: the ODE
        // solver could in principle break the linear invariant. It does not --
        // seulex and the Rosenbrock family form each stage as a linear
        // combination of f evaluations, and sum(q.f) = 0 holds for every one
        // of them -- but that is a property worth measuring rather than
        // asserting from the method's name.
        if (!chemCellsReported_)
        {
            // Positivity: the source is explicit, so a species whose loss over
            // dt exceeds its own density would be driven negative. That is not
            // a stability nuisance to be clipped away -- it means dt is longer
            // than the chemistry it is trying to represent, and the mean rate
            // is then not a meaningful description of the step.
            for (label s = 0; s < nSp; ++s)
            {
                if (chemN0_[s][celli] + chemRR_[s][celli]*dt < 0)
                {
                    ++chemNegativeCells_;
                    break;
                }
            }
        }
    }

    if (!chemCellsReported_)
    {
        chemCellsReported_ = true;
        reduce(nActive, sumOp<label>());
        label nTot = neI.size();
        reduce(nTot, sumOp<label>());

        if (nActive == 0)
        {
            // Not a performance note: with no active cells the run evolves no
            // chemistry whatsoever, and every species simply stays where it
            // started. That looks like a converged solution.
            WarningInFunction
                << "chemActivityThreshold = " << chemActivityThreshold_
                << " 1/m3 excludes EVERY cell, so no chemistry is being"
                << " integrated at all." << nl
                << "    The largest electron density is "
                << gMax(neI) << " 1/m3. Lower the threshold, or remove it."
                << endl;
        }
        else
        {
            reduce(chemResidualMax_, maxOp<scalar>());
            reduce(chemSourceChargeMax_, maxOp<scalar>());
            reduce(chemNegativeCells_, sumOp<label>());
            reduce(chemUnprojectedCells_, sumOp<label>());
            reduce(chemSourceChargeAfter_, maxOp<scalar>());
            Info<< "plasmaChemistry: integrating " << nActive << " of " << nTot
                << " cells (" << 100.0*nActive/max(nTot, 1) << "%) above the"
                << " activity threshold" << nl
                << "    charge residual of the RHS, worst cell:    "
                << chemResidualMax_ << nl
                << "    charge residual of the SOURCE, worst cell: "
                << chemSourceChargeMax_
                << "  (sum q_s RR_s / sum |q_s RR_s|, BEFORE projection;"
                << ")" << nl
                << "    charge residual of the SOURCE, after projection:  "
                << chemSourceChargeAfter_ << endl;

            if (chemUnprojectedCells_ > 0)
            {
                WarningInFunction
                    << chemUnprojectedCells_ << " cell(s) had a charge"
                    << " residual above 1e-3 and were left UNPROJECTED." << nl
                    << "    That is too large to be solver arithmetic and"
                    << " points at the mechanism, so it is reported rather"
                    << " than corrected away." << endl;
            }

            if (chemNegativeCells_ > 0)
            {
                WarningInFunction
                    << chemNegativeCells_ << " cell(s) would be driven"
                    << " negative by the explicit chemistry source over this"
                    << " timestep." << nl
                    << "    The source is a MEAN rate over dt, so this means"
                    << " dt is longer than the chemistry it represents and the"
                    << " mean is not a meaningful description of the step."
                    << " Reduce deltaT." << endl;
            }
        }
    }
}


bool Foam::plasmaTransport::mechanismSourceTerms
(
    List<autoPtr<fvScalarMatrix>>& eqns,
    const volScalarField& ne,
    const volScalarField& Emag
)
{
    if (!rates_)
    {
        return false;   // caller falls back to the legacy fits
    }

    // `off` means OFF: no chemistry source of any kind reaches the species
    // equations, leaving drift-diffusion and Poisson alone. Returning TRUE is
    // the whole point -- false would send the caller into its legacy
    // Townsend-fit branch and quietly reinstate ionisation.
    if (chemReactions_ == rxNone)
    {
        return true;
    }

    // With operator-split chemistry the sources are applied by
    // applyChemistry(), not here. Returning true keeps the caller from adding
    // its legacy fits on top; returning false would do exactly that.
    // INSTANTANEOUS rate at the current iterate, loss implicit.
    //
    // BDF2 asks for the right-hand side at t^{n+1}. A mean rate over the step
    // is not that, and supplying one is a Lie composition that caps the scheme
    // at first order however accurately the mean itself is computed. Evaluating
    // the rate at the current iterate makes the source implicit as the outer
    // loop converges, which is exactly why the legacy electron-impact path
    // reaches second order -- generalised here to the whole mechanism.
    //
    // The loss goes in through fvm::Sp so a fast sink cannot drive a density
    // negative and cannot break the Picard iteration.
    if (chemistrySolver_ == csImplicitRate || chemistrySolver_ == csAdaptive)
    {
        const label nSp = species_.nSpecies();
        const label nTab = rates_->size();

        chemP_.setSize(nSp);
        chemL_.setSize(nSp);
        forAll(chemP_, s)
        {
            chemP_[s].setSize(mesh_.nCells(), Zero); chemP_[s] = Zero;
            chemL_[s].setSize(mesh_.nCells(), Zero); chemL_[s] = Zero;
        }



        scalarField n(nSp, Zero), kTab(max(nTab, 1), Zero);
        scalarField P(nSp, Zero), L(nSp, Zero);
        const scalar dtNow = mesh_.time().deltaTValue();
        const scalarField& neI =
            species_.numberDensity(species_.electronSpeciesID()).primitiveField();

        forAll(neI, celli)
        {
            if (neI[celli] < chemActivityThreshold_) continue;

            for (label sp = 0; sp < nSp; ++sp)
            {
                n[sp] = species_.numberDensities()[sp].primitiveField()[celli];
            }
            for (label j = 0; j < nTab; ++j)
            {
                kTab[j] = rates_->k(j).primitiveField()[celli];
            }

            chem_->productionLoss(n, kTab, TgasCell(celli), P, L);

            // Stiffness, as the solver actually experiences it: L*dt is the
            // number of loss timescales crossed in one step. Large is not
            // unstable -- the implicit sink handles that -- but it is where
            // linearising around the current state stops being enough.
            scalar cellStiff = 0;
            for (label sp = 0; sp < nSp; ++sp)
            {
                cellStiff = max(cellStiff, L[sp]*dtNow);
            }
            chemStiffness_ = max(chemStiffness_, cellStiff);

            // What the linearised step would actually DO to this cell. L*dt
            // above is a diagonal estimate and misses the off-diagonal
            // instability entirely -- see chemChangeLimit_ in the header. This
            // predicts the step and measures it, which catches the mode the
            // estimate cannot see.
            scalar cellChange = 0;
            for (label sp = 0; sp < nSp; ++sp)
            {
                const scalar n0s = chemN0_[sp][celli];
                if (n0s <= chemChangeFloor_) continue;
                const scalar nPred =
                    (n0s + P[sp]*dtNow)/(1.0 + L[sp]*dtNow);
                cellChange = max(cellChange, mag(nPred - n0s)/n0s);
            }

            // ADAPTIVE: integrate this cell instead of linearising it, when
            // the linearisation cannot be expected to hold. The result is
            // returned to the SAME P/L form, so the equation assembly does not
            // know or care which path a cell took -- and a stiff cell keeps
            // the implicit sink, hence positivity, rather than reverting to a
            // bare explicit source.
            const bool tooStiff  = (cellStiff  > chemStiffLimit_);
            const bool tooBigStep = (cellChange > chemChangeLimit_);

            if (chemistrySolver_ == csAdaptive && (tooStiff || tooBigStep))
            {
                ++chemNstiff_;
                if (tooStiff)   ++chemNstiffByL_;
                if (tooBigStep) ++chemNstiffByChange_;

                scalarField n0(nSp), nEnd(nSp);
                for (label sp = 0; sp < nSp; ++sp)
                {
                    n0[sp] = chemN0_[sp][celli];
                }
                nEnd = n0;
                chem_->integrate(nEnd, kTab, TgasCell(celli), dtNow);

                // An ODE solver can fail on a step far longer than the
                // chemistry it is integrating, and it does not always say so:
                // it returns non-finite values. Feeding those into the matrix
                // raises SIGFPE inside the linear solver, several layers away
                // from the cause -- which is how this was found.
                //
                // Fall back to the linearised P/L already computed for this
                // cell. That is exactly the right fallback: it is stable at
                // any L*dt, only less accurate, which is the trade the cell
                // was going to make anyway.
                bool ok = true;
                for (label sp = 0; sp < nSp; ++sp)
                {
                    if (!std::isfinite(nEnd[sp]) || nEnd[sp] < 0)
                    {
                        ok = false;
                        break;
                    }
                }
                if (!ok)
                {
                    ++chemODEFailures_;
                }
                else

                {
                // Mean rate over the step, then projected onto exact charge
                // conservation -- the integrator does not preserve the linear
                // invariant exactly in floating point.
                scalarField rr(nSp);
                scalar netQ = 0, traffic = 0;
                for (label sp = 0; sp < nSp; ++sp)
                {
                    rr[sp] = (nEnd[sp] - n0[sp])/dtNow;
                    const scalar q = species_.speciesChargeNumbers()[sp];
                    netQ    += q*rr[sp];
                    traffic += mag(q*rr[sp]);
                }
                if (traffic > VSMALL && mag(netQ)/traffic < 1e-3)
                {
                    for (label sp = 0; sp < nSp; ++sp)
                    {
                        const scalar q = species_.speciesChargeNumbers()[sp];
                        if (mag(q) < SMALL) continue;
                        rr[sp] -= (netQ/traffic)*mag(rr[sp])*sign(q);
                    }
                }

                // Back into production / loss-coefficient form, so the sink
                // stays implicit. Exact at the current state by construction.
                for (label sp = 0; sp < nSp; ++sp)
                {
                    if (rr[sp] >= 0)
                    {
                        P[sp] = rr[sp];
                        L[sp] = 0;
                    }
                    else if (n[sp] > 1.0)
                    {
                        // Loss coefficient, capped so that a species being
                        // consumed to exhaustion cannot put an unbounded
                        // entry on the matrix diagonal. L*dt = 1e3 already
                        // removes all but 1e-3 of it in one step, so the cap
                        // is physically indistinguishable from complete
                        // consumption and numerically finite.
                        P[sp] = 0;
                        L[sp] = min(-rr[sp]/n[sp], 1.0e3/dtNow);
                    }
                    else
                    {
                        // Below one particle per cubic metre there is nothing
                        // to lose.
                        P[sp] = 0;
                        L[sp] = 0;
                    }
                }
                }
            }

            // Charge conservation needs NO projection in this mode. The source
            // is the instantaneous rate, so sum q_s (P_s - L_s n_s) is exactly
            // the right-hand-side residual -- zero for a balanced mechanism, no
            // integration involved to spoil it. Measured rather than assumed.
            if (!chemCellsReported_)
            {
                scalar netQ = 0, traffic = 0;
                for (label sp = 0; sp < nSp; ++sp)
                {
                    const scalar q = species_.speciesChargeNumbers()[sp];
                    const scalar rr = P[sp] - L[sp]*n[sp];
                    netQ    += q*rr;
                    traffic += mag(q*rr);
                }
                if (traffic > VSMALL)
                {
                    chemSourceChargeAfter_ =
                        max(chemSourceChargeAfter_, mag(netQ)/traffic);
                }
                ++chemActiveCount_;
            }

            for (label sp = 0; sp < nSp; ++sp)
            {
                chemP_[sp][celli] = P[sp];
                chemL_[sp][celli] = L[sp];
            }
        }

        // Peak-hold between reports, so a stiff step is not missed just
        // because it fell between two samples.
        chemStiffnessPeak_ = max(chemStiffnessPeak_, chemStiffness_);
        chemNstiffPeak_    = max(chemNstiffPeak_, chemNstiff_);
        chemChargePeak_    = max(chemChargePeak_, chemSourceChargeAfter_);


        const label ti = mesh_.time().timeIndex();
        const bool due =
            !chemCellsReported_
         || (chemReportInterval_ > 0 && ti - chemLastReport_ >= chemReportInterval_);

        // ---- transport-choice advisory, once ---------------------------
        //
        // Whether a species needs to be transported is a PHYSICAL question --
        // how far it travels within its own chemical lifetime, against the
        // scale it has to resolve -- and until now the user had to answer it
        // from memory. The loss coefficients L are already assembled here, so
        // the lifetime is free; report it beside the diffusion length and the
        // cell size so the choice can be checked rather than assumed.
        //
        // Printed after the first chemistry evaluation rather than at start-up,
        // because before it there is no L to report and a table of estimates
        // would be worth less than no table.
        // Peak-hold the loss coefficients on EVERY call, so that whenever the
        // advisory prints it reports the shortest lifetime reached so far
        // rather than the one that happened to hold on the sampled step.
        if (advisoryLpeak_.size() != chemL_.size())
        {
            advisoryLpeak_.setSize(chemL_.size(), 0.0);
        }
        forAll(chemL_, sp)
        {
            if (!chemL_[sp].empty())
            {
                advisoryLpeak_[sp] = max(advisoryLpeak_[sp], gMax(chemL_[sp]));
            }
        }

        // `due` is true on the very first call, because chemCellsReported_
        // starts false -- so requiring it to be ALREADY set defers the advisory
        // to the second report, i.e. a full chemReportInterval of steps in.
        if (!transportAdvisoryDone_ && due && chemCellsReported_)
        {
            transportAdvisoryDone_ = true;

            // Reference diffusivity for a neutral in air at 1 bar, 300 K.
            // Species carrying their own diffusivity model use it instead.
            const scalar Dref = advisoryD_;

            const scalar hCell =
                Foam::cbrt(gMin(mesh_.V().field()));

            Info<< nl
                << "  Transport advisory -- is each species' motion resolvable"
                << " within its own lifetime?" << nl
                << "  (diffusion length over the chemical lifetime, against the"
                << " smallest cell)" << nl
                << "  " << string(78, '-').c_str() << nl
                << "    species      transport        tau_chem [s]"
                << "   L_diff [m]    vs cell   verdict" << nl;

            forAll(chemL_, sp)
            {
                if (chemL_[sp].empty()) continue;

                const scalar Lmax = advisoryLpeak_[sp];
                const scalar tau = (Lmax > VSMALL) ? 1.0/Lmax : GREAT;

                // Per-species D from the mechanism where mechc supplied it,
                // scaled to the gas temperature; the single reference value is
                // only the fallback for older mechanisms.
                scalar D = Dref;
                bool haveOwnD = false;
                const word& sn = species_.speciesNames()[sp];
                if (species_.mechanismDiffusivity().found(sn))
                {
                    const scalar Tnow =
                        TgasField_ ? gAverage(TgasField_->primitiveField())
                                   : TgasConst_;
                    D = species_.mechanismDiffusivity()[sn]
                      * Foam::pow
                        (
                            Tnow/species_.diffusivityTref(),
                            species_.diffusivityExponent(sn)
                        );
                    haveOwnD = true;
                }

                const scalar Ldiff =
                    (tau < GREAT) ? Foam::sqrt(D*tau) : GREAT;
                const scalar ratio = Ldiff/max(hCell, VSMALL);

                const word tmName = transportModels_[sp].type();
                const bool charged =
                    mag(species_.speciesChargeNumber(sp)) > SMALL;

                // Charged species are not advised on: they drift, and drift is
                // not optional whatever the diffusion length says.
                word verdict;
                if (charged)              verdict = "charged: transport it";
                else if (ratio < 0.1)     verdict = "local is fine";
                else if (ratio < 1.0)     verdict = "local ok, marginal";
                else                      verdict = "CONSIDER diffusion";

                if (!charged && tmName == "diffusion" && ratio < 0.1)
                {
                    verdict = "diffusion not needed";
                }

                char row[256];
                snprintf
                (
                    row, sizeof(row),
                    "    %-12s %-15s %10.3g   %10.3g%s %8.3g   %s",
                    species_.speciesNames()[sp].c_str(),
                    tmName.c_str(), tau, Ldiff,
                    haveOwnD ? "  " : " *", ratio, verdict.c_str()
                );
                Info<< row << nl;
            }

            Info<< "  " << string(78, '-').c_str() << nl
                << "  smallest cell " << hCell << " m; reference D = "
                << Dref << " m^2/s used only where the mechanism"
                << " carries no D (override: chemistry/advisoryDiffusivity)"
                << nl
                << "  Convection is not assessed: this solver carries no"
                << " momentum equation." << endl;
        }

        if (due)
        {
            chemCellsReported_ = true;
            chemLastReport_ = ti;

            // Report the PEAK since the last report, not this step's value.
            chemSourceChargeAfter_ = chemChargePeak_;
            chemStiffness_ = chemStiffnessPeak_;
            chemNstiff_ = chemNstiffPeak_;
            chemStiffnessPeak_ = 0;
            chemNstiffPeak_ = 0;
            chemChargePeak_ = 0;

            reduce(chemSourceChargeAfter_, maxOp<scalar>());
            reduce(chemActiveCount_, sumOp<label>());
            reduce(chemStiffness_, maxOp<scalar>());
            reduce(chemNstiff_, sumOp<label>());
            // Cumulative like chemNstiff_ itself, and reduced the same way, so
            // the three numbers in the report are directly comparable.
            reduce(chemNstiffByL_, sumOp<label>());
            reduce(chemNstiffByChange_, sumOp<label>());
            reduce(chemODEFailures_, sumOp<label>());
            if (chemODEFailures_ > 0)
            {
                WarningInFunction
                    << chemODEFailures_ << " cell(s) had the stiff integration"
                    << " return non-finite values and fell back to the"
                    << " linearised source." << nl
                    << "    deltaT is far longer than the chemistry in those"
                    << " cells. The fallback is stable but less accurate;"
                    << " reduce deltaT if they matter." << endl;
            }
            if (chemistrySolver_ == csAdaptive)
            {
                // Split by CRITERION, because the two catch different things
                // and the difference is the diagnostic. Cells caught only by
                // the change test are ones L*dt called safe and the linearised
                // step would have mangled -- the off-diagonal mode. If that
                // count is ever large, the timestep is past the linearisation's
                // limit for this mechanism, whatever L*dt says.
                Info<< "plasmaChemistry: adaptive, " << chemNstiff_ << " of "
                    << chemActiveCount_ << " cells integrated, the rest"
                    << " linearised" << nl
                    << "    by L*dt > " << chemStiffLimit_ << ": "
                    << chemNstiffByL_
                    << ",  by step change > " << chemChangeLimit_ << ": "
                    << chemNstiffByChange_
                    << "  (overlap counted in both)" << endl;
            }
            Info<< "plasmaChemistry: " << (chemistrySolver_ == csAdaptive
                                           ? "adaptive" : "implicitRate")
                << ", " << chemActiveCount_
                << " active cells, max(L*dt) = " << chemStiffness_ << nl
                << "    charge residual of the SOURCE, worst cell: "
                << chemSourceChargeAfter_
                << " (exact by construction: no integration)" << endl;
        }

        // Has the Picard iteration converged inside the outer loop?
        //
        // This is the question that decides implicitRate versus ode, and it is
        // measurable rather than a matter of judgement: compare this
        // iteration's source with the previous one. A contraction means the
        // linearisation is resolving the chemistry; no contraction means the
        // chemistry is stiffer than dt and the source is being extrapolated
        // rather than converged.
        {
            scalar num = 0, den = 0;
            forAll(chemP_, sp)
            {
                const scalarField& nsp =
                    species_.numberDensities()[sp].primitiveField();
                forAll(chemP_[sp], c)
                {
                    const scalar srcNow = chemP_[sp][c] - chemL_[sp][c]*nsp[c];
                    const scalar srcOld =
                        (chemSrcPrev_.size() > sp && chemSrcPrev_[sp].size() > c)
                      ? chemSrcPrev_[sp][c] : 0.0;
                    num += sqr(srcNow - srcOld);
                    den += sqr(srcNow);
                }
            }
            reduce(num, sumOp<scalar>());
            reduce(den, sumOp<scalar>());

            const scalar change = (den > VSMALL) ? Foam::sqrt(num/den) : 0.0;
            if (chemOuterCount_ > 0)
            {
                chemPicardChange_ = change;
            }
            ++chemOuterCount_;

            chemSrcPrev_.setSize(chemP_.size());
            forAll(chemP_, sp)
            {
                const scalarField& nsp =
                    species_.numberDensities()[sp].primitiveField();
                chemSrcPrev_[sp].setSize(chemP_[sp].size());
                forAll(chemP_[sp], c)
                {
                    chemSrcPrev_[sp][c] = chemP_[sp][c] - chemL_[sp][c]*nsp[c];
                }
            }
        }

        forAll(chemP_, sp)
        {
            if (sp >= eqns.size() || !eqns[sp]) continue;

            volScalarField Pf
            (
                IOobject("chemP", mesh_.time().timeName(), mesh_,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh_, dimensionedScalar(dimensionSet(0,-3,-1,0,0,0,0), Zero)
            );
            volScalarField Lf
            (
                IOobject("chemL", mesh_.time().timeName(), mesh_,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh_, dimensionedScalar(dimensionSet(0,0,-1,0,0,0,0), Zero)
            );
            Pf.primitiveFieldRef() = chemP_[sp];
            Lf.primitiveFieldRef() = chemL_[sp];

            *eqns[sp] -= Pf;
            *eqns[sp] += fvm::Sp(Lf, species_.numberDensity(sp));
        }
        return true;
    }

    // With the source formulation the chemistry is delivered HERE, as an
    // explicit source per species, rather than by advancing the fields.
    if (chemistrySolver_ == csODE)
    {
        forAll(chemRR_, s)
        {
            if (s >= eqns.size() || !eqns[s] || chemRR_[s].empty()) continue;

            volScalarField src
            (
                IOobject
                (
                    "chemRR_" + species_.speciesNames()[s],
                    mesh_.time().timeName(), mesh_,
                    IOobject::NO_READ, IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), Zero)
            );
            src.primitiveFieldRef() = chemRR_[s];

            *eqns[s] -= src;
        }
        return true;
    }

    // Interpolation only; cheap, and never throttled -- skipping it would leave
    // k_j frozen at its t=0 value while the discharge moves. The EEDF re-solve
    // is the expensive thing, and that is governed separately.
    rates_->correct();

    const scalarField& neI = ne.primitiveField();
    const label nCells = mesh_.nCells();

    // Drift speed from the registered mobility field the mobility model already
    // maintains -- so it cannot disagree with the species dictionary, which is
    // precisely what the legacy hard-coded mu_ref = 2.398 could do.
    const label eIdx = species_.electronSpeciesID();
    const word muName = "mu_" + species_.speciesNames()[eIdx];
    scalarField vDrift(nCells, Zero);
    if (mesh_.foundObject<volScalarField>(muName))
    {
        const scalarField& mu =
            mesh_.lookupObject<volScalarField>(muName).primitiveField();
        const scalarField& E = Emag.primitiveField();
        forAll(vDrift, c) { vDrift[c] = mu[c]*E[c]; }
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find mobility field '" << muName << "'." << nl
            << "The electron transport model must be driftDiffusion for the"
            << " mechanism source terms to have a drift velocity." << nl
            << exit(FatalError);
    }

    // Background reservoir left, recomputed each step so it depletes as the
    // discharge burns through it rather than being held fixed.
    scalarField nBg(nCells, species_.backgroundDensity().value());
    for (const label n : species_.neutralSpeciesIDs())
    {
        nBg -= species_.numberDensities()[n].primitiveField();
    }
    forAll(nBg, c) { nBg[c] = max(nBg[c], scalar(0)); }

    if (rates_->eedfRefreshDue())
    {
        // Composition weighted by ELECTRON DENSITY, not by volume.
        //
        // The rate coefficients describe what electrons do, so the mixture that
        // matters is the one the electrons are actually in. A volume average
        // over a domain that is mostly undisturbed gas would report very nearly
        // the initial composition no matter what the streamer head had done to
        // the gas inside it -- so the refresh would faithfully re-solve for a
        // mixture that exists nowhere near the chemistry, at full cost.
        //
        // Mole fractions are taken over the NEUTRALS only, which is what an
        // EEDF solve means by composition: the electron-impact cross sections
        // are per neutral target, and ions are a ~1e-6 fraction whose own
        // cross sections this mechanism does not carry.
        HashTable<scalar> comp;
        {
            scalar wSum = 0;
            forAll(neI, c) { wSum += neI[c]*mesh_.V()[c]; }
            reduce(wSum, sumOp<scalar>());

            if (wSum > VSMALL)
            {
                // Transported neutrals first, then whatever of the background
                // the mechanism's reference composition names. nBg already has
                // the transported neutrals subtracted, so the two do not
                // double-count.
                scalar nTot = 0;
                HashTable<scalar> weighted;

                forAll(species_.speciesNames(), i)
                {
                    if (mag(species_.speciesChargeNumbers()[i]) > SMALL) continue;

                    const scalarField& ni =
                        species_.numberDensities()[i].primitiveField();
                    scalar acc = 0;
                    forAll(ni, c) { acc += neI[c]*mesh_.V()[c]*ni[c]; }
                    reduce(acc, sumOp<scalar>());
                    weighted.set(species_.speciesNames()[i], acc/wSum);
                    nTot += acc/wSum;
                }

                scalar bgAcc = 0;
                forAll(nBg, c) { bgAcc += neI[c]*mesh_.V()[c]*nBg[c]; }
                reduce(bgAcc, sumOp<scalar>());
                const scalar bgAvg = bgAcc/wSum;

                for (const entry& e : mechComposition_)
                {
                    const scalar x = readScalar(e.stream());
                    weighted.set(e.keyword(), weighted.lookup(e.keyword(), 0.0)
                                            + x*bgAvg);
                    nTot += x*bgAvg;
                }

                if (nTot > VSMALL)
                {
                    forAllConstIters(weighted, it)
                    {
                        comp.set(it.key(), it.val()/nTot);
                    }
                }
            }
        }

        // Gas temperature: the background value, until an energy equation for
        // the heavy species exists. Passing -1 would keep the dictionary's,
        // which is the same thing today but stops being so the moment gas
        // heating is solved.
        // One representative temperature for a global table rebuild. The MEAN
        // over active cells, not the peak: a single hot cell should not set
        // the EEDF the whole domain then uses.
        scalar Tref = TgasConst_;
        if (gasHeating_)
        {
            if (TgasField_) Tref = gAverage(TgasField_->primitiveField());
        }
        rates_->refreshEEDF(comp, Tref);
    }

    List<scalarField> src(species_.nSpecies(), scalarField(nCells, Zero));
    scalarField sIon(nCells, Zero), sAtt(nCells, Zero);

    // speciesID() looks up a HashTable with at(), which THROWS on a missing
    // key rather than returning -1. Mechanism species that the case does not
    // transport are entirely normal -- N2 and O2 are background gas, and most
    // excited states will not be carried -- so membership has to be tested
    // without calling it.
    HashSet<word> transported;
    forAll(species_.speciesNames(), i)
    {
        transported.insert(species_.speciesNames()[i]);
    }
    auto idOf = [&](const word& nm) -> label
    {
        return transported.found(nm) ? species_.speciesID(nm) : -1;
    };

    forAll(rates_->reactions(), r)
    {
        const auto& rx = rates_->reactions()[r];
        const label tgt = idOf(rx.target);

        scalarField Rj(nCells);
        // Three-body processes always take the k path. alpha_j/N is a
        // SECOND-order construction -- nu_j/N divided by drift speed -- so it
        // has neither the units nor the collider density a third-order rate
        // needs, and applying the collider factor on top of it would count the
        // density twice.
        if (sourceModel_ == smTownsend && rates_->hasTownsend(r)
         && !rx.threeBody())
        {
            // alpha_j/N is already mole-fraction weighted (it is nu_j/N over
            // v_drift), so neither a composition factor nor a target density
            // belongs here -- both are inside it.
            const scalarField& aj = rates_->townsend(r).primitiveField();
            forAll(Rj, c) { Rj[c] = aj[c]*nBg[c]*vDrift[c]*neI[c]; }
        }
        else
        {
            // k_j is UNWEIGHTED -- per unit target density. Do NOT multiply by
            // a mole fraction; the tables already exclude it.
            const scalarField& kj = rates_->k(r).primitiveField();
            scalarField nTarget(nCells);
            if (tgt >= 0)
            {
                nTarget = species_.numberDensities()[tgt].primitiveField();
            }
            else
            {
                const scalar x =
                    mechComposition_.getOrDefault<scalar>(rx.target, 0.0);
                if (x <= 0)
                {
                    FatalErrorInFunction
                        << "Reaction " << rx.id << " target '" << rx.target
                        << "' is neither transported nor in the mechanism"
                        << " composition" << nl << exit(FatalError);
                }
                forAll(nTarget, c) { nTarget[c] = x*nBg[c]; }
            }
            forAll(Rj, c) { Rj[c] = kj[c]*neI[c]*nTarget[c]; }

            // Third body. e + O2 + M -> O2- + M is the dominant electron loss
            // channel in atmospheric air at low E/N, and it is third order:
            // without the collider density the rate is not merely inaccurate,
            // it has the wrong units and the wrong pressure scaling. rateScale
            // carries the cm^3 -> m^3 conversion for the density-normalised
            // curve LXCat tabulates, so that k*rateScale is in m^6/s.
            if (rx.threeBody())
            {
                const label cid = idOf(rx.collider);
                scalarField nCol(nCells);
                if (cid >= 0)
                {
                    nCol = species_.numberDensities()[cid].primitiveField();
                }
                else
                {
                    const scalar xc =
                        mechComposition_.getOrDefault<scalar>(rx.collider, 0.0);
                    if (xc <= 0)
                    {
                        FatalErrorInFunction
                            << "Reaction " << rx.id << " has collider '"
                            << rx.collider << "', which is neither transported"
                            << " nor in the mechanism composition" << nl
                            << exit(FatalError);
                    }
                    forAll(nCol, c) { nCol[c] = xc*nBg[c]; }
                }
                forAll(Rj, c) { Rj[c] *= rx.rateScale*nCol[c]; }
            }
        }

        // Both sides are applied by the SAME loop, with only the sign
        // differing. The left-hand side used to be handled by hand -- subtract
        // the target, and remember that an electron is consumed too -- and the
        // "remember" is what failed: the electron was never subtracted, so
        // every excitation created one out of nothing and every ionisation
        // created two instead of one. The mechanism now states its reactants
        // explicitly, and the code below cannot treat the two sides
        // inconsistently because it does not distinguish them.
        //
        // Counted WITH MULTIPLICITY throughout: (N2+ e e) really does add two
        // electrons, and e + N2 really does consume one.
        label dNe = 0;
        for (label side = 0; side < 2; ++side)
        {
            const wordList& names = side ? rx.products : rx.reactants;
            const scalar sign = side ? 1.0 : -1.0;

            forAll(names, i)
            {
                const word& nm = names[i];

                // The electron is matched through the mechanism's declared
                // name, not a hard-coded "e": a case is free to call its own
                // electron species anything, and a literal comparison would
                // silently match nothing for a case that did.
                const bool isElectron = (nm == rates_->electronSpecies());
                if (isElectron) { dNe += side ? 1 : -1; }

                const label pid = isElectron ? eIdx : idOf(nm);
                if (pid >= 0)
                {
                    if (side) { src[pid] += Rj; } else { src[pid] -= Rj; }
                    continue;
                }
                if (!side)
                {
                    // An untransported REACTANT is normal and benign: N2 and O2
                    // are background gas, held fixed by construction. Its
                    // depletion is not tracked, which is the whole point of
                    // calling it background.
                    continue;
                }
                {
                // Discarding an untransported product is NOT uniformly benign,
                // and treating it so is how a run looks healthy while being
                // wrong. Two different situations:
                //
                //   neutral  losing an excited state is an approximation. The
                //            species is absent; the field is unaffected.
                //
                //   CHARGED  breaks charge conservation. Ionisation then adds
                //            electrons with no counter-charge, net charge grows
                //            out of nothing, and the Poisson solve is silently
                //            corrupted. Measured on this very case: +1.45e18
                //            electrons per step against +0 positive ions.
                //
                // So the charged case is fatal by default. `pIon`-style lumped
                // models are still available through allowChargeNonConservation,
                // but they have to be asked for.
                const label q = rates_->chargeOf(nm);
                if (q != 0 && !allowChargeNonConservation_)
                {
                    FatalErrorInFunction
                        << "Reaction " << rx.id << " produces '" << nm
                        << "' (charge " << q << "), which is not a transported"
                        << " species." << nl
                        << "    Discarding a CHARGED product breaks charge"
                        << " conservation: electrons would be created without"
                        << " their counter-charge and the Poisson solve would be"
                        << " wrong." << nl
                        << "    Either add '" << nm << "' to activeSpecies, or"
                        << " set allowChargeNonConservation true in the"
                        << " chemistry dictionary to accept a lumped-ion model."
                        << nl << exit(FatalError);
                }
                if (!warnedMissing_.found(nm))
                {
                    warnedMissing_.insert(nm);
                    WarningInFunction
                        << "Reaction " << rx.id << " produces '" << nm
                        << "' (neutral), not a transported species; discarded."
                        << " Its population is lost, but charge is unaffected."
                        << endl;
                }
                }
            }
        }

        // Electron balance from the reaction's own stoichiometry, both sides:
        // +1 ionisation, -1 attachment (real attachment, its own cross
        // section), 0 excitation/dissociation.
        const label dE = dNe;
        if (dE > 0)      { sIon += dE*Rj; }
        else if (dE < 0) { sAtt += (-dE)*Rj; }
    }

    // Charge-balance diagnostic. Sum q_i S_i over the transported species: for
    // a charge-conserving mechanism this is identically zero, so any drift is a
    // species that is being created or destroyed without its counter-charge.
    // Reported as a fraction of the ionisation source so it is scale-free --
    // an absolute number means nothing without knowing how much chemistry is
    // happening.
    {
        scalar netQ = 0, refQ = 0;
        forAll(src, sIdx)
        {
            // The CASE's charge number, not the mechanism's. The two name the
            // electron differently -- `e` here, `Electron` in the mechanism --
            // so looking the electron up in the mechanism table returns 0 and
            // drops the one species the diagnostic exists to track, turning it
            // into a permanent false alarm. rates_->chargeOf() is for products
            // that are NOT transported, where the case has no opinion.
            const scalar q = species_.speciesChargeNumbers()[sIdx];
            if (mag(q) < SMALL) continue;
            forAll(src[sIdx], c)
            {
                netQ += q*src[sIdx][c];

                // Normalise against the total charge TRAFFIC, not against the
                // ionisation source. netQ is a near-total cancellation, so its
                // roundoff floor is set by the largest term being cancelled --
                // and once three-body attachment is present that is attachment,
                // running at ~6e7 s^-1, not ionisation. Dividing by the
                // ionisation source made pure roundoff read as 2.9e-05 and trip
                // a threshold it had passed at 2.5e-09 the run before.
                refQ += mag(q*src[sIdx][c]);
            }
        }
        reduce(netQ, sumOp<scalar>());
        reduce(refQ, sumOp<scalar>());

        // Tolerance is 1e-6, not machine epsilon: netQ is a near-total
        // cancellation over ~1e6 cells, so a few orders above double precision
        // is arithmetic rather than physics. A genuinely missing charged
        // product is an O(1) violation of this ratio -- every reaction that
        // creates a charge without its counter-charge contributes with the
        // same sign -- so nothing real hides underneath.
        // Checked EVERY step, and the running maximum is what gets reported.
        // Reporting only the first step would grade the run on its quietest
        // moment, before the discharge has produced any chemistry to unbalance.
        if (mag(refQ) > VSMALL)
        {
            const scalar rel = mag(netQ)/mag(refQ);
            const bool worse = rel > 10*chargeRelMax_;
            chargeRelMax_ = max(chargeRelMax_, rel);

            if (rel > 1e-6 && !chargeWarned_)
            {
                WarningInFunction
                    << "charge is not conserved by the mechanism source terms: "
                    << "sum(q_i S_i) = " << netQ << ", which is " << rel
                    << " of the total charge traffic sum|q_i S_i|." << nl
                    << "    Every charged product must be a transported species."
                    << endl;
                chargeWarned_ = true;
            }
            else if (worse)
            {
                // Reported even when it passes: a check whose only evidence of
                // success is the absence of output cannot be distinguished from
                // a check that never ran. Rate-limited to order-of-magnitude
                // increases so it does not print every step.
                Info<< "plasmaTransport: charge balance OK, |sum(q_i S_i)| = "
                    << rel << " of sum|q_i S_i| (roundoff)" << endl;
            }
        }
    }

    // eqns is dense over species ids -- List<autoPtr<fvScalarMatrix>> sized
    // nSpecies() with every entry filled. Checked in the source.
    forAll(src, s)
    {
        if (eqns.size() > s && eqns[s])
        {
            volScalarField S
            (
                IOobject
                (
                    "S_" + species_.speciesNames()[s],
                    mesh_.time().timeName(), mesh_,
                    IOobject::NO_READ, IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), Zero)
            );
            S.primitiveFieldRef() = src[s];
            *eqns[s] -= S;
        }
    }

    // S_iz_, k_eff_ and alpha_ from the REAL reactions -- no Townsend fit and
    // no effective eta anywhere. alpha_ becomes its own definition,
    // S_iz/(n_e v_drift), rather than a fit to it, so the AMR criterion and the
    // photoionization model keep reading it unchanged.
    {
        scalarField& siz  = S_iz_.primitiveFieldRef();
        scalarField& keff = k_eff_.primitiveFieldRef();
        scalarField& a    = alpha_.primitiveFieldRef();
        forAll(siz, c)
        {
            siz[c]  = sIon[c];
            const scalar n = max(neI[c], SMALL);
            keff[c] = (sIon[c] - sAtt[c])/n;
            a[c]    = sIon[c]/(n*max(vDrift[c], SMALL));
        }
    }
    S_iz_.correctBoundaryConditions();
    k_eff_.correctBoundaryConditions();
    alpha_.correctBoundaryConditions();

    return true;
}


// ************************************************************************* //












