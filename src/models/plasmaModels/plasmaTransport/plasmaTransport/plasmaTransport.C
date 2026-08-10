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
#include "IFstream.H"
#include "plasmaTransportModel.H"
#include "plasmaWallBC.H"
#include "photoionizationModel.H"

// Remove these headers later
#include "interpolationTable.H"
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

plasmaTransport::~plasmaTransport() = default;

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
void plasmaTransport::solve()
{
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

    {
        IFstream is(cd.get<fileName>("mechanism"));
        dictionary mech(is);
        if (mech.found("composition"))
        {
            mechComposition_ = mech.subDict("composition");
        }
    }

    Info<< "plasmaTransport: sourceModel " << model
        << ", " << rates_->size() << " electron-impact reactions" << nl
        << "    NOTE this path evaluates ELECTRON-IMPACT reactions only."
        << " Heavy chemistry (recombination," << nl
        << "    quenching, detachment) lives in the Cantera mechanism and is"
        << " NOT applied here, so ion and" << nl
        << "    excited-state densities have no loss channel in this mode."
        << " Adequate for a short pulse;" << nl
        << "    not for long transients. Stage 2 hands the chemistry to Cantera."
        << endl;
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
        const scalar Tgas =
            species_.backgroundDict().subOrEmptyDict("energy")
                .getOrDefault<scalar>("T", -1);

        rates_->refreshEEDF(comp, Tgas);
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












