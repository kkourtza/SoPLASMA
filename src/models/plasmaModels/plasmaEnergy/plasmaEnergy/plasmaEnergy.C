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
#include "DynamicList.H"
#include "FlatOutput.H"
#include "Pair.H"
#include "OSspecific.H"
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
                // NOT gated on active() -- contribute() only records the
                // iterate so the joint residual (and rho) can be formed. See
                // plasmaOuterRelaxation::enrol().
                if (r)
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

    // Or is the electron carrying its own energy equation? Only LMEA does.
    //
    // ASKED, not re-derived: plasmaSpecies resolved the closure once, from a
    // single top-level key, and this is the one condition that follows from it.
    // Note the two are deliberately OR-ed, never coupled -- gas heating without
    // LMEA and LMEA without gas heating are both normal, and an earlier version
    // that keyed LMEA off `energy { solve }` made enabling one switch on the
    // other silently.
    return species.isLMEA();
}


namespace
{

// Validate the leaves the CASE wrote, and fill in a lookupVariable that
// `fromMechanism` would otherwise default to the WRONG field.
//
// The split follows the one already drawn in plasmaTransport for the chemistry
// key: `lookupVariable` is a FIELD NAME and is unambiguously fixed by the
// closure, so a contradiction is derived-and-validated (fatal); `quantity`,
// `file` and `tableDir` name TABLE FILES, a different namespace, where a case
// may legitimately differ, so those are defaulted and never validated.
void validateElectronCoeffs(Foam::dictionary& c, const Foam::word& key)
{
    using namespace Foam;

    // The six leaves that must follow the closure, and the one that must not.
    const List<Pair<word>> leaves
    ({
        Pair<word>("mobility",           key),
        Pair<word>("diffusivity",        key),
        Pair<word>("energyMobility",     key),
        Pair<word>("energyDiffusivity",  key),
        Pair<word>("powerLoss/elastic",   key),
        Pair<word>("powerLoss/inelastic", key),
        // STRUCTURAL carve-out, fatal in the OPPOSITE direction.
        Pair<word>("initialMeanEnergy",  word("reducedE"))
    });

    for (const Pair<word>& leaf : leaves)
    {
        const word& path = leaf.first();
        const word& want = leaf.second();

        const auto slash = path.find('/');
        const bool nested = (slash != std::string::npos);

        word parent;
        word name(path);

        if (nested)
        {
            parent = word(path.substr(0, slash), false);
            name   = word(path.substr(slash + 1), false);
        }

        dictionary* owner = &c;

        if (nested)
        {
            owner = c.found(parent) ? &c.subDict(parent) : nullptr;
        }

        if (!owner || !owner->found(name)) continue;

        dictionary& pd = owner->subDict(name);
        const word type = pd.getOrDefault<word>("type", "fromMechanism");

        // The density-division trap. PelasticN/PinelasticN are already per
        // density and meanEnergy is an absolute eV value, so fromMechanism
        // divides by the gas density a second time.
        if
        (
            type == "fromMechanism"
         && (path == "initialMeanEnergy" || nested)
        )
        {
            FatalErrorInFunction
                << "energyModelCoeffs/" << path << " is `type fromMechanism`,"
                << " which DIVIDES BY THE GAS NUMBER DENSITY." << nl
                << "    That is correct for the reduced quantities (muN, DLN,"
                << " muEpsN, DEpsN) but wrong here: PelasticN and PinelasticN"
                << " are ALREADY per density, and meanEnergy is an absolute"
                << " energy in eV." << nl
                << "    Measured consequence: a loss term 2.4e25 times too"
                << " small, so the electrons heat with no sink and the mean"
                << " energy climbs to the clamp. The seed comes out at"
                << " ~4e-26 eV." << nl
                << "    Use `type tabulated1D` with an explicit `file`." << nl
                << exit(FatalError);
        }

        if (pd.found("lookupVariable"))
        {
            const word given = pd.get<word>("lookupVariable");

            if (given == want) continue;

            if (path == "initialMeanEnergy")
            {
                FatalErrorInFunction
                    << "energyModelCoeffs/initialMeanEnergy/lookupVariable is `"
                    << given << "`, but it must be `reducedE`." << nl
                    << "    TWO independent reasons, both structural:" << nl
                    << "      1. the Boltzmann sweep excludes `meanEnergy` from"
                    << " its mean-energy loop, so meanEnergy_vs_meanE is NEVER"
                    << " written -- it would be the identity;" << nl
                    << "      2. this entry SEEDS the mean energy, so it is"
                    << " evaluated before `meanE` exists." << nl
                    << "    This is the one entry that does not follow"
                    << " electronEnergyModel. Delete it and it is derived."
                    << nl
                    << exit(FatalError);
            }

            FatalErrorInFunction
                << "energyModelCoeffs/" << path << "/lookupVariable is `"
                << given << "` but electronEnergyModel requires `" << want
                << "`." << nl
                << "    That combination is the HALF-LMEA: transport and"
                << " losses responding to one variable while the reaction"
                << " rates respond to another. It is a measured runaway and it"
                << " produces no error of its own." << nl
                << "    The key is DERIVED from electronEnergyModel. Delete it."
                << nl
                << exit(FatalError);
        }
        else if (type == "fromMechanism")
        {
            // MechanismProperty's OWN default is `reducedE`, so an omitted key
            // inside an LMEA block silently read the LFA table -- a working
            // half-LMEA with no warning anywhere. Insert the derived key.
            pd.add("lookupVariable", want);
        }
    }
}

} // End anonymous namespace


Foam::dictionary plasmaEnergy::resolveEnergyModelCoeffs
(
    const word& modelName,
    const label speciesIndex,
    const dictionary& userCoeffs
) const
{
    // Only the electron under LMEA is configured from the model. Everything
    // else gets its own block back, byte for byte.
    //
    // SCOPED TO THE ELECTRON on purpose: meanE, muEpsN and PelasticN are all
    // electron quantities, and ion coefficients are functions of E/N or |E|
    // and live in each species' own driftDiffusionCoeffs. Nothing here can
    // reach them.
    if
    (
        modelName != "localEnergy"
     || speciesIndex != species_.electronSpeciesID()
    )
    {
        return userCoeffs;
    }

    // The table directory. Defaulted from the chemistry dictionary so the
    // energy model reads the same tables the sweep just wrote; a case that
    // keeps its energy tables elsewhere overrides per entry, which is why this
    // is DEFAULTED and never validated.
    fileName tableDir("constant/plasmaTables");
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

        tableDir = transportDict.subOrEmptyDict("chemistry")
            .getOrDefault<fileName>("tableDir", tableDir);
    }

    // PRE-FLIGHT, before any evaluator is built.
    //
    // TabulatedProperty1D constructs its table in its member-init list, so its
    // own nicer error is unreachable and a missing defaulted file surfaces as
    // interpolationTable's raw "cannot open file". Worse, MechanismProperty's
    // missing-table error names a DIFFERENT likely cause -- "the case has no
    // chemistry dictionary" -- which is false here and actively misleading.
    //
    // Only run on the defaulting path: a case that writes the block itself may
    // point anywhere it likes, and its own evaluators will complain.
    const bool haveMeanE = isFile(tableDir/"muN_vs_meanE");
    const bool haveReducedE = isFile(tableDir/"muN_vs_reducedE");

    if (!haveMeanE)
    {
        if (haveReducedE)
        {
            // The sweep RAN and DECLINED. This is the real diagnosis.
            FatalErrorInFunction
                << "`" << tableDir
                << "` has muN_vs_reducedE but no muN_vs_meanE, so this"
                << " mechanism cannot be run under LMEA." << nl
                << "    The Boltzmann sweep writes mean-energy-keyed tables"
                << " only where the mean electron energy is STRICTLY"
                << " INCREASING with E/N. It is not, for this mixture." << nl
                << "    That happens in ATTACHMENT-DOMINATED mixtures:"
                << " attachment removes electrons energy-selectively, so"
                << " <eps> can fall as E/N rises, and a mean-energy lookup"
                << " would then be ambiguous -- two fields, two different"
                << " coefficients, one key. The sweep says so in its own log"
                << " (\"mean electron energy is not strictly increasing\")."
                << nl
                << "    Note the k_*_vs_meanE rate tables are missing for the"
                << " same reason, so this is not a limitation of the energy"
                << " equation alone: where the key is not invertible, the LMEA"
                << " closure itself does not apply." << nl
                << "    Remedies, in order:" << nl
                << "      1. run this case under LFA -- set"
                << " `electronEnergyModel LFA;` and everything else follows;"
                << nl
                << "      2. change the sweep's E/N range so the"
                << " non-monotonic stretch falls outside it;" << nl
                << "      3. if you have validated mean-energy tables from"
                << " elsewhere, write `energyModelCoeffs` explicitly and point"
                << " at them." << nl
                << exit(FatalError);
        }

        FatalErrorInFunction
            << "`electronEnergyModel LMEA` configures itself from the"
            << " Boltzmann sweep's tables, and `" << tableDir
            << "` holds none." << nl
            << "    The sweep runs only when the case has a `chemistry`"
            << " dictionary in plasmaTransportProperties. Add one, or write"
            << " `energyModelCoeffs` yourself (mobility, diffusivity,"
            << " energyMobility, energyDiffusivity, powerLoss/elastic,"
            << " powerLoss/inelastic, initialMeanEnergy)." << nl
            << exit(FatalError);
    }

    // Start from the case's own block and fill in whole LEAVES it omitted.
    //
    // NOT dictionary::merge: OpenFOAM's merge RECURSES into sub-dictionaries,
    // so a user writing `mobility { lookupVariable reducedE; }` would inherit
    // type/quantity from the default and get a WORKING half-LMEA. The rule is
    // therefore: if you write a property sub-dictionary it is yours entirely;
    // if you omit it, it is derived entirely. There is no half-inherited
    // sub-dictionary.
    dictionary c(userCoeffs);

    const word& key = species_.electronLookupKey();   // `meanE` under LMEA

    DynamicList<word> defaulted;
    DynamicList<word> fromCase;

    // --- the four REDUCED, table-keyed transport coefficients -------------
    //
    // `fromMechanism` divides by the gas number density, which is exactly
    // right for muN / DLN / muEpsN / DEpsN because those tables hold reduced
    // quantities. tableDir is omitted deliberately: MechanismProperty already
    // defaults it to constant/plasmaTables, and writing it here would be a
    // second copy of the same constant.
    auto addMech = [&](const word& name, const word& quantity)
    {
        if (c.found(name))
        {
            fromCase.append(name);
            return;
        }

        dictionary d;
        d.add("type", word("fromMechanism"));
        d.add("quantity", quantity);
        d.add("tableDir", tableDir);
        d.add("lookupVariable", key);
        c.add(name, d);
        defaulted.append(name + " (" + quantity + "_vs_" + key + ")");
    };

    addMech("mobility",          "muN");
    addMech("diffusivity",       "DLN");
    addMech("energyMobility",    "muEpsN");
    addMech("energyDiffusivity", "DEpsN");

    // --- the ALREADY-PER-DENSITY tables ----------------------------------
    //
    // tabulated1D, NEVER fromMechanism. PelasticN / PinelasticN are already
    // per density and meanEnergy is an absolute energy in eV, so running them
    // through fromMechanism divides by the gas density a second time. Measured
    // consequence: a loss term 2.4e25 times too small, the electrons heat with
    // no sink, and the mean energy climbs to the clamp. Same trap gives
    // ~4e-26 eV for the seed.
    auto addTable = [&](const word& name, const word& file, const word& lookup)
    {
        if (c.found(name))
        {
            fromCase.append(name);
            return;
        }

        dictionary d;
        d.add("type", word("tabulated1D"));
        d.add("file", fileName(tableDir/file));
        d.add("lookupVariable", lookup);
        c.add(name, d);
        defaulted.append(name + " (" + file + ")");
    };

    // initialMeanEnergy is keyed on reducedE, and that is STRUCTURAL, not a
    // convention: the sweep excludes `meanEnergy` from its mean-energy loop, so
    // meanEnergy_vs_meanE is never written -- it would be the identity -- and
    // the seed is needed before `meanE` exists at all.
    addTable("initialMeanEnergy", "meanEnergy_vs_reducedE", "reducedE");

    // powerLoss defaults PER LEAF, not as a block: a case writing only
    // `elastic` still gets `inelastic`. Omitting a sink can only ever make the
    // equation wrong in the runaway direction, and a user who genuinely wants
    // none can say so explicitly with `inelastic { type constant; value 0; }`.
    dictionary& pl = c.subDictOrAdd("powerLoss");
    auto addLoss = [&](const word& name, const word& file)
    {
        if (pl.found(name))
        {
            fromCase.append("powerLoss/" + name);
            return;
        }

        dictionary d;
        d.add("type", word("tabulated1D"));
        d.add("file", fileName(tableDir/file));
        d.add("lookupVariable", key);
        pl.add(name, d);
        defaulted.append("powerLoss/" + name + " (" + file + ")");
    };

    addLoss("elastic",   "PelasticN_vs_meanE");
    addLoss("inelastic", "PinelasticN_vs_meanE");

    // --- validate what the CASE wrote ------------------------------------
    validateElectronCoeffs(c, key);

    // --- report, per the rule that a fallback is never silent ------------
    Info<< "plasmaEnergy: energyModelCoeffs for `"
        << species_.speciesNames()[speciesIndex]
        << "` resolved from electronEnergyModel LMEA." << nl
        << "    tables: " << tableDir << nl;

    if (defaulted.size())
    {
        Info<< "    DERIVED from the model (the case omitted these): "
            << flatOutput(defaulted) << nl
            << "      -- every one is the Boltzmann sweep's own EEDF moment;"
            << " nothing here is a fit or an assumed ratio." << nl;
    }

    if (fromCase.size())
    {
        Info<< "    READ from the case (used whole, the default dropped): "
            << flatOutput(fromCase) << nl;
    }

    Info<< "    initialMeanEnergy is keyed on `reducedE`, NOT `" << key
        << "`: meanEnergy_vs_meanE is the identity and is never written, and"
        << " the seed precedes " << key << "." << nl;

    // reportInterval: 25 ON THE DEFAULTING PATH ONLY.
    //
    // The model's own default is 1, and that is left alone for any case that
    // writes its own block -- silently changing an existing default is the
    // thing this change exists to make unnecessary.
    //
    // But 1 is the wrong value to HAND a new user, because the per-term dump
    // evaluates ddt/div/laplacian over the whole mesh once per CORRECTOR:
    //   - measured on the 1.15M-cell 2 ns streamer: ~44 s -> ~23 s per step
    //     when set to 25, i.e. reporting every step nearly DOUBLES the run;
    //   - measured 2026-09-01 on the 130x130 bed: only 1.7% (76.94 vs
    //     75.68 s), because at that size the dump is small next to the step.
    // The cost is therefore mesh-dependent and severe exactly where it hurts,
    // so a one-switch case would conclude "LMEA is slow" from a diagnostics
    // default. 25 is what every hand-written block in the repository uses.
    if (!c.found("reportInterval"))
    {
        c.add("reportInterval", 25);
        Info<< "    reportInterval defaulted to 25. The per-term dump is a"
            << " whole-mesh evaluation per corrector: measured ~44 s -> ~23 s"
            << " per step on the 1.15M-cell streamer between 1 and 25." << nl;
    }

    Info<< endl;

    return c;
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
        // THE OUTER `backgroundGas/pressure`, not `backgroundGas/energy/pressure`.
        // The nested key duplicated the outer one in the same file; the outer is
        // what plasmaSpecies closes the gas density from, so it is the owner.
        pGasAtm_    = constant::plasma::atmFromPa
        (
            bgGasEnergyDict.getOrDefault<scalar>
            (
                "pressure",
                species_.backgroundDict().getOrDefault<scalar>
                (
                    "pressure", constant::plasma::PaPerAtm
                )
            )
        );
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

        // ASKED, not derived. plasmaSpecies::resolveElectronEnergyModel() is
        // the one owner of the electron closure; re-reading `energyModel` here
        // is how three sites came to derive the same condition independently,
        // and how a run came to be half-LMEA.
        //
        // Heavy species always follow the gas temperature. That is not a
        // narrowing of anything real: the per-species temperature vocabulary
        // that used to be readable here filled five labelLists with no readers
        // (see the removal note in plasmaSpecies.C), and `gasTemperature` is
        // precisely "use T_gas", which is the solved field when
        // backgroundGas/energy/solve is true and the fixed dictionary value
        // when it is not. One switch, one place.
        const bool isElectron = (i == species_.electronSpeciesID());

        const word modelName =
            isElectron
          ? (species_.isLMEA() ? word("localEnergy") : word("gasTemperature"))
          : word("gasTemperature");

        // Optional for the same reason energyModel is: the default model
        // needs no coefficients, so demanding the sub-dictionary would be a
        // second migration tax behind the first.
        //
        // For the electron under LMEA the block is also RESOLVED: whatever
        // leaves the case omitted are filled in from the model, so the switch
        // alone is a complete configuration. Stored in modelDicts_ because a
        // synthesised dictionary has no home in the case file.
        modelDicts_[i] = resolveEnergyModelCoeffs
        (
            modelName,
            i,
            sDict.subOrEmptyDict("energyModelCoeffs")
        );

        const dictionary& modelDict = modelDicts_[i];

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
    energyModels_(species.nSpecies()),
    modelDicts_(species.nSpecies())
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

        // NOT gated on active() -- enrolling only records the field so the
        // joint Picard residual (and rho) can be formed. See
        // plasmaOuterRelaxation::enrol().
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

