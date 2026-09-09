/*---------------------------------------------------------------------------*\
  File: plasmaSpecies.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaSpecies.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaStepAudit.H"
#include "plasmaSpecies.H"
#include "IFstream.H"
#include "HashSet.H"
#include "plasmaConstants.H"
#include "DynamicList.H"
#include "FlatOutput.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaSpecies, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::plasmaSpecies::readMechanismSpecies()
{
    fromMechanism_ = true;
    speciesNames_ = speciesFromMechanism(*this, &mechCharge_, &mechMass_,
                                        &ionTransport_, &ionFluxScheme_,
                                        &neutralTransport_, &mechDiff_,
                                        &diffTref_, &diffPref_,
                                        &mechDiffExp_);
}


Foam::wordList Foam::plasmaSpecies::speciesFromMechanism
(
    const dictionary& speciesDict,
    HashTable<scalar>* chargeOut,
    HashTable<scalar>* massOut,
    word* ionTrOut,
    word* fluxOut,
    word* neutralTrOut,
    HashTable<scalar>* diffOut,
    scalar* diffTrefOut,
    scalar* diffPrefOut,
    HashTable<scalar>* diffExpOut
)
{
    // Static, and taking only the dictionary, because two callers need it: this
    // class, and plasmaCreateSpeciesFields, which runs before a mesh exists and
    // must create exactly the fields the solver will later look for. When they
    // derived the list separately they disagreed -- the utility did not
    // understand `fromMechanism` at all, so it created nothing and the solver
    // then failed on a missing 0/n_e.
    wordList names;
    HashTable<scalar> mechCharge_;
    HashTable<scalar> mechMass_;

    const dictionary md = speciesDict.subOrEmptyDict("mechanismSpecies");
    const fileName mechFile = md.getOrDefault<fileName>
    (
        "mechanism", "constant/air_plasma.foam"
    );

    IFstream is(mechFile);
    if (!is.good())
    {
        FatalIOErrorInFunction(speciesDict)
            << "activeSpecies is `fromMechanism` but the mechanism dictionary "
            << mechFile << " cannot be read." << nl
            << "    It is written by mechc alongside the .mech.json manifest."
            << " Set `mechanism` in the `mechanismSpecies` dictionary if it"
            << " lives elsewhere." << nl << exit(FatalIOError);
    }
    const dictionary mech(is);

    const dictionary& chargeDict = mech.subDict("speciesCharge");
    const dictionary& massDict   = mech.subDict("speciesMass");

    for (const entry& e : chargeDict)
    {
        mechCharge_.insert(e.keyword(), readScalar(e.stream()));
    }
    for (const entry& e : massDict)
    {
        mechMass_.insert(e.keyword(), readScalar(e.stream()));
    }

    // Which species are SOLVED -- carried as fields and evolved by the
    // chemistry. NOT which are transported: that is `transportModel` per
    // species, with `ionTransport` and `neutralTransport` as bulk defaults.
    // The two questions are separate, and a species can perfectly well be in
    // the chemistry while sitting still.
    //
    //   charged             electron + ions. The minimum for a self-consistent
    //                       discharge: these carry the space charge that drives
    //                       the field, so omitting one is not an approximation,
    //                       it breaks charge conservation.
    //
    //   chargedAndExcited   the above plus excited states and radicals. Needed
    //                       once excited-state chemistry matters -- stepwise
    //                       ionisation, quenching, associative processes.
    //
    //                       A species left OUT is not merely untransported, it
    //                       has DENSITY ZERO: every reaction consuming it is
    //                       dead. With `charged`, the N2(A3,B3,C3,a1) + O2
    //                       quenching that carries fast gas heating in air
    //                       cannot fire at all, however complete the mechanism.
    //
    //                       The cost is one field each plus their chemistry;
    //                       neutrals default to `immobile`, so it is a diagonal
    //                       solve rather than a transport equation unless
    //                       `neutralTransport diffusion` asks for more.
    //
    //   all                 everything the mechanism names, background gas
    //                       included. Rarely wanted: the background is held
    //                       fixed by construction, so transporting it both
    //                       doubles its cost and lets it drift from the density
    //                       the rate tables were built at.
    // How derived IONS are transported.
    //
    //   driftDiffusion   (DEFAULT since 2026-09-04) they drift and diffuse,
    //                    using the mu*N and D*N tables ionmob writes from
    //                    LXCat measurements. Ion motion sets the timescale of
    //                    afterglow, and of a DBD's memory between pulses.
    //   immobile         they carry space charge but do not move. Correct on
    //                    nanosecond timescales -- an ion drifts a few microns
    //                    while an electron crosses the domain -- and what the
    //                    published streamer benchmark assumes.
    //
    // WHY THE DEFAULT CHANGED. `immobile` is right for a nanosecond single
    // pulse and wrong for everything else, and "everything else" is what a
    // user arriving with a DBD or an afterglow problem has. A default that is
    // correct only for the shortest case in the suite reads, to anyone who
    // does not already know, as "ions are handled" -- and the failure is
    // silent: the space charge is simply missing its slower half.
    //
    // The four `positiveStreamer_*` tutorials now PIN `immobile` explicitly,
    // so the validation baseline did not move with this default. A validation
    // case must not inherit a default it depends on.
    const word ionTr = md.getOrDefault<word>("ionTransport", "driftDiffusion");
    if (ionTr != "immobile" && ionTr != "driftDiffusion")
    {
        FatalIOErrorInFunction(speciesDict)
            << "Unknown `ionTransport` '" << ionTr << "' in mechanismSpecies"
            << nl << "Valid: immobile | driftDiffusion" << nl
            << exit(FatalIOError);
    }

    // Default transport for derived NEUTRAL species, in the SAME vocabulary as
    // ionTransport. `immobile` is right for anything whose chemical lifetime is
    // short enough that it is quenched before it moves -- which is every N2
    // electronic state at atmospheric pressure -- and `diffusion` for radicals
    // and metastables that live long enough to spread. The start-up transport
    // advisory reports the diffusion length against the cell size per species,
    // so the choice can be checked rather than assumed.
    const word neutralTr = md.getOrDefault<word>("neutralTransport", "immobile");
    if (neutralTr != "immobile" && neutralTr != "diffusion")
    {
        FatalIOErrorInFunction(speciesDict)
            << "Unknown `neutralTransport` '" << neutralTr
            << "' in mechanismSpecies" << nl
            << "Valid: immobile | diffusion" << nl
            << "    (`driftDiffusion` is not offered: a neutral has no drift.)"
            << nl << exit(FatalIOError);
    }
    if (neutralTrOut) *neutralTrOut = neutralTr;

    // Neutral diffusivities, written by mechc from Lennard-Jones data. Absent
    // in mechanisms compiled before that existed, which is not an error unless
    // a case actually asks for `neutralTransport diffusion`.
    {
        const dictionary& dd = mech.subOrEmptyDict("speciesDiffusivity");
        if (diffTrefOut) *diffTrefOut = dd.getOrDefault<scalar>("Tref", 300.0);
        if (diffPrefOut) *diffPrefOut = dd.getOrDefault<scalar>("pref", 1.0e5);

        // One sub-dictionary per species, each with its OWN exponent: the
        // Chapman-Enskog temperature dependence is T^1.5/Omega_D(T*), and
        // Omega_D differs between species, so a single shared exponent would
        // be wrong for all but one of them.
        for (const entry& e : dd)
        {
            if (!e.isDict()) continue;
            const dictionary& sd = e.dict();
            if (diffOut) diffOut->insert(e.keyword(), sd.get<scalar>("D0"));
            if (diffExpOut)
            {
                diffExpOut->insert
                (
                    e.keyword(), sd.getOrDefault<scalar>("exponent", 1.68)
                );
            }
        }
    }

    const word select = md.getOrDefault<word>("include", "charged");
    if (select != "charged" && select != "chargedAndExcited" && select != "all")
    {
        FatalIOErrorInFunction(speciesDict)
            << "Unknown `include` '" << select << "' in mechanismSpecies" << nl
            << "Valid: charged | chargedAndExcited | all" << nl
            << exit(FatalIOError);
    }

    // Background-gas species are named by the mechanism's own reference
    // composition, so this needs no list of "things that are air".
    wordHashSet background;
    if (mech.found("composition"))
    {
        for (const entry& e : mech.subDict("composition"))
        {
            background.insert(e.keyword());
        }
    }

    const word electronName = mech.getOrDefault<word>("electronSpecies", "Electron");
    const wordList exclude  = md.getOrDefault<wordList>("exclude", wordList());
    wordHashSet excluded(exclude);

    // The electron goes first, because its index is the electron species ID
    // everything else looks up.
    const word caseElectron = md.getOrDefault<word>("electronName", "e");
    names.append(caseElectron);
    mechCharge_.set(caseElectron, -1);
    if (mechMass_.found(electronName))
    {
        mechMass_.set(caseElectron, mechMass_[electronName]);
    }

    for (const word& name : chargeDict.toc())
    {
        if (name == electronName || name == caseElectron) continue;
        if (excluded.found(name)) continue;

        const scalar q = mechCharge_[name];
        const bool charged = (mag(q) > SMALL);
        const bool isBackground = background.found(name);

        bool take = false;
        if (select == "all")                    take = true;
        else if (select == "chargedAndExcited") take = charged || !isBackground;
        else                                    take = charged;

        if (take && !(select != "all" && isBackground && !charged))
        {
            names.append(name);
        }
    }

    Info<< "plasmaSpecies: activeSpecies from " << mechFile
        << " (include " << select << "): " << names << endl;

    if (chargeOut) *chargeOut = mechCharge_;
    if (massOut)   *massOut   = mechMass_;
    if (ionTrOut)  *ionTrOut  = ionTr;
    if (fluxOut)
    {
        // THE ION DEFAULT FOLLOWS THE ELECTRON, it is not "standard".
        //
        // A flux scheme is a statement about how the drift-diffusion flux is
        // discretised, and running the electron on ScharfetterGummel while
        // every ion stays on the split div/laplacian is a mixed
        // discretisation nobody asked for: `fluxScheme ScharfetterGummel` in
        // the case reads as "use SG", and before this it silently switched
        // ONLY the electron, because the derived ions are synthesised from
        // this key and it was hard-defaulted. An explicit `ionFluxScheme` in
        // mechanismSpecies still wins, for the case that means it.
        word eScheme("standard");

        if (speciesDict.found("speciesProperties"))
        {
            const dictionary& sp = speciesDict.subDict("speciesProperties");

            if (sp.found(caseElectron))
            {
                const dictionary& ed = sp.subDict(caseElectron);

                if (ed.found("driftDiffusionCoeffs"))
                {
                    eScheme = ed.subDict("driftDiffusionCoeffs")
                        .lookupOrDefault<word>("fluxScheme", "standard");
                }
            }
        }

        *fluxOut = md.getOrDefault<word>("ionFluxScheme", eScheme);

        if (*fluxOut != eScheme)
        {
            Info<< "plasmaSpecies: ion fluxScheme `" << *fluxOut
                << "` OVERRIDES the electron's `" << eScheme
                << "` (ionFluxScheme set explicitly)." << endl;
        }
    }
    return names;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaSpecies::plasmaSpecies
(
    const fvMesh& mesh,
    electromagneticsModel& em
)
:
    IOdictionary
    (
        IOobject
        (
            "plasmaSpeciesProperties",
            mesh.time().constant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    em_(em),
    nSpecies_(0),
    speciesNames_(),
    speciesIDs_(),
    speciesMasses_(),
    speciesCharges_(),
    speciesChargeNumbers_(),
    numberDensities_(),
    speciesMinNumberDensities_(),
    speciesDicts_(),
    defaultSpeciesDict_(),
    backgroundName_("none"),
    backgroundDensity_
    (
        "backgroundDensity", 
        dimensionSet(0, -3, 0, 0, 0, 0, 0),
        0.0
    ),
    totalNeutralDensity_
    (
        IOobject
        (
            "totalNeutralDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero", 
            dimensionSet(0, -3, 0, 0, 0, 0, 0), 
            0.0
        )
    ),
    backgroundDict_(),
    // Overwritten by resolveElectronEnergyModel(); LFA is the safe value to
    // hold in the window before it runs, since it demands nothing extra.
    electronEnergyModel_("LFA"),
    electronSpeciesID_(-1),
    ionSpeciesIDs_(),
    positiveIonSpeciesIDs_(),
    negativeIonSpeciesIDs_(),
    neutralSpeciesIDs_(),
    chargedSpeciesIDs_(),
    mobileSpeciesIDs_(),
    immobileSpeciesIDs_(),
    reactingSpeciesIDs_(),
    nonReactingSpeciesIDs_()
{
    if (!found("backgroundGas"))
    {
        FatalIOErrorInFunction(*this)
            << "Missing required dictionary 'backgroundGas' in "
            << objectPath() << nl << exit(FatalIOError);
    }
    
    const dictionary& bgDict = subDict("backgroundGas");
    backgroundName_ = bgDict.get<word>("name");
    // Background number density: either stated, or closed by the ideal gas law
    // from pressure and temperature.
    //
    // The tutorial hard-coded 2.4463e25, which IS p/(k_B T) at 1 atm and 300 K
    // -- but written as a literal it silently stops being that the moment the
    // case temperature or pressure changes, and nothing connects it to the
    // `energy { T ... }` sitting three lines below in the same dictionary.
    // Deriving it means the mechanism tables, the gas temperature and the
    // background density cannot disagree about what gas is being modelled.
    //
    // `numberDensity` still wins when given, because a case may be
    // deliberately non-ideal or matching a reference calculation.
    if (bgDict.found("numberDensity"))
    {
        backgroundDensity_.value() = bgDict.get<scalar>("numberDensity");
    }
    else
    {
        const scalar pAbs = bgDict.getOrDefault<scalar>("pressure", 101325.0);
        const dictionary eDict = bgDict.subOrEmptyDict("energy");
        const scalar Tgas = eDict.getOrDefault<scalar>("T", 300.0);

        if (pAbs <= 0 || Tgas <= 0)
        {
            FatalIOErrorInFunction(*this)
                << "backgroundGas: pressure and temperature must be positive"
                << " to close the ideal gas law (got p = " << pAbs
                << " Pa, T = " << Tgas << " K)" << nl << exit(FatalIOError);
        }

        backgroundDensity_.value() = pAbs/(constant::plasma::kappaBoltzmann.value()*Tgas);

        Info<< "    background gas: N = p/(k_B T) = "
            << backgroundDensity_.value() << " 1/m3"
            << "  (p = " << pAbs << " Pa, T = " << Tgas << " K)" << endl;
    }
    backgroundDict_ = bgDict;

    totalNeutralDensity_ == backgroundDensity_;

    // One owner for the gas density. The electromagnetics model needs it to
    // form reducedE = |E|/N, and reading its own copy is how E/N came to be
    // evaluated at 1 atm while the case was at 1 bar.
    em_.setBackgroundDensity(backgroundDensity_);

    // backgroundGas/energy is NOT read here any more. It is a global property
    // and plasmaEnergy is its one owner (see the removal note further down);
    // reading `solve`/`T` here as well is how a second, disagreeing copy of a
    // setting gets established.

    if (!found("activeSpecies"))
    {
        FatalIOErrorInFunction(*this)
            << "Required entry 'activeSpecies' is missing in dictionary "
            << objectPath() << nl << exit(FatalIOError);
    }

    // Read species list
    // `activeSpecies fromMechanism;` derives the list from the compiled
    // mechanism instead of repeating it by hand. The mechanism already knows
    // which species its reactions create; writing them out again is how a case
    // ends up transporting an ion the chemistry does not produce, or -- far
    // worse and the reason this exists -- NOT transporting one that it does,
    // which silently breaks charge conservation.
    //
    // An explicit list still works and still wins. Deriving is a convenience,
    // not a policy: a case may deliberately carry a subset.
    {
        ITstream& is = lookup("activeSpecies");
        token firstToken(is);
        is.rewind();

        if (firstToken.isWord() && firstToken.wordToken() == "fromMechanism")
        {
            readMechanismSpecies();
        }
        else
        {
            is >> speciesNames_;
        }
    }
    nSpecies_ = speciesNames_.size();

    speciesChargeNumbers_.setSize(nSpecies_);
    speciesCharges_.setSize(nSpecies_);
    speciesMasses_.setSize(nSpecies_);
    numberDensities_.setSize(nSpecies_);
    speciesMinNumberDensities_.setSize(nSpecies_);
    speciesDicts_.resize(nSpecies_);

    // Read species properties dictionary
    if (!found("speciesProperties"))
    {
        FatalIOErrorInFunction(*this)
            << "Missing required dictionary 'speciesProperties' in "
            << objectPath() << nl << exit(FatalIOError);
    }

    const dictionary& propsDict = subDict("speciesProperties");

    // Read defaultProperties if present
    defaultSpeciesDict_ = propsDict.subOrEmptyDict("defaultProperties");

    // Loop over all species
    for (label i = 0; i < nSpecies_; ++i)
    {
        const word& sName = speciesNames_[i];
        speciesIDs_.insert(sName, i);

        // A derived species need not have its own sub-dictionary: the
        // defaults plus the mechanism's charge and mass are enough to carry
        // it. An explicitly listed species still must, so a typo in a
        // hand-written activeSpecies list is still caught.
        if (!propsDict.found(sName) && !fromMechanism_)
        {
            FatalIOErrorInFunction(*this)
                << "Species '" << sName << "' is listed in 'activeSpecies' but "
                << "has no sub-dictionary in " << objectPath() << nl
                << exit(FatalIOError);
        }

        // Build merged properties (defaults + overrides)
        dictionary mergedDict(defaultSpeciesDict_);
        if (propsDict.found(sName))
        {
            mergedDict.merge(propsDict.subDict(sName));
        }

        // `ionTransport driftDiffusion` gives every derived ION a
        // drift-diffusion model reading the tables ionmob generated, without
        // the case naming a single one of them. The point of deriving species
        // from the mechanism is that the case does not restate what the
        // mechanism already knows; making the user then hand-write a
        // driftDiffusionCoeffs block per ion would give that back.
        //
        // A species with its own sub-dictionary still wins, so one ion can be
        // treated differently without opting the rest out.
        // Gated on `transportModel` specifically, not on whether the species
        // has a block at all: the ions have blocks that set only a floor
        // density, and treating any block as an override made
        // `ionTransport driftDiffusion` a silent no-op.
        // Derived species need a transport model, and `ionTransport` is what
        // decides it. Handling only the driftDiffusion branch left the default
        // (`immobile`) case with no transportModel at all, which fails at
        // construction -- the switch has to answer for both of its values.
        // Derived NEUTRALS take neutralTransport. Done before the ion branch
        // so that `ionTransport driftDiffusion` cannot reach a neutral, which
        // has no mobility to drift with.
        const bool isNeutral =
            mag(mechCharge_.lookup(sName, 0.0)) <= SMALL
         && sName != speciesNames_[0];

        if (fromMechanism_
         && !mergedDict.found("transportModel")
         && isNeutral
         && neutralTransport_ == "diffusion")
        {
            if (!mechDiff_.found(sName))
            {
                FatalErrorInFunction
                    << "`neutralTransport diffusion` needs a diffusivity for '"
                    << sName << "', and the mechanism does not carry one."
                    << nl << nl
                    << "    mechc writes `speciesDiffusivity` from"
                    << " Lennard-Jones data; a mechanism compiled" << nl
                    << "    before that existed has no such block."
                    << " Recompile it, or give the species its own" << nl
                    << "    `transportModel diffusion` with an explicit"
                    << " `diffusivity` sub-dictionary." << nl
                    << exit(FatalError);
            }

            // D(T) = D0 (T/Tref)^e, written in the powerLaw evaluator's form,
            // amplitude*var^exponent, so amplitude absorbs Tref^-e.
            const scalar D0 = mechDiff_[sName];
            const scalar expo = mechDiffExp_.lookup(sName, 1.68);

            // T_gas exists as a FIELD only when the energy equation is solved.
            // With heating off there is nothing for a powerLaw to look up, and
            // asking for one dies at construction on a missing registry entry
            // -- so the temperature dependence is folded into a constant at the
            // dictionary temperature instead. Same D either way; the difference
            // is only whether it can follow a temperature that moves.
            const dictionary bgD = subDict("backgroundGas");
            const dictionary eD = bgD.subOrEmptyDict("energy");
            const bool solvesT = eD.getOrDefault<bool>("solve", false);
            const scalar Tfix = eD.getOrDefault<scalar>("T", 300.0);

            dictionary dif;
            if (solvesT)
            {
                dif.add("type", word("powerLaw"));
                dif.add("amplitude", D0/Foam::pow(diffTref_, expo));
                dif.add("exponent", expo);
                dif.add("lookupVariable", word("T_gas"));
            }
            else
            {
                dif.add("type", word("constant"));
                dif.add
                (
                    "value",
                    D0*Foam::pow(Tfix/diffTref_, expo)
                );
            }

            dictionary dc;
            dc.add("diffusivity", dif);

            mergedDict.add("transportModel", word("diffusion"));
            mergedDict.add("diffusionCoeffs", dc);
        }

        // THE `immobile` VALUE OF EACH SWITCH, kept SEPARATE, because they are
        // separate switches over separate sets of species.
        //
        // This was ONE blanket default gated on `ionTransport_ !=
        // "driftDiffusion"`, which meant that asking for MOBILE IONS silently
        // switched off the default for NEUTRALS -- and since the driftDiffusion
        // branch below only ever reaches charged species, every neutral was
        // left with no transportModel at all. Measured 2026-09-04 on needleDBD,
        // the first case to want `ionTransport driftDiffusion` together with
        // `neutralTransport immobile`:
        //     "Species 'N' is missing required entry 'transportModel'"
        //
        // The comment above already records fixing exactly this shape of bug
        // for the ion switch. The same mistake was still here for the neutral
        // one: a switch must answer for both of ITS values, and must not
        // answer for another switch's.
        if (fromMechanism_
         && !mergedDict.found("transportModel")
         && isNeutral
         && neutralTransport_ == "immobile")
        {
            mergedDict.add("transportModel", word("immobile"));
        }

        if (fromMechanism_
         && !mergedDict.found("transportModel")
         && !isNeutral
         && sName != speciesNames_[0]
         && ionTransport_ == "immobile")
        {
            mergedDict.add("transportModel", word("immobile"));
        }

        if (fromMechanism_
         && ionTransport_ == "driftDiffusion"
         && mag(mechCharge_.lookup(sName, 0.0)) > SMALL
         && sName != speciesNames_[0]
         && !mergedDict.found("transportModel"))
        {
            dictionary dd;
            dd.add("fluxScheme", ionFluxScheme_);

            // WHERE THE ION TABLES LIVE.
            //
            // Not the same directory as the electron ones, and that is
            // structural rather than a preference: the electron tables come
            // from the Boltzmann sweep into `constant/plasmaTables`, while ion
            // mobilities come from a different chain entirely (LXCat Viehland
            // data through tools/ionmob.py) into `constant/ionTables`.
            //
            // The derived block used to carry NO tableDir, so it fell back to
            // the electron directory and died with
            //   `fromMechanism needs "constant/plasmaTables/muN_N2p_vs_reducedE"`
            // -- measured 2026-09-04 on needleDBD, the first case to derive its
            // ion transport at all. `ionTableDir` overrides it for a case that
            // keeps them elsewhere.
            const word ionTabDir
            (
                subDict("mechanismSpecies").getOrDefault<word>
                (
                    "ionTableDir", "constant/ionTables"
                )
            );

            dictionary mu;
            mu.add("type", word("fromMechanism"));
            mu.add("quantity", word("muN_" + sName));
            mu.add("tableDir", ionTabDir);
            dd.add("mobility", mu);

            dictionary dif;
            dif.add("type", word("fromMechanism"));
            dif.add("quantity", word("DLN_" + sName));
            dif.add("tableDir", ionTabDir);
            dd.add("diffusivity", dif);

            mergedDict.add("transportModel", word("driftDiffusion"));
            mergedDict.add("driftDiffusionCoeffs", dd);
        }
        else if (fromMechanism_
              && ionTransport_ == "driftDiffusion"
              && mag(mechCharge_.lookup(sName, 0.0)) > SMALL
              && sName != speciesNames_[0])
        {
            // Asked for mobile ions and got an override. Legitimate -- one ion
            // may need different treatment -- but silence here would let a
            // case believe its ions move when they do not.
            WarningInFunction
                << "`ionTransport driftDiffusion` is set, but species '"
                << sName << "' declares its own transportModel ("
                << mergedDict.get<word>("transportModel")
                << "), which takes precedence." << endl;
        }
        speciesDicts_.insert(sName, mergedDict);

        // Charge and mass: from the case if stated, otherwise from the
        // mechanism. Both are properties of the SPECIES, not of the
        // simulation, and the mechanism computes them from the elemental
        // composition -- so a case that omits them cannot get them wrong.
        //
        // The tutorial that motivated this carried `mass 1.67e-26` for its
        // generic positive ion. That is the proton mass; N2+ is 4.65e-26, so
        // it was out by a factor of 2.8. Nothing detected it, because nothing
        // else in the case knew what the ion was supposed to be.
        if (!mergedDict.found("charge") && mechCharge_.found(sName))
        {
            mergedDict.add("charge", mechCharge_[sName]);
        }
        if (!mergedDict.found("mass") && mechMass_.found(sName))
        {
            mergedDict.add("mass", mechMass_[sName]);
        }

        if (!mergedDict.found("charge"))
        {
            FatalIOErrorInFunction(*this)
                << "Species '" << sName << "' is missing required entry "
                << "'charge' in " << objectPath() << nl
                << "    It can also come from the mechanism: give"
                << " `mechanism` in the `mechanismSpecies` dictionary." << nl
                << exit(FatalIOError);
        }

        if (!mergedDict.found("mass"))
        {
            FatalIOErrorInFunction(*this)
                << "Species '" << sName << "' is missing required entry "
                << "'mass' in " << objectPath() << nl
                << "    It can also come from the mechanism: give"
                << " `mechanism` in the `mechanismSpecies` dictionary." << nl
                << exit(FatalIOError);
        }

        speciesChargeNumbers_[i] = readScalar(mergedDict.lookup("charge"));
        scalar massValue = readScalar(mergedDict.lookup("mass"));

        speciesCharges_.set
        (
            i,
            new dimensionedScalar
            (
                "q_" + sName,
                constant::plasma::eCharge.dimensions(),
                speciesChargeNumbers_[i] * constant::plasma::eCharge.value()
            )
        );

        speciesMasses_.set
        (
            i,
            new dimensionedScalar
            (
                "m_" + sName,
                constant::plasma::eMass.dimensions(),
                massValue
            )
        );

        numberDensities_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "n_" + sName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

        speciesMinNumberDensities_[i] = 
                       mergedDict.getOrDefault<scalar>("minNumberDensity", 0.0);

        // Groups
        scalar Z = speciesChargeNumbers_[i];
        if (Z == -1 && massValue < 1e-29)
        {
            electronSpeciesID_ = i;
            chargedSpeciesIDs_.append(i);
        }

        // Ions (Charged, but not electrons)
        else if (Z != 0)
        {
            chargedSpeciesIDs_.append(i);
            ionSpeciesIDs_.append(i); 
            if (Z > 0)
                positiveIonSpeciesIDs_.append(i);
            else
                negativeIonSpeciesIDs_.append(i);
        }
        // Active neutrals
        else
        {
            neutralSpeciesIDs_.append(i);
        }
            
        // Mobile vs immobile species
        const word transport =
            mergedDict.getOrDefault<word>("transportModel", "immobile");

        if (transport == "immobile")
        {
            immobileSpeciesIDs_.append(i);
        }
        else
        {
            mobileSpeciesIDs_.append(i);
        }

        // NO per-species temperature grouping here.
        //
        // REMOVED 2026-09-01, with the five labelLists it filled. It read
        // `energyModel` with its own vocabulary -- isothermal (default) /
        // backgroundGas / localField / solveEnergy -- while plasmaEnergy and
        // plasmaTransport read the SAME key as {gasTemperature (default),
        // isothermal, localEnergy, localField}. Two enums, one key name: only
        // two values were shared, `gasTemperature` matched nothing here, and
        // `backgroundGas`/`solveEnergy` would have been fatal as
        // plasmaEnergyModel type names.
        //
        // It was also inert -- all five accessors had zero callers -- so the
        // whole block decided nothing while presenting the user with four
        // keywords to choose between. The gas temperature is a GLOBAL
        // property: `backgroundGas/energy/solve` selects whether it is solved
        // and `backgroundGas/energy/T` fixes its value, in one place. The
        // electron closure is `electronEnergyModel` (LFA|LMEA), which is
        // orthogonal to gas heating -- all four combinations are supported.
    }

    // SECOND PASS, after the loop: electronSpeciesID_ is only assigned inside
    // it, and the legacy per-species spelling lives on the electron's own
    // sub-dictionary.
    resolveElectronEnergyModel();
}


void Foam::plasmaSpecies::resolveElectronEnergyModel()
{
    // ONE OWNER for the electron energy closure.
    //
    // WHY A GLOBAL KEY RATHER THAN A PER-SPECIES ONE (user, 2026-09-01):
    // "gasTemperature for electrons does not mean anything to me and the user
    // won't know if that's LFA". Two independent problems with the old shape:
    //
    //   1. LMEA is meaningful for EXACTLY ONE species -- meanE, muEpsN and
    //      PelasticN are all electron quantities -- so a per-species key
    //      invited `energyModel localEnergy` on an ion, which means nothing.
    //   2. It spelled the closure in a vocabulary that never mentions LFA or
    //      LMEA, so the switch could not be found by the name it is known by.
    //
    // The key is also read in ONE place because it was previously derived
    // independently in three (plasmaEnergy twice, plasmaTransport once), which
    // is how a run came to be half-LMEA: transport keyed on the mean energy
    // while ionisation still followed the field, a measured runaway.
    const bool haveGlobal = found("electronEnergyModel");

    // NO DEFAULT. The key is REQUIRED (user's decision, 2026-09-01).
    //
    // WHY THIS ONE IS NOT DEFAULTED, when so much else in this change is:
    // LFA and LMEA are DIFFERENT PHYSICS that give different answers, so there
    // is no value that is merely "the historical one". The migration-tax
    // argument that justifies defaulting `energyModel` or `energyModelCoeffs`
    // -- an entry whose only sensible value is the historical one is a tax, not
    // a safety feature -- turns on there being ONE sensible value. Here there
    // are two, and picking either silently chooses the closure for the user.
    //
    // That is the same failure this whole change exists to remove: the
    // half-LMEA is fatal rather than warned precisely because quietly-wrong
    // physics, with plausible output and no error, is the worst outcome.
    //
    // The one-switch win is that ONE key configures ~80 lines of coefficients,
    // not that the key itself can be omitted. So the error below TEACHES the
    // choice rather than making it.
    word resolved;

    if (haveGlobal)
    {
        resolved = get<word>("electronEnergyModel");

        if (resolved != "LFA" && resolved != "LMEA")
        {
            FatalIOErrorInFunction(*this)
                << "electronEnergyModel is `" << resolved
                << "`, which is not a closure this solver has." << nl
                << "    The only values are:" << nl
                << "      LFA   -- local field approximation: electron"
                << " coefficients and rates follow the local E/N." << nl
                << "      LMEA  -- local mean energy approximation: the"
                << " electron energy-density equation is solved and"
                << " coefficients and rates follow the local mean energy."
                << nl
                << exit(FatalIOError);
        }
    }

    // LEGACY per-species spelling, accepted for one release.
    //
    // Mapped, not ignored: silently dropping it would turn every existing LMEA
    // case into an LFA case that still runs and still looks plausible, which is
    // precisely the failure this whole change exists to remove.
    word legacy(word::null);

    if (electronSpeciesID_ >= 0)
    {
        const dictionary& eDict = speciesDict(electronSpeciesID_);

        if (eDict.found("energyModel"))
        {
            const word em = eDict.get<word>("energyModel");

            if (em == "localEnergy")
            {
                legacy = "LMEA";
            }
            else if (em == "gasTemperature" || em == "isothermal")
            {
                legacy = "LFA";
            }
            else
            {
                FatalIOErrorInFunction(*this)
                    << "The electron species carries the legacy key"
                    << " `energyModel " << em << "`, which does not map onto"
                    << " an electron energy closure." << nl
                    << "    Replace it with the top-level key"
                    << " `electronEnergyModel LFA;` or"
                    << " `electronEnergyModel LMEA;`." << nl
                    << exit(FatalIOError);
            }
        }
    }

    // THE PER-SPECIES SPELLING IS REJECTED, not translated.
    //
    // Keeping it readable would leave two live spellings of one setting, which
    // is how the two-vocabulary mess this change removes came about in the
    // first place. The translation is unambiguous and printed, so the fix is a
    // one-line edit.
    if (!legacy.empty())
    {
        FatalIOErrorInFunction(*this)
            << "The electron species carries the REMOVED per-species key"
            << " `energyModel`." << nl
            << "    The electron energy closure is now a single TOP-LEVEL key."
            << nl
            << nl
            << "    Replace it with:" << nl
            << nl
            << "        electronEnergyModel " << legacy << ";" << nl
            << nl
            << "    at the top level of plasmaSpeciesProperties, and delete the"
            << " `energyModel` entry from the electron's block."
            << (haveGlobal
                  ? " (You already set electronEnergyModel; the two spellings"
                    " of one setting must not coexist.)"
                  : "")
            << nl
            << exit(FatalIOError);
    }

    // The key is REQUIRED, so there is no defaulted-LMEA case to rescue and no
    // fallback here. An explicit `LMEA` on a mechanism with no
    // mean-energy-keyed tables reaches plasmaEnergy's pre-flight, which gives
    // the full diagnosis (non-invertible mean energy in an
    // attachment-dominated mixture) instead of quietly running a different
    // closure than the one that was asked for.
    if (resolved.empty())
    {
        FatalIOErrorInFunction(*this)
            << "`electronEnergyModel` is not set, and it has no default." << nl
            << "    It selects the ELECTRON ENERGY CLOSURE, which is a physics"
            << " choice: LFA and LMEA give different answers, so neither can be"
            << " assumed on your behalf." << nl
            << nl
            << "    Add ONE of these at the TOP LEVEL of"
            << " constant/plasmaSpeciesProperties:" << nl
            << nl
            << "      electronEnergyModel LFA;   // local FIELD approximation."
            << " Electron coefficients and reaction rates are read at the local"
            << " reduced field E/N. Assumes the electron energy distribution is"
            << " in equilibrium with the local field. Cheaper, and the standard"
            << " choice for streamer work where the field varies slowly"
            << " compared with the energy relaxation length." << nl
            << nl
            << "      electronEnergyModel LMEA;  // local MEAN ENERGY"
            << " approximation. Solves an electron energy-density equation and"
            << " reads the coefficients at the local mean energy instead, so the"
            << " energy is transported and may lag the field. More accurate"
            << " where LFA breaks -- in a streamer head and near boundaries --"
            << " at the cost of one more equation." << nl
            << nl
            << "    If you are unsure, LMEA is the better physics and needs no"
            << " further configuration: it derives the electron transport key,"
            << " the chemistry rate key and the whole energyModelCoeffs block"
            << " from the mechanism's own tables." << nl
            << "    See docs/models/energy/lmea.md." << nl
            << nl
            << "    NOTE this is separate from gas heating"
            << " (backgroundGas/energy/solve). The two are independent and all"
            << " four combinations are supported." << nl
            << exit(FatalIOError);
    }

    electronEnergyModel_ = resolved;

    // Heavy species: report that a legacy energy key does nothing, rather than
    // dropping it silently. The five temperature groupings it used to fill had
    // no readers at all (see the removal note above), so a case carrying one
    // has always been configuring nothing -- and a user is entitled to know
    // that the entry they wrote is not doing what they assumed.
    DynamicList<word> staleKeys;

    forAll(speciesNames_, i)
    {
        if (i == electronSpeciesID_) continue;

        if (speciesDict(i).found("energyModel"))
        {
            staleKeys.append(speciesNames_[i]);
        }
    }

    if (staleKeys.size())
    {
        FatalIOErrorInFunction(*this)
            << "The REMOVED key `energyModel` is set on " << staleKeys.size()
            << " heavy species: " << flatOutput(staleKeys) << "." << nl
            << "    Delete these entries. Heavy-species temperature is not a"
            << " per-species setting: `backgroundGas/energy/solve` decides"
            << " whether the gas temperature is solved, and"
            << " `backgroundGas/energy/T` fixes its value." << nl
            << "    Rejected rather than ignored because these entries NEVER"
            << " had a reader -- they filled five index lists nothing consumed"
            << " -- so leaving them in place would let a case go on believing"
            << " it had configured something." << nl
            << exit(FatalIOError);
    }

    Info<< "plasmaSpecies: electronEnergyModel " << electronEnergyModel_
        << " -- electron coefficients and reaction rates are keyed on `"
        << electronLookupKey() << "`." << nl
        << "    Independent of gas heating: `backgroundGas/energy/solve`"
        << " is a separate switch and all four combinations are supported."
        << endl;

    deriveElectronTransportKey();
}


void Foam::plasmaSpecies::deriveElectronTransportKey()
{
    if (electronSpeciesID_ < 0)
    {
        return;
    }

    const word& eName = speciesNames_[electronSpeciesID_];

    // Nothing to point anywhere. A drift-diffusion electron with no
    // coefficients fails loudly in driftDiffusion itself, which is the right
    // place for that error; inventing a block here would hide it.
    if (!speciesDicts_[eName].found("driftDiffusionCoeffs"))
    {
        return;
    }

    dictionary& dd = speciesDicts_[eName].subDict("driftDiffusionCoeffs");
    const word& want = electronLookupKey();

    const wordList props({word("mobility"), word("diffusivity")});

    DynamicList<word> derived;
    DynamicList<word> ownKey;

    for (const word& prop : props)
    {
        if (!dd.found(prop))
        {
            continue;
        }

        dictionary& pd = dd.subDict(prop);
        const word type = pd.getOrDefault<word>("type", "fromMechanism");

        // ONLY `fromMechanism` is derived-and-validated.
        //
        // For that type the key names the TABLE FILE the sweep wrote
        // (<quantity>_vs_<key>), so LFA-versus-LMEA is precisely what the entry
        // selects, and a mismatch is the measured half-LMEA.
        //
        // For a powerLaw / function1 / tabulated1D fit the key is part of the
        // FIT's own definition -- the hand-fitted electron mobility kept in the
        // streamer tutorial is `powerLaw` in |E|, which is a legitimate model
        // and has no meanE-keyed counterpart. Forcing `meanE` onto it would
        // reject valid physics, so those are reported and left alone. Same
        // reasoning as `tableKey` in plasmaTransport: derive what is
        // unambiguously implied, default or leave what is a modelling choice.
        if (type != "fromMechanism")
        {
            ownKey.append
            (
                prop + " (" + type + ", keyed on "
              + pd.getOrDefault<word>("lookupVariable", "<unset>") + ")"
            );
            continue;
        }

        if (pd.found("lookupVariable"))
        {
            const word given = pd.get<word>("lookupVariable");

            if (given != want)
            {
                FatalIOErrorInFunction(*this)
                    << "Electron driftDiffusionCoeffs/" << prop
                    << "/lookupVariable is `" << given
                    << "` but electronEnergyModel is `" << electronEnergyModel_
                    << "`, which requires `" << want << "`." << nl
                    << "    That combination is the HALF-LMEA: the electron"
                    << " flux responding to one variable while the reaction"
                    << " rates respond to another. It is a measured runaway,"
                    << " and it produces no error of its own." << nl
                    << "    The key is DERIVED from electronEnergyModel."
                    << " Delete it from the electron's driftDiffusionCoeffs."
                    << nl
                    << exit(FatalIOError);
            }
        }
        else
        {
            // MechanismProperty's own default is `reducedE`, so an omitted key
            // under LMEA silently read the LFA table. Insert it.
            pd.add("lookupVariable", want);
            derived.append(prop);
        }
    }

    if (derived.size())
    {
        Info<< "plasmaSpecies: electron driftDiffusionCoeffs "
            << flatOutput(derived) << " keyed on `" << want
            << "`, DERIVED from electronEnergyModel." << endl;
    }

    if (ownKey.size())
    {
        Info<< "plasmaSpecies: electron " << flatOutput(ownKey)
            << " is not table-based, so its lookup variable is part of the fit"
            << " and was left as written -- it does NOT follow"
            << " electronEnergyModel." << endl;
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void Foam::plasmaSpecies::verifyChargeDensity(const label everyN) const
{
    if (everyN <= 0) return;
    if (mesh_.time().timeIndex() % everyN != 0) return;

    // Recompute the invariant's right-hand side from the CURRENT densities.
    scalarField rhs(mesh_.nCells(), Zero);

    forAll(chargedSpeciesIDs_, i)
    {
        const label id = chargedSpeciesIDs_[i];
        rhs += numberDensities_[id].primitiveField()
             * speciesCharges_[id].value();
    }

    const scalarField& rho = em_.chargeDensity().primitiveField();

    // A RELATIVE test against the field's own scale, with an absolute floor.
    //
    // The floor matters: SI charge densities here are ~1e-8 C/m^3, which is
    // far below any sensible absolute tolerance, so a purely absolute test
    // passes on ANY input -- the trap that has caught a diagnostic here
    // before. The scale is the max |rhs| over the domain, so the test means
    // the same thing before and after breakdown.
    const scalar scale = max(gMax(mag(rhs)), SMALL);
    const scalar tol = 1e-8;

    scalar worst = 0;
    label worstCell = -1;

    forAll(rho, celli)
    {
        const scalar err = mag(rho[celli] - rhs[celli])/scale;
        if (err > worst) { worst = err; worstCell = celli; }
    }

    if (returnReduce(worst, maxOp<scalar>()) > tol)
    {
        const label nCharged = chargedSpeciesIDs_.size();

        FatalErrorInFunction
            << "CHARGE DENSITY IS NOT sum_i q_i n_i." << nl << nl
            << "    Poisson's source disagrees with the species densities it"
            << " is supposed to be built" << nl
            << "    from, by " << worst << " of the domain scale ("
            << scale << " C/m^3) at cell " << worstCell << ":" << nl << nl
            << "        chargeDensity   = "
            << (worstCell >= 0 ? rho[worstCell] : 0.0) << nl
            << "        sum_i q_i n_i   = "
            << (worstCell >= 0 ? rhs[worstCell] : 0.0) << nl
            << "        charged species = " << nCharged << nl << nl
            << "    A STALE OR WRONG SOURCE MAKES THE FIELD LOOK LIKE"
            << " PHYSICS. Measured 2026-09-06:" << nl
            << "    with the source frozen, E/N was 1104 Td against V/L ="
            << " 1106 Td -- the VACUUM field" << nl
            << "    to 0.2% -- while the ion density grew three decades, and"
            << " nothing screened." << nl << nl
            << "    Candidate causes, in the order worth checking:" << nl
            << "      * updateChargeDensity() skipped this step -- see the"
            << " plasmaStepAudit report" << nl
            << "      * a density CLAMPED or floored after the charge update"
            << nl
            << "      * a charged species missing from chargedSpeciesIDs_"
            << " (there are " << nCharged << ")" << nl
            << "      * a sign or unit error in speciesCharges_" << nl << nl
            << "    See CLAUDE.md rule 27 and doc/case-monitoring-plan.md."
            << nl << exit(FatalError);
    }
}


void Foam::plasmaSpecies::updateChargeDensity()
{
    // Audited: this ran on 0.04% of steps on 2026-09-06 and froze Poisson's
    // source. See plasmaStepAudit.H and CLAUDE.md rule 27.
    plasmaStepAudit::record("updateChargeDensity");

    em_.chargeDensity() == dimensionedScalar
                                        (em_.chargeDensity().dimensions(), 0.0);

    forAll(chargedSpeciesIDs_, i)
    {
        const label id = chargedSpeciesIDs_[i];
        
        em_.chargeDensity() += numberDensities_[id] * speciesCharges_[id];
    }

    em_.chargeDensity().correctBoundaryConditions();

    Info << "Charge density updated." << endl;
}

void plasmaSpecies::clampNumberDensities()
{
    forAll(numberDensities_, i)
    {
        clampNumberDensity(i);
        numberDensities_[i].correctBoundaryConditions();
    }

    Info << "Species' number densities clamped." << endl;
}

void plasmaSpecies::clampNumberDensity(const label i)
{
    if (speciesMinNumberDensities_[i] > 0.0)
    {
        volScalarField& n = numberDensities_[i];
        const dimensionedScalar nMin
        (
            "nMin",
            n.dimensions(),
            speciesMinNumberDensities_[i]
        );

        // ---- POSITIVITY DIAGNOSTIC, measured BEFORE the clamp -------------
        //
        // Scharfetter-Gummel and CompleteFlux carry no limiter, so both CAN
        // produce a negative density in principle: SG's matrix is an M-matrix
        // and is provably positive, but CFS adds an explicit inhomogeneous
        // flux Gamma^i that can be negative. Whether that ACTUALLY happens is
        // an empirical question this reports on, rather than being inferred
        // from the timestep-floor counter (which is a different thing).
        //
        // THE REDUCTIONS ARE UNCONDITIONAL -- a rank whose subdomain happens
        // to be clean must still take part or the collective deadlocks
        // (rule 31).
        // FILE-LOCAL, NOT CLASS MEMBERS -- and deliberately so.
        //
        // Adding members to plasmaSpecies.H changes the CLASS LAYOUT, and
        // plasmaSpecies is included across the tree (plasmaEnergy::required
        // takes a plasmaSpecies const&). On 2026-09-08 doing exactly that
        // segfaulted every case: wmake rebuilt the solver but left
        // libplasmaEnergy compiled against the old layout, so it read garbage.
        // A DIAGNOSTIC MUST NOT CHANGE A WIDELY-INCLUDED CLASS'S ABI. There is
        // one plasmaSpecies per process, so file-local state is equivalent
        // here and costs nothing.
        static List<scalar> negWorst_;
        static List<label>  negCount_;

        if (negWorst_.size() != numberDensities_.size())
        {
            negWorst_.setSize(numberDensities_.size(), 0.0);
            negCount_.setSize(numberDensities_.size(), 0);
        }

        const scalarField& nI = n.primitiveField();

        scalar preMin = GREAT;
        label nNeg = 0;

        forAll(nI, c)
        {
            preMin = Foam::min(preMin, nI[c]);
            if (nI[c] < 0.0) ++nNeg;
        }

        reduce(preMin, minOp<scalar>());
        reduce(nNeg, sumOp<label>());

        if (nNeg > 0)
        {
            negCount_[i] += nNeg;

            // Report only on a new worst excursion, so a persistent problem
            // does not flood the log.
            if (preMin < negWorst_[i])
            {
                negWorst_[i] = preMin;

                WarningInFunction
                    << "NEGATIVE number density clipped for species '"
                    << speciesNames_[i] << "': min = " << preMin
                    << " in " << nNeg << " cell(s), floor is "
                    << speciesMinNumberDensities_[i] << "." << nl
                    << "    POSITIVITY IS VIOLATED. The flux scheme in use"
                       " carries no limiter, so this is possible by"
                    << nl
                    << "    construction. Cumulative negative clips for this"
                       " species: " << negCount_[i] << "." << endl;
            }
        }

        n = Foam::max(n, nMin);

        n.correctBoundaryConditions();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
