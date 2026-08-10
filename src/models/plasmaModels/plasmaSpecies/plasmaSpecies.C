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

#include "plasmaSpecies.H"
#include "IFstream.H"
#include "HashSet.H"
#include "plasmaConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaSpecies, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::plasmaSpecies::readMechanismSpecies()
{
    fromMechanism_ = true;
    speciesNames_ = speciesFromMechanism(*this, &mechCharge_, &mechMass_);
}


Foam::wordList Foam::plasmaSpecies::speciesFromMechanism
(
    const dictionary& speciesDict,
    HashTable<scalar>* chargeOut,
    HashTable<scalar>* massOut
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

    // Which species to transport.
    //
    //   charged             electron + ions. The minimum for a self-consistent
    //                       discharge: these carry the space charge that drives
    //                       the field, so omitting one is not an approximation,
    //                       it breaks charge conservation.
    //
    //   chargedAndExcited   the above plus excited states and radicals. Needed
    //                       once excited-state chemistry matters -- stepwise
    //                       ionisation, quenching, associative processes. Costs
    //                       one transport equation each.
    //
    //   all                 everything the mechanism names, background gas
    //                       included. Rarely wanted: the background is held
    //                       fixed by construction, so transporting it both
    //                       doubles its cost and lets it drift from the density
    //                       the rate tables were built at.
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
    electronSpeciesID_(-1),
    ionSpeciesIDs_(),
    positiveIonSpeciesIDs_(),
    negativeIonSpeciesIDs_(),
    neutralSpeciesIDs_(),
    chargedSpeciesIDs_(),
    mobileSpeciesIDs_(),
    immobileSpeciesIDs_(),
    constantTemperatureSpeciesIDs_(),
    dynamicTemperatureSpeciesIDs_(),
    followBackgroundTempSpeciesIDs_(),
    solveEnergySpeciesIDs_(),
    fieldTemperatureSpeciesIDs_(),
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

    const dictionary bgEnergyDict = backgroundDict_.subOrEmptyDict("energy");
    bool solveBg = bgEnergyDict.getOrDefault<bool>("solve", false);
    bool bgIsField = solveBg || !bgEnergyDict.found("T");

    constantTemperatureSpeciesIDs_.clear();
    dynamicTemperatureSpeciesIDs_.clear();
    solveEnergySpeciesIDs_.clear();
    fieldTemperatureSpeciesIDs_.clear();

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

        // Energy groups
        const word energy =
            mergedDict.getOrDefault<word>("energyModel", "isothermal");

        // Constant vs Dynamic Temperature
        if (energy == "isothermal")
        {
            constantTemperatureSpeciesIDs_.append(i);
        }
        else if (energy == "backgroundGas")
        {
            followBackgroundTempSpeciesIDs_.append(i);
            if (bgIsField)
            {
                fieldTemperatureSpeciesIDs_.append(i);
                dynamicTemperatureSpeciesIDs_.append(i);
            }
            else
            {
                constantTemperatureSpeciesIDs_.append(i);
            }
        }
        else if (energy == "localField")
        {
            fieldTemperatureSpeciesIDs_.append(i);
            dynamicTemperatureSpeciesIDs_.append(i);
        }
        else if (energy == "solveEnergy")
        {
            solveEnergySpeciesIDs_.append(i);
            fieldTemperatureSpeciesIDs_.append(i);
            dynamicTemperatureSpeciesIDs_.append(i);
        }
    }
}  

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void Foam::plasmaSpecies::updateChargeDensity()
{
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
        n = Foam::max(n, nMin);

        n.correctBoundaryConditions();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
