/*---------------------------------------------------------------------------*\
  File: electromagneticsModel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::electromagneticsModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "materialLibrary.H"
#include "multiRegionPoisson.H"
#include "mappedPatchBase.H"
#include "fixedValueFvPatchFields.H"
#include "electromagneticsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(electromagneticsModel, 0);
defineRunTimeSelectionTable(electromagneticsModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from meshes
electromagneticsModel::electromagneticsModel
(
    const fvMesh& mesh,
    const UPtrList<fvMesh>& dielectricMeshes
)
:
    IOdictionary
    (
        IOobject
        (
            "electromagneticsProperties",
            mesh.time().constant(),
            mesh.time(),
            // READ_IF_PRESENT, not MUST_READ. The file is now OPTIONAL: the
            // model is derived from the region topology (see New()), per-region
            // permittivity lives in constant/<region>/electricalProperties
            // following OpenFOAM's own per-region convention, and the global
            // Poisson numerics live in system/plasmaSimulationControls. What
            // remains here is a legacy fallback, read with a notice.
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),
    mesh_(mesh),
    dielectricMeshes_(dielectricMeshes),
    ePotential_
    (
        IOobject
        (
            "ePotential",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -fvc::grad(ePotential_)
    ),
    Emag_
    (
        IOobject
        (
            "Emag",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(E_)
    ),
    phiE_
    (
        IOobject
        (
            "phiE",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -fvc::snGrad(ePotential_) * mesh.magSf()
    ),
    reducedE_
    (
        IOobject
        (
            "reducedE",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, 4, -3, 0, 0, -1, 0), 0.0)
    ),
    chargeDensity_
    (
        IOobject
        (
            "chargeDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, -3, 1, 0, 0, 1, 0), 0.0)
    ),
    surfCharge_
    (
        IOobject
        (
            "surfCharge",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, -2, 1, 0, 0, 1, 0), 0.0)
    ),
    epsilon_
    (
        "epsilon",
        dimensionSet(-1, -3, 4, 0, 0, 2, 0),
        0.0
    ),
    epsilonR_(1.0)
{}

// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

//- Largest face non-orthogonality angle in the mesh, in degrees: the angle
//  between the face area vector and the owner->neighbour centre vector.
//
//  THE REDUCTION IS UNCONDITIONAL. A rank holding no internal faces must still
//  take part or the collective deadlocks -- rule 31, and the same defect as the
//  gAverage behind a local guard in plasmaExternalCircuit.
static Foam::scalar maxNonOrthogonalityDeg(const Foam::fvMesh& mesh)
{
    const Foam::vectorField& Sf = mesh.Sf().primitiveField();
    const Foam::scalarField& magSf = mesh.magSf().primitiveField();
    const Foam::vectorField& CC = mesh.C().primitiveField();
    const Foam::labelUList& own = mesh.owner();
    const Foam::labelUList& nei = mesh.neighbour();

    Foam::scalar minCos = 1.0;

    forAll(own, facei)
    {
        const Foam::vector d(CC[nei[facei]] - CC[own[facei]]);
        const Foam::scalar denom = magSf[facei]*Foam::mag(d);

        if (denom > Foam::VSMALL)
        {
            minCos = Foam::min(minCos, (Sf[facei] & d)/denom);
        }
    }

    Foam::reduce(minCos, Foam::minOp<Foam::scalar>());

    const Foam::scalar ang =
        Foam::acos
        (
            Foam::min(Foam::max(minCos, Foam::scalar(-1)), Foam::scalar(1))
        );

    return ang*180.0/Foam::constant::mathematical::pi;
}


electromagneticsModel::poissonNumerics
electromagneticsModel::readPoissonNumerics
(
    const fvMesh& mesh,
    const dictionary& legacyCoeffs
)
{
    poissonNumerics n;      // members carry the defaults

    // TOP-LEVEL ONLY -- read with `mesh.time()` as the database, deliberately
    // NOT region-scoped. plasmaTransport has to try system/<region>/ first
    // because outerCoupling can legitimately differ per region; these cannot.
    // One model owns the whole coupled Poisson solve and reads them once, and
    // forming E by a different scheme on either side of an interface would make
    // the coupled flux inconsistent.
    IOdictionary controls
    (
        IOobject
        (
            "plasmaSimulationControls",
            mesh.time().system(),
            mesh.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    // By value: subOrEmptyDict returns a copy, not a reference.
    const dictionary p(controls.subOrEmptyDict("poisson"));

    // Two-level lookup: the `poisson` block wins, the legacy coeffs fill gaps,
    // the struct's own value is the default. `used*` records whether anything
    // actually came from the legacy file, so the notice is only printed when it
    // is doing work.
    bool usedLegacy = false;

    auto pick = [&](const word& newKey, const word& oldKey, auto fallback)
    {
        using T = decltype(fallback);
        if (p.found(newKey))
        {
            return p.get<T>(newKey);
        }
        if (legacyCoeffs.found(oldKey))
        {
            usedLegacy = true;
            return legacyCoeffs.get<T>(oldKey);
        }
        return fallback;
    };

    n.scheme  = pick("scheme",  "PoissonScheme", n.scheme);
    n.EScheme = pick("EScheme", "EScheme",       n.EScheme);
    // NON-ORTHOGONAL CORRECTORS: DERIVED FROM THE MESH unless the case states
    // a number.
    //
    // The default was 0, which on a skewed mesh SILENTLY DROPS the
    // non-orthogonal correction: the Poisson solution is then wrong and
    // nothing says so. Same class of defect as the one found in the
    // Scharfetter-Gummel operator on 2026-09-08, where orthogonal deltaCoeffs
    // left the scheme on a fixed error floor with ZERO convergence under
    // refinement at 11.3 deg non-orthogonality.
    //
    // A corrector count is a property of the MESH, not of the physics, so the
    // user should not have to supply it. Thresholds are the usual practice; an
    // explicit entry in the case still wins.
    {
        const scalar maxNonOrth = maxNonOrthogonalityDeg(mesh);

        const label autoCorr =
            (maxNonOrth <  5.0) ? 0
          : (maxNonOrth < 35.0) ? 2
          : (maxNonOrth < 60.0) ? 3
          :                       4;

        const label asked =
            pick("nNonOrthogonalCorrectors", "nNonOrthogonalCorrectors",
                 autoCorr);

        // THE MESH-DERIVED VALUE IS A FLOOR, not merely a default.
        //
        // A default alone would never fire: every generated case writes the
        // key, so a stale 0 in a case would keep silently dropping the
        // correction on a skewed mesh. Too FEW correctors gives a wrong
        // answer with no symptom; too many only costs time. So the floor is
        // enforced and the case can still ask for MORE.
        n.nNonOrthogonalCorrectors = Foam::max(asked, autoCorr);

        Info<< "  Poisson: mesh max non-orthogonality " << maxNonOrth
            << " deg -> nNonOrthogonalCorrectors "
            << n.nNonOrthogonalCorrectors
            << (n.nNonOrthogonalCorrectors == asked
                  ? (p.found("nNonOrthogonalCorrectors")
                        ? "  (SET BY CASE)" : "  (DERIVED from the mesh)")
                  : "  (RAISED from the case value)")
            << endl;

        if (n.nNonOrthogonalCorrectors > asked)
        {
            WarningInFunction
                << "the case asks for " << asked << " non-orthogonal"
                   " corrector(s) but this mesh has " << maxNonOrth
                << " deg of non-orthogonality," << nl
                << "    which needs at least " << autoCorr
                << ". Using " << autoCorr << "." << nl
                << "    Too few correctors does not fail -- it silently"
                   " returns a wrong potential, which is why this is"
                << nl
                << "    raised rather than obeyed. Refine or improve the mesh"
                   " to remove the cost." << endl;
        }
    }

    // The segregated Picard loop keeps its own sub-dictionary in both places.
    {
        // By value: subOrEmptyDict returns a copy, not a reference.
        const dictionary newRC(p.subOrEmptyDict("nonCoupledResidualControl"));
        const dictionary oldRC
        (
            legacyCoeffs.subOrEmptyDict("nonCoupledResidualControl")
        );

        if (newRC.found("maxIter"))
        {
            n.maxNonCoupledIterations = newRC.get<label>("maxIter");
        }
        else if (oldRC.found("maxIter"))
        {
            usedLegacy = true;
            n.maxNonCoupledIterations = oldRC.get<label>("maxIter");
        }

        if (newRC.found("tolerance"))
        {
            n.nonCoupledTolerance = newRC.get<scalar>("tolerance");
        }
        else if (oldRC.found("tolerance"))
        {
            usedLegacy = true;
            n.nonCoupledTolerance = oldRC.get<scalar>("tolerance");
        }
    }

    if (n.EScheme != "grad" && n.EScheme != "reconstruct")
    {
        FatalErrorInFunction
            << "Unknown EScheme `" << n.EScheme << "`." << nl
            << "    Valid options are: (grad | reconstruct)" << nl
            << "    Set it in system/plasmaSimulationControls under `poisson`."
            << nl << exit(FatalError);
    }

    if (n.scheme != "explicit" && n.scheme != "semiImplicit")
    {
        FatalErrorInFunction
            << "Unknown Poisson `scheme` `" << n.scheme << "`." << nl
            << "    Valid options are: (explicit | semiImplicit)" << nl
            << "    Set it in system/plasmaSimulationControls under `poisson`."
            << nl << exit(FatalError);
    }

    if (usedLegacy)
    {
        Info<< "poisson: some numerics were read from the DEPRECATED"
            << " `<model>Coeffs` in constant/electromagneticsProperties." << nl
            << "    They have moved to system/plasmaSimulationControls:" << nl
            << "        poisson" << nl
            << "        {" << nl
            << "            scheme                    " << n.scheme << ";" << nl
            << "            EScheme                   " << n.EScheme << ";" << nl
            << "            nNonOrthogonalCorrectors  "
            << n.nNonOrthogonalCorrectors << ";" << nl
            << "        }" << nl
            << "    Note `PoissonScheme` is now `scheme`: it is a SCHEME, not a"
            << " linear solver -- the linear solver lives in fvSolution." << nl
            << endl;
    }

    return n;
}


const wordList& electromagneticsModel::knownRegionKinds()
{
    // `farField` is deliberately a KIND rather than a convention.
    //
    // It was expressible before as `dielectric` + epsilonR 1.0 + three correct
    // boundary-condition choices, and that works -- verified exactly: at
    // epsilonR 1.0 the coupled interface reduces to plain continuity of V and
    // dV/dn, and the two-region analytic bed returns the single-medium answer
    // to a relative error of 0.00e+00. But all three BC choices were the
    // author's to remember every time, and getting one wrong is silent.
    //
    // Naming the kind lets the solver default the permittivity and CHECK the
    // boundary conditions instead of hoping.
    static const wordList kinds({"gas", "dielectric", "farField"});
    return kinds;
}


word electromagneticsModel::regionKind(const fvMesh& regionMesh)
{
    IOdictionary rpDict
    (
        IOobject
        (
            "regionProperties",
            regionMesh.time().constant(),
            regionMesh.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    // No regionProperties at all -> a single-region case, which is a gas.
    if (!rpDict.found("regions"))
    {
        return "gas";
    }

    HashTable<wordList> regions;
    rpDict.readEntry("regions", regions);

    forAllConstIters(regions, iter)
    {
        if (iter.val().found(regionMesh.name()))
        {
            return iter.key();
        }
    }

    // Named nowhere. The gas mesh of a case whose regionProperties lists only
    // dielectrics lands here, as does the default region.
    return "gas";
}


scalar electromagneticsModel::readEpsilonR
(
    const fvMesh& regionMesh,
    const dictionary& legacyCoeffs
)
{
    const word kind(regionKind(regionMesh));
    const bool mayDefault = (kind == "gas" || kind == "farField");

    // constant/<region>/electricalProperties -- registered on the REGION mesh,
    // which is what puts it in the per-region directory.
    IOdictionary props
    (
        IOobject
        (
            "electricalProperties",
            regionMesh.time().constant(),
            regionMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    if (props.found("epsilonR"))
    {
        const scalar e = props.get<scalar>("epsilonR");

        // epsilonR == 1 IS ALLOWED AND IS SOMETIMES THE POINT.
        //
        // A "dielectric" region with epsilonR 1.0 is a FICTITIOUS region: it
        // solves Laplace in vacuum/air and nothing else. That is the standard
        // extended-electrostatic-domain trick -- surround a small plasma region
        // with a large Poisson-only region so the far field is resolved without
        // paying for species transport, chemistry or the electron energy
        // equation there, none of which a dielectric region solves.
        //
        // So this is NOT rejected. What is rejected is a value that cannot be a
        // permittivity at all: at epsilonR <= 0 the operator
        // fvm::laplacian(epsilon, ePotential) loses positive-definiteness and
        // the solve returns garbage rather than failing.
        if (e <= 0)
        {
            FatalErrorInFunction
                << "epsilonR = " << e << " in " << props.objectPath() << nl
                << "    Relative permittivity must be POSITIVE." << nl
                << "    At epsilonR <= 0 the Poisson operator is no longer"
                << " positive-definite and the linear solve returns" << nl
                << "    garbage instead of failing." << nl
                << nl
                << "    Note epsilonR = 1 IS valid and is deliberately allowed:"
                << " it makes the region a vacuum/air Laplace" << nl
                << "    domain, which is how you extend the electrostatic"
                << " domain around a small plasma region without" << nl
                << "    solving species there." << nl
                << exit(FatalError);
        }

        return e;
    }

    // NAMED MATERIAL (G2): one line in the case, cited numbers from the
    // library. Reached only when no explicit `epsilonR` was given, so an
    // explicit value always WINS -- see materialLibrary.H for the precedence.
    if (props.found("material"))
    {
        const word mat(props.get<word>("material"));

        const scalar e = materialLibrary::get
        (
            mat,
            "epsilonR",
            "constant/" + regionMesh.dbDir() + "/electricalProperties"
                " (region `" + regionMesh.name() + "`)"
        );

        Info<< "electricalProperties: region `" << regionMesh.name()
            << "` epsilonR = " << e << "  (from MATERIAL `" << mat
            << "`, " << materialLibrary::path() << ")" << endl;

        return e;
    }

    // The path to quote back to the user. dbDir() is EMPTY for the default
    // region, where the file is constant/electricalProperties with no region
    // level at all -- quoting `constant/region0/...` would send a single-region
    // user to a directory that must not exist.
    const fileName where
    (
        fileName("constant")/regionMesh.dbDir()/fileName("electricalProperties")
    );

    // Legacy fallback: the gas read its value from the coeffs dictionary
    // directly; a dielectric from a sub-dictionary named after itself.
    const dictionary legacy
    (
        (kind == "gas") ? legacyCoeffs
                        : legacyCoeffs.subOrEmptyDict(regionMesh.name())
    );

    if (legacy.found("dielectricConstant"))
    {
        const scalar e = legacy.get<scalar>("dielectricConstant");

        Info<< "electricalProperties: region `" << regionMesh.name()
            << "` has no " << where << ", so epsilonR = " << e
            << " was taken from the DEPRECATED `dielectricConstant` in"
            << " constant/electromagneticsProperties." << nl
            << "    Move it, one small file per region, as OpenFOAM does for"
            << " thermophysicalProperties:" << nl
            << "        " << where << ":   epsilonR " << e << ";" << nl
            << endl;

        return e;
    }

    if (mayDefault)
    {
        // 1.0, and for two different reasons.
        //
        //   gas       -- a gas at atmospheric density is a vacuum to within a
        //                few parts in 10^4, so this is physics, not a guess.
        //   farField  -- 1.0 IS the definition of the kind: a fictitious
        //                air/vacuum region that solves Laplace and nothing
        //                else. Requiring the user to state it would be asking
        //                them to repeat the name they already chose.
        if (kind == "farField")
        {
            Info<< "electricalProperties: farField region `"
                << regionMesh.name() << "` uses epsilonR = 1 (air/vacuum)."
                << endl;
        }
        return 1.0;
    }

    // Distinguish "no file" from "file there, value not set". The second is the
    // normal state of a freshly generated template from
    // tools/plasmaSetupRegions.sh, which writes the file with `epsilonR`
    // COMMENTED OUT on purpose -- telling that user to "create" a file they are
    // looking at would be actively misleading.
    const bool haveFile = props.headerOk();

    FatalErrorInFunction
        << "Dielectric region `" << regionMesh.name()
        << "` has no relative permittivity." << nl
        << nl;

    // NOT a ternary: `where` is a fileName and the two branches would build
    // std::string, which has no Ostream inserter here.
    if (haveFile)
    {
        FatalError
            << "    " << where << " EXISTS but sets no `epsilonR`." << nl
            << "    If it was generated by tools/plasmaSetupRegions.sh the"
            << " entry is commented out and waiting for you:" << nl
            << "        epsilonR   4.6;      // uncomment, and use YOUR"
            << " material's value" << nl;
    }
    else
    {
        FatalError
            << "    Create " << where << " containing:" << nl
            << "        epsilonR   4.6;      // your material's value" << nl;
    }

    FatalError
        << nl
        << "    IT IS NOT DEFAULTED ON PURPOSE. A DEFAULT of 1 would be"
        << " vacuum, so a real barrier would go electrically" << nl
        << "    invisible while the run converged and looked fine -- and in a"
        << " DBD the barrier is the whole point, since surface" << nl
        << "    charge on it shields the gap and makes the discharge"
        << " self-limiting. Nothing in the log would say it had gone." << nl
        << "    A guessed value is no better: a run that completes on a"
        << " fabricated permittivity is the kind of result nobody" << nl
        << "    re-checks." << nl
        << nl
        << "    epsilonR 1.0 IS VALID if you MEAN it -- a `dielectric` region"
        << " at 1.0 is a fictitious vacuum/air region that" << nl
        << "    solves only Laplace, which is how you extend the electrostatic"
        << " domain around a small plasma region. State it" << nl
        << "    explicitly and it is accepted." << nl
        << exit(FatalError);

    return 0;
}


void electromagneticsModel::validateRegionKinds(const wordList& kinds)
{
    // constant/regionProperties has the form
    //     regions ( <kind> ( <meshRegionName> ... ) ... );
    // where the KIND is a keyword and the parenthesised names are free. The
    // kinds are looked up by literal string, so a misspelt or invented kind
    // matches nothing and its regions are simply never created.
    //
    // THE FAILURE IS SILENT AND THE ANSWER IS WRONG. Measured 2026-09-02 on the
    // needle-DBD bed: writing
    //     regions ( gas (gas)  plexiglass (plexiglass) );
    // instead of `dielectric (plexiglass)` leaves the dielectric list EMPTY, so
    // New() below derives `singleRegionPoisson`, the dielectric mesh is never
    // built, and the potential is solved in the gas alone -- a converged,
    // plausible, wrong DBD. `foamListRegions` still prints both names, because
    // it flattens every list regardless of kind, so that is no check either.
    //
    // singleRegionPoisson's own guard cannot catch it: that fires when
    // dielectric meshes ARE present, and here there are none.
    //
    // An unknown kind is always a mistake -- there is no use for a region the
    // solver will not create -- so this is fatal rather than a warning.
    const wordList& knownKinds = knownRegionKinds();

    wordList unknown;
    for (const word& kind : kinds)
    {
        if (!knownKinds.found(kind))
        {
            unknown.append(kind);
        }
    }

    if (!unknown.empty())
    {
        FatalErrorInFunction
            << "constant/regionProperties declares region kind(s) this solver"
            << " does not know: " << unknown << nl << nl
            << "    Valid kinds are exactly: " << knownKinds << nl << nl
            << "    The format is" << nl
            << "        regions ( <kind> ( <meshRegionName> ... ) ... );" << nl
            << "    where the KIND is one of the keywords above and the names"
            << " in parentheses are MESH region names." << nl
            << "    Those names are FREE, so" << nl
            << "        dielectric (plexiglass)     is correct" << nl
            << "        plexiglass (plexiglass)     is not" << nl << nl
            << "    WHAT EACH KIND MEANS:" << nl
            << "        gas         the plasma -- species, chemistry, electron"
            << " energy, Poisson" << nl
            << "        dielectric  a solid barrier -- Poisson only,"
            << " `epsilonR` REQUIRED in constant/<region>/" << nl
            << "                    electricalProperties" << nl
            << "        farField    a FICTITIOUS air region -- Poisson only,"
            << " epsilonR defaults to 1.0. Use it to extend" << nl
            << "                    the ELECTROSTATIC domain without extending"
            << " the plasma one." << nl << nl
            << "    This is fatal because it would otherwise be SILENT: an"
            << " unrecognised kind matches nothing, its regions are" << nl
            << "    never created, and a case with a dielectric would run as a"
            << " single-region gas case and give a plausible" << nl
            << "    wrong answer." << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<electromagneticsModel> electromagneticsModel::New
(
    const fvMesh& mesh,
    const UPtrList<fvMesh>& dielectricMeshes
)
{
    // THE MODEL IS DERIVED FROM THE MESH TOPOLOGY, not stated by the case.
    //
    // `singleRegionPoisson` and `multiRegionPoisson` are not two physical
    // models: they solve the same Poisson equation, and which one applies is
    // decided entirely by whether the case has any dielectric region. The
    // solver already knows -- `dielectricMeshes` is right there in the
    // signature -- so asking the user to state it is asking them to maintain a
    // derived quantity by hand.
    //
    // It was also a real footgun in both directions: `multiRegionPoisson` with
    // no dielectric is fatal ("Switch to 'singleRegionPoisson'"), and
    // `singleRegionPoisson` with dielectrics present silently ignores them --
    // a case that looks like a DBD and has no barrier in the field solution.
    //
    // CONTRAST `electronEnergyModel`, which is deliberately REQUIRED with no
    // default: LFA and LMEA are different physics giving different answers, so
    // neither can be assumed. Nothing analogous is being chosen here.
    const bool haveDielectrics = !dielectricMeshes.empty();

    const word modelName
    (
        haveDielectrics ? "multiRegionPoisson" : "singleRegionPoisson"
    );

    // The file is OPTIONAL now, and holds nothing this function needs. Read it
    // only to catch a case still stating the model, so an old case is told
    // where the setting went rather than silently ignored.
    {
        IOdictionary legacyDict
        (
            IOobject
            (
                "electromagneticsProperties",
                mesh.time().constant(),
                mesh.time(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        if (legacyDict.found("electromagneticsModel"))
        {
            const word given(legacyDict.get<word>("electromagneticsModel"));

            if (given != modelName)
            {
                FatalIOErrorInFunction(legacyDict)
                    << "electromagneticsModel is stated as `" << given
                    << "` but this case has "
                    << dielectricMeshes.size() << " dielectric region(s), which"
                    << " requires `" << modelName << "`." << nl
                    << "    The model is DERIVED from constant/regionProperties"
                    << " -- it is not a choice, because both models solve the"
                    << " same equation and only the mesh topology decides which"
                    << " applies." << nl
                    << "    Delete the entry." << nl
                    << exit(FatalIOError);
            }

            Info<< "electromagneticsModel: `" << given << "` in"
                << " constant/electromagneticsProperties is no longer read --"
                << " the model is DERIVED from constant/regionProperties." << nl
                << "    It agrees with the topology here, so the run continues."
                << " Delete the entry; per-region permittivity now lives in"
                << " constant/<region>/electricalProperties and the global"
                << " Poisson numerics in system/plasmaSimulationControls."
                << endl;
        }
    }

    Info<< "electromagneticsModel: " << modelName << " (DERIVED -- "
        << dielectricMeshes.size() << " dielectric region(s) in"
        << " constant/regionProperties)" << endl;

    // Look up the constructor in the table
    auto* ctorPtr = dictionaryConstructorTable(modelName);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "No electromagneticsModel registered as '" << modelName
            << "'. This is a build problem, not a case problem." << nl
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<electromagneticsModel>
    (
        ctorPtr(mesh, dielectricMeshes)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam



void Foam::electromagneticsModel::classifyElectrodePatches
(
    wordList& driven,
    wordList& grounded,
    wordList& floating
) const
{
    DynamicList<word> dr, gr, fl;

    auto classify = [&](const volScalarField& ePot)
    {
        const volScalarField::Boundary& bf = ePot.boundaryField();

        forAll(bf, patchi)
        {
            const fvPatchScalarField& pf = bf[patchi];

            // Not an imposed-value condition: not a conductor.
            if (!isA<fixedValueFvPatchScalarField>(pf)) continue;

            // A region interface is never an electrode, whatever sits on it.
            if (isA<mappedPatchBase>(pf.patch().patch())) continue;

            const word& pname = pf.patch().name();

            // FLOATING first, because it IS a fixedValue by inheritance and
            // would otherwise be classified as a driven or grounded electrode
            // -- its stored value is simply the last equipotential it reached.
            //
            // Compared BY TYPE NAME rather than with isA<>, deliberately: the
            // condition lives in libplasmaBcs, which depends on THIS library,
            // so including its header here would be circular. The name is part
            // of the condition's public interface (its TypeName), so this is a
            // stable contract rather than a guess.
            if (pf.type() == "floatingElectrodePotential")
            {
                fl.append(pname);
                continue;
            }

            const bool timeVarying =
                (pf.type() != fixedValueFvPatchScalarField::typeName);

            const scalar peak = gMax(mag(pf));

            if (!timeVarying && peak < SMALL)
            {
                gr.append(pname);
            }
            else
            {
                dr.append(pname);
            }
        }
    };

    classify(ePotential_);

    if (isA<multiRegionPoisson>(*this))
    {
        const multiRegionPoisson& mrp =
            refCast<const multiRegionPoisson>(*this);

        for (label i = 0; i < mrp.nDielectrics(); ++i)
        {
            classify(mrp.dielectric(i).ePotential());
        }
    }

    driven.transfer(dr);
    grounded.transfer(gr);
    floating.transfer(fl);
}




void Foam::electromagneticsModel::correctFloatingElectrode
(
    const volScalarField* effEpsGas
)
{
    if (!floatingChecked_)
    {
        floatingChecked_ = true;

        if (floatingElectrode::present(*this))
        {
            floating_.reset(new floatingElectrode(mesh_, *this));
        }
    }

    if (floating_)
    {
        floating_->correct(*this, effEpsGas);
    }
}


// ************************************************************************* //
