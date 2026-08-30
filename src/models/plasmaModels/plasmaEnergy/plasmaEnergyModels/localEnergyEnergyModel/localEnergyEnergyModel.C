/*---------------------------------------------------------------------------*\
  File: localEnergyEnergyModel.C
  Part of: SoPLASMA
  Copyright (C) 2026
  License: GNU General Public License v3 or later
\*---------------------------------------------------------------------------*/

#include "zeroGradientFvPatchFields.H"
#include "localEnergyEnergyModel.H"
#include "plasmaSpecies.H"
#include "plasmaConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvc.H"
#include "IFstream.H"

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
        dimensionedScalar("zero", dimless/dimVolume, 0.0),

        // SEED WITH zeroGradient, not the default `calculated`.
        //
        // READ_IF_PRESENT means the field is created here when the case does
        // not supply 0/nEps_e -- which is the normal way to enable LMEA. The
        // default patch type for this constructor is `calculated`, and
        // fvMatrix REFUSES to solve a field with calculated patches:
        //
        //     cannot be called for a calculatedFvPatchField
        //     on patch <name> of field nEps_e
        //
        // so the solver died at the first energy solve and the user had to
        // hand-write 0/nEps_e to get past it. Enabling a model with one
        // keyword should not require also knowing which field file to author.
        //
        // Constraint patches (wedge, empty, processor, cyclic) still receive
        // their own types -- GeometricField honours constraintType over a
        // requested default -- so this is safe on the wedge meshes these
        // cases use. A case that wants something other than zeroGradient
        // (e.g. energyDDWallFluxMixed at the electrodes) still supplies its
        // own 0/nEps_e and that is read instead.
        zeroGradientFvPatchScalarField::typeName
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
            // BARE name, not per-species, and that is required rather than
            // stylistic: fromMechanism builds its path as
            //     quantity + "_vs_" + lookupVariable
            // so `lookupVariable meanE` resolves to muN_vs_meanE, the table
            // the Boltzmann solver actually writes. `meanE_e` would look for
            // muN_vs_meanE_e and fail. It also matches `reducedE`, which is
            // likewise global -- mean energy is an electron property and there
            // is one electron species.
            "meanE",
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
        dimensionedScalar("zero", dimensionSet(-1, 0, 2, 0, 0, 1, 0), 0.0)
    ),
    muEpsf_
    (
        IOobject("muEps_lmea", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, 0, 2, 0, 0, 1, 0), 0.0)
    ),
    DEf_
    (
        IOobject("DE_lmea", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
    ),
    DEpsEff_
    (
        IOobject("DEps", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
    ),
    muEpsEff_
    (
        IOobject("muEps", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, 0, 2, 0, 0, 1, 0), 0.0)
    ),
    PlossN_
    (
        IOobject("PlossN_lmea", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0, 0, 0), 0.0)
    ),
    // WRITTEN, because a claim about it cannot otherwise be checked.
    //
    // This is the outer loop's self-damping coefficient, -dS/d(n_eps). The
    // Option 4 Jacobian notes record that it "was MEASURED going to ~0 in cold
    // cells, because P_loss sits below every inelastic threshold there" -- i.e.
    // exactly where the Picard coupling would be left undamped. That is a
    // testable statement about a field, and it was NO_WRITE, so it could not
    // be tested on any run that had already happened.
    dSdEps_
    (
        IOobject("dSdEps_lmea", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    )
{
    // `dict` IS energyModelCoeffs -- plasmaEnergy passes that subdict
    // directly, as localFieldEnergyModel's `temperature` lookup shows. A
    // nested `localEnergyCoeffs` would break the convention every other
    // energy model follows.
    const dictionary& c = dict;

    // See nEfloor_: eps_bar = n_eps/n_e is 0/0 in the quiescent far field.
    nEfloor_  = c.getOrDefault<scalar>("electronDensityFloor", 1.0);
    // FLOOR FROM COLD ELECTRONS, not an arbitrary small number.
    //
    // Electrons cannot be colder than the gas they collide with: elastic
    // exchange drives them to the gas temperature, and the tabulated elastic
    // power loss changes sign there. So the physical floor on the mean energy
    // is the thermal value
    //
    //     eps_min = (3/2) k_B T_gas    [eV]
    //
    // which is 0.0388 eV at 300 K. This is RIGID physics -- unlike, say, the
    // secondary-electron emission energy, which has no such derivation and is
    // therefore a user default. Overridable all the same, for a case that
    // wants to probe below it.
    {
        const dictionary& bg = species.backgroundDict();
        const scalar Tgas =
            bg.subOrEmptyDict("energy").getOrDefault<scalar>("T", 300.0);

        const scalar kB_eV = constant::plasma::kappaBoltzmann.value()
                           / constant::plasma::eCharge.value();

        meanEmin_ = c.getOrDefault<scalar>("meanEnergyMin", 1.5*kB_eV*Tgas);
    }
    meanEmax_ = c.getOrDefault<scalar>("meanEnergyMax", 100.0);
    sourceSensitivity_ =
        c.getOrDefault<Switch>("sourceSensitivity", true);
    transportEnergy_ = c.getOrDefault<Switch>("transportEnergy", true);
    reportInterval_ = c.getOrDefault<label>("reportInterval", 1);
    lastReport_ = -1;
    meanEprev_ = 0.0;
    clampedSteps_ = 0;

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

    // ENERGY transport coefficients.
    //
    // PREFER the tabulated muEpsN / DEpsN, which are the real energy-flux
    // moments of the EEDF. Fall back to (5/3) x the particle values only when
    // a case has no such tables -- the Maxwellian relation, and MEASURED on
    // this air set to be wrong by up to 35% in mobility and 60% in diffusion:
    //
    //     E/N [Td]    mu_eps/mu    D_eps/D_T     (Maxwellian: 1.6667 both)
    //         13        1.557        1.062
    //        106        1.500        1.505
    //        428        1.492        1.305
    //       1723        1.545        1.467
    //
    // Which fallback applied is REPORTED at start-up, never silent: a reader
    // must be able to tell whether the coefficients came from the Boltzmann
    // solution or from an assumption about it.
    energyCoeffsTabulated_ =
        c.found("energyMobility") && c.found("energyDiffusivity");

    const dictionary& muDict =
        energyCoeffsTabulated_ ? c.subDict("energyMobility")
                               : c.subDict("mobility");

    const dictionary& DDict =
        energyCoeffsTabulated_ ? c.subDict("energyDiffusivity")
                               : c.subDict("diffusivity");

    // PARTICLE mobility, always from `mobility`: the Joule term acts on the
    // electrons, not on their energy.
    muE_ = plasmaPropertyEvaluator::New
    (
        c.subDict("mobility").get<word>("type"), c.subDict("mobility"),
        mesh, "ElectronMobility", dimensionSet(-1, 0, 2, 0, 0, 1, 0)
    );

    // ENERGY mobility, for the energy flux only.
    muEps_ = plasmaPropertyEvaluator::New
    (
        muDict.get<word>("type"), muDict,
        mesh, "ElectronEnergyMobility", dimensionSet(-1, 0, 2, 0, 0, 1, 0)
    );

    DE_ = plasmaPropertyEvaluator::New
    (
        DDict.get<word>("type"), DDict,
        mesh, "ElectronEnergyDiffusivity", dimensionSet(0, 2, -1, 0, 0, 0, 0)
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

    // SEED. Without one the run starts at meanEmin_, and the transport
    // coefficients looked up at ~1e-3 eV are nowhere near equilibrium -- a
    // pathological first step in a quantity that relaxes in ~5 ps. Seeding a
    // plausible mean energy costs nothing and is overwritten within a few
    // steps by the physics.
    // Seed on a FRESH START even though the field was read from 0/.
    //
    // The 0/ file exists to carry the BOUNDARY CONDITIONS -- a transported
    // field needs them -- but its internalField is a placeholder. Letting it
    // win means starting from a uniform guess inconsistent with the field
    // that already exists at t=0, which is what produced the first-step
    // excursion to 1499 eV. A RESTART is different: there the stored field is
    // a real solution and must not be touched. The two are told apart by the
    // time index -- at startTime nothing has been solved yet.
    const bool freshStart =
        mesh.time().timeIndex() == mesh.time().startTimeIndex();

    if
    (
        gMax(nEps_.primitiveField()) <= 0
     || (freshStart && c.found("initialMeanEnergy"))
    )
    {
        // SEED FROM THE LFA EQUILIBRIUM, not from a uniform guess.
        //
        // At t=0 the field is ALREADY established -- the Poisson solve has the
        // applied voltage and the seed space charge -- so a uniform cold start
        // is inconsistent with it, and the first energy solve has to travel
        // from cold to field-consistent in one iteration. MEASURED with a
        // uniform 1 eV seed: `Initial residual 0.9997` on the first nEps
        // solve, and an intermediate iterate reaching 1499 eV, which pins the
        // mean energy at the clamp and feeds EXTRAPOLATED rate coefficients --
        // exponential in energy -- to that corrector. It recovered here; that
        // is luck, not robustness.
        //
        // meanEnergy_vs_reducedE is exactly the state LMEA relaxes to without
        // transport, so seeding it starts the run in equilibrium and leaves
        // the energy equation to supply only the NON-local correction, which
        // is what it is for.
        // The LFA seed is DEFERRED to the first correct(): at construction
        // the Poisson solve has not run, reducedE is zero everywhere, and the
        // lookup returns the BOTTOM of the table. Measured exactly that:
        // "meanE in [0.0834, 0.0834] eV", the table's first entry.
        if (c.found("initialMeanEnergy"))
        {
            seedFromLFA_ = true;

            seedEval_ = plasmaPropertyEvaluator::New
            (
                c.subDict("initialMeanEnergy").get<word>("type"),
                c.subDict("initialMeanEnergy"), mesh,
                "InitialMeanEnergy", dimless
            );

            Info<< "    n_eps will be seeded from the LFA equilibrium after"
                << " the first field solve" << endl;
        }
        else
        {
            const scalar e0 = c.getOrDefault<scalar>("meanEnergy0", 1.0);

            nEps_ == e0*species_.numberDensity(specieIndex_);

            Info<< "    seeded n_eps at a UNIFORM " << e0 << " eV -- supply"
                << " `initialMeanEnergy` (meanEnergy vs reducedE) to start"
                << " from the LFA equilibrium instead and avoid a first-step"
                << " transient" << endl;
        }
    }

    // Top of the tabulated range, read from the mobility table this model
    // actually uses, so the extrapolation warning refers to the real limit
    // rather than a hardcoded guess.
    {
        const dictionary& mdict = c.subDict("mobility");
        const fileName td
        (
            mdict.getOrDefault<fileName>("tableDir", "constant/plasmaTables")
        );
        fileName tf
        (
            mdict.found("file")
          ? mdict.get<fileName>("file")
          : td/(mdict.get<word>("quantity") + "_vs_"
                + mdict.get<word>("lookupVariable"))
        );
        tf.expand();

        IFstream is(tf);
        if (is.good())
        {
            scalar xLast = 0;
            string line;
            while (is.getLine(line))
            {
                const auto o = line.find('(');
                const auto e = line.find(' ', o == string::npos ? 0 : o + 1);
                if (o != string::npos && e != string::npos)
                {
                    const string num = line.substr(o + 1, e - o - 1);
                    if (num.size() && (isdigit(num[0]) || num[0] == '-'))
                    {
                        xLast = max(xLast, readScalar(IStringStream(num)()));
                    }
                }
            }
            tableMaxMeanE_ = xLast;
        }
    }

    Info<< "  energyModel localEnergy (LMEA) for `"
        << species.speciesNames()[specieIndex] << "`:" << nl
        << "    transported n_eps, mean energy published as `"
        << meanE_.name() << "`" << nl
        << "    energy transport: "
        << (energyCoeffsTabulated_
              ? "muEpsN / DEpsN from the EEDF (exact moments)"
              : "5/3 x particle values (MAXWELLIAN approximation -- provide"
                " energyMobility/energyDiffusivity to use the real ones)")
        << nl
        << "    n_e floor " << nEfloor_ << " 1/m^3, mean energy clamped to ["
        << meanEmin_ << ", " << meanEmax_ << "] eV" << endl;

    updateDerived();
}


// One volt, as a dimensioned quantity.
//
// n_eps is carried in eV/m^3 with the eV treated as DIMENSIONLESS, while E
// carries full SI dimensions. So mu ne |E|^2 -- a power density in W/m^3 --
// exceeds the eV/(m^3 s) this equation needs by exactly one factor of volts.
// Dividing by this converts joules to electron-volts and is the only place
// that conversion happens.
static const dimensionedScalar oneVolt
(
    "oneVolt", dimensionSet(1, 2, -3, 0, 0, -1, 0), 1.0
);



// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void localEnergyEnergyModel::updateDerived()
{
    const volScalarField& ne = species_.numberDensity(specieIndex_);

    scalarField& me = meanE_.primitiveFieldRef();
    scalarField& Tf = T_.primitiveFieldRef();

    // WRITABLE: the upper clamp must be written back into the state, not
    // merely into the reported mean energy. See the write-back below.
    scalarField& nEpsI = nEps_.primitiveFieldRef();
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

        const scalar raw = nEpsI[c]/n;

        // CAPTURE the state of the first cell that needs the upper clamp.
        // Counting clamp hits says a cell is out of range; it does not say
        // WHY, and the two causes want opposite fixes: a genuinely large
        // energy density is a physics/step problem, while a collapsed n_e is
        // a FLOOR problem -- eps = n_eps/n_e blowing up because the
        // denominator vanished, not because the numerator grew.
        if (raw > meanEmax_)
        {
            // GLOBAL quantities, not a rank-local cell index. Info<< prints
            // on master only, while the offending cell lives on whichever
            // rank owns it -- so storing a local index means master has
            // nothing to report and the diagnostic prints nothing at all.
            // (Third time this locality trap has bitten in one session.)
            clampRaw_ = max(clampRaw_, raw);
            clampNe_ = min(clampNe_, neI[c]);
            clampCell_ = 1;
        }

        me[c] = min(max(raw, meanEmin_), meanEmax_);
        Tf[c] = me[c]*eVtoK;

        // WRITE THE UPPER CLAMP BACK INTO THE STATE.
        //
        // Clamping only meanE_ leaves nEps_ free to grow without bound, and
        // the loss reaches the matrix as a RATE, L = P_loss/n_eps (see
        // eEqn()). P_loss is looked up at the CLAMPED energy so its magnitude
        // is right, but dividing it by an UNCLAMPED n_eps makes L too small
        // by exactly the overshoot factor. That is a positive feedback --
        // bigger n_eps -> weaker damping -> bigger n_eps -- and a cell that
        // once crosses the clamp can never cool back down.
        //
        // MEASURED, before this fix: raw eps 5.3e7 eV against a 100 eV clamp,
        // so L was 530,000x too weak. The step diagnostic reported J/L =
        // 0.0028, a sink 360x the source, while the matrix received almost no
        // sink at all -- the reported balance and the solved balance had come
        // apart. 67 of 68 steps were clamped and the run died on SIGFPE.
        //
        // The LOWER bound is already enforced on the state further down
        // ("COLD-ELECTRON floor", v[c] = floorV). This restores the symmetry:
        // both bounds now limit n_eps itself, not just what is reported.
        //
        // Non-conservative by construction, exactly like the species density
        // clip -- energy is discarded here. clampRaw_ records how far out the
        // worst cell was, so a run leaning on this says so in the log.
        if (raw > meanEmax_)
        {
            nEpsI[c] = me[c]*n;
        }
    }

    meanE_.correctBoundaryConditions();
    T_.correctBoundaryConditions();

    // NON-FINITE / NEGATIVE CHECK, with the state that produced it.
    //
    // Without this the run dies on SIGFPE with a stack trace that names the
    // linear solver and tells you nothing about WHICH cell went bad or WHY.
    // Reporting the local state turns an opaque crash into a located one --
    // the same reasoning as the minDeltaT diagnostic.
    {
        // ORDER MATTERS, and the first version had it backwards: it
        // detected, THEN clipped, THEN reported the clipped cell as fatal --
        // killing the run over a cell it had just successfully repaired.
        //
        // Negatives are RECOVERABLE and are auto-corrected; only a non-finite
        // value is unrecoverable, because nothing sensible can be substituted
        // for a NaN.

        // 1. Clip negatives. An energy density is n_e times a mean energy,
        //    both non-negative, so a negative value is meaningless -- exactly
        //    the clip the species densities already apply.
        //
        //    MEASURED cause: at the offending cell the source was strongly
        //    POSITIVE (Joule 1.09e25 against loss 9.2e22), so the source did
        //    not do it. The convective term did: energy drifts at 5/3 the
        //    electron speed, so its Courant number is 5/3 of the species one.
        label nClipped = 0;
        {
            scalarField& v = nEps_.primitiveFieldRef();
            const scalarField& neI =
                species_.numberDensity(specieIndex_).primitiveField();
            forAll(v, c)
            {
                // Floor at the COLD-ELECTRON energy density, not at zero.
                // Zero would make eps_bar = 0, colder than the gas, which is
                // unphysical and leaves the next lookup at the bottom edge of
                // every table.
                const scalar floorV = meanEmin_*max(neI[c], scalar(0));

                if (v[c] < floorV)
                {
                    v[c] = floorV;
                    ++nClipped;
                }
            }
        }
        reduce(nClipped, sumOp<label>());

        if (nClipped > 0 && !negativeWarned_)
        {
            negativeWarned_ = true;

            WarningInFunction
                << "electron energy density fell below the COLD-ELECTRON"
                << " floor in " << nClipped << " cell(s) and was clipped."
                << nl
                << "    Floor: eps_min = (3/2) k_B T_gas = " << meanEmin_
                << " eV, the thermal energy of electrons in equilibrium with"
                << " the gas." << nl
                << "    The energy convects at 5/3 the electron drift speed,"
                << " so its Courant number is 5/3 of the species one: a case"
                << " at maxSpeciesConvectiveCo 1.5 runs the ENERGY equation"
                << " at about 2.5." << nl
                << "    Lower maxSpeciesConvectiveCo if this persists."
                << endl;
        }

        // 2. Non-finite is fatal, and reported WITH the state that produced
        //    it. A cell index is not meaningful across ranks, so each rank
        //    reports its own finding and the exit is collective.
        label bad = -1;
        forAll(nEps_.primitiveField(), c)
        {
            if (!std::isfinite(nEps_.primitiveField()[c]))
            {
                bad = c;
                break;
            }
        }

        bool anyBad = (bad >= 0);
        reduce(anyBad, orOp<bool>());

        if (bad >= 0)
        {
            const volScalarField& neD =
                species_.numberDensity(specieIndex_);

            FatalErrorInFunction
                << "electron energy density is NOT FINITE." << nl
                << "    cell        " << mesh_.C()[bad] << nl
                << "    n_eps       " << nEps_.primitiveField()[bad] << nl
                << "    n_e         " << neD[bad] << " 1/m^3" << nl
                << "    meanE       " << meanE_.primitiveField()[bad]
                << " eV" << nl
                << "    |E|         " << mag(E_[bad]) << " V/m  (E/N = "
                << mag(E_[bad])/species_.backgroundDensity().value()/1e-21
                << " Td)" << nl
                << "    mu          " << muEf_[bad] << nl
                << "    PlossN      " << PlossN_[bad] << nl
                << "    dS/d(n_eps) " << dSdEps_[bad] << " 1/s" << nl << nl
                << "    The source or a transport coefficient diverged."
                << " Check the step against the energy relaxation time"
                << " eps/(N PlossN), and the mean energy against the"
                << " tabulated range." << nl
                << abort(FatalError);
        }

        if (anyBad)
        {
            FatalErrorInFunction
                << "electron energy density diverged on another rank; see the"
                << " located report above." << nl
                << exit(FatalError);
        }
    }

    // EXTRAPOLATION WARNING. Extrapolating is deliberate -- clamping to the
    // table would be less robust, since a transient excursion would silently
    // freeze the physics at the table edge -- but the user must be told the
    // coefficients are no longer interpolated data.
    if (!extrapWarned_)
    {
        const scalar hi = gMax(meanE_.primitiveField());

        if (hi > tableMaxMeanE_ && tableMaxMeanE_ > 0)
        {
            extrapWarned_ = true;

            WarningInFunction
                << "mean energy reached " << hi << " eV, beyond the tabulated"
                << " range (" << tableMaxMeanE_ << " eV)." << nl
                << "    Transport and rate coefficients there are"
                << " EXTRAPOLATED, not interpolated: they are no longer"
                << " supported by the Boltzmann solution." << nl
                << "    Extend the EEDF sweep past this energy if the region"
                << " matters, or check whether the excursion is itself a"
                << " symptom." << endl;
        }
    }
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void localEnergyEnergyModel::correct()
{
    // FIRST call only. The field now exists, so the LFA equilibrium can be
    // evaluated at the REAL reducedE -- which is the whole point of deferring
    // it out of the constructor.
    if (seedFromLFA_)
    {
        seedFromLFA_ = false;

        volScalarField meanE0(meanE_);
        seedEval_->correct(meanE0);

        nEps_ == meanE0*species_.numberDensity(specieIndex_);

        Info<< "  LMEA: seeded n_eps from the LFA equilibrium at the SOLVED"
            << " field, meanE in [" << gMin(meanE0.primitiveField()) << ", "
            << gMax(meanE0.primitiveField()) << "] eV" << endl;
    }

    updateDerived();

    // MEAN-ENERGY RANGE, the diagnostic that tells runaway from oscillation.
    muE_->correct(muEf_);
    muEps_->correct(muEpsf_);
    DE_->correct(DEf_);

    // Published for boundary conditions, from the SAME expression the energy
    // equation uses below, so the two cannot drift apart.
    {
        const dimensionedScalar energyFactor
        (
            "energyFactor", dimless, energyCoeffsTabulated_ ? 1.0 : 5.0/3.0
        );
        DEpsEff_  == energyFactor*DEf_;
        muEpsEff_ == energyFactor*muEpsf_;
    }

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

    // The report runs AFTER the coefficients are refreshed.
    //
    // It used to run before, so it printed the NEW mean energy beside the
    // PREVIOUS iteration's mu and PlossN -- inconsistent state that read as
    // 'extrapolated nonsense' when the coefficients were simply stale. A
    // diagnostic that mixes two iterates is worse than none.

    if (reportInterval_ > 0)
    {
        const label ti = mesh_.time().timeIndex();

        const scalar lo = gMin(meanE_.primitiveField());
        const scalar hi = gMax(meanE_.primitiveField());

        // A cell sitting ON the upper clamp is not a diagnostic detail: the
        // transport and rate lookups there are pinned at the top of the
        // tables, so the physics has stopped responding to the solution.
        // Counted per CALL, and correct() runs once per outer iteration --
        // so this counts intermediate iterates, not converged solutions.
        // Those need opposite responses: an intermediate overshoot that later
        // correctors pull back is the outer loop doing its job, while a
        // CONVERGED value beyond the table means the solution genuinely left
        // the tabulated range.
        //
        // MEASURED here: raw eps reached 1499 eV on the FIRST outer iteration
        // of a step (the nEps solve starts from residual ~1, the seed being
        // replaced wholesale by a field-driven solution), while the converged
        // per-step values were 1 -> 3.33 -> 3.24 eV. So the excursions seen so
        // far are transients, NOT divergence -- which is why the distinction
        // is drawn rather than the clamp simply being tightened.
        if (hi >= meanEmax_*(1.0 - 1e-6))
        {
            ++clampedSteps_;

            if (mesh_.time().timeIndex() == lastClampTimeIndex_)
            {
                ++clampedItersThisStep_;
            }
            else
            {
                lastClampTimeIndex_ = mesh_.time().timeIndex();
                clampedItersThisStep_ = 1;
            }
        }

        if (ti != lastReport_ && (ti % reportInterval_) == 0)
        {
            lastReport_ = ti;

            const scalar d = (meanEprev_ > 0) ? (hi - meanEprev_)/meanEprev_
                                              : 0.0;

            // SOURCE BALANCE. The single most diagnostic number this model
            // has, and it is checkable on STEP 1 rather than after a run:
            // in quasi-equilibrium Joule and P_loss must be the SAME ORDER,
            // differing only by the growth term. Orders of magnitude apart
            // means a term is wrong -- which is how the power-loss tables
            // being run through a density-dividing evaluator was found, after
            // watching a mean energy climb for seven steps instead.
            const volScalarField& neD = species_.numberDensity(specieIndex_);

            const scalar jouleAvg = gSum
            (
                mesh_.V()*(muEf_*neD*magSqr(E_)/oneVolt)().primitiveField()
            );
            const scalar plossAvg = gSum
            (
                mesh_.V()
               *(PlossN_*neD*species_.backgroundDensity())().primitiveField()
            );

            // WHERE the maximum sits, and the local field there. Without
            // this a high mean energy is uninterpretable: 20 eV at a sharp
            // electrode corner where E/N is thousands of Td is expected, the
            // same value in the bulk is a bug. Reported rather than guessed.
            label iMax = -1;
            scalar hiLocal = -GREAT;
            forAll(meanE_.primitiveField(), c)
            {
                if (meanE_.primitiveField()[c] > hiLocal)
                {
                    hiLocal = meanE_.primitiveField()[c];
                    iMax = c;
                }
            }

            // iMax/hiLocal are RANK-LOCAL. `Info<<` prints on master only, so
            // reporting them directly prints MASTER'S local maximum -- which
            // in parallel is some ordinary bulk cell -- immediately below the
            // GLOBAL `hi` from gMax(), which lives on a different rank.
            //
            // MEASURED, before this was fixed: the same report read "[0.128,
            // 100] eV" and then "at that cell: meanE 0.6817 eV, J/L 0.66", a
            // healthy cell at local equilibrium. Two unrelated cells printed
            // as one. Every diagnosis drawn from that block was wrong.
            //
            // Reduce the VALUE to find the owning rank (ties broken by the
            // lowest rank), then have ONLY that rank contribute the payload
            // to a sum-reduce, so master prints the real cell. This is
            // outside the per-cell loop, so the collectives cannot deadlock
            // on ranks holding different cell counts.
            const scalar hiGlobal = returnReduce(hiLocal, maxOp<scalar>());

            // Lowest rank holding the global maximum owns it, so a tie
            // contributes exactly once to the sum-reduce below.
            const label ownerProc = returnReduce
            (
                (iMax >= 0 && hiLocal >= hiGlobal)
              ? Pstream::myProcNo() : labelMax,
                minOp<label>()
            );
            const bool iAmOwner = (Pstream::myProcNo() == ownerProc);

            // 0 hiLocal, 1-3 x, 4 E/N, 5 n_e, 6 mu, 7 PlossN, 8 |E|
            scalarList pay(9, Zero);
            if (iAmOwner)
            {
                const scalar NL = species_.backgroundDensity().value();
                pay[0] = hiLocal;
                pay[1] = mesh_.C()[iMax].x();
                pay[2] = mesh_.C()[iMax].y();
                pay[3] = mesh_.C()[iMax].z();
                pay[4] = mag(E_[iMax])/NL/1e-21;
                pay[5] = species_.numberDensity(specieIndex_)[iMax];
                pay[6] = muEf_[iMax];
                pay[7] = PlossN_[iMax];
                pay[8] = mag(E_[iMax]);
            }
            forAll(pay, i) reduce(pay[i], sumOp<scalar>());

            const vector xMax(pay[1], pay[2], pay[3]);
            const scalar ENmax = pay[4];

            Info<< "  LMEA mean energy:         [" << lo << ", " << hi
                << "] eV   max " << (d >= 0 ? "+" : "") << 100*d << "% since"
                << " last report" << nl
                << "    source balance:         Joule " << jouleAvg
                << ", loss " << plossAvg << ", ratio "
                << jouleAvg/max(plossAvg, VSMALL)
                << "  (order 1 expected in quasi-equilibrium)" << nl
                << "    max at:                 " << xMax
                << "  local E/N = " << ENmax << " Td";

            // LOCAL terms at that cell. Domain integrals hide this: the
            // balance can look order-1 overall while one cell is wildly out,
            // and it was the per-cell numbers that were needed all along.
            {
                // From the reduced payload, so this is the GLOBAL maximum's
                // cell whichever rank owns it -- not master's local one.
                const scalar neL = pay[5];
                const scalar muL = pay[6];
                const scalar plL = pay[7];
                const scalar EL  = pay[8];
                const scalar NL  = species_.backgroundDensity().value();

                hiLocal = pay[0];

                const scalar jL = muL*neL*EL*EL;
                const scalar lL = plL*neL*NL;

                Info<< nl
                    << "    at that cell:           meanE " << hiLocal
                    << " eV, n_e " << neL << " 1/m3" << nl
                    << "      mu " << muL << ", PlossN " << plL
                    << ", |E| " << EL << " V/m" << nl
                    << "      Joule " << jL << ", loss " << lL
                    << ", J/L " << jL/max(lL, VSMALL)
                    << "   (must be ~1 at the local equilibrium)";
            }

            if (clampedSteps_ > 0)
            {
                Info<< nl << "    WARNING: max pinned at the " << meanEmax_
                    << " eV clamp on " << clampedSteps_ << " step(s) -- the"
                    << " lookups are at the top of the tables and no longer"
                    << " respond to the solution." << nl
                    << "      n_eps is RESET to eps_max*n_e in those cells, so"
                    << " the loss rate L = P_loss/n_eps stays consistent with"
                    << " the clamped energy. Energy is DISCARDED there: this"
                    << " is a limiter, not physics.";

                // Reduced across ranks, so master reports the real extremes.
                scalar rawG = clampRaw_;
                scalar neG = clampNe_;
                reduce(rawG, maxOp<scalar>());
                reduce(neG, minOp<scalar>());

                // Only when something clamped THIS step. The counter is
                // cumulative, so printing on it alone reported `worst raw
                // eps 0` on steps where nothing happened -- a stale line that
                // reads like a measurement.
                if (rawG > 0)
                {
                // Compare against the floor that ACTUALLY BINDS. `nEfloor_`
                // defaults to 1.0, thirteen decades below the species
                // `minNumberDensity` (1e13), so max(n_e, nEfloor_) NEVER
                // fires and quoting it made cells sitting exactly ON the
                // species floor read as "far from the floor, so healthy" --
                // the exact opposite of the truth. Report the domain minimum
                // n_e and let the comparison be visible rather than asserted.
                scalar neDomMin =
                    gMin(species_.numberDensity(specieIndex_).primitiveField());

                Info<< nl << "      worst raw eps " << rawG
                    << " eV, smallest n_e among clamped cells " << neG
                    << " 1/m3" << nl
                    << "      (domain min n_e " << neDomMin
                    << ", model floor " << nEfloor_ << ";"
                    << " if the clamped n_e sits AT the domain minimum the"
                    << " denominator has collapsed to the species floor --"
                    << " the cell is empty, not hot)";
                }
            }
            Info<< endl;

            meanEprev_ = hi;
            clampRaw_ = 0;
            clampNe_ = GREAT;
            clampCell_ = -1;
        }
    }
    // SOURCE SENSITIVITY, by finite difference on the tabulated coefficients:
    //
    //     D = -dS/deps / n_e = N [ dPlossN/deps - (E/N)^2 dmuN/deps ]
    //
    // Evaluated by re-running the evaluators on a perturbed mean energy. That
    // is why meanE_ is a REGISTERED field: the evaluators read their lookup
    // coordinate from the registry, so perturbing it and re-correcting is the
    // only way to differentiate a table this code does not own.
    dSdEps_ == dimensionedScalar("zero", dSdEps_.dimensions(), 0.0);

    if (sourceSensitivity_)
    {
        const volScalarField meanE0(meanE_);

        // Step proportional to the local energy, floored so a near-zero cell
        // does not divide by a vanishing h.
        const volScalarField h
        (
            max(1e-3*meanE0, dimensionedScalar("hmin", dimless, 1e-4))
        );

        volScalarField muP(muEf_), muM(muEf_);
        volScalarField PlP(PlossN_), PlM(PlossN_);

        meanE_ == meanE0 + h;
        muE_->correct(muP);
        {
            PlP == dimensionedScalar("z", PlP.dimensions(), 0.0);
            if (Pelastic_)   { volScalarField t(PlP); Pelastic_->correct(t);   PlP += t; }
            if (Pinelastic_) { volScalarField t(PlP); Pinelastic_->correct(t); PlP += t; }
        }

        meanE_ == max(meanE0 - h, dimensionedScalar("z", dimless, 0.0));
        muE_->correct(muM);
        {
            PlM == dimensionedScalar("z", PlM.dimensions(), 0.0);
            if (Pelastic_)   { volScalarField t(PlM); Pelastic_->correct(t);   PlM += t; }
            if (Pinelastic_) { volScalarField t(PlM); Pinelastic_->correct(t); PlM += t; }
        }

        meanE_ == meanE0;      // restore before anything else reads it

        const volScalarField dMu((muP - muM)/(2.0*h));
        const volScalarField dPl((PlP - PlM)/(2.0*h));

        const volScalarField& ne = species_.numberDensity(specieIndex_);

        // D is built from THE SAME EXPRESSIONS the source uses, differentiated
        // term by term, rather than from a formula re-derived in reduced
        // variables.
        //
        // Writing it as N[dPlossN/deps - (E/N)^2 dmuN/deps] was WRONG here and
        // failed on dimensions: the evaluator returns mu ALREADY divided by
        // the gas density, not muN, so the reduced form does not match the
        // units this code carries. Differentiating the source expressions
        // keeps the dimensions right by construction, which is the point --
        // there is no second derivation to get out of step with the first.
        //
        //     S      =  mu ne |E|^2  -  ne N PlossN
    //
    // dSdEps_ IS the damping coefficient -dS/d(n_eps) = -(dS/deps)/ne, so it
    // is used DIRECTLY in fvm::SuSp. Dividing by ne there as well -- which is
    // what the first version did -- double-counts and leaves [0 0 -1] where
    // the equation needs [0 -3 -1].
        //     dS/deps = dMu ne |E|^2  -  ne N dPlossN
        //     D       = -(dS/deps)/ne
        const volScalarField dS
        (
            dMu*ne*magSqr(E_)/oneVolt - ne*species_.backgroundDensity()*dPl
        );

        // NOT clamped positive here: fvm::SuSp does that job per cell, putting
        // it implicit where it damps and explicit where it would anti-damp.
        // This is the genuinely sign-indefinite coefficient the codebase had
        // no use for until now -- see the note on why Sp is used everywhere
        // else.
        dSdEps_ = -dS/max(ne, dimensionedScalar("f", ne.dimensions(), nEfloor_));
    }
}


bool Foam::localEnergyEnergyModel::canExportEnergyCoefficients() const
{
    return muE_.valid() && muE_->hasScalarLookup()
        && Pelastic_.valid() && Pelastic_->hasScalarLookup()
        && Pinelastic_.valid() && Pinelastic_->hasScalarLookup();
}


std::function<Foam::scalar(Foam::scalar)>
Foam::localEnergyEnergyModel::muNofEps() const
{
    // Captures the EVALUATOR, not a value: the table is owned by this model
    // for the life of the run, and re-reading it per call is what lets the
    // integrator see a coefficient that responds to the eps it is currently
    // at. Returning mu*N raw -- see plasmaPropertyEvaluator::atValue().
    const plasmaPropertyEvaluator* mu = muE_.get();
    return [mu](const scalar e) { return mu->atValue(e); };
}


std::function<Foam::scalar(Foam::scalar)>
Foam::localEnergyEnergyModel::PlossNofEps() const
{
    const plasmaPropertyEvaluator* pe = Pelastic_.get();
    const plasmaPropertyEvaluator* pi = Pinelastic_.get();
    return [pe, pi](const scalar e)
    {
        return pe->atValue(e) + pi->atValue(e);
    };
}


Foam::scalar Foam::localEnergyEnergyModel::maxEnergyRelaxationRate() const
{
    // 1/tau_eps = P_loss/n_eps -- the SAME expression the equation puts on the
    // diagonal as Lrate, rebuilt here rather than cached so the limiter and
    // the matrix cannot disagree about how fast the energy relaxes.
    //
    // Before the first correct() the coefficient fields are still zero, so
    // there is no rate to report and the limiter must not be handed one --
    // the same guard as maxEnergyRate().
    if (gMax(mag(PlossN_.primitiveField())) <= 0)
    {
        return 0;
    }

    const scalarField& Pl = PlossN_.primitiveField();
    const scalarField& nE = nEps_.primitiveField();
    const scalarField& ne =
        species_.numberDensity(specieIndex_).primitiveField();
    const scalar N = species_.backgroundDensity().value();

    // EXPLICIT LOOP with a guarded denominator, for the same reason
    // maxEnergyRate() uses one: a field expression divides by n_eps, which is
    // legitimately zero in the quiescent far field before the seed spreads.
    scalar rmax = 0;
    forAll(Pl, c)
    {
        if (nE[c] > VSMALL)
        {
            rmax = max(rmax, (Pl[c]*ne[c]*N)/nE[c]);
        }
    }

    reduce(rmax, maxOp<scalar>());

    return rmax;
}


Foam::scalar Foam::localEnergyEnergyModel::maxEnergyRate() const
{
    // The SAME flux the energy equation convects with, so the limiter and the
    // equation cannot disagree: rebuilt here rather than cached, because a
    // cached copy would silently go stale against a changed mobility.
    // Before the first correct() the coefficient fields are still zero, so
    // there is no rate to report and the limiter must not be handed one.
    if (gMax(mag(muEpsf_.primitiveField())) <= 0)
    {
        return 0;
    }

    const scalar Z = species_.speciesChargeNumber(specieIndex_);

    const dimensionedScalar energyFactor
    (
        "energyFactor", dimless, energyCoeffsTabulated_ ? 1.0 : 5.0/3.0
    );

    const surfaceScalarField phiEps
    (
        energyFactor*Z
       *fvc::interpolate(muEpsf_, "interpolate(mu_e)")
       *species_.em().phiE()
    );

    // Volumetric flux, so no density division -- unlike the species rate,
    // whose particleFlux carries n and must be divided by it.
    //
    // EXPLICIT LOOP, not a field expression. `sum/(mesh_.V() + VSMALL)`
    // resolved to a scalar/field divide and raised SIGFPE at construction,
    // where this is first called: muEpsf_ is still zero there because
    // correct() has not run, so the whole expression degenerates. An explicit
    // loop has no operator ambiguity and guards the denominator directly.
    const scalarField sumFlux
    (
        fvc::surfaceSum(mag(phiEps))().primitiveField()
    );

    const scalarField& Vc = mesh_.V().field();

    scalar rmax = 0;
    forAll(sumFlux, c)
    {
        rmax = max(rmax, 0.5*sumFlux[c]/max(Vc[c], VSMALL));
    }

    reduce(rmax, maxOp<scalar>());

    return rmax;
}


tmp<fvScalarMatrix> localEnergyEnergyModel::eEqn() const
{
    const volScalarField& ne = species_.numberDensity(specieIndex_);

    // Energy flux moments of the two-term expansion: the energy drifts and
    // diffuses at 5/3 of the particle rates (Hagelaar & Kroesen Eqs. 10-13).
    // 5/3 ONLY when falling back to the particle coefficients. With the
    // tabulated energy coefficients the factor is already inside them, and
    // applying it again would double-count.
    const dimensionedScalar energyFactor
    (
        "energyFactor", dimless, energyCoeffsTabulated_ ? 1.0 : 5.0/3.0
    );

    // Built the SAME way driftDiffusion builds its convective flux:
    //
    //     Z * fvc::interpolate(mu) * em().phiE()
    //
    // rather than fvc::flux(-mu*E). Two reasons, both learned the hard way
    // here: fvc::flux of a product auto-generates the scheme name
    // `flux((-muE_lmea*E))`, which no case declares and which aborted the run;
    // and interpolating mu ALONE against the precomputed face flux phiE is
    // what the species equations already do, so the energy rides on exactly
    // the same discrete velocity field as the electrons it belongs to.
    //
    // The scheme name is borrowed from the electron mobility for the same
    // reason the laplacian is: this IS that interpolation, and a case should
    // not have to declare a scheme for a field it did not ask for.
    const scalar Z = species_.speciesChargeNumber(specieIndex_);

    const surfaceScalarField phiEps
    (
        "phiEps",
        energyFactor*Z
       *fvc::interpolate(muEpsf_, "interpolate(mu_e)")
       *species_.em().phiE()
    );

    const volScalarField DEps("DEps", energyFactor*DEf_);

    // TRANSPORT RATE, straight from the operators, at the current iterate.
    //
    // Sign follows the equation below: ddt(nEps) + div - laplacian + ... = P,
    // so transport contributes -div + laplacian to d(nEps)/dt -- the same
    // convention the species use for chemExt_, so the integrator can add it
    // to the energy row exactly as it adds theirs.
    //
    // The scheme names are BORROWED from the electron equation for the same
    // reason the implicit operators borrow them: this case declares
    // div(phi_e,n_e) and laplacian(D_e,n_e) and nothing else, so asking for
    // div(phiEps,nEps_e) here would abort a case that runs today.
    if (transportEnergy_)
    {
        energyTransportRate_ =
        (
            -fvc::div(phiEps, nEps_, "div(phi_e,n_e)")()
           + fvc::laplacian(DEps, nEps_, "laplacian(D_e,n_e)")()
        )().primitiveField();
    }
    else
    {
        energyTransportRate_.setSize(mesh_.nCells(), Zero);
    }

    // JOULE HEATING, explicit in phase A: -e Gamma_e . E, which for electrons
    // drifting against the field is a positive power input. Written from the
    // drift flux rather than from the full Gamma_e, because the diffusive part
    // carries E.grad(n_e) -- the sign-indefinite term that makes Hagelaar's
    // bracket not unconditionally damping.
    const volScalarField jouleHeating
    (
        "jouleHeating",
        muEf_*ne*magSqr(E_)/oneVolt
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

    // ---- WHERE THE SOURCE COMES FROM -----------------------------------
    //
    // OPTION 4: when `energySource chemistry` is set, the per-cell stiff
    // solver has ALREADY integrated n_eps alongside the species and returned
    // the end-state production and loss RATE. Using them here instead of the
    // locally built jouleHeating/Lrate is not an optimisation -- applying both
    // would count the source twice.
    //
    // The Newton sensitivity dSdEps_ is dropped in that mode as well: it
    // exists to linearise a source this equation is no longer evaluating, and
    // the integrator resolved the nonlinearity properly rather than about a
    // tangent that was MEASURED going to zero in cold cells.
    // BY NAME, not by type: plasmaTransport links plasmaEnergy, so depending
    // on it here would be a circular library dependency. The two fields are
    // registered by plasmaTransport exactly when it is integrating the energy,
    // so their presence IS the signal.
    const bool chemOwnsSource =
        mesh_.foundObject<volScalarField>("chemPeps")
     && mesh_.foundObject<volScalarField>("chemLeps");

    volScalarField Psrc(jouleHeating);
    volScalarField Lsrc(Lrate);

    if (chemOwnsSource)
    {
        Psrc.primitiveFieldRef() =
            mesh_.lookupObject<volScalarField>("chemPeps").primitiveField();
        Lsrc.primitiveFieldRef() =
            mesh_.lookupObject<volScalarField>("chemLeps").primitiveField();

        if (!chemSourceReported_)
        {
            chemSourceReported_ = true;
            Info<< "  LMEA: the electron energy SOURCE comes from the"
                << " chemistry ODE (energySource chemistry)." << nl
                << "    Joule and P_loss are integrated per cell with the"
                << " species and returned as an end-state P/L pair; this"
                << " model contributes transport only." << endl;
        }
    }

    tmp<fvScalarMatrix> tEqn
    (
        new fvScalarMatrix
        (
            fvm::ddt(nEps_)

            // Scheme names BORROWED from the electron equation, checked
            // against the case before running rather than discovered by
            // crashing: this case declares div(phi_e,n_e) and
            // laplacian(D_e,n_e) with no `default`, so the auto-generated
            // div(phiEps,nEps_e) and laplacian(DEps,nEps_e) abort at start-up.
            //
            // Borrowing is also the right physics: the energy is convected by
            // the same velocity field and diffuses on the same mesh as the
            // electrons, at 5/3 of the rates. A case should not have to
            // declare schemes for an equation it enabled with one keyword.
          + (transportEnergy_
              ? fvm::div(phiEps, nEps_, "div(phi_e,n_e)")
              : fvm::Sp(dimensionedScalar("z", dimless/dimTime, 0.0), nEps_))
          - (transportEnergy_
              ? fvm::laplacian(DEps, nEps_, "laplacian(D_e,n_e)")
              : fvm::Sp(dimensionedScalar("z", dimless/dimTime, 0.0), nEps_))
          + fvm::Sp(Lsrc, nEps_)
          + (chemOwnsSource
              ? fvm::Sp(dimensionedScalar("z", dimless/dimTime, 0.0), nEps_)
              : fvm::SuSp(dSdEps_, nEps_))
         ==
            Psrc
          + (chemOwnsSource
              ? 0.0*dSdEps_*nEps_
              : dSdEps_*nEps_)
        )
    );

    // ---- PER-TERM DUMP AT THE WORST CELL --------------------------------
    //
    // Every wrong diagnosis in this model came from inferring a mechanism
    // from aggregate numbers. These are the actual terms the matrix receives,
    // at the cell with the largest eps = n_eps/n_e, so the candidates are
    // separated by measurement rather than by argument:
    //
    //   * a near-zero DIAGONAL explains the SIGFPE directly;
    //   * a large div/laplacian against a small Joule means the energy is
    //     being DELIVERED to the cell (boundary or neighbour), which is the
    //     zeroGradient-at-the-electrode story;
    //   * a large Joule against a small transport means it is generated
    //     locally, and the wall condition is not to blame.
    //
    // All reduces are OUTSIDE any per-cell loop, and the payload comes from
    // the owning rank only -- a rank-local index would report the wrong cell,
    // which is exactly the bug this diagnostic replaces.
    if (reportInterval_ > 0 && mesh_.time().timeIndex() % reportInterval_ == 0)
    {
        const volScalarField& neR = species_.numberDensity(specieIndex_);
        const scalarField& nER = nEps_.primitiveField();
        const scalarField& neI2 = neR.primitiveField();

        label iW = -1;
        scalar epsW = -GREAT;
        forAll(nER, c)
        {
            const scalar e = nER[c]/max(neI2[c], nEfloor_);
            if (e > epsW) { epsW = e; iW = c; }
        }

        const scalar epsG = returnReduce(epsW, maxOp<scalar>());
        const label own = returnReduce
        (
            (iW >= 0 && epsW >= epsG) ? Pstream::myProcNo() : labelMax,
            minOp<label>()
        );

        // Explicit evaluations of the SAME expressions the matrix was built
        // from, so the two cannot drift apart.
        const volScalarField ddtT(fvc::ddt(nEps_));

        volScalarField divT
        (
            IOobject("divT_lmea", mesh_.time().timeName(), mesh_,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh_,
            dimensionedScalar(nEps_.dimensions()/dimTime, Zero)
        );
        volScalarField lapT(divT);

        if (transportEnergy_)
        {
            divT = fvc::div(phiEps, nEps_, "div(phi_e,n_e)");
            lapT = fvc::laplacian(DEps, nEps_, "laplacian(D_e,n_e)");
        }
        const scalarField& diagF = tEqn->diag();
        const scalarField& VF = mesh_.V();

        // REPORT WHAT THE MATRIX RECEIVED, not what this model would have
        // built. Under `energySource chemistry` the source terms are REPLACED
        // by Psrc/Lsrc above, so reporting Lrate/jouleHeating here described a
        // matrix that was never assembled -- the exact failure the comment at
        // the top of this block warns against.
        //
        // The LOSS-RATE GAP (q[9]..q[11]) is the quantity under test:
        // chemLeps is L/nEps_end, but the equation applies it to the nEps it
        // SOLVES for, so the realised loss is L*(nEps_new/nEps_end). If that
        // ratio departs from 1 by O(dt) -- halving as dt halves -- the loss
        // term is the first-order defect. If it does not scale, the candidate
        // is refuted.
        const scalarField* nEpsEndF = nullptr;
        if (mesh_.foundObject<volScalarField>("chemNEpsEnd"))
        {
            nEpsEndF =
               &mesh_.lookupObject<volScalarField>("chemNEpsEnd")
                .primitiveField();
        }
        const scalarField* epsEndF = nullptr;
        if (mesh_.foundObject<volScalarField>("chemEpsEnd"))
        {
            epsEndF =
               &mesh_.lookupObject<volScalarField>("chemEpsEnd")
                .primitiveField();
        }

        // TERM-BY-TERM DIFF OF THE TWO SOURCE PATHS.
        //
        // `energySource model` measures p = 2.7 and `chemistry` p = 1.05 from
        // algebra that is IDENTICAL on paper. Rather than guess which term
        // differs -- eight hypotheses have been tried, four refuted by
        // measurement -- both paths are evaluated at the same cell and the
        // same instant, and printed side by side:
        //
        //   eps    : the model's meanE FIELD vs the ODE's integrated end state
        //   mu     : the SAME table, looked up at those two eps
        //   J, L   : the assembled terms
        //
        // Whichever line shows the O(dt) discrepancy IS the defect; the ones
        // that agree are eliminated without further argument.
        scalarList q(18, Zero);
        if (Pstream::myProcNo() == own && iW >= 0)
        {
            q[0] = epsW;
            q[1] = ddtT[iW];
            q[2] = divT[iW];
            q[3] = lapT[iW];
            q[4] = Lsrc[iW]*nEps_[iW];
            q[5] = chemOwnsSource ? 0.0 : dSdEps_[iW]*nEps_[iW];
            q[6] = Psrc[iW];
            q[7] = diagF[iW]/max(VF[iW], VSMALL);  // per unit volume
            q[8] = neI2[iW];
            q[9] = nEpsEndF ? (*nEpsEndF)[iW] : 0.0;
            q[10] = nEps_[iW];
            // L itself = rate * the state it was normalised BY.
            q[11] = nEpsEndF ? Lsrc[iW]*(*nEpsEndF)[iW] : 0.0;

            // --- MODEL path: evaluated at the meanE FIELD ---
            q[12] = epsW;                       // eps_model = nEps/n_e
            q[13] = muEf_[iW];                  // mu at the field eps
            q[14] = jouleHeating[iW];           // J_model
            q[15] = Lrate[iW]*nEps_[iW];        // L_model (absolute power)

            // --- ODE path: evaluated at the integrated end state ---
            q[16] = epsEndF ? (*epsEndF)[iW] : 0.0;
            // The SAME table the ODE used, looked up here at the ODE's eps,
            // so a difference in q[13] vs q[17] is the eps, not the table.
            q[17] = (epsEndF && (*epsEndF)[iW] > 0)
                  ? muNofEps()((*epsEndF)[iW])
                     /species_.backgroundDensity().value()
                  : 0.0;
        }
        forAll(q, i) reduce(q[i], sumOp<scalar>());

        Info<< "    LMEA terms at worst cell (eps " << q[0]
            << " eV, n_e " << q[8] << "):" << nl
            << "      ddt " << q[1] << "  div " << q[2]
            << "  lap " << q[3] << nl
            << "      Sp(L)*nEps " << q[4] << "  SuSp*nEps " << q[5]
            << "  Joule " << q[6] << nl
            << "      diag/V " << q[7]
            << "   (a diag/V near zero is the SIGFPE; transport >> Joule"
               " means the energy is DELIVERED, not generated)" << endl;

        if (chemOwnsSource && q[9] > VSMALL)
        {
            const scalar ratio = q[10]/q[9];
            Info<< "      LOSS-RATE GAP: nEps_end " << q[9]
                << "  nEps_solved " << q[10] << nl
                << "        ratio " << ratio
                << "  rel.gap " << (ratio - 1.0)
                << "   L " << q[11] << "  realised " << q[4] << nl
                << "        (rel.gap ~ O(dt) and halving with dt =>"
                   " the loss term is the first-order defect)" << endl;

            const scalar dEps = q[16] > 0 ? (q[12] - q[16])/q[16] : 0.0;
            const scalar dMu  = q[17] > 0 ? (q[13] - q[17])/q[17] : 0.0;
            const scalar dJ   = q[6]  != 0 ? (q[14] - q[6])/q[6]  : 0.0;
            const scalar dL   = q[11] != 0 ? (q[15] - q[11])/q[11] : 0.0;

            Info<< "      PATH DIFF (model - ODE)/ODE at the same cell:" << nl
                << "        eps  model " << q[12] << "  ode " << q[16]
                << "   rel " << dEps << nl
                << "        mu   model " << q[13] << "  ode " << q[17]
                << "   rel " << dMu << nl
                << "        J    model " << q[14] << "  ode " << q[6]
                << "   rel " << dJ << nl
                << "        L    model " << q[15] << "  ode " << q[11]
                << "   rel " << dL << nl
                << "        (the line whose rel is O(dt) IS the defect;"
                   " lines that agree are eliminated)" << endl;
        }
    }

    return tEqn;
}


void localEnergyEnergyModel::discardStep()
{
    // THE STEP IS BEING THROWN AWAY, so n_eps must go back with the species.
    //
    // plasmaTransport::discardStep() restores every number density to
    // oldTime(), but n_eps lives here, and nothing restored it. A retried
    // step therefore re-solved with n_e at t^n and n_eps at the DISCARDED
    // t^{n+1}, and the ratio meanE = n_eps/n_e was formed from two different
    // instants. Measured: SIGFPE on the first discard of an LMEA run
    // (`retry-iso`, step 6, one discard, zero non-convergences).
    //
    // retryStep was developed and validated on LFA cases only -- none of the
    // four co1.5-retry* beds declares an energyModel -- so this path had never
    // been exercised with an electron-energy equation.
    nEps_ == nEps_.oldTime();
    nEps_.correctBoundaryConditions();

    // Everything downstream of n_eps -- meanE, T_e, the clamp counters -- is
    // derived, so it must be rebuilt rather than left describing the attempt
    // that was discarded.
    updateDerived();

    // Per-pass diagnostics, for the same reason plasmaTransport zeroes its
    // own: the discarded attempt must not be counted alongside the kept one.
    clampedSteps_ = 0;
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
