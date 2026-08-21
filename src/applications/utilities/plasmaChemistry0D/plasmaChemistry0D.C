/*---------------------------------------------------------------------------*\
Application
    plasmaChemistry0D

Description
    Integrate the compiled mechanism in a single well-mixed cell and write the
    history as CSV.

    Exists to validate the chemistry ODE against Cantera on the identical
    mechanism, in a setting where a disagreement can only come from the
    chemistry. There is no mesh, no transport, no field: everything that a
    streamer run would confound with the chemistry is absent.

    Rate coefficients for the electron-impact reactions are read from the same
    tables the CFD run interpolates, at one fixed E/N, so the EEDF is identical
    on both sides by construction and only the integration differs.

Usage
    plasmaChemistry0D -mechanism <f.foam> -tables <dir> -EN <Td> -endTime <s>

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "OFstream.H"
#include "dictionary.H"
#include "plasmaChemistry.H"

// SoEEDF, for re-solving the EEDF in process.
#undef Log            // OpenFOAM defines this; SoEEDF's headers do not want it
#include "Mechanism.H"
#include "BoltzmannSolver.H"
#include "TransportCoefficients.H"
#include "DynamicList.H"
#include "wallLoss.H"
#include "vibRelax.H"
#include "janafMixture.H"
#include <cstdlib>
#include <cmath>

using namespace Foam;


// Vibrational-translational relaxation and the mixture thermodynamics both
// live in headers so they can be unit-tested without a reactor: see
// vibRelax.H (tested by testVibRelax) and janafMixture below.
using Foam::vibRelax::tauVT_N2;

static Foam::scalar cvAirFallback(const Foam::scalar T)
{
    const scalar t = min(max(T, scalar(200)), scalar(3500));
    return 6.165129e+02 + 3.212567e-01*t
         - 8.907555e-05*t*t + 8.870573e-09*t*t*t;
}






// Read an OpenFOAM ((x y) ...) table and interpolate linearly.
static scalar tableAt(const fileName& path, const scalar x)
{
    IFstream is(path);
    if (!is.good())
    {
        FatalErrorInFunction << "cannot open " << path << exit(FatalError);
    }
    // Parsed with a tolerant scan rather than a stream operator: the file
    // carries a `//` provenance header, and the stream operator raises a
    // FatalIOError on a token it cannot read instead of skipping it.
    DynamicList<scalar> xs, ys;
    string line;
    while (is.good())
    {
        is.getLine(line);
        const auto o = line.find('(');
        if (o == std::string::npos) continue;
        const char* p0 = line.c_str() + o + 1;
        char* end = nullptr;
        const double a = std::strtod(p0, &end);
        if (end == p0) continue;
        const char* p1 = end;
        const double b = std::strtod(p1, &end);
        if (end == p1) continue;
        xs.append(a); ys.append(b);
    }
    if (xs.size() < 2) return 0;
    if (x <= xs.first()) return ys.first();
    if (x >= xs.last())  return ys.last();
    for (label i = 1; i < xs.size(); ++i)
    {
        if (xs[i] >= x)
        {
            const scalar t = (x - xs[i-1])/(xs[i] - xs[i-1]);
            return ys[i-1] + t*(ys[i] - ys[i-1]);
        }
    }
    return ys.last();
}


// Integrate one substep, carrying the electron energy density in the state
// vector when Option 4 is active.
//
// `n` stays the SPECIES vector everywhere else in this file; the padded copy
// lives only for the duration of the call. Resizing `n` itself would mean
// every output loop, every chemistry call and every energy-budget sum had to
// know about the extra slot.
static void integrateWithEnergy
(
    const Foam::plasmaChemistry& chem,
    Foam::scalarField& n,
    Foam::scalar& nEps,
    const Foam::scalarField& kTab,
    const Foam::scalar T,
    const Foam::scalar dt,
    const bool withEnergy,
    const Foam::scalar Emag,
    const Foam::scalar Ngas
)
{
    if (!withEnergy)
    {
        chem.integrate(n, kTab, T, dt);
        return;
    }

    chem.setEnergyCell(Emag, Ngas);

    const Foam::label nSp = n.size();
    Foam::scalarField y(nSp + 1);
    for (Foam::label i = 0; i < nSp; ++i) y[i] = n[i];
    y[nSp] = nEps;

    chem.integrate(y, kTab, T, dt);

    for (Foam::label i = 0; i < nSp; ++i) n[i] = y[i];
    nEps = Foam::max(y[nSp], Foam::scalar(0));
}

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("mechanism", "file", "compiled mechanism .foam");
    argList::addOption("tables", "dir", "rate table directory");
    argList::addOption("EN", "Td", "reduced field, held fixed");
    argList::addOption("endTime", "s", "integration time");
    argList::addOption("nOut", "n", "output samples");
    argList::addOption("Tgas", "K", "gas temperature");
    argList::addOption("pressure", "Pa", "gas pressure");
    argList::addOption("ne0", "m^-3", "initial electron density");
    argList::addOption("out", "file", "CSV output");

    // ---- gas heating -------------------------------------------------------
    argList::addBoolOption("gasHeating",
        "solve the gas temperature from the deposited power");
    argList::addOption("ENpulse", "peak:centre_ns:fwhm_ns",
        "Gaussian E/N pulse instead of a fixed field");
    argList::addOption("manifest", "file",
        "<mech>.mech.json. Supplying it enables DYNAMIC EEDF: the Boltzmann"
        " equation is re-solved at the current composition and gas temperature"
        " instead of interpolating tables frozen at the initial ones");
    argList::addOption("eedfTol", "frac",
        "re-solve when any mole fraction or T has drifted by this relative"
        " amount since the last solve (default 0.05)");
    argList::addOption("X", "N2=0.774 O2=0.186 O=0.04",
        "initial mole fractions, overriding the mechanism's own composition."
        " Needed to reproduce an experiment that starts partly dissociated by"
        " earlier pulses");
    argList::addOption("profile", "file",
        "CSV of t_ns,EN_Td,ne_cm3 -- measured discharge history, which is how a"
        " 0-D reactor is compared with an experiment");
    argList::addOption("wallLoss",
        "none|ambipolar|effectiveDiffusion|QGM|hFactorLP|hFactorHPEN",
        "charged-particle loss to the walls, for a BOUNDED discharge"
        " (default none). 0-D ONLY: the coupled CFD resolves this transport"
        " explicitly, so using both would count the same loss twice.");
    argList::addOption("wallR", "m", "discharge radius for -wallLoss");
    argList::addOption("wallL", "m", "discharge length for -wallLoss");
    argList::addOption("wallMuiN", "1/(V m s)",
        "reduced ion mobility for -wallLoss (default 5e21, ~O2+ in O2)");
    // ---- LMEA (transported electron energy) --------------------------------
    argList::addBoolOption("lmea",
        "solve the LMEA electron energy equation instead of taking the mean "
        "energy from the LFA table");
    argList::addOption("lmeaSource", "explicit|implicit|newton|ode",
        "how the source enters the energy update (default implicit). "
        "`newton` adds the SENSITIVITY dS/deps -- Hagelaar's linearisation -- "
        "not merely the loss magnitude.");
    argList::addOption("lmeaDt", "s",
        "fixed energy sub-step; default follows the chemistry step");

    argList::addOption("chemistryBackend", "native|cantera",
        "who evaluates the heavy reactions. Same keyword as the CFD.");
    argList::addOption("chemistrySource",
        "ode|implicitRate|adaptive|adaptiveError",
        "how the chemistry is advanced. Same keyword as the CFD (default ode)."
        " In 0-D there is no transport and therefore no operator splitting, so"
        " `ode` here is the converged REFERENCE rather than a first-order"
        " alternative as it is in the CFD.");
    argList::addOption("nOuterCorrectors", "n",
        "Picard iterations per implicitRate step. Mirrors the CFD keyword;"
        " 1 is explicit in the production term and can diverge (default 4).");
    argList::addOption("changeTol", "x",
        "adaptive: relative change of any significant species allowed in one"
        " implicitRate step before it is rejected for the stiff substep."
        " Default 0.1: 0.5 was measured to be too loose, letting a single"
        " accepted step leave an afterglow 190x off.");
    argList::addOption("stiffTol", "x",
        "adaptive threshold on max(L)*dt; above it the stiff substep is used"
        " (default 1).");
    argList::addOption("errRelTol", "x",
        "adaptiveError: relative part of the local-error tolerance"
        " (default 1e-2).");
    argList::addOption("errAbsFrac", "x",
        "adaptiveError: absolute part of the local-error tolerance, as a"
        " fraction of the cell's electron density (default 1e-4). Together:"
        " scale_s = errAbsFrac*n_e + errRelTol*n_s. Without the absolute part"
        " a trace species moving orders of magnitude outranks a dominant"
        " species moving a few percent.");
    argList::addBoolOption("constPressure",
        "hold pressure instead of volume: the gas expands as it heats, so the"
        " temperature rises against c_p rather than c_v. Right past the"
        " acoustic time; wrong within a nanosecond pulse, which is isochoric.");
    argList::addOption("tauVT", "s",
        "vibrational-translational relaxation time (default 1e-5 s)");
    argList args(argc, argv, false, false, false);

    const fileName mechFile = args.getOrDefault<fileName>
        ("mechanism", "constant/air_plasma.foam");
    const fileName tableDir = args.getOrDefault<fileName>
        ("tables", "constant/plasmaTables");
    const scalar EN_Td   = args.getOrDefault<scalar>("EN", 150.0);
    const bool   heating = args.found("gasHeating");

    // LMEA. The point of doing this in 0-D is that the answer is KNOWN: at
    // fixed E/N the steady mean energy must return meanEnergy_vs_reducedE,
    // the LFA table, because that is the local-equilibrium limit LMEA reduces
    // to. Getting it wrong here costs seconds; getting it wrong at 1.15M
    // cells costs a day.
    const bool lmea = args.found("lmea");
    const word lmeaSrc = args.getOrDefault<word>("lmeaSource", "implicit");
    const bool lmeaImplicit = (lmeaSrc == "implicit");
    const bool lmeaNewton   = (lmeaSrc == "newton");

    // OPTION 4: integrate n_eps INSIDE the stiff chemistry ODE rather than
    // updating it by hand after the chemistry. See docs/lmea-option4-plan.md.
    const bool lmeaOde      = (lmeaSrc == "ode");
    if (lmeaOde && !lmea)
    {
        FatalErrorInFunction
            << "-lmeaSource ode integrates the electron energy in the"
               " chemistry ODE, but -lmea was not given."
            << exit(FatalError);
    }

    // -lmeaDt PINS the integration step, overriding the pulse-derived dtFine
    // and the afterglow growth. It exists to map the LMEA stability boundary
    // in dt, which is impossible while the step is chosen for you.
    //
    // It was declared here and NEVER READ for its first several commits: an
    // 18-run dt scan came back bit-identical because every run silently used
    // dtFine = fwhm/50. A dead option is worse than a missing one -- it
    // reports success while measuring nothing.
    const scalar lmeaDtFix = args.getOrDefault<scalar>("lmeaDt", 0.0);
    if (lmeaDtFix > 0 && !lmea)
    {
        FatalErrorInFunction
            << "-lmeaDt pins the step for the LMEA energy integration but"
               " -lmea was not given, so there is nothing to pin."
            << exit(FatalError);
    }

    // A PRESCRIBED relaxation time, not a Landau-Teller model, and that is a
    // deliberate first step. Within a nanosecond pulse tau_VT is microseconds,
    // so the reservoir only fills and its emptying rate barely matters; what
    // matters is that the energy is REMOVED from the gas budget rather than
    // counted as prompt heat. A physics-based tau (Millikan-White, with the
    // atomic-oxygen correction that moves it by decades once O2 dissociates)
    // is the tier-2 model in docs/gas-heating-plan.md.
    // Negative means COMPUTE it per step from Millikan-White plus the
    // atomic-oxygen channel; a positive value forces a fixed one, which is
    // still useful for isolating the reservoir from the relaxation model.
    const scalar tauVTfixed = args.getOrDefault<scalar>("tauVT", -1.0);

    // Gaussian pulse "peak:centre_ns:fwhm_ns", or a fixed field.
    scalar pkEN = EN_Td, tc_ns = 0.0, fwhm_ns = 0.0;
    const bool pulsed = args.found("ENpulse");
    if (pulsed)
    {
        const std::string spec(args.get<string>("ENpulse"));
        DynamicList<std::string> f;
        std::size_t a = 0, b;
        while ((b = spec.find(':', a)) != std::string::npos)
        {
            f.append(spec.substr(a, b - a)); a = b + 1;
        }
        f.append(spec.substr(a));
        if (f.size() != 3)
        {
            FatalErrorInFunction << "-ENpulse wants peak:centre_ns:fwhm_ns"
                                 << exit(FatalError);
        }
        pkEN = std::atof(f[0].c_str());
        tc_ns = std::atof(f[1].c_str());
        fwhm_ns = std::atof(f[2].c_str());
    }
    // Measured discharge history, if given. A pin-to-pin experiment supplies
    // E/N(t) from voltage and gap, and n_e(t) from Stark broadening; imposing
    // both is what makes the comparison a test of the CHEMISTRY rather than of
    // a discharge model we do not have in 0-D.
    DynamicList<scalar> pT, pEN, pNE;
    const bool profiled = args.found("profile");
    if (profiled)
    {
        IFstream pf(args.get<fileName>("profile"));
        if (!pf.good())
        {
            FatalErrorInFunction << "cannot open profile" << exit(FatalError);
        }
        string line;
        while (pf.good())
        {
            pf.getLine(line);
            if (line.empty() || line[0] == '#' || line[0] == 't') continue;
            std::string l(line);
            for (auto& c : l) if (c == ',') c = ' ';
            IStringStream is(l);
            scalar a, b, c;
            is >> a >> b >> c;
            pT.append(a*1e-9); pEN.append(b); pNE.append(c*1e6);
        }
        Info<< "profile: " << pT.size() << " points, "
            << pT[0]*1e9 << " to " << pT[pT.size()-1]*1e9 << " ns" << endl;
    }

    auto interpAt = [&](const DynamicList<scalar>& y, const scalar t) -> scalar
    {
        if (t <= pT[0]) return y[0];
        if (t >= pT[pT.size()-1]) return y[y.size()-1];
        label i = 0;
        while (i + 1 < pT.size() && pT[i+1] < t) ++i;
        const scalar w = (t - pT[i])/(pT[i+1] - pT[i]);
        return y[i] + w*(y[i+1] - y[i]);
    };

    auto fieldAt = [&](const scalar t) -> scalar
    {
        if (profiled) return interpAt(pEN, t);
        if (!pulsed) return EN_Td;
        const scalar sig = fwhm_ns*1e-9/(2.0*Foam::sqrt(2.0*Foam::log(2.0)));
        const scalar x = (t - tc_ns*1e-9)/sig;
        const scalar en = pkEN*Foam::exp(-0.5*x*x);

        // The field goes to ZERO after the pulse, and it matters that it does.
        // An earlier version floored it at 0.1 Td to keep table lookups inside
        // range -- and that floor kept depositing into the vibrational
        // reservoir for the whole 100 us afterglow, because at low E/N nearly
        // all electron energy goes to vibration. The apparent vibrational
        // share climbed from 47% to 66% on nothing but the floor.
        return (en < 1e-4*pkEN) ? 0.0 : en;
    };
    const scalar endTime = args.getOrDefault<scalar>("endTime", 1e-8);
    const label  nOut    = args.getOrDefault<label>("nOut", 100);
    const scalar Tgas    = args.getOrDefault<scalar>("Tgas", 300.0);
    const scalar pres    = args.getOrDefault<scalar>("pressure", 1e5);
    const scalar ne0     = args.getOrDefault<scalar>("ne0", 1e16);
    const fileName out   = args.getOrDefault<fileName>("out", "chem0d.csv");

    const scalar nGas = pres/(1.380649e-23*Tgas);

    // Every species the mechanism names is carried, so nothing is held fixed
    // and the comparison with Cantera is like for like.
    dictionary mech((IFstream(mechFile))());
    const dictionary& chargeDict = mech.subDict("speciesCharge");
    const word electron = mech.getOrDefault<word>("electronSpecies", "Electron");

    wordList species;
    DynamicList<scalar> charge;
    for (const entry& e : chargeDict)
    {
        species.append(e.keyword());
        charge.append(readScalar(e.stream()));
    }

    dictionary cfg;
    cfg.add("electronName", electron);
    cfg.add("backgroundDensity", nGas);

    // MIRRORS THE CFD. plasmaTransport selects who evaluates the heavy
    // reactions with exactly this keyword, and a 0-D tool that could not
    // reproduce the production configuration would be a debugging tool you
    // cannot debug with: a backend disagreement would only ever show up in
    // the run that is hardest to interrogate.
    const word backend = args.getOrDefault<word>("chemistryBackend", "native");
    cfg.add("chemistryBackend", backend);

    plasmaChemistry chem(mechFile, species, scalarField(charge), cfg);

    // Electron-impact rate coefficients at this fixed E/N, from the same
    // tables the CFD run reads.
    const scalar EN_SI = EN_Td*1e-21;
    scalarField kTab(chem.nTabulated(), Zero);
    forAll(chem.tabulatedIds(), i)
    {
        kTab[i] = tableAt(tableDir/("k_" + chem.tabulatedIds()[i]
                                    + "_vs_reducedE"), EN_SI);
    }

    // ---- dynamic EEDF ------------------------------------------------------
    // With a manifest, the Boltzmann equation is re-solved during the run at
    // the CURRENT composition and gas temperature. Tables are a functional of
    // the whole mixture, and a discharge changes the mixture: dissociating O2
    // removes vibrational and Herzberg energy sinks, the EEDF runs hotter, and
    // every threshold rate accelerates. Measured on the Rusterholtz benchmark,
    // N2 ionisation at 100 Td is 84% higher at the end of the run than the
    // frozen table says. That is a feedback loop -- dissociation begets faster
    // ionisation -- and no table can carry it.
    //
    // Re-solved on a TOLERANCE rather than every step: one solve costs
    // milliseconds, which is far more than one RHS evaluation, and the EEDF
    // responds to composition that moves slowly compared with the timestep.
    // The trigger is the same idea Shao et al. use in ChemPlasKin.
    std::unique_ptr<Boltzmann::Mechanism> bmech;
    std::unique_ptr<Boltzmann::BoltzmannSolver> bsolver;
    std::unique_ptr<Boltzmann::TransportCalculator> bcalc;
    List<label> idToRate(chem.nTabulated(), -1);

    // Unit conversion per process, and it is NOT cosmetic. A `density-cm3`
    // process -- LXCat's three-body attachment -- is tabulated as a true
    // third-order k3 [m^6/s], with the collider density divided back out,
    // because the chemistry multiplies by that density itself. The solver
    // returns a second-order k [m^3/s] with the density already inside the
    // scaled cross section.
    //
    // Taking the solver's value raw made the attachment rate 1.03e24 times too
    // large, and the electron density was annihilated inside one timestep. The
    // table writer applies exactly this factor (MechTables.C); the dynamic path
    // has to apply it too, or the two disagree by twenty-four orders of
    // magnitude while looking like the same quantity.
    scalarField rateConv(chem.nTabulated(), 1.0);
    const bool dynamicEEDF = args.found("manifest");
    const scalar eedfTol = args.getOrDefault<scalar>("eedfTol", 0.05);
    const bool isobaric = args.found("constPressure");

    enum chemSourceType { csODE, csImplicitRate, csAdaptive, csAdaptiveError };
    const word csName =
        args.getOrDefault<word>("chemistrySource", "ode");
    const chemSourceType chemSource =
        (csName == "ode")           ? csODE
      : (csName == "implicitRate")  ? csImplicitRate
      : (csName == "adaptive")      ? csAdaptive
      : (csName == "adaptiveError") ? csAdaptiveError
      : (FatalErrorInFunction
            << "unknown -chemistrySource " << csName << nl
            << "Valid: ode | implicitRate | adaptive | adaptiveError"
            << exit(FatalError), csODE);
    const scalar stiffTol = args.getOrDefault<scalar>("stiffTol", 1.0);
    const scalar changeTol = args.getOrDefault<scalar>("changeTol", 0.1);

    // `adaptiveError`: the CFD's local-error switch (chemistry/chemErrorRelTol
    // and chemErrorAbsFrac in plasmaTransport), brought here because 0-D is
    // where it can be checked against TRUTH. There is no transport in 0-D and
    // therefore no operator splitting, so `-chemistrySource ode` carries no
    // splitting error and IS the converged reference -- something the CFD has
    // no equivalent of, where every reference is itself an approximation.
    //
    // Same mixed scale as the CFD. There the absolute part is a fraction of the
    // DOMAIN-maximum electron density; a 0-D cell is the whole domain, so it
    // degenerates to a fraction of this cell's n_e. That is the same
    // discrimination -- a trace species measured against the dominant one --
    // rather than a different criterion.
    const scalar errRelTol  = args.getOrDefault<scalar>("errRelTol", 1.0e-2);
    const scalar errAbsFrac = args.getOrDefault<scalar>("errAbsFrac", 1.0e-4);

    // WALL LOSS. Default none, and that is not laziness: every case validated
    // here is a nanosecond atmospheric discharge where the wall timescale is
    // ~8.5 us against a 100 ns run, so switching a model on by default would
    // move validated results by ~1% for no physical reason. See wallLoss.H.
    const word wlName = args.getOrDefault<word>("wallLoss", "none");
    const wallLoss::model wlModel = wallLoss::modelFromWord(wlName);
    if (wlName != "none" && wlModel == wallLoss::wlNone)
    {
        FatalErrorInFunction
            << "unknown -wallLoss " << wlName << nl
            << "Valid: none | ambipolar | effectiveDiffusion | QGM |"
            << " hFactorLP | hFactorHPEN" << exit(FatalError);
    }
    wallLoss::state wlState;
    wlState.R = args.getOrDefault<scalar>("wallR", 0.01);
    wlState.L = args.getOrDefault<scalar>("wallL", 0.1);
    wlState.muiN = args.getOrDefault<scalar>("wallMuiN", 5.0e21);
    if (wlModel != wallLoss::wlNone && !args.found("wallR"))
    {
        WarningInFunction
            << "-wallLoss " << wlName << " with no -wallR/-wallL given;"
            << " using R = " << wlState.R << " m, L = " << wlState.L
            << " m." << nl
            << "    The loss frequency scales as 1/Lambda^2, so the geometry"
            << " is not a detail." << endl;
    }
    // Densities below this are numerical dust: a trace species going from
    // 1e-6 to 1e-3 m^-3 is a 1000x relative change and no physics at all.
    const scalar nFloorRel = 1.0e6;
    const label nOuterCorr =
        max(label(1), args.getOrDefault<label>("nOuterCorrectors", 4));
    label nODEsteps = 0, nRateSteps = 0, nRejected = 0;
    label nSolves = 0, nUnconverged = 0;

    // Transport and power channels from the most recent dynamic solve. Negative
    // means "no solve yet" -- the tables are used until the first one lands.
    // Multiplicative corrections to the tabulated channels, from the most
    // recent dynamic solve. 1 means "table as written" -- which is exactly
    // right before the first solve, and for a run with no dynamic EEDF at all.
    scalar muNcorr = 1.0, PelCorr = 1.0, PgsCorr = 1.0, PvbCorr = 1.0;
    scalarField kCorr(chem.nTabulated(), 1.0);

    if (dynamicEEDF)
    {
        // Gas density in cm^-3, for the density-scaled (three-body) processes.
        // It must be p/(k_B T) at the CASE temperature, not Loschmidt's number
        // scaled by pressure: at 1500 K those differ by 5.5x, and three-body
        // attachment scales linearly with it. Getting this wrong made
        // attachment 5.5x too strong and annihilated the electron density
        // mid-pulse -- a collapse that looked like a stiff-chemistry failure
        // and was a units error.
        //
        // Constant through the run because the reactor is isochoric.
        bmech.reset(new Boltzmann::Mechanism(
            args.get<fileName>("manifest"), nGas*1.0e-6));

        // Map our tabulated ids onto the solver's rate-coefficient vector.
        // Done ONCE, by string id: the two orderings are different and pairing
        // them by position is exactly the class of error the mechanism hash
        // exists to prevent.
        forAll(chem.tabulatedIds(), i)
        {
            const int mi = bmech->indexOf(chem.tabulatedIds()[i]);
            idToRate[i] = (mi >= 0) ? bmech->mechToRate(mi) : -1;
            if (mi >= 0)
            {
                const auto& mp = bmech->processes()[mi];
                if (mp.scaling == "density-cm3" && mp.sigmaScale > 0)
                {
                    rateConv[i] = 1.0e-6/mp.sigmaScale;   // CM3_TO_M3/sigmaScale
                }
            }
        }
        label nMapped = 0;
        forAll(idToRate, i) if (idToRate[i] >= 0) ++nMapped;
        Info<< "dynamic EEDF: " << nMapped << " of " << chem.nTabulated()
            << " tabulated processes mapped, tolerance " << eedfTol << endl;
    }

    // Re-solve and refill kTab at this (E/N, T, composition).
    auto solveEEDF = [&](const scalar en_Td, const scalar Tg,
                         const scalarField& dens)
    {
        std::map<std::string, double> X;
        scalar tot = 0;
        forAll(species, si) if (dens[si] > 0) tot += dens[si];
        forAll(species, si)
        {
            // Only species the Boltzmann database actually carries; the rest
            // (ions, excited states without cross sections) are not part of
            // the mixture it can represent.
            bool known = false;
            for (const auto& gs : bmech->db().species())
            {
                if (gs.name == std::string(species[si])) { known = true; break; }
            }
            if (dens[si] > 0 && known)
            {
                X[species[si]] = dens[si]/tot;
            }
        }
        if (X.empty()) return;

        bmech->setMixture(X);
        Boltzmann::SolverConfig cfg;
        cfg.T_gas_K = Tg;
        cfg.T_exc_K = Tg;
        cfg.growth  = Boltzmann::GrowthModel::Temporal;
        Boltzmann::BoltzmannSolver sol(bmech->db(), cfg);
        const auto r = sol.solve(en_Td);
        ++nSolves;

        // An unconverged EEDF is not merely inaccurate -- its normalisation is
        // wrong, so every rate taken from it is scaled by an unknown factor.
        // Counted rather than silently used.
        if (!r.converged) ++nUnconverged;

        // The POWER channels must come from the same EEDF as the rates. Left on
        // the tables they are evaluated at the sweep's reference composition and
        // temperature (here N2/O2/O = 0.774/0.186/0.04 at 1500 K) while the gas
        // is actually approaching 2600 K and 17% atomic O -- so the elastic and
        // vibrational shares, which depend on WHICH species the electrons hit,
        // would be taken from a state the gas left long ago. Refreshing the
        // rates but not the partition is not a consistent EEDF.
        // Stored as a RATIO to the table at the same E/N, not as an absolute
        // value. Held absolute, the channel is frozen between solves while E/N
        // keeps moving -- through the rising edge of a nanosecond pulse that is
        // a 4.5% staircase error in the deposited energy, because the tables
        // resolve E/N continuously and a piecewise-constant value does not.
        //
        // As a ratio it carries only what the table CANNOT know -- the change
        // in composition and gas temperature -- while the sharp E/N dependence
        // stays with the table where it is resolved properly. Both effects,
        // neither approximated.
        {
            const scalar enSI = en_Td*1e-21;
            const scalar mT = tableAt(tableDir/"muN_vs_reducedE", enSI);
            const scalar eT = tableAt(tableDir/"PelasticN_vs_reducedE", enSI);
            const scalar gT = tableAt(tableDir/"PgasN_vs_reducedE", enSI);
            const scalar vT = tableAt(tableDir/"PvibN_vs_reducedE", enSI);
            muNcorr = (mT > 0) ? r.transport.muN/mT : 1.0;
            PelCorr = (eT > 0) ? r.transport.PelasticN/eT : 1.0;
            PgsCorr = (gT > 0) ? r.transport.PgasN/gT : 1.0;
            PvbCorr = (vT > 0) ? r.transport.PvibN/vT : 1.0;
        }

        static bool dbg = Foam::getEnv("SOEEDF_DEBUG_EEDF").size();
        forAll(idToRate, i)
        {
            if (idToRate[i] >= 0
             && idToRate[i] < label(r.transport.rateCoeffs.size()))
            {
                const scalar kNew =
                    rateConv[i]*r.transport.rateCoeffs[idToRate[i]];
                // Same ratio argument as the power channels above: the table
                // keeps the E/N dependence, the solve supplies the composition
                // and temperature correction.
                const scalar kRaw = tableAt(
                    tableDir/("k_" + chem.tabulatedIds()[i] + "_vs_reducedE"),
                    en_Td*1e-21);
                kCorr[i] = (kRaw > 0) ? kNew/kRaw : 1.0;
                if (dbg && nSolves <= 3)
                {
                    Info<< "  [eedf " << nSolves << "] EN=" << en_Td
                        << " T=" << Tg << "  " << chem.tabulatedIds()[i]
                        << "  table=" << kRaw << "  solved=" << kNew
                        << "  ratio=" << kCorr[i]
                        << endl;
                }
            }
        }
    };

    // Initial state: dry air at the reference composition, lightly ionised.
    scalarField n(species.size(), Zero);
    const dictionary& comp = mech.subDict("composition");
    forAll(species, s)
    {
        if (comp.found(species[s]))
        {
            n[s] = readScalar(comp.lookup(species[s]))*nGas;
        }
    }
    if (args.found("X"))
    {
        n = Zero;
        std::string spec(args.get<string>("X"));
        for (auto& c : spec) if (c == '=' || c == ',') c = ' ';
        IStringStream is(spec);
        while (true)
        {
            word nm; scalar x;
            is >> nm;
            if (!is.good() && nm.empty()) break;
            is >> x;
            const label si = species.find(nm);
            if (si < 0)
            {
                FatalErrorInFunction << "-X names " << nm
                    << ", which the mechanism does not carry" << exit(FatalError);
            }
            n[si] = x*nGas;
            if (!is.good()) break;
        }
    }

    const label ie = species.find(electron);
    if (ie >= 0) n[ie] = ne0;
    // charge-neutral start: the electrons are balanced by N2+
    const label iN2p = species.find("N2p");
    if (iN2p >= 0) n[iN2p] = ne0;

    Info<< "plasmaChemistry0D: E/N = " << EN_Td << " Td, T = " << Tgas
        << " K, N = " << nGas << " 1/m3, n_e0 = " << ne0 << " 1/m3" << nl
        << "  charge residual of the RHS at t=0: "
        << chem.chargeResidual(n, kTab, Tgas) << endl;

    // ---- Jacobian verification -------------------------------------------
    //
    // A linearly-implicit method preserves a linear invariant EXACTLY only if
    // the Jacobian is exact: q.f = 0 for all states implies q.J = 0, and the
    // preservation follows from that. So a charge leak in the integrated
    // result is evidence about the JACOBIAN, not about the integrator -- which
    // is why this check exists rather than a tolerance study.
    {
        const label nEq = chem.nSpecie();
        scalarField y(n);
        scalarField f0(nEq), fp(nEq), dfdx(nEq);
        scalarSquareMatrix J(nEq, Zero);

        chem.derivatives(y, kTab, Tgas, f0);
        chem.jacobian(y, kTab, Tgas, dfdx, J);

        // Scaled comparison. A naive relative error is meaningless here:
        // d(N2)/dt is ~1e20 while perturbing a trace species changes it by
        // ~1e-3, so the finite difference is pure cancellation noise and
        // reports an "error" of 1 for a perfectly correct entry. Only entries
        // whose contribution rises above the cancellation floor of the column
        // are compared.
        scalar worst = 0; label wi = -1, wj = -1; label nCmp = 0;
        for (label j = 0; j < nEq; ++j)
        {
            const scalar h = 1e-6*max(mag(y[j]), scalar(1));
            scalarField yp(y); yp[j] += h;
            chem.derivatives(yp, kTab, Tgas, fp);

            // Cancellation floor: double precision on the largest |f| in this
            // column, times a safety margin.
            scalar fmax = 0;
            for (label i = 0; i < nEq; ++i) fmax = max(fmax, mag(f0[i]));
            const scalar floor = 1e4*SMALL*fmax/h;

            for (label i = 0; i < nEq; ++i)
            {
                const scalar fd = (fp[i] - f0[i])/h;
                const scalar sc = max(mag(fd), mag(J(i, j)));
                if (sc < floor) continue;          // unresolvable, not wrong
                ++nCmp;
                const scalar e = mag(fd - J(i, j))/sc;
                if (e > worst) { worst = e; wi = i; wj = j; }
            }
        }
        Info<< "  Jacobian: " << nCmp << " entries above the cancellation floor"
            << endl;
        Info<< "  Jacobian vs finite difference: worst relative error "
            << worst;
        if (wi >= 0)
        {
            Info<< "  at d(" << chem.species()[wi] << ")/d("
                << chem.species()[wj] << ")";
        }
        Info<< endl;
    }

    OFstream os(out);
    os << "t,EN_Td,Tgas,e_vib,E_dep,E_gas,E_vib";
    if (lmea) os << ",meanE_lmea,meanE_lfa,nEps";
    forAll(species, s) os << "," << species[s];
    os << nl;

    // ISOCHORIC, and that is the physics rather than a simplification. A
    // nanosecond pulse deposits its energy far faster than the gas can expand:
    // the acoustic time across a streamer channel of ~100 um is ~300 ns, so on
    // a 10 ns pulse the density is frozen and the PRESSURE rises. That pressure
    // rise is precisely what launches the blast wave, and it is the quantity
    // G3/G4 exist to propagate.
    //
    // So c_v, not c_p, DURING THE PULSE. Using c_p -- the natural habit, since
    // Cantera's is a constant-pressure reactor -- under-predicts the
    // temperature rise by 29% (dT scales as 1/c, and c_v/c_p = 1/gamma =
    // 0.714). Worse than the 29%, it loses the pressure rise ENTIRELY: at
    // constant pressure there is no pressure perturbation, and therefore no
    // blast wave to propagate.
    //
    // AFTER the pulse the opposite holds. Once the acoustic wave has crossed
    // the kernel -- hundreds of ns -- the pressure equalises with the
    // surroundings and the gas expands, which is a c_p process AND a genuine
    // energy loss from the kernel. Shao et al., Appl. Energy Combust. Sci. 19
    // (2024) 100280 handle exactly this: constant volume during the pulse
    // (their eq. 2a), an isentropic expansion at pulse end (eq. 12), then
    // constant pressure between pulses (eq. 14). This reactor implements only
    // the first, so it is right within a single pulse and over-predicts the
    // temperature of anything longer -- which matters most for repetitive
    // pulsing, where the error accumulates pulse on pulse.

    const scalar EVJ = 1.602176634e-19;

    // Thermodynamics from the mechanism, not from a hard-coded air fit. See
    // janafMixture. The fallback keeps older mechanisms running, but it is
    // air-specific, so it says so rather than quietly producing a temperature.
    const janafMixture thermo(mech, species, electron);
    if (!thermo.valid() && heating)
    {
        WarningInFunction
            << "this mechanism carries no `speciesThermo` block, so the gas" << nl
            << "    heat capacity falls back to a polynomial fitted to AIR." << nl
            << "    That is wrong for any other gas. Recompile with mechc to" << nl
            << "    emit thermo from the mechanism's own NASA7 data." << endl;
    }

    scalar T = Tgas;
    scalar eVib = 0.0;                        // J/m^3 in the reservoir
    scalar Edep = 0.0, Egas = 0.0, Evib = 0.0;
    scalar tauVTend = 0.0;

    // TWO-PHASE STEPPING, because a nanosecond pulse followed by a millisecond
    // afterglow spans six decades and a uniform step cannot serve both. At
    // dt = endTime/nOut a 1 ms run steps at 250 ns and NEVER SAMPLES a 10 ns
    // pulse: the run reported zero deposited energy and a temperature rise of
    // 0.01 K, which looks like a physics result and is a sampling failure.
    //
    // So: fixed fine steps through the pulse, then geometric growth. The
    // afterglow is smooth on a log axis -- V-T relaxation and recombination
    // are exponentials -- so geometric steps resolve it at a few dozen points.
    // A PROFILED run must resolve its own profile. With -profile the field is
    // not a Gaussian this code knows the width of, so the fine-step window runs
    // until the profile's field has decayed to 1% of its peak, and the step is
    // set from the profile's own sampling.
    //
    // Without this the step grew to 2 ns while the pulse was 10 ns wide. The
    // static-table run absorbed it; the dynamic-EEDF run, whose rates are
    // higher, went unstable and oscillated the electron density between zero
    // and 8.6e21 m^-3 -- which reads as a chemistry failure and is a
    // resolution failure.
    scalar tPulseEnd = pulsed ? (tc_ns + 5.0*fwhm_ns)*1e-9 : 0.0;
    scalar dtFine = pulsed ? fwhm_ns*1e-9/50.0 : endTime/nOut;
    if (profiled)
    {
        scalar enMax = 0;
        forAll(pEN, i) enMax = max(enMax, pEN[i]);
        tPulseEnd = pT[pT.size()-1];
        for (label i = pT.size() - 1; i >= 0; --i)
        {
            if (pEN[i] > 0.01*enMax) { tPulseEnd = pT[i]; break; }
        }
        scalar dtMin = GREAT;
        for (label i = 1; i < pT.size(); ++i)
        {
            if (pT[i] <= tPulseEnd) dtMin = min(dtMin, pT[i] - pT[i-1]);
        }
        dtFine = min(dtMin/10.0, tPulseEnd/200.0);
        Info<< "profile: fine step " << dtFine << " s until " << tPulseEnd
            << " s, then growing" << endl;
    }

    if (lmeaDtFix > 0)
    {
        dtFine = lmeaDtFix;
        Info<< "lmeaDt: step PINNED at " << dtFine
            << " s (pulse-derived step and afterglow growth both overridden)"
            << endl;
    }

    scalar t = 0.0, dt = min(dtFine, endTime/10.0);
    scalar tNextOut = 0.0;                    // next sample time, see below
    scalarField Pchem(chem.nSpecie(), 0.0), Lchem(chem.nSpecie(), 0.0);
    scalarField n0chem(chem.nSpecie(), 0.0);
    scalar wallLossPeak = 0.0;

    // LMEA state: the ENERGY DENSITY n_eps [eV/m^3], not the mean energy.
    // Conservative by choice -- see the model header. Seeded from the LFA
    // table so the run starts in equilibrium and any drift away from it is
    // the model's own doing rather than a transient from a cold start.
    scalar nEps = 0.0;
    if (lmea)
    {
        const scalar en0 = fieldAt(0.0);
        const scalar e0 = (en0 > 0)
            ? tableAt(tableDir/"meanEnergy_vs_reducedE", en0*1e-21) : 0.0;
        nEps = e0*((ie >= 0) ? max(n[ie], scalar(0)) : 0.0);
    }
    scalar meanELmea = 0.0, meanELfa = 0.0;

    // OPTION 4 set-up. Done ONCE, before any integration: enableEnergyEquation
    // changes nEqns() and rebuilds the ODE solver, whose workspace is sized at
    // construction.
    //
    // The lambdas are pure table reads with NO captured mutable state -- the
    // gas density they would otherwise need is passed per cell through
    // setEnergyCell(), because it changes with the gas temperature.
    if (lmeaOde)
    {
        if (ie < 0)
        {
            FatalErrorInFunction
                << "-lmeaSource ode needs the electron in the transported set."
                << exit(FatalError);
        }

        const fileName td = tableDir;
        chem.enableEnergyEquation
        (
            ie,
            1.0e-3,            // eps floor, as the hand-rolled path uses
            100.0,             // eps ceiling, as meanEnergyMax
            1.0,               // n_e floor inside the ratio
            [td](const scalar e) { return tableAt(td/"muN_vs_meanE", e); },
            [td](const scalar e)
            {
                return tableAt(td/"PelasticN_vs_meanE",   e)
                     + tableAt(td/"PinelasticN_vs_meanE", e);
            }
        );
    }

    label nWritten = 0;
    label k = 0;
    while (t < endTime)
    {
        if (pulsed && t > tPulseEnd)
        {
            // Grow, but CAPPED. Left to reach endTime/50 the step reaches
            // ~20 us in the afterglow, where recombination is still fast
            // enough that the stiff solver gives up -- it failed with species
            // driven to denormals (1e-323). The cap is a property of the
            // chemistry rather than of the output cadence, hence absolute.
            dt = min(min(dt*1.05, endTime/50.0), 1.0e-6);
            if (lmeaDtFix > 0) dt = lmeaDtFix;   // pinned: no growth at all
        }
        dt = min(dt, endTime - t);
        const scalar en = fieldAt(t);
        ++k;

        // ---- LMEA electron energy ------------------------------------------
        //
        //   d(n_eps)/dt = Joule - P_loss
        //   Joule  = mu_e n_e E^2   = (muN) n_e N (E/N)^2
        //   P_loss = n_e N (PelasticN + PinelasticN)(eps_bar)
        //
        // NOTE what is NOT here: a term for the energy carried by newly
        // created electrons. It is not missing -- n_eps is conservative and
        // n_e grows through the chemistry, so eps = n_eps/n_e carries that
        // term automatically via d(n_e eps) = eps dn_e + n_e deps. The 0-D
        // table gate measured exactly this: Joule - P_loss = eps nu_eff in
        // steady state, to 0.3%.
        if (lmea)
        {
            const scalar ne = (ie >= 0) ? max(n[ie], scalar(0)) : 0.0;
            const scalar Ngas = pres/(1.380649e-23*T);

            meanELfa = (en > 0)
                ? tableAt(tableDir/"meanEnergy_vs_reducedE", en*1e-21) : 0.0;

            // Floor then divide, never divide then clamp.
            meanELmea = nEps/max(ne, scalar(1.0));
            meanELmea = min(max(meanELmea, scalar(1e-3)), scalar(100.0));

            const scalar muN = tableAt(tableDir/"muN_vs_meanE", meanELmea);
            const scalar PlN =
                tableAt(tableDir/"PelasticN_vs_meanE",   meanELmea)
              + tableAt(tableDir/"PinelasticN_vs_meanE", meanELmea);

            const scalar enSI = en*1e-21;
            const scalar joule = muN*ne*Ngas*enSI*enSI;
            const scalar ploss = ne*Ngas*PlN;

            if (lmeaOde)
            {
                // NOTHING TO DO. The energy density is a component of the
                // integrated state vector and was advanced by the stiff
                // solver alongside the species -- see integrateWithEnergy().
                // Updating it here as well would apply the source twice.
                //
                // Note the reported meanE therefore lags by one step: this
                // block runs BEFORE the chemistry, so `nEps` here is the
                // start-of-step value. Harmless for a diagnostic, and stated
                // rather than left to be discovered.
            }
            else if (lmeaNewton)
            {
                // HAGELAAR'S LINEARISATION: implicit in the SENSITIVITY of the
                // source to the mean energy, not merely in the magnitude of
                // the loss.
                //
                //     S(eps),  eps = n_eps/n_e   =>   dS/d(n_eps) = S'(eps)/n_e
                //     D = -S'(eps)/n_e = N [ dPlossN/deps - (E/N)^2 dmuN/deps ]
                //
                // Both terms are POSITIVE for air -- PlossN rises with eps,
                // muN falls -- so D damps. That is a property of the gas, not
                // a theorem, which is why D is clamped at zero: a negative D
                // is anti-damping, the same shape that broke the Rosenbrock
                // controller and the transport cross-term sink.
                const scalar h = max(1e-3*meanELmea, 1e-4);

                const scalar dMu =
                    (tableAt(tableDir/"muN_vs_meanE", meanELmea + h)
                   - tableAt(tableDir/"muN_vs_meanE", max(meanELmea - h, 0.0)))
                  / (2.0*h);

                const scalar dPl =
                    ( tableAt(tableDir/"PelasticN_vs_meanE",   meanELmea + h)
                    + tableAt(tableDir/"PinelasticN_vs_meanE", meanELmea + h)
                    - tableAt(tableDir/"PelasticN_vs_meanE",
                              max(meanELmea - h, 0.0))
                    - tableAt(tableDir/"PinelasticN_vs_meanE",
                              max(meanELmea - h, 0.0)) ) / (2.0*h);

                const scalar D = max(Ngas*(dPl - enSI*enSI*dMu), 0.0);

                // Implicit in the increment: (1 + D dt) dn = S dt
                nEps = max(nEps + (joule - ploss)*dt/(1.0 + D*dt), 0.0);
            }
            else if (lmeaImplicit)
            {
                // n^{k+1} = (n + P dt)/(1 + L dt): unconditionally stable and
                // positivity-preserving for L >= 0, which is the whole reason
                // to write the loss as a RATE rather than as a source.
                const scalar L = (nEps > VSMALL) ? max(ploss/nEps, 0.0) : 0.0;
                nEps = (nEps + joule*dt)/(1.0 + L*dt);
            }
            else
            {
                nEps = max(nEps + (joule - ploss)*dt, 0.0);
            }
        }
        {
            // Rates follow the field, and the field follows the pulse. With
            // the field off there is no electron-impact chemistry at all --
            // the heavy reactions carry the afterglow on their own.
            // Table x correction: E/N from the table, where it is resolved;
            // composition and T_gas from the last Boltzmann solve. Refreshed
            // whenever any of the three can have moved.
            if (pulsed || heating || dynamicEEDF)
            {
                forAll(chem.tabulatedIds(), i)
                {
                    kTab[i] = (en > 0)
                        ? kCorr[i]*tableAt(
                              tableDir/("k_" + chem.tabulatedIds()[i]
                                        + "_vs_reducedE"), en*1e-21)
                        : 0.0;
                }
            }

            // Re-solve the EEDF when the state has drifted past tolerance.
            if (dynamicEEDF)
            {
                static scalarField Xlast; static scalar Tlast = -1, ENlast = -1;
                scalar tot = 0;
                forAll(species, si) if (n[si] > 0) tot += n[si];
                if (Xlast.size() != species.size())
                {
                    Xlast.setSize(species.size(), Zero); Tlast = -1;
                }
                bool due = (Tlast < 0)
                        || (mag(T - Tlast) > eedfTol*Tlast)
                        || (ENlast > 0 && mag(en - ENlast) > eedfTol*ENlast);
                if (!due && tot > 0)
                {
                    forAll(species, si)
                    {
                        const scalar x = n[si]/tot;
                        // Absolute floor as well as relative: a trace species
                        // doubling from 1e-12 is not a change in the mixture.
                        if (mag(x - Xlast[si]) > eedfTol*max(Xlast[si], 1e-3))
                        {
                            due = true; break;
                        }
                    }
                }
                if (due && en > 0 && tot > 0)
                {
                    solveEEDF(en, T, n);
                    forAll(species, si) Xlast[si] = n[si]/tot;
                    Tlast = T; ENlast = en;
                }
            }

            // The measured electron density is IMPOSED, not integrated. In 0-D
            // there is no space charge to stop ionisation running away, so the
            // experiment's own n_e(t) is the only honest driver -- the same
            // choice Cheng et al. (Combust. Flame 240 (2022) 111990) and the
            // simulations they compare against make.
            if (profiled && ie >= 0) n[ie] = interpAt(pNE, t);

            // Floor the state before integrating. Species pushed to denormal
            // values (~1e-323) destabilise the stiff solver without carrying
            // any physics: below one particle per cubic metre there is nothing
            // there. Same floor the CFD source term uses.
            forAll(n, si) if (n[si] < 1.0) n[si] = 0.0;

            // CHEMISTRY SOURCE, mirroring plasmaTransport's `chemistrySource`.
            //
            //   ode           stiff substep over dt (rodas23). Robust, and the
            //                 default, because plasma chemistry is stiff.
            //   implicitRate  one semi-implicit Euler step on the instantaneous
            //                 linearisation dn/dt = P - L n, which is exactly
            //                 what fvm::Sp(L, n) does in the CFD:
            //                     n <- (n + P dt)/(1 + L dt)
            //   adaptive      implicitRate while the linearisation resolves the
            //                 step, the stiff substep when it does not.
            //
            // The adaptive CRITERION is an analogue, not the CFD's. There the
            // test is the Picard change between outer correctors -- a quantity
            // that does not exist without outer correctors. Here it is the
            // stiffness the linearisation actually sees, max(L) dt: below the
            // threshold one implicit step resolves the fastest loss, above it
            // the step is being asked to leap over that timescale.
            if (chemSource == csODE)
            {
                integrateWithEnergy
                (
                    chem, n, nEps, kTab, T, dt,
                    lmeaOde, en*1e-21*(pres/(1.380649e-23*T)), pres/(1.380649e-23*T)
                );
                ++nODEsteps;
            }
            else
            {
                chem.productionLoss(n, kTab, T, Pchem, Lchem);
                scalar Ldtmax = 0.0;
                forAll(Lchem, si) Ldtmax = max(Ldtmax, Lchem[si]*dt);

                // The L*dt screen is a STABILITY floor and applies to
                // BOTH adaptive variants. `adaptiveError` adds its error
                // test on top rather than replacing this: measured
                // without it, the error estimate admitted steps 25-359%
                // wrong, because Richardson compares two estimates that
                // share the same instability and therefore agree while
                // both are wrong.
                if ((chemSource == csAdaptive
                  || chemSource == csAdaptiveError)
                 && Ldtmax > stiffTol)
                {
                    integrateWithEnergy
                    (
                        chem, n, nEps, kTab, T, dt,
                        lmeaOde, en*1e-21*(pres/(1.380649e-23*T)), pres/(1.380649e-23*T)
                    );
                    ++nODEsteps;
                }
                else
                {
                    // PICARD-ITERATED, which is what makes it implicit.
                    //
                    // Only the LOSS is implicit in one pass; the production
                    // term is evaluated at the old iterate, so a pair of
                    // species that feed each other is coupled explicitly and
                    // the step can diverge -- measured at 1e65 on an afterglow
                    // before the outer loop was added. The CFD gets its outer
                    // iterations from PIMPLE's correctors and warns when they
                    // are set to 1; without them the scheme is merely first
                    // order there and outright unstable here, where there is
                    // no transport to damp it.
                    //
                    // Always integrated FROM the start-of-step state, never
                    // from the partially-updated iterate, so repeating the
                    // sweep is idempotent rather than cumulative. Same reason
                    // the CFD keeps chemN0_.
                    n0chem = n;
                    for (label it = 0; it < nOuterCorr; ++it)
                    {
                        if (it > 0) chem.productionLoss(n, kTab, T, Pchem, Lchem);
                        forAll(n, si)
                        {
                            n[si] = (n0chem[si] + Pchem[si]*dt)
                                  / (1.0 + Lchem[si]*dt);
                        }
                    }

                    // VERIFY THE STEP, do not merely predict it.
                    //
                    // max(L) dt screens the DIAGONAL only, and the instability
                    // of this scheme lives off the diagonal: a pair of species
                    // that feed each other -- attachment and detachment, here
                    // -- each sees the other at the old iterate, and the pair
                    // can grow without bound while every individual L dt stays
                    // small. Screened on L alone, an afterglow reached 1e85 m^-3
                    // with the criterion reporting the step as safe.
                    //
                    // So the cheap step is TAKEN and then checked, which costs
                    // nothing extra when it succeeds and is exact when it does
                    // not. A step that moves a significant species by more than
                    // changeTol, or produces anything non-finite, is discarded
                    // and redone with the stiff substep from the same state.
                    bool bad = false;
                    forAll(n, si)
                    {
                        if (!std::isfinite(n[si])) { bad = true; break; }
                        // Relative-change test: csAdaptive only. Its
                        // normalisation (each species against itself) is the
                        // one csAdaptiveError replaces.
                        if (chemSource != csAdaptive) continue;
                        const scalar ref = max(n0chem[si], nFloorRel);
                        if (n0chem[si] > nFloorRel
                         && mag(n[si] - n0chem[si]) > changeTol*ref)
                        {
                            bad = true; break;
                        }
                    }

                    // adaptiveError: replace the relative-change verdict with a
                    // Richardson estimate of the LOCAL ERROR. `n` already holds
                    // the full step from n0chem; redo the same Picard update as
                    // two half-steps and compare. Non-finite still rejects
                    // outright -- an error estimate on a NaN means nothing.
                    if (chemSource == csAdaptiveError && !bad)
                    {
                        scalarField nHalf(n0chem);
                        scalarField Ph(Pchem.size()), Lh(Lchem.size());

                        for (label half = 0; half < 2; ++half)
                        {
                            const scalarField nStart(nHalf);
                            for (label it = 0; it < nOuterCorr; ++it)
                            {
                                chem.productionLoss(nHalf, kTab, T, Ph, Lh);
                                forAll(nHalf, si)
                                {
                                    nHalf[si] = (nStart[si] + Ph[si]*0.5*dt)
                                              / (1.0 + Lh[si]*0.5*dt);
                                }
                            }
                        }

                        const scalar absPart =
                            errAbsFrac*max(n0chem[ie], scalar(0));

                        scalar errNorm = 0;
                        forAll(n, si)
                        {
                            if (!std::isfinite(nHalf[si])) { errNorm = GREAT; break; }
                            const scalar scale = absPart + errRelTol*n0chem[si];
                            errNorm = max
                            (
                                errNorm,
                                mag(nHalf[si] - n[si])/max(scale, VSMALL)
                            );
                        }
                        bad = (errNorm > 1.0);
                    }

                    if (bad
                     && (chemSource == csAdaptive
                      || chemSource == csAdaptiveError))
                    {
                        n = n0chem;
                        integrateWithEnergy
                        (
                            chem, n, nEps, kTab, T, dt,
                            lmeaOde, en*1e-21*(pres/(1.380649e-23*T)), pres/(1.380649e-23*T)
                        );
                        ++nODEsteps;
                        ++nRejected;
                    }
                    else
                    {
                        ++nRateSteps;
                    }
                }
            }

            // WALL LOSS, applied to every CHARGED species at the same
            // fractional rate.
            //
            // That uniformity is the point: a quasineutral plasma losing the
            // same FRACTION of every charged species loses exactly zero net
            // charge, so the charge residual stays machine-zero and the models
            // cannot smuggle in a space-charge error. The alternative -- a
            // per-species nu from each ion's own mobility -- is more faithful
            // to a multi-ion plasma and does NOT conserve charge on its own,
            // which is why Alves et al. cluster the ions into single + and -
            // components (their eq (6)) before applying it.
            //
            // Exponential rather than n -= nu n dt, so a step longer than
            // 1/nu decays instead of going negative.
            if (wlModel != wallLoss::wlNone)
            {
                scalar nePlus = 0.0, neMinus = 0.0, neE = 0.0;
                forAll(n, si)
                {
                    if (si == ie) { neE = max(n[si], scalar(0)); continue; }
                    if (charge[si] > 0) nePlus  += max(n[si], scalar(0));
                    if (charge[si] < 0) neMinus += max(n[si], scalar(0));
                }
                // T_e = (2/3)<eps>, their footnote to eq (2e). Taken from the
                // same sweep the rates come from, so it is the temperature of
                // the EEDF actually driving this chemistry.
                const scalar meanE = (en > 0)
                    ? tableAt(tableDir/"meanEnergy_vs_reducedE", en*1e-21)
                    : 0.03;
                wlState.Te = max(meanE*(2.0/3.0), scalar(0.02));
                wlState.Tg = T;
                wlState.N  = nGas;
                wlState.alpha  = (neE > 0) ? neMinus/neE : 0.0;
                wlState.nPlus  = nePlus;
                wlState.nMinus = neMinus;

                const scalar nuW = wallLoss::nuTransport(wlModel, wlState);
                if (nuW > 0)
                {
                    const scalar f = Foam::exp(-nuW*dt);
                    forAll(n, si) if (charge[si] != 0) n[si] *= f;
                    wallLossPeak = max(wallLossPeak, nuW);
                }
            }

            if (heating)
            {
                const scalar ne = (ie >= 0) ? max(n[ie], scalar(0)) : 0.0;
                // Exact: sum of n_i m_i over the mechanism's own species. Constant
                // in an isochoric reactor by mass conservation, but computed
                // each step so it stays right when that assumption is relaxed.
                const scalar rho = thermo.valid() ? thermo.rho(n)
                                                  : nGas*(28.96e-3/6.02214076e23);

                // LIVE heavy-particle density, not the initial one. The power
                // channels are tabulated per electron per unit GAS density, so
                // the multiplier has to be the density now -- and dissociation
                // raises it, by 7.5% over the Rusterholtz case as O2 -> 2O.
                // Frozen, the field power and every loss channel are under-
                // counted by that much at late time. In the isobaric reactor it
                // moves far more, because the gas expands as it heats.
                scalar nHeavy = 0.0;
                forAll(n, si)
                {
                    if (si != ie) nHeavy += max(n[si], scalar(0));
                }
                if (nHeavy <= 0.0) nHeavy = nGas;

                // Per electron per unit gas density. From the live EEDF when
                // one is being solved, otherwise from the sweep. Mixing the two
                // is what has to be avoided: the rates and the power channels
                // are moments of the SAME f0, and taking them from different
                // states breaks the energy budget they are supposed to close.
                const scalar Pel = (en <= 0) ? 0.0 : PelCorr*
                    tableAt(tableDir/"PelasticN_vs_reducedE", en*1e-21);
                const scalar Pgs = (en <= 0) ? 0.0 : PgsCorr*
                    tableAt(tableDir/"PgasN_vs_reducedE",     en*1e-21);
                const scalar Pvb = (en <= 0) ? 0.0 : PvbCorr*
                    tableAt(tableDir/"PvibN_vs_reducedE",     en*1e-21);
                const scalar muN = (en <= 0) ? 0.0 : muNcorr*
                    tableAt(tableDir/"muN_vs_reducedE",       en*1e-21);

                const scalar EN_SI2 = en*1e-21;
                const scalar Pdep = muN*EN_SI2*EN_SI2*ne*nHeavy*EVJ;   // W/m^3

                // Prompt heat: elastic/rotational plus the gas share of the
                // inelastic defect. The heavy reactions add fast gas heating
                // on top, from their own enthalpies.
                const scalar Qprompt = (Pel + Pgs)*ne*nHeavy*EVJ;
                const scalar Qheavy  = chem.heavyHeatRelease(n, T)*EVJ;

                // Vibrational reservoir. It FILLS during the pulse and empties
                // on tau_VT, which is microseconds -- so on this timescale the
                // point of tracking it is that the energy is withheld from the
                // gas, not that it comes back.
                const scalar Pvib = Pvb*ne*nHeavy*EVJ;

                // e_vib here is the NON-EQUILIBRIUM (excess) vibrational
                // energy, so it relaxes towards zero rather than towards
                // e_vib^eq(T). Tracking the excess is what lets the sink be a
                // plain e/tau: with the TOTAL energy as the variable the
                // equilibrium term is mandatory, and omitting it drives the
                // reservoir negative as soon as the gas is warm.
                const label iO = species.find("O");
                const scalar nO = (iO >= 0) ? max(n[iO], scalar(0)) : 0.0;
                const label iN2 = species.find("N2");
                const label iO2 = species.find("O2");
                const scalar ntot = max(sum(n), scalar(1));
                const scalar xN2 = (iN2 >= 0) ? max(n[iN2], scalar(0))/ntot : 0.0;
                const scalar xO2 = (iO2 >= 0) ? max(n[iO2], scalar(0))/ntot : 0.0;

                const scalar tauVT = (tauVTfixed > 0)
                    ? tauVTfixed
                    : tauVT_N2(T, pres/101325.0, xN2, xO2, nO);

                const scalar Qvt  = eVib/tauVT;
                tauVTend = tauVT;
                eVib += (Pvib - Qvt)*dt;

                // ISOCHORIC heats against c_v, ISOBARIC against c_p. The gas
                // cannot expand within a nanosecond pulse -- the acoustic time
                // across a 450 um channel is ~0.3 us -- so c_v is right there,
                // and it is what drives the pressure rise and hence the blast
                // wave. Past the acoustic time the gas has expanded, the work
                // has been done, and c_p is right. Neither is right for both,
                // which is why this is a mode and not a default.
                const scalar Qgas = Qprompt + Qheavy + Qvt;
                const scalar cGas = thermo.valid()
                    ? (isobaric ? thermo.cp(n, T) : thermo.cv(n, T))
                    : cvAirFallback(T)*(isobaric ? 1.4 : 1.0);
                T += Qgas*dt/(rho*cGas);

                if (isobaric)
                {
                    // Hold p by letting the gas EXPAND: every density is scaled
                    // so that n_heavy k T = p again. This absorbs both effects
                    // at once -- thermal expansion, and the extra moles the
                    // chemistry makes when it dissociates something.
                    //
                    // The electron partial pressure is neglected: at n_e/N ~
                    // 1e-3 and T_e ~ 15 T_gas it is ~1% of the total, well
                    // inside everything else here, and including it would mean
                    // committing to a T_e the reactor does not track.
                    scalar nh = 0.0;
                    forAll(n, si) if (si != ie) nh += max(n[si], scalar(0));
                    if (nh > 0.0)
                    {
                        const scalar f = (pres/(1.380649e-23*T))/nh;
                        forAll(n, si) n[si] *= f;
                    }
                }

                Edep += Pdep*dt;
                Egas += Qgas*dt;
                Evib += Pvib*dt;
            }
        }

        // THROTTLED to ~nOut rows. -nOut used to set only the step size, so a
        // profiled run -- which picks its own step from the profile -- wrote a
        // row per step: 357,143 rows and 80 MB for a 100 ns case, from a flag
        // that said 2000. Sampled on TIME rather than step count so the rows
        // stay evenly spaced through the two-phase stepping, and the final
        // state is always written whatever the sampling lands on.
        if (t >= tNextOut || t + dt >= endTime)
        {
            os << t << ',' << en << ',' << T << ',' << eVib
               << ',' << Edep << ',' << Egas << ',' << Evib;
            if (lmea)
            {
                os << ',' << meanELmea << ',' << meanELfa << ',' << nEps;
            }
            forAll(n, s) os << "," << n[s];
            os << nl;
            ++nWritten;
            tNextOut = t + endTime/scalar(max(nOut, label(1)));
        }

        t += dt;
    }
    Info<< "  steps taken: " << k << ", rows written: " << nWritten << endl;
    Info<< "  chemistry: " << csName << " (backend " << chem.backend() << ")"
        << ", stiff substeps " << nODEsteps
        << ", implicit-rate steps " << nRateSteps
        << (nRejected > 0
            ? " (" + Foam::name(nRejected) + " rejected -> stiff)" : "")
        << endl;

    // ALWAYS report which EEDF path ran. Omitting -manifest silently falls
    // back to the frozen tables and the run still completes normally, so the
    // only evidence was the absence of a line -- and absence is not something
    // a reader notices. A whole validation run was lost to exactly this.
    if (!dynamicEEDF)
    {
        Info<< "EEDF: FROZEN tables from " << tableDir << nl
            << "      (no -manifest, so the Boltzmann equation is NOT re-solved."
            << " Rates are" << nl
            << "       taken at the composition and T_gas the sweep was"
            << " generated at.)" << endl;
    }
    else
    {
        Info<< "dynamic EEDF: " << nSolves << " Boltzmann solves over the run, "
            << nUnconverged << " unconverged" << endl;
        if (nUnconverged > 0)
        {
            WarningInFunction
                << nUnconverged << " of " << nSolves << " EEDF solves did not"
                << " converge. Their normalisation is wrong, so every rate"
                << " taken from them is scaled by an unknown factor." << endl;
        }
    }

    if (wlModel != wallLoss::wlNone)
    {
        Info<< "wall loss: " << wlName << ", R = " << wlState.R
            << " m, L = " << wlState.L << " m, Lambda = "
            << wallLoss::diffusionLength(wlState.R, wlState.L) << " m" << nl
            << "    peak nu_transp = " << wallLossPeak << " 1/s"
            << "  (tau = " << (wallLossPeak > 0 ? 1.0/wallLossPeak : 0.0)
            << " s)" << endl;
    }

    if (heating)
    {
        // A SEEDED ion is not free. Starting with n_e = n_i means starting
        // with n_i x 15.6 eV of ionisation energy already in the box, put
        // there by the initial condition rather than by the field. It comes
        // back out through recombination and is indistinguishable, in the
        // output, from gas heating that the discharge actually paid for.
        //
        // At a 3e18 seed that is 7.5 J/m^3 -- which was a hundred times the
        // deposited energy in the first 100 Torr case run here, and produced a
        // "gas heating fraction" of 10000%. The fraction is meaningless
        // whenever the seed is comparable to the deposition, so say so rather
        // than print it.
        const scalar Eseed = ne0*15.6*EVJ;      // ~ionisation energy of the seed

        Info<< nl << "gas heating:" << nl
            << "  seed energy    " << Eseed << " J/m^3"
            << "   (ionisation energy of the initial n_e, NOT deposited)" << nl
            << "  deposited      " << Edep << " J/m^3" << nl
            << "  to the gas     " << Egas << " J/m^3"
            << "   (" << (Edep > 0 ? 100.0*Egas/Edep : 0.0) << " %)" << nl
            << "  to vibration   " << Evib << " J/m^3"
            << "   (" << (Edep > 0 ? 100.0*Evib/Edep : 0.0) << " %)" << nl
            << "  dT             " << T - Tgas << " K" << nl
            << "  tau_VT(end)    " << tauVTend << " s" << nl;

        if (isobaric)
        {
            Info<< "  reactor        ISOBARIC (p held, gas expands, heats"
                   " against c_p)" << nl
                << "  V/V0           " << T/Tgas
                << "   (thermal expansion; the extra moles from dissociation"
                   " add to this)" << nl
                << "  p/p0           1   (by construction -- so this mode"
                   " produces NO blast wave)" << endl;
        }
        else
        {
            Info<< "  reactor        ISOCHORIC (V held, heats against c_v)" << nl
                << "  p/p0           " << T/Tgas
                << "   (the density is frozen on this timescale, so the"
                   " pressure rises in proportion to T -- this is what drives"
                   " the blast wave)" << endl;
        }

        if (Eseed > 0.05*Edep)
        {
            WarningInFunction
                << "the seed carries " << Eseed << " J/m^3 of ionisation"
                << " energy against " << Edep << " J/m^3 deposited." << nl
                << "    Recombination returns the seed energy to the gas, so"
                << " the heating FRACTIONS above are not attributable to the"
                << " discharge. Lower -ne0 until the seed is negligible, or"
                << " read the absolute energies only." << endl;
        }
    }

    Info<< "  charge residual of the RHS at t=end: "
        << chem.chargeResidual(n, kTab, Tgas) << nl
        << "wrote " << out << endl;
    return 0;
}
