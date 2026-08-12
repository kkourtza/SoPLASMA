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
#include <cstdlib>

using namespace Foam;


// ---------------------------------------------------------------------------
// Vibrational-translational relaxation time of N2(v) in air [s].
//
// Landau-Teller relaxation needs a timescale, and a single constant is wrong by
// DECADES in a discharge: V-T relaxation of N2 by atomic OXYGEN is ~1e3 times
// faster than by N2 or O2, so tau collapses as the discharge dissociates O2.
// That is the whole reason this is a function of composition and not a number.
//
// Millikan & White, J. Chem. Phys. 39 (1963) 3209, for the molecular partners:
//
//     p0 tau_{k,m} = exp[ a (T^-1/3 - b) - 18.42 ]      [atm s]
//
// with (a,b) = (220, 0.03) for N2-N2 and (162, 0.03) for N2-O2 [Shao et al.,
// Appl. Energy Combust. Sci. 19 (2024) 100280, from Colgan & Levitt 1967].
//
// The N2-O channel is taken from POPOV's rate constant directly,
//
//     k = 4.5e-21 (T/300)^2.1  m^3/s     [Popov, J. Phys. D 44 (2011) 285201]
//     tau_{N2,O} = 1/(k n_O)
//
// rather than from the fitted form 488.5/(p0 T^1.1) quoted alongside it in
// Shao et al. eq. (7). Those two disagree by a factor of 101325 -- exactly Pa
// per atm -- so that constant requires p0 in PASCALS while the surrounding text
// says atm. Implemented as written with atm, V-T relaxation comes out 1e5 times
// too slow and the reservoir never empties. Deriving from the rate constant
// avoids the trap entirely, which is the general lesson: prefer the underlying
// quantity to a fitted restatement of it.
//
// Mixing rule, Millikan-White extended to several partners:
//     1/tau_k = SUM_m X_m / tau_{k,m}
// Constant-volume heat capacity of air [J/kg/K], cubic in T over 300-3500 K.
//
// NOT the 300 K value. c_v rises from 723 to 996 J/kg/K between 300 and 2400 K
// as the vibrational modes of N2 and O2 become active, and a discharge that
// heats air to 2000 K spends most of its time in the range where the constant
// is 30-40% wrong. Rusterholtz et al. (J. Phys. D 46 (2013) 464010) implicitly
// use ~1038 J/kg/K when converting their measured 900 K rise into 140 uJ.
//
// Fitted here against Cantera's own NASA polynomials for 79/21 N2/O2, max
// error 18 J/kg/K over the range. Composition dependence is neglected: a
// discharge that dissociates half the O2 changes c_v by a few percent, which
// is far below the uncertainty in the energy partitions feeding it.
static Foam::scalar cvAir(const Foam::scalar T)
{
    const scalar t = min(max(T, scalar(200)), scalar(3500));
    return 6.165129e+02 + 3.212567e-01*t
         - 8.907555e-05*t*t + 8.870573e-09*t*t*t;
}


static Foam::scalar tauVT_N2
(
    const Foam::scalar T,
    const Foam::scalar p_atm,
    const Foam::scalar xN2,
    const Foam::scalar xO2,
    const Foam::scalar nO           // atomic oxygen number density [m^-3]
)
{
    const scalar Tm13 = Foam::pow(T, -1.0/3.0);
    auto mw = [&](const scalar a, const scalar b)
    {
        return Foam::exp(a*(Tm13 - b) - 18.42)/max(p_atm, SMALL);   // [s]
    };

    scalar inv = 0.0;
    if (xN2 > 0) inv += xN2/mw(220.0, 0.03);       // N2 - N2
    if (xO2 > 0) inv += xO2/mw(162.0, 0.03);       // N2 - O2

    // N2 - O, from the rate constant. Written as a frequency because that is
    // what it is: k n_O, with no mole fraction needed.
    if (nO > 0)
    {
        inv += 4.5e-21*Foam::pow(T/300.0, 2.1)*nO;
    }

    return (inv > 0) ? 1.0/inv : GREAT;
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
    argList::addOption("tauVT", "s",
        "vibrational-translational relaxation time (default 1e-5 s)");
    argList args(argc, argv, false, false, false);

    const fileName mechFile = args.getOrDefault<fileName>
        ("mechanism", "constant/air_plasma.foam");
    const fileName tableDir = args.getOrDefault<fileName>
        ("tables", "constant/plasmaTables");
    const scalar EN_Td   = args.getOrDefault<scalar>("EN", 150.0);
    const bool   heating = args.found("gasHeating");

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
    label nSolves = 0, nUnconverged = 0;

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

        static bool dbg = Foam::getEnv("SOEEDF_DEBUG_EEDF").size();
        forAll(idToRate, i)
        {
            if (idToRate[i] >= 0
             && idToRate[i] < label(r.transport.rateCoeffs.size()))
            {
                const scalar kNew =
                    rateConv[i]*r.transport.rateCoeffs[idToRate[i]];
                if (dbg && nSolves <= 3)
                {
                    Info<< "  [eedf " << nSolves << "] EN=" << en_Td
                        << " T=" << Tg << "  " << chem.tabulatedIds()[i]
                        << "  table=" << kTab[i] << "  solved=" << kNew
                        << "  ratio=" << (kTab[i] > 0 ? kNew/kTab[i] : -1)
                        << endl;
                }
                kTab[i] = kNew;
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

    const scalar Mair = 28.96e-3/6.02214076e23;   // kg per particle
    const scalar EVJ = 1.602176634e-19;

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

    scalar t = 0.0, dt = min(dtFine, endTime/10.0);
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
        }
        dt = min(dt, endTime - t);
        const scalar en = fieldAt(t);
        ++k;
        {
            // Rates follow the field, and the field follows the pulse. With
            // the field off there is no electron-impact chemistry at all --
            // the heavy reactions carry the afterglow on their own.
            if (pulsed || heating)
            {
                forAll(chem.tabulatedIds(), i)
                {
                    kTab[i] = (en > 0)
                        ? tableAt(tableDir/("k_" + chem.tabulatedIds()[i]
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

            chem.integrate(n, kTab, T, dt);

            if (heating)
            {
                const scalar ne = (ie >= 0) ? max(n[ie], scalar(0)) : 0.0;
                const scalar rho = nGas*Mair;

                // Per electron per unit gas density, straight from the sweep.
                const scalar Pel = (en > 0) ? tableAt(tableDir/"PelasticN_vs_reducedE", en*1e-21) : 0.0;
                const scalar Pgs = (en > 0) ? tableAt(tableDir/"PgasN_vs_reducedE",     en*1e-21) : 0.0;
                const scalar Pvb = (en > 0) ? tableAt(tableDir/"PvibN_vs_reducedE",     en*1e-21) : 0.0;
                const scalar muN = (en > 0) ? tableAt(tableDir/"muN_vs_reducedE",       en*1e-21) : 0.0;

                const scalar EN_SI2 = en*1e-21;
                const scalar Pdep = muN*EN_SI2*EN_SI2*ne*nGas*EVJ;   // W/m^3

                // Prompt heat: elastic/rotational plus the gas share of the
                // inelastic defect. The heavy reactions add fast gas heating
                // on top, from their own enthalpies.
                const scalar Qprompt = (Pel + Pgs)*ne*nGas*EVJ;
                const scalar Qheavy  = chem.heavyHeatRelease(n, T)*EVJ;

                // Vibrational reservoir. It FILLS during the pulse and empties
                // on tau_VT, which is microseconds -- so on this timescale the
                // point of tracking it is that the energy is withheld from the
                // gas, not that it comes back.
                const scalar Pvib = Pvb*ne*nGas*EVJ;

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

                const scalar Qgas = Qprompt + Qheavy + Qvt;
                T += Qgas*dt/(rho*cvAir(T));

                Edep += Pdep*dt;
                Egas += Qgas*dt;
                Evib += Pvib*dt;
            }
        }

        os << t << ',' << en << ',' << T << ',' << eVib
           << ',' << Edep << ',' << Egas << ',' << Evib;
        forAll(n, s) os << "," << n[s];
        os << nl;

        t += dt;
    }
    Info<< "  steps taken: " << k << endl;

    if (dynamicEEDF)
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
            << "  tau_VT(end)    " << tauVTend << " s" << nl
            << "  p/p0           " << T/Tgas
            << "   (isochoric: the density is frozen on this timescale, so the"
               " pressure rises in proportion to T -- this is what drives the"
               " blast wave)" << endl;

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
