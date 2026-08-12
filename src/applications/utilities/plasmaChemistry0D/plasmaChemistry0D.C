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
#include "DynamicList.H"
#include <cstdlib>

using namespace Foam;

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
    const scalar tauVT = args.getOrDefault<scalar>("tauVT", 1.0e-5);

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
    auto fieldAt = [&](const scalar t) -> scalar
    {
        if (!pulsed) return EN_Td;
        const scalar sig = fwhm_ns*1e-9/(2.0*Foam::sqrt(2.0*Foam::log(2.0)));
        const scalar x = (t - tc_ns*1e-9)/sig;
        return max(pkEN*Foam::exp(-0.5*x*x), scalar(0.1));
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
    // So c_v, not c_p. Using c_p here -- the obvious mistake, since "constant
    // pressure reactor" is the Cantera habit -- would under-predict the
    // temperature rise by 40% and the pressure rise with it.
    const scalar cv = 718.0;                  // J/kg/K, air at 300 K
    const scalar Mair = 28.96e-3/6.02214076e23;   // kg per particle
    const scalar EVJ = 1.602176634e-19;

    scalar T = Tgas;
    scalar eVib = 0.0;                        // J/m^3 in the reservoir
    scalar Edep = 0.0, Egas = 0.0, Evib = 0.0;

    const scalar dt = endTime/nOut;
    for (label k = 0; k <= nOut; ++k)
    {
        const scalar t = k*dt;
        const scalar en = fieldAt(t);

        if (k > 0)
        {
            // Rates follow the field, and the field follows the pulse.
            if (pulsed || heating)
            {
                forAll(chem.tabulatedIds(), i)
                {
                    kTab[i] = tableAt(tableDir/("k_" + chem.tabulatedIds()[i]
                                                + "_vs_reducedE"), en*1e-21);
                }
            }

            chem.integrate(n, kTab, T, dt);

            if (heating)
            {
                const scalar ne = (ie >= 0) ? max(n[ie], scalar(0)) : 0.0;
                const scalar rho = nGas*Mair;

                // Per electron per unit gas density, straight from the sweep.
                const scalar Pel = tableAt(tableDir/"PelasticN_vs_reducedE", en*1e-21);
                const scalar Pgs = tableAt(tableDir/"PgasN_vs_reducedE",     en*1e-21);
                const scalar Pvb = tableAt(tableDir/"PvibN_vs_reducedE",     en*1e-21);
                const scalar muN = tableAt(tableDir/"muN_vs_reducedE",       en*1e-21);

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
                const scalar Qvt  = eVib/tauVT;
                eVib += (Pvib - Qvt)*dt;

                const scalar Qgas = Qprompt + Qheavy + Qvt;
                T += Qgas*dt/(rho*cv);

                Edep += Pdep*dt;
                Egas += Qgas*dt;
                Evib += Pvib*dt;
            }
        }

        os << t << ',' << en << ',' << T << ',' << eVib
           << ',' << Edep << ',' << Egas << ',' << Evib;
        forAll(n, s) os << "," << n[s];
        os << nl;
    }

    if (heating)
    {
        Info<< nl << "gas heating:" << nl
            << "  deposited      " << Edep << " J/m^3" << nl
            << "  to the gas     " << Egas << " J/m^3"
            << "   (" << (Edep > 0 ? 100.0*Egas/Edep : 0.0) << " %)" << nl
            << "  to vibration   " << Evib << " J/m^3"
            << "   (" << (Edep > 0 ? 100.0*Evib/Edep : 0.0) << " %)" << nl
            << "  dT             " << T - Tgas << " K" << nl
            << "  p/p0           " << T/Tgas
            << "   (isochoric: the density is frozen on this timescale, so the"
               " pressure rises in proportion to T -- this is what drives the"
               " blast wave)" << endl;
    }

    Info<< "  charge residual of the RHS at t=end: "
        << chem.chargeResidual(n, kTab, Tgas) << nl
        << "wrote " << out << endl;
    return 0;
}
