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
    argList args(argc, argv, false, false, false);

    const fileName mechFile = args.getOrDefault<fileName>
        ("mechanism", "constant/air_plasma.foam");
    const fileName tableDir = args.getOrDefault<fileName>
        ("tables", "constant/plasmaTables");
    const scalar EN_Td   = args.getOrDefault<scalar>("EN", 150.0);
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

    OFstream os(out);
    os << "t";
    forAll(species, s) os << "," << species[s];
    os << nl;

    const scalar dt = endTime/nOut;
    for (label k = 0; k <= nOut; ++k)
    {
        if (k > 0) chem.integrate(n, kTab, Tgas, dt);
        os << k*dt;
        forAll(n, s) os << "," << n[s];
        os << nl;
    }

    Info<< "  charge residual of the RHS at t=end: "
        << chem.chargeResidual(n, kTab, Tgas) << nl
        << "wrote " << out << endl;
    return 0;
}
