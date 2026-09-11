/*---------------------------------------------------------------------------*\
Application
    testVibRelax

Description
    Checks the V-T relaxation model against its own sources and against the
    physics it is supposed to represent.

    This closure was the least tested thing in the 0-D reactor and the one the
    documentation flagged as weakest: a single lumped vibrational reservoir
    drained on a composition-dependent tau. It sets how much of the deposited
    energy reaches the gas over microseconds, so getting it wrong moves the
    slow heating without touching the fast heating that is already validated.

    There is no experiment in this file. These are checks against Millikan &
    White's correlation, Popov's rate constant, and the limits both must obey
    -- which is what can be verified without reference data, and is exactly the
    part that would break silently on a units slip.

\*---------------------------------------------------------------------------*/

#include "vibRelax.H"

#include <cmath>
#include <cstdio>
#include <string>

using namespace Foam;
using namespace Foam::vibRelax;

static int nFail = 0;

static void check(const char* name, bool ok, const std::string& detail)
{
    std::printf("  [%s] %-50s (%s)\n", ok ? "ok  " : "FAIL", name,
                detail.c_str());
    if (!ok) ++nFail;
}

static std::string f2(const char* f, double a, double b = 0, double c = 0)
{
    char buf[256]; std::snprintf(buf, sizeof(buf), f, a, b, c); return buf;
}

static bool close(double a, double b, double tol)
{
    return std::fabs(a - b) <= tol*std::fabs(b);
}


int main()
{
    std::printf("testVibRelax: Millikan-White + Popov, and their limits\n\n");

    // ---- the correlation is implemented as written -----------------------
    std::printf("Millikan-White, eq. as published\n");
    {
        // Independent evaluation of exp[a(T^-1/3 - b) - 18.42]/p.
        auto ref = [](double T, double p, double a, double b)
        {
            return std::exp(a*(std::pow(T, -1.0/3.0) - b) - 18.42)/p;
        };
        bool ok = true;
        for (double T : {300.0, 1000.0, 2000.0, 3000.0})
        {
            if (!close(millikanWhite(T, 1.0, 220.0, 0.03), ref(T, 1.0, 220.0, 0.03), 1e-12))
                ok = false;
        }
        check("N2-N2 matches the published form", ok, "300-3000 K");

        // Pressure scaling is exactly 1/p: V-T is a two-body process.
        check("tau scales as 1/p",
              close(millikanWhite(1000.0, 2.0, 220.0, 0.03),
                    0.5*millikanWhite(1000.0, 1.0, 220.0, 0.03), 1e-12),
              "doubling p halves tau");

        // Landau-Teller: relaxation accelerates steeply with temperature.
        const double t300 = millikanWhite(300.0, 1.0, 220.0, 0.03);
        const double t3000 = millikanWhite(3000.0, 1.0, 220.0, 0.03);
        check("tau falls steeply with T (Landau-Teller)",
              t3000 < t300/1e6,
              f2("%.3e s at 300 K -> %.3e s at 3000 K", t300, t3000));

        // N2 at room temperature relaxes on SECONDS at 1 atm. This is the
        // well-known fact that makes nitrogen a vibrational energy reservoir
        // at all, and a model that got it merely "slow" would still be wrong.
        check("N2-N2 at 300 K is of order 1e3 s at 1 atm",
              t300 > 1e2 && t300 < 1e5, f2("%.3e s", t300));

        // O2 is the more efficient molecular partner at every temperature.
        bool o2Faster = true;
        for (double T : {300.0, 1000.0, 3000.0})
        {
            if (millikanWhite(T, 1.0, 162.0, 0.03)
              >= millikanWhite(T, 1.0, 220.0, 0.03)) o2Faster = false;
        }
        check("O2 relaxes N2 faster than N2 does", o2Faster, "300-3000 K");
    }

    // ---- the mixing rule --------------------------------------------------
    std::printf("\nmixing rule, 1/tau = SUM x_m/tau_m\n");
    {
        const double T = 1000.0, p = 1.0;
        // Pure N2 must reproduce the single-partner value exactly.
        check("x_N2 = 1 reproduces the N2-N2 time",
              close(tauVT_N2(T, p, 1.0, 0.0, 0.0),
                    millikanWhite(T, p, 220.0, 0.03), 1e-12),
              "pure N2");

        // A mixture must relax faster than either pure component alone --
        // frequencies add, so tau can only fall as partners are added.
        const double tMix = tauVT_N2(T, p, 0.8, 0.2, 0.0);
        check("a mixture relaxes faster than pure N2",
              tMix < tauVT_N2(T, p, 1.0, 0.0, 0.0),
              f2("%.3e s vs %.3e s", tMix, tauVT_N2(T, p, 1.0, 0.0, 0.0)));

        // Explicit sum, evaluated independently.
        const double inv = 0.8/millikanWhite(T, p, 220.0, 0.03)
                         + 0.2/millikanWhite(T, p, 162.0, 0.03);
        check("mixture equals the explicit frequency sum",
              close(tMix, 1.0/inv, 1e-12), "80/20 N2/O2");
    }

    // ---- the atomic-oxygen channel ---------------------------------------
    std::printf("\natomic oxygen, Popov's rate constant\n");
    {
        // tau = 1/(k n_O) with k = 4.5e-21 (T/300)^2.1 m^3/s.
        const double T = 2700.0, nO = 8.7e23;
        check("nu_O = k n_O as published",
              close(nuVT_N2_O(T, nO), 4.5e-21*std::pow(T/300.0, 2.1)*nO, 1e-12),
              f2("%.4e 1/s", nuVT_N2_O(T, nO)));

        check("no O means no O channel", nuVT_N2_O(T, 0.0) == 0.0, "");

        // THE UNITS TRAP. Shao et al. also quote tau = 488.5/(p0 T^1.1), which
        // requires p0 in PASCALS although the text says atm -- a factor of
        // 101325. Implemented from the rate constant instead, so the two must
        // agree only when that form is read with p0 in Pa.
        const double p_Pa = 101325.0;
        const double tauFit = 488.5/(p_Pa*std::pow(T, 1.1));
        const double tauOurs = 1.0/nuVT_N2_O(T, nO);
        check("agrees with the Pa-reading of the fitted form to O(1)",
              tauOurs/tauFit > 0.2 && tauOurs/tauFit < 5.0,
              f2("ours %.3e s, fit(Pa) %.3e s", tauOurs, tauFit));
        std::printf("        (the atm reading would give %.3e s -- 1e5 too "
                    "slow, and the reservoir would never empty)\n",
                    488.5/(1.0*std::pow(T, 1.1)));
    }

    // ---- the state that actually matters ---------------------------------
    std::printf("\nat the Rusterholtz 20 ns state (2700 K, 1 atm, 17%% O)\n");
    {
        const double T = 2700.0, nO = 8.7e23;
        const double tauMol = tauVT_N2(T, 1.0, 0.72, 0.11, 0.0);
        const double tauAll = tauVT_N2(T, 1.0, 0.72, 0.11, nO);

        // The reactor reports tau_VT(end) ~ 1.8-2.4 us on this benchmark, so
        // this is the number the slow heating actually rides on.
        check("tau is microseconds, not milliseconds",
              tauAll > 1e-7 && tauAll < 1e-5, f2("%.3e s", tauAll));

        check("atomic O dominates the relaxation here",
              tauAll < tauMol/10.0,
              f2("%.3e s -> %.3e s, %.0fx", tauMol, tauAll, tauMol/tauAll));

        // AND the crucial timescale separation: tau must be far longer than a
        // nanosecond pulse, or the reservoir would drain while it fills and
        // the fast/slow heating split would be meaningless.
        check("tau >> 10 ns pulse (the reservoir fills before it drains)",
              tauAll > 100e-9, f2("%.3e s vs 1e-8 s", tauAll));
    }

    std::printf("\n%s\n", nFail ? "FAILED" : "all checks passed");
    return nFail ? 1 : 0;
}
