/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.
    Copyright (C) 2026
    GNU General Public License v3 or later.

Application
    testCoulombHeating

Description
    Unit bed for coulombHeatingClosure::epsilonEff (eq. 14-15 of Eliseev,
    Bogdanov & Kudryavtsev, Phys. Plasmas 24, 093503 (2017)), exercised
    directly against analytic ground truth -- no mesh, no dictionary, no
    plasmaPropertyEvaluator machinery. Same pattern as testWallFlux.

    Checks, in order:
      1) nu_ee reproduces a HAND-VERIFIED reference value (2026-09-08,
         confirmed against the standard NRL Plasma Formulary to 1.7% and
         against the user's own reading of eq. (14) from the actual PDF).
      2) epsilon_eff -> 0 as n_e -> 0 (bounded below, no ionisation-degree
         gate needed -- see the closure header).
      3) epsilon_eff -> epsilon_s = secondaryEnergyFraction*epsilon_1 as
         n_e -> large (nu_ee dominates the denominator).
      4) epsilon_eff is monotonically increasing in n_e at fixed Te.
      5) A SMALLER diffusion length gives a SMALLER (never larger)
         epsilon_eff -- the direction that matters: getting L_diff wrong
         on the small side silently weakens the very correction this term
         exists to supply (see doc/coulomb-heating-term-plan.md).
      6) A sweep at Eliseev et al.'s own Fig. 1(a) conditions (Ar, p = 1
         Torr, R = 3.85 cm used as the confinement length, epsilon_1 =
         11.55 eV), reported for visual comparison against the published
         curve -- NOT a numeric pass/fail, since the paper does not state
         the Te held fixed for that figure. The one number that IS exact
         ground truth from the formula itself, independent of Te, is the
         saturation plateau: epsilon_s = 0.5*11.55 = 5.775 eV, visually
         consistent with Fig. 1(a)'s ~5.5-6 eV plateau.

    Written 2026-09-08, alongside the Coulomb-heating energy-equation term
    (localEnergyEnergyModel, genericPlasmaPropertyTemplates.H).

\*---------------------------------------------------------------------------*/

#include "coulombHeatingClosure.H"

#include <cstdio>
#include <cmath>
#include <string>

using namespace Foam;

static int nFail = 0;
static int nRun  = 0;

static void check(const char* name, bool ok, const std::string& detail)
{
    std::printf("  [%s] %-58s (%s)\n", ok ? "ok  " : "FAIL", name,
                detail.c_str());
    ++nRun;
    if (!ok) ++nFail;
}

static std::string fmt(const char* f, double a, double b = 0)
{
    char buf[256];
    std::snprintf(buf, sizeof(buf), f, a, b);
    return buf;
}

static bool close(double a, double b, double tol)
{
    return std::fabs(a - b) <= tol*std::max(std::fabs(a), std::fabs(b));
}

int main()
{
    using namespace coulombHeatingClosure;

    std::printf("\ntestCoulombHeating -- unit bed for eq. (14)-(15) of"
                " Eliseev, Bogdanov & Kudryavtsev,\n"
                "Phys. Plasmas 24, 093503 (2017)\n\n");

    // ------------------------------------------------------------------
    // 1) nu_ee regression, against the hand-verified reference value.
    //
    //    n_e = 1e11 cm^-3 = 1e17 m^-3, Te = 0.2 eV, Lambda = 10.
    //    Confirmed 2026-09-08 with the user, reading eq. (14) directly
    //    from the PDF: "3.7 n_e[cm^-3] Lambda / ((e/k_b) Te[eV])^1.5",
    //    i.e. the denominator is Te IN KELVIN. Expected ~3.31e7 s^-1
    //    (1.7% from the standard NRL Plasma Formulary value of 3.25e7).
    // ------------------------------------------------------------------
    {
        const scalar ne = 1e17;   // m^-3
        const scalar Te = 0.2;    // eV
        const scalar lnLambda = 10.0;

        const scalar nuee = nuEE(ne, Te, lnLambda);
        const scalar expected = 3.309e7; // s^-1, hand-computed reference

        check
        (
            "nu_ee regression (n_e=1e11cm^-3, Te=0.2eV, Lambda=10)",
            close(nuee, expected, 0.02),
            fmt("got %.4e s^-1, expected ~%.4e s^-1", nuee, expected)
        );
    }

    // Shared "generic" parameters for the tests below: argon.
    const scalar coulombLog = 10.0;
    const scalar epsilon1 = 11.55;              // eV, Ar first excitation
    const scalar backgroundAtomMass = 39.948*1.66053906660e-27; // kg
    const scalar secondaryEnergyFraction = 0.5;  // Eliseev's Ref. 51 default
    const scalar muE = 10.0;                     // m^2/(V.s), representative
    const scalar meanE_typical = 0.3;            // eV -> Te = 0.2 eV

    // ------------------------------------------------------------------
    // 2) epsilon_eff -> 0 as n_e -> 0.
    // ------------------------------------------------------------------
    {
        const scalar epsEffLow = epsilonEff
        (
            1e6, meanE_typical, muE, coulombLog, epsilon1,
            backgroundAtomMass, secondaryEnergyFraction, GREAT
        );
        check
        (
            "epsilon_eff -> 0 as n_e -> 0",
            epsEffLow < 1e-6,
            fmt("epsilon_eff = %.3e eV at n_e = 1e6 m^-3", epsEffLow)
        );
    }

    // ------------------------------------------------------------------
    // 3) epsilon_eff -> epsilon_s at large n_e (nu_ee dominates).
    // ------------------------------------------------------------------
    {
        const scalar epsilonS = secondaryEnergyFraction*epsilon1;
        const scalar epsEffHigh = epsilonEff
        (
            1e22, meanE_typical, muE, coulombLog, epsilon1,
            backgroundAtomMass, secondaryEnergyFraction, GREAT
        );
        check
        (
            "epsilon_eff -> epsilon_s at large n_e",
            close(epsEffHigh, epsilonS, 0.01),
            fmt("epsilon_eff = %.4f eV, epsilon_s = %.4f eV",
                epsEffHigh, epsilonS)
        );
    }

    // ------------------------------------------------------------------
    // 4) Monotonicity in n_e at fixed Te.
    // ------------------------------------------------------------------
    {
        bool monotonic = true;
        scalar prev = 0.0;
        for (int i = 8; i <= 22; ++i)
        {
            const scalar ne = std::pow(10.0, scalar(i));
            const scalar e = epsilonEff
            (
                ne, meanE_typical, muE, coulombLog, epsilon1,
                backgroundAtomMass, secondaryEnergyFraction, GREAT
            );
            if (e < prev - SMALL) monotonic = false;
            prev = e;
        }
        check
        (
            "epsilon_eff monotonically increasing in n_e",
            monotonic,
            "swept n_e = 1e8 .. 1e22 m^-3"
        );
    }

    // ------------------------------------------------------------------
    // 5) Smaller diffusion length must NEVER give a LARGER epsilon_eff --
    //    the direction that matters, since getting L_diff wrong small
    //    silently weakens the correction this term exists to supply.
    // ------------------------------------------------------------------
    {
        const scalar ne = 1e18; // an intermediate density, off both limits
        const scalar eLarge = epsilonEff
        (
            ne, meanE_typical, muE, coulombLog, epsilon1,
            backgroundAtomMass, secondaryEnergyFraction, GREAT
        );
        const scalar eSmall = epsilonEff
        (
            ne, meanE_typical, muE, coulombLog, epsilon1,
            backgroundAtomMass, secondaryEnergyFraction, 0.01 // 1 cm
        );
        check
        (
            "smaller diffusion length -> smaller (never larger) epsilon_eff",
            eSmall <= eLarge + SMALL,
            fmt("epsilon_eff(L=GREAT) = %.4f eV, epsilon_eff(L=1cm) = %.4f eV",
                eLarge, eSmall)
        );
    }

    // ------------------------------------------------------------------
    // 6) Sweep at Eliseev et al.'s own Fig. 1(a) conditions: Ar, p = 1
    //    Torr (N_a from the ideal gas law at 300 K), R = 3.85 cm as the
    //    confinement length. Reported for visual comparison -- NOT a
    //    numeric pass/fail (Te held fixed for that figure is not stated
    //    in the paper), except for the plateau, which IS exact ground
    //    truth from the formula regardless of Te.
    // ------------------------------------------------------------------
    {
        const scalar torr = 133.322;      // Pa
        const scalar kB = 1.380649e-23;
        const scalar Na = (1.0*torr)/(kB*300.0); // m^-3, ideal gas law

        std::printf("\n  Fig. 1(a) sweep (Ar, p=1 Torr, R=3.85 cm,"
                    " Te=0.2 eV held fixed):\n");
        std::printf("    ionisation degree   epsilon_eff [eV]\n");

        for (int i = -8; i <= -3; ++i)
        {
            const scalar ionisationDegree = std::pow(10.0, scalar(i));
            const scalar ne = ionisationDegree*Na;

            const scalar e = epsilonEff
            (
                ne, meanE_typical, muE, coulombLog, epsilon1,
                backgroundAtomMass, secondaryEnergyFraction, 0.0385
            );
            std::printf("    %.0e              %.4f\n",
                        double(ionisationDegree), double(e));
        }

        const scalar epsilonS = secondaryEnergyFraction*epsilon1;
        std::printf("    (plateau = epsilon_s = %.3f eV; Fig. 1(a) shows"
                    " ~5.5-6 eV -- visual, not numeric, agreement)\n\n",
                    epsilonS);
    }

    // ------------------------------------------------------------------
    // 7) ORDER-OF-MAGNITUDE cross-check against Eliseev et al.'s OWN
    //    reported H_secondary, 2026-09-09. Not a numeric pass/fail: the
    //    inputs are visual reads off log-scale figures (Figs. 2, 5a, 7),
    //    not stated coefficients, so a factor-of-2-3 spread is expected
    //    from read error alone, not from formula disagreement.
    //
    //    Case: Ar, p = 107 Pa, I = 5 mA, R = 3.85 cm (the "standard"
    //    tube radius used throughout the paper for the Ref. 11
    //    comparison cases; Fig. 2/7 do not restate it but nothing in the
    //    text suggests a different geometry for this case).
    //
    //    Fig. 7 (p=107 Pa curve) peaks at H_se ~ 1.5-2e21 eV/(m^3 s).
    //    Fig. 2(c)'s source function peaks at the SAME location, at
    //    S_fast ~ 1-2e21 m^-3/s (log-scale read, gridlines at
    //    1e16/1e18/1e20/1e22). Fig. 5(a) gives the discharge's peak n_e
    //    as ~2e16 m^-3 at 107 Pa (a cleaner linear-scale read than
    //    Fig. 2(b)'s log axis, though not guaranteed to be AT the exact
    //    same x as the S_fast peak). Fig. 2(d)'s Te near that same
    //    steeply-falling region (not yet at the deep-bulk plateau of
    //    ~0.35 eV) is read as ~0.5 eV.
    //
    //    This checks: does epsilonEff(), fed THEIR numbers, reproduce
    //    THEIR H_secondary magnitude -- i.e. is the closure itself
    //    faithful to eq. (16)/(17), independent of whether SoPLASMA's
    //    own S_iz can ever reach a comparable value in a real run.
    // ------------------------------------------------------------------
    {
        const scalar ne_eliseev = 2e16;      // m^-3, Fig. 5(a) peak, 107 Pa
        const scalar Te_eliseev = 0.5;       // eV, Fig. 2(d), near S_fast peak
        const scalar Sfast_eliseev = 1.5e21; // m^-3/s, Fig. 2(c) peak
        const scalar R_eliseev = 0.0385;     // m, tube radius throughout

        const scalar e = epsilonEff
        (
            ne_eliseev, Te_eliseev, muE, coulombLog, epsilon1,
            backgroundAtomMass, secondaryEnergyFraction, R_eliseev
        );
        const scalar Hsecondary = e*Sfast_eliseev;

        std::printf("\n  Cross-check against Eliseev Fig. 2/7 (Ar, p=107 Pa,"
                    " I=5 mA):\n");
        std::printf("    inputs (read off figures): n_e=%.1e m^-3,"
                    " Te=%.2f eV, S_fast=%.1e m^-3/s\n",
                    double(ne_eliseev), double(Te_eliseev),
                    double(Sfast_eliseev));
        std::printf("    epsilon_eff (this closure) = %.4f eV\n", double(e));
        std::printf("    H_secondary = epsilon_eff * S_fast = %.3e"
                    " eV/(m^3 s)\n", double(Hsecondary));
        std::printf("    Fig. 7 reports ~1.5-2e21 eV/(m^3 s) at this"
                    " condition -- %s\n\n",
                    (Hsecondary > 1e20 && Hsecondary < 1e22)
                  ? "same order of magnitude"
                  : "DIFFERENT order of magnitude -- investigate");
    }

    std::printf("\n%d/%d checks passed.\n\n", nRun - nFail, nRun);
    return nFail == 0 ? 0 : 1;
}

// ************************************************************************* //
