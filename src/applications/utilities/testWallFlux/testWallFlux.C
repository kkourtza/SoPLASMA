/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.
    Copyright (C) 2026
    GNU General Public License v3 or later.

Application
    testWallFlux

Description
    Unit bed for the WALL-FLUX CLOSURE of Hagelaar's HDR chapter 6
    (Literature/), eqs (6.1), (6.2), (6.3), (6.6), (6.7), (6.8), (6.14),
    (6.15).

    It exercises the SHIPPED algebra -- `electronDDWallFluxMixed`'s static
    `hagelaarClosure()` and `hagelaarEnergyWeight()`, the same functions the
    boundary condition calls per face -- against ANALYTIC GROUND TRUTH. That
    is the point of a unit bed over a discharge case: eq. (6.7) and the
    r -> 1 energy limit are exact numbers the source states, where a CFD run
    could only show that two diagnostics agree with each other.

    Written 2026-09-04, with the reflection closure it verifies. Extended
    2026-09-06 with the MIXED CONDITION'S f for both flux schemes -- the
    closure algebra and the condition that imposes it are separate questions,
    and only the first was covered.

\*---------------------------------------------------------------------------*/

#include "electronDDWallFluxMixedFvPatchScalarField.H"
#include "plasmaConstants.H"

#include <cstdio>
#include <cmath>
#include <string>

using namespace Foam;

typedef electronDDWallFluxMixedFvPatchScalarField wallBC;

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

//- Relative closeness, with an absolute floor so exact zeros compare.
static bool close(double a, double b, double tol = 1e-12)
{
    return std::fabs(a - b) <= tol*std::max(1.0, std::fabs(b));
}

//- Convenience wrapper around the shipped kernel.
static double W(double A, double Dd, double GcOverN, double r)
{
    scalar w = 0, g = 0;
    wallBC::hagelaarClosure(A, Dd, GcOverN, r, w, g);
    return w;
}

static double gOf(double A, double Dd, double GcOverN, double r)
{
    scalar w = 0, g = 0;
    wallBC::hagelaarClosure(A, Dd, GcOverN, r, w, g);
    return g;
}


int main()
{
    const double pi = constant::mathematical::pi;
    const double sq = std::sqrt(pi);

    std::printf
    (
        "testWallFlux: Hagelaar HDR ch. 6 wall-flux closure,"
        " analytic ground truth\n\n"
    );

    // A representative thermal base: 2 eV electrons.
    // vT = sqrt(2 e Te/m_e); A = vT/sqrt(pi) = eq. (6.6).
    // vT = sqrt(2 k T/m) -- the HDR writes it as sqrt(2 e T/m) with T in
    // volts, which is the same quantity.
    const double TeK = 23208.6;                        // ~ 2 eV
    const double vT = std::sqrt
    (
        2.0*constant::plasma::kappaBoltzmann.value()*TeK
       /constant::plasma::eMass.value()
    );
    const double A = vT/sq;

    std::printf("  thermal base A = vT/sqrt(pi) = %.6e m/s  (Te = %.0f K"
                " ~ 2 eV)\n\n", A, TeK);

    // ---------------------------------------------------------------- (6.3)
    std::printf("=== eqs (6.3) and (6.6): the thermal base\n");
    {
        // (6.3) = vT/(2 sqrt(pi)) = (1/4) sqrt(8 e Te / pi m)
        const double e63 = 0.25*std::sqrt
        (
            8.0*constant::plasma::kappaBoltzmann.value()*TeK
           /(pi*constant::plasma::eMass.value())
        );
        check("eq. (6.6) is EXACTLY twice eq. (6.3)",
              close(A, 2.0*e63), fmt("A/e63 = %.15f", A/e63));

        check("(6.6) = 0.5*sqrt(8 k T / pi m), the shipped form",
              close(A, 0.5*std::sqrt
              (
                  8.0*constant::plasma::kappaBoltzmann.value()*TeK
                 /(pi*constant::plasma::eMass.value())
              )),
              "matches calcThermalVelocity");
    }

    // ---------------------------------------------------------------- (6.7)
    std::printf("\n=== eq. (6.7): reflection-only wall flux -- THE ground truth\n");
    {
        // Gamma.n = (1/sqrt(pi)) (1-r)/(1+r) n vT, i.e. W = A (1-r)/(1+r).
        bool all = true;
        double worst = 0;
        for (double r : {0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.99})
        {
            const double got = W(A, 0.0, 0.0, r);
            const double ref = (1.0/sq)*(1.0 - r)/(1.0 + r)*vT;
            const double rel = std::fabs(got - ref)/ref;
            worst = std::max(worst, rel);
            if (rel > 1e-14) all = false;
        }
        check("W = (1/sqrt(pi))(1-r)/(1+r) n vT for r in [0, 0.99]",
              all, fmt("worst rel err = %.3e over 8 values", worst));

        check("r = 0 recovers the bare (6.6) base",
              close(W(A, 0.0, 0.0, 0.0), A),
              fmt("W = %.9e", W(A, 0.0, 0.0, 0.0)));

        check("W -> 0 as r -> 1 (a perfectly reflecting wall loses nothing)",
              W(A, 0.0, 0.0, 1.0) == 0.0,
              fmt("W(r=1) = %.3e", W(A, 0.0, 0.0, 1.0)));

        // Monotone in r, which (6.7) requires and a sign slip would break.
        bool mono = true;
        for (double r = 0.0; r < 0.95; r += 0.05)
        {
            if (!(W(A, 0.0, 0.0, r + 0.05) < W(A, 0.0, 0.0, r))) mono = false;
        }
        check("W falls monotonically with r", mono, "20 steps");
    }

    // ---------------------------------------------------- the r = 0 identity
    std::printf("\n=== r = 0 must reproduce the pre-2026-09-04 form\n");
    {
        // old: max(A - Gc/n, 0) + Dd    new: (1-r) max((A + Dd - Gc/n)/(1+r), Dd)
        // The claim is the algebraic identity
        //     max(A + Dd - Gc/n, Dd) == Dd + max(A - Gc/n, 0).
        // ALGEBRAIC, NOT BITWISE: the two forms associate the additions
        // differently, so they agree to rounding and not to the last bit.
        // Asserting bit equality was this bed's own first failure, 2026-09-04.
        double worstAbs = 0, worstRel = 0;
        int n = 0;
        for (double dd = 0.0; dd <= 3.0*A; dd += 0.25*A)
        {
            for (double gc = 0.0; gc <= 3.0*A; gc += 0.1*A)
            {
                const double got = W(A, dd, gc, 0.0);
                const double old = std::max(A - gc, 0.0) + dd;
                worstAbs = std::max(worstAbs, std::fabs(got - old));
                worstRel = std::max
                (
                    worstRel,
                    std::fabs(got - old)/std::max(1.0, std::fabs(old))
                );
                ++n;
            }
        }
        check("agrees with the old form to floating-point rounding",
              worstRel < 1e-14,
              fmt("%.0f points, max rel = %.2e", double(n), worstRel));

        check("and the absolute error is far below any physical scale",
              worstAbs < 1e-8*A,
              fmt("max |diff| = %.2e m/s vs A = %.2e m/s", worstAbs, A));
    }

    // ------------------------------------------------------- (6.1) + (6.6)
    std::printf("\n=== eqs (6.1)+(6.6): the FACTOR 2 on emission\n");
    {
        // Gamma = n w_w - Gamma_w = n(A - Gamma_w/n) - Gamma_w = n A - 2 Gamma_w
        const double gc = 0.2*A;                     // Gc/n, well below A
        const double got = W(A, 0.0, gc, 0.0) - gc;  // Gamma/n = W - Gc/n
        check("net Gamma/n = A - 2 Gc/n at r = 0",
              close(got, A - 2.0*gc),
              fmt("got %.9e, want %.9e", got, A - 2.0*gc));

        // With reflection the exact closure gives the creation factor 2/(1+r),
        // which is what the Implicit family now carries.
        bool all = true;
        double worst = 0;
        for (double r : {0.0, 0.2, 0.5, 0.8})
        {
            const double net = W(A, 0.0, gc, r) - gc;
            const double ref = (1.0 - r)/(1.0 + r)*A - 2.0/(1.0 + r)*gc;
            const double rel = std::fabs(net - ref)/std::fabs(ref);
            worst = std::max(worst, rel);
            if (rel > 1e-13) all = false;
        }
        check("net Gamma/n = [(1-r)/(1+r)]A - [2/(1+r)]Gc/n",
              all, fmt("worst rel err = %.3e", worst));
    }

    // ---------------------------------------------------- the (6.8) clamp
    std::printf("\n=== eq. (6.8): the clamp, and its floor at Dd\n");
    {
        check("strong creation clamps W to zero when Dd = 0",
              W(A, 0.0, 5.0*A, 0.0) == 0.0,
              fmt("W = %.3e", W(A, 0.0, 5.0*A, 0.0)));

        check("the clamp FLOOR is Dd, not 0 -- drift is outside the max",
              close(W(A, 0.7*A, 5.0*A, 0.0), 0.7*A),
              fmt("W = %.6e, Dd = %.6e", W(A, 0.7*A, 5.0*A, 0.0), 0.7*A));

        check("floor carries the (1-r) factor too",
              close(W(A, 0.7*A, 9.0*A, 0.5), 0.5*0.7*A),
              fmt("W = %.6e, (1-r)Dd = %.6e",
                  W(A, 0.7*A, 9.0*A, 0.5), 0.5*0.7*A));

        // THE DERIVATION ITSELF, against an INDEPENDENT numerical solve.
        //
        // The load-bearing step is that the closed form solves the implicit
        // relation of eqs (6.2)+(6.8),
        //     w_w = max(A - (r*w_w + Gc/n), 0) + Dd
        // whose right-hand side is a contraction in w_w for r < 1. Iterating
        // it to convergence gives w_w without using any of the algebra under
        // test, so this compares the shipped closure against the SOURCE
        // EQUATIONS rather than against my rearrangement of them.
        //
        // This replaced a continuity probe that could not fail: max(.,c) is
        // continuous for ANY floor c, so the "jump" it measured was only the
        // finite-difference slope. Rejected 2026-09-04 by its own liveness
        // control.
        auto fixedPoint = [&](double Dd, double gc, double r)
        {
            double w = A + Dd;
            for (int it = 0; it < 20000; ++it)
            {
                const double next = std::max(A - (r*w + gc), 0.0) + Dd;
                if (std::fabs(next - w) < 1e-15*std::max(1.0, A)) return next;
                w = next;
            }
            return w;
        };

        bool solves = true;
        double worstRel = 0;
        int nPts = 0;
        for (double r : {0.0, 0.1, 0.3, 0.5, 0.7, 0.9})
        {
            for (double ddf : {0.0, 0.2, 0.4, 1.0, 2.5})
            {
                for (double gcf = 0.0; gcf <= 4.0; gcf += 0.05)
                {
                    const double Dd = ddf*A, gc = gcf*A;
                    const double wRef = fixedPoint(Dd, gc, r);
                    const double got  = W(A, Dd, gc, r)/(1.0 - r);  // w_w
                    const double rel  =
                        std::fabs(got - wRef)/std::max(1.0, std::fabs(wRef));
                    worstRel = std::max(worstRel, rel);
                    if (rel > 1e-10) solves = false;
                    ++nPts;
                }
            }
        }
        check("the closed form SOLVES the implicit eq. (6.2)+(6.8) relation",
              solves, fmt("%.0f (r, Dd, Gc) points, max rel = %.2e",
                          double(nPts), worstRel));

        // LIVENESS. A comparison that cannot fail has tested nothing. The two
        // errors the derivation had to rule out must both be CONVICTED, and
        // over the SAME domain as the correctness claim -- a single hand-
        // picked point lands in the clamped branch, where every variant
        // correctly returns Dd and nothing can be told apart. That is how the
        // first version of this control passed the wrong form, 2026-09-04.
        {
            double worstFloor = 0, worstScale = 0;
            for (double r : {0.0, 0.1, 0.3, 0.5, 0.7, 0.9})
            {
                for (double ddf : {0.0, 0.2, 0.4, 1.0, 2.5})
                {
                    for (double gcf = 0.0; gcf <= 4.0; gcf += 0.05)
                    {
                        const double Dd = ddf*A, gc = gcf*A;
                        const double wRef = fixedPoint(Dd, gc, r);

                        // (a) floor at 0 instead of Dd
                        const double wrongFloor =
                            std::max((A + Dd - gc)/(1.0 + r), 0.0);

                        // (b) divide only A by (1+r) -- the pre-2026-09-04
                        //     form, which also subtracted Gc AFTER scaling
                        const double wrongScale =
                            std::max(A*(1.0 - r)/(1.0 + r) - gc, 0.0) + Dd;

                        worstFloor =
                            std::max(worstFloor, std::fabs(wrongFloor - wRef));
                        worstScale =
                            std::max(worstScale, std::fabs(wrongScale - wRef));
                    }
                }
            }

            check("...and it CONVICTS a floor of 0 instead of Dd",
                  worstFloor > 1e-2*A,
                  fmt("max deviation %.3e m/s = %.2f A", worstFloor,
                      worstFloor/A));

            check("...and it CONVICTS scaling A alone by (1-r)/(1+r)",
                  worstScale > 1e-2*A,
                  fmt("max deviation %.3e m/s = %.2f A", worstScale,
                      worstScale/A));
        }

        bool mono = true;
        for (double gc = 0.0; gc < 3.0*A; gc += 0.05*A)
        {
            if (W(A, 0.0, gc + 0.05*A, 0.3) > W(A, 0.0, gc, 0.3) + 1e-30)
                mono = false;
        }
        check("W never RISES with creation", mono, "60 steps");
    }

    // --------------------------------------------------------------- (6.15)
    std::printf("\n=== eq. (6.15): the electron-energy weight\n");
    {
        const double f0 = 5.0/3.0;

        check("eps_w/eps = 5/3 at zero creation and zero reflection",
              close(wallBC::hagelaarEnergyWeight(f0, gOf(A, 0.0, 0.0, 0.0)),
                    5.0/3.0),
              fmt("got %.15f",
                  wallBC::hagelaarEnergyWeight(f0, gOf(A, 0.0, 0.0, 0.0))));

        check("that is EXACTLY 5/4 x the old (6.14) weight of 4/3",
              close(5.0/3.0, 1.25*4.0/3.0),
              "5/3 = (5/4)(4/3)");

        // THE PUBLISHED LIMIT. Hagelaar: substituting the reflection-only
        // Gamma_w and letting r -> 1 yields eps_w -> 2 Te, "as for a centred
        // Maxwellian". With eps = (3/2)Te that is eps_w/eps -> 4/3.
        const double wLim = wallBC::hagelaarEnergyWeight(f0, gOf(A, 0.0, 0.0, 1.0));
        check("r -> 1 gives eps_w -> 2 Te (eps_w/eps -> 4/3), as the text asserts",
              close(wLim, 4.0/3.0),
              fmt("got %.15f, want %.15f", wLim, 4.0/3.0));

        // gRatio for reflection only must be r/(1+r) -- the quantity the text
        // substitutes into (6.15).
        bool all = true;
        double worst = 0;
        for (double r : {0.0, 0.25, 0.5, 0.75, 1.0})
        {
            const double got = gOf(A, 0.0, 0.0, r);
            const double ref = r/(1.0 + r);
            worst = std::max(worst, std::fabs(got - ref));
            if (std::fabs(got - ref) > 1e-14) all = false;
        }
        check("gRatio = r/(1+r) for reflection only",
              all, fmt("max |diff| = %.3e", worst));

        check("eps_w falls monotonically with reflection",
              wallBC::hagelaarEnergyWeight(f0, gOf(A, 0, 0, 0.2))
            > wallBC::hagelaarEnergyWeight(f0, gOf(A, 0, 0, 0.8)),
              "r = 0.2 vs 0.8");

        check("eps_w/eps is CLAMPED at zero under extreme creation",
              wallBC::hagelaarEnergyWeight(f0, gOf(A, 0.0, 50.0*A, 0.0)) == 0.0,
              "Gc/n = 50 A");

        // A negative weight would reverse the sign of the wall energy flux --
        // the failure the clamp exists to prevent. Sweep for it.
        bool nonNeg = true;
        for (double r = 0.0; r <= 0.9; r += 0.1)
        {
            for (double gc = 0.0; gc <= 20.0*A; gc += 0.1*A)
            {
                if (wallBC::hagelaarEnergyWeight(f0, gOf(A, 0.3*A, gc, r)) < 0.0)
                    nonNeg = false;
            }
        }
        check("eps_w/eps >= 0 everywhere in (r, Gc/n)", nonNeg,
              "10 x 200 sweep");
    }

    // -------------------------------------------------------- degenerate A
    std::printf("\n=== degenerate inputs\n");
    {
        scalar w = 0, g = 0;
        wallBC::hagelaarClosure(0.0, 0.0, 0.0, 0.0, w, g);
        check("A = 0 gives W = 0 and gRatio = 0, not a division by zero",
              w == 0.0 && g == 0.0, fmt("W = %.1e, g = %.1e", w, g));

        wallBC::hagelaarClosure(0.0, 0.0, 1.0e6, 0.0, w, g);
        check("A = 0 with creation: gRatio guarded to 0",
              g == 0.0, fmt("g = %.1e", g));
    }

    // ------------------------------- the mixed condition's f, BOTH schemes
    //
    // THE GAP THIS CLOSES. Everything above verifies the closure ALGEBRA. It
    // says nothing about whether the assembled boundary condition actually
    // IMPOSES that flux -- and that is a separate question, because the mixed
    // condition sets a face VALUE (n_p = (1-f) n_c) and the flux is whatever
    // the discretisation's own face-flux formula makes of it. f must therefore
    // be chosen PER SCHEME, and each choice has to be checked against
    // Gamma = n_p*W under ITS OWN extraction formula.
    //
    // Added 2026-09-06 after reporting a defect that was not one: I had
    // compared the two branches' f against a SINGLE extraction formula, which
    // is meaningless. This is the check that would have caught that
    // immediately -- both branches pass it exactly.
    std::printf("\n=== the mixed condition's f imposes the closure flux\n");
    {
        // Bernoulli, and the identity the SG derivation turns on.
        auto Bern = [](scalar x) -> scalar
        {
            const scalar ax = std::fabs(x);
            if (ax < 1e-4) return 1.0 - 0.5*x + (x*x)/12.0;
            if (x > 100.0)  return 0.0;
            if (x < -100.0) return -x;
            return x/(std::exp(x) - 1.0);
        };

        bool idOK = true;
        for (scalar x = -50.0; x <= 50.0; x += 0.37)
        {
            if (std::fabs((Bern(-x) - Bern(x)) - x) > 1e-9*std::max(1.0, std::fabs(x)))
            {
                idOK = false;
            }
        }
        check("Bernoulli identity Bern(-x) - Bern(x) = x, which the SG f needs",
              idOK, "271 points over x = -50..50");

        const scalar g = 1.0;          // D/delta, in velocity units
        const scalar A = 1.0;          // thermal wall speed
        scalar worstStd = 0, worstSG = 0;
        scalar peWorstStd = 0, peWorstSG = 0;
        int nPt = 0;

        const scalar peList[9] =
            {0.0, 0.01, 0.1, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0};
        const scalar rList[4] = {0.0, 0.1, 0.36, 0.8};

        for (int ip = 0; ip < 9; ++ip)
        {
            for (int ir = 0; ir < 4; ++ir)
            {
                const scalar Pe = peList[ip];
                const scalar r  = rList[ir];
                const scalar ud = Pe*g;

                // The closure's total loss speed, with the drift flux in.
                const scalar W = (1.0 - r)/(1.0 + r)*A + std::max(scalar(0), ud);
                const scalar uEff = W - ud;    // the bookkeeping identity
                const scalar nc = 1.0;
                ++nPt;

                // STANDARD: f = uEff/(D/delta + uEff); extracted as
                // diffusion + drift.
                {
                    const scalar f = uEff/(g + uEff);
                    const scalar np = (1.0 - f)*nc;
                    const scalar Gamma = g*(nc - np) + ud*np;
                    const scalar e = std::fabs(Gamma - np*W)
                                   / std::max(std::fabs(np*W), scalar(1e-300));
                    if (e > worstStd) { worstStd = e; peWorstStd = Pe; }
                }

                // SCHARFETTER-GUMMEL: f = uEff/(D/delta*Bern(Pe) + W);
                // extracted by ScharfetterGummel.H's boundary coefficients,
                // Gamma = (D/delta)[Bern(-Pe) n_c - Bern(Pe) n_p].
                {
                    const scalar f = uEff/(g*Bern(Pe) + W);
                    const scalar np = (1.0 - f)*nc;
                    const scalar Gamma = g*(Bern(-Pe)*nc - Bern(Pe)*np);
                    const scalar e = std::fabs(Gamma - np*W)
                                   / std::max(std::fabs(np*W), scalar(1e-300));
                    if (e > worstSG) { worstSG = e; peWorstSG = Pe; }
                }
            }
        }

        const std::string span =
            fmt(" over %g (Pe, r) points", scalar(nPt));

        check("STANDARD f imposes Gamma = n_p*W exactly",
              worstStd < 1e-12,
              fmt("worst rel. err %.2e at Pe = %g", worstStd, peWorstStd)
                  + span);

        check("SCHARFETTER-GUMMEL f imposes Gamma = n_p*W exactly, under ITS"
              " own extraction",
              worstSG < 1e-12,
              fmt("worst rel. err %.2e at Pe = %g", worstSG, peWorstSG)
                  + span);

        // LIVENESS, both ways: the check must CONVICT each branch's f when it
        // is used with the other branch's extraction formula. Without this,
        // the two checks above could both be passing on an identity that
        // holds for any f.
        {
            const scalar Pe = 100.0, r = 0.0, ud = Pe*g;
            const scalar W = (1.0 - r)/(1.0 + r)*A + std::max(scalar(0), ud);
            const scalar uEff = W - ud, nc = 1.0;

            const scalar fStd = uEff/(g + uEff);
            const scalar npStd = (1.0 - fStd)*nc;
            const scalar GstdViaSG =
                g*(Bern(-Pe)*nc - Bern(Pe)*npStd);
            check("...and it CONVICTS the standard f under SG extraction",
                  std::fabs(GstdViaSG - npStd*W)
                      > 0.01*std::fabs(npStd*W),
                  fmt("ratio %.3f at Pe = 100", GstdViaSG/(npStd*W)));

            const scalar fSG = uEff/(g*Bern(Pe) + W);
            const scalar npSG = (1.0 - fSG)*nc;
            const scalar GsgViaStd = g*(nc - npSG) + ud*npSG;
            check("...and it CONVICTS the SG f under standard extraction",
                  std::fabs(GsgViaStd - npSG*W)
                      > 0.01*std::fabs(npSG*W),
                  fmt("ratio %.3f at Pe = 100 -- THIS is the number I"
                      " misreported as a defect", GsgViaStd/(npSG*W)));
        }

        // And that the two f agree where they must: as Pe -> 0 the drift
        // vanishes, W -> uEff and Bern(0) = 1, so the schemes coincide.
        {
            const scalar Pe = 1e-8, r = 0.36, ud = Pe*g;
            const scalar W = (1.0 - r)/(1.0 + r)*A + std::max(scalar(0), ud);
            const scalar uEff = W - ud;
            const scalar fStd = uEff/(g + uEff);
            const scalar fSG  = uEff/(g*Bern(Pe) + W);
            check("the two f COINCIDE as Pe -> 0, as they must",
                  std::fabs(fStd - fSG) < 1e-7,
                  fmt("f_std = %.10f, f_SG = %.10f", fStd, fSG));
        }
    }

    std::printf("\n%s: %d checks, %d failed\n",
                nFail ? "FAILURES" : "ALL PASS", nRun, nFail);

    return nFail ? 1 : 0;
}

// ************************************************************************* //
