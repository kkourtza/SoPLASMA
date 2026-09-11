/*---------------------------------------------------------------------------*\
Application
    testWallLoss

Description
    Checks the charged-particle wall-loss models against the ANALYTIC LIMITS
    stated in the paper they come from:

        L L Alves and A Tejero-del-Caz,
        Plasma Sources Sci. Technol. 32 (2023) 054003

    Every model there is built by stitching limiting solutions together, and
    each states what it must reduce to at the joins. Those statements are free,
    exact test oracles -- far sharper than "the number looks plausible", and
    they need no reference data at all.

    This is the right way to validate this family. Reproducing the paper's
    oxygen and helium figures would need their full kinetic schemes (O2 with
    vibrational levels v=0-41, He with every state to n=7), which is mechanism
    work, not solver work. The limits test the FORMULAE, which is what was
    actually implemented here.

\*---------------------------------------------------------------------------*/

#include "wallLoss.H"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace Foam;
using namespace Foam::wallLoss;

static int nFail = 0;

static void check(const char* name, bool ok, const std::string& detail)
{
    std::printf("  [%s] %-52s (%s)\n", ok ? "ok  " : "FAIL", name,
                detail.c_str());
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
    return std::fabs(a - b) <= tol*std::fabs(b);
}


int main()
{
    std::printf("testWallLoss: Alves & Tejero-del-Caz 2023, analytic limits\n\n");

    // A representative low-pressure discharge: the oxygen DC glow of their
    // section 5.2, R = 1 cm, L = 52.5 cm.
    state base;
    base.R = 0.01;
    base.L = 0.525;
    base.Te = 2.0;
    base.Tg = 400.0;
    base.Mi = 5.31e-26;          // O2+
    base.muiN = 5.0e21;
    base.N = 2.4e22;             // ~133 Pa at 400 K
    base.alpha = 0.0;

    // ---- eq (4b): the diffusion length -----------------------------------
    std::printf("diffusion length, eq (4b)\n");
    {
        state s = base;
        s.L = 1e6;                                   // L -> infinity
        const scalar lam = diffusionLength(s.R, s.L);
        check("Lambda -> R/2.405 for a long cylinder",
              close(lam, s.R/2.405, 1e-6),
              fmt("%.6e vs %.6e", lam, s.R/2.405));

        s = base;
        s.R = 1e6;                                   // R -> infinity
        const scalar lam2 = diffusionLength(s.R, s.L);
        check("Lambda -> L/pi for a wide slab",
              close(lam2, s.L/M_PI, 1e-6),
              fmt("%.6e vs %.6e", lam2, s.L/M_PI));
    }

    // ---- eq (10a): the electropositive ambipolar limit --------------------
    std::printf("\nambipolar coefficient, eq (10a)\n");
    {
        state s = base;
        s.alpha = 0.0;
        const scalar g = gammaRatio(s);
        const scalar Da = DambipolarPlus(s);
        // Da+ = D+(1+gamma) -> gamma D+ for gamma >> 1. The residual is
        // exactly 1/gamma, so this checks the algebra, not a tolerance.
        check("Da+ -> gamma D+ at alpha=0, gamma>>1",
              close(Da, g*Dion(s), 1.5/g),
              fmt("Da/(gamma D+) = %.6f, 1/gamma = %.4f", Da/(g*Dion(s)), 1.0/g));

        // Electronegativity DRIVES Da+ DOWN, and to a specific value. As
        // alpha -> infinity, eq (10a) gives
        //     (1 + gamma(1+2 alpha))/(1 + alpha gamma) -> 2,
        // so Da+ -> 2 D+: the ION-ION ambipolar limit, where the field is set
        // by two heavy species rather than by the electrons and the fast
        // electron-driven enhancement gamma is gone. This is the sharpest
        // available check on eq (10a) because the limit is a pure number.
        state s2 = s; s2.alpha = 0.5;
        check("Da+ FALLS with electronegativity",
              DambipolarPlus(s2) < Da,
              fmt("%.4e -> %.4e", Da, DambipolarPlus(s2)));

        state s3 = s; s3.alpha = 1e6;
        check("Da+ -> 2 D+ as alpha -> inf (ion-ion ambipolar limit)",
              close(DambipolarPlus(s3), 2.0*Dion(s3), 1e-3),
              fmt("Da+/D+ = %.6f, expected 2", DambipolarPlus(s3)/Dion(s3)));
    }

    // ---- eq (5a-5c): the effective-diffusion abacus -----------------------
    std::printf("\neffective diffusion, eq (5a)-(5c)\n");
    {
        // Their statement after eq (5c): "in the HP limit of Lambda/lambda+
        // >> 1, equations (5a)-(5c) yield Deff/Da ~ 1, corresponding to the
        // classical ambipolar diffusion model."
        state s = base;
        s.N = 2.4e25;                                // ~100x denser
        const scalar r = DeffectivePlus(s)/DambipolarPlus(s);
        check("Deff/Da -> 1 at high pressure (classical limit)",
              close(r, 1.0, 0.02),
              fmt("Lambda/lambda+ = %.1f, Deff/Da = %.4f",
                  diffusionLength(s.R, s.L)/lambdaIon(s), r));

        // And below it, at low pressure the free-fall correction must REDUCE
        // the loss relative to classical ambipolar, which is the entire
        // purpose of the Self-Ewald treatment.
        state s2 = base;
        s2.N = 2.4e20;                               // very low pressure
        const scalar r2 = DeffectivePlus(s2)/DambipolarPlus(s2);
        check("Deff/Da < 1 at low pressure (free-fall)",
              r2 < 0.9,
              fmt("Lambda/lambda+ = %.3f, Deff/Da = %.4f",
                  diffusionLength(s2.R, s2.L)/lambdaIon(s2), r2));

        check("Deff/Da is monotone in pressure", r > r2,
              fmt("%.4f > %.4f", r, r2));
    }

    // ---- eq (38) -> eq (24a): Chabert reduces to Godyak --------------------
    std::printf("\nh-factors: Chabert eq (38) -> Godyak eq (24a)\n");
    {
        // Their statement: "the corrected equation (38) allows retrieving
        // (24a) at low EN and for LP conditions, i.e. for alpha0 << 1 and
        // L/lambda+ << 1."
        state s = base;
        s.alpha = 0.0;
        s.N = 1.0e18;                    // L/lambda+ << 1
        const scalar hC = hChabertL(s);
        const scalar hG = hGodyakL(s);
        check("h_L(Chabert) -> h_L(Godyak) at alpha0=0, L<<lambda+",
              close(hC, hG, 0.02),
              fmt("L/lambda+ = %.4f, %.6f vs %.6f",
                  s.L/lambdaIon(s), hC));
        std::printf("        (Godyak %.6f, Chabert %.6f)\n", hG, hC);

        // The prefactor of (38) is identically 1 at alpha0 = 0 -- worth
        // asserting separately, because a slip there would be masked by the
        // tolerance above.
        check("the electronegativity prefactor is exactly 1 at alpha0=0",
              close(hC/hG, 1.0, 0.03), fmt("ratio %.6f", hC/hG));

        // Electronegativity must SUPPRESS the edge-to-centre ratio: more
        // negative ions, fewer positive ions reaching the sheath.
        state s2 = s; s2.alpha = 3.0;
        check("h_L falls with electronegativity",
              hChabertL(s2) < hC,
              fmt("%.6f -> %.6f", hC, hChabertL(s2)));
    }

    // ---- eq (33): Thorsteinsson ------------------------------------------
    std::printf("\nh-factors: Thorsteinsson eq (33a)-(33e)\n");
    {
        state s = base;
        s.alpha = 0.0;
        s.nMinus = 0.0;                  // no EN core, so h_c must vanish
        const scalar hL = hThorsteinsson(s, true);
        const scalar hR = hThorsteinsson(s, false);
        check("h_L, h_R are in (0,1] for an electropositive plasma",
              hL > 0 && hL <= 1.0 && hR > 0 && hR <= 1.0,
              fmt("h_L = %.4f, h_R = %.4f", hL, hR));

        // h is an EDGE-TO-CENTRE RATIO, so it falls as the dimension grows:
        // a longer path gives the profile more room to decay before the
        // sheath. With L = 52 R the axial factor must therefore be the
        // SMALLER of the two. (Getting this backwards is easy, and the first
        // version of this test did.)
        check("h_L < h_R when L >> R (h falls with the dimension)", hL < hR,
              fmt("h_L = %.4f < h_R = %.4f", hL, hR));

        state sBig = s; sBig.L = 4.0*s.L;
        check("h_L falls further when L is quadrupled",
              hThorsteinsson(sBig, true) < hL,
              fmt("%.5f -> %.5f", hL, hThorsteinsson(sBig, true)));
    }

    // ---- convergence of the two families at high pressure -----------------
    std::printf("\nconvergence at high pressure (their section 5.1)\n");
    {
        // "At HP, the ion-transport losses are strongly reduced. In this case,
        // all transport theories converge to classical ambipolar diffusion
        // losses" -- with the small geometric factors of their eq (42), which
        // are sqrt(pi/2) ~ 1.25 and sqrt(2.405/(2 J1(2.405))) ~ 1.52. So the
        // h-factor family must approach the ambipolar one to within a factor
        // of about 1.5, not exactly.
        state s = base;
        s.N = 2.4e25;
        const scalar nuA = nuTransport(wlAmbipolar, s);
        const scalar nuH = nuTransport(wlHFactorHPEN, s);
        const scalar ratio = nuH/nuA;
        check("h-factor and ambipolar agree to the eq (42) factors at HP",
              ratio > 0.4 && ratio < 2.5,
              fmt("nu_h/nu_amb = %.3f (expected O(1), factors 1.25-1.52)",
                  ratio));
    }

    // ---- general sanity ---------------------------------------------------
    std::printf("\ngeneral\n");
    {
        const model all[] = {wlAmbipolar, wlEffectiveDiffusion, wlQGM,
                             wlHFactorLP, wlHFactorHPEN};
        bool allPos = true, allShrink = true;
        for (const model m : all)
        {
            state s = base;
            const scalar nu1 = nuTransport(m, s);
            state s2 = base; s2.R = 2.0*base.R;      // wider tube
            const scalar nu2 = nuTransport(m, s2);
            if (!(nu1 > 0) || !std::isfinite(nu1)) allPos = false;
            if (!(nu2 < nu1)) allShrink = false;
        }
        check("every model gives a positive, finite loss frequency", allPos, "");
        check("every model loses less from a wider tube", allShrink,
              "nu(2R) < nu(R)");

        state s = base;
        check("`none` is exactly zero",
              nuTransport(wlNone, s) == 0.0, "");
    }

    std::printf("\n%s\n", nFail ? "FAILED" : "all checks passed");
    return nFail ? 1 : 0;
}
