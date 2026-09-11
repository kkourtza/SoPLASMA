/*---------------------------------------------------------------------------*\
    Charged-particle wall-loss models for a 0-D (global) reactor.

    A bounded discharge loses charge to its walls, and in a global model that
    loss is an algebraic closure rather than a solved flux:

        dn_s/dt = S_s - nu_transp,s n_s

    This header supplies nu_transp for the families reviewed by

        L L Alves and A Tejero-del-Caz,
        "Charged-particle transport models for global models",
        Plasma Sources Sci. Technol. 32 (2023) 054003

    Equation numbers below are theirs.

    WHY THIS IS 0-D ONLY, AND MUST STAY THAT WAY
    --------------------------------------------
    An h-factor is a SURROGATE for spatial resolution: it stands in for the
    density profile between the plasma centre and the sheath edge. The coupled
    CFD solves that profile explicitly with drift-diffusion. Applying both
    would count the same loss twice, and it would not be obvious in the output
    -- the discharge would simply be harder to sustain than it should be.

    This is the one place where the 0-D reactor deliberately does NOT mirror
    the CFD. Everything else about the two is kept in step on purpose.

    WHY THE DEFAULT IS `none`
    -------------------------
    Every case validated in this framework so far is a nanosecond atmospheric
    discharge, where the wall-loss timescale is Lambda^2/D_amb ~ 8.5 us against
    a 100 ns run: about 1% of the charge, well inside the other uncertainties.
    Switching a model on by default would change validated results by a
    percent for no physical reason. It is the low-pressure, steady-state cases
    -- ICP, CCP, microwave, DC glow -- where transport is 50-90% of all charge
    loss and a global model without it has no operating point at all.

    WHAT THE PAPER'S OWN CONCLUSION IS
    ----------------------------------
    That no single formulation is right everywhere: the choice of model moves
    E/N by 20-60% in oxygen and 8-115% in helium, and species densities by up
    to 150%. So the models are offered as a SET with their validity domains
    stated, not ranked. Choosing one is a modelling decision the user makes and
    should report, exactly like choosing a cross-section database.

\*---------------------------------------------------------------------------*/

#ifndef wallLoss_H
#define wallLoss_H

#include "scalar.H"
#include "word.H"
#include "List.H"
#include <cmath>

namespace Foam
{
namespace wallLoss
{

// Physical constants, SI.
constexpr scalar kB_ = 1.380649e-23;
constexpr scalar qe_ = 1.602176634e-19;

//- Which formulation. `none` is the default and returns zero.
enum model
{
    wlNone,
    wlAmbipolar,          //- classical ambipolar, eq (10a) + (4a)
    wlEffectiveDiffusion, //- Self-Ewald/Ferreira-Ricard abacus, eq (5a-5c)
    wlQGM,                //- Quantemol/GlobalKin, ambipolar + thermal, eq (11d)
    wlHFactorLP,          //- low-pressure h-factor, eq (41a)
    wlHFactorHPEN         //- high-pressure / electronegative, eq (41b)
};


inline model modelFromWord(const word& w)
{
    if (w == "none")                return wlNone;
    if (w == "ambipolar")           return wlAmbipolar;
    if (w == "effectiveDiffusion")  return wlEffectiveDiffusion;
    if (w == "QGM")                 return wlQGM;
    if (w == "hFactorLP")           return wlHFactorLP;
    if (w == "hFactorHPEN")         return wlHFactorHPEN;
    return wlNone;
}


inline word wordFromModel(const model m)
{
    switch (m)
    {
        case wlAmbipolar:           return "ambipolar";
        case wlEffectiveDiffusion:  return "effectiveDiffusion";
        case wlQGM:                 return "QGM";
        case wlHFactorLP:           return "hFactorLP";
        case wlHFactorHPEN:         return "hFactorHPEN";
        default:                    return "none";
    }
}


//- The state a wall-loss model needs. Deliberately a plain struct of scalars:
//  every model here is a pure function of these, which is what makes them
//  testable against the paper's analytic limits without a reactor.
struct state
{
    scalar R = 0.01;        //- discharge radius [m]
    scalar L = 0.1;         //- discharge length [m]

    scalar Te = 2.0;        //- electron temperature [eV], = (2/3)<eps>
    scalar Tg = 300.0;      //- gas (and ion) temperature [K]

    scalar Mi = 5.31e-26;   //- representative positive-ion mass [kg]
    scalar muiN = 5.0e21;   //- reduced ion mobility mu*N [1/(V m s)]
    scalar N = 2.45e25;     //- gas number density [m^-3]

    scalar alpha = 0.0;     //- electronegativity n_-/n_e at the centre
    scalar krec = 1.0e-13;  //- ion-ion recombination coefficient [m^3/s]
    scalar nMinus = 0.0;    //- negative-ion density [m^-3]
    scalar nPlus = 0.0;     //- positive-ion density [m^-3]
};


//- Diffusion length, eq (4b). Zero-density boundary at the wall, separation of
//  variables in r and z.
inline scalar diffusionLength(const scalar R, const scalar L)
{
    const scalar iL = (L > 0) ? (M_PI/L) : 0.0;
    const scalar iR = (R > 0) ? (2.405/R) : 0.0;
    const scalar s = iL*iL + iR*iR;
    return (s > 0) ? 1.0/std::sqrt(s) : 0.0;
}


//- Ion free-diffusion coefficient from its mobility, eq (2a) with the Einstein
//  relation. D/mu = kT/e, with the ION temperature -- taken equal to the gas
//  temperature here, which is the usual global-model assumption and is what
//  Alves et al. adopt for both their cases (their section 5.1).
inline scalar Dion(const state& s)
{
    const scalar mui = s.muiN/s.N;
    return mui*kB_*s.Tg/qe_;
}


//- Ion mean free path, eq (2e), with the RMS thermal speed used by the
//  ambipolar family, v_th = sqrt(3 kT/M) [their ref. 22].
inline scalar lambdaIon(const state& s)
{
    const scalar mui = s.muiN/s.N;
    const scalar nu = qe_/(s.Mi*mui);          // ion-neutral collision freq.
    const scalar vth = std::sqrt(3.0*kB_*s.Tg/s.Mi);
    return (nu > 0) ? vth/nu : 0.0;
}


//- Bohm speed, eq (20).
inline scalar bohmSpeed(const state& s)
{
    return std::sqrt(qe_*s.Te/s.Mi);
}


//- gamma = Te/T+, eq (2c). Te is in eV, Tg in K.
inline scalar gammaRatio(const state& s)
{
    return (s.Tg > 0) ? (s.Te*qe_/(kB_*s.Tg)) : 0.0;
}


//- Classical ambipolar diffusion coefficient for the positive ions, eq (10a).
//  Reduces to Da+ = gamma D+ for gamma >> 1 and alpha -> 0, which is the
//  textbook electropositive result.
inline scalar DambipolarPlus(const state& s)
{
    const scalar g = gammaRatio(s);
    return Dion(s)*(1.0 + g*(1.0 + 2.0*s.alpha))/(1.0 + s.alpha*g);
}


//- Effective (Self-Ewald) diffusion coefficient, eq (5a)-(5c).
//
//  Bridges the classical ambipolar limit at high pressure and the free-fall
//  limit at low pressure, where ion inertia matters and the ambipolar
//  description over-predicts the loss. The abacus is written in terms of
//  Lambda/lambda+, and the electronegativity correction of Coche et al. maps
//  the EP result onto a weakly EN plasma.
inline scalar DeffectivePlus(const state& s)
{
    const scalar Da = DambipolarPlus(s);
    const scalar lam = lambdaIon(s);
    if (lam <= 0) return Da;

    const scalar g = gammaRatio(s);
    const scalar Lam = diffusionLength(s.R, s.L);

    // eq (5b), solved for the Self coordinate
    const scalar ratio = (Lam/lam)
        *std::sqrt((1.0 + s.alpha*g)/(1.0 + s.alpha));

    // eq (5c)
    const scalar f = std::pow(std::tanh(std::pow(ratio, 0.36)), 2.9);

    // eq (5a)
    return Da*f*(1.0 + s.alpha)*(1.0 + s.alpha)/(1.0 + s.alpha*g);
}


//- Godyak's low-pressure h-factors for an electropositive discharge,
//  eq (24a) and (24b). Retained because every h-factor model below reduces to
//  these in the LP/EP corner, which is the cheapest available check on them.
inline scalar hGodyakL(const state& s)
{
    const scalar lam = lambdaIon(s);
    return 0.86/std::sqrt(3.0 + s.L/(2.0*lam));
}

inline scalar hGodyakR(const state& s)
{
    const scalar lam = lambdaIon(s);
    return 0.8/std::sqrt(4.0 + s.R/lam);
}


//- Chabert's corrected axial h-factor for LP and electronegative plasmas,
//  eq (38). Reduces to Godyak (24a) at alpha0 << 1 and L/lambda+ << 1.
inline scalar hChabertL(const state& s)
{
    const scalar lam = lambdaIon(s);
    const scalar g = gammaRatio(s);
    const scalar gm = g;                 // T_- = T_+ assumed, their section 5.1
    const scalar a0 = s.alpha;

    const scalar pre = std::sqrt(
        (gm - 1.0)/(gm*(1.0 + a0)*(1.0 + a0)) + 1.0/gm);

    const scalar den = std::sqrt(
        3.0 + s.L/(2.0*lam)
      + std::sqrt(1.0 + a0)*(s.L/lam)*(s.L/lam)/(5.0*g));

    return pre*0.86/den;
}


//- Thorsteinsson-Gudmundsson h-factors, eq (33a)-(33e), which join the LP
//  parabolic-EN-core solution to the HP flat-topped one.
inline scalar hThorsteinsson
(
    const state& s,
    const bool axial
)
{
    const scalar lam = lambdaIon(s);
    const scalar g = gammaRatio(s);
    const scalar gm = g;
    const scalar a0 = s.alpha;
    const scalar Da = DambipolarPlus(s);
    const scalar uB = bohmSpeed(s);

    // eq (32a): eta = 2T+/(T+ + T-) = 1 when the two are equal
    const scalar eta = 1.0;

    scalar h0;
    if (axial)
    {
        // eq (33c)
        const scalar t = 0.86*uB*s.L/(M_PI*Da);
        h0 = 0.86/std::sqrt(3.0 + eta*s.L/(2.0*lam) + t*t);
    }
    else
    {
        // eq (33d). J1(2.405) = 0.5191474
        const scalar t = 0.8*uB*s.R/(2.405*0.5191474*Da);
        h0 = 0.8/std::sqrt(4.0 + eta*s.R/lam + t*t);
    }

    // eq (32b) with the mean speed (12c), and eq (33e)
    const scalar vth = std::sqrt(8.0*kB_*s.Tg/(M_PI*s.Mi));
    const scalar nStar = (15.0/56.0)*vth*eta*eta/(s.krec*lam);

    scalar hc = 0.0;
    if (s.nMinus > 0 && s.nPlus > 0)
    {
        hc = 1.0/(std::sqrt(gm)
                + std::sqrt(g)*std::sqrt(nStar)*s.nPlus
                  /std::pow(s.nMinus, 1.5));
    }

    const scalar ha = h0/(1.0 + a0);
    return std::sqrt(ha*ha + hc*hc);
}


//- Loss frequency [1/s] for the positive ions.
//
//  The ambipolar family returns D/Lambda^2, eq (4a). The h-factor family
//  returns 2 u_B (h_L/L + h_R/R), eq (21). They are different constructions of
//  the same quantity and are NOT expected to agree except in the limits -- the
//  spread between them is the modelling uncertainty the paper is about.
inline scalar nuTransport(const model m, const state& s)
{
    if (m == wlNone) return 0.0;

    const scalar Lam = diffusionLength(s.R, s.L);
    if (Lam <= 0) return 0.0;

    switch (m)
    {
        case wlAmbipolar:
            return DambipolarPlus(s)/(Lam*Lam);

        case wlEffectiveDiffusion:
            return DeffectivePlus(s)/(Lam*Lam);

        case wlQGM:
        {
            // eq (11d): ambipolar transport and thermal wall loss combined as
            // characteristic times, over the volume-averaged collection area.
            const scalar A = 2.0*M_PI*s.R*s.L + 2.0*M_PI*s.R*s.R;
            const scalar V = M_PI*s.R*s.R*s.L;
            const scalar Da = DambipolarPlus(s);
            const scalar vth = std::sqrt(8.0*kB_*s.Tg/(M_PI*s.Mi));
            const scalar gr = 1.0;              // recombination probability
            const scalar den = Lam*gr + 4.0*Da/vth;
            return (den > 0) ? (A/V)*Da*gr/den : 0.0;
        }

        case wlHFactorLP:
        {
            // eq (41a): Chabert axially, Thorsteinsson radially
            const scalar uB = bohmSpeed(s);
            return 2.0*uB*(hChabertL(s)/s.L
                         + hThorsteinsson(s, false)/s.R);
        }

        case wlHFactorHPEN:
        {
            // eq (41b): Thorsteinsson in both directions
            const scalar uB = bohmSpeed(s);
            return 2.0*uB*(hThorsteinsson(s, true)/s.L
                         + hThorsteinsson(s, false)/s.R);
        }

        default:
            return 0.0;
    }
}

} // namespace wallLoss
} // namespace Foam

#endif
