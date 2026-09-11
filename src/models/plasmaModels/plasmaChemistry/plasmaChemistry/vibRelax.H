/*---------------------------------------------------------------------------*\
    Vibrational-translational relaxation time of N2(v) in air.

    Landau-Teller relaxation needs a timescale, and a single constant is wrong
    by DECADES in a discharge: V-T relaxation of N2 by atomic OXYGEN is far
    faster than by N2 or O2, so tau collapses as the discharge dissociates its
    oxygen. That is why this is a function of composition and not a number, and
    it is the whole reason the reservoir cannot be drained with a fitted
    constant.

    HOW MUCH faster is itself strongly temperature dependent, and quoting a
    single figure for it -- an earlier version of this comment said "~1e3
    times" -- is misleading. Per collision partner, against N2:

        T [K]      300     500    1000    1500    2000    2700    3000
        O / N2    2.8e8   2.9e6   2.0e4   2.0e3   4.6e2   1.2e2   8.0e1

    So "~1e3" happens to be right near 1500 K and is wrong by six orders at
    room temperature. The enhancement collapses as T rises because the
    molecular channels are the ones with the steep Landau-Teller temperature
    dependence, while Popov's atomic rate goes only as T^2.1.

    MOLECULAR PARTNERS -- Millikan & White, J. Chem. Phys. 39 (1963) 3209:

        p0 tau_{k,m} = exp[ a (T^-1/3 - b) - 18.42 ]      [atm s]

    with (a,b) = (220, 0.03) for N2-N2 and (162, 0.03) for N2-O2 [Shao et al.,
    Appl. Energy Combust. Sci. 19 (2024) 100280, from Colgan & Levitt 1967].

    THE ATOMIC-OXYGEN CHANNEL is taken from Popov's rate constant directly,

        k = 4.5e-21 (T/300)^2.1  m^3/s     [Popov, J. Phys. D 44 (2011) 285201]
        tau_{N2,O} = 1/(k n_O)

    rather than from the fitted form 488.5/(p0 T^1.1) quoted alongside it in
    Shao et al. eq. (7). Those two disagree by a factor of 101325 -- exactly Pa
    per atm -- so that constant requires p0 in PASCALS while the surrounding
    text says atm. Implemented as written with atm, V-T relaxation comes out
    1e5 times too slow and the reservoir never empties.

    The general lesson, and the reason this is written the way it is: prefer
    the underlying physical quantity to a fitted restatement of it.

    MIXING RULE, Millikan-White extended to several partners:

        1/tau_k = SUM_m X_m / tau_{k,m}

\*---------------------------------------------------------------------------*/

#ifndef vibRelax_H
#define vibRelax_H

#include "scalar.H"
#include <cmath>

namespace Foam
{
namespace vibRelax
{

//- Millikan-White relaxation time of N2 against one molecular partner [s].
//  p_atm is the TOTAL pressure in atmospheres; the partner's mole fraction is
//  applied by the mixing rule, not here.
inline scalar millikanWhite
(
    const scalar T,
    const scalar p_atm,
    const scalar a,
    const scalar b
)
{
    const scalar Tm13 = std::pow(T, -1.0/3.0);
    return std::exp(a*(Tm13 - b) - 18.42)/std::max(p_atm, scalar(1e-30));
}


//- Frequency of the N2-O V-T channel [1/s], from Popov's rate constant.
inline scalar nuVT_N2_O(const scalar T, const scalar nO)
{
    return (nO > 0) ? 4.5e-21*std::pow(T/300.0, 2.1)*nO : 0.0;
}


//- V-T relaxation time of N2(v) in an N2/O2/O mixture [s].
inline scalar tauVT_N2
(
    const scalar T,
    const scalar p_atm,
    const scalar xN2,
    const scalar xO2,
    const scalar nO            //- atomic oxygen number density [m^-3]
)
{
    scalar inv = 0.0;
    if (xN2 > 0) inv += xN2/millikanWhite(T, p_atm, 220.0, 0.03);   // N2-N2
    if (xO2 > 0) inv += xO2/millikanWhite(T, p_atm, 162.0, 0.03);   // N2-O2
    inv += nuVT_N2_O(T, nO);
    return (inv > 0) ? 1.0/inv : GREAT;
}

} // namespace vibRelax
} // namespace Foam

#endif
