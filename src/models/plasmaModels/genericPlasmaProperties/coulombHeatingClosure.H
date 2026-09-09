/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Namespace
    Foam::coulombHeatingClosure

Description
    Free-standing physics core of `CoulombHeatingProperty`
    (genericPlasmaPropertyTemplates.H), pulled out on purpose so it can be
    unit-tested directly against analytic ground truth -- see
    src/applications/utilities/testCoulombHeating -- without constructing
    an fvMesh, a dictionary, or the plasmaPropertyEvaluator machinery at
    all. Same pattern as electronDDWallFluxMixed's static hagelaarClosure(),
    tested by testWallFlux.

    epsilon_eff, the order-of-magnitude estimate of eq. (15) in Eliseev,
    Bogdanov & Kudryavtsev, "Slow electron energy balance for hybrid models
    of direct-current glow discharges", Phys. Plasmas 24, 093503 (2017):

        epsilon_eff = [ 2 nu_ee / (2 nu_ee + (2 me/ma) nu_ea + 1/tau_diff) ]
                      * epsilon_s

    nu_ee is the electron-electron (Coulomb) collision frequency, eq. (14)
    of the same paper, verified 2026-09-08 against the standard NRL Plasma
    Formulary value (2.91e-6 n_e[cm^-3] lnLambda / Te[eV]^1.5): the two
    agree to 1.7% once eq. (14)'s "(e/k_b) Te[eV]" factor is read as the
    standard eV-to-Kelvin conversion (e/k_b ~ 11604 K/eV), i.e. that
    denominator is Te IN KELVIN, recovered from the eV number the case
    actually supplies:

        nu_ee = 3.7 n_e[cm^-3] lnLambda / ( (e/k_b) Te[eV] )^1.5

    nu_ea is the electron-neutral elastic (momentum-transfer) collision
    frequency, recovered from the ALREADY-TABULATED electron mobility via
    the standard relation mu_e = e/(me nu_ea) -- NO new cross-section data
    is needed. tau_diff = L_diff^2/D_e (eq. 13), with D_e likewise
    recovered from mu_e via the Einstein relation D_e = mu_e Te[eV], rather
    than looked up as a fourth registered field.

    WHAT THIS RETURNS: epsilon_eff alone, an ENERGY [eV] -- eq. (12)'s
    H_secondary = epsilon_eff * S_fast is NOT formed here. See the class-
    level comment on CoulombHeatingProperty in
    genericPlasmaPropertyTemplates.H for why the S_iz multiplication is
    deliberately done OUTSIDE this closure and outside the evaluator, and
    why its output must never be routed through the same path as
    Pelastic_/Pinelastic_.

    BOUNDED BY CONSTRUCTION, NO IONISATION-DEGREE GATE NEEDED: the ratio
    2 nu_ee / (2 nu_ee + positive terms) lies in [0, 1) for any n_e, Te, so
    epsilon_eff is bounded in [0, epsilon_s) regardless of how n_e or Te
    behave -- it cannot diverge even at a floored, near-vacuum density; it
    simply goes to zero (nu_ee -> 0 as n_e -> 0).

SourceFiles
    (header-only)
\*---------------------------------------------------------------------------*/

#ifndef coulombHeatingClosure_H
#define coulombHeatingClosure_H

#include "plasmaConstants.H"
#include "scalar.H"
#include <algorithm>
#include <cmath>

namespace Foam
{
namespace coulombHeatingClosure
{

//- eq. (14): the electron-electron (Coulomb) collision frequency [s^-1].
//  ne_SI: electron density [m^-3]. Te_eV: electron temperature [eV].
//  coulombLog: the Coulomb logarithm (~10 for typical DC-glow conditions).
inline scalar nuEE
(
    const scalar ne_SI,
    const scalar Te_eV,
    const scalar coulombLog
)
{
    const scalar e  = constant::plasma::eCharge.value();
    const scalar kB = constant::plasma::kappaBoltzmann.value();

    const scalar TeEvFloored = std::max(Te_eV, SMALL);
    const scalar neCC = std::max(ne_SI*1e-6, 0.0); // m^-3 -> cm^-3

    // Te IN KELVIN, recovered from the eV number via e/k_b ~ 11604 K/eV --
    // see the file header for why this is not a bare Te[eV]^1.5 formula.
    const scalar TeK = (e/kB)*TeEvFloored;

    return 3.7*neCC*coulombLog/std::pow(TeK, 1.5);
}

//- eq. (15): epsilon_eff [eV], the order-of-magnitude estimate.
//
//  ne_SI:     electron density [m^-3]
//  meanE_eV:  mean energy [eV] (LMEA's own n_eps/n_e; Te = (2/3) meanE)
//  muE_SI:    electron mobility [m^2/(V.s)]
//  coulombLog:          Coulomb logarithm, Lambda
//  epsilon1_eV:         first inelastic (excitation) threshold [eV]
//  backgroundAtomMass:  neutral gas atom/molecule mass [kg], m_a
//  secondaryEnergyFraction: epsilon_s / epsilon_1 (default 0.5 upstream)
//  diffusionLength:     characteristic confinement length [m], L_diff
//                       (GREAT => no confining wall)
inline scalar epsilonEff
(
    const scalar ne_SI,
    const scalar meanE_eV,
    const scalar muE_SI,
    const scalar coulombLog,
    const scalar epsilon1_eV,
    const scalar backgroundAtomMass,
    const scalar secondaryEnergyFraction,
    const scalar diffusionLength
)
{
    const scalar e  = constant::plasma::eCharge.value();
    const scalar me = constant::plasma::eMass.value();

    // Te[eV] from the mean energy, same relation this model's own T_ field
    // uses (localEnergyEnergyModel.H: "(2/3) eps_bar").
    const scalar TeEv = std::max((2.0/3.0)*meanE_eV, SMALL);

    const scalar nuEe = nuEE(ne_SI, TeEv, coulombLog);

    // nu_ea from the already-tabulated mobility: mu_e = e/(me nu_ea).
    const scalar muE = std::max(muE_SI, SMALL);
    const scalar nuEa = e/(me*muE);

    // tau_diff = L_diff^2/D_e (eq. 13), D_e from the Einstein relation
    // D_e = mu_e Te[eV] rather than a fourth registered-field lookup.
    const scalar De = muE*TeEv;
    const scalar tauDiff = diffusionLength*diffusionLength/std::max(De, SMALL);

    const scalar epsilonS = secondaryEnergyFraction*epsilon1_eV;

    const scalar numerator = 2.0*nuEe;
    const scalar denominator =
        2.0*nuEe + (2.0*me/backgroundAtomMass)*nuEa + 1.0/tauDiff;

    return (numerator/std::max(denominator, SMALL))*epsilonS;
}

} // End namespace coulombHeatingClosure
} // End namespace Foam

#endif

// ************************************************************************* //
