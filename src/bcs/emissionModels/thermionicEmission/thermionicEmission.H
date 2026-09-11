/*---------------------------------------------------------------------------*\
  File: thermionicEmission.H
  Part of: SoPLASMA

Class
    Foam::emissionModels::thermionicEmission

Description
    Richardson-Dushman thermionic emission: electrons with enough thermal
    energy pass OVER the surface barrier.

        J = lambda_R A0 T^2 exp( - phi_eff / (k_B T) )        [A/m^2]

        A0 = 4 pi m_e e k_B^2 / h^3 = 1.201735e6 A m^-2 K^-2

    A0 is an exact combination of fundamental constants. lambda_R is a
    material-dependent correction of order 0.5 that absorbs quantum reflection
    at the surface and the true density of states; it is NOT derivable from
    anything the solver knows, so it defaults to 1 (the free-electron value)
    and is reported as such. For tungsten it is ~0.5, so a default of 1
    overestimates by 2x -- which is small against the exponential and is
    stated rather than hidden.

    THE SCHOTTKY EFFECT couples this to the field. An applied field lowers the
    barrier by

        dphi = sqrt( e^3 F / (4 pi eps0) ) = 3.794686e-5 sqrt(F[V/m])  [eV]

    so phi_eff = phi - dphi. At F = 1e7 V/m that is 0.12 eV, which at 2000 K
    (k_B T = 0.172 eV) multiplies the current by exp(0.12/0.172) = 2.0. NOT a
    correction to neglect at a discharge electrode. `schottky` (default true)
    applies it; it requires `Emag`, and is the one coupling between this model
    and the field.

    WHERE THE TEMPERATURE COMES FROM. The BOUNDARY VALUE of `T_gas` at this
    patch IS the wall temperature. So a heated electrode is configured by
    setting T_gas's boundary condition on that patch -- a fixedValue for a
    controlled temperature, or a gradient/mixed condition for conduction into
    the wall. Nothing extra is needed here, and nothing is invented: if T_gas
    does not exist, because the gas energy equation is not being solved, the
    model refuses to run rather than assuming a temperature.

SourceFiles
    thermionicEmission.C

\*---------------------------------------------------------------------------*/

#ifndef thermionicEmission_H
#define thermionicEmission_H

#include "emissionModel.H"
#include "Switch.H"

namespace Foam
{
namespace emissionModels
{

/*---------------------------------------------------------------------------*\
                      Class thermionicEmission Declaration
\*---------------------------------------------------------------------------*/

class thermionicEmission
:
    public emissionModel
{
    //- Richardson correction lambda_R [-]
    const scalar lambdaR_;

    //- Whether to apply the Schottky barrier lowering
    const Switch schottky_;

    //- Work function [eV], resolved once
    const scalar phi_;


public:

    TypeName("thermionicEmission");

    thermionicEmission(const fvPatch& p, const dictionary& dict);

    virtual ~thermionicEmission() = default;

    virtual tmp<scalarField> emittedFlux() const;

    virtual void report(Ostream& os) const;
};


} // End namespace emissionModels
} // End namespace Foam

#endif

// ************************************************************************* //
