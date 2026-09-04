/*---------------------------------------------------------------------------*\
  File: fieldEmission.H
  Part of: SoPLASMA

Class
    Foam::emissionModels::fieldEmission

Description
    Fowler-Nordheim field emission: electrons TUNNEL through the surface
    barrier when the local field is high enough.

    In Forbes' elementary form (Forbes & Deane, Proc. R. Soc. A 463 (2007)
    2907; Forbes, Appl. Phys. Lett. 89 (2006) 113122),

        J = (a / phi) F^2 exp( - v_F b phi^{3/2} / F )      [A/m^2]

        a = e^3 / (8 pi h)         = 1.541434e-6  A eV V^-2
        b = (8 pi / 3) sqrt(2m)/(e h) = 6.830890e9 eV^-3/2 V m^-1

    with phi the work function [eV] and F the barrier field [V/m]. Both
    constants are exact combinations of fundamental constants, not fits.

    THE BARRIER CORRECTION v_F MATTERS BY ORDERS OF MAGNITUDE. The elementary
    form takes a triangular barrier, v_F = 1. The real barrier is lowered by the
    image charge of the escaping electron (Schottky-Nordheim), which for the
    standard approximation gives

        f = (e^3 / 4 pi eps0) F / phi^2,   v_F ~ 1 - f + (f/6) ln f

    and RAISES the current, by a factor of ~100 at f = 0.3. `barrier
    schottkyNordheim` (the default) applies it; `barrier elementary` does not.
    The default is the physical one, and the switch exists so the difference can
    be seen rather than assumed.

    FIELD ENHANCEMENT. Emission is exponential in F, so surface roughness
    dominates it: a protrusion of aspect ratio ~100 raises the local field by
    ~100 and the current by many orders. `fieldEnhancement` (beta, default 1)
    multiplies the field the model sees. THIS IS THE LEAST KNOWN NUMBER IN THE
    MODEL, and it is not derivable from anything the solver knows -- it is a
    property of a real surface's microscopic topography. A case that needs field
    emission quantitatively must fit beta, and should say so.

SourceFiles
    fieldEmission.C

\*---------------------------------------------------------------------------*/

#ifndef fieldEmission_H
#define fieldEmission_H

#include "emissionModel.H"

namespace Foam
{
namespace emissionModels
{

/*---------------------------------------------------------------------------*\
                         Class fieldEmission Declaration
\*---------------------------------------------------------------------------*/

class fieldEmission
:
    public emissionModel
{
    //- Field enhancement factor beta [-]
    const scalar beta_;

    //- Barrier treatment: elementary | schottkyNordheim
    const word barrier_;

    //- Work function [eV], resolved once
    const scalar phi_;


public:

    TypeName("fieldEmission");

    fieldEmission(const fvPatch& p, const dictionary& dict);

    virtual ~fieldEmission() = default;

    virtual tmp<scalarField> emittedFlux() const;

    virtual void report(Ostream& os) const;
};


} // End namespace emissionModels
} // End namespace Foam

#endif

// ************************************************************************* //
