/*---------------------------------------------------------------------------*\
  File: floatingPotentialInterface.H
  Part of: SoPLASMA

Class
    Foam::floatingPotentialInterface

Description
    The narrow contract between Foam::floatingElectrode, which owns the
    charge-conservation constraint, and the boundary condition that carries the
    resulting equipotential value.

    WHY AN INTERFACE RATHER THAN THE CONCRETE TYPE.
    floatingElectrode lives in libplasmaElectromagnetics; the boundary
    condition lives in libplasmaBcs, which DEPENDS on that library because the
    coupled-potential condition reads the model's permittivity. Including the
    condition's header from the model side would therefore be circular.

    So the model side declares what it needs, the condition inherits it, and
    the model reaches it by dynamic_cast. Three functions is the whole surface:
    the charge in, the potential out, and the potential back in.

\*---------------------------------------------------------------------------*/

#ifndef floatingPotentialInterface_H
#define floatingPotentialInterface_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class floatingPotentialInterface Declaration
\*---------------------------------------------------------------------------*/

class floatingPotentialInterface
{
public:

    //- Destructor
    virtual ~floatingPotentialInterface() = default;

    //- The initial total charge Q0 [C], stated by the case
    virtual scalar initialCharge() const = 0;

    //- The floating potential [V] the condition currently carries
    virtual scalar floatingPotential() const = 0;

    //- Set the equipotential value [V]
    virtual void setFloatingPotential(const scalar Vf) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
