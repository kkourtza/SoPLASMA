/*---------------------------------------------------------------------------*\
  File: emissionModel.H
  Part of: SoPLASMA

Class
    Foam::emissionModel

Description
    Abstract base for a mechanism that EMITS ELECTRONS FROM A WALL.

    Every such mechanism answers one question -- how many electrons per unit
    area per second leave this surface -- and differs only in what drives it:

        ionInducedSEE        an arriving ion, neutralised at the surface
        electronInducedSEE   an arriving electron, above the yield threshold
        fieldEmission        the local field, tunnelling through the barrier
        thermionic           the wall temperature, over the barrier
        photoEmission        an arriving photon

    So they are ONE runtime-selectable family summed into the electron wall
    flux, not a growing set of flags on a boundary condition. A mechanism
    contributes its flux exactly the way secondary emission already does:
    added to the mixed condition's refGradient as Gamma/D.

    EVERY ACTIVE MECHANISM REPORTS ITSELF AT START-UP, and that is not
    decoration. Emission decides whether a discharge is SELF-SUSTAINING rather
    than a single avalanche, so a mechanism that is silently on -- or silently
    off -- produces a plausible, healthy-looking, wrong discharge. report() is
    part of the interface for that reason.

SourceFiles
    emissionModel.C

\*---------------------------------------------------------------------------*/

#ifndef emissionModel_H
#define emissionModel_H

#include "fvPatch.H"
#include "scalarField.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class emissionModel Declaration
\*---------------------------------------------------------------------------*/

class emissionModel
{
protected:

    //- The patch this mechanism acts on
    const fvPatch& patch_;

    //- Its own dictionary
    const dictionary dict_;


    // Protected Member Functions

        //- The field magnitude at the patch [V/m], from `Emag`.
        //
        //  FATAL when absent rather than defaulted to zero: a field-driven
        //  mechanism with no field is not a mechanism contributing nothing, it
        //  is a mechanism that cannot be evaluated, and silently returning
        //  zero would switch it off invisibly.
        tmp<scalarField> patchField() const;

        //- The wall temperature at the patch [K], from `T_gas`.
        //
        //  The BOUNDARY value of the gas temperature IS the wall temperature:
        //  set that field's boundary condition to run a heated electrode. FATAL
        //  when T_gas does not exist, which is the case when the gas energy
        //  equation is not solved -- a thermionic cathode with no temperature
        //  field is a configuration error, not a zero.
        tmp<scalarField> patchTemperature() const;

        //- The work function [eV] of this surface, from the material library
        //  when a `material` is named, else from an explicit `workFunction`.
        scalar workFunction() const;


public:

    //- Runtime type information
    TypeName("emissionModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        emissionModel,
        dictionary,
        (
            const fvPatch& p,
            const dictionary& dict
        ),
        (p, dict)
    );


    // Constructors

        emissionModel(const fvPatch& p, const dictionary& dict);

        //- Selector
        static autoPtr<emissionModel> New
        (
            const word& modelType,
            const fvPatch& p,
            const dictionary& dict
        );


    //- Destructor
    virtual ~emissionModel() = default;


    // Member Functions

        //- Electron flux leaving the surface [1/m2/s], per face. Always
        //  POSITIVE or zero: it is a magnitude, and the caller decides the
        //  sign convention of the flux it is added to.
        virtual tmp<scalarField> emittedFlux() const = 0;

        //- One line naming the mechanism and the numbers it will use, printed
        //  once at start-up. See the class description for why this is part of
        //  the interface.
        virtual void report(Ostream& os) const = 0;

        //- The patch
        inline const fvPatch& patch() const { return patch_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
