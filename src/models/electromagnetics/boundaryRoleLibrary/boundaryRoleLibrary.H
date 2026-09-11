/*---------------------------------------------------------------------------*\
  File: boundaryRoleLibrary.H
  Part of: SoPLASMA

Class
    Foam::boundaryRoleLibrary

Description
    Resolves a boundary KIND to what it means, from the shipped library at
    $SoPLASMA_ETC/boundaryRoles, and reads the case's own semantic description
    from configuration/boundaries.

    G2: the user says WHAT a surface IS; each module decides HOW it treats it.
    G3: the KIND names are the shared vocabulary, and a new module adds a
        sub-dictionary to the library rather than a name of its own.

    ONE OWNER, deliberately. Two utilities need this -- plasmaSetupBoundaries
    for the potential's conditions and plasmaCreateSpeciesFields for the
    species' -- and a second reader would be a second library the moment one of
    them grew a special case. That is not hypothetical: before this class
    existed, the species generator INFERRED the kind by reading the
    ELECTROSTATICS boundary condition type and mapping it back to a semantic
    class, which is deriving one module's answer from another's rather than both
    from the description.

SourceFiles
    boundaryRoleLibrary.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryRoleLibrary_H
#define boundaryRoleLibrary_H

#include "dictionary.H"
#include "wordList.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;

/*---------------------------------------------------------------------------*\
                     Class boundaryRoleLibrary Declaration
\*---------------------------------------------------------------------------*/

class boundaryRoleLibrary
{
    //- The library, loaded once on first use
    static const dictionary& library();


public:

    // Static Member Functions

        //- The library's file path, for provenance and error messages
        static fileName path();

        //- The kinds the library defines
        static wordList kinds();

        //- The entry for `kind`. FATAL if absent, listing what is available.
        static const dictionary& lookup
        (
            const word& kind,
            const string& context
        );

        //- What SPECIES do on a surface of this kind:
        //      mechanical | open | conductor | chargingSurface
        //  See the library's own header for what each means. FATAL if the kind
        //  does not declare a plasmaTransport/surface -- a kind that cannot say
        //  what species do on it is a library defect, not a case defect.
        static word surfaceClass
        (
            const word& kind,
            const string& context
        );

        //- The case's semantic boundary description, configuration/boundaries.
        //  FATAL if absent, because there is nothing to derive from.
        static dictionary caseDeclaration(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
