/*---------------------------------------------------------------------------*\
  File: materialLibrary.H
  Part of: SoPLASMA

Class
    Foam::materialLibrary

Description
    Resolves a MATERIAL NAME to its properties, from the shipped library at
    $SoPLASMA_ETC/materials/dielectrics.

    Governing directive G2: anything with literature values becomes a library
    with cited numbers, not a number typed into a case. So a case writes

        material    alumina96;

    in constant/<region>/electricalProperties and both epsilonR and gammaSEE
    follow, with a citation, from one place.

    PRECEDENCE, per G1 (derive, else default, else require):
      1. an explicit numeric value in the case WINS -- it is an override, and
         is reported as one;
      2. otherwise the named material supplies it;
      3. otherwise the region KIND's default applies.

    Naming a material that is not in the library is FATAL and lists what is
    available, rather than falling back to a default nobody chose.

    ONE OWNER for the lookup. epsilonR is resolved by
    electromagneticsModel::readEpsilonR and gammaSEE by dielectricRegion, and
    both come through here: two lookups would be two libraries the moment one
    of them grew a special case.

SourceFiles
    materialLibrary.C

\*---------------------------------------------------------------------------*/

#ifndef materialLibrary_H
#define materialLibrary_H

#include "dictionary.H"
#include "wordList.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class materialLibrary Declaration
\*---------------------------------------------------------------------------*/

class materialLibrary
{
    // Private Member Functions

        //- The library dictionary, loaded once on first use
        static const dictionary& library();


public:

    // Static Member Functions

        //- Whether the library could be found at all
        static bool available();

        //- The materials the library defines
        static wordList names();

        //- The entry for `name`. FATAL if absent, listing what is available.
        //  `context` names who asked, so the error points at the case file
        //  rather than at this class.
        static const dictionary& lookup
        (
            const word& name,
            const string& context
        );

        //- One property of a material [dimensionless]. FATAL if the material
        //  exists but does not define it -- a material that cannot answer for
        //  a property it is being asked for is a library defect, not a case
        //  defect, and should say so.
        static scalar get
        (
            const word& name,
            const word& property,
            const string& context
        );

        //- The library's own file path, for provenance reporting
        static fileName path();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
