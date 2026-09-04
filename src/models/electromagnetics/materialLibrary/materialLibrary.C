/*---------------------------------------------------------------------------*\
  File: materialLibrary.C
  Part of: SoPLASMA
\*---------------------------------------------------------------------------*/

#include "materialLibrary.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Static Data  * * * * * * * * * * * * * * //

namespace Foam
{
namespace
{
    //- Resolved once: the library is read-only and the same for every region.
    bool libLoaded_ = false;
    dictionary lib_;
    fileName libPath_;
}
}


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

const Foam::dictionary& Foam::materialLibrary::library()
{
    if (libLoaded_) return lib_;

    libLoaded_ = true;

    // ANCHORED ON $SoPLASMA_ETC, which etc/bashrc exports, rather than searched
    // for. A search path would make WHICH library was used depend on the
    // working directory -- and the whole point of a library is that two cases
    // asking for `alumina96` get the same number.
    const fileName etc(Foam::getEnv("SoPLASMA_ETC"));

    if (etc.empty()) return lib_;      // reported by available()/lookup()

    libPath_ = etc/"materials"/"dielectrics";

    IFstream is(libPath_);

    if (!is.good()) return lib_;

    lib_ = dictionary(is);

    // The FoamFile header is not a material.
    lib_.remove("FoamFile");

    return lib_;
}


// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

bool Foam::materialLibrary::available()
{
    return !library().empty();
}


Foam::fileName Foam::materialLibrary::path()
{
    library();
    return libPath_;
}


Foam::wordList Foam::materialLibrary::names()
{
    return library().toc();
}


const Foam::dictionary& Foam::materialLibrary::lookup
(
    const word& name,
    const string& context
)
{
    const dictionary& lib = library();

    if (lib.empty())
    {
        FatalErrorInFunction
            << "A material (`" << name << "`) was named by " << context
            << ", but the material library could not be read." << nl << nl
            << "    Expected at  $SoPLASMA_ETC/materials/dielectrics" << nl
            << "    SoPLASMA_ETC = "
            << (Foam::getEnv("SoPLASMA_ETC").empty()
                    ? "(NOT SET -- source the project's etc/bashrc)"
                    : Foam::getEnv("SoPLASMA_ETC"))
            << nl << nl
            << "    State the numbers directly in the case instead, or fix the"
               " environment. A named" << nl
            << "    material must not silently become a default nobody chose."
            << nl
            << exit(FatalError);
    }

    if (!lib.found(name))
    {
        FatalErrorInFunction
            << "Unknown material `" << name << "`, named by " << context
            << "." << nl << nl
            << "    Available in " << libPath_ << ":" << nl
            << "        " << lib.toc() << nl << nl
            << "    Add it to that file -- WITH a `reference` -- or state"
               " epsilonR and gammaSEE" << nl
            << "    directly in the case. A misspelt material must not fall"
               " back to a default." << nl
            << exit(FatalError);
    }

    return lib.subDict(name);
}


Foam::scalar Foam::materialLibrary::get
(
    const word& name,
    const word& property,
    const string& context
)
{
    const dictionary& m = lookup(name, context);

    if (!m.found(property))
    {
        // A LIBRARY defect, not a case defect, and the message says so: the
        // user did nothing wrong by naming a material that exists.
        FatalErrorInFunction
            << "Material `" << name << "` does not define `" << property
            << "`, which " << context << " needs." << nl << nl
            << "    This is a defect in the LIBRARY, not in the case:" << nl
            << "        " << libPath_ << nl << nl
            << "    Add `" << property << "` to that entry, with a reference."
            << nl
            << exit(FatalError);
    }

    return m.get<scalar>(property);
}


// ************************************************************************* //
