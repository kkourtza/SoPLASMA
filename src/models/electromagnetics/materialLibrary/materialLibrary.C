/*---------------------------------------------------------------------------*\
  File: materialLibrary.C
  Part of: SoPLASMA
\*---------------------------------------------------------------------------*/

#include "materialLibrary.H"
#include "HashTable.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Static Data  * * * * * * * * * * * * * * //

namespace Foam
{
namespace
{
    //- Resolved once PER LIBRARY: they are read-only and the same everywhere.
    HashTable<dictionary> libs_;
    HashTable<fileName> libPaths_;
}
}


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

const Foam::dictionary& Foam::materialLibrary::library(const word& lib)
{
    auto iter = libs_.find(lib);
    if (iter.good()) return iter.val();

    // ANCHORED ON $SoPLASMA_ETC, which etc/bashrc exports, rather than searched
    // for. A search path would make WHICH library was used depend on the
    // working directory -- and the whole point of a library is that two cases
    // asking for `alumina96` get the same number.
    const fileName etc(Foam::getEnv("SoPLASMA_ETC"));

    // Cache an empty dictionary on failure too, so a missing library is
    // reported once by lookup() rather than re-read on every call.
    libs_.insert(lib, dictionary());
    libPaths_.insert(lib, fileName::null);

    if (etc.empty()) return libs_[lib];   // reported by available()/lookup()

    const fileName p(etc/"materials"/lib);

    libPaths_.set(lib, p);

    IFstream is(p);

    if (!is.good()) return libs_[lib];

    dictionary d(is);

    // The FoamFile header is not an entry.
    d.remove("FoamFile");

    libs_.set(lib, d);

    return libs_[lib];
}


// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

bool Foam::materialLibrary::available(const word& lib)
{
    return !library(lib).empty();
}


Foam::fileName Foam::materialLibrary::path(const word& lib)
{
    library(lib);
    return libPaths_[lib];
}


Foam::wordList Foam::materialLibrary::names(const word& lib)
{
    return library(lib).toc();
}


const Foam::dictionary& Foam::materialLibrary::lookup
(
    const word& name,
    const string& context,
    const word& libName
)
{
    const dictionary& lib = library(libName);
    const fileName& libPath_ = libPaths_[libName];

    if (lib.empty())
    {
        FatalErrorInFunction
            << "A material (`" << name << "`) was named by " << context
            << ", but the material library could not be read." << nl << nl
            << "    Expected at  $SoPLASMA_ETC/materials/" << libName << nl
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
    const string& context,
    const word& libName
)
{
    const dictionary& m = lookup(name, context, libName);
    const fileName& libPath_ = libPaths_[libName];

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
