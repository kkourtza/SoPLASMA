/*---------------------------------------------------------------------------*\
  File: boundaryRoleLibrary.C
  Part of: SoPLASMA
\*---------------------------------------------------------------------------*/

#include "boundaryRoleLibrary.H"
#include "Time.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Static Data  * * * * * * * * * * * * * * //

namespace Foam
{
namespace
{
    bool loaded_ = false;
    dictionary lib_;
    fileName path_;
}
}


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

const Foam::dictionary& Foam::boundaryRoleLibrary::library()
{
    if (loaded_) return lib_;

    loaded_ = true;

    // ANCHORED on $SoPLASMA_ETC, not searched for: a search path would make
    // WHICH library answered depend on the working directory, and the point of
    // a library is that two cases naming the same kind agree.
    const fileName etc(Foam::getEnv("SoPLASMA_ETC"));

    if (etc.empty())
    {
        FatalErrorInFunction
            << "SoPLASMA_ETC is not set, so the boundary-role library cannot be"
            << " found." << nl
            << "    Source the project's etc/bashrc." << nl
            << exit(FatalError);
    }

    path_ = etc/"boundaryRoles";

    IFstream is(path_);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read the boundary-role library:" << nl
            << "    " << path_ << nl << exit(FatalError);
    }

    lib_ = dictionary(is);
    lib_.remove("FoamFile");

    return lib_;
}


// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

Foam::fileName Foam::boundaryRoleLibrary::path()
{
    library();
    return path_;
}


Foam::wordList Foam::boundaryRoleLibrary::kinds()
{
    return library().toc();
}


const Foam::dictionary& Foam::boundaryRoleLibrary::lookup
(
    const word& kind,
    const string& context
)
{
    const dictionary& lib = library();

    if (!lib.found(kind))
    {
        FatalErrorInFunction
            << "Unknown boundary kind `" << kind << "`, named by " << context
            << "." << nl << nl
            << "    Available: " << lib.toc() << nl
            << "    Defined in " << path_ << nl
            << "    `plasmaSetupBoundaries -listKinds` describes them." << nl
            << exit(FatalError);
    }

    return lib.subDict(kind);
}


Foam::word Foam::boundaryRoleLibrary::surfaceClass
(
    const word& kind,
    const string& context
)
{
    const dictionary& r = lookup(kind, context);

    if (!r.found("plasmaTransport"))
    {
        FatalErrorInFunction
            << "Boundary kind `" << kind << "` does not say what SPECIES do on"
            << " it." << nl << nl
            << "    This is a defect in the LIBRARY, not in the case:" << nl
            << "        " << path_ << nl << nl
            << "    Add a `plasmaTransport { surface <class>; }` block, where"
               " <class> is one of" << nl
            << "    mechanical | open | conductor | chargingSurface." << nl
            << exit(FatalError);
    }

    const word s(r.subDict("plasmaTransport").get<word>("surface"));

    if
    (
        s != "mechanical" && s != "open"
     && s != "conductor" && s != "chargingSurface"
    )
    {
        FatalErrorInFunction
            << "Boundary kind `" << kind << "` declares surface class `" << s
            << "`, which is not one of" << nl
            << "    mechanical | open | conductor | chargingSurface." << nl
            << "    In " << path_ << nl
            << exit(FatalError);
    }

    return s;
}


Foam::dictionary Foam::boundaryRoleLibrary::caseDeclaration
(
    const Time& runTime
)
{
    const fileName p(runTime.path()/"configuration"/"boundaries");

    IFstream is(p);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read the case's boundary description:" << nl
            << "    " << p << nl << nl
            << "    This file is LAYER 1: one block per patch saying WHAT that"
               " surface is." << nl
            << "    Every boundary condition in the case is derived from it,"
               " so there is nothing" << nl
            << "    to fall back on." << nl << nl
            << "    `plasmaSetupBoundaries -listKinds` lists the available"
               " kinds." << nl
            << exit(FatalError);
    }

    dictionary d(is);
    d.remove("FoamFile");
    return d;
}


// ************************************************************************* //
