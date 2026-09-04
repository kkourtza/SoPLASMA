/*---------------------------------------------------------------------------*\
  File: plasmaSetupBoundaries.C
  Part of: SoPLASMA

Application
    plasmaSetupBoundaries

Description
    Generate the ePotential and surfCharge boundary conditions for every region
    from the case's SEMANTIC boundary description.

    Governing directive G2: the user says WHAT a surface IS, in
    configuration/boundaries; this reads $SoPLASMA_ETC/boundaryRoles to learn
    what each KIND means to each module, and writes the dictionaries.

    THREE THINGS ARE DERIVED AND MUST NOT BE DECLARED BY A CASE:

      * MECHANICAL patches -- empty, wedge, symmetry, cyclic, processor. These
        are read from the MESH. A mechanical constraint is not a physical
        description of a surface, and a case that restated it could disagree
        with the mesh.

      * REGION INTERFACES -- `<this>_to_<other>`, created by splitMeshRegions.
        The topology already says what they are. A gas/dielectric pair gets
        coupledElectricPotential with the surface charge owned by the GAS side
        ONLY (exactly one owner, or it is double-counted); a gas/farField pair
        gets zeroGradient and no charging, because a farField region is
        fictitious air and not a surface at all.

      * The permittivity and gamma behind a named MATERIAL, which come from the
        material library.

    Run from the case directory, AFTER splitMeshRegions:

        plasmaSetupBoundaries              # generate
        plasmaSetupBoundaries -listKinds   # list the kinds and their parameters

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "IFstream.H"
#include "volFields.H"
#include "OSspecific.H"

#include "materialLibrary.H"

#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * Helpers * * * * * * * * * * * * * * * * //

//- The shipped role library, anchored on $SoPLASMA_ETC like the material one:
//  a searched path would make WHICH library answered depend on the working
//  directory, and the point of a library is that two cases agree.
static fileName rolesPath()
{
    const fileName etc(Foam::getEnv("SoPLASMA_ETC"));

    if (etc.empty())
    {
        FatalErrorInFunction
            << "SoPLASMA_ETC is not set, so the boundary-role library cannot"
            << " be found." << nl
            << "    Source the project's etc/bashrc." << nl
            << exit(FatalError);
    }

    return etc/"boundaryRoles";
}


static dictionary readRoles()
{
    const fileName p(rolesPath());

    IFstream is(p);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read the boundary-role library:" << nl
            << "    " << p << nl << exit(FatalError);
    }

    dictionary d(is);
    d.remove("FoamFile");
    return d;
}


//- Whether a patch is a MECHANICAL constraint rather than a physical surface.
static bool mechanicalType(const polyPatch& pp, word& type)
{
    if (isA<emptyPolyPatch>(pp))              { type = "empty";         return true; }
    if (isA<wedgePolyPatch>(pp))              { type = "wedge";         return true; }
    if (isA<symmetryPolyPatch>(pp))           { type = "symmetry";      return true; }
    if (isA<symmetryPlanePolyPatch>(pp))      { type = "symmetryPlane"; return true; }
    if (isA<cyclicPolyPatch>(pp))             { type = "cyclic";        return true; }
    if (isA<processorPolyPatch>(pp))           { type = "processor";     return true; }
    return false;
}


//- The value tokens of an entry, written back verbatim. Used for `$param`
//  substitution, where the parameter may be a whole Function1 such as
//  `table ((0 0) (1e-7 8e3))` -- which cannot be reduced to a scalar and must
//  survive as written.
static string entryValueText(const dictionary& d, const word& key)
{
    // READ THE TOKEN LIST, not the stream.
    //
    // The first version drove `ITstream&` with `is >> t` in a loop and got
    // EMPTY strings: an ITstream returned by lookup() carries its own read
    // position, and the eof/good handling around it is fragile enough that the
    // whole value came back blank -- generating `type ;` for every patch.
    // ITstream IS a tokenList, so indexing it needs no stream state at all.
    const ITstream& is = d.lookup(key);

    OStringStream os;

    forAll(is, i)
    {
        if (i) os << ' ';
        os << is[i];
    }

    return os.str();
}


// * * * * * * * * * * * * * * * * * * Main  * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate ePotential/surfCharge boundary conditions from"
        " configuration/boundaries and the shipped boundary-role library"
    );
    argList::addBoolOption
    (
        "listKinds",
        "list the boundary kinds the library defines, and exit"
    );

    #include "setRootCase.H"

    const dictionary roles(readRoles());

    if (args.found("listKinds"))
    {
        Info<< nl << "Boundary kinds in " << rolesPath() << nl << nl;

        for (const word& k : roles.toc())
        {
            const dictionary& r = roles.subDict(k);

            Info<< "  " << k << nl
                << "      " << r.getOrDefault<string>("description", "") << nl;

            if (r.found("requires"))
            {
                Info<< "      requires: " << r.get<wordList>("requires") << nl;
            }
            if (r.found("optional"))
            {
                Info<< "      optional: " << r.subDict("optional").toc() << nl;
            }
            Info<< nl;
        }

        Info<< "  Region INTERFACES and MECHANICAL patches are DERIVED and must"
               " not be declared." << nl << endl;

        return 0;
    }

    #include "createTime.H"

    // ---- which regions exist, and of what kind ----------------------------
    HashTable<word> kindOfRegion;      // region name -> gas | dielectric | farField
    {
        IOdictionary rp
        (
            IOobject
            (
                "regionProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        HashTable<wordList> regions;
        rp.readEntry("regions", regions);

        forAllConstIters(regions, iter)
        {
            for (const word& r : iter.val())
            {
                kindOfRegion.insert(r, iter.key());
            }
        }
    }

    // ---- the case's semantic description (LAYER 1) ------------------------
    dictionary decl;
    {
        const fileName p(runTime.path()/"configuration"/"boundaries");

        IFstream is(p);

        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot read the case's boundary description:" << nl
                << "    " << p << nl << nl
                << "    This file is LAYER 1: one block per patch saying WHAT"
                   " that surface is." << nl
                << "    Run `plasmaSetupBoundaries -listKinds` to see the"
                   " available kinds." << nl
                << exit(FatalError);
        }

        decl = dictionary(is);
        decl.remove("FoamFile");
    }

    label nWritten = 0;

    forAllConstIters(kindOfRegion, regIter)
    {
        const word& regionName = regIter.key();
        const word& regionKind = regIter.val();

        fvMesh mesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // A dielectric or farField region carries the potential only.
        const bool isGas = (regionKind == "gas");

        IOobject io
        (
            "ePotential",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        );

        mkDir(io.path());

        OFstream os(io.objectPath());
        io.writeHeader(os, volScalarField::typeName);

        os  << "dimensions      [1 2 -3 0 0 -1 0];" << nl << nl
            << "internalField   uniform 0;" << nl << nl
            << "// GENERATED by plasmaSetupBoundaries from"
               " configuration/boundaries" << nl
            << "// and " << rolesPath() << "." << nl
            << "// Edit the SEMANTIC description, not this file." << nl << nl
            << "boundaryField" << nl << "{" << nl;

        forAll(mesh.boundaryMesh(), patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];
            const word& name = pp.name();

            os << "    " << name << nl << "    {" << nl;

            // 1. MECHANICAL -- from the mesh.
            word mech;
            if (mechanicalType(pp, mech))
            {
                os  << "        type            " << mech << ";" << nl
                    << "    }" << nl << nl;
                continue;
            }

            // 2. A REGION INTERFACE -- from the topology.
            bool handled = false;
            forAllConstIters(kindOfRegion, other)
            {
                if (other.key() == regionName) continue;

                if (name != regionName + "_to_" + other.key()) continue;

                handled = true;

                if (other.val() == "farField")
                {
                    // Fictitious air, not a surface: nothing accumulates.
                    os  << "        // DERIVED: interface to the FICTITIOUS"
                           " farField region `" << other.key() << "`." << nl
                        << "        type            zeroGradient;" << nl;
                }
                else
                {
                    // Surface charge has EXACTLY ONE OWNER: the gas side,
                    // because the plasma is what deposits it.
                    os  << "        // DERIVED: meshed interface to `"
                        << other.key() << "`. The surface charge is owned by"
                           " the GAS side only." << nl
                        << "        type            coupledElectricPotential;"
                        << nl
                        << "        value           uniform 0;" << nl
                        << "        surfCharge      "
                        << (isGas ? "surfCharge" : "none") << ";" << nl
                        << "        surfChargeNbr   "
                        << (isGas ? "none" : "surfCharge") << ";" << nl
                        << "        useImplicit     true;" << nl;
                }

                os << "    }" << nl << nl;
                break;
            }
            if (handled) continue;

            // 3. DECLARED by the case.
            if (!decl.found(name))
            {
                FatalErrorInFunction
                    << "Patch `" << name << "` of region `" << regionName
                    << "` is not described." << nl << nl
                    << "    Add a block to configuration/boundaries saying what"
                       " it IS:" << nl << nl
                    << "        " << name << nl
                    << "        {" << nl
                    << "            kind    <one of the kinds below>;" << nl
                    << "        }" << nl << nl
                    << "    Available kinds: " << roles.toc() << nl
                    << "    `plasmaSetupBoundaries -listKinds` describes them"
                       " and lists their parameters." << nl << nl
                    << "    A patch is NOT declared when it is a mechanical"
                       " constraint (empty, wedge," << nl
                    << "    symmetry, processor) or a region interface -- both"
                       " are derived. This one is" << nl
                    << "    neither, so it needs a physical description." << nl
                    << exit(FatalError);
            }

            const dictionary& pd = decl.subDict(name);

            const word kind(pd.get<word>("kind"));

            if (!roles.found(kind))
            {
                FatalErrorInFunction
                    << "Unknown boundary kind `" << kind << "` on patch `"
                    << name << "`." << nl << nl
                    << "    Available: " << roles.toc() << nl
                    << "    Defined in " << rolesPath() << nl
                    << exit(FatalError);
            }

            const dictionary& role = roles.subDict(kind);

            // Required parameters, checked HERE so the error names the patch
            // and the kind rather than surfacing later as a missing key in a
            // generated file nobody wrote.
            if (role.found("requires"))
            {
                for (const word& need : role.get<wordList>("requires"))
                {
                    if (!pd.found(need))
                    {
                        FatalErrorInFunction
                            << "Patch `" << name << "` is `" << kind
                            << "`, which REQUIRES `" << need << "`." << nl
                            << nl
                            << "    " << role.getOrDefault<string>
                                          ("description", "") << nl << nl
                            << "    Add it to that block in"
                               " configuration/boundaries." << nl
                            << exit(FatalError);
                    }
                }
            }

            os  << "        // " << kind << ": "
                << role.getOrDefault<string>("description", "") << nl;

            const dictionary& es = role.subDict("electrostatics");

            for (const word& key : es.toc())
            {
                string val(entryValueText(es, key));

                // "<param>"          -> the patch's own entry, verbatim
                // "<material:prop>"   -> through the material library
                //
                // Quoted-and-bracketed rather than `$param`: OpenFOAM's own
                // dictionary reader resolves `$name` at PARSE time and aborts
                // on one it cannot find, so a `$` placeholder never reaches
                // here. A quoted string is inert.
                if (val.size() > 4 && val.starts_with("\"<") && val.ends_with(">\""))
                {
                    const string ref(val.substr(2, val.size() - 4));

                    const auto colon = ref.find(':');

                    if (colon != std::string::npos)
                    {
                        const word mkey(ref.substr(0, colon));
                        const word prop(ref.substr(colon + 1));

                        const word mat(pd.get<word>(mkey));

                        OStringStream mo;
                        mo << materialLibrary::get
                        (
                            mat, prop,
                            "configuration/boundaries, patch `" + name + "`"
                        );
                        val = mo.str();
                    }
                    else if (ref == "internalField")
                    {
                        val = "uniform 0";
                    }
                    else if (pd.found(ref))
                    {
                        val = entryValueText(pd, ref);
                    }
                    else if (role.found("optional")
                          && role.subDict("optional").found(ref))
                    {
                        val = entryValueText(role.subDict("optional"), ref);
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Kind `" << kind << "` needs `" << ref
                            << "`, which patch `" << name << "` does not give"
                            << " and the library does not default." << nl
                            << exit(FatalError);
                    }
                }

                os  << "        " << key;
                for (label i = key.size(); i < 20; ++i) os << ' ';
                os << val.c_str() << ';' << nl;
            }

            os << "    }" << nl << nl;
        }

        os << "}" << nl << nl;
        io.writeEndDivider(os);

        Info<< "  wrote " << io.objectPath() << endl;
        ++nWritten;
    }

    Info<< nl << "plasmaSetupBoundaries: " << nWritten
        << " ePotential field(s) generated." << nl << endl;

    Info<< "End" << nl << endl;
    return 0;
}

// ************************************************************************* //
