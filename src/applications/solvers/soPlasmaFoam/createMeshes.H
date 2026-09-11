/*---------------------------------------------------------------------------*\
  File: createMeshes.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createMeshes.H in OpenFOAM v2412.
      Creates and initializes the region meshes.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

regionProperties rp(runTime, IOobjectOption::READ_IF_PRESENT);

// Reject any region kind this solver does not implement, BEFORE anything is
// derived from the kinds. Shared with multiRegionElectrostaticFoam and
// implemented in electromagneticsModel, which is the class that would otherwise
// degrade silently to singleRegionPoisson -- see the reasoning there.
electromagneticsModel::validateRegionKinds(rp.sortedToc());

wordList gasNames;
if (rp.found("gas"))
{
    gasNames = rp["gas"];
}

// ONE `gas`, ANY NUMBER OF `dielectric`.
//
// This is not a limitation on the geometry. A stack such as
// electrode|dielectric|gas|dielectric|gas|dielectric|electrode is set up as ONE
// gas region whose cellZone happens to be geometrically disconnected -- an
// OpenFOAM mesh region may consist of several disjoint parts -- plus three
// dielectric regions, each with its own constant/<region>/electricalProperties.
// Every gas gap then shares one species set, one chemistry and one Boltzmann
// solution, which is what you want when the gaps hold the same gas.
//
// What is genuinely NOT supported is two gas regions with DIFFERENT gas
// mixtures (air in one gap, CO2 in the other): that needs a second species
// list, a second mechanism and a second EEDF, none of which are per-region
// today.
if (gasNames.size() > 1)
{
    FatalErrorInFunction
        << "Only one `gas` region is supported by this solver, but"
        << " constant/regionProperties names " << gasNames.size()
        << ": " << gasNames << nl << nl
        << "    If these are SEPARATE GAPS OF THE SAME GAS, you do not need"
        << " more than one region: put every gas cell in a single cellZone."
        << " It may be geometrically disconnected -- the dielectrics between"
        << " the gaps simply are not part of it -- and all of it is then solved"
        << " with one species set and one chemistry." << nl << nl
        << "    If they hold DIFFERENT GASES, that is not supported: the"
        << " species list, mechanism and EEDF are global, not per-region." << nl
        << exit(FatalError);
}

word gasRegionName;

if (gasNames.empty())
{
    gasRegionName = polyMesh::defaultRegion;
    Info<< "No gas region specified in constant/regionProperties." << nl
        << "Falling back to default mesh: " << gasRegionName << nl << endl;
}
else
{
    gasRegionName = gasNames[0];
}

wordList dielectricNames;
if (rp.found("dielectric"))
{
    dielectricNames = rp["dielectric"];
}

// const label nTotalRegions = 
//     (gasNames.empty() ? 1 : gasNames.size()) + dielectricNames.size();

// FICTITIOUS FAR-FIELD REGIONS join the dielectric list.
//
// A `farField` region is an air/vacuum region that solves ONLY Poisson, exactly
// as a dielectric does -- so it needs no new machinery here, just a name. What
// the name buys is that `epsilonR` DEFAULTS to 1.0 instead of being demanded
// (see electromagneticsModel::readEpsilonR) and that the intent is declared, so
// the solver can check the interface rather than hope.
//
// WHY IT EXISTS. Poisson is long-range, so a correct potential distribution can
// need a domain far larger than the plasma. Species transport, chemistry and
// the electron energy equation are not long-range, and photoionization is
// exponentially screened at ~1 mm in atmospheric air. Wrapping a small gas
// region in a large farField region therefore buys the electrostatics a big
// domain at almost no cost -- none of the expensive equations are solved there.
//
// THE ONE PHYSICAL CONSTRAINT: a farField region carries NO SPACE CHARGE, since
// Poisson in a Poisson-only region has no rho on the right-hand side. Put the
// interface where the charge density is negligible, and check that it stays so.
wordList farFieldNames;
if (rp.found("farField"))
{
    farFieldNames = rp["farField"];

    Info<< "Fictitious far-field region(s): " << farFieldNames << nl
        << "    Poisson only, epsilonR defaults to 1.0 (air/vacuum)." << nl
        << "    No species, chemistry, electron energy or photoionization is"
        << " solved there," << nl
        << "    and no space charge exists there -- keep the interface where"
        << " rho is negligible." << nl << endl;

    dielectricNames.append(farFieldNames);
}

#include "createGasMesh.H"
#include "createDielectricMeshes.H"
