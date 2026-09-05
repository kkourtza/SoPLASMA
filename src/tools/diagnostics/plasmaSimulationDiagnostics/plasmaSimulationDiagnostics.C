/*---------------------------------------------------------------------------*\
  File: plasmaSimulationDiagnostics.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaSimulationDiagnostics.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaSimulationDiagnostics.H"
#include "plasmaTransport.H"
#include "plasmaConstants.H"
#include "plasmaRateTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaSimulationDiagnostics::plasmaSimulationDiagnostics
(
    const Time& runTime,
    const plasmaTransport& transport
)
:
    runTime_(runTime),
    transport_(transport),
    dict_(dictionary::null),
    printSpecies_(true),
    printElectromagnetics_(true),
    printLocalityValidity_(true),
    localityMargin_(10.0),
    localityReportFraction_(0.1),
    localityFieldFraction_(0.1),
    reducedEOld_(nullptr),
    reducedEOldTime_(-1),
    lfaCriterionUnavailableReported_(false),
    localityWarmup_(5),
    localityReports_(0),
    lastPctLFA_(-1),
    lastPctLEA_(-1),
    meanEofEN_(nullptr),
    meanEofENTried_(false)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// OUT OF LINE so autoPtr<plasmaRateTable> is destroyed where that type is
// complete. Keeping it inline in the header would force plasmaRateTable.H on
// every consumer of the tools library.
plasmaSimulationDiagnostics::~plasmaSimulationDiagnostics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void plasmaSimulationDiagnostics::read()
{
    IOdictionary plasmaDict
        (
            IOobject
            (
                "plasmaSimulationControls",
                runTime_.system(),
                runTime_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

    if (plasmaDict.found("plasmaSimulationDiagnostics"))
    {
        dict_ = plasmaDict.subDict("plasmaSimulationDiagnostics");

        printSpecies_ =
            dict_.lookupOrDefault<Switch>("printSpecies", true);

        printElectromagnetics_ =
            dict_.lookupOrDefault<Switch>("printElectromagnetics", true);

        printLocalityValidity_ =
            dict_.lookupOrDefault<Switch>("printLocalityValidity", true);

        // What ">>" means. Dias & Guerra state the conditions as strong
        // inequalities without a number; 10 is this solver's reading of that,
        // stated as a setting rather than buried as a literal.
        localityMargin_ =
            dict_.lookupOrDefault<scalar>("localityMargin", 10.0);

        localityReportFraction_ =
            dict_.lookupOrDefault<scalar>("localityReportFraction", 0.1);

        localityFieldFraction_ =
            dict_.lookupOrDefault<scalar>("localityFieldFraction", 0.1);

        localityWarmup_ =
            dict_.lookupOrDefault<label>("localityWarmupReports", 5);
    }
}

void plasmaSimulationDiagnostics::reportLocalityValidity()
{
    const fvMesh& mesh = transport_.species().mesh();

    // reducedE and mu_e are registered under BOTH closures. Without them there
    // is nothing to say, so say nothing.
    if
    (
        !mesh.foundObject<volScalarField>("reducedE")
     || !mesh.foundObject<volScalarField>("mu_e")
    )
    {
        return;
    }

    const volScalarField& EN = mesh.lookupObject<volScalarField>("reducedE");
    const volScalarField& mue = mesh.lookupObject<volScalarField>("mu_e");

    // d(E/N)/dt against the PREVIOUS REPORT, not the previous step: the
    // advisory is a report-interval quantity and differencing against a step
    // the caller may not have taken would be a different number than the one
    // named.
    if (!reducedEOld_.valid())
    {
        reducedEOld_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "reducedE0_locality", runTime_.timeName(), mesh,
                    IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
                ),
                EN
            )
        );
        reducedEOldTime_ = runTime_.value();
        return;                        // no interval yet
    }

    const scalar dt = runTime_.value() - reducedEOldTime_;
    if (dt <= VSMALL) return;

    // nu_eps needs a mean energy, and only LMEA transports one. Under LFA the
    // LEA criterion is still evaluated; the LFA criterion is not, and the
    // advisory says so ONCE rather than pretending to a number it cannot form.
    const bool haveMeanE = mesh.foundObject<volScalarField>("meanE");

    // Under LFA there is no meanE field; fall back to the sweep's own
    // E/N -> <eps> map, which is what the LFA assumes by definition.
    if (!haveMeanE && !meanEofENTried_)
    {
        meanEofENTried_ = true;
        const fileName path =
            transport_.tableDir()/"meanEnergy_vs_reducedE";

        if (isFile(path))
        {
            meanEofEN_.reset
            (
                new plasmaRateTable(path, plasmaRateTable::bhClamp)
            );
        }
    }
    const bool haveMeanEofEN = meanEofEN_.valid();
    const bool havePower =
        mesh.foundObject<volScalarField>("PelasticN")
     && mesh.foundObject<volScalarField>("PgasN")
     && mesh.foundObject<volScalarField>("PvibN");

    const scalarField& en   = EN.primitiveField();
    const scalarField& en0  = reducedEOld_().primitiveField();
    const scalarField& mu   = mue.primitiveField();

    const scalar eOverMe =
        constant::plasma::eCharge.value()/constant::plasma::eMass.value();

    // The domain's own peak field sets the floor: a cell at 1% of the peak is
    // not where the discharge is, and its relative field rate is noise.
    scalar enPeak = 0;
    forAll(en, c) enPeak = max(enPeak, en[c]);
    reduce(enPeak, maxOp<scalar>());
    const scalar enFloor = max(localityFieldFraction_*enPeak, 1.0e-24);

    label nActive = 0, nFailLEA = 0, nFailLFA = 0;
    scalar worstLEA = GREAT, worstLFA = GREAT;

    forAll(en, c)
    {
        // A cell with no field has no forcing to be slow against. Excluded
        // rather than counted as passing, so the fraction is a fraction of the
        // cells the criterion actually applies to.
        //
        // THE FLOOR IS PHYSICAL, NOT `SMALL`. reducedE is SI (V m^2), so a
        // discharge sits at 1e-22 to 1e-19 -- entirely BELOW OpenFOAM's
        // SMALL = 1e-15. Gating on SMALL rejected every cell in the domain and
        // the whole advisory was silently dead: it entered, found every field
        // it needed, and reported nothing. Measured 2026-09-05, and it took a
        // margin of 1e12 (which must fire) to tell "passing" from "dead".
        //
        // The floor is now RELATIVE to the domain's own peak field, because
        // nu_EN is a relative rate and diverges wherever E/N is near zero.
        // See localityFieldFraction_ for the calibration.
        if (en[c] <= enFloor) continue;
        ++nActive;

        const scalar nuEN = mag(en[c] - en0[c])/(dt*en[c]);
        if (nuEN <= VSMALL) continue;          // field steady here: both valid

        // eq. (5): momentum transfer against the forcing rate.
        if (mu[c] > VSMALL)
        {
            const scalar nu_m = eOverMe/mu[c];
            const scalar r = nu_m/nuEN;
            worstLEA = min(worstLEA, r);
            if (r < localityMargin_) ++nFailLEA;
        }

        // eq. (2): energy relaxation against the forcing rate.
        if ((haveMeanE || haveMeanEofEN) && havePower)
        {
            const scalar eps =
                haveMeanE
              ? mesh.lookupObject<volScalarField>("meanE")[c]
              : meanEofEN_().value(en[c]);

            if (eps > SMALL)
            {
                // N from the definition of the reduced field, so it cannot
                // disagree with the field the criterion is built on.
                const scalar Emag =
                    mesh.foundObject<volVectorField>("E")
                  ? mag(mesh.lookupObject<volVectorField>("E")[c]) : 0.0;
                const scalar N = (Emag > SMALL) ? Emag/en[c] : 0.0;

                const scalar P =
                    ( mesh.lookupObject<volScalarField>("PelasticN")[c]
                    + mesh.lookupObject<volScalarField>("PgasN")[c]
                    + mesh.lookupObject<volScalarField>("PvibN")[c] )*N;

                if (P > VSMALL)
                {
                    const scalar r = (P/eps)/nuEN;
                    worstLFA = min(worstLFA, r);
                    if (r < localityMargin_) ++nFailLFA;
                }
            }
        }
    }

    reduce(nActive,  sumOp<label>());
    reduce(nFailLEA, sumOp<label>());
    reduce(nFailLFA, sumOp<label>());
    reduce(worstLEA, minOp<scalar>());
    reduce(worstLFA, minOp<scalar>());

    reducedEOld_() = EN;
    reducedEOldTime_ = runTime_.value();

    if (nActive == 0) return;

    const scalar fLEA = scalar(nFailLEA)/scalar(nActive);
    const scalar fLFA = scalar(nFailLFA)/scalar(nActive);

    const word closure
    (
        haveMeanE ? word("LMEA") : word("LFA")
    );

    // Warm-up: the field establishing from zero is not a locality failure.
    if (++localityReports_ <= localityWarmup_) return;

    const label pctLFA = label(100.0*fLFA + 0.5);
    const label pctLEA = label(100.0*fLEA + 0.5);

    const bool changed = (pctLFA != lastPctLFA_) || (pctLEA != lastPctLEA_);
    lastPctLFA_ = pctLFA;
    lastPctLEA_ = pctLEA;

    if
    (
        changed
     && (fLFA > localityReportFraction_ || fLEA > localityReportFraction_)
    )
    {
        Info<< "  TIME-LOCALITY: the field is changing faster than the"
            << " electrons can follow." << nl
            << "    Dias & Guerra 2025 eqs (2), (5); \">>\" taken as a ratio"
            << " above " << localityMargin_ << "," << nl
            << "    judged only where E/N exceeds "
            << 100.0*localityFieldFraction_ << "% of the domain peak." << nl;

        if (fLFA > localityReportFraction_)
        {
            Info<< "    LFA  invalid in " << pctLFA << "% of "
                << nActive << " cells (worst nu_eps/nu_EN = " << worstLFA
                << ")." << nl;
        }
        if (fLEA > localityReportFraction_)
        {
            Info<< "    LMEA invalid in " << pctLEA << "% of "
                << nActive << " cells (worst nu_m/nu_EN = " << worstLEA
                << ")." << nl
                << "      This is the WEAKER condition -- failing it means"
                << " NEITHER closure applies here, and only a kinetic or"
                << nl
                << "      time-dependent treatment does." << nl;
        }

        Info<< "    This case runs " << closure << "." << endl;
    }

    if (!haveMeanE && !haveMeanEofEN && !lfaCriterionUnavailableReported_)
    {
        lfaCriterionUnavailableReported_ = true;
        Info<< "  TIME-LOCALITY: the LFA criterion (nu_eps/nu_EN) cannot be"
            << " evaluated -- this case transports no" << nl
            << "    mean energy and `"
            << (transport_.tableDir()/"meanEnergy_vs_reducedE")
            << "` is absent." << nl
            << "    Only the LMEA criterion is reported. See"
            << " Projects/SoEEDF/docs/lfa-vs-lmea.md." << endl;
    }
}


void plasmaSimulationDiagnostics::report()
{
    read();

    if (printLocalityValidity_)
    {
        reportLocalityValidity();
    }

    if (!printSpecies_ && !printElectromagnetics_)
        return;

    auto fmtRow = [](
        const std::string& label,
        const std::string& unit,
        scalar             minVal,
        scalar             maxVal,
        int                lw = 24) -> std::string
    {
        std::string line = "  " + label + " [" + unit + "]:";
        while (static_cast<int>(line.size()) < lw) line += ' ';
        line += "min = " + std::string(Foam::name(minVal).c_str());
        while (static_cast<int>(line.size()) < lw + 24) line += ' ';
        line += "max = " + std::string(Foam::name(maxVal).c_str());
        return line;
    };

    Info<< nl << "  Plasma Diagnostics" << nl
        << "  " << std::string(52, '-').c_str() << nl;

    if (printSpecies_)
    {
        Info<< "  Species number densities" << nl;
        for (label i = 0; i < transport_.species().nSpecies(); ++i)
        {
            const volScalarField& n = transport_.species().numberDensity(i);
            const word& name = transport_.species().speciesName(i);

            Info<< "  " << fmtRow
                (
                    name.c_str(), "m^-3",
                    gMin(n), gMax(n)
                ).c_str() << nl;
        }
    }

    if (printElectromagnetics_)
    {
        
        const volScalarField& Emag = transport_.species().em().Emag();
        const volScalarField& surfCharge = 
                                         transport_.species().em().surfCharge();
        
        scalar sigmaMin =  GREAT;
        scalar sigmaMax = -GREAT;
        forAll(surfCharge.boundaryField(), patchi)
        {
            const scalarField& sp = surfCharge.boundaryField()[patchi];
            if (sp.size())
            {
                sigmaMin = min(sigmaMin, min(sp));
                sigmaMax = max(sigmaMax, max(sp));
            }
        }
        reduce(sigmaMin, minOp<scalar>());
        reduce(sigmaMax, maxOp<scalar>());
        if (sigmaMin > sigmaMax) { sigmaMin = sigmaMax = 0; }

        Info<< "  Electromagnetics" << nl;

        Info<< "  " << fmtRow
            (
                "Emag",       "V/m",
                gMin(Emag),   gMax(Emag)
            ).c_str() << nl;

        Info<< "  " << fmtRow
            (
                "surfCharge",      "C/m^2",
                sigmaMin,  sigmaMax
            ).c_str() << nl;
    }

    Info<< "  " << std::string(52, '-').c_str() << nl << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
