/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::plasmaSimulationDiagnostics

Description
    Class for calculating and reporting runtime diagnostics in plasma
    simulations.

    This class reads the 'plasmaSimulationDiagnostics' dictionary to determine
    which fields (e.g., species number densities, electric field, surface
    charge) should have their minimum and maximum values reported to the console
    at each time step.

SourceFiles
    plasmaSimulationDiagnostics.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaSimulationDiagnostics_H
#define plasmaSimulationDiagnostics_H

#include "fvCFD.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations. plasmaRateTable is held by pointer and its
// definition is needed only in the .C -- this header is PUBLIC, and
// including a model header here forces the include path on every
// consumer of the tools library (it broke soPlasmaFoam, 2026-09-05).
class plasmaTransport;
class plasmaRateTable;
/*---------------------------------------------------------------------------*\
                 Class plasmaSimulationDiagnostics Declaration
\*---------------------------------------------------------------------------*/

class plasmaSimulationDiagnostics
{
private:
// Private data

        //- Reference to the Time object
        const Time& runTime_;

        //- Reference to the plasma transport model
        const plasmaTransport& transport_;

        //- Dictionary of control parameters
        dictionary dict_;

    // Diagnostic switches

        //- Switch to enable/disable printing species min/max values
        Switch printSpecies_;

        //- Switch to enable/disable printing electromagnetic fields min/max
        Switch printElectromagnetics_;

        //- TIME-LOCALITY VALIDITY ADVISORY. Default ON.
        //
        //  Dias & Guerra, J. Phys. D 58 185204 (2025), give two computable
        //  conditions that say whether a closure is inside its domain of
        //  validity, rather than leaving it to judgement. With
        //
        //      nu_EN  = |d(E/N)/dt| / (E/N)        the forcing rate
        //
        //  their eq. (2) is  nu_eps/nu_EN >> 1  for the LFA, and their eq. (5)
        //  is  nu_m/nu_EN >> 1  for the LEA. Since nu_m >> nu_eps under
        //  virtually all conditions, the LEA's condition is the far weaker one
        //  -- which is the formal statement of why it survives transients that
        //  break the LFA.
        //
        //  Reported as the FRACTION OF CELLS failing each condition, and only
        //  when that fraction is non-negligible: an advisory that fires on
        //  every run is one nobody reads.
        //
        //  See Projects/SoEEDF/docs/lfa-vs-lmea.md for the evidence and for
        //  the verification cases these criteria come from.
        Switch printLocalityValidity_;

        //- What ">>" means here. A ratio below this counts as a failure.
        scalar localityMargin_;

        //- Report only when this fraction of cells fails.
        //
        //  0.1, CALIBRATED not chosen -- see localityFieldFraction_.
        scalar localityReportFraction_;

        //- Only judge cells whose E/N is at least this fraction of the
        //  domain's own maximum.
        //
        //  WITHOUT THIS THE CRITERION IS MEANINGLESS, and Dias & Guerra say so
        //  themselves: the comparison is useful only "if the numerical
        //  precision is sufficient to prevent ARTIFICIALLY HIGH VALUES of
        //  nu_EN at low reduced electric fields". nu_EN is a RELATIVE rate,
        //  |d(E/N)/dt|/(E/N), so it diverges wherever the field is near zero
        //  -- which is most of a discharge domain, and none of the physics.
        //
        //  CALIBRATED 2026-09-05 against a published answer rather than taste.
        //  Dias & Guerra report, for their nanosecond pulse, LFA errors of
        //  40-50% for air at 10 Torr and "minimal" at 100 Torr. Evaluating
        //  their eq. (2) on the same tables, over the pulse:
        //
        //    restriction        10 Torr frac below 10   100 Torr
        //    none (whole tail)        95.3%               46.7%   <- flags both
        //    E/N > 10% of peak        90.2%                3.0%   <- correct
        //
        //  So margin 10 with reportFraction 0.1 and this restriction FLAGS
        //  10 Torr and PASSES 100 Torr, which is the published verdict. The
        //  median ratios are 5.5 and 55 -- a factor of ten, the density ratio.
        //  Script: Projects/SoEEDF/validation/dias2025/calibrate_locality.py
        scalar localityFieldFraction_;

        //- E/N at the previous report, for d(E/N)/dt. Sized on first use.
        autoPtr<volScalarField> reducedEOld_;

        //- Time at which reducedEOld_ was taken.
        scalar reducedEOldTime_;

        //- One-shot note that the LFA criterion cannot be evaluated.
        bool lfaCriterionUnavailableReported_;

        //- `meanEnergy_vs_reducedE`, for the mean energy under LFA.
        //
        //  The LFA criterion needs nu_eps = P_loss/<eps>, and LFA transports
        //  no mean energy -- which left the criterion that matters most to an
        //  LFA user unevaluable in the first version of this advisory. The
        //  table is the sweep's own E/N -> <eps> map, i.e. exactly the
        //  local-equilibrium mean energy the LFA assumes, so reading it here
        //  is not an approximation but a statement of what LFA already means.
        //
        //  Loaded lazily: absent for a mechanism compiled before the column
        //  existed, which is not an error, only a criterion that cannot be
        //  reported.
        mutable autoPtr<plasmaRateTable> meanEofEN_;
        mutable bool meanEofENTried_;

        //- Reports to skip before the advisory may speak.
        //
        //  A case starting from zero field changes E/N by O(1) PER STEP while
        //  the field establishes, so nu_EN is enormous and both criteria fail
        //  everywhere. That is start-up, not a statement about the discharge,
        //  and an advisory that fires on every run's first steps is one nobody
        //  reads. Measured 2026-09-05 on a 5-step needle-DBD bed: 100% of
        //  71715 cells, every step.
        label localityWarmup_;

        //- Reports seen so far, against localityWarmup_.
        label localityReports_;

        //- The last message's rounded percentages, so a STEADY verdict is
        //  stated once rather than every step. Repeating an unchanged warning
        //  is how a real one gets skipped.
        label lastPctLFA_;
        label lastPctLEA_;

        //- Evaluate and report the two locality criteria.
        void reportLocalityValidity();

    // Private Member Functions

        //- Read controls from dictionary
        void read();

public:

    // Constructors

        //- Construct from Time and plasmaTransport
        plasmaSimulationDiagnostics
        (
            const Time& runTime,
            const plasmaTransport& transport
        );

    //- Destructor. OUT OF LINE: autoPtr<plasmaRateTable> needs the
    //  complete type to destroy, and that type is deliberately not
    //  visible in this header.
    virtual ~plasmaSimulationDiagnostics();

    // Member Functions

        //- Report the requested diagnostics
        void report();

    // Accessors (Getters)

        //- Is species reporting enabled?
        inline bool printSpecies() const
        {
            return printSpecies_;
        }

        //- Is electromagnetics reporting enabled?
        inline bool printElectromagnetics() const
        {
            return printElectromagnetics_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
