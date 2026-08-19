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
    Foam::plasmaEnergy

Description
    High-level manager for plasma energy models.

    This class owns and manages one energy model per plasma species (mainly 
    for electrons). It does not implement any energy physics itself. Instead it
    serves as a controller that:

        - Retrieves species-level energy configuration
        - Constructs the appropriate energy model (LFA, LEA, etc.)
        - Provides temperature fields to the transport models

    The plasmaEnergy object is typically created once in the solver and used
    during the solution loop to update the energy behavior of the plasma
    species.

SourceFiles
    plasmaEnergyI.H
    plasmaEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaEnergy_H
#define plasmaEnergy_H

#include "fvCFD.H"
#include "fvMesh.H"
#include "PtrList.H"

#include "plasmaSpecies.H"
#include "plasmaEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class plasmaEnergy Declaration
\*---------------------------------------------------------------------------*/

class plasmaEnergy
:
    public regIOobject
{
private:
    //- Flag to solve background gas energy equation
    Switch solveGasEnergy_;

    //- Flag to track if gas temperature is a field (solved or fixed)
    Switch isGasTempField_;

    //- The scalar value if isGasTempField_ == false
    dimensionedScalar TgasValue_;

    //- The volScalarField if isGasTempField_ == true
    autoPtr<volScalarField> TgasFieldPtr_;

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Reference to species database
    plasmaSpecies& species_;

    //- Reference to the electric field
    const volVectorField& E_;

    //- Per-species energy models
    PtrList<plasmaEnergyModel> energyModels_;

    //- Non-equilibrium vibrational energy density [J/m^3].
    //
    //  Carried separately from the gas temperature because it does not reach
    //  the gas within a nanosecond pulse: tau_VT is microseconds. Folding it
    //  into T would deliver the energy far too early and roughly double the
    //  prompt temperature rise, which is exactly the error the fast/slow
    //  heating split exists to avoid.
    autoPtr<volScalarField> eVibPtr_;

    //- Start-of-timestep vibrational energy, and the index it belongs to.
    //
    //  e_vib is the ONE piece of solution state in this solver that is not
    //  recomputed from oldTime(): it carries no fvm::ddt, so it has no stored
    //  old time level, and it used to be advanced in place with `+=`. That is
    //  correct exactly once per timestep -- and it is called once, under
    //  `finalIter`, so there was never a bug -- but it makes the field the
    //  only thing here that is not IDEMPOTENT within a step.
    //
    //  Everything else in this solver can be re-run within a timestep and land
    //  on the same answer: the species and T_gas are recomputed from
    //  oldTime(), surfCharge resets from its own oldTime() before
    //  accumulating, and the chemistry caches are keyed on timeIndex. Only
    //  this one had to be advanced exactly once, by convention rather than by
    //  construction.
    //
    //  Baselining it restores that property: the update becomes
    //  ev = ev_start + (Pvib - Qvt) dt rather than ev += ..., which is
    //  BIT-EXACT for a single call (Qvt is evaluated at the same start-of-step
    //  state either way) and simply repeatable for any number of calls. That
    //  is what lets a non-converged step be discarded and re-run.
    mutable scalarField eVibStart_;
    mutable label eVibTimeIndex_ = -1;

    //- Thermal conductivity of the gas [W/m/K]. Constant for now; the
    //  mixture-averaged form needs Lennard-Jones data the mechanism does not
    //  yet carry.
    scalar kappaGas_ = 0.026;

    //- Fixed V-T time [s] if positive, else the composition-dependent
    //  Millikan-White + Popov form.
    scalar tauVTfixed_ = -1;

    //- Gas pressure [atm], for those relaxation times.
    scalar pGasAtm_ = 1.0;


    //- Disallow copying
    plasmaEnergy(const plasmaEnergy&) = delete;
    void operator=(const plasmaEnergy&) = delete;

    //- Construct per-species models
    void constructModels();

public:

    //- Runtime type information
    TypeName("plasmaEnergy");

    // Constructors

        //- Construct from plasmaSpecies and mesh
        plasmaEnergy
        (
            plasmaSpecies& species,
            const fvMesh& mesh,
            const volVectorField& E
        );

    //- Destructor
    virtual ~plasmaEnergy() {}

    // Member functions

        //- Correct all models
        void correct();

        //- Solve the background gas energy equation over dt.
        //
        //  Q_gas and rho*c_v are supplied by the caller because they belong to
        //  the chemistry, not to the energy model: only plasmaTransport knows
        //  the reaction enthalpies, the EEDF power channels and the mechanism
        //  thermodynamics. This class owns the TEMPERATURE and the equation;
        //  it does not own the physics that heats it.
        //
        //  Pvib is the power into the vibrational reservoir, drained here on
        //  tauVT and returned to the gas -- so the reservoir and the
        //  temperature stay on one timeline.
        void solveGasEnergy
        (
            const volScalarField& Qgas,
            const volScalarField& rhoCv,
            const volScalarField& Pvib,
            const volScalarField& tauVT,
            const scalar dt
        );

        //- The gas temperature FIELD, by reference.
        //
        //  Tgas() returns a tmp, and taking the address of a tmp is a dangling
        //  pointer waiting to happen -- it is only safe today because the
        //  solved path happens to return a non-owning tmp. A consumer that
        //  wants to hold on to the field should hold this instead.
        //  Valid only when solvesGasEnergy() is true.
        const volScalarField& TgasField() const { return *TgasFieldPtr_; }

        //- Throw away the step just solved: T_gas back to its own oldTime().
        //  e_vib needs nothing -- it integrates from eVibStart_, which is
        //  keyed on the time index a retry holds fixed, so re-running simply
        //  recomputes it with the shorter deltaT.
        void discardStep()
        {
            if (!solveGasEnergy_ || !TgasFieldPtr_) return;
            volScalarField& T = *TgasFieldPtr_;
            T == T.oldTime();
            T.correctBoundaryConditions();
        }

        //- The vibrational reservoir [J/m^3].
        const volScalarField& eVib() const { return eVibPtr_(); }

        //- Is the gas energy equation being solved?
        bool solvesGasEnergy() const noexcept { return solveGasEnergy_; }

        // //- Return total total Joule heating [W/m^3]
        // tmp<volScalarField> totalJouleHeating() const;

        virtual bool writeData(Ostream& os) const;

    // Access

        //- Mesh and species metadata access
        inline const fvMesh& mesh() const { return mesh_; }
        inline const plasmaSpecies& species() const { return species_; }

        //- Model access
inline plasmaEnergyModel& model(const label i)
        {
            return energyModels_[i];
        }

        inline const plasmaEnergyModel& model(const label i) const
        {
            return energyModels_[i];
        }

        //- Temperature field access for species i
inline tmp<volScalarField> T(const label i) const
        {
            return energyModels_[i].T();
        }

        //- Return the background gas temperature as a tmp field
        tmp<volScalarField> Tgas() const;

        //- Return the background gas temperature as a dimensionedScalar
        const dimensionedScalar& TgasValue() const;

//- Return true if gas temperature is a single uniform value (not a field)
inline bool isGasUniform() const { return !isGasTempField_; }

//- Return true if gas temperature is NOT being solved (remains constant in time)
inline bool isGasConstant() const { return !solveGasEnergy_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaEnergyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
