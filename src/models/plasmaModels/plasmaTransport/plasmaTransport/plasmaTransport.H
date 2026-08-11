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
    Foam::plasmaTransport

Description
    Manager class for per-species plasma transport models.

    Owns one plasmaTransportModel per active species. Delegates all
    transport physics to the individual models and aggregates their
    outputs for use by the electromagnetics solver.

    Provides:
        - Per-species particle flux (for number density equations)
        - Total electrical conductivity (for semi-implicit Poisson)
        - Total diffusive charge source (for semi-implicit Poisson)
        - Wall flux list (for surface charge update via plasmaSpecies)

    Has no knowledge of reactions or energy models.

SourceFiles
    plasmaTransportI.H
    plasmaTransport.C
    
\*---------------------------------------------------------------------------*/

#ifndef plasmaTransport_H
#define plasmaTransport_H

#include "fvCFD.H"
#include "fvMesh.H"
#include "PtrList.H"
#include "regIOobject.H"

#include "plasmaSpecies.H"
#include "HashSet.H"
#include "plasmaTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration keeps plasmaReactionRates.H out of this public
// header; the include lives in plasmaTransport.C. Including it here forces
// the include path onto bcs, tools and the solver, all of which then fail
// to build.
class plasmaReactionRates;
class plasmaChemistry;
    
class photoionizationModel;
/*---------------------------------------------------------------------------*\
                    Class plasmaTransport Declaration
\*---------------------------------------------------------------------------*/

class plasmaTransport
:
    public regIOobject
{
private:

    //- Disallow copy
    plasmaTransport(const plasmaTransport&) = delete;
    void operator=(const plasmaTransport&) = delete;

    //- Construct per-species transport models
    void constructTransportModels();

    //- Build rates_ from the `chemistry` sub-dictionary, if present
    void readChemistry(const dictionary& dict);

    //- Integrate the chemistry over `dt` in every active cell.
    //
    //  Called twice per timestep with dt/2, either side of the transport
    //  solve. Strang splitting is second order in time where applying the
    //  chemistry once, before or after, is only first order -- and splitting
    //  this way, by halving the CHEMISTRY substep rather than the transport
    //  one, leaves the transport time discretisation untouched.
    void applyChemistry(const scalar dt);

    //- Mechanism-driven species sources. Returns false when no mechanism is
    //  configured, so the caller falls back to the legacy hard-coded fits and
    //  existing cases keep running unchanged.
    bool mechanismSourceTerms
    (
        List<autoPtr<fvScalarMatrix>>& eqns,
        const volScalarField& ne,
        const volScalarField& Emag
    );

protected:

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Reference to species database
    plasmaSpecies& species_;

    //- Per-species transport models
    PtrList<plasmaTransportModel> transportModels_;

    //- Convective particle flux [1/s]
    PtrList<surfaceScalarField> convectiveFlux_;

    //- Diffusive particle flux [1/s]
    PtrList<surfaceScalarField> diffusiveFlux_;

    //- Total particle surface flux [1/s]
    PtrList<surfaceScalarField> particleFlux_;

    //- Mechanism-driven electron-impact rates. Null when no `chemistry`
    //  sub-dictionary is present, in which case the legacy hard-coded Townsend
    //  fits are used -- so this is additive and existing cases keep running.
    autoPtr<plasmaReactionRates> rates_;

    //- How the source terms are formed. smReactions sums k_j n_e n_target;
    //  smTownsend uses per-reaction (alpha_j/N) N v_drift n_e. Algebraically
    //  identical; they differ only if the solver's drift velocity differs from
    //  the one the tables were built with.
    enum sourceModelType { smReactions, smTownsend };
    sourceModelType sourceModel_ = smReactions;

    //- Reference composition from the mechanism, for background-gas targets
    dictionary mechComposition_;

    //- Products already warned about as untransported, so the warning is
    //  emitted once rather than every cell of every timestep
    mutable HashSet<word> warnedMissing_;

    //- Accept a lumped-ion model in which charged products may be discarded.
    //  Off by default: discarding a charged product breaks charge conservation
    //  silently, which is not something to warn about.
    bool allowChargeNonConservation_ = false;

    //- Operator-split chemistry. Null unless `chemistrySolver` asks for it,
    //  in which case it replaces the in-equation source terms entirely --
    //  applying both would count the chemistry twice.
    autoPtr<plasmaChemistry> chem_;

    //- none: sources go into the transport equations, as before.
    //  ode:   Strang splitting with a stiff ODE substep.
    enum chemSolverType { csNone, csODE };
    chemSolverType chemistrySolver_ = csNone;

    //- Cells below this electron density are skipped by the chemistry.
    //  Most of a streamer domain is undisturbed gas at the floor, where the
    //  ODE integrates a state that does not change; skipping it is what makes
    //  a per-cell stiff solve affordable at all.
    scalar chemActivityThreshold_ = 0;

    //- Active-cell count reported once rather than every step
    bool chemCellsReported_ = false;

    //- Worst relative charge residual seen in the first step
    scalar chemResidualMax_ = 0;

    //- Charge-imbalance warning emitted once
    mutable bool chargeWarned_ = false;

    //- Largest relative charge imbalance reported so far, so the diagnostic
    //  tracks the whole run rather than whichever step happened to be first
    mutable scalar chargeRelMax_ = 0;

    //- Species that carries the ion produced by photoionization. In air the
    //  VUV emitted by excited N2 ionizes O2 (Zheleznyak et al. 1982), so the
    //  photoelectron is paired with O2+. When this species is not carried the
    //  source falls back to the first positive ion, reproducing the previous
    //  behaviour exactly.
    word photoIonSpecies_ = "O2p";

    volScalarField k_eff_;
    volScalarField alpha_;
     volScalarField alphaDx_;
//- Pure ionization source α·μ·E·ne [1/(m3 s)] (for photoionization)
volScalarField S_iz_;

autoPtr<photoionizationModel> photoionization_;
public:

    //- Runtime type information
    TypeName("plasmaTransport");

    // Constructors

        //- Construct from mesh and species
        plasmaTransport
        (
            const fvMesh& mesh,
            plasmaSpecies& species
        );

    //- Destructor
    virtual ~plasmaTransport();

    // Member functions

        //- Correct all transport models
        void correctTransportModels();

        //- Solve all species transport equations
        void solve();

        //- Update surface charge
        void updateSurfaceCharge();

        //- Return the LHS matrix for a specific species
        tmp<fvScalarMatrix> transportMatrix(const label i);

        //- Return total electrical conductivity [S/m]
        tmp<volScalarField> electricalConductivity() const;

        //- Return total diffusive charge source [C/m^3/s]
        tmp<volScalarField> diffusiveChargeSource() const;

        virtual bool writeData(Ostream& os) const;

    // Access

        //- Metadata access
        inline const fvMesh& mesh() const;
        inline plasmaSpecies& species();
        inline const plasmaSpecies& species() const;

        //- Access the EM model through species
        inline const electromagneticsModel& em() const;

        //- Model access
        inline plasmaTransportModel& model(const label i);
        inline const plasmaTransportModel& model(const label i) const;

        //- Access fluxes
        inline surfaceScalarField& particleFlux(const label i);
        inline const surfaceScalarField& particleFlux(const label i) const;
        inline PtrList<surfaceScalarField>& particleFlux();
        inline const PtrList<surfaceScalarField>& particleFlux() const;
        inline surfaceScalarField& convectiveFlux(const label i);
        inline const surfaceScalarField& convectiveFlux(const label i) const;
        inline PtrList<surfaceScalarField>& convectiveFlux();
        inline const PtrList<surfaceScalarField>& convectiveFlux() const;
        inline surfaceScalarField& diffusiveFlux(const label i);
        inline const surfaceScalarField& diffusiveFlux(const label i) const;
        inline PtrList<surfaceScalarField>& diffusiveFlux();
        inline const PtrList<surfaceScalarField>& diffusiveFlux() const;

        inline const volScalarField& k_eff() const { return k_eff_; };
        const volScalarField& alpha() const noexcept { return alpha_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaTransportI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
