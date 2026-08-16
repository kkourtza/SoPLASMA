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

// Forward declaration rather than an include: janafMixture lives in
// libplasmaChemistry, and pulling it through this PUBLIC header would force
// every consumer of plasmaTransport.H -- the boundary conditions, the tools,
// the solver -- to carry that include path for a type none of them name.
class janafMixture;
class plasmaEnergy;

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
    //- Integrate the chemistry over `dt` from the start-of-step state and
    //  store the mean rate in chemRR_. Does NOT modify the species fields.
    //
    //  Modifying them directly is what the first implementation did, and it
    //  is wrong in OpenFOAM: fvm::ddt(n) differences against n.oldTime(),
    //  which is stored before the chemistry runs and does not follow it.
    //  The transport solve then sees the chemical change as a transport
    //  change and partially undoes it -- measured directly, oldTime kept the
    //  pre-chemistry value at every step. Delivering the chemistry as a
    //  SOURCE keeps ddt and oldTime consistent by construction, which is
    //  what reactingFoam does with its RR field.
    void computeChemistrySources(const scalar dt);

    //- True on the last PIMPLE outer iteration of this timestep. Retained
    //  for reporting; the source formulation no longer needs it, because
    //  recomputing the source each outer iteration is idempotent.
    //
    //  solve() is called once per OUTER ITERATION, not once per timestep, and
    //  Strang splitting has to bracket the timestep. plasmaTransport holds no
    //  pimpleControl, so the count is tracked here against nOuterCorrectors
    //  read from fvSolution -- with a flush on timestep change, so that a
    //  PIMPLE that exits early on residual control cannot drop the trailing
    //  half-step.
    bool finalOuterIteration();

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
    //  none          sources in the equations, electron-impact only (legacy)
    //  ode           whole mechanism, MEAN rate from a stiff ODE over the
    //                step. Lie composition: robust, first order.
    //  implicitRate  whole mechanism, INSTANTANEOUS rate at the current
    //                iterate, loss implicit. Second order with BDF2, because
    //                that is the right-hand side BDF2 asks for.
    //  adaptive      per CELL: implicitRate where the linearisation resolves
    //                the chemistry, the stiff ODE substep where it does not.
    //                The criterion is local because the stiffness is -- a
    //                streamer head can be stiff while the bulk is not, so
    //                switching globally would either pay the ODE cost
    //                everywhere or lose robustness where it is needed.
    // TWO ORTHOGONAL AXES, deliberately on two keys.
    //
    // `chemistry/reactions` says WHICH reactions exist -- physics.
    // `chemistry/solver`    says HOW they are integrated -- numerics.
    //
    // They used to share one key, `chemistrySolver`, whose value `none`
    // actually meant "electron-impact reactions, integrated explicitly". A
    // physics choice was therefore hidden inside something named after a
    // solver, which is how the shipped tutorial came to stop reproducing the
    // streamer validation without anything in the case saying so.

    //- WHICH reactions are evaluated
    //  rxAll             electron-impact + heavy (recombination, ion-ion
    //                    neutralisation, detachment, quenching)
    //  rxElectronImpact  electron-impact only. IONS NEVER DECAY: every loss
    //                    channel for ions lives in the heavy set.
    //  rxNone            no chemistry sources of any kind
    enum chemReactionSet { rxAll, rxElectronImpact, rxNone };
    chemReactionSet chemReactions_ = rxAll;

    //- HOW those reactions are integrated. No bearing on which ones exist.
    //  csAdaptive        per cell: implicitRate where valid, ode where stiff
    //  csImplicitRate    linearised at the current iterate, loss via fvm::Sp
    //  csODE             stiff substep integrated through the step
    //  csExplicitSource  legacy in-equation explicit source
    enum chemSolverType
    {
        csAdaptive,
        csImplicitRate,
        csODE,
        csExplicitSource
    };
    chemSolverType chemistrySolver_ = csAdaptive;

    //- Cells below this electron density are skipped by the chemistry.
    //  Most of a streamer domain is undisturbed gas at the floor, where the
    //  ODE integrates a state that does not change; skipping it is what makes
    //  a per-cell stiff solve affordable at all.
    scalar chemActivityThreshold_ = 0;

    //- Diagnostics are reported every `chemReportInterval_` timesteps, not
    //  once. Reporting only at start-up samples the run at t=0, where the
    //  densities are at their floor and the chemistry is at its LEAST stiff --
    //  the least informative moment there is. A discharge becomes stiff when
    //  it ignites, which is after the one sample a start-up-only report takes.
    //
    //  0 disables periodic reporting (start-up only, the old behaviour).
    label chemReportInterval_ = 0;

    //- Timestep index of the last report, per report family
    label chemLastReport_ = -1;
    label chemPicardLastReport_ = -1;

    //- Running maxima since the last report. Peak values matter more than
    //  instantaneous ones: a single step whose chemistry was unresolved is a
    //  defect in the solution whether or not the step that gets printed
    //  happens to show it.
    scalar chemStiffnessPeak_ = 0;
    label  chemNstiffPeak_ = 0;
    scalar chemPicardPeak_ = 0;
    scalar chemChargePeak_ = 0;

    //- Active-cell count reported once rather than every step
    bool chemCellsReported_ = false;

    //- The transport advisory is a one-shot report, deferred to the first
    //  chemistry REPORT rather than the first chemistry evaluation. At t=0 the
    //  densities sit at their floor and the chemistry is at its least active,
    //  so lifetimes sampled there are not the ones that matter -- the same
    //  reason chemReportInterval does not sample at start-up either.
    bool transportAdvisoryDone_ = false;

    //- Peak loss coefficient per species since the run began, so the advisory
    //  reports the SHORTEST lifetime the species actually reached.
    scalarField advisoryLpeak_;

    //- Reference diffusivity for the advisory, m^2/s. A neutral in air at
    //  1 bar, 300 K; species with their own diffusivity model use that instead.
    scalar advisoryD_ = 2.0e-5;

    //- Worst relative charge residual seen in the first step
    scalar chemResidualMax_ = 0;

    //- Worst relative charge residual of the INTEGRATED source, and how many
    //  cells the explicit source would drive negative
    scalar chemSourceChargeMax_ = 0;
    label  chemNegativeCells_ = 0;
    label  chemUnprojectedCells_ = 0;
    label  chemActiveCount_ = 0;

    //- Stiffness as the solver meets it, max over cells and species of L*dt
    scalar chemStiffness_ = 0;

    //- L*dt above which a cell is integrated rather than linearised
    scalar chemStiffLimit_ = 1.0;

    //- Relative change of a significant species that one LINEARISED step is
    //  allowed to make before the cell is integrated instead.
    //
    //  L*dt alone is not sufficient, and this is the lesson from the 0-D
    //  solver comparison (docs/validation-solvers.md in SoEEDF). L*dt tests
    //  the DIAGONAL -- the species' own loss. The instability of the
    //  linearised step lives OFF the diagonal: a pair such as attachment and
    //  detachment each see the other at the old iterate, so the pair can grow
    //  without bound while every individual L*dt stays comfortably below the
    //  limit. Screened on L*dt alone a 0-D afterglow reached 1e85 m^-3 with
    //  every step reported safe.
    //
    //  So the linearised outcome is PREDICTED and checked, not merely
    //  estimated. Costs one extra pass over the species in a cell that was
    //  already having P and L evaluated.
    scalar chemChangeLimit_ = 0.1;

    //- Density below which a relative change is meaningless [m^-3]
    scalar chemChangeFloor_ = 1.0e6;

    //- Cells sent to the stiff path by each criterion, for reporting
    label chemNstiffByL_ = 0;
    label chemNstiffByChange_ = 0;

    // ---- G2: gas energy ---------------------------------------------------
    //
    //  The gas temperature was a CONSTANT read from a dictionary. Every rate
    //  that depends on it -- three-body recombination, V-T relaxation, the
    //  EEDF's rotational and superelastic channels -- was therefore frozen at
    //  a value the discharge immediately leaves behind. Solving it closes the
    //  thermal-ionisation feedback: heating lowers the gas density, which
    //  raises E/N at fixed E, which raises ionisation, which heats further.
    //
    //  Solved as a proper fvScalarMatrix rather than a per-cell update, so it
    //  gets OpenFOAM's ddt schemes, its boundary conditions and its linear
    //  solvers -- and so conduction is available the moment a gradient exists.
    bool gasHeating_ = false;

    //- Borrowed view of plasmaEnergy's temperature field. Not owned: the
    //  energy model owns T_gas, and duplicating it here is how the two would
    //  come to disagree.
    const volScalarField* TgasField_ = nullptr;

    //- Thermal conductivity of the gas [W/m/K]. Constant by default; the
    //  mixture-averaged form needs Lennard-Jones data the mechanism does not
    //  yet carry (deferred as G3t).
    scalar kappaGas_ = 0.026;

    //- Fixed V-T time [s] if positive; otherwise the composition-dependent
    //  Millikan-White + Popov form of vibRelax.H.
    scalar tauVTfixed_ = -1;

    //- Diagnostics, peak-held between reports like the chemistry ones.
    scalar TgasMax_ = 0, QgasMax_ = 0, eVibMax_ = 0;

    //- Where the sweep lives, for the power channels.
    fileName tableDir_;

    //- Gas pressure [atm], for the Millikan-White relaxation times.
    scalar pGasAtm_ = 1.0;

    //- Mixture thermodynamics from the mechanism's NASA7 block. Shared with
    //  the 0-D reactor so the two cannot disagree about a heat capacity.
    autoPtr<janafMixture> thermo_;

    //- Thermodynamics of the BACKGROUND gas, over the mechanism's reference
    //  composition. A separate object because the background species have no
    //  index in plasmaSpecies at all -- they are a reservoir, not transported
    //  fields -- so they cannot be folded into the transported state vector.
    //
    //  rho*c_v is additive over disjoint subsets of the gas (it is
    //  SUM_i n_i m_i c_v,i), so the two contributions simply add.
    autoPtr<janafMixture> bgThermo_;
    wordList bgNames_;
    scalarField bgX_;

        //- Solve the gas energy equation and the vibrational reservoir over dt.
        void solveGasEnergy(const scalar dt);

        //- Gas temperature in a cell: the solved field when gas heating is on,
        //  the dictionary constant when it is not. One accessor so the two
        //  cannot diverge, which they did while three call sites each read the
        //  dictionary themselves.
        inline scalar TgasCell(const label celli) const
        {
            return (gasHeating_ && TgasField_)
                 ? (*TgasField_)[celli]
                 : TgasConst_;
        }

    scalar TgasConst_ = 300.0;

    //- Cells taking each path, for reporting
    label chemNstiff_ = 0;
    label chemODEFailures_ = 0;

    //- Relative change of the chemistry source between consecutive outer
    //  iterations. Small means the Picard iteration is converging and
    //  implicitRate is resolving the chemistry; not small means it is not.
    scalar chemPicardChange_ = 0;
    List<scalarField> chemSrcPrev_;
    bool   chemPicardReported_ = false;
    scalar chemSourceChargeAfter_ = 0;

    //- Timestep whose start-of-step state is held in chemN0_
    label chemTimeIndex_ = -1;

    //- Species densities at the START of the timestep, n^n.
    //
    //  The chemistry is always integrated FROM this state, never from the
    //  partially-updated field, so recomputing the source on every outer
    //  iteration is idempotent rather than cumulative.
    List<scalarField> chemN0_;

    //- Mean chemical production rate over the timestep, added to the
    //  transport equations as an explicit source.
    List<scalarField> chemRR_;

    //- Mean TRANSPORT rate over the timestep, from the previous outer
    //  iteration. Fed into the chemistry ODE as a constant source so that the
    //  trajectory it integrates is the cell's actual one rather than a
    //  closed-cell one -- which is what lifts the coupling to second order
    //  once the outer loop converges. Zero on the first iteration, where the
    //  scheme degrades gracefully to the closed-cell (first-order) form.
    List<scalarField> chemExt_;

    //- Instantaneous production and loss coefficient, for csImplicitRate
    List<scalarField> chemP_;
    List<scalarField> chemL_;

    //- Outer iterations taken in the current timestep, and how many are
    //  configured
    label chemOuterCount_ = 0;
    label nOuterCorrectors_ = 1;

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
        //- Advance the species over one timestep.
        //
        //  `finalIter` must be PIMPLE's finalIter() for this outer iteration.
        //  solve() is called once per outer iteration, and Strang splitting
        //  has to bracket the TIMESTEP: the leading half-step goes on the
        //  first call of a timestep and the trailing one on the last, and
        //  only the caller knows which that is. Counting iterations here
        //  instead is wrong as soon as residual control lets PIMPLE exit
        //  early. Defaults to true so a single-corrector caller is correct
        //  without knowing any of this.
        void solve(const bool finalIter = true);

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
