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

#include "plasmaOuterRelaxation.H"
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

class localEnergyEnergyModel;

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

    //- Enable the electron energy in the chemistry ODE (once) and return the
    //  LMEA model, or nullptr. Called from BOTH chemistry paths.
    const localEnergyEnergyModel* enableChemistryEnergy();

    //- True on the last PIMPLE outer iteration of this timestep. Retained
    //  for reporting; the source formulation no longer needs it, because
    //  recomputing the source each outer iteration is idempotent.
    //
    //  It mattered when the chemistry was an in-place update that had to
    //  bracket the TIMESTEP rather than the outer iteration; that scheme is
    //  gone -- see computeChemistrySources(). The counter is kept because
    //  plasmaTransport holds no pimpleControl, and is still flushed on a
    //  timestep change so a PIMPLE exiting early on residual control cannot
    //  leave it stale.
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
    //  csAdaptive        per cell: implicitRate where valid, ode where stiff,
    //                    judged by two THRESHOLD tests (L*dt, predicted change)
    //  csAdaptiveError   same switch, judged by an estimated LOCAL ERROR
    //                    against a tolerance -- one step vs two half-steps
    //  csImplicitRate    linearised at the current iterate, loss via fvm::Sp
    //  csODE             stiff substep integrated through the step
    //  csExplicitSource  legacy in-equation explicit source
    enum chemSolverType
    {
        csAdaptive,
        csAdaptiveError,
        csImplicitRate,
        csODE,
        csExplicitSource
    };
    chemSolverType chemistrySolver_ = csAdaptive;

    //- Either per-cell switch. They differ only in HOW a cell is judged, so
    //  everything else -- the P/L assembly, the outer-corrector requirement,
    //  the Picard reporting, the order banner -- treats them identically.
    inline bool isAdaptive() const
    {
        return chemistrySolver_ == csAdaptive
            || chemistrySolver_ == csAdaptiveError;
    }

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

    //- Volume-integrated gas-heating power [W] in each channel, this step.
    //  Split rather than summed so the two-step (Popov) share can be read off
    //  at one plasma state instead of inferred by differencing two runs.
    mutable scalar gasHeatPrompt_ = 0.0;
    mutable scalar gasHeatHeavy_ = 0.0;
    mutable label gasHeatLastReport_ = 0;

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

    //- Pass the transport cross term to the per-cell stiff integration.
    //  False integrates every stiff cell as CLOSED. See the measured numbers
    //  at the read site in plasmaTransport.C -- the term currently flips the
    //  sign of the Co 1.5 streamer error, so this is a live question, not a
    //  legacy switch.
    bool chemCrossTerm_ = true;

    //- Cap on the cross term's sink coefficient, in units of 1/deltaT.
    //  Zero disables the sink, which is the default; see
    //  plasmaChemistryODE::extSinkLimit_ for the measurements behind that.
    scalar chemCrossSinkLimit_ = 0.0;

    //- Density below which a relative change is meaningless [m^-3]
    //
    //  This is the knob that decides how much of the mesh the change test
    //  claims. At 1e6 a trace species sitting at 1e7 and driven to 1e9 in one
    //  step registers a relative change of 100 and flags its cell, even where
    //  n_e is 1e20 and that species cannot affect anything. Measured on the
    //  1.15M-cell streamer: the change test flagged 100% of cells from step 3
    //  onward while `L*dt` flagged 5-7%, so `adaptive` was running the stiff
    //  ODE path over the whole domain -- first order everywhere, which is what
    //  the temporal-order sweep was seeing.
    //
    //  A floor is a blunt instrument: the right value depends on the case's own
    //  densities, which is why the error-estimating variant exists. Exposed so
    //  it can at least be set; the default is unchanged. Note the real fix for
    //  the 100%-flagged problem is chemChangeScale_ below, not this.
    scalar chemChangeFloor_ = 1.0e6;

    //- Fraction of the cell's DOMINANT density used as the absolute part of
    //  the change test's scale: scale_s = n0_s + f*max_s(n0_s).
    //
    //  Without it the change test is purely relative, normalising every
    //  species by its own density, so a trace species moving orders of
    //  magnitude outranks a dominant species moving a few percent. That is
    //  what made the test flag every cell of the streamer case.
    //
    //  f = 0 is the purely relative test, and is the DEFAULT: it is what every
    //  existing case and every validated result was produced with, so the
    //  default does not move under them.
    //
    //  MEASURED, on the 1.15M-cell streamer at Co 1.5 (2026-08-18):
    //
    //    f = 0     100% of cells flagged from step 3 on -- `adaptive` is the
    //              first-order stiff path over the whole mesh
    //    f = 1e-3  ~83% flagged, and only 0.03% at step 2
    //
    //  So 1e-3 is a strict improvement and the right value to reach for first,
    //  but it does NOT fix the problem: it plateaus because in the quiescent
    //  far field the cell's own dominant density is itself tiny (N2 and O2
    //  ground states are the background, not transported), so the absolute
    //  part of the scale shrinks with it. Any purely LOCAL normalisation hits
    //  this -- whether a cell's chemistry matters is a question about its
    //  contribution to the global solution. That is the argument for the
    //  error-estimating variant rather than a tuned threshold.
    scalar chemChangeScale_ = 0.0;

        //- Significance floor for Co_chem, as a fraction of the largest
        //  density anywhere: scale_s = max(n_s, chemCoScale_*max_all(n)).
        //  Same idea as chemChangeScale_ and deliberately the same vocabulary,
        //  but a separate knob because it steers the TIME STEP rather than the
        //  per-cell solver choice. 1e-3 keeps five decades of dynamic range
        //  below the peak, which covers the streamer channel and excludes the
        //  quiescent far field. See maxChemStateRate().
        scalar chemCoScale_ = 1e-3;

    // ---- csAdaptiveError: the local-error switch --------------------------
    //
    //  The threshold tests above ask "is this cell stiff / does it move a lot".
    //  Neither is the question that matters, which is "would linearising this
    //  cell put an error into the solution worth paying an ODE solve to avoid".
    //  This variant asks that directly: take the linearised step over dt, take
    //  two linearised half-steps, and use the difference as a Richardson
    //  estimate of the local error. Route to the stiff path only when that
    //  estimate exceeds a tolerance.
    //
    //  The scale is the part that matters, and is where the threshold tests
    //  fail. It carries a GLOBAL absolute part:
    //
    //      scale = absFrac * max_domain(n_e) + relTol * n_cell
    //
    //  so a cell can only claim the stiff path by making an error comparable to
    //  a fraction of the discharge's own peak. A purely local scale cannot do
    //  this: in the quiescent far field every local density is small, so any
    //  relative measure saturates there -- which is exactly why
    //  chemChangeScale_ plateaus at 83% of cells.

    //- Which species the error norm is measured over
    //    neElectron  the electron density alone (default)
    //    neCharged   every charged species
    //    neAll       every transported species
    enum chemErrorNormType
    {
        neElectron,
        neCharged,
        neAll
    };
    chemErrorNormType chemErrorNorm_ = neElectron;

    //- Relative tolerance, against the cell's own density
    scalar chemErrorRelTol_ = 1.0e-2;

    //- Absolute tolerance, as a fraction of the DOMAIN-maximum electron
    //  density. Not a literal density: a fixed one would have to be re-chosen
    //  per case, which is the trap chemChangeFloor_ documents.
    scalar chemErrorAbsFrac_ = 1.0e-4;

    //- max_domain(n_e) for this pass, one gMax per outer iteration
    scalar chemErrorRefDensity_ = 0.0;

    //- Cells sent to the stiff path by the error test, per pass and peak-held
    label chemNstiffByError_ = 0;
    label chemNstiffByEnergy_ = 0;
    label chemNstiffByErrorPeak_ = 0;
    label chemNstiffByEnergyPeak_ = 0;

    //- Cells sent to the stiff path by each criterion, for reporting.
    //
    //  Per PASS: computeChemistrySources() runs once per outer iteration, so
    //  these are zeroed at the top of every pass and peak-held into the
    //  companions below. They used to accumulate instead, which is what made
    //  the report read "169443621 of 9171728 cells" -- a count summed over
    //  every iteration of every step, against a denominator that the in-place
    //  reduce() had already multiplied by the rank count.
    label chemNstiffByL_ = 0;
    label chemNstiffByChange_ = 0;

    //- Peak of the three stiff-cell counts since the last report, so the
    //  numbers printed together describe the SAME pass and stay comparable.
    label chemNstiffByLPeak_ = 0;
    label chemNstiffByChangePeak_ = 0;

    //- Stiff integrations that returned non-finite values since the last
    //  report. Summed rather than peak-held: for a failure count, "how many
    //  happened" is the useful number, not "how many at once".
    label chemODEFailuresSinceReport_ = 0;

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

    //- The previous step's chemExt_, and the slope derived from the pair.
    //
    //  A CONSTANT cross term caps the coupling at first order: measured
    //  p = 0.988 for `solver ode` on the 40x40 study with the outer loop
    //  converged to 1e-12 and chemExt_ at its fixed point. The transport rate
    //  varies across the step, and integrating a nonlinear chemistry against
    //  a constant stand-in perturbs the trajectory at O(dt) however exact
    //  that constant is.
    //
    //  chemExtSlope_ = (chemExt_ - chemExtPrev_)/dt is lagged by one step and
    //  therefore only O(dt) accurate -- which is enough, because it enters
    //  the substep multiplied by (x - dt/2) ~ dt, leaving T correct to
    //  O(dt^2). Zero on the first step, which is first order there alone.
    List<scalarField> chemExtPrev_;
    List<scalarField> chemExtSlope_;

    //- THE ENERGY'S OWN TRANSPORT RATE, and the state it was measured at.
    //
    //  The species hand the stiff integrator their transport rate through
    //  chemExt_, so the ODE evolves chemistry AND transport together and its
    //  end state is the real t^{n+1}. The ENERGY row did not: its transport
    //  was left to the field equation alone, so `yE[nSpecie_]` came back as a
    //  chemistry-only state while `yE[eIndex]` carried transport.
    //
    //  The reconstruction then formed eps = yE[nSpecie_]/yE[eIndex] -- a ratio
    //  of two states that do not correspond to the same physics -- and
    //  evaluated the Joule and loss terms there. That is wrong by O(dt), which
    //  caps the scheme at first order. MEASURED on the order-study bed:
    //  p = 0.79 for the energy, and 0.67 for the SPECIES too once
    //  `lookupVariable meanE` makes the rate coefficients depend on the same
    //  broken eps. With the chemistry insulated (lookupVariable reducedE) the
    //  species recover to 1.97 while the energy stays at 0.79 -- one cause,
    //  both symptoms.
    scalarField chemExtEps_;
    scalarField chemExtEpsPrev_;
    scalarField chemExtEpsSlope_;
    scalarField chemNEps0_;

    //- Instantaneous production and loss coefficient, for csImplicitRate
    List<scalarField> chemP_;

    //- ENERGY source from the chemistry ODE (Option 4), in the same
    //  production/loss-rate form the species use, so the energy equation can
    //  put the sink on its diagonal exactly as before.
    //
    //  Written only when the ODE integrates the energy; the LMEA model checks
    //  energyInChemistry() and falls back to its own hand-rolled source when
    //  it does not. Applying BOTH would count the source twice.
    //  REGISTERED FIELDS, looked up BY NAME by localEnergyEnergyModel.
    //
    //  Not passed by type: plasmaTransport links plasmaEnergy, so a
    //  plasmaEnergy -> plasmaTransport dependency would be circular. The
    //  registry is the seam that already exists between them.
    autoPtr<volScalarField> chemPeps_;
    autoPtr<volScalarField> chemLeps_;
    bool energyInChemistry_ = false;

    //- `energySource chemistry` in the chemistry dict: integrate the electron
    //  energy with the species instead of sourcing it in the energy model.
    //  Default OFF, so existing cases are unchanged.
    bool chemEnergyRequested_ = false;
    List<scalarField> chemL_;

    //- Outer iterations taken in the current timestep, and how many are
    //  configured
    label chemOuterCount_ = 0;
    label nOuterCorrectors_ = 1;

    //- ADAPTIVE OUTER-LOOP RELAXATION (opt-in:
    //  plasmaSimulationControls/outerCoupling/adaptiveRelaxation).
    //
    //  The species and electron-energy equations couple through a Picard
    //  iteration inside the outer loop, and above a critical deltaT that
    //  iteration has gain > 1 and enters a period-2 limit cycle. A FIXED
    //  under-relaxation factor damps it but the usable value is case-specific
    //  AND time-varying -- measured on the LMEA benchmark, 0.9 was stable at
    //  6.19 ps and unstable at 8.03 ps in the same run. Aitken derives the
    //  factor from the iteration itself, per corrector.
    bool adaptiveRelax_ = false;

    //- The per-mesh JOINT relaxation coordinator, shared with plasmaEnergy.
    //  Not owned -- it lives on the registry.
    plasmaOuterRelaxation* relax_ = nullptr;

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

    //- Whether the electron energy is integrated by the chemistry ODE rather
    //  than sourced by the energy model itself.
    bool energyInChemistry() const { return energyInChemistry_; }

    //- End-state energy production [eV/m^3/s] and loss RATE [1/s]. Valid
    //  only when energyInChemistry().
    const volScalarField& chemPeps() const { return chemPeps_(); }
    const volScalarField& chemLeps() const { return chemLeps_(); }

    //- Runtime type information
    TypeName("plasmaTransport");

        //- Throw away the step just solved, restoring the start-of-step state.
        //
        //  Every transported density goes back to its own oldTime(), which is
        //  the last ACCEPTED state: oldTime() only rotates when the time INDEX
        //  advances, and a retry deliberately holds the index fixed, so it is
        //  still the right level to come back to -- and oldTime().oldTime(),
        //  which `backward` (BDF2) also needs, is untouched.
        //
        //  Nothing else needs restoring, which is a property of how this
        //  solver was already written: surfCharge resets from its own
        //  oldTime() before accumulating, chemN0_/chemExtPrev_ are keyed on
        //  the time index and so still hold start-of-step values, the
        //  photoionization and rate-table caches are index-keyed too, and the
        //  chemistry source is rebuilt from n^n on every outer iteration
        //  rather than from the partially updated field. e_vib was the one
        //  exception and is now baselined per index for exactly this reason.
        void discardStep();

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

        //- Advance the species over one timestep.
        //
        //  `finalIter` must be PIMPLE's finalIter() for this outer iteration.
        //  It no longer gates the chemistry: the source is recomputed from the
        //  start-of-step state on every outer iteration and is idempotent, so
        //  finalIter is used only for end-of-step reporting. Defaults to true
        //  so a single-corrector caller is correct without knowing any of this.
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

        //- Smallest Aitken relaxation factor used over the last time step,
        //  across every species. This is the outer coupling's STABILITY
        //  MARGIN: 1 means nothing needed damping at this deltaT, while a
        //  small value means the Picard iteration is only stable because it
        //  is being damped hard, i.e. deltaT is near the coupling's limit.
        //
        //  plasmaTimeControl consumes it to govern deltaT growth. Without
        //  that feedback a damped run simply walks deltaT back up to the new
        //  stability boundary and sits on it -- which is exactly what a fixed
        //  relaxation factor was measured doing.
        //
        //  Returns 1 when adaptive relaxation is off, so a case that does not
        //  use it imposes no limit.
        scalar relaxationMargin() const
        {
            return relax_ ? relax_->margin() : 1.0;
        }

        //- Largest joint relaxation factor over the last step. Reported beside
        //  the minimum so "omega sat near 0.5 all step" is distinguishable
        //  from "omega swung from 1.0 to 0.18" -- those describe different
        //  couplings and a single number conflates them.
        scalar relaxationOmegaMax() const
        {
            return relax_ ? relax_->omegaMaxStep() : 1.0;
        }
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

        //- Fastest FRACTIONAL RATE OF CHANGE OF THE STATE [1/s], over all
        //  species and cells. The source for Co_chem.
        //
        //      scale_s = max( n_s , chemCoScale_ * max_over_all(n) )
        //      rate    = max over (s, cells) of |P_s - L_s n_s| / scale_s
        //
        //  For electrons this reduces EXACTLY to |k_eff| = |S_ion - S_att|/n_e,
        //  which is what the legacy air-fitted limiter measured -- so the
        //  existing 0.9 calibration carries over rather than needing a fresh
        //  sweep. It is the mechanism-derived generalisation of it: any gas,
        //  any mechanism, no mu_ref and no constant eta.
        //
        //  Why the NET rate and not max_s(L_s), which was tried and reverted:
        //  MEASURED 2026-08-20, limiting on L pinned Co_chem at its cap every
        //  step, took the 0.5 ns bed from 52 to 113 steps and moved the answer
        //  55%. The driver was the ELECTRON (tau_chem = 1.76e-12 s), because
        //  P ~ L n -- quasi-equilibrium, where the channels are enormous and
        //  the STATE barely moves. That is what the stiff integrator already
        //  handles under its own error control, so limiting on it double-counts.
        //  Co_chem must measure how far the state MOVES.
        //
        //  The scale floor does two jobs at once. It stops a trace species at
        //  1e7 driven to 1e9 outranking a dominant species at 1e20 moving 5%
        //  (the documented chemChangeScale trap), and it suppresses the
        //  quiescent far field, where n is negligible but its fractional rate
        //  need not be -- which is what made the old mag(k_eff) throttle the
        //  step on benign attachment decay.
        //
        //  Returns 0 before the first chemistry evaluation and with
        //  `reactions none`; the caller must cope with 0.
        scalar maxChemStateRate() const;

        //- Fastest chemical loss rate over all species and cells [1/s], or 0
        //  when the chemistry path does not assemble P/L.
        //
        //  DIAGNOSTIC ONLY -- reported so stiffness stays visible, never used
        //  to steer the step. See maxChemStateRate() for why.
        //
        //  This exists because k_eff() is a NET rate, (S_ion - S_att)/n_e:
        //  where ionisation and attachment nearly cancel it approaches zero
        //  while BOTH processes are fast, so a limiter reading it alone is
        //  blind to exactly the regime it should constrain. The per-species
        //  loss coefficient has no such cancellation.
        //
        //  chemL_ is sized in TWO places: the in-equation path guarded by
        //  `csImplicitRate || isAdaptive()`, and computeChemistrySources(),
        //  which is guarded only by `!chem_ || dt <= 0` and therefore runs for
        //  every solver including pure `ode`. So this is normally live in all
        //  chemistry modes. It still returns 0 before the first chemistry
        //  evaluation, and with `reactions none`, so the caller must cope with
        //  0 rather than assume a value.
        scalar maxChemLossRate() const;

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
