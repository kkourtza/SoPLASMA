/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::plasmaChemistry

Description
    Reads the compiled mechanism and integrates the chemistry of a cell.

    Builds one reaction list from BOTH families in the mechanism dictionary --
    the electron-impact `processes` and the heavy `heavyReactions` -- so a
    single source-term assembly applies both and they cannot drift apart in
    convention. That has been the recurring failure in this work: an implicit
    electron reactant, a collider density counted twice, a units block that
    reached one consumer and not another.

    Species that the case does not transport (background N2, O2) are not state
    variables. Their densities enter as fixed parameters, folded in once at
    setup rather than re-multiplied on every right-hand-side evaluation.

    Integration is by OpenFOAM's ODESolver, selected in the dictionary:

        chemistrySolver ode;
        odeCoeffs { solver seulex; absTol 1e-8; relTol 1e-4; }

    `seulex` and `rodas23` are the stiff choices and the sensible defaults for
    plasma chemistry, whose rates span twenty orders of magnitude. An explicit
    solver will appear to work and then quietly demand timesteps far below the
    transport limit.

SourceFiles
    plasmaChemistry.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaChemistry_H
#define plasmaChemistry_H

#include "plasmaChemistryODE.H"
#include "ODESolver.H"
#include "dictionary.H"
#include "plasmaChemistryCantera.H"
#include "fileName.H"
#include "HashTable.H"
#include "autoPtr.H"
#include "scalarMatrices.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class plasmaChemistry Declaration
\*---------------------------------------------------------------------------*/

class plasmaChemistry
{
    // Private Data

        //- Species carried as state variables, in state order
        wordList species_;

        //- Charge number per state entry
        scalarField charge_;

        //- Fixed densities for species the case does not transport [m^-3]
        HashTable<scalar> background_;

        List<plasmaReactionSpec> reactions_;

        //- How many of the above are electron-impact, i.e. table-driven
        label nTabulated_ = 0;

        //- Mechanism process ids in the order their rate coefficients are
        //  expected, so the caller can fill the array unambiguously
        wordList tabulatedIds_;

        mutable scalarField kTab_;
        mutable autoPtr<plasmaChemistryODE> ode_;
        mutable autoPtr<ODESolver> solver_;

        dictionary odeDict_;

        //- `chemistryBackend`: who evaluates the HEAVY reactions.
        //
        //    native   the Arrhenius forms in the mechanism dictionary,
        //             evaluated by plasmaChemistryODE. No dependencies.
        //    cantera  the same reactions from <mech>.heavy.yaml, evaluated by
        //             libcantera. Adds falloff, PLOG, Chebyshev, reversible
        //             reactions with equilibrium constants from thermo data --
        //             forms the native parser does not implement.
        //
        //  Electron-impact rates come from the Boltzmann tables in BOTH cases.
        //  Letting Cantera evaluate those as well would double-count them and
        //  would need the EEDF injected per cell, which is the cost the tables
        //  exist to avoid.
        word backend_;

        autoPtr<plasmaChemistryCantera> cantera_;

        //- Cantera-side names of the carried species (`N2+` where OpenFOAM
        //  must say `N2p`), from the mechanism's speciesAlias table
        wordList canteraNames_;

        //- The same mapping for EVERY species the mechanism mentions, carried
        //  or not: an untransported background species needs translating too.
        HashTable<word> alias_;

        //- Heavy production and loss from Cantera, sized once
        mutable scalarField Pc_, Lc_;

        //- Largest relative disagreement between the native and Cantera heavy
        //  rates at the reference state, measured once at construction
        scalar backendMismatch_ = 0.0;


    // Private Member Functions

        label stateIndex(const word& name) const;

        //- Compare native and Cantera heavy P/L at a reference state.
        //  Returns the largest relative difference; see the constructor for
        //  why this is checked rather than assumed.
        scalar checkBackends(const scalar Tgas) const;
        void readMechanism(const fileName& path, const dictionary& dict);


public:

    // Constructors

        //- Construct from the mechanism dictionary and the species to carry
        plasmaChemistry
        (
            const fileName& mechanismDict,
            const wordList& species,
            const scalarField& charge,
            const dictionary& dict
        );


    // Member Functions

        //- Cap on the transport cross term's sink coefficient, in units of
        //  1/deltaT. See plasmaChemistryODE::extRef_ for what it trades.
        void setCrossSinkLimit(const scalar c)
        {
            ode_->setCrossSinkLimit(c);
        }

        label nSpecie() const { return species_.size(); }
        const wordList& species() const { return species_; }

        //- Mechanism ids of the table-driven reactions, in the order
        //  integrate() expects their rate coefficients
        const wordList& tabulatedIds() const { return tabulatedIds_; }

        //- Warn once if any rate coefficient evaluated negative and was
        //  clamped. COLLECTIVE -- call once per step from outside the cell
        //  loop, never from integrate(), which is per cell and would deadlock
        //  when ranks hold different cell counts.
        void reportRateGuard() const
        {
            ode_->reportNegativeRateGuard(tabulatedIds_);
        }

        //- Which backend evaluates the heavy reactions
        const word& backend() const { return backend_; }

        //- True when the Cantera backend supplies its analytic Jacobian
        bool canteraAnalyticJacobian() const
        {
            return cantera_ && cantera_->analyticJacobian();
        }

        //- Largest relative native-vs-Cantera disagreement found at startup
        scalar backendMismatch() const { return backendMismatch_; }

        label nReactions() const { return reactions_.size(); }
        label nTabulated() const { return nTabulated_; }

        //- Fixed background density of an untransported species [m^-3]
        void setBackground(const word& name, const scalar n);

        //- Integrate one cell from t to t+dt in place.
        //
        //  `kTab` holds the electron-impact rate coefficients at this cell,
        //  already interpolated and in tabulatedIds() order. They are held
        //  fixed across the substep: under operator splitting the field does
        //  not change while the chemistry integrates.
        //  `ext`, when given, is the constant TRANSPORT rate for this cell
        //  over the step. Supplying it makes the integrated trajectory the
        //  cell's actual one rather than a closed-cell one, which is what
        //  makes the integrated trajectory the cell's actual one.
        //
        //  `extSlope` is d(ext)/dt. Supplying it makes the cross term LINEAR
        //  across the substep instead of constant, which is what actually
        //  reaches second order: a constant ext leaves the scheme at p = 0.99
        //  even with the outer loop converged to 1e-12 (measured 2026-08-18,
        //  40x40 study). Centred so the mean stays `ext`.
        void integrate
        (
            scalarField& n,
            const scalarField& kTab,
            const scalar Tgas,
            const scalar dt,
            const scalarField* ext = nullptr,
            const scalarField* extSlope = nullptr
        ) const;

        //- Right-hand side and Jacobian at a state, for verification.
        void derivatives(const scalarField& n, const scalarField& kTab,
                         const scalar Tgas, scalarField& dndt) const;
        void jacobian(const scalarField& n, const scalarField& kTab,
                      const scalar Tgas, scalarField& dfdx,
                      scalarSquareMatrix& dfdy) const;

        //- Instantaneous production and loss coefficient at this state,
        //  dn_s/dt = P_s - L_s n_s. See plasmaChemistryODE::productionLoss.
        void productionLoss(const scalarField& n, const scalarField& kTab,
                            const scalar Tgas, scalarField& P,
                            scalarField& L) const;

        //- Heat released to the gas by the HEAVY reactions [eV m^-3 s^-1].
        //  Fast gas heating. Excludes the electron-impact channels, which are
        //  covered by the tabulated PgasN -- see plasmaChemistryODE.
        scalar heavyHeatRelease(const scalarField& n, const scalar Tgas) const;

        //- Relative charge residual of the right-hand side at this state.
        //  Identically zero for a charge-conserving mechanism, whatever the
        //  state, so a non-zero value is an unbalanced reaction and nothing
        //  else.
        scalar chargeResidual(const scalarField& n, const scalarField& kTab,
                              const scalar Tgas) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
