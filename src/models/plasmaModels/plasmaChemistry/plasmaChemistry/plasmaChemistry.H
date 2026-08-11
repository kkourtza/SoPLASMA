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
#include "fileName.H"
#include "HashTable.H"
#include "autoPtr.H"

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


    // Private Member Functions

        label stateIndex(const word& name) const;
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

        label nSpecie() const { return species_.size(); }
        const wordList& species() const { return species_; }

        //- Mechanism ids of the table-driven reactions, in the order
        //  integrate() expects their rate coefficients
        const wordList& tabulatedIds() const { return tabulatedIds_; }

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
        void integrate
        (
            scalarField& n,
            const scalarField& kTab,
            const scalar Tgas,
            const scalar dt
        ) const;

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
