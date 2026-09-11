/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
\*---------------------------------------------------------------------------*/

#include "plasmaChemistryCantera.H"
#include "error.H"

#ifndef HAVE_CANTERA

// Built without Cantera. The class still exists so the selector can name it
// and fail with a message that says what to do, rather than the backend simply
// not appearing in the list of valid options.
Foam::plasmaChemistryCantera::plasmaChemistryCantera
(
    const fileName&, const wordList&, const wordList&
)
{
    FatalErrorInFunction
        << "`chemistryBackend cantera` was selected, but libplasmaChemistry"
        << " was built without Cantera." << nl
        << "    Rebuild with CANTERA_DIR pointing at a Cantera install"
        << " prefix (one containing include/cantera/core.h), or use"
        << " `chemistryBackend native`, which evaluates the same heavy"
        << " reactions from the mechanism dictionary." << nl
        << exit(FatalError);
}

Foam::plasmaChemistryCantera::~plasmaChemistryCantera() = default;
void Foam::plasmaChemistryCantera::setBackground
(
    const word& canteraName,
    const scalar n
)
{
    const size_t k = sol_->thermo()->speciesIndex(std::string(canteraName), false);
    if (k != Cantera::npos) background_[label(k)] = n;
}


Foam::wordList Foam::plasmaChemistryCantera::unmatched() const { return wordList(); }
void Foam::plasmaChemistryCantera::productionLoss
(
    const scalarField&, const scalar, scalarField&, scalarField&
) const {}
void Foam::plasmaChemistryCantera::setBackground(const word&, const scalar) {}
void Foam::plasmaChemistryCantera::jacobian
(
    const scalarField&, const scalar, scalarSquareMatrix&
) const {}
bool Foam::plasmaChemistryCantera::setState(const scalarField&, const scalar) const
{ return false; }

#else

#include "cantera/core.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/numerics/eigen_sparse.h"
#include "cantera/base/AnyMap.h"

#include <vector>

/* Cantera 3 and Cantera 4 differ in the shape of the few calls used here:
   4 passes std::span where 3 passed a raw pointer. The two are wrapped rather
   than the file being forked, because the difference is purely a calling
   convention -- there is no version-specific PHYSICS below this point, and a
   fork would be an invitation for the two copies to drift.

   CANTERA_MAJOR comes from Make/options, which reads it out of the installed
   config.h. Cantera 4 additionally requires C++20 for std::span; that flag is
   set there too. */
#if defined(CANTERA_MAJOR) && CANTERA_MAJOR >= 4
    #include <span>
    #define CT_OUT(v)  std::span<double>((v).data(), std::size_t((v).size()))
    #define CT_IN(v)   std::span<const double>((v).data(), std::size_t((v).size()))
#else
    #define CT_OUT(v)  (v).data()
    #define CT_IN(v)   (v).data()
#endif

// Avogadro's number in Cantera's units: particles per kmol. Cantera works in
// kmol/m^3 and this code in m^-3, and the conversion is applied here so the
// caller sees the same units whichever backend is selected.
static constexpr Foam::scalar NA_PER_KMOL = 6.02214076e26;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasmaChemistryCantera::plasmaChemistryCantera
(
    const fileName& heavyMechanism,
    const wordList& species,
    const wordList& canteraNames
)
:
    species_(species)
{
    try
    {
        sol_ = Cantera::newSolution(heavyMechanism, "", "none");
    }
    catch (const std::exception& e)
    {
        FatalErrorInFunction
            << "Cantera could not load " << heavyMechanism << ":" << nl
            << e.what() << nl
            << "    This file is written by mechc alongside the .foam"
            << " dictionary and holds the HEAVY reactions only." << nl
            << exit(FatalError);
    }

    const auto& thermo = sol_->thermo();

    // Map the state vector onto Cantera's species order. The mechanism names
    // are used, not the OpenFOAM-safe ones: `N2p` means nothing to Cantera,
    // `N2+` does.
    toCantera_.setSize(species_.size(), -1);
    forAll(species_, i)
    {
        const std::string nm(canteraNames[i]);

        // `false` suppresses the throw on an unknown name: Cantera 4 raises by
        // default where 3 returned npos, and an unmatched species is expected
        // here -- the electron takes no part in heavy chemistry and is absent
        // from the heavy-only mechanism by construction. The overload taking
        // the flag exists in both series, so this is not version-dependent.
        const size_t k = thermo->speciesIndex(nm, false);
        if (k != Cantera::npos)
        {
            toCantera_[i] = label(k);
        }
    }

    fromCantera_.setSize(label(thermo->nSpecies()), -1);
    forAll(toCantera_, i)
    {
        if (toCantera_[i] >= 0) fromCantera_[toCantera_[i]] = i;
    }

    // Sized on the KINETICS, not the thermo: getCreationRates writes one
    // entry per species the kinetics manager knows, which for a multi-phase
    // manager exceeds the phase's own count. Equal here, but a buffer sized
    // from the wrong object overflows the heap silently.
    const label nK = label(max(sol_->kinetics()->nTotalSpecies(), thermo->nSpecies()));
    X_.setSize(nK, Zero);
    creation_.setSize(nK, Zero);
    destruction_.setSize(nK, Zero);
    background_.setSize(nK, Zero);

    // Ask for a COMPLETE composition derivative. Cantera can be told to skip
    // the third-body and falloff contributions, which is a legitimate speed /
    // accuracy trade for a preconditioner but not for a Rosenbrock Jacobian:
    // those are exactly the terms that make a mechanism stiff.
    try
    {
        Cantera::AnyMap settings;
        settings["skip-third-bodies"] = false;
        settings["skip-falloff"] = false;
        sol_->kinetics()->setDerivativeSettings(settings);
    }
    catch (const std::exception&)
    {
        // Not every kinetics manager accepts settings; the defaults are then
        // whatever it implements, and the probe below is what decides.
    }

    // Probe the analytic derivative once, at a state where it is defined.
    // Calling it is the only test: the method exists on the base class and
    // throws NotImplementedError for managers that do not provide it.
    try
    {
        X_ = 1.0;
        thermo->setMoleFractions(CT_IN(X_));
        thermo->setState_TD(300.0, thermo->meanMolecularWeight()*1e-3);
        sol_->kinetics()->netProductionRates_ddCi();
        analyticJacobian_ = true;
    }
    catch (const std::exception&)
    {
        analyticJacobian_ = false;
    }
    X_ = Zero;
}


bool Foam::plasmaChemistryCantera::setState
(
    const scalarField& n,
    const scalar Tgas
) const
{
    auto thermo = sol_->thermo();

    // Cantera's state is set from the NUMBER DENSITIES directly, never from a
    // fixed pressure. setState_TPX would impose the ideal-gas density implied
    // by (T, p) and rescale every concentration to match -- the plasma's own
    // composition is the state here, and a cell in a streamer head is not at
    // the nominal pressure of the case anyway.
    //
    // So: mole fractions plus the total molar density that the densities
    // themselves define. Cantera normalises the mole fractions, and the
    // separately-set density restores the absolute scale that normalisation
    // removed, leaving every concentration exactly as supplied.
    //
    // The working arrays are members, sized once at construction: this runs
    // per cell per timestep, and allocating them each time would cost more
    // than the kinetics evaluation it feeds.
    X_ = background_;
    forAll(toCantera_, i)
    {
        if (toCantera_[i] >= 0) X_[toCantera_[i]] += max(n[i], scalar(0));
    }

    const scalar nSum = sum(X_);
    if (nSum <= 0) return false;

    thermo->setMoleFractions(CT_IN(X_));
    thermo->setState_TD(Tgas, (nSum/NA_PER_KMOL)*thermo->meanMolecularWeight());
    return true;
}


void Foam::plasmaChemistryCantera::jacobian
(
    const scalarField& n,
    const scalar Tgas,
    scalarSquareMatrix& dfdy
) const
{
    if (analyticJacobian_)
    {
        // The state has to be pushed here rather than assumed from the last
        // productionLoss(): the ODE solver asks for the Jacobian at states the
        // rates were not evaluated at.
        if (!setState(n, Tgas)) return;

        const Eigen::SparseMatrix<double> J =
            sol_->kinetics()->netProductionRates_ddCi();

        // Scatter the Cantera block into the state-ordered matrix. Iterating
        // the sparse structure rather than indexing it: a mechanism of any
        // size is overwhelmingly zero, and coeff() on a compressed matrix is a
        // binary search per entry.
        for (label k = 0; k < J.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(J, k); it; ++it)
            {
                const label i = fromCantera_[label(it.row())];
                const label j = fromCantera_[label(it.col())];
                if (i >= 0 && j >= 0) dfdy(i, j) += it.value();
            }
        }
        return;
    }

    // Finite differences on the Cantera rates. Correct but nSpecie+1 kinetics
    // evaluations, so it is a fallback and reported as one -- never a silent
    // substitution of the other backend's Jacobian, which would be wrong in
    // precisely the cases that made Cantera necessary.
    //
    productionLoss(n, Tgas, P0_, L0_);

    forAll(n, i) P0_[i] -= L0_[i]*max(n[i], scalar(0));   // P0_ now holds f(n)

    nPert_ = n;
    forAll(n, j)
    {
        const scalar nj = max(n[j], scalar(0));

        // Step relative to the species' own density, floored at one particle
        // per cubic metre: the densities span twenty orders of magnitude, so a
        // single absolute step is either roundoff for the small ones or a
        // finite jump for the large ones.
        //
        // 1e-4 rather than the more obvious sqrt(eps) ~ 1e-8. The quotient
        // differences two rates of the size of the LARGEST rate in the column,
        // not of the entry being computed, so its roundoff is eps*f_max/h --
        // and for a structurally zero entry that noise is the entire answer.
        // At 1e-6 it was 7% of the column scale; each decade of h buys a
        // decade of it back, against a truncation error that only grows
        // linearly. Measured: 1e-4 puts both at ~1e-3 of the column scale.
        const scalar h = max(1e-4*nj, 1.0);

        nPert_[j] = nj + h;
        productionLoss(nPert_, Tgas, Pp_, Lp_);

        forAll(n, i)
        {
            const scalar fp = Pp_[i] - Lp_[i]*max(nPert_[i], scalar(0));
            dfdy(i, j) += (fp - P0_[i])/h;
        }

        nPert_[j] = n[j];
    }
}


Foam::plasmaChemistryCantera::~plasmaChemistryCantera() = default;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plasmaChemistryCantera::setBackground
(
    const word& canteraName,
    const scalar n
)
{
    const size_t k = sol_->thermo()->speciesIndex(std::string(canteraName), false);
    if (k != Cantera::npos) background_[label(k)] = n;
}


Foam::wordList Foam::plasmaChemistryCantera::unmatched() const
{
    DynamicList<word> out;
    forAll(toCantera_, i)
    {
        if (toCantera_[i] < 0) out.append(species_[i]);
    }
    return wordList(out);
}


void Foam::plasmaChemistryCantera::productionLoss
(
    const scalarField& n,
    const scalar Tgas,
    scalarField& P,
    scalarField& L
) const
{
    P.setSize(species_.size()); P = Zero;
    L.setSize(species_.size()); L = Zero;

    if (!setState(n, Tgas)) return;

    auto kin = sol_->kinetics();
    kin->getCreationRates(CT_OUT(creation_));
    kin->getDestructionRates(CT_OUT(destruction_));

    // Cantera returns kmol/m^3/s. Production is used as is; destruction is
    // converted to a loss COEFFICIENT by dividing out the species' own density,
    // so the sink stays implicit exactly as on the native path -- the caller
    // must not be able to tell which backend produced these.
    //
    // The floor is the same one the native path uses: below one particle per
    // cubic metre there is nothing to lose, and dividing by a vanishing density
    // produces a coefficient that overflows the matrix.
    forAll(toCantera_, i)
    {
        const label k = toCantera_[i];
        if (k < 0) continue;

        P[i] += creation_[k]*NA_PER_KMOL;

        const scalar d = destruction_[k]*NA_PER_KMOL;
        const scalar ns = max(n[i], scalar(0));
        if (ns > 1.0)
        {
            L[i] += d/ns;
        }
        else
        {
            P[i] -= d;
        }
    }
}


// ************************************************************************* //

#endif // HAVE_CANTERA

// ************************************************************************* //
