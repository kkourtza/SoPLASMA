/*---------------------------------------------------------------------------*\
    Mixture thermodynamics from the mechanism's own NASA7 coefficients.

    Shared by the 0-D reactor and the coupled solver's energy equation, which
    is the point: a heat capacity that disagreed between the two would make the
    0-D validation of the CFD meaningless.

\*---------------------------------------------------------------------------*/

#ifndef janafMixture_H
#define janafMixture_H

#include "scalar.H"
#include "scalarField.H"
#include "wordList.H"
#include "dictionary.H"
#include "DynamicList.H"
#include "FixedList.H"
#include "error.H"

namespace Foam
{

//- Mixture thermodynamics from the mechanism's NASA7 coefficients.
//
//  The heat capacity and the density both used to be hard-coded for air: a
//  polynomial fit to 79/21 N2/O2, and a fixed mean molecular mass of
//  28.96 g/mol. Both are exact here instead, and for ANY mechanism:
//
//      rho    = SUM_i n_i m_i                     (mass, so exactly conserved)
//      c_v    = SUM_i Y_i (c_p,i(T) - R/W_i)
//      c_p,i  = (R/W_i) (a0 + a1 T + a2 T^2 + a3 T^3 + a4 T^4)
//
//  Electrons are excluded: they are at T_e, not T_gas, so they have no place
//  in a mixture heat capacity (and carry ~1e-5 of the mass regardless).
class janafMixture
{
    static constexpr scalar RU_ = 8314.462618;   // J/(kmol K)

    DynamicList<scalar> W_;                      // kg/kmol
    DynamicList<scalar> Tlow_, Thigh_, Tcommon_;
    DynamicList<FixedList<scalar, 7>> hi_, lo_;
    DynamicList<label> idx_;                     // index into the state vector
    bool ok_ = false;

public:

    //- Read the `speciesThermo` block, pairing it with the state vector.
    //  Species present in the state but absent from the block are a hard
    //  error rather than a silent omission: a missing N2 would quietly drop
    //  three quarters of the heat capacity and the run would still finish.
    janafMixture
    (
        const dictionary& mech,
        const wordList& species,
        const word& electron
    )
    {
        if (!mech.found("speciesThermo")) return;
        const dictionary& td = mech.subDict("speciesThermo");

        forAll(species, s)
        {
            if (species[s] == electron) continue;
            if (!td.found(species[s]))
            {
                FatalErrorInFunction
                    << "species " << species[s] << " is integrated by the"
                    << " reactor but has no entry in `speciesThermo`." << nl
                    << "Recompile the mechanism with mechc." << exit(FatalError);
            }
            const dictionary& e = td.subDict(species[s]);
            W_.append(e.get<scalar>("molWeight"));
            Tlow_.append(e.get<scalar>("Tlow"));
            Thigh_.append(e.get<scalar>("Thigh"));
            Tcommon_.append(e.get<scalar>("Tcommon"));
            hi_.append(e.get<FixedList<scalar, 7>>("highCpCoeffs"));
            lo_.append(e.get<FixedList<scalar, 7>>("lowCpCoeffs"));
            idx_.append(s);
        }
        ok_ = true;
    }

    bool valid() const { return ok_; }

    //- Mass density [kg/m^3] from number densities [1/m^3].
    scalar rho(const scalarField& n) const
    {
        scalar r = 0.0;
        forAll(idx_, i)
        {
            r += max(n[idx_[i]], scalar(0))*W_[i]/6.02214076e26;
        }
        return r;
    }

    //- Mixture c_p [J/(kg K)]. The isobaric reactor heats against this rather
    //  than c_v, because the work done expanding against the ambient pressure
    //  comes out of the same deposited energy. For air the two differ by
    //  gamma = 1.4, so choosing the wrong one is a 29% error in dT -- which is
    //  larger than most of the physics being argued about.
    scalar cp(const scalarField& n, const scalar T) const
    {
        scalar num = 0.0, den = 0.0;
        forAll(idx_, i)
        {
            const scalar mi = max(n[idx_[i]], scalar(0))*W_[i];
            if (mi <= 0.0) continue;
            const scalar t = min(max(T, Tlow_[i]), Thigh_[i]);
            const FixedList<scalar, 7>& a = (t < Tcommon_[i]) ? lo_[i] : hi_[i];
            const scalar cpR = a[0] + t*(a[1] + t*(a[2] + t*(a[3] + t*a[4])));
            num += mi*cpR*RU_/W_[i];
            den += mi;
        }
        return (den > 0.0) ? num/den : 0.0;
    }

    //- Mixture c_v [J/(kg K)] at the current composition and temperature.
    scalar cv(const scalarField& n, const scalar T) const
    {
        scalar num = 0.0, den = 0.0;
        forAll(idx_, i)
        {
            const scalar mi = max(n[idx_[i]], scalar(0))*W_[i];   // ~mass
            if (mi <= 0.0) continue;
            // Clamped to the polynomial's own validity range. A NASA7 fit
            // extrapolates catastrophically -- the quartic term dominates
            // and c_p can go negative, which would drive T the wrong way.
            const scalar t = min(max(T, Tlow_[i]), Thigh_[i]);
            const FixedList<scalar, 7>& a = (t < Tcommon_[i]) ? lo_[i] : hi_[i];
            const scalar cpR = a[0] + t*(a[1] + t*(a[2] + t*(a[3] + t*a[4])));
            num += mi*(cpR - 1.0)*RU_/W_[i];      // c_v = c_p - R/W
            den += mi;
        }
        return (den > 0.0) ? num/den : 0.0;
    }
};


} // End namespace Foam

#endif
