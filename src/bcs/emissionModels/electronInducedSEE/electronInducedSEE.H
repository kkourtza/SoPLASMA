/*---------------------------------------------------------------------------*\
  File: electronInducedSEE.H
  Part of: SoPLASMA

Class
    Foam::emissionModels::electronInducedSEE

Description
    Secondary electron emission by ELECTRON impact:

        Gamma_emit = delta(E) max(0, Gamma_e . n)

    with delta the yield per incident electron. Unlike ion-induced emission,
    delta can EXCEED ONE, and where it does the surface emits more electrons
    than it collects -- which is what makes an insulator charge POSITIVE.

    VAUGHAN'S UNIVERSAL CURVE, which is why this mechanism is generalizable
    where gamma_eff(E/N) is not: the shape is material-independent and only two
    numbers move it (Vaughan, IEEE Trans. Electron Devices 36 (1989) 1963, with
    the correction in 40 (1993) 830):

        v = (E - E0) / (Emax - E0),            E0 = 12.5 eV
        delta/delta_max = (v exp(1-v))^k       k = 0.56  for v < 1
                                               k = 0.25  for 1 <= v <= 3.6
        delta/delta_max = 1.125 v^-0.35        for v > 3.6

    delta rises from the threshold E0, peaks at Emax, and falls beyond. It
    crosses delta = 1 twice, and those two crossings decide the sign of the
    charge a dielectric accumulates -- the reason this matters in a barrier
    discharge at all.

    delta_max AND Emax ARE REQUIRED INPUTS, from the material library when the
    entry carries them or explicitly. THEY ARE NOT DEFAULTED, and for
    insulators they are deliberately NOT SHIPPED: the published spread for an
    insulator (MgO spans roughly 4-25) is far larger than the difference
    between materials, which is exactly the case G2 says must not be given a
    point value. Metals are better characterised and a few are in the library.

    WHERE THE INCIDENT ENERGY COMES FROM, and its limitation. The model uses the
    MEAN electron energy at the wall (`meanE` under LMEA). That is an
    approximation: delta is a strongly non-linear function of energy, so the
    correct treatment integrates delta(E) over the electron energy
    DISTRIBUTION at the wall, not evaluates it at the mean. The two agree only
    where delta is locally linear. This is stated in the report, not buried.

    AND WHERE IT MATTERS LEAST. In an atmospheric-pressure discharge electrons
    reach a wall with a few eV -- far below Vaughan's 12.5 eV threshold -- so
    delta is ZERO and this mechanism contributes nothing. It matters in a
    low-pressure sheath, in a plasma display panel, and anywhere the wall sees
    energetic electrons. A case that switches it on at atmospheric pressure and
    sees no effect is seeing the right answer.

SourceFiles
    electronInducedSEE.C

\*---------------------------------------------------------------------------*/

#ifndef electronInducedSEE_H
#define electronInducedSEE_H

#include "emissionModel.H"

namespace Foam
{
namespace emissionModels
{

/*---------------------------------------------------------------------------*\
                      Class electronInducedSEE Declaration
\*---------------------------------------------------------------------------*/

class electronInducedSEE
:
    public emissionModel
{
    //- Peak yield [-]
    const scalar deltaMax_;

    //- Energy of the peak [eV]
    const scalar energyDeltaMax_;

    //- Threshold energy [eV]; Vaughan's 12.5 eV unless overridden
    const scalar E0_;

    //- Whether the mean-energy caveat has been reported
    mutable bool caveatReported_;


    // Private Member Functions

        //- One of the two required parameters, from the material library when
        //  the entry has it, else explicit. FATAL when neither.
        scalar required(const word& key, const word& libKey) const;


public:

    TypeName("electronInducedSEE");

    electronInducedSEE(const fvPatch& p, const dictionary& dict);

    virtual ~electronInducedSEE() = default;

    //- Vaughan's universal curve. Public and static so the unit test can
    //  check the CLASS against the same closed form it checks the curve with.
    static scalar vaughanYield
    (
        const scalar E,
        const scalar deltaMax,
        const scalar Emax,
        const scalar E0
    );

    virtual tmp<scalarField> emittedFlux() const;

    virtual void report(Ostream& os) const;
};


} // End namespace emissionModels
} // End namespace Foam

#endif

// ************************************************************************* //
