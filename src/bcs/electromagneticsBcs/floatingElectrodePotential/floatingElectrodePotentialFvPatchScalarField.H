/*---------------------------------------------------------------------------*\
  File: floatingElectrodePotentialFvPatchScalarField.H
  Part of: SoPLASMA

Class
    Foam::floatingElectrodePotentialFvPatchScalarField

Description
    A conductor connected to NOTHING -- no supply, no ground. A probe, an
    isolated pin, a metal fragment, a floating guard ring, or one disconnected
    segment of a segmented electrode.

    The usual boundary condition is INVERTED. For a driven or grounded
    electrode the potential is known and the charge is whatever the circuit
    supplies; here the CHARGE is known (it is conserved, changing only by what
    the plasma delivers) and the POTENTIAL is the unknown.

    Two conditions define it:

      1. IT IS AN EQUIPOTENTIAL. Charge in a conductor redistributes over its
         own relaxation time, eps/sigma ~ 1e-18 s for a metal, instantaneous
         against every plasma timescale. So V = V_f, ONE unknown scalar shared
         by every face -- a CONSTRAINT, not a Dirichlet value and not
         zeroGradient.

      2. ITS TOTAL CHARGE IS CONSERVED:  Q(t) = Q0 + INT I_plasma dt'.
         On a conductor there is NO LOCAL sigma: charge redistributes, so only
         the total is meaningful. That is the essential difference from
         thinDielectricPotential, where sigma is local and stays where it
         lands.

    This class is only the CARRIER of the equipotential value. It is a uniform
    Dirichlet condition whose value is set each step by
    Foam::floatingElectrode, which owns the constraint: it is that class, not
    this one, that closes Gauss's law and computes V_f.

    THERE IS DELIBERATELY NO POTENTIAL ENTRY. Supplying one would make this a
    driven electrode, and accepting one here would let a case silently ask for
    a contradiction.

    Theory, the closed-form solution and the tests:
    docs/models/poisson_equation/floating-electrode.md

SourceFiles
    floatingElectrodePotentialFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef floatingElectrodePotentialFvPatchScalarField_H
#define floatingElectrodePotentialFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "floatingPotentialInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class floatingElectrodePotentialFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class floatingElectrodePotentialFvPatchScalarField
:
    public fixedValueFvPatchScalarField,
    public floatingPotentialInterface
{
    // Private data

        //- Initial total charge Q0 [C]. Zero is the honest default: an
        //  electrode that has never been connected to anything starts neutral.
        scalar initialCharge_;

        //- The floating potential V_f [V]. A RESULT, written out so the
        //  settling behaviour can be read rather than assumed.
        scalar floatingPotential_;


public:

    //- Runtime type information
    TypeName("floatingElectrodePotential");


    // Constructors

        floatingElectrodePotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        floatingElectrodePotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        floatingElectrodePotentialFvPatchScalarField
        (
            const floatingElectrodePotentialFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        floatingElectrodePotentialFvPatchScalarField
        (
            const floatingElectrodePotentialFvPatchScalarField&
        );

        floatingElectrodePotentialFvPatchScalarField
        (
            const floatingElectrodePotentialFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new floatingElectrodePotentialFvPatchScalarField(*this)
            );
        }

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new floatingElectrodePotentialFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- The initial total charge Q0 [C]
        virtual scalar initialCharge() const { return initialCharge_; }

        //- The floating potential [V], a result
        virtual scalar floatingPotential() const { return floatingPotential_; }

        //- Set the equipotential value. Called by Foam::floatingElectrode,
        //  which owns the constraint that determines it.
        virtual void setFloatingPotential(const scalar Vf);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
