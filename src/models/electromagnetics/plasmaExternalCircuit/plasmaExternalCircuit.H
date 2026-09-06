/*---------------------------------------------------------------------------*\
  SoPLASMA -- external circuit coupled to a driven electrode.

Class
    Foam::plasmaExternalCircuit

Description
    THE DISCHARGE DOES NOT SET ITS OWN CURRENT; THE CIRCUIT DOES.

    A gap held at a FIXED applied voltage above its breakdown value has nothing
    to limit the current it draws. Measured 2026-09-06 on the Grubert dc glow:
    the gap broke down at -179 V, correctly (Paschen for argon at
    pd = 0.75 Torr cm sits near its minimum), and the current then went
    0.51 -> 6.7 -> 40 -> 59 mA/cm^2 and was still climbing when the timestep
    collapsed. It passes STRAIGHT THROUGH the 0.511 mA/cm^2 the reference
    reports, because nothing in a fixed-voltage transient selects an operating
    point. That is physics, not a numerical defect: a normal glow has negative
    differential resistance and is UNSTABLE at fixed voltage, which is why real
    glow discharges are run through a ballast.

    Grubert avoid the problem by solving the TIME-INDEPENDENT equations, so
    their operating point is computed rather than integrated towards. A
    transient solver needs the circuit that the experiment actually has.

    With a series ballast the electrode potential becomes an OUTPUT,

        V_electrode(t) = V_source(t) - R * I_circuit,

    and the operating point is self-selecting: draw more current, lose more
    voltage across the ballast, draw less current.

    I_circuit is Sato's I_total = I_cond + I_disp from plasmaDischargeCurrent,
    which is by construction the current an external circuit carries. It is
    PASSED IN rather than looked up, so this class has no dependency on the
    diagnostic library and cannot grow a second, independently-derived current
    -- two expressions for one quantity is how they come to disagree.

    THIS IS THE FIRST AND SIMPLEST MODEL OF A CAPABILITY THAT SHOULD GROW:
    see doc/external-circuit-plan.md for the staged design (RC and RLC
    ballasts, current sources, matched RF networks, several electrodes on one
    netlist). `type` is the extension point; today it accepts `seriesResistor`
    only, and rejects anything else rather than silently doing something else.

\*---------------------------------------------------------------------------*/

#ifndef plasmaExternalCircuit_H
#define plasmaExternalCircuit_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "Function1.H"
#include "autoPtr.H"
#include "OFstream.H"

namespace Foam
{

class plasmaExternalCircuit
{
    // Private data

        const fvMesh& mesh_;

        //- Off unless an `externalCircuit` sub-dictionary asks for it.
        bool enabled_;

        //- Circuit topology. Only `seriesResistor` today.
        word type_;

        //- The driven electrode patch this circuit feeds.
        word electrode_;

        //- Patch index of `electrode_` on this mesh.
        label patchi_;

        //- Open-circuit source voltage V_source(t) [V].
        autoPtr<Function1<scalar>> source_;

        //- Ballast resistance [Ohm].
        scalar R_;

        //- OPTIONAL extra damping. The implicit update below does not need
        //  it -- it is unconditionally stable -- so this defaults to 1 (none)
        //  and exists only for a case that wants to slow the approach to the
        //  operating point deliberately.
        scalar relax_;

        //- The electrode potential [V]: a RESULT, not an input.
        scalar V_;

        //- Whether V_ has been set from a measured current yet.
        bool started_;

        //- postProcessing/externalCircuit/circuit.csv
        autoPtr<OFstream> file_;


    // Private Member Functions

        plasmaExternalCircuit(const plasmaExternalCircuit&) = delete;
        void operator=(const plasmaExternalCircuit&) = delete;


public:

    //- Construct from the mesh and the controls dictionary
    plasmaExternalCircuit(const fvMesh& mesh, const dictionary& dict);

    //- Is a circuit active?
    inline bool enabled() const { return enabled_; }

    //- The electrode potential this circuit is imposing [V]
    inline scalar V() const { return V_; }

    //- The patch this circuit drives
    inline const word& electrode() const { return electrode_; }

    //- Apply the circuit: set the electrode potential from the measured
    //  CONDUCTION current and the gap capacitance, and write the trace.
    //  Does nothing when disabled.
    //
    //  IMPLICIT, UNCONDITIONALLY. The electrode sees a ballast in series with
    //  the gap capacitance, so V = V_source - R*(I_cond + C dV/dt) is an RC
    //  differential equation, not an algebraic relation. Evaluated explicitly
    //  it amplifies by R*C/dt per step -- MEASURED 885 at dt = 1e-10 s on the
    //  Grubert case, which took the electrode from 0 to -1.25e17 V in five
    //  steps. Backward Euler on the same relation,
    //
    //      V^{n+1} = (V_src + (RC/dt) V^n - R I_cond) / (1 + RC/dt),
    //
    //  is unconditionally stable and costs nothing extra. There is no explicit
    //  option: it was tried, it is unusable, and keeping it would only offer a
    //  way to reproduce the divergence.
    void update
    (
        const scalar Icond,
        const scalar Cgap,
        volScalarField& ePotential
    );
};

} // End namespace Foam

#endif
