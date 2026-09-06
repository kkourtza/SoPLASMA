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

        //- Circuit topology: seriesResistor, seriesRC or currentSource.
        word type_;

        //- The driven electrode patch this circuit feeds.
        word electrode_;

        //- Patch index of `electrode_` on this mesh.
        label patchi_;

        //- Open-circuit source voltage V_source(t) [V].
        //  Unset for `currentSource`, which imposes a current instead.
        autoPtr<Function1<scalar>> source_;

        //- Regulated current MAGNITUDE I_set(t) [A], for `currentSource`.
        //  A magnitude: the polarity comes from the sign of compliance_, so
        //  the user never has to know the internal sign of I_cond.
        autoPtr<Function1<scalar>> setCurrent_;

        //- Compliance (rail) voltage [V], SIGNED, for `currentSource`.
        //  Every real supply has one, and it is also what makes the
        //  pre-ignition phase well posed: with no plasma there is no voltage
        //  at which I_set flows.
        scalar compliance_;

        //- Ballast resistance [Ohm].
        scalar R_;

        //- Shunt capacitance across the electrode [F], in PARALLEL with the
        //  gap's own. Zero for `seriesResistor`.
        //
        //  This is what makes a real dc rig light gently: ~100 pF/m of coax
        //  gives R*C in the millisecond range, so the voltage creeps up and
        //  the gas breaks down at the lowest voltage that sustains it. With
        //  the gap capacitance alone (1.8e-16 F here) R*C is 17.7 ns and the
        //  electrode is slammed far past Paschen before any plasma exists.
        scalar Cext_;

        //- OPTIONAL extra damping. The implicit update below does not need
        //  it -- it is unconditionally stable -- so this defaults to 1 (none)
        //  and exists only for a case that wants to slow the approach to the
        //  operating point deliberately.
        scalar relax_;

        //- The electrode potential [V]: a RESULT, not an input.
        scalar V_;

        //- Whether V_ has been set from a measured current yet.
        bool started_;

        //- Previous accepted (V, I) pair, for the secant estimate of the
        //  plasma's differential conductance dI/dV.
        scalar Vprev_;
        scalar Iprev_;
        bool   havePrev_;

        //- Last conductance estimate [S], kept so a step with a negligible
        //  voltage change reuses it rather than discarding the linearisation.
        scalar g_;

        //- postProcessing/externalCircuit/circuit.csv
        autoPtr<OFstream> file_;


    // Private Member Functions

        plasmaExternalCircuit(const plasmaExternalCircuit&) = delete;
        void operator=(const plasmaExternalCircuit&) = delete;

        //- Secant estimate of the plasma's differential conductance
        //  g = |dI/dV| [S] from the last two accepted steps. Shared by every
        //  topology, deliberately: two copies of this would drift.
        void updateSecantConductance(const scalar Icond);


public:

    //- Construct from the mesh alone.
    //
    //  The circuit is read from configuration/boundaries, on the patch
    //  declared `kind circuitDrivenElectrode`. There is no separate circuit
    //  dictionary and no `electrode` key: the electrode a circuit feeds is
    //  the patch it was written on, so the name cannot be stated twice and
    //  cannot drift.
    explicit plasmaExternalCircuit(const fvMesh& mesh);

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
    //  is unconditionally stable AGAINST THE CAPACITIVE TERM -- but not
    //  against the plasma, whose current was still explicit. That failed
    //  post-breakdown: the discharge drew 8.5x the current the circuit could
    //  supply at that instant, the ballast demanded a gap voltage the source
    //  could not provide, and the CATHODE WAS DRIVEN POSITIVE (+43.5 V) with
    //  dt -> 0. Measured 2026-09-06.
    //
    //  So the plasma is linearised implicitly too. With
    //  I(V) = I^n + g (V - V^n) and g = dI/dV estimated by secant from the
    //  last two accepted steps,
    //
    //      V^{n+1} = (V_src - R I^n + (Rg + a) V^n) / (1 + Rg + a).
    //
    //  This damps itself exactly where the explicit form diverged: after
    //  breakdown g is large, so Rg dominates the denominator and the step
    //  shrinks. Before breakdown g -> 0 and it reduces to the RC form above.
    //  g is clamped non-negative -- a discharge cannot have negative
    //  differential conductance to the CIRCUIT without the circuit being the
    //  thing that stabilises it, which is the case a ballast exists for.
    //
    //  There is no explicit option: it was tried, it is unusable, and keeping
    //  it would only offer a way to reproduce the divergence.
    void update
    (
        const scalar Icond,
        const scalar Cgap,
        volScalarField& ePotential
    );
};

} // End namespace Foam

#endif
