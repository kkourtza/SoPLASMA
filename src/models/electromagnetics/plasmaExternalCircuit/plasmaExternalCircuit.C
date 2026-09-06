/*---------------------------------------------------------------------------*\
  SoPLASMA -- external circuit coupled to a driven electrode.
\*---------------------------------------------------------------------------*/

#include "plasmaExternalCircuit.H"
#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

Foam::plasmaExternalCircuit::plasmaExternalCircuit
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    enabled_(false),
    type_("none"),
    electrode_(word::null),
    patchi_(-1),
    R_(0.0),
    relax_(1.0),
    V_(0.0),
    started_(false)
{
    if (!dict.found("externalCircuit"))
    {
        return;
    }

    const dictionary& cd = dict.subDict("externalCircuit");

    enabled_ = cd.getOrDefault<Switch>("enabled", true);
    if (!enabled_)
    {
        return;
    }

    // ONE TOPOLOGY TODAY, AND UNKNOWN NAMES ARE FATAL.
    //
    // `type` is the extension point for the staged design in
    // doc/external-circuit-plan.md. Accepting an unrecognised name and
    // defaulting to a series resistor would mean a case asking for an RC
    // ballast silently got a resistive one -- the class of failure where the
    // run completes and the answer is for a different circuit.
    type_ = cd.get<word>("type");
    if (type_ != "seriesResistor")
    {
        FatalIOErrorInFunction(cd)
            << "externalCircuit/type is `" << type_
            << "`, which is not implemented." << nl
            << "    Available: seriesResistor" << nl
            << "    Planned (doc/external-circuit-plan.md): seriesRC,"
            << " seriesRLC, currentSource, matchedRF, netlist." << nl
            << exit(FatalIOError);
    }

    electrode_ = cd.get<word>("electrode");
    R_ = cd.get<scalar>("resistance");
    relax_ = cd.getOrDefault<scalar>("relaxation", 1.0);

    if (R_ <= 0)
    {
        FatalIOErrorInFunction(cd)
            << "externalCircuit/resistance must be positive; got " << R_
            << "." << nl
            << "    A zero ballast is a fixed-voltage electrode, which is what"
            << " a drivenElectrode already is -- use that instead of a circuit"
            << " that does nothing." << nl
            << exit(FatalIOError);
    }

    if (relax_ <= 0 || relax_ > 1)
    {
        FatalIOErrorInFunction(cd)
            << "externalCircuit/relaxation must be in (0, 1]; got " << relax_
            << "." << nl << exit(FatalIOError);
    }

    source_ = Function1<scalar>::New("sourceVoltage", cd);

    patchi_ = mesh_.boundaryMesh().findPatchID(electrode_);
    if (patchi_ < 0)
    {
        FatalIOErrorInFunction(cd)
            << "externalCircuit/electrode is `" << electrode_
            << "`, which is not a patch on this mesh." << nl
            << "    Patches: " << mesh_.boundaryMesh().names() << nl
            << exit(FatalIOError);
    }

    // The source at t = 0, so the first step starts from the open-circuit
    // value rather than from zero -- before any current is measured the
    // circuit carries none, and V = V_source is then exactly right.
    V_ = source_->value(mesh_.time().value());

    Info<< "plasmaExternalCircuit: seriesResistor on patch `" << electrode_
        << "`" << nl
        << "    V_electrode(t) = V_source(t) - R*I_circuit,  R = " << R_
        << " Ohm" << (relax_ < 1 ? ", relaxation " : "")
        << (relax_ < 1 ? name(relax_) : "") << nl
        << "    THE ELECTRODE POTENTIAL IS NOW AN OUTPUT. A fixed-voltage gap"
        << " above breakdown has nothing to limit its current; with a ballast"
        << " the operating point is self-selecting." << nl
        << "    I_circuit is Sato's I_total from the dischargeCurrent"
        << " diagnostic, which must therefore be enabled." << nl
        << "    Starting from the open-circuit value V = " << V_ << " V."
        << endl;

    if (Pstream::master())
    {
        const fileName dir(mesh_.time().globalPath()/"postProcessing"/"externalCircuit");
        mkDir(dir);
        file_.reset(new OFstream(dir/"circuit.csv"));
        *file_ << "# external circuit: V_electrode = V_source - R*I_circuit" << nl
               << "# R = " << R_ << " Ohm, electrode = " << electrode_ << nl
               << "# I_circuit is Sato's I_total (conduction + displacement),"
                  " LAGGED by one step:" << nl
               << "#   it is measured after the step it belongs to, so the"
                  " potential it sets applies to the NEXT step." << nl
               << "time,V_source,I_circuit,V_electrode" << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::plasmaExternalCircuit::update
(
    const scalar Icircuit,
    volScalarField& ePotential
)
{
    if (!enabled_)
    {
        return;
    }

    const scalar t = mesh_.time().value();
    const scalar Vsrc = source_->value(t);

    // Ohm's law across the ballast. The sign works out without a special case:
    // a cathode at negative potential draws negative current, and
    // -1011 - 5e8*(-1.02e-6) = -500 raises the electrode towards zero, which
    // is the direction a ballast must act.
    const scalar Vtarget = Vsrc - R_*Icircuit;

    // Under-relax, because the current is LAGGED: it was measured after the
    // step that produced it, so this is an explicit fixed-point iteration in
    // (V, I) with one step of delay. A large R makes it stiff.
    V_ = started_ ? (V_ + relax_*(Vtarget - V_)) : Vtarget;
    started_ = true;

    ePotential.boundaryFieldRef()[patchi_] == V_;

    if (file_.valid() && Pstream::master())
    {
        *file_ << t << ',' << Vsrc << ',' << Icircuit << ',' << V_ << endl;
    }
}


// ************************************************************************* //
