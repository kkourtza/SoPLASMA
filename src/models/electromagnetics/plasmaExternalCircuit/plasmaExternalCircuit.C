/*---------------------------------------------------------------------------*\
  SoPLASMA -- external circuit coupled to a driven electrode.
\*---------------------------------------------------------------------------*/

#include "plasmaExternalCircuit.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "boundaryRoleLibrary.H"

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

Foam::plasmaExternalCircuit::plasmaExternalCircuit(const fvMesh& mesh)
:
    mesh_(mesh),
    enabled_(false),
    type_("none"),
    electrode_(word::null),
    patchi_(-1),
    R_(0.0),
    Cext_(0.0),
    relax_(1.0),
    V_(0.0),
    started_(false),
    Vprev_(0.0),
    Iprev_(0.0),
    havePrev_(false),
    g_(0.0)
{
    // THE CIRCUIT IS DECLARED WHERE THE ELECTRODE IS, AND NOWHERE ELSE.
    //
    // It used to live in system/plasmaSimulationControls with an `electrode`
    // key naming the patch -- so the patch was named twice, once there and
    // once as `kind ballastedElectrode` in configuration/boundaries. Two
    // copies of one name is a defect even while they agree (CLAUDE.md G1),
    // and the failure mode is silent: rename the patch in one file and the
    // circuit quietly drives nothing.
    //
    // Now the user writes ONE thing, in the semantic layer:
    //
    //     cathode
    //     {
    //         kind        ballastedElectrode;
    //         circuit
    //         {
    //             type          seriesResistor;
    //             sourceVoltage table ((0 0) (5e-6 -1011));
    //             resistance    5e8;
    //         }
    //     }
    //
    // and the electrode this circuit feeds is the patch it was declared on.
    const dictionary decl(boundaryRoleLibrary::caseDeclaration(mesh_.time()));

    wordList ballasted;
    forAllConstIters(decl, iter)
    {
        if (!iter().isDict()) continue;
        const dictionary& pd = iter().dict();
        if (pd.getOrDefault<word>("kind", word::null) == "ballastedElectrode")
        {
            ballasted.append(iter().keyword());
        }
    }

    if (ballasted.empty())
    {
        return;                       // no ballasted electrode: no circuit
    }

    // ONE CIRCUIT TODAY. Several ballasted electrodes means a netlist, which
    // is stage 2 of doc/external-circuit-plan.md; refusing is better than
    // driving the first one and silently ignoring the rest.
    if (ballasted.size() > 1)
    {
        FatalErrorInFunction
            << "several patches are declared `ballastedElectrode`: "
            << ballasted << nl
            << "    Only ONE external circuit is supported today. Coupling"
            << " several electrodes through one network is a netlist -- stage"
            << " 2 of doc/external-circuit-plan.md." << nl
            << exit(FatalError);
    }

    electrode_ = ballasted[0];
    const dictionary& pd = decl.subDict(electrode_);

    if (!pd.found("circuit"))
    {
        FatalIOErrorInFunction(pd)
            << "patch `" << electrode_ << "` is declared"
            << " `kind ballastedElectrode` but carries no `circuit`"
            << " sub-dictionary." << nl
            << "    A ballasted electrode's potential is an OUTPUT of its"
            << " circuit, so without one there is nothing to compute it from."
            << nl
            << "    Add, inside the patch's block in configuration/boundaries:"
            << nl
            << "        circuit { type seriesResistor; sourceVoltage"
            << " <Function1>; resistance <Ohm>; }" << nl
            << exit(FatalIOError);
    }

    const dictionary& cd = pd.subDict("circuit");
    enabled_ = true;

    // ONE TOPOLOGY TODAY, AND UNKNOWN NAMES ARE FATAL.
    //
    // `type` is the extension point for the staged design in
    // doc/external-circuit-plan.md. Accepting an unrecognised name and
    // defaulting to a series resistor would mean a case asking for an RC
    // ballast silently got a resistive one -- the class of failure where the
    // run completes and the answer is for a different circuit.
    type_ = cd.get<word>("type");
    if (type_ != "seriesResistor" && type_ != "seriesRC")
    {
        FatalIOErrorInFunction(cd)
            << "circuit/type is `" << type_
            << "`, which is not implemented." << nl
            << "    Available: seriesResistor, seriesRC" << nl
            << "    Planned (doc/external-circuit-plan.md): seriesRLC,"
            << " currentSource, matchedRF." << nl
            << exit(FatalIOError);
    }

    R_ = cd.get<scalar>("resistance");

    // THE SHUNT CAPACITANCE, and why a real rig needs one.
    //
    // With the GAP's own capacitance alone -- 1.77e-16 F for this geometry --
    // R*C is 17.7 ns, so the electrode snaps to full voltage before any plasma
    // exists. The gap then sits at hundreds of volts with alpha*d enormous,
    // and the first avalanche overshoots by orders of magnitude before space
    // charge can screen anything. Measured 2026-09-06: n_e reached 6.3e19
    // m^-3, an ionisation degree of 2.6e-3, where the screening-arrest
    // estimate and the reference both say ~2e15.
    //
    // A REAL dc glow rig has far more capacitance across the gap than the gap
    // itself: ~100 pF/m of coax from the supply, which with R = 1e8 gives
    // tau = 10 ms. That is why a laboratory dc discharge lights GENTLY -- the
    // voltage creeps up and the gas breaks down at the lowest voltage that
    // sustains it, rather than being slammed far past Paschen.
    //
    // 10 ms is 200x longer than we can afford to simulate, so a case picks C
    // to make the rise slow against the IONISATION and SCREENING timescales
    // (both ~1.5e-10 s) while staying affordable: C = 5e-14 F gives 5 us,
    // gentle by 3e4 and comparable to the 4.5 us ion transit.
    //
    // It ADDS to the gap capacitance rather than replacing it: both are
    // physically in parallel across the electrode.
    Cext_ = (type_ == "seriesRC") ? cd.get<scalar>("capacitance") : 0.0;

    if (type_ == "seriesRC" && Cext_ <= 0)
    {
        FatalIOErrorInFunction(cd)
            << "circuit/capacitance must be positive for seriesRC; got "
            << Cext_ << "." << nl
            << "    A zero shunt capacitance IS `seriesResistor` -- use that"
            << " rather than an RC that is not one." << nl
            << exit(FatalIOError);
    }
    relax_ = cd.getOrDefault<scalar>("relaxation", 1.0);

    if (R_ <= 0)
    {
        FatalIOErrorInFunction(cd)
            << "circuit/resistance must be positive; got " << R_ << "." << nl
            << "    A zero ballast is a fixed-voltage electrode, which is what"
            << " `kind drivenElectrode` already is -- use that instead of a"
            << " circuit that does nothing." << nl
            << exit(FatalIOError);
    }

    if (relax_ <= 0 || relax_ > 1)
    {
        FatalIOErrorInFunction(cd)
            << "circuit/relaxation must be in (0, 1]; got " << relax_ << "."
            << nl
            << "    It is not needed for stability -- the update is implicit."
            << nl << exit(FatalIOError);
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

    // ON A RESTART, TAKE V FROM THE FIELD, NOT FROM THE SOURCE.
    //
    // The electrode potential is this class's own state, and it is NOT the
    // open-circuit value once current is flowing. Re-initialising from
    // V_source after a restart injects a step: at t = 1e-6 s on the Grubert
    // case, V_source is -202 V while the gap was actually at -163 V, so the
    // circuit would have jerked the electrode by 39 V at the resume.
    //
    // ePotential already carries the value, having been read back from the
    // restart directory, so the state is recoverable without checkpointing
    // anything of our own.
    V_ = source_->value(mesh_.time().value());

    if (mesh_.foundObject<volScalarField>("ePotential"))
    {
        const volScalarField& phi =
            mesh_.lookupObject<volScalarField>("ePotential");

        const fvPatchScalarField& pf = phi.boundaryField()[patchi_];

        if (pf.size())
        {
            const scalar Vfield = gAverage(pf);

            // A FRESH START has the generated `uniform 0`, which is not a
            // state to resume from; the open-circuit source value is right
            // there. Anything else is a genuine restart.
            if (mag(Vfield) > SMALL)
            {
                V_ = Vfield;
                started_ = true;

                Info<< "plasmaExternalCircuit: RESUMED from the field,"
                    << " V_electrode = " << V_ << " V" << nl
                    << "    (the open-circuit source value here would be "
                    << source_->value(mesh_.time().value())
                    << " V -- resuming from that would step the electrode)."
                    << endl;
            }
        }
    }

    Info<< "plasmaExternalCircuit: seriesResistor on patch `" << electrode_
        << "`" << nl
        << "    V_electrode(t) = V_source(t) - R*I_circuit,  R = " << R_
        << " Ohm" << (relax_ < 1 ? ", relaxation " : "")
        << (relax_ < 1 ? name(relax_) : "") << nl
        << "    THE ELECTRODE POTENTIAL IS NOW AN OUTPUT. A fixed-voltage gap"
        << " above breakdown has nothing to limit its current; with a ballast"
        << " the operating point is self-selecting." << nl
        << "    Coupling is IMPLICIT (backward Euler on the RC relation the"
        << " ballast forms with the gap), so it is unconditionally stable and"
        << " needs no relaxation." << nl
        << "    It uses the CONDUCTION current from the dischargeCurrent"
        << " diagnostic, which must therefore be enabled; the capacitive part"
        << " is carried by the RC term, not taken from I_total." << nl
        << "    Starting from the open-circuit value V = " << V_ << " V."
        << endl;

    if (Pstream::master())
    {
        const fileName dir(mesh_.time().globalPath()/"postProcessing"/"externalCircuit");
        mkDir(dir);
        file_.reset(new OFstream(dir/"circuit.csv"));
        *file_ << "# external circuit: V_electrode = V_source - R*I_circuit" << nl
               << "# R = " << R_ << " Ohm, electrode = " << electrode_ << nl
               << "# I_cond is the CONDUCTION current. The capacitive part is"
                  " carried implicitly by the" << nl
               << "# RC update, NOT taken from I_total -- using I_total here"
                  " makes V = V_src - R*I the RC" << nl
               << "# equation evaluated explicitly, which amplifies by R*C/dt"
                  " per step (885 measured)." << nl
               << "time,V_source,I_cond,V_electrode,g_dIdV" << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::plasmaExternalCircuit::update
(
    const scalar Icond,
    const scalar Cgap,
    volScalarField& ePotential
)
{
    if (!enabled_)
    {
        return;
    }

    const scalar t  = mesh_.time().value();
    const scalar dt = mesh_.time().deltaTValue();
    const scalar Vsrc = source_->value(t);

    // THE ELECTRODE SEES AN RC CIRCUIT, NOT A RESISTOR.
    //
    // The ballast is in series with the gap, and the gap has capacitance, so
    //     V = V_source - R*(I_cond + C dV/dt)
    // which is an ODE. Backward Euler on it gives the update below; it is
    // unconditionally stable, which the explicit form is emphatically not
    // (R*C/dt = 885 at dt = 1e-10 s here).
    //
    // C is the GAP capacitance the discharge-current diagnostic already
    // derives from its unit-potential solve, so it is a computed property of
    // this geometry rather than a number anybody types.
    // C_ext is in PARALLEL with the gap, so the capacitances add.
    const scalar tau = R_*(Cgap + Cext_);   // [s]
    const scalar a   = tau/max(dt, SMALL);

    // THE PLASMA IS LINEARISED IMPLICITLY, via a secant estimate of its
    // differential conductance g = dI/dV from the last two accepted steps.
    //
    // Without this the circuit is implicit only in its own capacitance, and
    // post-breakdown that is not enough: the current the circuit sees is one
    // step old, the discharge outruns it, and the ballast demands a voltage
    // the source cannot supply. Including g makes the update contract exactly
    // where the discharge is stiff, because Rg then dominates the denominator.
    if (havePrev_)
    {
        const scalar dV = V_ - Vprev_;

        // A negligible voltage change gives no information about dI/dV; reuse
        // the previous estimate rather than dividing by nothing. The threshold
        // is relative to the working voltage, not absolute, so it means the
        // same thing on a 100 V and a 1000 V discharge.
        if (mag(dV) > 1e-6*max(mag(V_), scalar(1)))
        {
            const scalar gNew = (Icond - Iprev_)/dV;

            // THE MAGNITUDE, NOT THE SIGNED VALUE, AND NOT CLAMPED TO ZERO.
            //
            // dI/dV IS GENUINELY NEGATIVE through breakdown, and that is the
            // physics rather than a sampling artefact: the current rises while
            // the ballast drags the gap voltage back up, so dI < 0 while
            // dV > 0. It is the NEGATIVE DIFFERENTIAL RESISTANCE of the glow --
            // the very thing a ballast exists to stabilise. Measured
            // 2026-09-06: g went negative exactly as j started to climb.
            //
            // Using it signed is wrong: 1 + Rg + a then shrinks and can change
            // sign, which inverts the update -- an amplifier. Clamping it to
            // zero is also wrong, and was the first attempt: it discards the
            // damping precisely where the discharge is stiffest, leaving the
            // RC-only form that already failed.
            //
            // The magnitude is right because of what the update then is:
            //
            //     V^{n+1} = V^n + (T - V^n)/(1 + R|g| + a),
            //     T = V_src - R I^n
            //
            // a damped step towards the explicit target, with a factor in
            // (0, 1] for ANY sign of g, and the fixed point V = T preserved
            // exactly. Steeper characteristic -> smaller step, which is the
            // behaviour wanted, and it is a relaxed Newton rather than a
            // heuristic.
            g_ = mag(gNew);
        }
    }

    Vprev_ = V_;
    Iprev_ = Icond;
    havePrev_ = true;

    const scalar Rg = R_*g_;   // g_ is already a magnitude

    const scalar Vtarget =
        (Vsrc - R_*Icond + (Rg + a)*V_)/(1.0 + Rg + a);

    // Optional extra damping. Not needed for stability; 1 by default.
    V_ = started_ ? (V_ + relax_*(Vtarget - V_)) : Vtarget;
    started_ = true;

    ePotential.boundaryFieldRef()[patchi_] == V_;

    if (file_.valid() && Pstream::master())
    {
        *file_ << t << ',' << Vsrc << ',' << Icond << ',' << V_
               << ',' << g_ << endl;
    }
}


// ************************************************************************* //
