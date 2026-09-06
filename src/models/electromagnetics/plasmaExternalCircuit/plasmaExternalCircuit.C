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
    compliance_(0.0),
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
    // once as `kind ballastedElectrode` (or `currentDrivenElectrode`) in
    // configuration/boundaries. Two
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
        // TWO ROLES, because the physical difference is WHICH QUANTITY THE
        // SUPPLY IMPOSES -- the same axis as driven / grounded / floating, not
        // a topology detail. `ballastedElectrode` sets a voltage and lets the
        // current find itself; `currentDrivenElectrode` does the reverse.
        // Which one a case wants is the user's choice and cannot be derived:
        // the operating current is usually the ANSWER, not an input.
        const word kind = pd.getOrDefault<word>("kind", word::null);

        if (kind == "ballastedElectrode" || kind == "currentDrivenElectrode")
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
            << "several patches are declared with a circuit: "
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
            << " `kind ballastedElectrode`/`currentDrivenElectrode` but carries no `circuit`"
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
    if (type_ != "seriesResistor" && type_ != "seriesRC"
     && type_ != "currentSource")
    {
        FatalIOErrorInFunction(cd)
            << "circuit/type is `" << type_
            << "`, which is not implemented." << nl
            << "    Available: seriesResistor, seriesRC, currentSource" << nl
            << "    Planned (doc/external-circuit-plan.md): seriesRLC,"
            << " matchedRF." << nl
            << exit(FatalIOError);
    }

    const bool isSource = (type_ == "currentSource");

    // THE ROLE AND THE TOPOLOGY MUST AGREE, both ways.
    //
    // `kind` says which quantity the supply imposes; `type` picks the topology
    // within that. A mismatch means the case has asked for two different
    // instruments, and guessing which one it meant is exactly the class of
    // silent wrong-answer this whole model exists to remove.
    const word kind = pd.getOrDefault<word>("kind", word::null);

    if (isSource && kind == "ballastedElectrode")
    {
        FatalIOErrorInFunction(cd)
            << "`type currentSource` on a `ballastedElectrode`." << nl
            << "    `ballastedElectrode` means the SUPPLY VOLTAGE is what you"
            << " set and the current is" << nl
            << "    an output. A current source imposes the current instead,"
            << " so it belongs on" << nl
            << "    `kind currentDrivenElectrode`." << nl
            << "    If what you want is a ballast, use `type seriesRC` with a"
            << " sourceVoltage." << nl
            << exit(FatalIOError);
    }

    if (!isSource && kind == "currentDrivenElectrode")
    {
        FatalIOErrorInFunction(cd)
            << "`type " << type_ << "` on a `currentDrivenElectrode`." << nl
            << "    `currentDrivenElectrode` means the CURRENT is what you"
            << " set. `" << type_ << "` sets a" << nl
            << "    voltage through a ballast and lets the current find"
            << " itself, which is" << nl
            << "    `kind ballastedElectrode`." << nl
            << exit(FatalIOError);
    }

    // A CURRENT SOURCE HAS NO BALLAST, so `resistance` is not merely optional
    // here -- supplying it means the case is asking for something this type
    // does not do, and silently ignoring it is how a case ends up believing it
    // has a ballast it does not have.
    if (isSource)
    {
        if (cd.found("resistance"))
        {
            FatalIOErrorInFunction(cd)
                << "circuit/resistance is set on a `currentSource`, which has"
                << " no ballast." << nl
                << "    A current source IS the R -> infinity limit of a"
                << " ballast: it regulates the" << nl
                << "    current directly, so R has no meaning. Either remove"
                << " `resistance`, or use" << nl
                << "    `type seriesRC` with that R and"
                << " `sourceVoltage = V_gap + R*I_set`." << nl
                << exit(FatalIOError);
        }
        R_ = 0.0;
    }
    else
    {
        R_ = cd.get<scalar>("resistance");
    }

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
    // For a current source C is OPTIONAL but strongly recommended: it is what
    // sets the pre-ignition ramp rate, dV/dt = I_set/C. See update().
    Cext_ = (type_ == "seriesRC")
          ? cd.get<scalar>("capacitance")
          : (isSource ? cd.getOrDefault<scalar>("capacitance", 0.0) : 0.0);

    if (type_ == "seriesRC" && Cext_ <= 0)
    {
        FatalIOErrorInFunction(cd)
            << "circuit/capacitance must be positive for seriesRC; got "
            << Cext_ << "." << nl
            << "    A zero shunt capacitance IS `seriesResistor` -- use that"
            << " rather than an RC that is not one." << nl
            << exit(FatalIOError);
    }
    if (isSource && Cext_ < 0)
    {
        FatalIOErrorInFunction(cd)
            << "circuit/capacitance must not be negative; got " << Cext_
            << "." << nl << exit(FatalIOError);
    }
    relax_ = cd.getOrDefault<scalar>("relaxation", 1.0);

    if (!isSource && R_ <= 0)
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

    if (isSource)
    {
        // THE COMPLIANCE VOLTAGE IS REQUIRED, and it is what fixes the sign.
        //
        // Every real current supply has a voltage rail it cannot exceed, and
        // that rail is also what makes the pre-ignition phase well posed: with
        // no plasma there is no voltage at which the demanded current flows,
        // so an unbounded regulator would run away. The rail is SIGNED,
        // because the user knows their supply's polarity, and `setCurrent` is
        // then a MAGNITUDE -- which removes the one sign convention nobody can
        // be expected to guess (I_cond is negative on a negative electrode,
        // verified against postProcessing/externalCircuit/circuit.csv on
        // 2026-09-06).
        compliance_ = cd.get<scalar>("compliance");

        if (compliance_ == 0)
        {
            FatalIOErrorInFunction(cd)
                << "circuit/compliance is zero, so the source can never drive"
                << " any current." << nl
                << "    It is the supply's voltage rail, SIGNED: negative for"
                << " a cathode." << nl
                << exit(FatalIOError);
        }

        setCurrent_ = Function1<scalar>::New("setCurrent", cd);

        // A magnitude, so a negative entry is a sign error the user made once
        // and would otherwise chase for an afternoon.
        const scalar I0 = setCurrent_->value(mesh_.time().value());
        if (I0 < 0)
        {
            FatalIOErrorInFunction(cd)
                << "circuit/setCurrent is negative (" << I0 << " A at t = "
                << mesh_.time().value() << ")." << nl
                << "    setCurrent is a MAGNITUDE in amperes; the POLARITY"
                << " comes from the sign of" << nl
                << "    `compliance` (" << compliance_ << " V here)." << nl
                << exit(FatalIOError);
        }
    }
    else
    {
        source_ = Function1<scalar>::New("sourceVoltage", cd);
    }

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
    // A CURRENT SOURCE HAS NO `sourceVoltage`, so there is no open-circuit
    // value to start from: it starts at ZERO and charges the electrode
    // capacitance at I_set. Dereferencing `source_` here aborted on the
    // model's FIRST EVER RUN, 2026-09-06 -- "unallocated autoPtr of type
    // Function1<double>" from this constructor. The guard below is not
    // cosmetic: three separate places assumed a voltage source exists.
    V_ = (type_ == "currentSource")
       ? 0.0
       : source_->value(mesh_.time().value());

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
                    << " V_electrode = " << V_ << " V";

                if (type_ != "currentSource")
                {
                    Info<< nl
                        << "    (the open-circuit source value here would be "
                        << source_->value(mesh_.time().value())
                        << " V -- resuming from that would step the"
                        << " electrode).";
                }
                Info<< endl;
            }
        }
    }

    if (type_ == "currentSource")
    {
        Info<< "plasmaExternalCircuit: currentSource on patch `" << electrode_
            << "`" << nl
            << "    I_set = "
            << setCurrent_->value(mesh_.time().value()) << " A (magnitude),"
            << " compliance " << compliance_ << " V, C_ext = " << Cext_
            << " F" << nl
            << "    THE ELECTRODE POTENTIAL IS AN OUTPUT, and so is the gap"
            << " voltage: the current is what is imposed. This is the"
            << " R -> infinity limit of a ballast, and it is the STRONGEST"
            << " pin on the operating point available." << nl
            << "    BEFORE IGNITION there is no voltage at which I_set flows,"
            << " so the source does what a real one does: it charges the"
            << " electrode capacitance at constant current, dV/dt = I_set/C"
            << " = "
            << (Cext_ > 0
                 ? name(setCurrent_->value(mesh_.time().value())/Cext_)
                 : word("(C_gap only -- set `capacitance` to control it)"))
            << " V/s, until the gas breaks down." << nl
            << "    AFTER IGNITION it regulates: the update is a damped Newton"
            << " step on I(V) using the secant conductance, whose fixed point"
            << " is I_cond = I_set exactly." << nl
            << "    THE COMPLIANCE RAIL BOUNDS IT BOTH WAYS. If the discharge"
            << " cannot pass I_set at any voltage up to the rail, the source"
            << " sits at the rail -- which is a real supply's behaviour and"
            << " NOT a solver failure. Watch for it in circuit.csv." << nl
            << "    It uses the CONDUCTION current from the dischargeCurrent"
            << " diagnostic, which must therefore be enabled." << nl
            << "    Starting from V = " << V_ << " V." << endl;
    }
    else
    {
    Info<< "plasmaExternalCircuit: " << type_ << " on patch `" << electrode_
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
    }

    if (Pstream::master())
    {
        const fileName dir(mesh_.time().globalPath()/"postProcessing"/"externalCircuit");
        mkDir(dir);
        file_.reset(new OFstream(dir/"circuit.csv"));
        if (type_ == "currentSource")
        {
            *file_
                << "# external circuit: currentSource, I_cond regulated to"
                   " I_set" << nl
                << "# I_set = "
                << setCurrent_->value(mesh_.time().value())
                << " A (magnitude), compliance = " << compliance_
                << " V, C_ext = " << Cext_ << " F, electrode = "
                << electrode_ << nl
                << "# V_source is the SET CURRENT here, signed by the"
                   " compliance rail, so the column" << nl
                << "# reads as the target and I_cond as what was achieved."
                   " V_electrode AT the rail means" << nl
                << "# the source ran out of compliance -- real behaviour, not"
                   " a failure." << nl
                << "time,I_set,I_cond,V_electrode,g_dIdV" << endl;
        }
        else
        {
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

    // ------------------------------------------------------------------
    // CURRENT SOURCE: the current is imposed, the voltage is the unknown.
    // ------------------------------------------------------------------
    //
    // A current source regulates the TERMINAL current, and the terminal
    // current is conduction plus the displacement drawn by the electrode
    // capacitance:
    //
    //     I_set = I_cond + C dV/dt                                     (1)
    //
    // Backward Euler on (1), with the plasma linearised by the same secant
    // conductance g = dI/dV the ballast uses, gives ONE update that is correct
    // in both regimes without a switch between them:
    //
    //     dV = dt (I_set - I_cond) / (C + |g| dt)                       (2)
    //
    // BEFORE IGNITION g -> 0 and (2) becomes dV = dt I_set / C: the source
    // charges the electrode capacitance at constant current, dV/dt = I_set/C.
    // That is EXACTLY what a real current-limited supply does into a
    // capacitor, and it means the pre-ignition voltage ramp is a consequence
    // of the circuit rather than a ramp anybody has to write. The ramp rate is
    // a design knob: it is set by C alone.
    //
    // AFTER IGNITION g is large, (2) becomes dV = (I_set - I_cond)/|g|, a
    // damped Newton step on I(V). Its fixed point is I_cond = I_set exactly,
    // for any g.
    //
    // |g| RATHER THAN g, for the reason documented at the secant estimate
    // below: dI/dV is genuinely negative through a glow's negative
    // differential resistance, and a signed g would inflate the step and
    // invert its direction precisely where the discharge is stiffest.
    //
    // WHY THIS IS THE STRONGEST AVAILABLE PIN ON THE OPERATING POINT.
    // A ballast's short-circuit current is I_sc/I_op = 1 + V_gap/(R I_op), so
    // it only pins the current in the limit R -> infinity with
    // V_src = V_gap + R I_op -- which IS this model. Measured 2026-09-06 on
    // the Grubert case: 5.89x the operating point at R = 1e8, 1.10x at
    // R = 5e9. A current source is the 1.00x end of that sweep.
    if (type_ == "currentSource")
    {
        const scalar Iset = sign(compliance_)*setCurrent_->value(t);
        const scalar C    = Cgap + Cext_;

        updateSecantConductance(Icond);

        Vprev_ = V_;
        Iprev_ = Icond;
        havePrev_ = true;

        const scalar dV = dt*(Iset - Icond)/max(C + g_*dt, VSMALL);

        V_ = started_ ? (V_ + relax_*dV) : (V_ + dV);
        started_ = true;

        // THE COMPLIANCE RAIL, clamped on BOTH sides.
        //
        // Above the rail the supply physically cannot go. Below zero it would
        // have to reverse polarity, which a single-quadrant supply cannot do
        // either -- and an unclamped regulator does try to, because
        // overshooting the set current asks for a voltage of the other sign.
        V_ = min(max(V_, min(compliance_, scalar(0))),
                 max(compliance_, scalar(0)));

        ePotential.boundaryFieldRef()[patchi_] == V_;

        if (file_.valid() && Pstream::master())
        {
            *file_ << t << ',' << Iset << ',' << Icond << ',' << V_
                   << ',' << g_ << endl;
        }

        return;
    }

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
    updateSecantConductance(Icond);

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


void Foam::plasmaExternalCircuit::updateSecantConductance(const scalar Icond)
{
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
}


// ************************************************************************* //
