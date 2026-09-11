/*---------------------------------------------------------------------------*\
  File: testEmission.C
  Part of: SoPLASMA

Application
    testEmission

Description
    Check every emission model against its CLOSED FORM, evaluated here
    independently of the model's own code.

    This is a UNIT test with ground truth, not a plausibility check: the
    Fowler-Nordheim and Richardson-Dushman currents are analytic functions of
    (phi, F, T), so agreement is exact to round-off and any disagreement is a
    defect. Both are EXPONENTIAL in their arguments, which is why a unit test
    is worth more here than anywhere else -- a sign slip or a units error does
    not shift the answer, it changes it by decades.

    TWO LEVELS, deliberately:

      1. THE CLOSED FORMS, evaluated in this file, checked against physical
         invariants and against literature orders of magnitude.
      2. THE MODEL CLASSES THEMSELVES, constructed through
         emissionModel::New on a real patch of a real mesh with prescribed
         Emag and T_gas, and compared face-by-face against (1).

    Level 2 is the one that matters: level 1 only says the formulae written
    HERE are right, which is no evidence at all about the code under test.
    Needs a case with a mesh.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "plasmaConstants.H"
#include "IOmanip.H"
#include "fvMesh.H"
#include "volFields.H"
#include "emissionModel.H"
#include "electronInducedSEE.H"

using namespace Foam;

// The closed forms, written out INDEPENDENTLY of the models. Constants from
// Forbes & Deane, Proc. R. Soc. A 463 (2007) 2907 (a, b) and CRC 95th (A0).
static scalar fnCurrent
(
    const scalar phi,
    const scalar F,
    const bool schottkyNordheim
)
{
    const scalar a = 1.541434e-6;
    const scalar b = 6.830890e9;
    const scalar c = 1.439964e-9;

    scalar v = 1.0;

    if (schottkyNordheim)
    {
        const scalar f = c*F/sqr(phi);
        v = (f < 1.0 && f > SMALL)
          ? max(0.0, 1.0 - f + (f/6.0)*Foam::log(f))
          : (f >= 1.0 ? 0.0 : 1.0);
    }

    return (a/phi)*sqr(F)*Foam::exp(-v*b*Foam::pow(phi, 1.5)/F);
}


static scalar vaughan
(
    const scalar E,
    const scalar dmax,
    const scalar Emax,
    const scalar E0
)
{
    if (E <= E0) return 0.0;
    const scalar v = (E - E0)/(Emax - E0);
    if (v <= 0) return 0.0;
    if (v > 3.6) return dmax*1.125*Foam::pow(v, -0.35);
    const scalar k = (v < 1.0) ? 0.56 : 0.25;
    return dmax*Foam::pow(v*Foam::exp(1.0 - v), k);
}


static scalar rdCurrent
(
    const scalar phi,
    const scalar T,
    const scalar lambdaR,
    const scalar F
)
{
    const scalar A0 = 1.201735e6;
    const scalar kB_eV = 8.617333262e-5;          // eV/K
    const scalar dphi = 3.794686e-5*Foam::sqrt(max(F, 0.0));

    return lambdaR*A0*sqr(T)*Foam::exp(-max(0.0, phi - dphi)/(kB_eV*T));
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addNote("Emission models against their closed forms");
    argList::addOption("patch", "name", "patch to evaluate on (default: the first wall)");
    #include "setRootCase.H"
    #include "createTime.H"

    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label nFail = 0;

    const scalar e = constant::plasma::eCharge.value();

    Info<< nl << "=== FOWLER-NORDHEIM  J = (a/phi) F^2 exp(-v b phi^3/2 / F)"
        << nl << nl
        << "    phi[eV]      F[V/m]   barrier            J[A/m2]"
           "     flux[1/m2/s]" << nl;

    for (const scalar phi : {4.28, 4.65, 5.65})
    {
        for (const scalar F : {1e9, 3e9, 5e9})
        {
            for (const bool sn : {false, true})
            {
                const scalar J = fnCurrent(phi, F, sn);

                Info<< "    " << setw(7) << phi
                    << setw(12) << F
                    << setw(20) << (sn ? "schottkyNordheim" : "elementary")
                    << setw(14) << J
                    << setw(16) << J/e << nl;
            }
        }
    }

    // THE PHYSICAL CHECKS, each of which a sign or units error breaks.
    Info<< nl << "    physical invariants:" << nl;
    {
        // 1. MONOTONIC AND STEEP in the field.
        const scalar J1 = fnCurrent(4.65, 1e9, true);
        const scalar J2 = fnCurrent(4.65, 3e9, true);
        const bool ok1 = (J2 > J1) && (J2/J1 > 1e3);
        Info<< "      steeply increasing with F      "
            << (ok1 ? "PASS" : "FAIL") << "   J(3e9)/J(1e9) = "
            << J2/max(J1, VSMALL) << nl;
        if (!ok1) ++nFail;

        // 2. DECREASING in the work function.
        const bool ok2 = fnCurrent(5.65, 3e9, true) < fnCurrent(4.28, 3e9, true);
        Info<< "      decreasing with phi            "
            << (ok2 ? "PASS" : "FAIL") << nl;
        if (!ok2) ++nFail;

        // 3. The image-charge correction RAISES the current, always.
        const bool ok3 = fnCurrent(4.65, 3e9, true) > fnCurrent(4.65, 3e9, false);
        Info<< "      Schottky-Nordheim > elementary "
            << (ok3 ? "PASS" : "FAIL") << "   ratio = "
            << fnCurrent(4.65, 3e9, true)/fnCurrent(4.65, 3e9, false) << nl;
        if (!ok3) ++nFail;

        // 4. ORDER OF MAGNITUDE against the literature. Field emission from a
        //    clean metal becomes measurable (~1-1e4 A/m2) in the low 1e9 V/m
        //    range; that is the classical FN onset and it is what makes the
        //    formula recognisable.
        const scalar Jref = fnCurrent(4.5, 3e9, false);
        const bool ok4 = (Jref > 1.0) && (Jref < 1e6);
        Info<< "      onset near 3e9 V/m is 1-1e6 A/m2  "
            << (ok4 ? "PASS" : "FAIL") << "   J = " << Jref << " A/m2" << nl;
        if (!ok4) ++nFail;
    }

    Info<< nl << "=== RICHARDSON-DUSHMAN  J = lambda A0 T^2 exp(-phi/kT)"
        << nl << nl
        << "    phi[eV]       T[K]     F[V/m]            J[A/m2]"
           "     flux[1/m2/s]" << nl;

    for (const scalar T : {1500.0, 2000.0, 2500.0})
    {
        for (const scalar F : {0.0, 1e7})
        {
            const scalar J = rdCurrent(4.55, T, 1.0, F);

            Info<< "    " << setw(7) << 4.55
                << setw(11) << T
                << setw(11) << F
                << setw(19) << J
                << setw(16) << J/e << nl;
        }
    }

    Info<< nl << "    physical invariants:" << nl;
    {
        // 5. Tungsten at 2500 K emits of order 0.1-10 A/cm2 -- the classical
        //    thermionic-cathode figure, and the reason W is used for one.
        const scalar J = rdCurrent(4.55, 2500.0, 1.0, 0.0);
        const scalar JcmSq = J/1e4;
        const bool ok5 = (JcmSq > 0.01) && (JcmSq < 100.0);
        Info<< "      W at 2500 K is 0.01-100 A/cm2  "
            << (ok5 ? "PASS" : "FAIL") << "   J = " << JcmSq << " A/cm2" << nl;
        if (!ok5) ++nFail;

        // 6. The Schottky effect RAISES it, and by a factor the exponent
        //    predicts: exp(dphi/kT) with dphi = 3.7947e-5 sqrt(F).
        const scalar Jno = rdCurrent(4.55, 2000.0, 1.0, 0.0);
        const scalar Jsc = rdCurrent(4.55, 2000.0, 1.0, 1e7);
        const scalar kT = 8.617333262e-5*2000.0;
        const scalar predicted = Foam::exp(3.794686e-5*Foam::sqrt(1e7)/kT);
        const scalar measured = Jsc/Jno;
        const bool ok6 = mag(measured - predicted)/predicted < 1e-10;
        Info<< "      Schottky factor = exp(dphi/kT) "
            << (ok6 ? "PASS" : "FAIL")
            << "   measured " << measured << " vs " << predicted << nl;
        if (!ok6) ++nFail;

        // 7. Cold walls emit NOTHING measurable: at 300 K, W gives ~1e-72.
        const scalar Jcold = rdCurrent(4.55, 300.0, 1.0, 0.0);
        const bool ok7 = Jcold < 1e-30;
        Info<< "      cold wall (300 K) is negligible "
            << (ok7 ? "PASS" : "FAIL") << "   J = " << Jcold << " A/m2" << nl;
        if (!ok7) ++nFail;
    }

    // ----------------------------------------------------------------------
    // LEVEL 2: THE CLASSES, on a real patch.
    // ----------------------------------------------------------------------
    Info<< nl << "=== THE MODEL CLASSES, on a real patch" << nl << nl;

    label patchi = -1;
    if (args.found("patch"))
    {
        patchi = mesh.boundaryMesh().findPatchID(args.get<word>("patch"));
    }
    else
    {
        forAll(mesh.boundary(), pi)
        {
            if (!mesh.boundary()[pi].coupled() && mesh.boundary()[pi].size())
            {
                patchi = pi;
                break;
            }
        }
    }

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "No usable patch. Available: " << mesh.boundaryMesh().names()
            << nl << exit(FatalError);
    }

    const fvPatch& p = mesh.boundary()[patchi];

    Info<< "    patch `" << p.name() << "`, " << p.size() << " faces" << nl
        << nl;

    // Prescribed drivers, uniform so the closed form is a single number the
    // model must reproduce on EVERY face.
    const scalar Etest = 3e9;      // V/m
    const scalar Ttest = 2000.0;   // K
    const scalar phiTest = 4.65;   // eV, copper

    volScalarField Emag
    (
        IOobject("Emag", runTime.timeName(), mesh, IOobject::NO_READ,
                 IOobject::NO_WRITE, IOobject::REGISTER),
        mesh,
        dimensionedScalar("Emag", dimensionSet(1, 1, -3, 0, 0, -1, 0), Etest)
    );
    Emag.boundaryFieldRef() == Etest;

    volScalarField Tgas
    (
        IOobject("T_gas", runTime.timeName(), mesh, IOobject::NO_READ,
                 IOobject::NO_WRITE, IOobject::REGISTER),
        mesh,
        dimensionedScalar("T_gas", dimTemperature, Ttest)
    );
    Tgas.boundaryFieldRef() == Ttest;

    auto check = [&](const word& type, const dictionary& d, const scalar want)
    {
        autoPtr<emissionModel> m(emissionModel::New(type, p, d));

        const scalarField f(m->emittedFlux());

        const scalar got = gMax(f);
        const scalar spread = gMax(f) - gMin(f);
        const scalar rel = mag(got - want)/max(mag(want), VSMALL);

        // TOLERANCE 1e-8, and the reason is worth stating because a tighter
        // one FAILED for a non-reason. The closed forms here carry their own
        // literal for k_B/e (8.617333262e-5) while the models compute it from
        // constant::plasma::kappaBoltzmann/eCharge; the two agree to 1.7e-11,
        // and the Richardson exponent phi/kT ~ 27 amplifies that to ~4.5e-10
        // in the current. That is round-off between two constant literals, not
        // a defect.
        //
        // The test stays discriminating at 1e-8: a sign slip, a units error or
        // a wrong constant in either of these models changes the answer by
        // DECADES, not by parts in 1e10 -- both are exponential in their
        // arguments, which is exactly why they are worth unit-testing.
        const bool ok = (rel < 1e-8) && (spread <= VSMALL*max(got, 1.0));

        Info<< "    " << setw(20) << type
            << "  model " << setw(14) << got
            << "  closed form " << setw(14) << want
            << "  rel " << setw(10) << rel
            << "  " << (ok ? "PASS" : "FAIL") << nl;

        if (!ok)
        {
            ++nFail;
            if (spread > VSMALL*max(got, 1.0))
            {
                Info<< "        NON-UNIFORM over the patch (spread "
                    << spread << ") on a uniform driver" << nl;
            }
        }
        m->report(Info);
    };

    {
        dictionary d;
        d.add("workFunction", phiTest);
        d.add("barrier", word("schottkyNordheim"));
        check("fieldEmission", d, fnCurrent(phiTest, Etest, true)/e);
    }
    {
        dictionary d;
        d.add("workFunction", phiTest);
        d.add("barrier", word("elementary"));
        check("fieldEmission", d, fnCurrent(phiTest, Etest, false)/e);
    }
    {
        dictionary d;
        d.add("workFunction", phiTest);
        d.add("schottky", Switch(false));
        check("thermionicEmission", d, rdCurrent(phiTest, Ttest, 1.0, 0.0)/e);
    }
    {
        dictionary d;
        d.add("workFunction", phiTest);
        d.add("schottky", Switch(true));
        check("thermionicEmission", d, rdCurrent(phiTest, Ttest, 1.0, Etest)/e);
    }
    {
        // THE MATERIAL ROUTE: `material copper` must give the same answer as
        // `workFunction 4.65`, since that is what the library says copper is.
        dictionary d;
        d.add("material", word("copper"));
        d.add("barrier", word("elementary"));
        check("fieldEmission", d, fnCurrent(4.65, Etest, false)/e);
    }

    // ----------------------------------------------------------------------
    Info<< nl << "=== VAUGHAN yield delta(E), delta_max = 1.3 at 600 eV" << nl
        << nl << "        E[eV]       delta   note" << nl;

    for (const scalar E : {5.0, 12.5, 50.0, 300.0, 600.0, 1200.0, 3000.0})
    {
        const scalar d = vaughan(E, 1.3, 600.0, 12.5);

        Info<< "    " << setw(9) << E << setw(12) << d << "   "
            << (
                   E <= 12.5      ? "below threshold: NOTHING emitted"
                 : mag(E - 600.0) < 1 ? "at the peak"
                 : (d > 1.0      ? "delta > 1: emits more than it collects"
                                 : "")
               ) << nl;
    }

    Info<< nl << "    physical invariants:" << nl;
    {
        const bool a = vaughan(5.0, 1.3, 600.0, 12.5) == 0.0;
        Info<< "      zero below the threshold        " << (a ? "PASS" : "FAIL")
            << "   <- the branch an ATMOSPHERIC discharge lives in" << nl;
        if (!a) ++nFail;

        const scalar peak = vaughan(600.0, 1.3, 600.0, 12.5);
        const bool b = mag(peak - 1.3) < 1e-12;
        Info<< "      peak equals delta_max exactly   " << (b ? "PASS" : "FAIL")
            << "   delta(Emax) = " << peak << nl;
        if (!b) ++nFail;

        const bool c = vaughan(3000.0, 1.3, 600.0, 12.5)
                     < vaughan(600.0, 1.3, 600.0, 12.5);
        Info<< "      falls beyond the peak           " << (c ? "PASS" : "FAIL")
            << nl;
        if (!c) ++nFail;

        // The two crossings of delta = 1 are what decide the SIGN of the charge
        // an insulator accumulates, so their existence is the physically
        // meaningful feature of the curve.
        scalar lo = 0, hi = 0;
        for (scalar E = 12.6; E < 5000.0; E *= 1.002)
        {
            const scalar d = vaughan(E, 1.3, 600.0, 12.5);
            if (!lo && d >= 1.0) lo = E;
            if (lo && d < 1.0 && E > 600.0) { hi = E; break; }
        }
        const bool dd = (lo > 12.5) && (hi > 600.0);
        Info<< "      two crossings of delta = 1      " << (dd ? "PASS" : "FAIL")
            << "   at " << lo << " and " << hi << " eV" << nl;
        if (!dd) ++nFail;
    }

    // THE CLASS, against the same closed form.
    {
        using namespace emissionModels;

        bool ok = true;
        for (const scalar E : {5.0, 50.0, 600.0, 3000.0})
        {
            const scalar want = vaughan(E, 1.3, 600.0, 12.5);
            const scalar got = electronInducedSEE::vaughanYield(E, 1.3, 600.0, 12.5);
            if (mag(got - want) > 1e-14) ok = false;
        }
        Info<< "      the CLASS matches the curve     " << (ok ? "PASS" : "FAIL")
            << nl;
        if (!ok) ++nFail;
    }

    Info<< nl << (nFail ? "FAILED" : "PASSED") << nl << endl;

    return nFail ? 1 : 0;
}

// ************************************************************************* //
