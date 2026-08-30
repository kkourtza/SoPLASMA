/*---------------------------------------------------------------------------*\
Application
    testAitken

Description
    Unit test for the outer-loop relaxation used by plasmaTransport, on a
    SYNTHETIC fixed-point iteration whose answer is known analytically.

    WHY A UNIT TEST. The relaxation was first tuned directly on a 1.15M-cell
    streamer case, where one run is minutes and every reading is confounded by
    the time-step controller, the chemistry and the Poisson solve. Two
    successive diagnoses were drawn from that setting and BOTH were wrong: an
    omega of 1 that turned out to be a governor that never ran, and an omega at
    its floor that was read as "Aitken cannot see local trouble". Here the
    iteration is linear, the fixed point is exactly zero, and the optimal
    factor is known in closed form, so a wrong implementation cannot hide.

    THE MODEL PROBLEM. Each cell c is an independent linear fixed-point map

        xTilde_c = -g_c x_c            (fixed point x_c = 0)

    so the Picard residual is r_c = xTilde_c - x_c = -(1 + g_c) x_c and the
    relaxed update x <- x + omega r gives x_c <- x_c (1 - omega (1 + g_c)).

    * g_c > 1 is the PERIOD-2 DIVERGENT case that kills the streamer run:
      undamped (omega = 1) the iterate alternates in sign and grows.
    * The optimal factor is omega* = 1/(1 + g_c), which converges in ONE step.
    * Convergence requires 0 < omega < 2/(1 + g_c).

    THE TWO QUESTIONS THIS ANSWERS.
      1. Does the scheme find omega ~ 1/(1+g) on a UNIFORM problem, or does it
         ratchet down to its floor? (Observed on the CFD case: it floored.)
      2. On a LOCALISED problem -- a few violent cells among many quiescent
         ones, which is what a streamer head is -- does a single GLOBAL omega
         still work, or is a per-cell factor needed?

\*---------------------------------------------------------------------------*/

#include "aitkenRelaxation.H"
#include "andersonRelaxation.H"
#include "IOstreams.H"
#include "scalarField.H"
#include <cmath>

using namespace Foam;

// Run the model problem with one GLOBAL omega from aitkenRelaxation.
// Returns the number of correctors to reach tol, or -1 if it diverged.
static label runGlobal
(
    const scalarField& g,
    const label maxIt,
    const scalar tol,
    scalar& omegaFinal,
    scalar& omegaMinSeen,
    const scalar omega0,
    const scalar omegaMin,
    const scalar omegaMax
)
{
    scalarField x(g.size(), 1.0);          // start away from the fixed point
    aitkenRelaxation relax("test", omega0, omegaMin, omegaMax);
    relax.reset();

    omegaFinal = omega0;
    omegaMinSeen = omegaMax;

    for (label k = 0; k < maxIt; ++k)
    {
        // "solve": xTilde = -g x, hence residual r = xTilde - x
        scalarField r(x.size());
        forAll(x, c)
        {
            r[c] = (-g[c]*x[c]) - x[c];
        }

        const scalar w = relax.omega(r);
        omegaFinal = w;
        omegaMinSeen = min(omegaMinSeen, w);

        // relaxed update
        forAll(x, c)
        {
            x[c] += w*r[c];
        }

        scalar nrm = 0;
        forAll(x, c) { nrm = max(nrm, mag(x[c])); }

        if (nrm < tol)   return k + 1;
        if (nrm > 1e12)  return -1;         // diverged
    }
    return -2;                              // did not converge in maxIt
}


// The same problem with a PER-CELL omega, computed pointwise by the same
// Aitken formula. No reductions at all.
static label runPerCell
(
    const scalarField& g,
    const label maxIt,
    const scalar tol,
    scalar& omegaMinSeen,
    const scalar omega0,
    const scalar omegaMin,
    const scalar omegaMax
)
{
    const label n = g.size();
    scalarField x(n, 1.0);
    scalarField rPrev(n, 0.0);
    scalarField w(n, omega0);
    bool haveHistory = false;

    omegaMinSeen = omegaMax;

    for (label k = 0; k < maxIt; ++k)
    {
        scalarField r(n);
        forAll(x, c)
        {
            r[c] = (-g[c]*x[c]) - x[c];
        }

        if (haveHistory)
        {
            forAll(x, c)
            {
                const scalar dr = r[c] - rPrev[c];

                // Guard: where the residual is not meaningfully changing there
                // is nothing to estimate from, and dividing by noise would
                // damp a cell that is behaving perfectly well.
                if (mag(dr) > 1e-300)
                {
                    scalar wc = -w[c]*rPrev[c]*dr/(dr*dr);
                    if (!(wc > omegaMin)) wc = omegaMin;
                    else if (wc > omegaMax) wc = omegaMax;
                    w[c] = wc;
                }
            }
        }
        haveHistory = true;
        rPrev = r;

        forAll(x, c)
        {
            omegaMinSeen = min(omegaMinSeen, w[c]);
            x[c] += w[c]*r[c];
        }

        scalar nrm = 0;
        forAll(x, c) { nrm = max(nrm, mag(x[c])); }

        if (nrm < tol)   return k + 1;
        if (nrm > 1e12)  return -1;
    }
    return -2;
}


static void report
(
    const word& label_,
    const label its,
    const scalar wFinal,
    const scalar wMin,
    const scalar wOptimal
)
{
    Info<< "  " << label_.c_str() << ": ";
    if (its == -1)      Info<< "DIVERGED";
    else if (its == -2) Info<< "not converged";
    else                Info<< its << " correctors";
    Info<< "   omega_final = " << wFinal
        << "  min(omega) = " << wMin
        << "   (optimal ~ " << wOptimal << ")" << endl;
}


int main()
{
    const label maxIt = 200;
    const scalar tol = 1e-10;
    const scalar wMinBound = 0.05, wMaxBound = 1.0;

    Info<< nl << "=== testAitken: synthetic Picard fixed point ===" << nl
        << "    map xTilde = -g x, fixed point 0, optimal omega = 1/(1+g)"
        << nl << endl;

    // ---- CASE 1: uniform, mildly divergent (the streamer's regime) --------
    {
        const scalar g = 1.2;
        scalarField gf(1000, g);
        scalar wF, wM;

        Info<< "CASE 1  uniform g = " << g
            << "  (undamped gain > 1: period-2 divergent)" << endl;

        const label a = runGlobal(gf, maxIt, tol, wF, wM, 1.0, wMinBound, wMaxBound);
        report("global Aitken ", a, wF, wM, 1.0/(1.0 + g));

        const label b = runPerCell(gf, maxIt, tol, wM, 1.0, wMinBound, wMaxBound);
        report("per-cell Aitken", b, 0.0, wM, 1.0/(1.0 + g));
    }

    // ---- CASE 2: LOCALISED -- a streamer head in a quiet domain ----------
    {
        scalarField gf(1000, 0.05);      // quiescent bulk, strongly contracting
        for (label c = 0; c < 10; ++c) { gf[c] = 1.2; }   // 1% violent cells
        scalar wF, wM;

        Info<< nl << "CASE 2  localised: 10 of 1000 cells at g = 1.2,"
            << " the rest at g = 0.05" << nl
            << "        (this is the streamer head in a quiet domain)" << endl;

        const label a = runGlobal(gf, maxIt, tol, wF, wM, 1.0, wMinBound, wMaxBound);
        report("global Aitken ", a, wF, wM, 1.0/2.2);

        const label b = runPerCell(gf, maxIt, tol, wM, 1.0, wMinBound, wMaxBound);
        report("per-cell Aitken", b, 0.0, wM, 1.0/2.2);
    }

    // ---- CASE 3: strongly divergent --------------------------------------
    {
        const scalar g = 4.0;
        scalarField gf(1000, g);
        scalar wF, wM;

        Info<< nl << "CASE 3  uniform g = " << g << "  (strongly divergent)"
            << endl;

        const label a = runGlobal(gf, maxIt, tol, wF, wM, 1.0, wMinBound, wMaxBound);
        report("global Aitken ", a, wF, wM, 1.0/(1.0 + g));

        const label b = runPerCell(gf, maxIt, tol, wM, 1.0, wMinBound, wMaxBound);
        report("per-cell Aitken", b, 0.0, wM, 1.0/(1.0 + g));
    }

    // ---- CASE 4: benign, must NOT be damped ------------------------------
    {
        const scalar g = 0.2;
        scalarField gf(1000, g);
        scalar wF, wM;

        Info<< nl << "CASE 4  uniform g = " << g
            << "  (already contracting: damping here only costs accuracy)"
            << endl;

        const label a = runGlobal(gf, maxIt, tol, wF, wM, 1.0, wMinBound, wMaxBound);
        report("global Aitken ", a, wF, wM, 1.0/(1.0 + g));

        const label b = runPerCell(gf, maxIt, tol, wM, 1.0, wMinBound, wMaxBound);
        report("per-cell Aitken", b, 0.0, wM, 1.0/(1.0 + g));
    }

    // ---- CASE 5: a CLAMP inside the iteration ----------------------------
    //
    // The real solver applies clampNumberDensities() to every species AFTER
    // the corrector's update, so the map Aitken actually sees is
    // x -> clamp(solve(x)), not the smooth linear map it assumes.
    //
    // The fixed point is x* = 1 and the clamp floor is BELOW it, so the clamp
    // bites only while the iterate overshoots downwards during the transient
    // and is inactive at the answer. (A first version of this case put the
    // fixed point AT the floor, so the clamp drove every cell onto it and the
    // test "converged" in one corrector while testing nothing.)
    {
        const scalar g = 1.2;
        const scalar xStar = 1.0;
        const scalar floorVal = 0.5;
        const label n = 1000;

        for (label mode = 0; mode < 2; ++mode)
        {
            const bool clamped = (mode == 1);

            scalarField x(n, 3.0);
            aitkenRelaxation relax("c", 1.0, wMinBound, wMaxBound);
            relax.reset();

            label its = -2; scalar wLast = 1.0, wMin = wMaxBound;

            for (label k = 0; k < maxIt; ++k)
            {
                // xTilde = x* - g (x - x*)   =>   r = (1+g)(x* - x)
                scalarField r(n);
                forAll(x, c) { r[c] = (1.0 + g)*(xStar - x[c]); }

                const scalar w = relax.omega(r);
                wLast = w; wMin = min(wMin, w);

                forAll(x, c)
                {
                    x[c] += w*r[c];
                    if (clamped && x[c] < floorVal) { x[c] = floorVal; }
                }

                scalar nrm = 0;
                forAll(x, c) { nrm = max(nrm, mag(x[c] - xStar)); }
                if (nrm < tol) { its = k + 1; break; }
                if (nrm > 1e12) { its = -1; break; }
            }

            if (!clamped)
            {
                Info<< nl << "CASE 5  fixed point x* = 1, g = " << g << endl;
                report("no clamp      ", its, wLast, wMin, 1.0/(1.0 + g));
            }
            else
            {
                report("WITH clamp    ", its, wLast, wMin, 1.0/(1.0 + g));
                Info<< "        -> a clamp inside the iteration is a"
                    << " NON-SMOOTH projection; Aitken's linear" << nl
                    << "           extrapolation is not valid across it."
                    << endl;
            }
        }
    }

    // ---- CASE 6: MAGNITUDE DISPARITY, as in a real density field ---------
    //
    // n_e spans ~1e5 in the background to ~1e19 in the streamer head. A GLOBAL
    // Aitken factor is built from inner products over the whole field, so it
    // is dominated by whichever cells are numerically largest -- which need
    // not be the cells carrying the unstable mode. A per-cell factor is
    // scale-free by construction: each cell's ratio has its own magnitude in
    // both numerator and denominator.
    {
        const label n = 1000;
        const scalar gOsc = 1.2, gCalm = 0.05;

        scalarField g(n, gCalm);
        scalarField x0(n);
        forAll(x0, c)
        {
            // background 1e5, rising to 1e19 in a small "head"
            x0[c] = (c < 50) ? 1e19 : 1e5;
            // the UNSTABLE cells are NOT the largest ones -- they sit at the
            // edge of the head, which is exactly where a streamer front is
            if (c >= 40 && c < 50) { g[c] = gOsc; }
        }

        Info<< nl << "CASE 6  magnitude spans 1e5..1e19; the 10 unstable cells"
            << " are NOT the largest" << nl
            << "        (a real n_e field: the front, not the peak, is what"
            << " oscillates)" << endl;

        // --- global
        {
            scalarField x(x0);
            aitkenRelaxation relax("g", 1.0, wMinBound, wMaxBound);
            relax.reset();
            label its = -2; scalar wLast = 1.0, wMin = wMaxBound;
            for (label k = 0; k < maxIt; ++k)
            {
                scalarField r(n);
                forAll(x, c) { r[c] = (-g[c]*x[c]) - x[c]; }
                const scalar w = relax.omega(r);
                wLast = w; wMin = min(wMin, w);
                forAll(x, c) { x[c] += w*r[c]; }
                scalar nrm = 0;
                forAll(x, c) { nrm = max(nrm, mag(x[c])/max(mag(x0[c]), SMALL)); }
                if (nrm < tol) { its = k + 1; break; }
                if (nrm > 1e12) { its = -1; break; }
            }
            report("global Aitken ", its, wLast, wMin, 1.0/(1.0 + gOsc));
        }

        // --- per cell
        {
            scalarField x(x0), rPrev(n, 0.0), w(n, 1.0);
            bool hist = false;
            label its = -2; scalar wMin = wMaxBound;
            for (label k = 0; k < maxIt; ++k)
            {
                scalarField r(n);
                forAll(x, c) { r[c] = (-g[c]*x[c]) - x[c]; }
                if (hist)
                {
                    forAll(x, c)
                    {
                        const scalar dr = r[c] - rPrev[c];
                        if (mag(dr) > 1e-300)
                        {
                            scalar wc = -w[c]*rPrev[c]*dr/(dr*dr);
                            if (!(wc > wMinBound)) wc = wMinBound;
                            else if (wc > wMaxBound) wc = wMaxBound;
                            w[c] = wc;
                        }
                    }
                }
                hist = true; rPrev = r;
                forAll(x, c) { wMin = min(wMin, w[c]); x[c] += w[c]*r[c]; }
                scalar nrm = 0;
                forAll(x, c) { nrm = max(nrm, mag(x[c])/max(mag(x0[c]), SMALL)); }
                if (nrm < tol) { its = k + 1; break; }
                if (nrm > 1e12) { its = -1; break; }
            }
            report("per-cell Aitken", its, 0.0, wMin, 1.0/(1.0 + gOsc));
        }
    }

    // ---- CASE 7: TWO COUPLED FIELDS, only one of them relaxed ------------
    //
    // This is what the CFD actually does today: plasmaTransport relaxes the
    // SPECIES, and nothing relaxes nEps_e. The instability is a period-2 cycle
    // BETWEEN the two, so damping one side leaves the loop gain almost
    // untouched while making the map Aitken sees non-stationary -- its
    // estimate then chases a target that moves every corrector.
    //
    // Model: aTilde = -gab*b, bTilde = -gba*a, fixed point (0,0), loop gain
    // gab*gba. Relaxing only a cannot change the product.
    {
        const label n = 200;
        // Composite gain must be NEGATIVE to be a PERIOD-2 cycle, which is
        // what the solver actually shows (n_e and nEps_e alternate in
        // antiphase). A first version of this case used two negative signs,
        // giving composite +1.43 -- a MONOTONE divergence that no
        // under-relaxation can fix, and which would have "proved" the method
        // useless. Composite here is -gab*gba = -1.43.
        const scalar gab = 1.3, gba = 1.1;
        Info<< nl << "CASE 7  two coupled fields, loop gain "
            << gab*gba << nl
            << "        (a <- +" << gab << " b,  b <- -" << gba
            << " a; the n_e / nEps_e cycle)" << endl;

        for (label mode = 0; mode < 2; ++mode)
        {
            const bool relaxBoth = (mode == 1);

            scalarField a(n, 1.0), b(n, 1.0);
            aitkenRelaxation rA("a", 1.0, wMinBound, wMaxBound);
            aitkenRelaxation rB("b", 1.0, wMinBound, wMaxBound);
            rA.reset(); rB.reset();

            label its = -2; scalar wMin = wMaxBound, wLast = 1.0;

            for (label k = 0; k < maxIt; ++k)
            {
                scalarField ra(n), rb(n);
                forAll(a, c) { ra[c] = ( gab*b[c]) - a[c]; }
                forAll(b, c) { rb[c] = (-gba*a[c]) - b[c]; }

                const scalar wa = rA.omega(ra);
                wLast = wa; wMin = min(wMin, wa);
                forAll(a, c) { a[c] += wa*ra[c]; }

                if (relaxBoth)
                {
                    const scalar wb = rB.omega(rb);
                    wMin = min(wMin, wb);
                    forAll(b, c) { b[c] += wb*rb[c]; }
                }
                else
                {
                    // b takes the full undamped update -- today's CFD
                    forAll(b, c) { b[c] += rb[c]; }
                }

                scalar nrm = 0;
                forAll(a, c) { nrm = max(nrm, max(mag(a[c]), mag(b[c]))); }
                if (nrm < tol) { its = k + 1; break; }
                if (nrm > 1e12) { its = -1; break; }
            }

            report
            (
                relaxBoth ? "BOTH relaxed  " : "only a relaxed",
                its, wLast, wMin, 0.0
            );
        }

        // --- JOINT omega: one factor from the CONCATENATED residual --------
        //
        // Aitken for a VECTOR fixed point takes the whole residual vector, not
        // one component at a time. The unstable mode here lives in the
        // a-b coupling, so no per-field factor can see it: the loop gain is
        // the PRODUCT gab*gba and neither field's own residual contains it.
        // Concatenating makes the coupled mode visible to the estimate.
        {
            scalarField a(n, 1.0), b(n, 1.0);
            aitkenRelaxation rJ("joint", 1.0, wMinBound, wMaxBound);
            rJ.reset();
            label its = -2; scalar wMin = wMaxBound, wLast = 1.0;

            for (label k = 0; k < maxIt; ++k)
            {
                scalarField rj(2*n);
                forAll(a, c) { rj[c]     = ( gab*b[c]) - a[c]; }
                forAll(b, c) { rj[n + c] = (-gba*a[c]) - b[c]; }

                const scalar w = rJ.omega(rj);
                wLast = w; wMin = min(wMin, w);

                forAll(a, c) { a[c] += w*rj[c];     }
                forAll(b, c) { b[c] += w*rj[n + c]; }

                scalar nrm = 0;
                forAll(a, c) { nrm = max(nrm, max(mag(a[c]), mag(b[c]))); }
                if (nrm < tol) { its = k + 1; break; }
                if (nrm > 1e12) { its = -1; break; }
            }
            report("JOINT omega   ", its, wLast, wMin, 0.0);
        }
        Info<< "        -> if 'only a relaxed' floors omega or diverges while"
            << " 'BOTH' converges," << nl
            << "           the missing nEps_e relaxation is the defect,"
            << " not the algorithm." << endl;
    }

    // ---- CASE 8: OUTLIER REJECTION (the descent rate limiter) ------------
    //
    // On the CFD case the omega trace read
    //     ... 0.561  0.184  0.05  0.05  0.539 ...
    // i.e. it fell onto its floor and came back. omega = 0.05 implies a loop
    // gain of ~19, while a CONSTANT factor of 0.8 runs that case stably (gain
    // < 0.25) -- so those floor hits were noise, not signal.
    //
    // Here a clean g = 1.2 iteration is interrupted by ONE corrector whose
    // residual is corrupted (a spike, as if the stiff chemistry had just moved
    // the state under the estimator). Without the limiter that single sample
    // drags omega down; with it, the estimate degrades by at most one factor
    // of descentLimit and recovers.
    {
        const scalar g = 1.2;
        const label n = 500;
        // The clean iteration converges in TWO correctors, so the spike must
        // land inside that window or it never fires -- a first version put it
        // at corrector 4 and the test silently exercised nothing.
        const label spikeAt = 1;

        for (label mode = 0; mode < 2; ++mode)
        {
            const scalar dl = (mode == 0) ? 0.0 : 0.5;   // 0 = limiter OFF

            scalarField x(n, 1.0);
            aitkenRelaxation relax("spike", 1.0, wMinBound, wMaxBound, dl);
            relax.reset();

            scalar wWorst = wMaxBound;
            scalarField wSeen(6, 0.0);
            label its = -2;

            for (label k = 0; k < maxIt; ++k)
            {
                scalarField r(n);
                forAll(x, c) { r[c] = (-g*x[c]) - x[c]; }

                // ONE corrupted corrector: the residual is not the linear
                // mode Aitken assumes.
                if (k == spikeAt)
                {
                    forAll(r, c)
                    {
                        r[c] *= ((c % 2) ? -37.0 : 41.0);
                    }
                }

                const scalar w = relax.omega(r);
                wWorst = min(wWorst, w);
                wSeen[min(k, 5)] = w;
                forAll(x, c) { x[c] += w*r[c]; }

                scalar nrm = 0;
                forAll(x, c) { nrm = max(nrm, mag(x[c])); }
                if (nrm < tol) { its = k + 1; break; }
                if (nrm > 1e12) { its = -1; break; }
            }

            Info<< (mode ? "" : "\nCASE 8  one corrupted corrector at index ")
                << (mode ? "" : Foam::name(spikeAt).c_str()) << endl;
            report
            (
                dl > 0 ? "limiter ON    " : "limiter OFF   ",
                its, wSeen[min(its > 0 ? its - 1 : 5, 5)], wWorst,
                1.0/(1.0 + g)
            );
            Info<< "                  omega trace:";
            for (label q = 0; q < min(its > 0 ? its : 6, 6); ++q)
            {
                Info<< " " << wSeen[q];
            }
            Info<< endl;
        }
        Info<< "        -> the limiter must keep min(omega) OFF the floor ("
            << wMinBound << ") without" << nl
            << "           blunting the genuine response seen in CASE 1."
            << endl;
    }

    // ---- CASE 9: DOES A VARYING omega DESTABILISE THE ITERATION? ---------
    //
    // The CFD evidence is ambiguous. A run with omega varying per corrector
    // (0.05 -> 0.539 -> 0.284 -> 0.05) crashed at a deltaT where a CONSTANT
    // 0.8 survived, even though the varying scheme was damping HARDER on
    // average. One explanation is that a per-corrector omega makes the
    // iteration matrix non-stationary and that this is itself destabilising.
    // That was a HYPOTHESIS; this case measures it.
    //
    // Three policies, same coupled two-field problem, same corrector budget:
    //   (a) CONSTANT omega for the whole run;
    //   (b) omega from joint Aitken, updated EVERY corrector (what we ship);
    //   (c) omega from joint Aitken but FROZEN for the duration of each step,
    //       updated only between steps -- stationary within a step.
    //
    // The gain CHANGES between steps, as it does in the real solver as the
    // streamer develops. A step "fails" if it does not converge in nCorrMax.
    {
        const label n = 100;
        const label nSteps = 40;
        const label nCorrMax = 20;
        // METRIC: does the step CONTRACT, not "did it reach 1e-8". For this
        // map the convergence factor is sqrt((1-w)^2 + w^2 g), which at
        // w = 0.45, g = 2 is 0.84 -- stable, but needing ~106 correctors to
        // reach 1e-8. A first version of this case asked for 1e-8 in 20 and
        // recorded 40/40 failures for ALL THREE policies, measuring the budget
        // rather than the stability it was supposed to discriminate.
        const scalar contractTarget = 1.0;   // final < initial => stable

        // gains sweeping through the marginal region, as deltaT grows
        auto gainAt = [&](const label st) -> scalar
        {
            return 0.6 + 1.4*scalar(st)/scalar(nSteps - 1);   // 0.6 .. 2.0
        };

        Info<< nl << "CASE 9  " << nSteps << " steps, loop gain sweeping"
            << " 0.6 -> 2.0, budget " << nCorrMax << " correctors/step" << nl
            << "        (a) constant  (b) omega per corrector"
            << "  (c) omega frozen within each step" << endl;

        for (label policy = 0; policy < 3; ++policy)
        {
            scalarField a(n, 1.0), b(n, 1.0);
            aitkenRelaxation rJ("j", 1.0, wMinBound, wMaxBound, 0.0);
            label failed = 0, totalCorr = 0;
            scalar wFrozen = 0.5;

            for (label st = 0; st < nSteps; ++st)
            {
                const scalar g = gainAt(st);
                const scalar gab = g, gba = 1.0;   // composite -g: period 2

                // a fresh fixed point each step, as in a real time step
                const scalar aStar = 1.0 + 0.1*st;
                forAll(a, c) { a[c] = aStar*1.5; b[c] = aStar*1.5; }

                rJ.reset();
                scalar wLastStep = wFrozen;
                scalar nrm0 = -1, nrm = 0;

                for (label k = 0; k < nCorrMax; ++k)
                {
                    scalarField rj(2*n);
                    forAll(a, c)
                    {
                        rj[c] = (aStar + gab*(b[c] - aStar)) - a[c];
                    }
                    forAll(b, c)
                    {
                        rj[n + c] = (aStar - gba*(a[c] - aStar)) - b[c];
                    }

                    scalar w = 0;
                    if (policy == 0)
                    {
                        w = 0.45;                        // constant
                    }
                    else if (policy == 1)
                    {
                        w = rJ.omega(rj);                // every corrector
                    }
                    else
                    {
                        // frozen: still measure, but apply last step's value
                        wLastStep = rJ.omega(rj);
                        w = wFrozen;
                    }

                    forAll(a, c) { a[c] += w*rj[c];     }
                    forAll(b, c) { b[c] += w*rj[n + c]; }

                    ++totalCorr;

                    nrm = 0;
                    forAll(a, c)
                    {
                        nrm = max(nrm, max(mag(a[c]-aStar), mag(b[c]-aStar)));
                    }
                    if (nrm0 < 0) { nrm0 = nrm; }
                    if (nrm > 1e12) break;
                }

                if (policy == 2) { wFrozen = wLastStep; }

                // DIVERGED if the step ended no closer than it started.
                if (!(nrm < contractTarget*nrm0)) { ++failed; }
            }

            const char* nm =
                (policy == 0) ? "constant 0.45 " :
                (policy == 1) ? "per corrector " : "frozen per step";

            Info<< "  " << nm << ": " << failed << " step(s) DIVERGED of "
                << nSteps << ",  " << totalCorr << " correctors total" << endl;
        }
        Info<< "        -> if (b) fails where (a) and (c) do not, a"
            << " per-corrector omega is" << nl
            << "           itself destabilising and should be frozen"
            << " within a step." << endl;
    }

    // ---- CASE 10: ANDERSON vs JOINT AITKEN on the CASE 7 cycle -----------
    //
    // CASE 7 is the regime that matters: the two-field period-2 cycle, which
    // is what n_e / nEps_e do above the critical deltaT. The best SCALAR
    // policy there is joint Aitken at 91 correctors -- damping buys stability
    // by slowing every mode, including the ones that were already fine.
    //
    // Anderson builds the next iterate from a least-squares combination of the
    // last m residuals, so it can CANCEL the offending mode instead. This
    //
    // HOW MANY VECTORS DOES IT NEED? The coupled map has Jacobian
    //     [  0    gab ]
    //     [ -gba   0  ]
    // whose eigenvalues are +/- i sqrt(gab gba) = +/- 1.196i -- a COMPLEX
    // CONJUGATE PAIR, not a single mode. A real Krylov method cannot represent
    // half of a conjugate pair, so m = 1 CANNOT converge here; m = 2 is the
    // true minimum, terminating in about m+1 correctors as GMRES does on a
    // 2-dimensional Krylov space.
    //
    // Stated precisely because the first version of this comment asserted "one
    // active mode, so m = 1 should terminate", and the very first run
    // contradicted it. ONE OSCILLATORY MODE IS TWO EIGENVALUES.
    //
    // Reported alongside: how often the safeguards fired. An acceleration that
    // only converges because it silently fell back to plain relaxation has not
    // accelerated anything, and the counters are the only way to tell.
    {
        const label n = 200;
        const scalar gab = 1.3, gba = 1.1;

        Info<< nl << "CASE 10  Anderson vs joint Aitken on the CASE 7 cycle"
            << nl << "        (a <- +" << gab << " b,  b <- -" << gba
            << " a; loop gain " << gab*gba
            << ", eigenvalues +/- 1.196i)" << endl;

        // --- reference: joint Aitken, as CASE 7 measures it
        {
            scalarField a(n, 1.0), b(n, 1.0);
            aitkenRelaxation rJ("joint", 1.0, wMinBound, wMaxBound);
            rJ.reset();
            label its = -2;

            for (label k = 0; k < maxIt; ++k)
            {
                scalarField rj(2*n);
                forAll(a, c) { rj[c]     = ( gab*b[c]) - a[c]; }
                forAll(b, c) { rj[n + c] = (-gba*a[c]) - b[c]; }

                const scalar w = rJ.omega(rj);
                forAll(a, c) { a[c] += w*rj[c];     }
                forAll(b, c) { b[c] += w*rj[n + c]; }

                scalar nrm = 0;
                forAll(a, c) { nrm = max(nrm, max(mag(a[c]), mag(b[c]))); }
                if (nrm < tol) { its = k + 1; break; }
                if (nrm > 1e12) { its = -1; break; }
            }
            Info<< "  joint Aitken   : ";
            if (its == -1)      Info<< "DIVERGED";
            else if (its == -2) Info<< "not converged in " << maxIt;
            else                Info<< its << " correctors";
            Info<< endl;
        }

        // --- Anderson, m = 1 .. 5
        for (label m = 1; m <= 5; ++m)
        {
            scalarField x(2*n, 1.0);
            andersonRelaxation aa("joint", m, 1.0);
            aa.reset();
            label its = -2;

            for (label k = 0; k < maxIt; ++k)
            {
                scalarField r(2*n);
                forAll(r, c)
                {
                    r[c] = (c < n)
                        ? ( gab*x[n + c]) - x[c]
                        : (-gba*x[c - n]) - x[c];
                }

                aa.correct(x, r);

                scalar nrm = 0;
                forAll(x, c) { nrm = max(nrm, mag(x[c])); }
                if (nrm < tol) { its = k + 1; break; }
                if (nrm > 1e12) { its = -1; break; }
            }

            Info<< "  Anderson m = " << m << "  : ";
            if (its == -1)      Info<< "DIVERGED";
            else if (its == -2) Info<< "not converged in " << maxIt;
            else                Info<< its << " correctors";
            Info<< "   fallbacks: step " << aa.nFallbackStep()
                << ", solve " << aa.nFallbackSolve()
                << "   max depth used " << aa.mUsedMax() << endl;
        }
        Info<< "        -> the oscillatory mode is a COMPLEX CONJUGATE PAIR,"
            << " so m = 1 cannot" << nl
            << "           represent it; m = 2 is the minimum and ~m+1"
            << " correctors is GMRES" << nl
            << "           behaviour. Fallbacks above zero would mean the"
            << " safeguards, not the" << nl
            << "           method, are doing the converging." << endl;
    }

    // ---- CASE 11: does Anderson survive a NON-SMOOTH iteration? ----------
    //
    // The honest counter-test. Anderson's extrapolation assumes the map is
    // (locally) linear, and this solver's is not: densities are clamped to a
    // floor, which is a projection, and CASE 5 already records that Aitken's
    // linear extrapolation is invalid across it. Anderson is MORE aggressive
    // than Aitken, so if the clamp is going to break something it breaks this.
    //
    // A method that wins CASE 10 and blows up here is not ready for the CFD,
    // and the step limiter is exactly what should catch it. This case measures
    // whether it does.
    {
        const label n = 200;
        const scalar gab = 1.3, gba = 1.1;
        const scalar floorVal = 0.35;    // the clamp, active near the answer

        Info<< nl << "CASE 11  the same cycle with a CLAMP inside the iteration"
            << nl << "        (x floored at " << floorVal
            << "; fixed point 1, so the clamp is active on the way in)" << endl;

        for (label mode = 0; mode < 6; ++mode)
        {
            const label m = mode;        // 0 = plain relaxation, else Anderson
            scalarField x(2*n, 3.0);
            andersonRelaxation aa("clamped", m, 1.0);
            aa.reset();
            label its = -2;

            for (label k = 0; k < maxIt; ++k)
            {
                // fixed point at x = 1 in both blocks
                scalarField r(2*n);
                forAll(r, c)
                {
                    r[c] = (c < n)
                        ? (1.0 + gab*(x[n + c] - 1.0)) - x[c]
                        : (1.0 - gba*(x[c - n] - 1.0)) - x[c];
                }

                aa.correct(x, r);

                // THE NON-SMOOTHNESS: a hard floor, as the solver applies to
                // every density.
                forAll(x, c) { x[c] = max(x[c], floorVal); }

                scalar nrm = 0;
                forAll(x, c) { nrm = max(nrm, mag(x[c] - 1.0)); }
                if (nrm < tol) { its = k + 1; break; }
                if (nrm > 1e12) { its = -1; break; }
            }

            if (m == 0) Info<< "  plain (m = 0)  : ";
            else        Info<< "  Anderson m = " << m << "  : ";
            if (its == -1)      Info<< "DIVERGED";
            else if (its == -2) Info<< "not converged in " << maxIt;
            else                Info<< its << " correctors";
            Info<< "   fallbacks: step " << aa.nFallbackStep()
                << ", solve " << aa.nFallbackSolve() << endl;
        }
        Info<< "        -> DIVERGED here would disqualify Anderson for the"
            << " solver regardless of" << nl
            << "           CASE 10. A rising `step` fallback count is the"
            << " limiter working." << endl;

        // --- THE HARDER VARIANT: the clamp is ACTIVE AT THE ANSWER ---------
        //
        // Above, the fixed point is 1 and the floor is 0.35, so the clamp is
        // crossed on the way in and INACTIVE once inside the basin -- which
        // flatters the method. The solver's situation is the opposite: most of
        // a streamer domain sits ON the density floor at convergence, so the
        // projection is active AT the fixed point for those cells, and the map
        // is non-differentiable exactly where the iteration has to stop.
        //
        // Half the cells are given a fixed point BELOW the floor, so their
        // converged value is the floor itself and the active set is part of the
        // answer.
        Info<< nl << "         clamp ACTIVE AT the fixed point"
            << " (half the cells converge ON the floor)" << endl;

        for (label mode = 0; mode < 6; ++mode)
        {
            const label m = mode;
            scalarField x(2*n, 3.0);

            // per-cell target: below the floor for half of them
            scalarField xStar(2*n, 1.0);
            forAll(xStar, c) { if ((c % 2) == 0) xStar[c] = 0.10; }

            andersonRelaxation aa("clampedActive", m, 1.0);
            aa.reset();
            label its = -2;

            for (label k = 0; k < maxIt; ++k)
            {
                scalarField r(2*n);
                forAll(r, c)
                {
                    r[c] = (c < n)
                        ? (xStar[c] + gab*(x[n + c] - xStar[n + c])) - x[c]
                        : (xStar[c] - gba*(x[c - n] - xStar[c - n])) - x[c];
                }

                // MEASURE THE STEP, NOT THE DISTANCE TO A GUESSED ANSWER.
                //
                // The first version of this case scored convergence against
                // max(xStar, floor). That is NOT a fixed point of the clamped
                // map -- P(G(x)) = x is a complementarity problem whose
                // solution is not the projection of the unclamped answer -- so
                // it recorded "not converged" for EVERY method including plain
                // relaxation, which is the tell that the target was wrong
                // rather than the methods.
                //
                // The stagnation of the iterate needs no analytic answer and is
                // exactly what the solver's own outer loop tests.
                const scalarField xOld(x);

                aa.correct(x, r);
                forAll(x, c) { x[c] = max(x[c], floorVal); }

                scalar nrm = 0;
                forAll(x, c) { nrm = max(nrm, mag(x[c] - xOld[c])); }
                if (nrm < tol) { its = k + 1; break; }
                if (nrm > 1e12) { its = -1; break; }
            }

            if (m == 0) Info<< "  plain (m = 0)  : ";
            else        Info<< "  Anderson m = " << m << "  : ";
            if (its == -1)      Info<< "DIVERGED";
            else if (its == -2) Info<< "not converged in " << maxIt;
            else                Info<< its << " correctors";
            Info<< "   fallbacks: step " << aa.nFallbackStep()
                << ", solve " << aa.nFallbackSolve() << endl;
        }
        Info<< "        -> this is the realistic one: the active set is part of"
            << " the answer, so the" << nl
            << "           map is non-differentiable AT the fixed point."
            << " Anderson failing HERE" << nl
            << "           but passing above would mean the earlier result was"
            << " an artefact of a" << nl
            << "           clamp that switched itself off." << endl;
    }

    // ---- CASE 12: Anderson on a LOCALISED / disparate-magnitude field ----
    //
    // CASES 10 and 11 are spatially uniform: every cell has the same gain, so a
    // single global correction fits every cell exactly. A streamer is not like
    // that. CASE 2 and CASE 6 exist because it matters -- global Aitken needs
    // 8-9 correctors on a localised problem against 2 on a uniform one.
    //
    // Anderson's correction is also global (one gamma for the whole vector),
    // so the same objection applies to it and has to be tested rather than
    // assumed. Two beds, both from the cases above:
    //   (a) 10 unstable cells in 1000 quiet ones;
    //   (b) the same, but with magnitudes spanning 1e5 .. 1e19, the unstable
    //       cells NOT being the largest -- a real n_e field, where the front
    //       oscillates and the peak does not.
    // (b) is the one that punishes an unscaled inner product: a residual sum
    // dominated by 1e19 cells cannot see a 1e5 cell misbehaving.
    {
        const label n = 1000;
        const label nBad = 10;

        for (label bed = 0; bed < 2; ++bed)
        {
            scalarField g(n, 0.05);
            for (label c = 0; c < nBad; ++c) { g[n/2 + c] = 1.2; }

            // (b) gives each cell its own magnitude; the unstable ones sit in
            // the middle of the range, not at the top.
            scalarField scaleF(n, 1.0);
            if (bed == 1)
            {
                forAll(scaleF, c)
                {
                    scaleF[c] = Foam::pow(10.0, 5.0 + 14.0*scalar(c)/scalar(n-1));
                }
                for (label c = 0; c < nBad; ++c) { scaleF[n/2 + c] = 1e9; }
            }

            Info<< nl << "CASE 12" << (bed ? "b  localised AND magnitudes 1e5..1e19"
                                            : "a  localised: 10 of 1000 cells at g = 1.2")
                << endl;

            // reference: global Aitken
            {
                scalar wF = 1, wMin = 1;
                const label its = runGlobal
                (
                    g, maxIt, tol, wF, wMin, 1.0, wMinBound, wMaxBound
                );
                report("global Aitken ", its, wF, wMin, 1.0/(1.0 + 1.2));
            }

            for (label m = 2; m <= 4; ++m)
            {
                scalarField x(n);
                forAll(x, c) { x[c] = scaleF[c]; }

                andersonRelaxation aa("loc", m, 1.0);
                aa.reset();
                label its = -2;

                for (label k = 0; k < maxIt; ++k)
                {
                    scalarField r(n);
                    forAll(r, c) { r[c] = (-g[c]*x[c]) - x[c]; }

                    aa.correct(x, r);

                    // RELATIVE convergence, per cell, against that cell's own
                    // scale: an absolute norm on a field spanning 1e14 in
                    // magnitude measures only the largest cells.
                    scalar nrm = 0;
                    forAll(x, c) { nrm = max(nrm, mag(x[c])/scaleF[c]); }
                    if (nrm < tol) { its = k + 1; break; }
                    if (nrm > 1e12) { its = -1; break; }
                }

                Info<< "  Anderson m = " << m << "  : ";
                if (its == -1)      Info<< "DIVERGED";
                else if (its == -2) Info<< "not converged in " << maxIt;
                else                Info<< its << " correctors";
                Info<< "   fallbacks: step " << aa.nFallbackStep()
                    << ", solve " << aa.nFallbackSolve() << endl;
            }
        }
        Info<< "        -> a global correction cannot fit cells with DIFFERENT"
            << " gains exactly, so" << nl
            << "           some iteration is expected here. Failing (b) while"
            << " passing (a) would" << nl
            << "           mean the inner products are dominated by magnitude"
            << " and the method is" << nl
            << "           blind to where the trouble actually is." << endl;
    }

    // ---- CASE 13: HOW FAR CAN THE GAIN GO? -------------------------------
    //
    // THE QUESTION THIS ANSWERS: if Anderson replaces Aitken, how much larger
    // can deltaT be? The Picard loop gain grows with deltaT, so the largest
    // gain each method can still converge at, within a realistic corrector
    // budget, IS the deltaT headroom in the coupling's own units.
    //
    // Measured on the two-field cycle, budget 20 correctors per step, which is
    // the shipped maxCorrectors. A method that "converges" only by taking 200
    // correctors is of no use to a solver that allows 20.
    {
        const label n = 100;
        const label budget = 20;

        Info<< nl << "CASE 13  largest loop gain still converged within "
            << budget << " correctors" << nl
            << "        (the coupling gain grows with deltaT, so this IS the"
            << " deltaT headroom)" << endl;

        // PARAMETRISE BY THE LOOP GAIN, KEEPING THE SHAPE OF CASE 10.
        //
        // A first version set gab = gain, gba = 1, which is a DEGENERATE
        // structure: it reported Anderson m = 1 converging at every gain up to
        // 200, flatly contradicting CASE 10 where m = 1 fails at gain 1.43, and
        // joint Aitken failing at gain 0.1 where CASE 7 has it converging at
        // 1.43. Both contradictions came from the bed, not the methods.
        //
        // Here gab*gba = G and gab/gba = 1.3/1.1, so the map is CASE 10's with
        // only its gain changed -- the one variable this sweep is about.
        auto converges = [&](const scalar G, const label m) -> bool
        {
            const scalar ratio = 1.3/1.1;
            const scalar gab = Foam::sqrt(G*ratio);
            const scalar gba = Foam::sqrt(G/ratio);
            scalarField x(2*n, 1.0);
            andersonRelaxation aa("sweep", m, 1.0);
            aitkenRelaxation rJ("sweep", 1.0, wMinBound, wMaxBound);
            aa.reset(); rJ.reset();

            for (label k = 0; k < budget; ++k)
            {
                scalarField r(2*n);
                forAll(r, c)
                {
                    r[c] = (c < n)
                        ? ( gab*x[n + c]) - x[c]
                        : (-gba*x[c - n]) - x[c];
                }

                if (m == 0)
                {
                    const scalar w = rJ.omega(r);
                    forAll(x, c) { x[c] += w*r[c]; }
                }
                else
                {
                    aa.correct(x, r);
                }

                scalar nrm = 0;
                forAll(x, c) { nrm = max(nrm, mag(x[c])); }
                if (nrm < tol) return true;
                if (nrm > 1e12) return false;
            }
            return false;
        };

        // The bisection assumes convergence is monotone in G: true at the
        // bottom, false at the top. Both ends are CHECKED, because a bisection
        // on a predicate that is false at the bottom or true at the top
        // returns the bracket end and looks like a result.
        // Wide enough to bracket BOTH methods: Aitken is slow enough that its
        // limit falls far below 1 within a 20-corrector budget, and Anderson's
        // is far above it. A bracket that fails to contain the answer reports
        // "nothing measured" rather than silently returning the bracket end.
        const scalar Glo = 1e-3, Ghi = 1e6;

        for (label m = 0; m <= 4; ++m)
        {
            const char* nm = (m == 0) ? "joint Aitken  " : "Anderson      ";

            Info<< "  " << nm;
            if (m) Info<< "m = " << m << "  : "; else Info<< "  : ";

            if (!converges(Glo, m))
            {
                Info<< "does not converge even at gain " << Glo
                    << " -- no bracket, nothing measured" << endl;
                continue;
            }
            if (converges(Ghi, m))
            {
                Info<< "still converges at gain " << Ghi
                    << " -- limit is above the bracket, nothing measured"
                    << endl;
                continue;
            }

            scalar lo = Glo, hi = Ghi;
            for (label it = 0; it < 40; ++it)
            {
                const scalar mid = 0.5*(lo + hi);
                if (converges(mid, m)) lo = mid; else hi = mid;
            }
            Info<< "largest loop gain " << lo << endl;
        }
        Info<< "        -> READ THIS WITH ITS CAVEAT. The map here is LINEAR,"
            << " and Anderson is" << nl
            << "           essentially exact on linear problems, so a limit"
            << " above the bracket is" << nl
            << "           EXPECTED and does NOT extrapolate to the solver."
            << " CASE 11's clamped" << nl
            << "           variant is the realistic bound: there the advantage"
            << " is real but small" << nl
            << "           (8 correctors, and non-monotone in m)." << nl
            << "           For deltaT: the coupling governed 23% of steps on"
            << " the streamer runs," << nl
            << "           so that is the addressable share; the transport"
            << " Courant number then" << nl
            << "           becomes binding, and its median headroom there was"
            << " about 4x." << endl;
    }

    // ---- CASE 14: CAN rho REPLACE omega AS THE deltaT GOVERNOR? ----------
    //
    // THE PROBLEM. plasmaTimeControl governs deltaT on `relaxMargin_`, the
    // Aitken min(omega). Anderson has no omega, so that governor would
    // silently stop governing the moment Anderson were enabled -- and the log
    // would still print a margin, which is worse than printing nothing.
    //
    // rho = max ||r_k||/||r_{k-1}|| over the step is defined for ANY scheme.
    // But "defined" is not "carries the signal", and that has to be measured.
    //
    // WHY HERE AND NOT ON THE STREAMER. On the CFD the only available check is
    // whether rho agrees with OMEGA -- one opinion against another. Here the
    // loop gain is set directly, so GROUND TRUTH is known: the step either
    // reached tol within the budget or it did not. That makes this the
    // stronger test as well as the one that runs in a second rather than an
    // hour.
    //
    // Loop gain stands in for deltaT: it is what grows as deltaT grows.
    {
        const label n = 100;
        const label budget = 20;
        const scalar ratio = 1.3/1.1;

        Info<< nl << "CASE 14  can rho govern deltaT where omega does?" << nl
            << "        loop gain stands in for deltaT; budget " << budget
            << " correctors/step" << endl;
        Info<< "        gain      omega_min     rho_max   converged?" << endl;

        DynamicList<scalar> omegas, rhos;
        DynamicList<bool> converged;

        for (label gi = 0; gi < 12; ++gi)
        {
            const scalar G = 0.02*Foam::pow(1.6, scalar(gi));
            const scalar gab = Foam::sqrt(G*ratio);
            const scalar gba = Foam::sqrt(G/ratio);

            scalarField a(n, 1.0), b(n, 1.0);
            aitkenRelaxation rJ("gov", 1.0, wMinBound, wMaxBound);
            rJ.reset();

            scalar wMin = 1.0, rhoMax = 0.0, rPrev = -1.0;
            bool conv = false;

            for (label k = 0; k < budget; ++k)
            {
                scalarField rj(2*n);
                forAll(a, c) { rj[c]     = ( gab*b[c]) - a[c]; }
                forAll(b, c) { rj[n + c] = (-gba*a[c]) - b[c]; }

                // rho FIRST, exactly as plasmaOuterRelaxation::applyJoint does
                scalar rN = 0;
                forAll(rj, c) { rN += rj[c]*rj[c]; }
                rN = Foam::sqrt(rN);
                if (rPrev > VSMALL) { rhoMax = max(rhoMax, rN/rPrev); }
                rPrev = rN;

                const scalar w = rJ.omega(rj);
                wMin = min(wMin, w);

                forAll(a, c) { a[c] += w*rj[c];     }
                forAll(b, c) { b[c] += w*rj[n + c]; }

                scalar nrm = 0;
                forAll(a, c) { nrm = max(nrm, max(mag(a[c]), mag(b[c]))); }
                if (nrm < tol) { conv = true; break; }
                if (nrm > 1e12) break;
            }

            omegas.append(wMin);
            rhos.append(rhoMax);
            converged.append(conv);

            Info<< "        " << G << "\t" << wMin << "\t" << rhoMax
                << "\t" << (conv ? "yes" : "NO") << endl;
        }

        // SEPARATION, not correlation. A governor needs a THRESHOLD that tells
        // the two populations apart; a correlation that does not separate them
        // is useless. Best achievable accuracy is reported for each, together
        // with the majority-class baseline -- the score a constant "always
        // says no" classifier would get, which is how a useless threshold can
        // still look good.
        const label N = omegas.size();
        label nConv = 0;
        forAll(converged, i) { if (converged[i]) ++nConv; }
        const scalar baseline =
            max(scalar(nConv), scalar(N - nConv))/scalar(N);

        auto bestAcc = [&](const DynamicList<scalar>& v, const bool below)
        {
            scalar best = 0;
            forAll(v, t)
            {
                label hit = 0;
                forAll(v, i)
                {
                    // `below`: predict converged when v < threshold (rho);
                    // otherwise predict converged when v >= threshold (omega).
                    const bool pred = below ? (v[i] < v[t]) : (v[i] >= v[t]);
                    if (pred == converged[i]) ++hit;
                }
                best = max(best, scalar(hit)/scalar(N));
            }
            return best;
        };

        Info<< nl << "        separating CONVERGED from NOT, over " << N
            << " gains:" << nl
            << "          majority-class baseline : " << baseline << nl
            << "          best omega threshold    : " << bestAcc(omegas, false)
            << nl
            << "          best rho threshold      : " << bestAcc(rhos, true)
            << endl;
        Info<< "        -> rho can govern only if it separates AT LEAST as well"
            << " as omega AND" << nl
            << "           beats the baseline. Equal-to-baseline means the"
            << " threshold is doing" << nl
            << "           nothing and the apparent skill is the class"
            << " imbalance." << endl;
    }

    // ---- CASE 15: the rho THRESHOLD on the beds that look like the solver -
    //
    // CASE 14 established that rho CAN carry the governor signal and put the
    // threshold at ~0.38 -- on a LINEAR, spatially uniform map. Neither
    // property holds in the solver: densities are clamped (a projection, so
    // the map is not differentiable at the active set) and the unstable cells
    // are a few in a quiet domain, spanning many decades in magnitude.
    //
    // A threshold fitted on the linear bed and shipped would be a number
    // derived from a problem the solver does not solve. This measures it on
    // the two beds that resemble it, and takes the MOST CONSERVATIVE.
    {
        const label budget = 20;
        const scalar ratio = 1.3/1.1;

        Info<< nl << "CASE 15  the rho threshold on non-smooth and localised"
            << " beds" << nl
            << "        (CASE 14's 0.38 came from a linear uniform map; the"
            << " solver is neither)" << endl;

        // bed 2 repeats bed 0 with rho computed over the cells NOT sitting on
        // the clamp -- the test of WHY bed 0 fails, and the candidate fix.
        for (label bed = 0; bed < 3; ++bed)
        {
            const label n = 100;
            DynamicList<scalar> omegas, rhos;
            DynamicList<bool> converged;

            const char* bedName =
                (bed == 0) ? "CLAMPED (projection active at the fixed point)"
              : (bed == 1) ? "LOCALISED (10 of 1000 cells unstable, 1e5..1e19)"
              :              "CLAMPED, rho over the UNCLAMPED cells only";

            Info<< nl << "        bed " << bedName << nl
                << "        stiffness  omega_min    rho_max   converged?"
                << endl;

            for (label gi = 0; gi < 12; ++gi)
            {
                const scalar G = 0.02*Foam::pow(1.6, scalar(gi));

                scalar wMin = 1.0, rhoMax = 0.0, rPrev = -1.0;
                bool conv = false;

                if (bed == 0 || bed == 2)
                {
                    // Clamped two-field cycle, floor active AT the answer for
                    // half the cells -- the realistic case, since most of a
                    // streamer domain sits on the density floor.
                    const scalar gab = Foam::sqrt(G*ratio);
                    const scalar gba = Foam::sqrt(G/ratio);
                    const scalar floorVal = 0.35;

                    scalarField x(2*n, 3.0);
                    scalarField xStar(2*n, 1.0);
                    forAll(xStar, c) { if ((c % 2) == 0) xStar[c] = 0.10; }

                    aitkenRelaxation rJ("g", 1.0, wMinBound, wMaxBound);
                    rJ.reset();

                    for (label k = 0; k < budget; ++k)
                    {
                        scalarField r(2*n);
                        forAll(r, c)
                        {
                            r[c] = (c < n)
                                ? (xStar[c] + gab*(x[n+c] - xStar[n+c])) - x[c]
                                : (xStar[c] - gba*(x[c-n] - xStar[c-n])) - x[c];
                        }

                        // THE HYPOTHESIS UNDER TEST. A cell pinned on the
                        // floor has a residual that does not change from one
                        // corrector to the next: x is frozen, so r = xStar - x
                        // is constant. Summed into the norm, those frozen
                        // terms swamp the cells that are still moving and pin
                        // the ratio at exactly 1. bed 2 excludes them.
                        scalar rN = 0;
                        forAll(r, c)
                        {
                            if (bed == 2 && x[c] <= floorVal*(1 + SMALL))
                            {
                                continue;           // on the clamp
                            }
                            rN += r[c]*r[c];
                        }
                        rN = Foam::sqrt(rN);
                        if (rPrev > VSMALL) { rhoMax = max(rhoMax, rN/rPrev); }
                        rPrev = rN;

                        const scalar w = rJ.omega(r);
                        wMin = min(wMin, w);

                        const scalarField xOld(x);
                        forAll(x, c) { x[c] += w*r[c]; }
                        forAll(x, c) { x[c] = max(x[c], floorVal); }

                        scalar nrm = 0;
                        forAll(x, c) { nrm = max(nrm, mag(x[c] - xOld[c])); }
                        if (nrm < tol) { conv = true; break; }
                        if (nrm > 1e12) break;
                    }
                }
                else
                {
                    // Localised single field: 10 unstable cells of 1000, the
                    // unstable ones NOT the largest in magnitude.
                    const label nc = 1000, nBad = 10;
                    scalarField g(nc, 0.05);
                    for (label c = 0; c < nBad; ++c) { g[nc/2 + c] = G; }

                    scalarField scaleF(nc, 1.0);
                    forAll(scaleF, c)
                    {
                        scaleF[c] =
                            Foam::pow(10.0, 5.0 + 14.0*scalar(c)/scalar(nc-1));
                    }
                    for (label c = 0; c < nBad; ++c) { scaleF[nc/2 + c] = 1e9; }

                    scalarField x(scaleF);
                    aitkenRelaxation rJ("g", 1.0, wMinBound, wMaxBound);
                    rJ.reset();

                    for (label k = 0; k < budget; ++k)
                    {
                        scalarField r(nc);
                        forAll(r, c) { r[c] = (-g[c]*x[c]) - x[c]; }

                        scalar rN = 0;
                        forAll(r, c) { rN += r[c]*r[c]; }
                        rN = Foam::sqrt(rN);
                        if (rPrev > VSMALL) { rhoMax = max(rhoMax, rN/rPrev); }
                        rPrev = rN;

                        const scalar w = rJ.omega(r);
                        wMin = min(wMin, w);
                        forAll(x, c) { x[c] += w*r[c]; }

                        scalar nrm = 0;
                        forAll(x, c) { nrm = max(nrm, mag(x[c])/scaleF[c]); }
                        if (nrm < tol) { conv = true; break; }
                        if (nrm > 1e12) break;
                    }
                }

                omegas.append(wMin);
                rhos.append(rhoMax);
                converged.append(conv);

                Info<< "        " << G << "\t" << wMin << "\t" << rhoMax
                    << "\t" << (conv ? "yes" : "NO") << endl;
            }

            // The LARGEST rho among converged steps and the SMALLEST among
            // non-converged. If the first is below the second the two
            // populations are separable and any threshold between them works;
            // if they overlap, rho cannot govern this bed and that is the
            // finding.
            scalar maxConv = -1, minFail = GREAT;
            label nConv = 0;
            forAll(converged, i)
            {
                if (converged[i]) { maxConv = max(maxConv, rhos[i]); ++nConv; }
                else              { minFail = min(minFail, rhos[i]); }
            }

            Info<< "          converged: " << nConv << " of " << rhos.size()
                << endl;
            if (nConv == 0 || nConv == rhos.size())
            {
                Info<< "          no contrast on this bed -- nothing measured"
                    << endl;
            }
            else if (maxConv < minFail)
            {
                Info<< "          SEPARABLE: rho <= " << maxConv
                    << " converged, rho >= " << minFail << " did not" << nl
                    << "          -> any threshold in (" << maxConv << ", "
                    << minFail << "); conservative choice " << maxConv << endl;
            }
            else
            {
                Info<< "          OVERLAP: converged up to rho " << maxConv
                    << " but failures start at " << minFail << nl
                    << "          -> rho alone cannot govern this bed" << endl;
            }
        }
        Info<< nl << "        -> ship the MOST CONSERVATIVE threshold across"
            << " beds, not CASE 14's," << nl
            << "           and if any bed OVERLAPS, rho needs a companion"
            << " signal." << endl;
    }

    // ---- CASE 16: IS `correctors used` A ROBUST GOVERNOR SIGNAL? ---------
    //
    // omega and rho are both PROXIES for "is this step about to fail to
    // converge". The corrector count IS that quantity, is already tracked by
    // plasmaTimeControl::noteOuterLoop(), and is method-agnostic by
    // construction -- which is what omega is not (Anderson has none) and rho is
    // not (CASE 15: blind on a clamped field).
    //
    // But it is only usable if it degrades SMOOTHLY. Today the controller reacts
    // at used == cap, a cliff: deltaT then oscillates across it, paying a
    // discarded step each time. A pre-emptive governor needs LEAD -- the count
    // must rise before the failure, not with it.
    //
    // FOUR WAYS IT COULD FAIL, each measured here rather than assumed:
    //   1. no lead: the count sits low then jumps to the cap;
    //   2. tolerance dependence: a user tightening `tolerance` raises the count
    //      everywhere, so a fixed used/cap threshold throttles deltaT for no
    //      stability reason;
    //   3. cap dependence: with maxCorrectors 5 the signal has five levels;
    //   4. scheme dependence: Anderson converges in ~3 until it does not, which
    //      would be a cliff by construction.
    {
        const label n = 100;
        const scalar ratio = 1.3/1.1;

        Info<< nl << "CASE 16  is `correctors used` a robust governor signal?"
            << endl;

        // Returns correctors used (== budget if it never converged), and sets
        // `conv`. Covers both schemes and both beds.
        auto runStep = [&]
        (
            const scalar G, const label m, const bool clamped,
            const label budget, const scalar ctol, bool& conv
        ) -> label
        {
            const scalar gab = Foam::sqrt(G*ratio);
            const scalar gba = Foam::sqrt(G/ratio);
            const scalar floorVal = 0.35;

            scalarField x(2*n, clamped ? 3.0 : 1.0);
            scalarField xStar(2*n, clamped ? 1.0 : 0.0);
            if (clamped)
            {
                forAll(xStar, c) { if ((c % 2) == 0) xStar[c] = 0.10; }
            }

            andersonRelaxation aa("g", m, 1.0);
            aitkenRelaxation rJ("g", 1.0, wMinBound, wMaxBound);
            aa.reset(); rJ.reset();
            conv = false;

            for (label k = 0; k < budget; ++k)
            {
                scalarField r(2*n);
                forAll(r, c)
                {
                    r[c] = (c < n)
                        ? (xStar[c] + gab*(x[n+c] - xStar[n+c])) - x[c]
                        : (xStar[c] - gba*(x[c-n] - xStar[c-n])) - x[c];
                }

                const scalarField xOld(x);

                if (m == 0)
                {
                    const scalar w = rJ.omega(r);
                    forAll(x, c) { x[c] += w*r[c]; }
                }
                else
                {
                    aa.correct(x, r);
                }

                if (clamped)
                {
                    forAll(x, c) { x[c] = max(x[c], floorVal); }
                }

                scalar nrm = 0;
                forAll(x, c)
                {
                    nrm = clamped ? max(nrm, mag(x[c] - xOld[c]))
                                  : max(nrm, mag(x[c]));
                }
                if (nrm < ctol) { conv = true; return k + 1; }
                if (nrm > 1e12) { return budget; }
            }
            return budget;
        };

        // --- 1 and 4: LEAD, for both schemes and both beds -----------------
        //
        // THE METRIC: correctors used at the LAST stiffness that still
        // converged, as a fraction of the budget. Near 1 means the count only
        // rises as it fails -- a cliff, no lead, and nothing to govern on.
        // Well below 1 means a threshold placed there catches the failure
        // before it happens.
        Info<< nl << "        LEAD: correctors at the last converging"
            << " stiffness, /budget 20" << endl;

        for (label variant = 0; variant < 4; ++variant)
        {
            const label m = (variant % 2) ? 2 : 0;      // Aitken | Anderson
            const bool clamped = (variant / 2) == 1;

            label lastUsed = -1;
            scalar lastG = 0;
            label firstFailUsed = -1;

            for (label gi = 0; gi < 16; ++gi)
            {
                const scalar G = 0.02*Foam::pow(1.5, scalar(gi));
                bool conv = false;
                const label used = runStep(G, m, clamped, 20, tol, conv);
                if (conv) { lastUsed = used; lastG = G; }
                else if (firstFailUsed < 0) { firstFailUsed = used; }
            }

            Info<< "          " << (m ? "Anderson m=2" : "Aitken      ")
                << (clamped ? " clamped: " : " linear : ");
            if (lastUsed < 0)
            {
                Info<< "never converged -- no lead measurable" << endl;
            }
            else
            {
                Info<< "last converged at G " << lastG << " using "
                    << lastUsed << "/20 (" << scalar(lastUsed)/20.0 << ")";
                if (firstFailUsed >= 0)
                {
                    Info<< ", first failure used " << firstFailUsed << "/20";
                }
                Info<< endl;
            }
        }

        // --- 2: TOLERANCE DEPENDENCE --------------------------------------
        //
        // If the count at the stability boundary moves with `tolerance`, a
        // FIXED used/cap threshold is not a stability signal at all -- it is
        // partly a measure of how tight the user set the criterion.
        Info<< nl << "        TOLERANCE DEPENDENCE (Aitken, linear bed):"
            << endl;
        for (label ti = 0; ti < 3; ++ti)
        {
            const scalar ctol = Foam::pow(10.0, -6.0 - 2.0*scalar(ti));
            label lastUsed = -1;
            for (label gi = 0; gi < 16; ++gi)
            {
                const scalar G = 0.02*Foam::pow(1.5, scalar(gi));
                bool conv = false;
                const label used = runStep(G, 0, false, 20, ctol, conv);
                if (conv) lastUsed = used;
            }
            Info<< "          tolerance " << ctol << " : last converging step"
                << " used " << lastUsed << "/20" << endl;
        }

        // --- 3: CAP DEPENDENCE --------------------------------------------
        Info<< nl << "        CAP DEPENDENCE (Aitken, linear bed):" << endl;
        for (label ci = 0; ci < 3; ++ci)
        {
            const label budget = (ci == 0) ? 5 : (ci == 1) ? 20 : 50;
            label lastUsed = -1;
            for (label gi = 0; gi < 16; ++gi)
            {
                const scalar G = 0.02*Foam::pow(1.5, scalar(gi));
                bool conv = false;
                const label used = runStep(G, 0, false, budget, tol, conv);
                if (conv) lastUsed = used;
            }
            Info<< "          cap " << budget << " : last converging step used "
                << lastUsed << "/" << budget << " ("
                << scalar(lastUsed)/scalar(budget) << ")" << endl;
        }

        Info<< nl << "        -> ROBUST requires: lead well below 1.0 on ALL"
            << " four variants, a" << nl
            << "           used/cap fraction that does NOT move with tolerance,"
            << " and one that" << nl
            << "           does not move with the cap. Any of those failing"
            << " means a fixed" << nl
            << "           threshold is not portable and the signal needs"
            << " normalising." << endl;
    }

    // ---- CASE 17: does ANY signal have lead UNDER ANDERSON? --------------
    //
    // CASE 16 killed `correctors used`: no lead on any variant, opposite
    // signatures for the two schemes, and it moves with both tolerance and the
    // cap. The proxies looked better because they are CONTINUOUS -- omega ran
    // 0.97 -> 0.12 and masked rho 0.15 -> 0.38 -> 35 across the sweep, which is
    // real lead.
    //
    // But that was measured under AITKEN. The worry is structural: a method
    // that either works or does not gives a cliff in every diagnostic, because
    // right up to the boundary it is converging comfortably. If masked rho is
    // flat under Anderson then NO pre-emptive signal exists, and the honest
    // architecture is that retryStep handles it -- the margin governor being a
    // device for Aitken's graceful-but-slow degradation, not a universal need.
    //
    // This measures rho (over unclamped cells) at the last converging and the
    // first failing stiffness, for both schemes.
    {
        const label n = 100;
        const scalar ratio = 1.3/1.1;
        const scalar floorVal = 0.35;

        Info<< nl << "CASE 17  does masked rho keep its lead under Anderson?"
            << endl;
        Info<< "        scheme          rho at last converged -> at first"
            << " failure" << endl;

        for (label variant = 0; variant < 4; ++variant)
        {
            const label m = (variant % 2) ? 2 : 0;
            const bool clamped = (variant / 2) == 1;

            scalar rhoLastOK = -1, rhoFirstBad = -1;

            for (label gi = 0; gi < 16; ++gi)
            {
                const scalar G = 0.02*Foam::pow(1.5, scalar(gi));
                const scalar gab = Foam::sqrt(G*ratio);
                const scalar gba = Foam::sqrt(G/ratio);

                scalarField x(2*n, clamped ? 3.0 : 1.0);
                scalarField xStar(2*n, clamped ? 1.0 : 0.0);
                if (clamped)
                {
                    forAll(xStar, c) { if ((c % 2) == 0) xStar[c] = 0.10; }
                }

                andersonRelaxation aa("g", m, 1.0);
                aitkenRelaxation rJ("g", 1.0, wMinBound, wMaxBound);
                aa.reset(); rJ.reset();

                scalar rhoMax = 0, rPrev = -1;
                bool conv = false;

                for (label k = 0; k < 20; ++k)
                {
                    scalarField r(2*n);
                    forAll(r, c)
                    {
                        r[c] = (c < n)
                            ? (xStar[c] + gab*(x[n+c] - xStar[n+c])) - x[c]
                            : (xStar[c] - gba*(x[c-n] - xStar[c-n])) - x[c];
                    }

                    // masked: skip cells sitting on the clamp (CASE 15)
                    scalar rN = 0;
                    forAll(r, c)
                    {
                        if (clamped && x[c] <= floorVal*(1 + SMALL)) continue;
                        rN += r[c]*r[c];
                    }
                    rN = Foam::sqrt(rN);
                    if (rPrev > VSMALL) { rhoMax = max(rhoMax, rN/rPrev); }
                    rPrev = rN;

                    const scalarField xOld(x);
                    if (m == 0)
                    {
                        const scalar w = rJ.omega(r);
                        forAll(x, c) { x[c] += w*r[c]; }
                    }
                    else
                    {
                        aa.correct(x, r);
                    }
                    if (clamped)
                    {
                        forAll(x, c) { x[c] = max(x[c], floorVal); }
                    }

                    scalar nrm = 0;
                    forAll(x, c)
                    {
                        nrm = clamped ? max(nrm, mag(x[c] - xOld[c]))
                                      : max(nrm, mag(x[c]));
                    }
                    if (nrm < tol) { conv = true; break; }
                    if (nrm > 1e12) break;
                }

                if (conv) { rhoLastOK = rhoMax; }
                else if (rhoFirstBad < 0) { rhoFirstBad = rhoMax; }
            }

            Info<< "          " << (m ? "Anderson m=2" : "Aitken      ")
                << (clamped ? " clamped: " : " linear : ")
                << rhoLastOK << "  ->  " << rhoFirstBad;
            if (rhoLastOK > 0 && rhoFirstBad > 0)
            {
                const scalar sep = rhoFirstBad/rhoLastOK;
                Info<< "   separation x" << sep
                    << (sep > 1.3 ? "  (LEAD)" : "  (no useful lead)");
            }
            Info<< endl;
        }
        Info<< "        -> a separation near 1 means the diagnostic looks the"
            << " same on the last" << nl
            << "           good step and the first bad one: no warning, and"
            << " retryStep -- not a" << nl
            << "           pre-emptive governor -- is the right mechanism for"
            << " that scheme." << endl;
    }

    Info<< nl << "=== end ===" << nl << endl;
    return 0;
}

// ************************************************************************* //
