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

    Info<< nl << "=== end ===" << nl << endl;
    return 0;
}

// ************************************************************************* //
