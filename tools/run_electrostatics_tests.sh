#!/bin/bash
#------------------------------------------------------------------------------
# Run every ELECTROSTATICS-ONLY test case and check it against its declared
# reference. Fast by construction -- all of these are 2D plate or 1D column
# beds and the whole suite is seconds, so there is no excuse not to run it.
#
# WHY THIS EXISTS
#   Every one of these cases was DEAD for three weeks and nothing noticed.
#   fe827ae (2026-08-11) made `backgroundDensity` fatal in the Poisson coeffs
#   but left `reducedE = Emag/backgroundDensity` unguarded in
#   singleRegionPoisson. The only caller of setBackgroundDensity() is
#   plasmaSpecies, which does not exist in an electrostatics-only solver, so
#   every case died with SIGFPE. Found 2026-09-03, by accident, while testing
#   an unrelated boundary condition. Fixed in 82fc587.
#
#   The lesson is build-all.sh's, one level up: A CASE THAT IS NOT RUN IS NOT
#   TESTED, AND IT ROTS SILENTLY. The plasma tutorials are hours long and
#   cannot be run casually; these cannot be an excuse, because they are cheap.
#
# USAGE
#   tools/run_electrostatics_tests.sh          # all of them
#   Exit status is 0 only if EVERY case passed.
#------------------------------------------------------------------------------
set -u

HERE="$(cd "${0%/*}/.." && pwd)"
ES="$HERE/tutorials/electrostatics"
CHECK="$HERE/tools/check_series_stack.py"

if [ -z "${WM_PROJECT_DIR:-}" ]; then
    echo "ERROR: OpenFOAM is not sourced. Source its bashrc first." >&2
    exit 2
fi

fail=0
pass=0

report () { # name verdict detail
    printf '  %-52s %-6s %s\n' "$1" "$2" "${3:-}"
}

# --- 1. the analytic series stack for thinDielectricPotential ---------------
# Its own sweep already prints per-row PASS/FAIL and exits nonzero on any
# failure, including the negative (guard) and round-trip (restart) tests.
name="singleRegionElectrostaticFoam/thinDielectricSeriesStack"
if ( cd "$ES/$name" && ./Allrun-sweep > log.suite 2>&1 ); then
    n=$(grep -c 'PASS' "$ES/$name/log.suite" || true)
    report "$name" PASS "$n checks"
    pass=$((pass+1))
else
    report "$name" FAIL "see $ES/$name/log.suite"
    fail=$((fail+1))
fi
( cd "$ES/$name" && ./Allclean > /dev/null 2>&1 )

# --- 2. the three series-stack cases with an ANALYTIC reference -------------
for name in \
    multiRegionElectrostaticFoam/plate2D_timeVaryingBC_explicitBoundary \
    multiRegionElectrostaticFoam/plate2D_timeVaryingBC_implicitBoundary \
    multiRegionElectrostaticFoam/plate2D_surfCharge_timeVaryingBC_implicitBoundary
do
    d="$ES/$name"
    ( cd "$d" && ./Allclean > /dev/null 2>&1; ./Allrun-serial > log.suite 2>&1 )
    # The reference and the tolerance live in the case's own COMPARE.md and are
    # applied by check_series_stack.py -- NOT restated here, so there is exactly
    # one copy of every reference number.
    if out=$( cd "$d" && "$CHECK" 2>&1 ); then
        report "$name" PASS "$(echo "$out" | grep -oE 'relative error [0-9.e-]+' | head -1)"
        pass=$((pass+1))
    else
        report "$name" FAIL "$(echo "$out" | tail -3 | tr '\n' ' ')"
        fail=$((fail+1))
    fi
    # Leave the tree clean -- run outputs are untracked and would otherwise
    # bury real changes in `git status`. The log is kept for a failure.
    [ "$fail" -eq 0 ] && ( cd "$d" && ./Allclean > /dev/null 2>&1 )
done

# --- 3. cases with no analytic reference: they must merely RUN --------------
# No physics check exists for these, so the only claim made is "it completed
# without a fatal error or a signal". Said explicitly rather than dressed up as
# a validation.
for name in singleRegionElectrostaticFoam/plate2D_timeVaryingBC ; do
    d="$ES/$name"
    ( cd "$d" && ./Allclean > /dev/null 2>&1; ./Allrun-serial > log.suite 2>&1 )
    if grep -qE 'FOAM FATAL|--> FOAM Warning.*sigFpe|sigFpe::sigHandler|Floating point exception \(core dumped\)|sigSegv' \
        "$d"/log.suite "$d"/logs/* 2>/dev/null; then
        report "$name" FAIL "fatal or signal in the log"
        fail=$((fail+1))
    else
        report "$name" RAN "no analytic reference -- completion only"
        pass=$((pass+1))
        ( cd "$d" && ./Allclean > /dev/null 2>&1 )
    fi
done

echo
echo "  $pass ok, $fail failed"
[ "$fail" -eq 0 ] || echo "  A FAILURE HERE MEANS AN ELECTROSTATICS CASE IS BROKEN."
exit $((fail > 0))
