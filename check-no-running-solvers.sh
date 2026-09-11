#!/usr/bin/env bash
#
# Refuse to build while a SoPLASMA executable is running.
#
# WHY THIS EXISTS
#
# wmake writes the rebuilt library straight over the path a running solver has
# mmap'd. In the common case wmake replaces the file by rename, so a running
# process keeps its old inode and survives -- but that is a property of how the
# linker happens to write, not a guarantee. When the file IS modified in place,
# the running process faults on a page it has not yet touched and dies with
# SIGBUS or a segfault, tens of minutes into a run, pointing at nothing that
# looks like the cause.
#
# It cost a scare on 2026-08-19: a rebuild was started while a 45-minute Co 2.5
# streamer was at step 230. It survived, but only by luck.
#
# For a USER the failure is worse than for a developer: they will not connect
# "I rebuilt in another terminal" with "my overnight run died at 3am", and the
# run is simply lost.
#
# Override with FORCE=1 when you really mean it (e.g. the only running solver
# is one you are about to kill).

set -u

# Executables this project installs. Matched by NAME, not by command line:
# `pkill -f soPlasmaFoam` does NOT match `soPlasmaFoam -parallel` reliably,
# which is its own trap.
_apps="soPlasmaFoam plasmaChemistry0D plasmaCreateSpeciesFields
       foamPlasmaCreateSpeciesFields testWallLoss testVibRelax
       testChemistryBackends"

_running=""
for _a in $_apps; do
    if pgrep -x "$_a" >/dev/null 2>&1; then
        _running="$_running $_a"
    fi
done

if [ -n "$_running" ]; then
    echo ""
    echo "  ============================================================"
    echo "  REFUSING TO BUILD: a SoPLASMA executable is running."
    echo "  ============================================================"
    echo ""
    echo "  Rebuilding overwrites the shared libraries these processes"
    echo "  have mapped into memory. They can die mid-run with SIGBUS or"
    echo "  a segfault that gives no hint of the real cause."
    echo ""
    for _a in $_running; do
        for _p in $(pgrep -x "$_a"); do
            _cwd=$(readlink "/proc/$_p/cwd" 2>/dev/null || echo "?")
            printf '    %-28s pid %-8s %s\n' "$_a" "$_p" "$_cwd"
        done
    done
    echo ""
    echo "  Wait for them to finish, or stop them first:"
    echo "      pkill -9 soPlasmaFoam; pkill -9 mpirun"
    echo ""
    echo "  To build anyway (you accept losing those runs):"
    echo "      FORCE=1 $0"
    echo ""
    if [ "${FORCE:-0}" != "1" ]; then
        exit 1
    fi
    echo "  FORCE=1 set -- building anyway."
    echo ""
fi

exit 0
