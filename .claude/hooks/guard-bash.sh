#!/usr/bin/env bash
# PreToolUse guard for Bash.
#
# TIER 1 (block): `git add -A` / `git add .` / `git add <dir>` on a tree that
#   holds case output. Rule D4 forbids these verbatim -- "never `git add -A`,
#   never `git add <dir>` on a tree holding case output; re-read
#   `git diff --cached --name-only` and unstage what does not belong."
#
# TIER 2 (warn): a bare `wmake` while a solver is running. build-all.sh already
#   refuses this through check-no-running-solvers.sh; only a hand-run wmake gets
#   past, and the cost is a stale .so under a running solver -- which presents
#   as a startup SEGV that looks like a physics bug.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
hooks_disabled && exit 0

payload="$(read_payload)"
cmd="$(json_field "$payload" tool_input.command)"
[ -z "${cmd:-}" ] && exit 0

# --- Tier 1: blanket git add ------------------------------------------------
#
# Matched on the git subcommand, not anywhere in the line, so that
# `echo "never git add -A"` and `git add -A -- specific/path` are not caught.
if printf '%s' "$cmd" | grep -qE '(^|[;&|]|&&)[[:space:]]*git[[:space:]]+add[[:space:]]+(-A|--all|\.)([[:space:]]|$)'; then
    hook_log BLOCK git-add-all "$cmd"
    block "BLOCKED by rule D4: never \`git add -A\` or \`git add .\` here.

This tree holds 71 GB of case output beside the source. A blanket add is how
run artefacts, 1.3 GB of locally built PETSc, or a half-finished experiment get
committed -- and the staging is what nobody re-reads afterwards.

Stage EXPLICITLY instead, then verify:
    git add <path> [<path> ...]
    git diff --cached --name-only     # re-read this, unstage what does not belong

If you want everything a pathspec covers, name the pathspec:
    git add src/models/plasmaModels/plasmaTransport/"
fi

# --- Tier 2: wmake while a solver runs --------------------------------------
if printf '%s' "$cmd" | grep -qE '(^|[;&|]|&&)[[:space:]]*(wmake|wclean)([[:space:]]|$)'; then
    if pgrep -x soPlasmaFoam >/dev/null 2>&1 ||
       pgrep -x singleRegionElectrostaticFoam >/dev/null 2>&1 ||
       pgrep -x multiRegionElectrostaticFoam >/dev/null 2>&1; then
        hook_log WARN wmake-while-solver-running "$cmd"
        echo "NOTE: a solver is RUNNING and this is a bare \`wmake\`.

Rebuilding a library a running solver has mapped is how a stale .so and an
ABI mismatch get in: adding a virtual to a shared header shifts every later
vtable slot, and the result presents as a startup SEGV that reads like a
physics bug.

  ./check-no-running-solvers.sh    # what build-all.sh runs first
  ./build-all.sh                   # dependency-ordered, per-component logs

Decide the kill point deliberately before rebuilding." >&2
    fi
fi

exit 0
