#!/usr/bin/env bash
# PreToolUse guard for Edit / Write / NotebookEdit.
#
# TIER 1 (block): the file is GENERATED. Editing it is silently undone on the
#   next run, and the edit never reaches the thing that produces it. Each block
#   names the file that should be edited instead.
# TIER 2 (warn, never block): the file is hand-authored but is one of the two
#   dictionaries where a LITERAL belongs in configuration/config as a $key.
#   Heuristic, so it must not be able to stop work.
#
# The Tier-1 list is DERIVED FROM .gitignore, which is where this project
# already records what is generated and why. Do not add a path here that
# .gitignore does not also treat as generated -- that is how the two copies
# start to disagree (G1).
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
hooks_disabled && exit 0

payload="$(read_payload)"
tool="$(json_field "$payload" tool_name)"
path="$(json_field "$payload" tool_input.file_path)"
[ -z "${path:-}" ] && exit 0

rel="${path#"$HOOK_DIR"/}"

# --- Tier 1 -----------------------------------------------------------------

# A system/fvSolution is generated ONLY where the thing that overwrites it
# exists: Allrun does `cp system/fvSolution-foam system/fvSolution` (or -petsc,
# chosen by `matrixSolver`). So test for the SOURCE, never for the path shape.
#
# Two families of hand-authored fvSolution would otherwise be blocked wrongly,
# and both were found by running this guard against the real tree:
#   * needleDBD -- no -foam/-petsc pair, because one ePotential block carries
#     both backends and $ePotentialSolver picks one. .gitignore un-ignores it.
#   * every verification/ bed -- unit beds with no Allrun copy step at all.
# "Block only when the overwriting mechanism demonstrably exists" is the rule;
# a guard that fires on a name rather than on a mechanism is the same defect as
# a diagnostic that cannot fail.
if [[ "$rel" == */system/fvSolution ]]; then
    d="$(dirname "$path")"
    if [ -e "$d/fvSolution-foam" ] || [ -e "$d/fvSolution-petsc" ]; then
        hook_log BLOCK generated-fvSolution "$rel"
        block "BLOCKED: $rel is GENERATED.
Allrun copies system/fvSolution-foam (or -petsc, chosen by \`matrixSolver\`)
over it, so this edit is discarded at the next run.
Edit instead:  $(dirname "$rel")/fvSolution-foam$([ -e "$d/fvSolution-petsc" ] && printf '\n               %s/fvSolution-petsc' "$(dirname "$rel")")"
    fi
fi

if [[ "$rel" == */system/changeDictionaryDict ]]; then
    hook_log BLOCK generated-changeDictionaryDict "$rel"
    block "BLOCKED: $rel is GENERATED.
Allrun copies etc/changeDictionary.<region> over it.
Edit instead:  the matching etc/changeDictionary.* in this case."
fi

# A case's 0/ is rewritten (`rm -rf 0 && mkdir -p 0`, then the two generators)
# ONLY where the generator inputs exist. Same principle as fvSolution above:
# test for the mechanism, not the path. A verification/ unit bed keeps its
# hand-authored 0/n and 0/phi and has neither configuration/ nor 0.orig/, so it
# must not be caught. 0.orig/ itself is the seed and is never generated.
if [[ "$rel" == */0/* ]] && [[ "$rel" != */0.orig/* ]]; then
    case_dir="${path%/0/*}"
    if [ -d "$case_dir/configuration" ] || [ -d "$case_dir/0.orig" ]; then
        hook_log BLOCK generated-initial-field "$rel"
        block "BLOCKED: $rel is a GENERATED initial field.
\`rm -rf 0 && mkdir -p 0\` then plasmaSetupBoundaries + plasmaCreateSpeciesFields
rewrite this whole directory, so this edit is discarded at the next run.
Edit instead, depending on what you are changing:
  configuration/boundaries   what a surface IS      (plasmaSetupBoundaries)
  configuration/config       a value, as \$key
  0.orig/<field>             the hand-authored seed, where the case has one"
    fi
fi

# Per-region stubs: ten-line #include shims that splitMeshRegions overwrites.
if [[ "$rel" == */constant/*/plasmaSpeciesProperties ]] ||
   [[ "$rel" == */constant/*/plasmaTransportProperties ]] ||
   [[ "$rel" == */constant/*/photoionizationProperties ]]; then
    hook_log BLOCK generated-region-stub "$rel"
    block "BLOCKED: $rel is a per-region STUB.
tools/plasmaSetupRegions.sh regenerates it, and splitMeshRegions overwrites it
with an empty dummy. Its only content is an #include of the parent.
Edit instead:  $(echo "$rel" | sed -E 's#/constant/[^/]+/#/constant/#')"
fi

# --- Tier 2 -----------------------------------------------------------------
#
# The 2026-09-10 defect was a LITERAL maxDeltaT typed into
# system/plasmaSimulationControls while configuration/config was the owner. The
# edit itself is legitimate -- these files are hand-authored and tracked -- so
# this only asks the question. It can never block.
if [[ "$rel" == */system/plasmaSimulationControls ]] ||
   [[ "$rel" == */system/controlDict ]]; then
    new="$(json_field "$payload" tool_input.new_string)"
    [ -z "${new:-}" ] && new="$(json_field "$payload" tool_input.content)"
    # a dictionary entry whose value is a bare number, with no $ indirection
    if printf '%s' "${new:-}" | grep -qE '^[[:space:]]*[A-Za-z_][A-Za-z0-9_]*[[:space:]]+-?[0-9][^;]*;' &&
       ! printf '%s' "${new:-}" | grep -q '\$'; then
        hook_log WARN literal-in-generated-consumer "$rel"
        echo "NOTE: $rel is gaining a bare numeric literal.
Per G1, a value belongs in configuration/config and this dictionary references
it as \$key -- a setting stated in two places is a defect even when both copies
agree, because they will not agree later. If this value genuinely has only one
home, carry on.
Then prove the key is not dangling:  grep -rlF '\$'\"\$key\" system/ constant/" >&2
    fi
fi

exit 0
