#!/usr/bin/env bash
# Shared helpers for the SoPlasma Claude Code hooks.
#
# WHY THESE EXIST AT ALL. This project's own measured conclusion is that
# "reminders do not fire, artifacts do" -- three separate memory entries covered
# the 2026-08-30 wrong-control failure and none of them prevented it; the
# on-disk COMPARE.md contract did. A hook is that same idea applied to the two
# mistakes that a rule in CLAUDE.md has repeatedly failed to stop.
#
# WHY THEY ARE NARROW. The prose rule says "NEVER hand-edit system/* or 0/*".
# Taken literally that would block 64 tracked controlDict, 60 fvSchemes, 55
# decomposeParDict, 41 plasmaSimulationControls, 28 blockMeshDict and 40
# fvSolution-foam -- every one hand-authored and legitimately edited. So these
# hooks do NOT implement the prose. They implement .gitignore, which already
# records precisely which files are GENERATED and why. One definition, per G1.

HOOK_DIR="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
HOOK_LOG="$HOOK_DIR/.claude/hook-events.log"

# THE KILL SWITCH. `touch .claude/hooks-off` disables every hook for as long as
# the file exists; `rm` re-enables. No JSON editing to get unstuck -- a guard
# you cannot turn off quickly is a guard that gets deleted permanently.
hooks_disabled() { [ -f "$HOOK_DIR/.claude/hooks-off" ]; }

# Every event is recorded, blocked or not, so the hook set can be tuned from
# evidence rather than from irritation. /save-state reports what fired.
hook_log() {  # tier  rule  detail
    printf '%s\t%s\t%s\t%s\n' "$(date -Is)" "$1" "$2" "$3" >> "$HOOK_LOG" 2>/dev/null || true
}

# Read the tool payload once; both guards need it.
read_payload() { cat; }

# Pull a field out of the PreToolUse JSON without assuming jq is installed.
json_field() {  # <json>  <dotted.path>
    python3 -c '
import json, sys
try:
    d = json.loads(sys.argv[1])
except Exception:
    sys.exit(0)
for k in sys.argv[2].split("."):
    if not isinstance(d, dict):
        sys.exit(0)
    d = d.get(k)
    if d is None:
        sys.exit(0)
print(d if isinstance(d, str) else json.dumps(d))
' "$1" "$2" 2>/dev/null
}

# Exit 2 makes Claude Code treat stderr as the blocking reason and show it to
# the model. The message must name the CORRECT action, not merely refuse: a
# guard that says only "no" gets worked around instead of obeyed.
block() {
    echo "$1" >&2
    exit 2
}
