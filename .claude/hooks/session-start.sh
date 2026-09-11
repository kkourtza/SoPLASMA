#!/usr/bin/env bash
# SessionStart hook. Never blocks anything; its stdout becomes session context.
#
# THE POINT. Rule R0 exists because a compaction summary preserves the NARRATIVE
# of what was done and loses the INVENTORY of what exists, and on 2026-09-10
# that cost a proposal for a capability already implemented in three places and
# a hand-rolled case cold start that re-earned three traps a script already
# encoded. R0 asks me to re-read the inventory. This hook removes the step where
# I have to remember to.
#
# It states FACTS ONLY -- where to read and what git currently says. It must
# never summarise the state itself, because a summary here would be a second
# copy of PROGRESS.md and would go stale (G1, D2).
set -uo pipefail
ROOT="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "$ROOT" 2>/dev/null || exit 0

echo "=== SoPlasma / SoPhy — session orientation (R0) ==="
echo
echo "READ FIRST, in this order:"
echo "  1. PROGRESS.md          the live state: next action, in flight, blocked,"
echo "                          and the failed approaches not to retry"
echo "  2. docs/CAPABILITIES.md the inventory: what exists, the tooling, what was refuted"
echo "  (CLAUDE.md is already in context. Read docs/rules-postmortems.md only when"
echo "   a rule looks arbitrary — it is the measured evidence, not the rule.)"

if [ -f PROGRESS.md ]; then
    # tolerate markdown emphasis around the stamp (**Last updated: ...**)
    stamp=$(grep -m1 -ioE 'Last updated:[^*]*' PROGRESS.md 2>/dev/null | sed 's/[[:space:]]*$//' || true)
    [ -n "$stamp" ] && echo && echo "PROGRESS.md — $stamp"
    # Section 1 is capped at five lines by contract, so echoing it is not a
    # second copy; it is the pointer with its first line attached.
    awk '/^## 1\./{f=1;next} /^## 2\./{f=0} f' PROGRESS.md 2>/dev/null |
        grep -v '^[[:space:]]*$' | head -6 | sed 's/^/    /'
else
    echo
    echo "  !! PROGRESS.md IS MISSING — it is the one file this project orients from."
fi

echo
n=$(git status --short 2>/dev/null | wc -l)
echo "git: $n uncommitted entries on $(git branch --show-current 2>/dev/null)"
if [ "$n" -gt 0 ]; then
    new=$(git status --short 2>/dev/null | grep -c '^??' || true)
    [ "${new:-0}" -gt 0 ] && echo "     $new UNTRACKED — D4 calls these the urgent case: they appear in no diff"
    echo "     (\`git status --short\`; stage explicitly, never \`git add -A\`)"
fi

# Anything the guards caught since the last session, so the hook set can be
# tuned from evidence rather than from irritation.
LOG="$ROOT/.claude/hook-events.log"
if [ -s "$LOG" ]; then
    recent=$(tail -20 "$LOG" | wc -l)
    blocked=$(tail -20 "$LOG" | grep -c '^[^	]*	BLOCK' || true)
    echo
    echo "hooks: $recent recent events in .claude/hook-events.log (${blocked:-0} blocked)."
    echo "       If a BLOCK was wrong, say so — that entry gets demoted, not argued with."
fi
[ -f "$ROOT/.claude/hooks-off" ] && echo && echo "hooks: DISABLED (.claude/hooks-off exists; \`rm\` it to re-enable)"

exit 0
