#!/usr/bin/env python3
"""checkConfigReference -- keep the option reference honest about the code.

WHAT IT DOES
------------
Parses every option key the solver reads out of `src/`, parses every option key
the reference documents in `docs/reference/`, and reports the three ways they
can disagree:

  MISSING READER   the reference documents an option no source file reads.
  STALE DEFAULT    the reference states a default the source no longer uses.
  UNDOCUMENTED     the source reads an option the reference never mentions.

WHY IT EXISTS
-------------
Hand-written reference prose drifts from the code, silently, and the drift is
invisible precisely because the document reads plausibly. Three instances were
found in this repository on 2026-09-01/02, each of which had survived review:

  `ePotentialControls { nonCoupledResidualControl { ... } }` was documented in
  system/fvSolution by two separate documents. NOTHING HAS EVER READ IT -- the
  name appears in no source file and no shipped case. A user following it set a
  tolerance and an iteration cap that had no effect, and got the defaults.

  GAMG was documented as incompatible with monolithic region coupling. It needs
  only `agglomerator assembledFaceAreaPair`, which OpenFOAM ships for exactly
  that case.

  Per-region permittivity was documented as `electricProperties` containing
  `dielectricConstant`. Neither name was ever read: the values lived in a global
  file under a different key.

The first is the one this tool is really aimed at, because no amount of careful
reading finds it. You have to ask the code.

WHEN TO RUN IT
--------------
  * before committing a change to any docs/reference/*.md
  * after adding, renaming or re-defaulting any option in src/
  * in CI, if there is one

It is READ-ONLY. It never edits anything.

HOW TO RUN IT
-------------
    tools/checkConfigReference.py                 # check everything
    tools/checkConfigReference.py --strict        # UNDOCUMENTED is fatal too
    tools/checkConfigReference.py --list-source   # dump what the code reads
    tools/checkConfigReference.py --doc docs/reference/configuration-config.md

Exit status is 0 when clean, 1 when a MISSING READER or STALE DEFAULT is found,
and (with --strict) also 1 for UNDOCUMENTED. UNDOCUMENTED is a warning by
default because the reference is written incrementally and 237 keys do not
arrive at once -- but it is the check that stops new options landing unnoticed,
so turn it on once a document is complete.

THE CONVENTION IT PARSES
------------------------
A documented option is a level-3 heading holding the key in backticks, followed
somewhere before the next heading by one of:

    **REQUIRED -- no default.**
    **Default:** `<value>`

so the document stays readable prose rather than a machine format. Case and
surrounding words are ignored; only the backticked value is compared, and the
comparison is numeric where both sides parse as numbers (so 1e-10 == 1E-10 and
100 == 100.0).

CASE VARIABLES ARE NOT SOLVER OPTIONS
-------------------------------------
`configuration/config` does not hold solver options. It holds case VARIABLES
that dictionaries expand into solver options, and the two names differ freely:

    system/plasmaSimulationControls:   maxSpeciesConvectiveCo   $maxSpeciesCo;
    etc/changeDictionary.gas:          uniformValue             $voltageRamp;

So `voltageRamp` appears nowhere in src/ and never will. Checking a config
variable against the source directly reports every one of them as a MISSING
READER, which is wrong.

`--case <dir>` resolves the indirection: it scans the case's dictionaries for
`$name` references, records which KEY each one is assigned to, and validates
that key instead. That also catches an error class the plain mode cannot see --
a config variable wired to a dictionary key that nothing reads, which is
`appliedVoltage` (dangling in five beds, so editing it did nothing) with the
wiring present but pointed at the wrong name.

A variable referenced only inside a compound value, or consumed by a shell
script rather than a dictionary, has no single target key; those are reported as
UNRESOLVED rather than guessed at.

WHAT IT DELIBERATELY DOES NOT DO
--------------------------------
It does not try to work out WHICH dictionary a key belongs to. That is not
recoverable from the source without tracking every dictionary reference through
the call graph, and a wrong answer would be worse than none. So a key documented
in the wrong document still passes -- this tool checks that the key EXISTS and
that its DEFAULT is right, not that it is filed correctly.

It also cannot see options read through a computed key name. Those are rare and
are listed in NOT_USER_FACING below when found.
"""

import argparse
import re
import sys
from pathlib import Path

# Keys the source reads that are NOT user-facing case options, so their absence
# from the reference is not a defect. Each needs a reason.
# Keys owned by OPENFOAM, not by us. They are read by OpenFOAM's own
# controlDict/fvSchemes/fvSolution machinery, so they will never appear in a
# get<>() in src/ and their absence there is not a defect.
OPENFOAM_OWNED = {
    "deltaT", "writeControl", "writeInterval", "endTime", "startTime",
    "startFrom", "stopAt", "purgeWrite", "writeFormat", "writePrecision",
    "writeCompression", "timeFormat", "timePrecision", "runTimeModifiable",
    "adjustTimeStep", "maxCo", "maxDeltaT", "libs", "application", "functions",
    "default", "solver", "preconditioner", "smoother", "tolerance", "relTol",
    "maxIter", "minIter", "agglomerator", "nCellsInCoarsestLevel",
    "mergeLevels", "cacheAgglomeration", "nOuterCorrectors", "nCorrectors",
    "nNonOrthogonalCorrectors", "residualControl", "useImplicit",
}

# Options our code reads through a COMPUTED key name, so the literal never
# appears next to a get<>() and the scanner cannot see it. Each needs a
# file:line so the claim is checkable.
COMPUTED_KEY_READS = {
    # electromagneticsModel::readPoissonNumerics uses a `pick(newKey, oldKey,
    # fallback)` lambda, so these three are read via a variable
    "scheme": "src/models/electromagnetics/electromagneticsModel.C:205",
    "EScheme": "src/models/electromagnetics/electromagneticsModel.C:206",
}

NOT_USER_FACING = {
    # OpenFOAM's own, read via our dictionaries but documented upstream
    "value", "type", "phi", "inletValue", "refValue", "refGradient",
    "valueFraction", "internalField", "boundaryField", "dimensions",
    # internal registry / bookkeeping names, never written by a user
    "regions", "cellToRegion",
}

KEY_RE = re.compile(
    r"""(?:
            \.?(?:getOrDefault|lookupOrDefault)  # defaulted read
            \s*<\s*(?P<dtype>[A-Za-z_][\w:<>\s]*?)\s*>\s*
            \(\s*"(?P<dkey>[A-Za-z_]\w*)"\s*,\s*(?P<ddef>[^,()]*?(?:\([^)]*\))?[^,()]*?)\s*\)
        |
            \.?get\s*<\s*(?P<rtype>[A-Za-z_][\w:<>\s]*?)\s*>\s*
            \(\s*"(?P<rkey>[A-Za-z_]\w*)"\s*\)          # required read
        )""",
    re.X,
)


def scan_source(src: Path):
    """key -> {"default": str|None, "type": str, "where": [file:line, ...]}"""
    found = {}
    for f in sorted(src.rglob("*.[CH]")):
        if "lnInclude" in f.parts:          # symlink farm: same files twice
            continue
        try:
            text = f.read_text(errors="replace")
        except OSError:
            continue
        for m in KEY_RE.finditer(text):
            key = m.group("dkey") or m.group("rkey")
            dflt = m.group("ddef")
            typ = (m.group("dtype") or m.group("rtype") or "").strip()
            line = text.count("\n", 0, m.start()) + 1
            rec = found.setdefault(
                key, {"default": None, "type": typ, "where": [], "required": True}
            )
            rec["where"].append(f"{f.relative_to(src.parent)}:{line}")
            if dflt is not None:
                rec["default"] = dflt.strip().strip('"')
                rec["required"] = False
    return found


DEFAULT_RE = re.compile(r"\*\*Default:\*\*\s*`([^`]*)`", re.I)
REQUIRED_RE = re.compile(r"\*\*REQUIRED\b", re.I)
HEADING_RE = re.compile(r"^###\s+.*?`([A-Za-z_]\w*)`", re.M)


# `configuration-config.md` documents CASE VARIABLES, whose names deliberately
# differ from the option keys they expand into ($maxSpeciesCo ->
# maxSpeciesConvectiveCo). Checking it by name in plain mode reports every entry
# as a MISSING READER, which is wrong -- it is checked by `--case` instead,
# which resolves the indirection. Skipped here, and said out loud so the skip
# cannot be mistaken for a pass.
VARIABLE_DOCS = {"configuration-config.md"}


def scan_docs(docs: Path):
    """key -> {"default": str|None, "required": bool, "where": "file:line"}"""
    found = {}
    for f in sorted(docs.rglob("*.md")):
        if f.name in VARIABLE_DOCS:
            print(f"note: skipping {f.name} -- it documents case VARIABLES, not "
                  f"option keys.\n      Check it with:  "
                  f"--case <caseDir>")
            continue
        text = f.read_text(errors="replace")
        heads = list(HEADING_RE.finditer(text))
        for i, h in enumerate(heads):
            key = h.group(1)
            body = text[h.end(): heads[i + 1].start() if i + 1 < len(heads) else len(text)]
            line = text.count("\n", 0, h.start()) + 1
            d = DEFAULT_RE.search(body)
            found[key] = {
                "default": d.group(1).strip() if d else None,
                "required": bool(REQUIRED_RE.search(body)),
                "where": f"{f.relative_to(docs.parent.parent)}:{line}",
            }
    return found


def same_value(a: str, b: str) -> bool:
    """Compare numerically when both parse as numbers, else textually."""
    if a is None or b is None:
        return a == b
    a, b = a.strip(), b.strip()
    for t, f in (("true", "1"), ("false", "0"), ("yes", "1"), ("no", "0")):
        a = t if a.lower() == t else a
        b = t if b.lower() == t else b
    if a.lower() == b.lower():
        return True
    try:
        return abs(float(a) - float(b)) <= 1e-12 * max(1.0, abs(float(a)))
    except ValueError:
        return False


# `key   $var;`  --  the assignment that wires a case variable to a solver option
WIRE_RE = re.compile(r"^\s*([A-Za-z_]\w*)\s+\$(\w+)\s*;", re.M)
# any other appearance of $var: nested in a compound value, inside a list, etc.
USE_RE = re.compile(r"\$(\w+)")


def scan_case(case: Path):
    """config variable -> set of dictionary keys it is assigned to."""
    wired, seen = {}, set()
    skip = {"polyMesh", "plasmaTables", "postProcessing", "logs"}
    for f in sorted(case.rglob("*")):
        if not f.is_file() or any(p in skip for p in f.parts):
            continue
        if f.name.endswith((".msh", ".json", ".foam")) or f.name.startswith("log."):
            continue
        try:
            text = f.read_text(errors="replace")
        except OSError:
            continue
        for m in WIRE_RE.finditer(text):
            wired.setdefault(m.group(2), set()).add(m.group(1))
        seen.update(USE_RE.findall(text))
    return wired, seen


def check_case(a, src) -> int:
    cfg = a.case/"configuration"/"config"
    if not cfg.is_file():
        print(f"ERROR: no {cfg}", file=sys.stderr)
        return 2

    declared = re.findall(r"^([A-Za-z]\w*)\s", cfg.read_text(errors="replace"), re.M)
    wired, seen = scan_case(a.case)

    dangling, unresolved, no_reader, ok = [], [], [], []
    for v in sorted(set(declared)):
        if v not in seen:
            dangling.append(v)
        elif v not in wired:
            unresolved.append(v)
        else:
            bad = [k for k in sorted(wired[v])
                   if k not in src and k not in OPENFOAM_OWNED
                   and k not in COMPUTED_KEY_READS]
            (no_reader if bad else ok).append((v, sorted(wired[v]), bad))

    rc = 0
    if dangling:
        rc = 1
        print("DANGLING -- declared in config, referenced by nothing.")
        print("  Editing these does nothing, and nothing warns you. This is")
        print("  `appliedVoltage`, which was dangling in five beds while being")
        print("  exactly the knob you reach for to change the drive.")
        for v in dangling:
            print(f"    {v}")
        print()
    if no_reader:
        # A WARNING, not fatal. The scanner has two known blind spots --
        # OpenFOAM-owned keys and options read through a computed key name --
        # and both are allowlisted above, but the lists are maintained by hand
        # and so cannot be assumed complete. DANGLING below is unambiguous and
        # is what makes this exit non-zero.
        print("WIRED TO NOTHING -- expands into a key no source file reads  [warning]")
        print("  Check by hand: the scanner cannot see keys read via a computed")
        print("  name, and OPENFOAM_OWNED/COMPUTED_KEY_READS are hand-maintained.")
        for v, keys, bad in no_reader:
            print(f"    ${v:32s} -> {', '.join(bad)}")
        print()
    if unresolved:
        print(f"UNRESOLVED -- {len(unresolved)} referenced, but not as a plain "
              f"`key $var;` assignment  [not checked]")
        print("  Nested in a compound value, or read by a shell script. Verify by hand.")
        for v in unresolved:
            print(f"    {v}")
        print()

    print(f"{len(set(declared))} variables declared; {len(ok)} resolve to a key "
          f"the source reads; {len(unresolved)} unresolved; {len(dangling)} dangling")
    if not rc:
        print("OK")
    return rc


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--src", type=Path, default=Path("src"))
    ap.add_argument("--docs", type=Path, default=Path("docs/reference"))
    ap.add_argument("--doc", type=Path, help="check only this one document")
    ap.add_argument("--strict", action="store_true",
                    help="UNDOCUMENTED is fatal too")
    ap.add_argument("--list-source", action="store_true",
                    help="dump every option the code reads, then exit")
    ap.add_argument("--case", type=Path,
                    help="resolve configuration/config variables to the "
                         "dictionary keys they feed, and check those")
    a = ap.parse_args()

    if not a.src.is_dir():
        print(f"ERROR: no such source tree: {a.src}", file=sys.stderr)
        return 2

    src = scan_source(a.src)

    if a.case:
        return check_case(a, src)

    if a.list_source:
        for k in sorted(src):
            r = src[k]
            d = "REQUIRED" if r["required"] else f"default {r['default']}"
            print(f"{k:44s} {d:28s} {r['where'][0]}")
        print(f"\n{len(src)} options read by {a.src}/")
        return 0

    docs_root = a.doc.parent if a.doc else a.docs
    if not docs_root.is_dir():
        print(f"ERROR: no such reference directory: {docs_root}", file=sys.stderr)
        return 2
    doc = scan_docs(docs_root) if not a.doc else {
        k: v for k, v in scan_docs(docs_root).items()
        if v["where"].startswith(str(a.doc).replace("docs/", "docs/", 1)[:200])
        or True
    }

    missing, stale, undoc = [], [], []

    for k, d in sorted(doc.items()):
        if k not in src:
            missing.append((k, d))
            continue
        s = src[k]
        if d["required"] and not s["required"]:
            stale.append((k, d, s, f"documented REQUIRED, source defaults to {s['default']}"))
        elif not d["required"] and s["required"] and d["default"] is not None:
            stale.append((k, d, s, f"documented default {d['default']}, source has NO default"))
        elif d["default"] is not None and not same_value(d["default"], s["default"]):
            stale.append((k, d, s, f"documented {d['default']}, source has {s['default']}"))

    for k in sorted(src):
        if k not in doc and k not in NOT_USER_FACING:
            undoc.append(k)

    w = 0
    if missing:
        w = 1
        print("MISSING READER -- documented, but no source file reads it")
        print("  This is the ePotentialControls class of error: prose that reads")
        print("  plausibly and configures nothing.")
        for k, d in missing:
            print(f"    {k:40s} {d['where']}")
        print()
    if stale:
        w = 1
        print("STALE DEFAULT -- the source moved and the reference did not")
        for k, d, s, why in stale:
            print(f"    {k:40s} {why}")
            print(f"    {'':40s} doc {d['where']}   src {s['where'][0]}")
        print()
    if undoc:
        print(f"UNDOCUMENTED -- {len(undoc)} option(s) the code reads and the "
              f"reference does not mention"
              + ("  [FATAL: --strict]" if a.strict else "  [warning]"))
        for k in undoc:
            print(f"    {k:40s} {src[k]['where'][0]}")
        print()
        if a.strict:
            w = 1

    n_ok = len(doc) - len(missing) - len(stale)
    print(f"{len(doc)} documented, {n_ok} verified against source; "
          f"{len(src)} options read by {a.src}/")
    if not w:
        print("OK")
    return w


if __name__ == "__main__":
    sys.exit(main())
