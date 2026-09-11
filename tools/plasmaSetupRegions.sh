#!/usr/bin/env bash
#
# plasmaSetupRegions.sh -- give every mesh region the per-region files it needs.
#
# Run this AFTER splitMeshRegions, from the case directory:
#
#     splitMeshRegions -cellZones -overwrite
#     tools/plasmaSetupRegions.sh
#
# It does two things:
#
#   1. NUMERICS -- writes system/<region>/{fvSchemes,fvSolution} as ten-line
#      stubs that #include the case-wide files. Nothing is duplicated.
#
#   2. MATERIALS -- creates constant/<region>/electricalProperties if it does
#      not exist, and NEVER touches one that does. For a `gas` the file is
#      complete (epsilonR 1.0 is genuinely right). For a `dielectric` the value
#      is left COMMENTED OUT and the run is listed as needing your attention --
#      see the reasoning at writeElectricalProperties below.
#
#
# WHY THIS EXISTS
#
# OpenFOAM reads fvSchemes and fvSolution PER REGION, from
# system/<region>/, and there is NO FALLBACK to the case-wide system/ copy.
# Measured 2026-09-02: with system/gas/fvSchemes absent, any tool that builds
# the gas mesh dies with
#     --> FOAM FATAL ERROR: cannot find file ".../system/gas/fvSchemes"
# It never looks at system/fvSchemes.
#
# splitMeshRegions knows this and writes a file -- but an EMPTY one, announcing
# `Writing dummy "gas/fvSchemes"`. Empty is worse than missing, because it fails
# late: `checkMesh -region gas` is perfectly happy with it, and the first actual
# scheme lookup then dies with
#     Entry 'grad(...)' not found in dictionary "system/gas/fvSchemes/gradSchemes"
#
# The usual answer -- what every other case in this repository does -- is to
# `cp` the case-wide files into each region directory and `sed` the #include
# depth from "../" to "../../". That means N copies of the same content, drifting
# apart the moment someone edits one of them, plus a path rewrite that fails
# silently if the sed pattern misses.
#
# WHAT THIS DOES INSTEAD: writes a ten-line STUB per region that #includes the
# case-wide file. One place to edit, no duplication, nothing to keep in sync,
# and no include-depth surgery -- a relative #include resolves against the
# directory of the file containing it, so "../fvSchemes" from system/<region>/
# lands on system/fvSchemes, and the "../configuration/config" inside THAT file
# still resolves from system/. Verified through both levels: the $var
# substitutions expand correctly in the region's view of the dictionary.
#
# MUST RUN AFTER THE SPLIT. splitMeshRegions OVERWRITES whatever is there with
# its dummy -- measured, the stub's md5 changes -- so writing the stubs first
# does not work.
#
# Idempotent: safe to run again.

set -eu

if [ ! -f system/controlDict ]; then
    echo "ERROR: run this from a case directory (no system/controlDict here)." >&2
    exit 1
fi

for f in fvSchemes fvSolution; do
    if [ ! -f "system/$f" ]; then
        echo "ERROR: system/$f does not exist -- the stubs would point at nothing." >&2
        exit 1
    fi
done

# The regions come from constant/regionProperties, the ONE place they are
# declared. foamListRegions reads exactly that file, so this cannot disagree
# with the rest of the case.
regions=$(foamListRegions 2>/dev/null || true)

if [ -z "$regions" ]; then
    echo "plasmaSetupRegions: no regions in constant/regionProperties -- nothing to do."
    echo "  (a single-region case reads system/fvSchemes directly; this is correct)"
    exit 0
fi

# WHICH KIND IS EACH REGION? foamListRegions flattens every list and loses the
# kind, and `foamDictionary -entry regions/dielectric` does not work because
# `regions` is a flat LIST, not a sub-dictionary:
#     regions   ( gas ( gas ) dielectric ( dielectric ) );
# So parse it: the tokens alternate <kind> ( <name> ... ) <kind> ( ... ).
kindOf() {
    foamDictionary constant/regionProperties -entry regions 2>/dev/null \
    | sed -e 's/^ *regions *//' -e 's/;$//' \
    | awk -v want="$1" '
        {
            # strip the outermost parentheses, then walk the token stream
            sub(/^[[:space:]]*\(/, ""); sub(/\)[[:space:]]*$/, "")
            kind = ""; depth = 0
            for (i = 1; i <= NF; i++) {
                t = $i
                if (t == "(")      { depth++ }
                else if (t == ")") { depth--; if (depth == 0) kind = "" }
                else if (depth == 0) { kind = t }
                else if (t == want) { print kind; exit }
            }
        }'
}

# Files the user must still fill in, reported together at the end.
needsAttention=""

writeElectricalProperties() {
    region="$1"
    kind="$2"
    target="constant/$region/electricalProperties"

    if [ -f "$target" ]; then
        echo "  $target  (exists, left alone)"
        return 0
    fi

    mkdir -p "constant/$region"

    if [ "$kind" = "farField" ]; then
        # A FICTITIOUS AIR/VACUUM REGION. epsilonR 1.0 is not a placeholder
        # here, it is the definition of the kind, and the solver defaults to it
        # anyway -- so this file is written for discoverability, and deleting it
        # changes nothing.
        cat > "$target" <<EOF
/*--------------------------------*- C++ -*----------------------------------*\\
| GENERATED by tools/plasmaSetupRegions.sh -- no action needed.               |
| Regenerated only if you DELETE it; an existing file is never overwritten.   |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/$region";
    object      electricalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Region \`$region\`, declared as kind \`farField\` in constant/regionProperties.
//
// A FICTITIOUS AIR/VACUUM REGION. It solves ONLY Poisson -- no species, no
// chemistry, no electron energy, no photoionization -- so it gives the
// electrostatics a large domain at almost no cost, letting the gas region be
// small.
//
// epsilonR 1.0 is the DEFINITION of this kind, not a placeholder, and the
// solver defaults to it whether or not this file exists. It is written out for
// discoverability; deleting it changes nothing.
//
// VERIFIED 2026-09-02 on the two-region analytic bed: at epsilonR 1.0 the
// coupled interface reduces to plain continuity of V and dV/dn, and the
// interface potential comes out at the SINGLE-MEDIUM value (0.5 of the applied
// voltage at the midpoint) to a relative error of 0.00e+00. The interface is
// exactly transparent.
epsilonR    1.0;

// WHAT YOU STILL HAVE TO GET RIGHT, on the GAS side of the interface, in
// etc/changeDictionary.<gasRegion>. The patch is named \`<gas>_to_$region\`.
//
//   1. ePotential: KEEP \`coupledElectricPotential\`, with
//      \`surfCharge none; surfChargeNbr none;\`.
//      Do NOT use zeroGradient/Neumann there: that forces E.n = 0 and
//      INSULATES the plasma region from the far field, destroying the reason
//      for having it.
//
//   2. Species: an OUTFLOW condition, e.g.
//          type inletOutlet; phi particleFlux_<sp>; inletValue uniform <bg>;
//      This also removes secondary emission, since SEE lives in the wall-BC
//      family. The solver CHECKS this one: a plasmaWallBC with
//      \`enableSurfaceCharging true\` on a farField interface is fatal, because
//      it would invent a dielectric surface in the middle of the gas and shield
//      the field the region was added to resolve.
//
//   3. THE ONE PHYSICAL CONSTRAINT: no space charge exists here, since Poisson
//      in a Poisson-only region has no rho on the right-hand side. Put the
//      interface where the charge density is negligible and CHECK it stays so.
//      Photoionization is truncated at the gas boundary too, but that error
//      decays as exp(-d/L) with L ~ 1 mm in atmospheric air, so a few mm of
//      clearance makes it negligible -- unlike Poisson, which is long-range and
//      is the whole reason this region exists.

// ************************************************************************* //
EOF
        echo "  $target  (farField: epsilonR 1.0, complete)"
        return 0
    fi

    if [ "$kind" = "gas" ]; then
        # A LIVE VALUE, because for a gas 1.0 is not a placeholder: at
        # atmospheric density a gas is a vacuum to within a few parts in 10^4.
        # The solver defaults to exactly this when the file is absent, so
        # writing it changes nothing -- it is here so that a multi-region case
        # does not have a properties file for only some of its regions.
        cat > "$target" <<EOF
/*--------------------------------*- C++ -*----------------------------------*\\
| GENERATED by tools/plasmaSetupRegions.sh -- but yours to edit.              |
| Regenerated only if you DELETE it; an existing file is never overwritten.   |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/$region";
    object      electricalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Region \`$region\`, declared as kind \`gas\` in constant/regionProperties.
//
// 1.0 is not a placeholder here -- at atmospheric density a gas is a vacuum to
// within a few parts in 10^4, and the solver defaults to this when the file is
// absent. Change it only if you know why.
epsilonR    1.0;

// ************************************************************************* //
EOF
        echo "  $target  (gas: epsilonR 1.0, complete)"
        return 0
    fi

    # A DIELECTRIC GETS NO VALUE, DELIBERATELY.
    #
    # It is tempting to write `epsilonR 1.0;` as a placeholder to be edited
    # later. Do not: for a dielectric 1.0 is not a placeholder, it is VACUUM.
    # The case would run, converge, and give a plausible WRONG answer with the
    # barrier electrically invisible -- and in a DBD the barrier is the whole
    # point, since it is the surface charge on it that shields the gap and makes
    # the discharge self-limiting. Nothing in the log would say so.
    #
    # A guessed-but-plausible value such as 3.0 is no better: it is a made-up
    # number for someone else's material, and a run that completes on a
    # fabricated permittivity is exactly the kind of result nobody re-checks.
    #
    # So the entry is COMMENTED OUT. The file exists, is easy to find, and
    # explains itself; and until you uncomment it the solver stops with an error
    # that names this exact path. Failing to start beats starting and lying.
    cat > "$target" <<EOF
/*--------------------------------*- C++ -*----------------------------------*\\
| GENERATED by tools/plasmaSetupRegions.sh -- ACTION REQUIRED.                |
| Regenerated only if you DELETE it; an existing file is never overwritten.   |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/$region";
    object      electricalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Region \`$region\`, declared as kind \`dielectric\` in constant/regionProperties.
//
// >>> UNCOMMENT THE LINE BELOW AND SET YOUR MATERIAL'S VALUE. <<<
//
// It is left commented rather than defaulted because there is no safe DEFAULT.
// A default of 1 would be vacuum, so a real barrier would go electrically
// invisible while the case ran, converged and looked entirely fine. In a DBD the
// barrier is the entire point -- surface charge on it shields the gap and makes
// the discharge self-limiting -- and nothing in the log would tell you it had
// gone.
//
// SETTING epsilonR 1.0 DELIBERATELY IS A DIFFERENT THING, AND IS VALID. A
// \`dielectric\` region at 1.0 is a FICTITIOUS vacuum/air region: it solves
// Laplace and nothing else. That is how you extend the electrostatic domain --
// wrap a small plasma region in a large Poisson-only one so the far field is
// resolved without paying for species transport, chemistry or the electron
// energy equation out there. See the note at the bottom of this file.
//
// A guessed-but-plausible number would be no better: a run that completes on a
// fabricated permittivity is exactly the kind of result nobody re-checks. So
// the solver stops instead, with an error naming this file.
//
// RELATIVE PERMITTIVITY OF SOME COMMON BARRIER MATERIALS (approximate, room
// temperature, low frequency -- check your own material and cite it):
//
//     quartz / fused silica     3.8
//     borosilicate glass        4.6
//     soda-lime glass           7.0
//     alumina (Al2O3, 96%)      9.0
//     PTFE (Teflon)             2.1
//     PMMA (plexiglass)         3.0
//     polyimide (Kapton)        3.4
//     mica                      6.0
//     barium titanate         >1000   (a ferroelectric -- also non-linear)
//
// WHY IT MATTERS TWICE. It sets how the applied voltage DIVIDES between gap and
// barrier, the two being capacitors in series, so a higher epsilonR puts more of
// the voltage across the gas. And with the barrier thickness it fixes the
// barrier capacitance, hence how much surface charge is needed to extinguish
// the discharge.
//
// ---------------------------------------------------------------------------
// USING THIS REGION AS A FICTITIOUS FAR FIELD (epsilonR 1.0)
// ---------------------------------------------------------------------------
// A dielectric region solves ONLY fvm::laplacian(epsilon, ePotential) -- no
// species, no chemistry, no electron energy. So a large region at epsilonR 1.0
// is an air/vacuum Poisson domain that costs almost nothing, and lets the
// plasma region be small. Three things to get right:
//
//   1. THE INTERFACE STAYS `coupledElectricPotential`. Do NOT make it
//      zeroGradient/Neumann: that forces E.n = 0 and INSULATES the plasma
//      region from the far field, destroying the reason for having it. The
//      coupled BC weights each side by epsilon*deltaCoeffs, so with equal
//      epsilon and zero surface charge it reduces to plain continuity of V and
//      dV/dn -- i.e. the interface becomes transparent, which is what you want.
//      Set `surfCharge none; surfChargeNbr none;` on it.
//
//   2. NO SURFACE CHARGE, and this is already the default. Surface charge
//      accumulates ONLY where a species patch field is a plasmaWallBC AND has
//      `enableSurfaceCharging true` (plasmaTransport::updateSurfaceCharge does
//      a dynamic_cast and skips everything else). Give the species an
//      `inletOutlet` at that patch and the cast fails, so nothing accumulates
//      and there is no secondary emission either -- SEE lives in the wall BC
//      family.
//
//   3. THE FAR FIELD HAS NO SPACE CHARGE. Poisson in a dielectric region has no
//      rho on the right-hand side. The trick is therefore only valid while the
//      plasma stays AWAY from the interface: put the interface where the charge
//      density is negligible, and check it stays so.

// epsilonR    4.6;

// ************************************************************************* //
EOF
    echo "  $target  (dielectric: NEEDS A VALUE)"
    needsAttention="$needsAttention $target"
    return 0
}

# The plasma description is read PER REGION too, from constant/<region>/, and
# only the gas has one. Same stub treatment, same reasons: one home, no
# duplication, and no include-depth surgery.
#
# What this replaces in a hand-written Allrun is three separate copy-and-sed
# passes, each of which has silently failed before -- the include resolved to a
# path that did not exist, every $var was left UNDEFINED, changeDictionary wrote
# them verbatim into the field, and the solver died on "Illegal dictionary entry
# or environment variable name voltageRamp" while reading 0/gas/ePotential. The
# lesson recorded there is worth keeping: three earlier "fixes" corrected include
# PATHS in files that no longer had an include at all.
GAS_DICTS="plasmaSpeciesProperties plasmaTransportProperties photoionizationProperties"

echo "plasmaSetupRegions: numerics stubs"
for r in $regions; do
    mkdir -p "system/$r"
    for f in fvSchemes fvSolution; do
        cat > "system/$r/$f" <<EOF
/*--------------------------------*- C++ -*----------------------------------*\\
| GENERATED by tools/plasmaSetupRegions.sh -- DO NOT EDIT.                    |
|                                                                             |
| A stub. OpenFOAM requires a per-region $f and does not fall back to  |
| the case-wide one, so this file exists only to point at it. Edit             |
| system/$f instead; every region sees the change.                     |
|                                                                             |
| Regenerate with:  tools/plasmaSetupRegions.sh                               |
| (must be re-run after splitMeshRegions, which overwrites this with an        |
|  empty dummy)                                                               |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system/$r";
    object      $f;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// The FoamFile header above is NOT optional: without it the reader fails with
// "problem while reading header for object $f" -- foamDictionary is
// lenient about a headerless file, the solver is not.
#include "../$f"

// ************************************************************************* //
EOF
    done
    echo "  system/$r/{fvSchemes,fvSolution}  ->  #include ../{fvSchemes,fvSolution}"
done

echo "plasmaSetupRegions: $(echo $regions | wc -w) region(s) point at the case-wide numerics."
echo
echo "plasmaSetupRegions: per-region material properties"
for r in $regions; do
    k=$(kindOf "$r")
    if [ -z "$k" ]; then
        echo "  WARNING: region \`$r\` has no kind in constant/regionProperties?" >&2
        continue
    fi
    writeElectricalProperties "$r" "$k"
done

# ---- per-region plasma dictionaries (gas only) -----------------------------
echo
echo "plasmaSetupRegions: per-region plasma dictionaries"
for r in $regions; do
    [ "$(kindOf "$r")" = "gas" ] || continue
    for f in $GAS_DICTS; do
        [ -f "constant/$f" ] || continue
        if [ -f "constant/$r/$f" ] && ! grep -q "^#include \"\.\./$f\"" "constant/$r/$f"; then
            echo "  constant/$r/$f  (exists and is not a stub, left alone)"
            continue
        fi
        cat > "constant/$r/$f" <<EOF
/*--------------------------------*- C++ -*----------------------------------*\\
| GENERATED by tools/plasmaSetupRegions.sh -- DO NOT EDIT.                    |
|                                                                             |
| A stub. The plasma dictionaries are read PER REGION from constant/<region>/, |
| so this file must exist; it contains nothing of its own. Edit               |
| constant/$f instead.                                    |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/$r";
    object      $f;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// The FoamFile header above is NOT optional: without it the reader fails with
// "problem while reading header for object $f".
//
// A relative #include resolves against the directory of the file CONTAINING
// it, so "../$f" lands on constant/$f, and the
// "../configuration/config" inside THAT file still resolves from constant/.
// Nothing needs its include depth rewritten.
#include "../$f"

// ************************************************************************* //
EOF
        echo "  constant/$r/$f  ->  #include ../$f"
    done
done

if [ -n "$needsAttention" ]; then
    echo
    echo "==========================================================================="
    echo " ACTION REQUIRED before this case will run"
    echo "==========================================================================="
    for f in $needsAttention; do
        echo "   $f"
    done
    echo
    echo " Each has a commented-out \`epsilonR\`. Uncomment it and set your"
    echo " material's value; the file lists common barrier materials."
    echo
    echo " Not defaulted on purpose: for a dielectric, epsilonR = 1 is not a"
    echo " placeholder but VACUUM, and the case would run, converge and give a"
    echo " plausible WRONG answer with the barrier electrically invisible."
    echo "==========================================================================="
fi
