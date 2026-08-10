#!/usr/bin/env bash
# Build the streamer case in a persistent location, with the mechanism wired in.
#   ./setup-case.sh <endTime> [sourceModel] [rampVoltsPerStep]
#
# rampVoltsPerStep: 0 (default) keeps the tutorial's hard 0->18.75 kV step, so
# the baseline stays comparable with the existing model. Non-zero ramps the
# electrode linearly at that rate.
#
# NOTE the rise-rate limiter in system/plasmaSimulationControls cannot do this.
# It only shrinks deltaT, and the Poisson equation is elliptic -- the potential
# is re-solved from scratch every step, so a smaller timestep does not reduce
# the spatial gradient at the electrode. A real ramp must be a time-dependent
# boundary condition.
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1 || true
set -e
export SoPLASMA=$HOME/soplasma-scratch SoPLASMA_SRC=$HOME/soplasma-scratch/src SoPLASMA_ETC=$HOME/soplasma-scratch/etc
END=${1:-1e-12}; MODEL=${2:-reactions}; RAMP=${3:-0}
B=$HOME/Projects/BoltzmannSolver
CASE=$HOME/streamer-case

rm -rf "$CASE"
cp -r "$SoPLASMA/tutorials/plasma/soPlasmaFoam/positiveStreamer/positiveStreamer_fixedMesh" "$CASE"
cd "$CASE"

mkdir -p constant/plasmaTables
cp "$B/build/mech/air_plasma.foam" constant/
cp "$B"/build/tables/{k_,alpha_,eta_}*_vs_reducedE constant/plasmaTables/ 2>/dev/null || true

cat > constant/plasmaTransportProperties <<EOF
FoamFile { version 2.0; format ascii; class dictionary; object plasmaTransportProperties; }

chemistry
{
    mechanism       "constant/air_plasma.foam";
    tableDir        "constant/plasmaTables";
    lookupVariable  reducedE;
    tableKey        reducedE;
    sourceModel     $MODEL;
    outOfBounds     extrapolate;
    eedfRefresh     never;
}
EOF

sed -i "s|^endTime .*|endTime                             $END;|" configuration/config
# Write often enough that a short run actually produces output: with the
# tutorial's 4e-9 interval a 2e-11 run writes nothing at all.
sed -i "s|^writeInterval .*|writeInterval                       $(~/ct-env/bin/python -c "print(f'{float('$END')/10:.6e}')");|" configuration/config
sed -i 's|^python initGaussianSeed.py|~/ct-env/bin/python initGaussianSeed.py|' Allrun-serial

# ---------------------------------------------------------------------------
# Replace the lumped `pIon` with the ions the mechanism actually produces.
#
# The tutorial carries one generic positive ion, which cannot represent a
# mechanism whose ionisation channels make N2+ and O2+ and whose attachment
# channels make O- and O2-. With only `pIon` present every charged product is
# discarded, and the run creates electrons with no counter-charge -- measured
# at +1.45e18 electrons against +0 ions per step before this change. The
# solver now refuses to start in that state, so the case has to carry them.
#
# All four are `immobile`, keeping the tutorial's approximation: over the ~1 ns
# of interest an ion drifts ~5 um against ~100 um for an electron, so ion
# transport is negligible while ion SPACE CHARGE -- which immobile ions still
# carry -- is what actually matters for the field.
# ---------------------------------------------------------------------------
~/ct-env/bin/python - <<'PYEOF'
import pathlib, re

u = 1.66053907e-27
IONS = [                      # name, charge, mass [kg]
    ("N2p",  +1, 28*u),
    ("O2p",  +1, 32*u),
    ("Om",   -1, 16*u),
    ("O2m",  -1, 32*u),
]

# --- constant/plasmaSpeciesProperties -------------------------------------
p = pathlib.Path('constant/plasmaSpeciesProperties'); s = p.read_text()
s = s.replace('activeSpecies (e pIon);',
              'activeSpecies (e ' + ' '.join(n for n, _, _ in IONS) + ');', 1)

blocks = "\n".join(
    f"""  {n}
  {{
    transportModel            immobile;
    charge                    {q:+d};
    mass                      {m:.4e};
    minNumberDensity          1e5;
  }}
""" for n, q, m in IONS)

old = re.search(r'\n  pIon\n  \{.*?\n  \}\n', s, re.S)
assert old, "pIon species block not found"
s = s[:old.start()] + "\n" + blocks + s[old.end():]
p.write_text(s)

# --- etc/changeDictionary --------------------------------------------------
# Initial densities: the seed is placed entirely in N2+, the dominant ion of an
# air discharge, and the other three sit at the floor. The split among positive
# ions is arbitrary at t=0 and is overwritten by the mechanism within a few
# hundred ps; what must be right is the TOTAL positive charge, which has to
# match the electron seed exactly.
p = pathlib.Path('etc/changeDictionary'); s = p.read_text()
m = re.search(r'\nn_pIon\n\{.*?\n\}\n(?=\n|\Z)', s, re.S)
assert m, "n_pIon changeDictionary block not found"
tmpl = m.group(0)
out = []
for n, _, _ in IONS:
    b = tmpl.replace('n_pIon', 'n_' + n, 1)
    if n != 'N2p':
        b = b.replace('internalField   uniform 1e13;',
                      'internalField   uniform 1e5;', 1)
    out.append(b)
s = s[:m.start()] + "".join(out) + s[m.end():]
p.write_text(s)

# --- system/fvSolution{,-foam,-petsc}, system/fvSchemes --------------------
for fn in ('system/fvSolution', 'system/fvSolution-foam', 'system/fvSolution-petsc'):
    p = pathlib.Path(fn)
    if not p.exists(): continue
    s = p.read_text()
    for pat, repl in (
        # solver block, and its ...Final twin
        (r'\n    n_pIon\n    \{.*?\n    \}\n    n_pIonFinal\n    \{.*?\n    \}\n',
         lambda t: "".join(t.replace('n_pIon', 'n_'+n) for n, _, _ in IONS)),
        # PIMPLE residualControl entry
        (r'\n        n_pIon\n        \{.*?\n        \}\n',
         lambda t: "".join(t.replace('n_pIon', 'n_'+n) for n, _, _ in IONS)),
    ):
        mm = re.search(pat, s, re.S)
        if mm:
            s = s[:mm.start()] + repl(mm.group(0)) + s[mm.end():]
    p.write_text(s)

p = pathlib.Path('system/fvSchemes'); s = p.read_text()
s = s.replace('    n_pIon          ;\n',
              "".join(f'    n_{n:<15s};\n' for n, _, _ in IONS), 1)
p.write_text(s)

# --- initGaussianSeed.py ---------------------------------------------------
# Two changes. The field is renamed, and the SAME profile is now written to the
# electrons as well.
#
# The tutorial seeded the Gaussian into the ion field alone and left n_e
# uniform, which is a net positive space charge of 5e18 e/m^3 in the blob
# rather than the neutral seed of the benchmark this case follows (Bagheri et
# al., Plasma Sources Sci. Technol. 27 (2018) 095002, case 1). A non-neutral
# seed launches its own field transient that has nothing to do with the applied
# voltage, so the seed is made neutral here.
p = pathlib.Path('initGaussianSeed.py'); s = p.read_text()
assert "'0', 'n_pIon'" in s, "seed target not found"
s = s.replace("field_path = os.path.join(sol_path, '0', 'n_pIon')",
              "field_paths = [os.path.join(sol_path, '0', f)\n"
              "               for f in ('n_N2p', 'n_e')]", 1)
s = s.replace("with open(field_path, 'r') as f:\n    lines = f.readlines()",
              "for field_path in field_paths:\n"
              "  with open(field_path, 'r') as f:\n    lines = f.readlines()", 1)
# n_e's `far` patch is inletOutlet with `value $internalField`. Harmless while
# the internal field is uniform; once it is a 1.1M-entry list the macro expands
# to a second copy of that list and OpenFOAM refuses to transfer the same
# compound token twice. The patch value is only an initial guess for an
# outflow boundary, so the background density is the right thing to pin it to.
s = s.replace("# 5. Write back to 0/n_pIon",
              "new_lines = [l.replace('$internalField', f'uniform {bg}')\n"
              "             for l in new_lines]\n\n"
              "# 5. Write back", 1)

# indent the reconstruct + write body into the new loop
head, sep, tail = s.partition("# 4. Filter and Reconstruct the file")
tail = "".join(('  ' + ln if ln.strip() else ln) for ln in tail.splitlines(True))
s = head + sep + tail
pathlib.Path('initGaussianSeed.py').write_text(s)

print("  species: e + " + ", ".join(n for n, _, _ in IONS) + " (ions immobile)")
print("  seed:    neutral -- Gaussian written to n_e and n_N2p")
PYEOF

if [ "$RAMP" != "0" ]; then
    ~/ct-env/bin/python - "$RAMP" <<'PYEOF'
import re, sys, pathlib
rate = float(sys.argv[1])
cfg = pathlib.Path('configuration/config').read_text()
dt = float(re.search(r'^deltaT\s+([0-9.eE+-]+)', cfg, re.M).group(1))
v  = float(re.search(r'^appliedVoltage\s+([0-9.eE+-]+)', cfg, re.M).group(1))
tr = v/rate*dt
p = pathlib.Path('etc/changeDictionary'); s = p.read_text()
old = """        active_electrode
        {
            type            fixedValue;
            value           uniform $appliedVoltage;
        }"""
new = ("        active_electrode\n        {\n"
       "            type            uniformFixedValue;\n"
       f"            uniformValue    table ((0 0) ({tr:.6e} {v:g}));\n"
       "            value           uniform 0;\n        }")
assert old in s, "electrode BC anchor not found"
p.write_text(s.replace(old, new, 1))
print(f"  electrode ramp: 0 -> {v:g} V over {tr:.4e} s  ({rate:g} V per {dt:g} s step)")
PYEOF
fi
echo "case ready at $CASE  (endTime=$END, sourceModel=$MODEL, ramp=$RAMP V/step)"
