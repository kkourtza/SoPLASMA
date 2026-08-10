#!/usr/bin/env bash
# Cheap health check on a streamer run. Reports only what would change a
# decision: did the mechanism load, which species are missing, are the fields
# finite and physically sized, and how fast is it going.
CASE=${1:-$HOME/streamer-case}
L=$CASE/log.soPlasmaFoam
[ -f "$L" ] || { echo "no solver log yet at $L"; exit 1; }

echo "== mechanism =="
grep -E "plasmaReactionRates:|plasmaTransport: sourceModel" "$L" | head -4

echo
echo "== untransported products (the to-do list for plasmaSpeciesProperties) =="
grep -oE "produces '[^']+'" "$L" | sort -u | sed "s/produces //" | tr '\n' ' '; echo
echo "  distinct: $(grep -oE "produces '[^']+'" "$L" | sort -u | wc -l)"

echo
echo "== progress =="
echo "  timesteps completed : $(grep -c '^Time = ' "$L")"
echo "  last time           : $(grep '^Time = ' "$L" | tail -1)"
echo "  ExecutionTime       : $(grep 'ExecutionTime' "$L" | tail -1)"

echo
echo "== field ranges (last occurrence) =="
for f in n_e alpha k_eff reducedE; do
  line=$(grep -iE "min\(.*$f.*\)|$f.*min" "$L" | tail -1)
  [ -n "$line" ] && echo "  $line"
done
grep -E "^(Min|Max|min|max).*(n_e|alpha|k_eff)" "$L" | tail -6

echo
echo "== trouble =="
grep -icE "FOAM FATAL|floating point|nan|inf\b" "$L" | sed 's/^/  suspicious lines: /'
grep -E "FOAM FATAL" -A3 "$L" | head -6
