#!/usr/bin/env bash
# KLUSTR rates for the conditions of the ACDC output (Sect. 5 of the paper).
#   ./klustr_rates.sh              all acdc/acdc_*_*.csv
#   ./klustr_rates.sh file.csv ...
# Writes klustr_<free energies>_<set>.csv next to each input. Needs gfortran (or FC=...) and awk.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
FC=${FC:-gfortran}
R=1.987204259e-3                     # kcal mol-1 K-1

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
(cd "$work" && $FC -O2 -o klustr_rates "$here/../fortran/klustr.F90" "$here/klustr_rates.F90")

if [ $# -gt 0 ]; then files=("$@"); else files=("$here"/acdc/acdc_*_*.csv); fi
for csv in "${files[@]}"; do
  name=$(basename "$csv" .csv); name=${name#acdc_}
  fe=${name%_*}; set=${name##*_}
  case $set in
    4x4) N=4 ;;
    band8) N=8 ;;
    *) echo "unknown cluster set in $csv" >&2; exit 1 ;;
  esac

  # free energy of the pair 1A1D (kcal/mol): H and S, or G per temperature
  case $fe in
    olenius2013) awk '$1 == "1A1D" {print "HS", $2, $3}' "$here/acdc/SA-DMA_olenius2013_HS.txt" ;;
    *) for g in "$here/acdc/$fe"/G*K.txt; do
         T=${g##*/G}; T=${T%K.txt}
         awk -v T="$T" '$1 == "1A1D" {print "G", T, $2}' "$g"
       done ;;
  esac > "$work/pair.txt"

  # one condition per line: T, CS, CA, CB, dG_1/kB T
  awk -v R="$R" '
    NR == FNR { if ($1 == "HS") { hs = 1; H = $2; S = $3 } else G[$2] = $3; next }
    FNR == 1 { next }
    { split($0, f, ","); T = f[1]
      if (hs) g = (H - T*S*1e-3)/(R*T)
      else if (T in G) g = G[T]/(R*T)
      else { print "no free energy of 1A1D at " T " K" > "/dev/stderr"; exit 1 }
      printf "%s %s %s %s %.12e\n", T, f[2], f[3], f[4], g }' "$work/pair.txt" "$csv" > "$work/in.txt"

  n=$(wc -l < "$work/in.txt")
  "$work/klustr_rates" "$N" "$n" < "$work/in.txt" > "$work/j.txt"
  out=$(dirname "$csv")/klustr_${fe}_${set}.csv
  { echo "T,CS,CA,CB,converged,J_ACDC,J_KLUSTR"
    paste -d, <(tail -n +2 "$csv" | cut -d, -f1-6) "$work/j.txt"; } > "$out"
  echo "$out: $n conditions, N = $N"
done
