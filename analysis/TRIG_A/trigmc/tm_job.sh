#!/bin/bash
# Węzeł DAG trigmc.  tm_job.sh gen <procid> | combine
# Scratch w katalogu condora; wyjścia tylko na EOS; na AFS tylko logi.
set -eo pipefail
TM=$(cd "$(dirname "$0")" && pwd); export TM
source "$TM/tm_common.sh"
step=$1; arg=${2:-0}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/tm_${step}_XXXXXX")
trap 'cd /; rm -rf "$WD"' EXIT
cd "$WD"
echo "[$(date -u +%FT%TZ)] TM $step $arg host $(hostname) WD=$WD TM_OUT=$TM_OUT"
q_env
case $step in
  gen)
    seed=$(( TM_SEED0 + arg + 1 )); out=$TM_OUT/gen/gen_$seed.tsv.gz
    if eos_exists "$out.json"; then echo "pomijam (istnieje) $out"; exit 0; fi
    "$TM/bin/pythia_dd_jet" "$TM_NEV" "$seed" "gen_$seed.tsv" A3 13600 nd > gen.log 2>&1 || { tail -40 gen.log; exit 1; }
    gzip -f "gen_$seed.tsv"
    put_output "gen_$seed.tsv.h3em" "$TM_OUT/gen/gen_$seed.tsv.h3em"
    put_output "gen_$seed.tsv.h3all" "$TM_OUT/gen/gen_$seed.tsv.h3all"
    put_output "gen_$seed.tsv.gz" "$out"
    put_output "gen_$seed.tsv.json" "$out.json"          # .json na końcu = znacznik kompletności
    ;;
  combine)
    python3 -u "$TM/tm_combine.py" --gen "$TM_OUT/gen" --out "$WD/res"
    for f in "$WD"/res/*; do put_output "$f" "$TM_OUT/res/$(basename "$f")"; done
    date -u +%FT%TZ > DONE; put_output DONE "$TM_OUT/DONE"
    ;;
  calib)
    python3 -u "$TM/tm_combine.py" --gen "$TM_OUT/gen" --out "$WD/res" --calib "$TM/j12eff.json"
    for f in "$WD"/res/*; do put_output "$f" "$TM_OUT/res_cal/$(basename "$f")"; done
    date -u +%FT%TZ > DONE; put_output DONE "$TM_OUT/res_cal/DONE"
    ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK TM $step $arg"
