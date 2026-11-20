#!/bin/bash
# Test lokalny MULT_charm (bez sieci): MC = prawdziwe pliki v1.1 mb_nd (tryb compat), „dane” = atrapy v2.0 z innych plików v1.1
# (make_synth_v20.py: emulacja triggerów/PS/GRL). Zamknięcie: dane ≡ MC → R≈1, Δα≈0; wstrzyknięcie --inject δ → Δα≈δ.
# Użycie: tests/run_local.sh [N_MC=16] [N_DATA=8] [J=4]
set -euo pipefail
cd "$(dirname "$0")/.."
P=${P:-$HOME/Documents/charm/.venv_h5/bin/python}
NMC=${1:-16}; NDA=${2:-8}; J=${3:-4}
SRC=$HOME/Documents/charm/mc23_minbiasND_Dpt3to150GeV_v9
T=tests/local; mkdir -p $T/mc $T/data $T/synth $T/logs
mapfile -t L < <(ls $SRC/*.root | sort)
jobs_mc() { for f in "${L[@]:0:$NMC}"; do t=$(basename $f | grep -o '\._[0-9]\{6\}\.' | tr -d '._'); echo "$P mult_extract.py --kind mc $f $T/mc/mb_nd_$t.npz --tag $t > $T/logs/mc_$t.log 2>&1"; done; }
jobs_da() { for f in "${L[@]:$NMC:$NDA}"; do t=$(basename $f | grep -o '\._[0-9]\{6\}\.' | tr -d '._')
  echo "[ -s $T/synth/data_$t.root ] || $P tests/make_synth_v20.py $f $T/synth/data_$t.root --data --max-events 1000000 --seed $((10#$t)) > $T/logs/syn_$t.log 2>&1; $P mult_extract.py --kind data $T/synth/data_$t.root $T/data/data_$t.npz --tag $t > $T/logs/da_$t.log 2>&1"; done; }
{ jobs_mc; jobs_da; } | while read -r c; do o=$(echo "$c" | grep -o "$T/[a-z]*/[a-z_]*[0-9]*\.npz" | tail -1); [ -s "$o" ] && continue; echo "$c"; done | xargs -P "$J" -I{} bash -c "{}"
for s in mbts sptrk2; do for e in bl blx opp trk500; do
  $P mult_analyze.py --data $T/data --mc $T/mc --scheme $s --est $e --out $T/res > $T/logs/an_${s}_$e.log 2>&1 || echo "BŁĄD $s $e"
done; done
$P mult_analyze.py --data $T/data --mc $T/mc --scheme mbts --est blx --inject 0.3 --label inject0.3_mbts_blx --out $T/res > $T/logs/an_inject.log 2>&1
$P mult_analyze.py --data $T/data --mc $T/mc --scheme mbts --est blx --w2 self --out $T/res > $T/logs/an_w2self.log 2>&1
$P mult_analyze.py --data $T/data --mc $T/mc --pairs-data --out $T/res > $T/logs/an_g5.log 2>&1 && echo "BŁĄD: G5 nie zatrzymał" || echo "G5 OK: $(cat $T/logs/an_g5.log)"
echo KONIEC
