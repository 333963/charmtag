#!/bin/bash
# Test lokalny gen_norm end-to-end (bez lxplus/EOS/condora). Użycie: tests/run_local.sh [KATALOG_ROBOCZY] [N_PLIKÓW_H5=3]
#  1) atrapy ROOT z prawdziwych lokalnych części h5 v1.2 mb_nd (FEAT_charmnet/data/mb_nd) + atrapy danych 428770 + atrapy scores (apply, losowy score);
#  2) ORYGINALNE skrypty etap1 (mb_truth.py, data_events.py, w2d.py, analyze2.py 2d; kopie ze zmienionymi ścieżkami) na tych atrapach → referencja;
#  3) łańcuch węzłów DAG (make_dag_gn.sh + gn_job.sh root/legs/combine) z lokalnym common.sh (RAW/EOS_OUT w katalogu roboczym);
#  4) porównanie: mb_rates_2d odtworzone vs etap1 (bit w bit), w2d.npy identyczne, test akceptacji z --ref = wynik etap1 → PASS;
#  5) tests/test_closure.py (zabawka z ε faktoryzowanym: truth/direct ≈ 1).
set -euo pipefail
G=$(cd "$(dirname "$0")/.." && pwd); Q=$(cd "$G/.." && pwd); W=$(cd "$Q/.." && pwd); CH=$(cd "$W/.." && pwd)
S=${1:-$(mktemp -d)}; N=${2:-3}; mkdir -p "$S"; S=$(cd "$S" && pwd)
PY=${PY:-$CH/.venv_h5/bin/python}; E1=$CH/etap1_definicje_masy_v20260925_01
export PYTHONDONTWRITEBYTECODE=1
H5=$(ls "$W"/FEAT_charmnet/data/mb_nd/mb_nd_*_none.h5 | head -"$N")
echo "== 1) atrapy w $S"
"$PY" "$G/tests/make_fake.py" $H5 --out "$S" | tail -2
PYTHONPATH=$G/src "$PY" "$G/tests/fake_scores.py" --dir "$S" | tail -1
echo "== 2) oryginalny etap1"
O=$S/etap1_orig; mkdir -p "$O/out2"; cp "$E1"/code/{mb_truth,data_events,w2d,analyze2}.py "$E1/weights_ntrk.csv" "$O/"
sed -i "s|/home/[^/]*/work/mbpairs/out2/|$O/out2/|; s|/home/[^/]*/work/etap1/weights_ntrk.csv|$O/weights_ntrk.csv|" "$O/w2d.py"
sed -i "s|/home/[^/]*/work/etap1/weights_ntrk.csv|$O/weights_ntrk.csv|; s|/home/[^/]*/work/trig/|$O/|g" "$O/analyze2.py"
( cd "$O"; for f in "$S"/root/mb_nd/*.root; do "$PY" mb_truth.py "$f" "out2/$(basename "$f" .root)" >/dev/null 2>&1; done
  for t in 000021 000030; do "$PY" data_events.py "$S"/root/data_428770/*_$t.*.root "ev_data_$t.parquet" >/dev/null; done
  "$PY" w2d.py >/dev/null; "$PY" analyze2.py 2d >/dev/null; echo "etap1: $(ls out2 | wc -l) plików, mb_rates_2d.json OK" )
echo "== 3) węzły DAG lokalnie"
MB=$(awk '$1=="mb_nd"{print $2}' "$Q/config/samples.txt"); DD=$(awk '$1=="data_428770"{print $2}' "$Q/config/samples.txt")
mkdir -p "$S/raw/$MB" "$S/raw/$DD"; ln -sf "$S"/root/mb_nd/*.root "$S/raw/$MB/"; ln -sf "$S"/root/data_428770/*.root "$S/raw/$DD/"
export GN_FMT=v11 Q RAW=$S/raw EOS_OUT=$S/eosout LCG_SETUP=/dev/null VENV=$CH/.venv_h5 BASE_CFG=$S/cfg.yaml NOMINAL=$S/run REF_RATES=
TMPG=$S/pkg; rm -rf "$TMPG"; mkdir -p "$TMPG"; cp -a "$G" "$TMPG/gen_norm"; cp "$Q/common.sh" "$TMPG/"; mkdir -p "$TMPG/config"; cp "$Q/config/samples.txt" "$TMPG/config/"
export Q=$TMPG; GG=$TMPG/gen_norm
"$GG/make_dag_gn.sh"
grep '^VARS' "$GG/gn.dag" | while read -r _ node rest; do
  v() { echo "$rest" | grep -o "$1=\"[^\"]*\"" | cut -d'"' -f2; }
  "$GG/gn_job.sh" "$(v step)" "$(v kind)" "$(v tag)" "$(v src)" > "$S/log_$node.txt" 2>&1 || { echo "BŁĄD węzła $node:"; tail -20 "$S/log_$node.txt"; exit 1; }
  echo "  $node: $(tail -1 "$S/log_$node.txt")"
done
OUT=$EOS_OUT/charm_q/gen_norm
echo "== 4) porównanie z etap1"
"$PY" "$GG/gn_combine.py" --mc "$OUT/parts/mc" --data "$OUT/parts/data" --legs "$OUT/legs" --out "$S/out_vs_etap1" \
      --ref "$O/mb_rates_2d.json" --ref-counts none 2>&1 | grep -E "akceptacja" | cut -c1-160
"$PY" - "$O" "$S/out_vs_etap1" <<'PY'
import json, sys, numpy as np
o, n = sys.argv[1:]
def fl(d, p=""):
    for k, v in d.items():
        yield from (fl(v, p + k + ".") if isinstance(v, dict) else [(p + k, v)])
A = dict(fl(json.load(open(f"{o}/mb_rates_2d.json")))); B = dict(fl(json.load(open(f"{n}/mb_rates_2d_repro.json"))))
bad = [k for k in A if A[k] != B.get(k)]
same_w = np.array_equal(np.load(f"{o}/w2d.npy"), np.load(f"{n}/w2d.npy"))
acc = json.load(open(f"{n}/acceptance.json"))["verdict"]
print(f"mb_rates_2d: {len(A)} pól, różnych {len(bad)} {bad[:5]}; w2d.npy identyczne: {same_w}; akceptacja vs etap1: {acc}")
sys.exit(0 if (not bad and same_w and acc == "PASS") else 1)
PY
echo "== 5) test domknięcia (zabawka)"
"$PY" "$G/tests/test_closure.py" | tail -6
echo "TEST LOKALNY OK ($S)"
