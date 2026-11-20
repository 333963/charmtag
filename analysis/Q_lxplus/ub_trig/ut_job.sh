#!/bin/bash
# Węzeł UBTRIG. ut_job.sh trg <g>  → $UT_OUT/trg/trg_<g>.npz (pliki ROOT 428770 o indeksie ≡ g mod UT_NG)
#                ut_job.sh fit      → $UT_OUT/{ubtrig.json,UBTRIG.md,DONE}
set -eo pipefail
UT=$(cd "$(dirname "$0")" && pwd); export UT
source "$UT/ut_common.sh"
step=$1; g=${2:--}
set -u
export CHARM_UNBLIND=1
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/ut_${step}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
echo "[$(date -u +%FT%TZ)] UBTRIG $step $g host $(hostname) UT_OUT=$UT_OUT"
uj_gate || exit 3
uj_env
export PYTHONPATH="$UT:$Q/ub_jk:${PYTHONPATH:-}"
mkdir -p "$UT_OUT/trg"
case $step in
  trg)
    mapfile -t ALL < <(ut_files)
    FILES=(); for i in "${!ALL[@]}"; do [ $((i % UT_NG)) -eq "$g" ] && FILES+=("${ALL[$i]}"); done
    [ ${#FILES[@]} -gt 0 ] || { echo "pusta grupa $g"; exit 2; }
    python3 -u "$UT/trg_root_extract.py" --out "$WD/trg_$g.npz" "${FILES[@]}"
    ut_put "$WD/trg_$g.npz" "$UT_OUT/trg/trg_$g.npz" ;;
  fit)
    n=$(ls "$UT_OUT"/trg/trg_*.npz | wc -l); [ "$n" -eq "$UT_NG" ] || { echo "STOP: trg $n / $UT_NG"; exit 4; }
    cp "$UT_OUT"/trg/trg_*.npz "$WD/"; cp "$BTD/bt_data_428770.npz" "$BTD/bt_fd_nd.npz" "$WD/"
    python3 -u "$UT/ubtrig_run.py" --pairs-dir "$UJ_PAIRS" --legs "$UJ_LEGS" --trg "$WD/trg_*.npz" --bt "$WD/bt_data_428770.npz" \
      --bt-mc "$WD/bt_fd_nd.npz" --out "$WD/res" 2>&1 | tee "$WD/ubtrig.log"
    for f in ubtrig.json UBTRIG.md; do ut_put "$WD/res/$f" "$UT_OUT/$f"; done
    ut_put "$WD/ubtrig.log" "$UT_OUT/ubtrig.log"
    date -u +%FT%TZ > "$WD/DONE"; ut_put "$WD/DONE" "$UT_OUT/DONE" ;;
  null)   # test zerowy: klasy losowe (ziarno g) → $UT_OUT/null/s<g>/
    cp "$UT_OUT"/trg/trg_*.npz "$WD/"; cp "$BTD/bt_data_428770.npz" "$BTD/bt_fd_nd.npz" "$WD/"
    python3 -u "$UT/ubtrig_run.py" --pairs-dir "$UJ_PAIRS" --legs "$UJ_LEGS" --trg "$WD/trg_*.npz" --bt "$WD/bt_data_428770.npz" \
      --bt-mc "$WD/bt_fd_nd.npz" --null "$g" --out "$WD/res" 2>&1 | tee "$WD/ubtrig.log"
    mkdir -p "$UT_OUT/null/s$g"
    for f in ubtrig.json UBTRIG.md; do ut_put "$WD/res/$f" "$UT_OUT/null/s$g/$f"; done ;;
  nullagg)
    n=$(ls "$UT_OUT"/null/s*/ubtrig.json 2>/dev/null | wc -l); [ "$n" -ge 10 ] || { echo "STOP: tylko $n ziaren"; exit 4; }
    [ -s "$UT_OUT/ubtrig.json" ] || { echo "STOP: brak prawdziwego ubtrig.json"; exit 4; }
    python3 -u "$UT/ubtrig_null_agg.py" --real "$UT_OUT/ubtrig.json" --null "$UT_OUT/null/s*/ubtrig.json" --out "$WD/agg"
    for f in ubtrig_null.json UBTRIG_NULL.md; do ut_put "$WD/agg/$f" "$UT_OUT/null/$f"; done
    date -u +%FT%TZ > "$WD/NULL_DONE"; ut_put "$WD/NULL_DONE" "$UT_OUT/null/NULL_DONE" ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK UBTRIG $step $g"
