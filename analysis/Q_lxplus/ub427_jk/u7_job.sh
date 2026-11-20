#!/bin/bash
# Węzeł UB427JK. u7_job.sh prep | fit <porcja> | finish   (ubjk_run.py z ../ub_jk, tryb fixed, grupy events:$U7_G)
# Scratch mktemp + trap; wyjścia tylko EOS ($U7_OUT); AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
U7=$(cd "$(dirname "$0")" && pwd); export U7
source "$U7/u7_common.sh"
step=$1; chunk=${2:--}
set -u
export CHARM_UNBLIND=1
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/u7_${step}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
put() { [ "${DRYRUN:-0}" = 1 ] || put_output "$1" "$2"; }
PY="python3 -u $UJ_RUN"
echo "[$(date -u +%FT%TZ)] UB427JK $step $chunk host $(hostname) WD=$WD U7_OUT=$U7_OUT G=$U7_G NCHUNK=$U7_NCHUNK"
if [ "${DRYRUN:-0}" != 1 ]; then
  u7_gate || exit 3
  u7_env
  [ "$step" = prep ] || [ -s "$U7_OUT/perfile.json" ] || { echo "STOP: brak $U7_OUT/perfile.json (prep)"; exit 4; }
fi
case $step in
  prep)
    # cfg i wejścia identyczne z węzłem `pairs` UB427 (u4_job.sh), out_dir nadpisuje ubjk (--rebuild-dir)
    run python3 "$Q/pf_pairs/tools/merge_cfg.py" "$WD/cfg427.yaml" "$U4_NOM_CFG" -s "paths.out_dir=$WD/run" -s "paths.scores_dir=$WD/scores" \
      -s "paths.samples.data_none=$U4_OUT/h5/data_427929_*_mb.h5" -s "paths.sample_info={}" -s "pairs.mb_rates_json=$WD/mb_rates_427929.json"
    if [ "${DRYRUN:-0}" != 1 ]; then
      cp "$U4_OUT/cw/mb_rates_427929.json" "$WD/"
      for sp in D0 Dplus; do mkdir -p "$WD/scores/${sp}_H"
        cp "$NOMINAL/scores/${sp}_H/signal_mc.h5" "$NOMINAL/scores/${sp}_H/incl_mc.h5" "$WD/scores/${sp}_H/"
        cp "$U4_OUT/scores/${sp}_H/data_none.h5" "$WD/scores/${sp}_H/"; done
      mkdir -p "$WD/pin"; cp "$U7_PAIRS"/pair_inputs.json "$U7_PAIRS"/pair_inputs.npz "$U7_PAIRS"/pairs_full.npz "$U7_PAIRS"/fit2d_unblind.json "$WD/pin/"
    fi
    LEGS=$U7_OUT/legs_data_rebuilt.parquet
    if [ "${DRYRUN:-0}" != 1 ] && [ -s "$LEGS" ] && [ "${U7_REBUILD:-0}" != 1 ]; then legs_arg=(--legs "$LEGS"); else legs_arg=(--legs "" --rebuild-dir "$WD/rebuild" --norm "$NORM"); fi
    run $PY prep --cfg "$WD/cfg427.yaml" --pairs-dir "$WD/pin" "${legs_arg[@]}" --groups "events:$U7_G" --out "$WD"
    if [ "${DRYRUN:-0}" != 1 ]; then
      [ -s "$WD/rebuild/pairs/legs_data.parquet" ] && put "$WD/rebuild/pairs/legs_data.parquet" "$LEGS"
      put "$WD/perfile.npz" "$U7_OUT/perfile.npz"; put "$WD/perfile.json" "$U7_OUT/perfile.json"
      [ -f "$WD/ubjk.log" ] && put "$WD/ubjk.log" "$U7_OUT/logs/prep.log"
    fi ;;
  fit)
    case $chunk in ''|*[!0-9]*) echo "fit: porcja musi być liczbą"; exit 2 ;; esac
    [ "${DRYRUN:-0}" = 1 ] || { mkdir -p "$WD/pin"; cp "$U7_PAIRS"/pair_inputs.json "$U7_PAIRS"/pair_inputs.npz "$U7_PAIRS"/pairs_full.npz "$U7_PAIRS"/fit2d_unblind.json "$WD/pin/"; }
    run $PY fit --pairs-dir "$WD/pin" --perfile "$U7_OUT" --chunk "$chunk" --nchunk "$U7_NCHUNK" --modes "$U7_MODES" --out "$WD"
    if [ "${DRYRUN:-0}" != 1 ]; then put "$WD/fit_${chunk}.json" "$U7_OUT/fit_${chunk}.json"; [ -f "$WD/ubjk.log" ] && put "$WD/ubjk.log" "$U7_OUT/logs/fit_${chunk}.log"; fi ;;
  finish)
    [ "${DRYRUN:-0}" = 1 ] || { mkdir -p "$WD/pin"; cp "$U7_PAIRS"/pair_inputs.json "$U7_PAIRS"/pair_inputs.npz "$U7_PAIRS"/pairs_full.npz "$U7_PAIRS"/fit2d_unblind.json "$WD/pin/"; }
    run $PY finish --pairs-dir "$WD/pin" --in "$U7_OUT" --out "$WD"
    if [ "${DRYRUN:-0}" != 1 ]; then
      put "$WD/ubjk.json" "$U7_OUT/ubjk.json"; put "$WD/ubjk.md" "$U7_OUT/ubjk.md"
      date -u +%FT%TZ > "$WD/DONE"; put "$WD/DONE" "$U7_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK UB427JK $step $chunk"
