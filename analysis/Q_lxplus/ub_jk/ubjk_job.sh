#!/bin/bash
# Węzeł UBJK. Użycie: ubjk_job.sh prep | fit <porcja> | ref | finish
# Jackknife μ po plikach danych na ODSŁONIĘTEJ próbce (decyzja 05.11, flaga $UB_OUT/UNBLIND_OK + hash cfg). CHARM_UNBLIND=1 ustawiane tu.
# Scratch mktemp + trap rm -rf; wyjścia tylko EOS ($UJ_OUT); AFS tylko logi. DRYRUN=1 → tylko komendy (bez blokad i env).
set -eo pipefail
UJ=$(cd "$(dirname "$0")" && pwd); export UJ
source "$UJ/ubjk_common.sh"
step=$1; chunk=${2:--}
set -u
export CHARM_UNBLIND=1
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/ubjk_${step}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
PY="python3 -u $UJ/ubjk_run.py"
echo "[$(date -u +%FT%TZ)] UBJK $step $chunk host $(hostname) WD=$WD UB_OUT=$UB_OUT UJ_OUT=$UJ_OUT NCHUNK=$UJ_NCHUNK MODES=$UJ_MODES"
if [ "${DRYRUN:-0}" != 1 ]; then
  uj_gate || exit 3
  uj_env
  mkdir -p "$UJ_OUT" 2>/dev/null || true
  [ "$step" = prep ] || [ -s "$UJ_OUT/perfile.json" ] || { echo "STOP: brak $UJ_OUT/perfile.json (prep)"; exit 4; }
fi
case $step in
  prep)
    legs_arg=(--legs "$UJ_LEGS")
    if [ "${UJ_REBUILD:-0}" = 1 ] || [ ! -s "$UJ_LEGS" ]; then
      legs_arg=(--legs "" --rebuild-dir "$WD/rebuild" --norm "$NORM")
    fi
    # shellcheck disable=SC2086
    run $PY prep --cfg "$UJ_CFG" --pairs-dir "$UJ_PAIRS" "${legs_arg[@]}" --files-json "$UJ_FILES_JSON" --groups "$UJ_GROUPS" --out "$WD"
    if [ "${DRYRUN:-0}" != 1 ]; then
      uj_put "$WD/perfile.npz" "$UJ_OUT/perfile.npz"; uj_put "$WD/perfile.json" "$UJ_OUT/perfile.json"
      if [ -s "$WD/rebuild/pairs/legs_data.parquet" ]; then   # przebudowane nogi → EOS (ref ich potrzebuje); ścieżka w perfile.json poprawiona
        uj_put "$WD/rebuild/pairs/legs_data.parquet" "$UJ_OUT/legs_data_rebuilt.parquet"
        python3 - "$UJ_OUT/perfile.json" "$UJ_OUT/legs_data_rebuilt.parquet" <<'PY'
import json, sys
p, legs = sys.argv[1], sys.argv[2]
j = json.load(open(p)); j["legs_src"] = legs; json.dump(j, open(p, "w"), indent=1)
PY
      fi
      [ -f "$WD/ubjk.log" ] && uj_put "$WD/ubjk.log" "$UJ_OUT/logs/prep.log"
    fi ;;
  fit)
    case $chunk in ''|*[!0-9]*) echo "fit: porcja musi być liczbą"; exit 2 ;; esac
    # shellcheck disable=SC2086
    run $PY fit --pairs-dir "$UJ_PAIRS" --perfile "$UJ_OUT" --chunk "$chunk" --nchunk "$UJ_NCHUNK" --modes "$UJ_MODES" --out "$WD"
    if [ "${DRYRUN:-0}" != 1 ]; then
      uj_put "$WD/fit_${chunk}.json" "$UJ_OUT/fit_${chunk}.json"
      [ -f "$WD/ubjk.log" ] && uj_put "$WD/ubjk.log" "$UJ_OUT/logs/fit_${chunk}.log"
    fi ;;
  ref)
    # shellcheck disable=SC2086
    run $PY ref --pairs-dir "$UJ_PAIRS" --perfile "$UJ_OUT" --out "$WD"
    [ "${DRYRUN:-0}" = 1 ] || uj_put "$WD/ref_unblinded_jackknife.json" "$UJ_OUT/ref_unblinded_jackknife.json" ;;
  finish)
    # shellcheck disable=SC2086
    run $PY finish --pairs-dir "$UJ_PAIRS" --in "$UJ_OUT" --out "$WD"
    if [ "${DRYRUN:-0}" != 1 ]; then
      uj_put "$WD/ubjk.json" "$UJ_OUT/ubjk.json"; uj_put "$WD/ubjk.md" "$UJ_OUT/ubjk.md"
      date -u +%FT%TZ > "$WD/DONE"; uj_put "$WD/DONE" "$UJ_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK UBJK $step $chunk"
