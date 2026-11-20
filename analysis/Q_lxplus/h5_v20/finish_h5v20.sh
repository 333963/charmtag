#!/bin/bash
# finish_h5v20.sh — węzeł FINAL DAG (też ręcznie): meta/<próbka>.json + meta/summary.json; parts/DONE tylko gdy komplet wszystkich części.
set -uo pipefail
HV=$(cd "$(dirname "$0")" && pwd); export HV
source "$HV/h5v20_common.sh"
hv_lists || echo "UWAGA: hv_lists zgłosił problem (brak DONE wejścia?)"
export PYTHONDONTWRITEBYTECODE=1
( q_env >/dev/null 2>&1; python3 "$HV/meta_h5v20.py" "$HV/lists/jobs.txt" "$PARTS" "$H5V20/meta" ) || { echo "BŁĄD meta"; exit 1; }
nmiss=$(wc -l < "$HV/lists/missing.txt")
if [ "$nmiss" -eq 0 ] && [ -z "$SAMPLES" ]; then
  echo "DONE $(date -u +%FT%TZ) $(wc -l < "$HV/lists/jobs.txt") plików ROOT v2.0 ${PTPROD}" > "$HV/DONE.tmp" && put_output "$HV/DONE.tmp" "$PARTS/DONE"; rm -f "$HV/DONE.tmp"
  echo "komplet → $PARTS/DONE"
else
  echo "brakujących części: $nmiss (SAMPLES=${SAMPLES:-wszystkie}) — bez DONE"; [ "$nmiss" -eq 0 ] || exit 1
fi
