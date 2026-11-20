#!/bin/bash
# make_dag_h5v20.sh — h5v20.dag: węzeł per (próbka, kawałek ≤ CHUNK plików) tylko z BRAKUJĄCYMI częściami + węzeł finish (meta, DONE). NIE wysyła.
set -euo pipefail
HV=$(cd "$(dirname "$0")" && pwd); export HV
source "$HV/h5v20_common.sh"
hv_lists >/dev/null || true
mkdir -p "$Q/logs/h5_v20" "$HV/lists/chunks"
rm -f "$HV/lists/chunks/"*.txt 2>/dev/null || true   # własne, generowane listy kawałków (nie dane)
PENV="H5V20=$H5V20 RAW20=$RAW20 PTPROD=$PTPROD STEP=$STEP SAMPLES=$SAMPLES FORCE=${FORCE:-0}"
[ -n "${MAX_EVENTS:-}" ] && PENV+=" MAX_EVENTS=$MAX_EVENTS"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych"; exit 1 ;; esac
D=$HV/h5v20.dag; nodes=(); nj=0
{
  echo "# wygenerowane $(date -u +%FT%TZ): RAW20=$RAW20 H5V20=$H5V20 PTPROD=$PTPROD CHUNK=$CHUNK SAMPLES=${SAMPLES:-wszystkie} MAX_EVENTS=${MAX_EVENTS:-}"
  while read -r k _ _ _; do
    grep "^$k " "$HV/lists/missing.txt" | awk '{gsub(",", "+", $4); print $1, $2, $3, $4}' > "$HV/lists/chunks/$k.all" || true
    [ -s "$HV/lists/chunks/$k.all" ] || { rm -f "$HV/lists/chunks/$k.all"; continue; }
    split -l "$CHUNK" -d -a 3 --additional-suffix=.txt "$HV/lists/chunks/$k.all" "$HV/lists/chunks/${k}_c"
    rm -f "$HV/lists/chunks/$k.all"
    for c in "$HV/lists/chunks/${k}"_c*.txt; do
      n=$(basename "$c" .txt); echo "JOB $n $HV/h5v20.sub"; echo "VARS $n list=\"$c\" penv=\"$PENV\""
      nodes+=("$n"); nj=$((nj + $(wc -l < "$c")))
    done
  done < <(hv_samples)
  echo "FINAL finish $HV/finish.sub"; echo "VARS finish penv=\"$PENV\""
  echo "RETRY ALL_NODES 1"
} > "$D"
echo "h5v20.dag: ${#nodes[@]} węzłów (kawałków), $nj procesów (plików ROOT) + finish → $PARTS"
