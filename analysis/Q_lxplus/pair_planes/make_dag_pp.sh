#!/bin/bash
# Generuje pp.dag: węzeł na (próbka, część) → finish. NIE wysyła.
set -euo pipefail
PX=$(cd "$(dirname "$0")" && pwd); export PX
source "$PX/pp_common.sh"
mkdir -p "$Q/logs/pair_planes"
S=$PX/pp.sub
PENV="PP_OUT=$PP_OUT H5V20=$H5V20 PP_SAMPLES=$(echo $PP_SAMPLES | tr ' ' ,) PP_NFILES=$(echo $PP_NFILES | tr ' ' ,)"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych"; exit 1 ;; esac
DAG=${PP_DAG:-$PX/pp.dag}
V="penv=\"$PENV\""
nodes=()
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV"
  for kv in $PP_SAMPLES; do
    s=${kv%%:*}; n=${kv##*:}
    for ((k = 0; k < n; k++)); do
      j=pp_${s}_$k; nodes+=("$j")
      echo "JOB $j $S"; echo "VARS $j step=\"$s\" part=\"$k\" np=\"$n\" cpus=\"1\" mem=\"6GB\" disk=\"$PP_DISK\" flavour=\"workday\" $V"; echo "RETRY $j 2"
    done
  done
  echo "JOB pp_finish $S"; echo "VARS pp_finish step=\"finish\" part=\"0\" np=\"1\" cpus=\"1\" mem=\"4GB\" disk=\"$PP_DISK\" flavour=\"longlunch\" $V"; echo "RETRY pp_finish 1"
  echo "PARENT ${nodes[*]} CHILD pp_finish"
} > "$DAG"
echo "pp.dag ($DAG): ${#nodes[@]} węzłów + finish → $PP_OUT"
