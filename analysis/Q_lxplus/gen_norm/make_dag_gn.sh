#!/bin/bash
# make_dag_gn.sh — gn.dag: root_mc × 80 (mb_nd) + root_data × 97 (v2.0; v1.1: 39|2) (428770) + legs (h5+scores) → combine. NIE wysyła.
set -euo pipefail
GN=$(cd "$(dirname "$0")" && pwd); export GN
source "$GN/gn_common.sh"
gn_lists
mkdir -p "$Q/logs/gen_norm"
S=$GN/gn.sub; D=$GN/gn.dag
PENV="GN_FMT=$GN_FMT PTPROD=$PTPROD RAW_GN=$RAW_GN GN_OUT=$GN_OUT BASE_CFG=$BASE_CFG NOMINAL=$NOMINAL DATA_REF_TAGS=$DATA_REF_TAGS REF_ANCHORS=$REF_ANCHORS REF_RATES=$REF_RATES FORCE=${FORCE:-0}"
[ -n "${GN_MB_EXCLUDE:-}" ] && PENV+=" GN_MB_EXCLUDE=$GN_MB_EXCLUDE"
[ -n "${MAX_EVENTS:-}" ] && PENV+=" MAX_EVENTS=$MAX_EVENTS"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
V="penv=\"$PENV\""
nodes=(); nm=0; nd=0
{
  echo "# wygenerowane $(date -u +%FT%TZ): GN_FMT=$GN_FMT PTPROD=$PTPROD GN_OUT=$GN_OUT DATA_SCOPE=$DATA_SCOPE DATA_REF_TAGS=$DATA_REF_TAGS NOMINAL=$NOMINAL MAX_EVENTS=${MAX_EVENTS:-}"
  while read -r f; do
    t=$(job_tag "$f"); [ -n "$t" ] || { echo "STOP: brak numeru joba w $f" >&2; exit 1; }
    echo "JOB mc_$t $S"; echo "VARS mc_$t step=\"root\" kind=\"mc\" tag=\"$t\" src=\"$f\" cpus=\"1\" mem=\"3GB\" flavour=\"microcentury\" $V"
    nodes+=("mc_$t"); nm=$((nm+1))
  done < "$GN/lists/mb_nd.txt"
  while read -r f; do
    t=$(job_tag "$f"); [ -n "$t" ] || { echo "STOP: brak numeru joba w $f" >&2; exit 1; }
    if [ "$GN_FMT" = v11 ] && [ "$DATA_SCOPE" = ref ] && ! echo ",$DATA_REF_TAGS," | grep -q ",$t,"; then continue; fi
    echo "JOB da_$t $S"; echo "VARS da_$t step=\"root\" kind=\"data\" tag=\"$t\" src=\"$f\" cpus=\"1\" mem=\"3GB\" flavour=\"microcentury\" $V"
    nodes+=("da_$t"); nd=$((nd+1))
  done < "$GN/lists/data_428770.txt"
  echo "JOB legs $S"; echo "VARS legs step=\"legs\" kind=\"-\" tag=\"incl_mc\" src=\"-\" cpus=\"1\" mem=\"8GB\" flavour=\"workday\" $V"
  nodes+=(legs)
  echo "JOB combine $S"; echo "VARS combine step=\"combine\" kind=\"-\" tag=\"all\" src=\"-\" cpus=\"1\" mem=\"12GB\" flavour=\"longlunch\" $V"
  echo "PARENT ${nodes[*]} CHILD combine"
  echo "RETRY ALL_NODES 1"
} > "$D"
echo "gn.dag: $nm × mc + $nd × data + legs + combine → $GN_OUT"
