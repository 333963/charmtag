#!/bin/bash
# make_dag_mu.sh — mu.dag: extract × (mb_nd 80 + data_427929 + data_428770 97) → analyze. NIE wysyła.
set -euo pipefail
MU=$(cd "$(dirname "$0")" && pwd); export MU
source "$MU/mu_common.sh"
mu_lists
mkdir -p "$Q/logs/mult_charm"
S=$MU/mu.sub; D=$MU/mu.dag
PENV="MU_OUT=$MU_OUT RAW20=$RAW20 PTPROD=$PTPROD FORCE=${FORCE:-0}"
[ -n "${MAX_EVENTS:-}" ] && PENV+=" MAX_EVENTS=$MAX_EVENTS"
V="penv=\"$PENV\""
nodes=(); n=0
{
  echo "# wygenerowane $(date -u +%FT%TZ): MU_OUT=$MU_OUT RAW20=$RAW20 PTPROD=$PTPROD SAMPLES=$SAMPLES"
  for s in $SAMPLES; do
    [ "${REANALYZE:-0}" = 1 ] && break                       # 19.11: tylko analyze na istniejących parts (poprawka fitów mult_fit/mult_analyze)
    k=data; case $s in mb_*|fd_*) k=mc ;; esac
    while read -r f; do
      t=$(job_tag "$f"); [ -n "$t" ] || { echo "STOP: brak numeru joba w $f" >&2; exit 1; }
      nm="x_${s}_$t"
      echo "JOB $nm $S"; echo "VARS $nm step=\"extract\" kind=\"$k\" smp=\"$s\" tag=\"$t\" src=\"$f\" cpus=\"1\" mem=\"4GB\" flavour=\"longlunch\" $V"
      nodes+=("$nm"); n=$((n+1))
    done < "$MU/lists/$s.txt"
  done
  echo "JOB analyze $S"; echo "VARS analyze step=\"analyze\" kind=\"-\" smp=\"all\" tag=\"all\" src=\"-\" cpus=\"1\" mem=\"16GB\" flavour=\"tomorrow\" $V"
  [ ${#nodes[@]} -gt 0 ] && echo "PARENT ${nodes[*]} CHILD analyze"
  echo "RETRY ALL_NODES 1"
} > "$D"
echo "mu.dag: $n × extract + analyze → $MU_OUT${REANALYZE:+ (REANALYZE=$REANALYZE)}"
