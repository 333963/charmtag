#!/bin/bash
# Generuje ubjk.dag: prep → fit_0..fit_{N-1} (+ ref) → finish. NIE wysyła.
# prep: gdy $UJ_LEGS (legs_data.parquet z odsłonięcia) jest na EOS → lekki (4 GB); gdy nie → przebudowa par w scratch (16 GB, 4 CPU, jak ub_pairs).
set -euo pipefail
UJ=$(cd "$(dirname "$0")" && pwd); export UJ
source "$UJ/ubjk_common.sh"
mkdir -p "$Q/logs/ubjk"
case "$UJ_NCHUNK" in ''|*[!0-9]*|0) echo "STOP: UJ_NCHUNK=$UJ_NCHUNK"; exit 1 ;; esac
REBUILD=${UJ_REBUILD:-}
if [ -z "$REBUILD" ]; then if [ -s "$UJ_LEGS" ]; then REBUILD=0; else REBUILD=1; fi; fi
S=$UJ/ubjk.sub
PENV="UB_OUT=$UB_OUT UJ_GROUPS=$UJ_GROUPS UJ_OUT=$UJ_OUT UJ_NCHUNK=$UJ_NCHUNK UJ_MODES=$UJ_MODES UJ_REBUILD=$REBUILD NORM=$NORM OMP_NUM_THREADS=1"
case "$PENV" in *\'*|*\"*) echo "STOP: cudzysłów w zmiennych biegu"; exit 1 ;; esac
V="penv=\"$PENV\""
if [ "$REBUILD" = 1 ]; then pres='cpus="4" mem="16GB" flavour="workday"'; else pres='cpus="1" mem="4GB" flavour="longlunch"'; fi
{
  echo "# wygenerowane $(date -u +%FT%TZ): $PENV UJ_REF=$UJ_REF (ODSŁONIĘTE — jackknife μ po plikach; blokada UNBLIND_OK w węźle)"
  echo "JOB ubjk_prep $S"; echo "VARS ubjk_prep step=\"prep\" chunk=\"-\" $pres $V"; echo "RETRY ubjk_prep 1"
  kids=""
  for ((i = 0; i < UJ_NCHUNK; i++)); do
    echo "JOB ubjk_fit_$i $S"; echo "VARS ubjk_fit_$i step=\"fit\" chunk=\"$i\" cpus=\"1\" mem=\"4GB\" flavour=\"workday\" $V"; echo "RETRY ubjk_fit_$i 1"
    kids="$kids ubjk_fit_$i"
  done
  if [ "$UJ_REF" = 1 ]; then
    echo "JOB ubjk_ref $S"; echo "VARS ubjk_ref step=\"ref\" chunk=\"-\" cpus=\"1\" mem=\"4GB\" flavour=\"workday\" $V"; echo "RETRY ubjk_ref 1"
    kids="$kids ubjk_ref"
  fi
  echo "JOB ubjk_finish $S"; echo "VARS ubjk_finish step=\"finish\" chunk=\"-\" cpus=\"1\" mem=\"2GB\" flavour=\"espresso\" $V"
  echo "PARENT ubjk_prep CHILD$kids"
  echo "PARENT$kids CHILD ubjk_finish"
} > "$UJ/ubjk.dag"
echo "ubjk.dag: prep (rebuild=$REBUILD) → $UJ_NCHUNK × fit$([ "$UJ_REF" = 1 ] && echo ' + ref') → finish; wyjście $UJ_OUT"
