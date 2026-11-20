#!/bin/bash
# submit_u7.sh check | submit | status | dag — DAG UB427JK: prep (przebudowa nóg 427929, 16 GB) → fit_0..fit_{N−1} → finish.
# Submit w czystym env (LCG psuje condor_dagman), condor_submit_dag -f. DRYRUN=1 → tylko DAG.
set -uo pipefail
U7=$(cd "$(dirname "$0")" && pwd); export U7
source "$U7/u7_common.sh"
mkdir -p "$Q/logs/ub427jk"
cmd=${1:-check}
check() {
  local ok=0
  u7_gate || ok=1
  case "$U7_NCHUNK" in ''|*[!0-9]*|0) echo "STOP: U7_NCHUNK=$U7_NCHUNK"; ok=1 ;; esac
  [ "$U7_G" -ge 20 ] 2>/dev/null || { echo "STOP: U7_G=$U7_G (< 20)"; ok=1; }
  [ -s "$U4_NOM_CFG" ] || { echo "STOP: brak $U4_NOM_CFG"; ok=1; }
  for f in "$NOMINAL/scores/D0_H/signal_mc.h5" "$NOMINAL/scores/Dplus_H/incl_mc.h5"; do [ -s "$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -e "$U7_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: $U7_OUT/DONE istnieje (FORCE=1 nadpisuje)"; ok=1; }
  echo "AFS: $(fs lq ~ 2>/dev/null | awk 'NR==2{print $4}')"
  echo "joby charm_ub427jk w kolejce: $(condor_q -batch 2>/dev/null | grep -c charm_ub427jk)"
  [ $ok = 0 ] && echo "OK — można wysłać (G=$U7_G, porcje $U7_NCHUNK, wyjście $U7_OUT)"; return $ok
}
dag() {
  local P="penv=\"U7_OUT=$U7_OUT U7_G=$U7_G U7_NCHUNK=$U7_NCHUNK U7_MODES=$U7_MODES NORM=$NORM OMP_NUM_THREADS=1\"" i kids=""
  echo "JOB prep $U7/u7.sub"; echo "VARS prep step=\"prep\" chunk=\"-\" cpus=\"4\" mem=\"16GB\" disk=\"20GB\" flavour=\"workday\" $P"; echo "RETRY prep 1"
  for ((i = 0; i < U7_NCHUNK; i++)); do
    echo "JOB fit_$i $U7/u7.sub"; echo "VARS fit_$i step=\"fit\" chunk=\"$i\" cpus=\"1\" mem=\"4GB\" disk=\"2GB\" flavour=\"workday\" $P"; echo "RETRY fit_$i 1"; kids="$kids fit_$i"
  done
  echo "JOB finish $U7/u7.sub"; echo "VARS finish step=\"finish\" chunk=\"-\" cpus=\"1\" mem=\"2GB\" disk=\"1GB\" flavour=\"espresso\" $P"
  echo "PARENT prep CHILD$kids"; echo "PARENT$kids CHILD finish"
}
case $cmd in
  check) check ;;
  submit)
    check || exit 1
    dag > "$U7/u7.dag"; echo "u7.dag: $(grep -c '^JOB' "$U7/u7.dag") węzłów"
    [ "${DRYRUN:-0}" = 1 ] && exit 0
    cd "$U7" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_ub427jk u7.dag ;;
  status) condor_q -batch 2>/dev/null | grep charm_ub427jk; ls "$U7_OUT" 2>/dev/null | head; ls "$U7_OUT" 2>/dev/null | grep -c '^fit_' ;;
  dag) dag ;;
  *) echo "użycie: $0 check|submit|status|dag"; exit 2 ;;
esac
