#!/bin/bash
# submit_pp.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/pair_planes)
set -uo pipefail
PX=$(cd "$(dirname "$0")" && pwd); export PX
source "$PX/pp_common.sh"
cmd=${1:-check}
check() {
  local ok=0 kv s q n e
  [ -f "$PX/pp_extract.py" ] && [ -f "$PX/pp_merge.py" ] || { echo "STOP: brak pp_extract.py/pp_merge.py"; ok=1; }
  for kv in $PP_SAMPLES; do
    s=${kv%%:*}; n=$(ls "$H5V20/$s" 2>/dev/null | grep -c '_none\.h5$'); e=$(pp_kv "$s" "$PP_NFILES")
    echo "  $s: $n plików (oczekiwane $e)"
    [ "$n" -gt 0 ] && [ "$n" = "$e" ] || { echo "STOP: $s niekompletne"; ok=1; }
  done
  condor_q -batch 2>/dev/null | grep -qE 'charm_pp|pp\.dag' && { echo "STOP: PAIR_PLANES w kolejce"; ok=1; }
  [ -s "$PP_OUT/DONE" ] && { echo "STOP: już policzone ($(cat "$PP_OUT/DONE"))"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  PP_OUT=$PP_OUT PP_SAMPLES=$PP_SAMPLES"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_pp.sh submit)" ;;
  submit) check || exit 1
          "$PX/make_dag_pp.sh" || exit 1
          chmod +x "$PX"/*.sh "$PX"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then DRYRUN=1 "$PX/pp_job.sh" fd_nd 0 3; DRYRUN=1 "$PX/pp_job.sh" finish; exit 0; fi
          ( pp_env ) || exit 1
          cd "$PX" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_pp pp.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_pp|Total for klyczak'; ls -la "$PP_OUT" 2>/dev/null | head -30; ls "$PP_OUT/parts" 2>/dev/null | wc -l ;;
  fetch)  find "$PP_OUT" -maxdepth 1 \( -name '*.npz' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
