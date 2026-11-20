#!/bin/bash
# submit_topo.sh check|submit|status   (lxplus, ~/private/charm_q/Q_lxplus/ub_topo). Batch charm_ubtopo; submit w czystym env (bez LCG).
set -uo pipefail
UT=$(cd "$(dirname "$0")" && pwd); export UT
source "$UT/topo_common.sh"
cmd=${1:-check}
check() {
  local ok=0
  [ -s "$UB_OUT/DONE" ] && [ -s "$UB_OUT/UNBLIND_OK" ] || { echo "STOP: brak odsłonięcia ($UB_OUT)"; ok=1; }
  [ -s "$PF/cfg/pf_unblind.yaml" ] || { echo "STOP: brak $PF/cfg/pf_unblind.yaml"; ok=1; }
  [ -f "$PF/src/charmtag/pairs/diff.py" ] || { echo "STOP: $PF/src bez charmtag/pairs/diff.py"; ok=1; }
  [ -s "$UT/ref/topo_gen.json" ] || { echo "STOP: brak ref/topo_gen.json (topo_gen.py lokalnie)"; ok=1; }
  [ -f "$UT/../pd_diff/tools/diff_fill.py" ] || { echo "STOP: brak pd_diff/tools/diff_fill.py"; ok=1; }
  condor_q -batch 2>/dev/null | grep -qE 'charm_ubtopo' && { echo "STOP: charm_ubtopo w kolejce"; ok=1; }
  [ -s "$UT_OUT/DONE" ] && { echo "STOP: już policzone ($UT_OUT/DONE)"; ok=1; }
  local q; q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "AFS ${q:-?}%  VAR=$UT_VAR CHANNELS=$UT_CHANNELS UT_OUT=$UT_OUT NORM=$NORM"; [ $ok = 0 ] && echo "OK — można wysłać"; return $ok
}
case $cmd in
  check)  check ;;
  submit) check || exit 1
          "$UT/make_dag_topo.sh" || exit 1
          chmod +x "$UT"/*.sh "$UT"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then for s in pairs fit2d report; do DRYRUN=1 "$UT/topo_job.sh" $s; done; exit 0; fi
          ( pd_env && python3 -c "import charmtag.pairs.diff; print('import OK')" ) || exit 1
          cd "$UT" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_ubtopo topo.dag ;;
  status) condor_q -batch 2>/dev/null | grep charm_ubtopo; ls -la "$UT_OUT" 2>/dev/null | head ;;
  *) echo "użycie: $0 check|submit|status"; exit 2 ;;
esac
