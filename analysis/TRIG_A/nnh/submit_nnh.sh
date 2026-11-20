#!/bin/bash
# submit_nnh.sh check|submit|status   (na lxplus, w ~/private/charm_q/Q_lxplus/nnh)
set -uo pipefail
NNH=$(cd "$(dirname "$0")" && pwd); Q=$(cd "$NNH/.." && pwd)
source "$Q/common.sh"
: "${NNH_IN:=$EOS_OUT/charm_q/b_trig}" "${NNH_OUT:=$EOS_OUT/charm_q/nnh}"
check() {
  local ok=0 q f s
  for f in nnh_train.py nnh_job.sh nnh.sub btrig_common.py src/charmtag/dstar_refit.py; do [ -f "$NNH/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  for s in fd_nd mb_nd data_428770 data_427929; do
    [ -s "$NNH_IN/bt_$s.npz" ] || [ "$(ls "$NNH_IN/parts" 2>/dev/null | grep -c "^bt_${s}_")" -gt 0 ] || { echo "STOP: brak bt_$s"; ok=1; }; done
  condor_q -batch 2>/dev/null | grep -q 'charm_nnh' && { echo "STOP: nnh w kolejce"; ok=1; }
  [ -s "$NNH_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone (FORCE=1)"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  NNH_IN=$NNH_IN NNH_OUT=$NNH_OUT"; return $ok
}
case ${1:-check} in
  check)  check && echo "OK — gotowe (./submit_nnh.sh submit)" ;;
  submit) check || exit 1; mkdir -p "$Q/logs/nnh"; chmod +x "$NNH"/*.sh "$NNH"/*.py
          cd "$NNH" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit nnh.sub ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_nnh|Total for klyczak' || true; exit 0 ;;
  *) echo "użycie: $0 check|submit|status"; exit 2 ;;
esac
