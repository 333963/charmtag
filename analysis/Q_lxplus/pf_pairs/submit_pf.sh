#!/bin/bash
# submit_pf.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/pf_pairs)
#  check  — główny DAG charmtag NIE w kolejce (IGNORE_MAIN=1 gdy zostały tylko toys), nominalne wyjścia są
#           (models/scores D0+D⁺, validation, pairs/fit2d_blind.json), PF nie w kolejce i brak DONE, kwota AFS < 85%, src z poprawką PF
#  submit — check → make_dag_pf.sh → test importu → condor_submit_dag pf.dag (batch charm_pf)
#  status — condor_q charm_pf + wyjścia EOS ; fetch — lista plików do pobrania (rsync robi START na Macu)
set -uo pipefail
PF=$(cd "$(dirname "$0")" && pwd); export PF
source "$PF/pf_common.sh"
cmd=${1:-check}
check() {
  local ok=0
  if condor_q -batch 2>/dev/null | grep -q 'charmtag.dag'; then
    if [ "${IGNORE_MAIN:-0}" = 1 ]; then echo "UWAGA: główny DAG w kolejce, IGNORE_MAIN=1"
    else echo "STOP: główny DAG charmtag.dag w kolejce (IGNORE_MAIN=1 gdy zostały tylko toys)"; ok=1; fi
  fi
  condor_q -batch 2>/dev/null | grep -q 'charm_pf\|pf.dag' && { echo "STOP: PF już w kolejce — nie dubluję"; ok=1; }
  [ -f "$PF_OUT/DONE" ] && { echo "PF skończony ($(cat "$PF_OUT/DONE")) — $PF_OUT/compare_pf.md; FORCE=1 → nowy bieg"; [ "${FORCE:-0}" = 1 ] || ok=1; }
  for f in models/D0_H/meta.json models/Dplus_H/meta.json scores/D0_H/data_none.h5 scores/Dplus_H/data_none.h5 \
           scores/D0_H/incl_mc.h5 scores/Dplus_H/incl_mc.h5 scores/D0_H/signal_mc.h5 scores/Dplus_H/signal_mc.h5 pairs/fit2d_blind.json; do
    [ -s "$NOMINAL/$f" ] || { echo "STOP: brak nominalnego $NOMINAL/$f"; ok=1; }
  done
  grep -q '"track_id": "d0z0"' "$PF/src/charmtag/pairs/build.py" || { echo "STOP: src bez poprawki PF (pairs.track_id d0z0)"; ok=1; }
  [ "$NORM" = truth ] && { grep -q truth_norm_json "$BASE_CFG" || echo "UWAGA: NORM=truth, ale brak pairs.truth_norm_json w $BASE_CFG"; }
  local q; q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}% ≥ 85%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  NORM=$NORM  PF_OUT=$PF_OUT"; return $ok
}
case $cmd in
  check)  check && echo "OK — można wysłać (./submit_pf.sh submit)" ;;
  submit) check || exit 1
          chmod +x "$PF"/*.sh "$PF"/tools/*.py
          "$PF/make_dag_pf.sh" || exit 1
          ( pf_env && python3 -c "import charmtag.pairs.build, charmtag.pairs.mixctl, charmtag.syst; print('import OK')" ) || exit 1
          cd "$PF" && condor_submit_dag -batch-name charm_pf -maxjobs 20 pf.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_pf|Total for klyczak'; ls -la "$PF_OUT" "$PF_OUT/pairs" 2>/dev/null | head -40
          [ -f "$PF_OUT/compare_pf.md" ] && head -60 "$PF_OUT/compare_pf.md" ;;
  fetch)  find "$PF_OUT" -maxdepth 3 \( -name '*.md' -o -name '*.json' \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
