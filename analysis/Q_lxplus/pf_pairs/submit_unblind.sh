#!/bin/bash
# submit_unblind.sh check|arm|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/pf_pairs)
#  check  — pf DONE na EOS, pf/unblind nie w kolejce, brak unblind/DONE, kwota AFS <85 %, cfg/pf.yaml jest, flaga UNBLIND_OK (lokalna kopia z Maca)
#  arm    — make_dag_unblind.sh + zapis $UB_OUT/UNBLIND_OK (treść flagi z Maca + cfg_sha256) — DOPIERO PO DECYZJI
#  submit — check + arm → condor_submit_dag unblind.dag (batch charm_unblind)
#  status / fetch — kolejka + lista plików unblind/*.md|json
set -uo pipefail
PF=$(cd "$(dirname "$0")" && pwd); export PF
source "$PF/pf_common.sh"
UB_OUT=${UB_OUT:-$PF_OUT/unblind}; FLAG=$PF/UNBLIND_OK
cmd=${1:-check}
check() {
  local ok=0
  [ -s "$PF_OUT/DONE" ] || { echo "STOP: pf_pairs nie skończony (brak $PF_OUT/DONE)"; ok=1; }
  [ -s "$PF_OUT/compare_pf.json" ] || { echo "STOP: brak compare_pf.json"; ok=1; }
  condor_q -batch 2>/dev/null | grep -qE 'charm_pf|charm_unblind|unblind.dag' && { echo "STOP: pf/unblind w kolejce"; ok=1; }
  [ -s "$UB_OUT/DONE" ] && { echo "STOP: odsłonięcie już wykonane ($(cat "$UB_OUT/DONE")) — nie powtarzam (FORCE nie istnieje celowo)"; ok=1; }
  [ -s "$PF/cfg/pf.yaml" ] || { echo "STOP: brak cfg/pf.yaml"; ok=1; }
  [ -s "$FLAG" ] || { echo "STOP: brak flagi $FLAG (decyzja o odsłonięciu: touch charm_work/ORCH/UNBLIND_OK na Macu → rsync)"; ok=1; }
  grep -q '"track_id": "d0z0"' "$PF/src/charmtag/pairs/build.py" || { echo "STOP: src bez PF"; ok=1; }
  local q; q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  NORM=$NORM  UB_OUT=$UB_OUT"; return $ok
}
arm() {
  "$PF/make_dag_unblind.sh" || return 1
  mkdir -p "$UB_OUT"
  { cat "$FLAG"; echo "armed $(date -u +%FT%TZ) on $(hostname)"; echo "cfg_sha256 $(sha256sum "$PF/cfg/pf_unblind.yaml" | cut -c1-64)"
    echo "pairs_src_sha256 $(cat $(ls "$PF"/src/charmtag/pairs/*.py | sort) | sha256sum | cut -c1-64)"
    echo "compare_pf_sha256 $(sha256sum "$PF_OUT/compare_pf.json" | cut -c1-64)"; } > "$UB_OUT/UNBLIND_OK"
  echo "uzbrojone: $UB_OUT/UNBLIND_OK"; cat "$UB_OUT/UNBLIND_OK"
}
case $cmd in
  check)  check && echo "OK — gotowe do odsłonięcia (./submit_unblind.sh submit)" ;;
  arm)    check || exit 1; arm ;;
  submit) check || exit 1; arm || exit 1
          chmod +x "$PF"/*.sh "$PF"/tools/*.py
          ( pf_env && python3 -c "import charmtag.pairs.fit2d; print('import OK')" ) || exit 1
          cd "$PF" && condor_submit_dag -batch-name charm_unblind unblind.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_unblind|Total for klyczak'; ls -la "$UB_OUT" "$UB_OUT/pairs" 2>/dev/null | head -30
          [ -f "$UB_OUT/tfin.md" ] && head -80 "$UB_OUT/tfin.md" ;;
  fetch)  find "$UB_OUT" -maxdepth 3 \( -name '*.md' -o -name '*.json' -o -name 'UNBLIND_OK' \) 2>/dev/null ;;
  *) echo "użycie: $0 check|arm|submit|status|fetch"; exit 2 ;;
esac
