#!/bin/bash
# submit_gn.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/gen_norm)
#  check  — ROOT mb_nd (80) i 428770 (v2.0: 97 *dpt3to150*; v1.1: 39, w tym DATA_REF_TAGS) widoczne, v2.0: format (nSelBL) i kotwice h5 v1.1 (REF_ANCHORS), h5 v1.1 + nominalne scores incl_mc (D0_H, Dplus_H),
#           BASE_CFG, src, nic w kolejce (charm_gn), brak $GN_OUT/DONE (FORCE=1 → nowy bieg), kwota AFS < 85 %
#  submit — check + make_dag_gn.sh → condor_submit_dag -f (batch charm_gn) w CZYSTYM env (LCG psuje condor_dagman); DRYRUN=1 → DAG + komendy
#  status / fetch — kolejka + pliki $GN_OUT
# Zmienne: GN_FMT=v20|v11 PTPROD=dpt3to150 RAW_GN GN_OUT DATA_SCOPE=all|ref DATA_REF_TAGS=000021,000030 MAX_EVENTS (test) FORCE MAXJOBS=40 NOMINAL BASE_CFG
set -uo pipefail
GN=$(cd "$(dirname "$0")" && pwd); export GN
source "$GN/gn_common.sh"
cmd=${1:-check}
check() {
  local ok=0 n t sp
  gn_lists || ok=1
  n=$(wc -l < "$GN/lists/mb_nd.txt" 2>/dev/null || echo 0)
  [ "$n" -eq "$N_MB_EXP" ] || { echo "STOP: mb_nd ROOT: $n plików (oczekiwane $N_MB_EXP jak etap1) — B1 się nie domknie"; ok=1; }
  [ "$n" -gt 0 ] || { echo "STOP: brak ROOT mb_nd w $GN_RAW/$MB_DIR"; ok=1; }
  for t in ${GN_MB_EXCLUDE//,/ }; do   # §3b: wykluczony tag nie może mieć kompletnego ROOT; h5 v1.1 tego tagu musi istnieć (nogi do usunięcia)
    grep -q "\._$t\." "$GN/lists/mb_nd.txt" && { echo "STOP: GN_MB_EXCLUDE=$t, a ROOT _$t jest na liście"; ok=1; }
    [ -s "/eos/user/k/klyczak/charm_h5_v11/parts/mb_nd/mb_nd_${t}_none.h5" ] || { echo "STOP: brak h5 v1.1 mb_nd_${t} (wykluczenie po tagu)"; ok=1; }
  done
  [ -n "${GN_MB_EXCLUDE:-}" ] && echo "  §3b: wykluczone mb_nd $GN_MB_EXCLUDE, oczekiwane $N_MB_EXP plików ROOT"
  n=$(wc -l < "$GN/lists/data_428770.txt" 2>/dev/null || echo 0); echo "  ROOT: mb_nd $(wc -l < "$GN/lists/mb_nd.txt") plików, data_428770 $n plików"
  [ "$n" -eq "$N_DATA_EXP" ] || { echo "STOP: data_428770 ROOT: $n plików (oczekiwane $N_DATA_EXP)"; ok=1; }
  if [ "$GN_FMT" = v20 ]; then
    for f in ${REF_ANCHORS//,/ }; do [ -s "$f" ] || { echo "STOP: brak kotwic referencji w2 $f (h5 v1.1)"; ok=1; }; done
    grep -q dpt2to3 "$GN/lists/"*.txt && { echo "STOP: listy zawierają dpt2to3"; ok=1; }
    ( q_env >/dev/null 2>&1; python3 - "$(head -1 "$GN/lists/mb_nd.txt")" "$(head -1 "$GN/lists/data_428770.txt")" <<'PY' ) || { echo "STOP: pliki nie wyglądają na v2.0"; ok=1; }
import sys, uproot
for f in sys.argv[1:]:
    k = set(uproot.open(f)["CharmAnalysis"].keys())
    assert "EventInfo_nSelBL_NOSYS" in k and "EventInfo_leadTrackPtBL_NOSYS" in k, f"brak nSelBL/leadTrackPtBL: {f}"
print("  v2.0 OK (nSelBL, leadTrackPtBL)")
PY
  fi
  for t in ${DATA_REF_TAGS//,/ }; do grep -q "\._$t\." "$GN/lists/data_428770.txt" || { echo "STOP: brak pliku danych _$t (referencja w2)"; ok=1; }; done
  [ -s "$BASE_CFG" ] || { echo "STOP: brak $BASE_CFG"; ok=1; }
  for sp in D0_H Dplus_H; do [ -s "$NOMINAL/scores/$sp/incl_mc.h5" ] || { echo "STOP: brak $NOMINAL/scores/$sp/incl_mc.h5"; ok=1; }; done
  ls /eos/user/k/klyczak/charm_h5_v11/parts/mb_nd/mb_nd_*_none.h5 >/dev/null 2>&1 || { echo "STOP: brak h5 v1.1 mb_nd (scores liczone na v1.1)"; ok=1; }
  [ -f "$GN/src/charmtag/pairs/build.py" ] && [ -f "$GN/ref/mb_rates_2d.json" ] || { echo "STOP: paczka niekompletna (src/ lub ref/)"; ok=1; }
  condor_q -batch 2>/dev/null | grep -qE 'charm_gn|gn.dag' && { echo "STOP: gen_norm w kolejce"; ok=1; }
  [ -s "$GN_OUT/DONE" ] && { echo "gen_norm już policzony ($(cat "$GN_OUT/DONE")) — FORCE=1 → nowy bieg (DONE → DONE.old_*)"; [ "${FORCE:-0}" = 1 ] || [ "${RECOMBINE:-0}" = 1 ] || ok=1; }
  local q; q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  GN_FMT=$GN_FMT PTPROD=$PTPROD GN_OUT=$GN_OUT DATA_SCOPE=$DATA_SCOPE DATA_REF_TAGS=$DATA_REF_TAGS NOMINAL=$NOMINAL MAX_EVENTS=${MAX_EVENTS:-}"; return $ok
}
case $cmd in
  check)  check && echo "OK — można wysłać (./submit_gn.sh submit)" ;;
  submit) check || exit 1
          chmod +x "$GN"/*.sh "$GN"/*.py
          "$GN/make_dag_gn.sh" || exit 1
          if [ "${DRYRUN:-0}" = 1 ]; then
            f=$(head -1 "$GN/lists/mb_nd.txt"); d=$(head -1 "$GN/lists/data_428770.txt")
            DRYRUN=1 "$GN/gn_job.sh" root mc "$(job_tag "$f")" "$f"; DRYRUN=1 "$GN/gn_job.sh" root data "$(job_tag "$d")" "$d"
            DRYRUN=1 "$GN/gn_job.sh" legs - incl_mc -; DRYRUN=1 "$GN/gn_job.sh" combine - all -; exit 0
          fi
          ( gn_env ) || exit 1   # test importów w podpowłoce: LCG view psuje condor_dagman (GLIBCXX_3.4.31)
          { [ "${FORCE:-0}" = 1 ] || [ "${RECOMBINE:-0}" = 1 ]; } && [ -f "$GN_OUT/DONE" ] && mv -f "$GN_OUT/DONE" "$GN_OUT/DONE.old_$(date +%s)"
          cd "$GN" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_gn -maxjobs "$MAXJOBS" gn.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_gn|Total for klyczak'
          echo "parts: mc $(ls "$GN_OUT/parts/mc" 2>/dev/null | wc -l), data $(ls "$GN_OUT/parts/data" 2>/dev/null | wc -l); legs: $(ls "$GN_OUT/legs" 2>/dev/null | tr '\n' ' ')"
          ls -la "$GN_OUT" 2>/dev/null | head -30; [ -f "$GN_OUT/acceptance.json" ] && grep -m1 verdict "$GN_OUT/acceptance.json" ;;
  fetch)  find "$GN_OUT" -maxdepth 1 -type f 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
