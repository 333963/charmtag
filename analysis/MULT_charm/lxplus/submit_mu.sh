#!/bin/bash
# submit_mu.sh check|submit|status   (na lxplus, w ~/private/charm_q/Q_lxplus/mult_charm; normalnie przez ORCH, etap `mult_charm` — domyślnie ZABLOKOWANY)
#  check  — ROOT v2.0 widoczne (mb_nd 80, data_428770 97, data_427929 ~583 przy PTPROD=dpt3to150), format v2.0 (nSelBL) w 1. pliku każdej próbki,
#           nic w kolejce (charm_mu), brak $MU_OUT/DONE (FORCE=1), kwota AFS < 85 %
#  submit — check + make_dag_mu.sh → condor_submit_dag -f w CZYSTYM env; DRYRUN=1 → DAG + komendy
set -uo pipefail
MU=$(cd "$(dirname "$0")" && pwd); export MU
source "$MU/mu_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s n f
  mu_lists
  for s in $SAMPLES; do
    n=$(wc -l < "$MU/lists/$s.txt"); echo "  $s: $n plików ($PTPROD)"
    case $s in mb_nd) [ "$n" -eq 80 ] || echo "UWAGA: mb_nd $n ≠ 80" ;; data_428770) [ "$n" -eq 97 ] || echo "UWAGA: 428770 $n ≠ 97" ;; esac
    [ "$n" -gt 0 ] || { if [ "$s" = mb_nd ]; then echo "STOP: brak plików mb_nd"; ok=1; else echo "UWAGA: brak plików $s (pominięte w analizie)"; fi; }
    f=$(head -1 "$MU/lists/$s.txt")
    [ -n "$f" ] && ( mu_env >/dev/null 2>&1 && python3 - "$f" <<'PY' ) || { [ -n "$f" ] && echo "STOP: $s nie wygląda na v2.0"; [ -n "$f" ] && ok=1; }
import sys, uproot
t = uproot.open(sys.argv[1])["CharmAnalysis"]; k = set(t.keys())
assert "EventInfo_nSelBL_NOSYS" in k and "EventInfo_leadTrackPtBL_NOSYS" in k, "brak nSelBL/leadTrackPtBL"
print("   v2.0 OK:", sys.argv[1].split("/")[-1], t.num_entries, "wpisów")
PY
  done
  condor_q -batch 2>/dev/null | grep -qE 'charm_mu|mu.dag' && { echo "STOP: mult_charm w kolejce"; ok=1; }
  [ -s "$MU_OUT/DONE" ] && { echo "już policzone ($(cat "$MU_OUT/DONE")) — FORCE=1"; [ "${FORCE:-0}" = 1 ] || [ "${REANALYZE:-0}" = 1 ] || ok=1; }
  local q; q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  MU_OUT=$MU_OUT"; return $ok
}
case $cmd in
  check)  check && echo "OK — można wysłać (./submit_mu.sh submit)" ;;
  submit) check || exit 1
          chmod +x "$MU"/*.sh "$MU"/*.py
          "$MU/make_dag_mu.sh" || exit 1
          if [ "${DRYRUN:-0}" = 1 ]; then f=$(head -1 "$MU/lists/mb_nd.txt"); DRYRUN=1 "$MU/mu_job.sh" extract mc mb_nd "$(job_tag "$f")" "$f"; DRYRUN=1 "$MU/mu_job.sh" analyze - - all -; exit 0; fi
          ( mu_env ) || exit 1
          if [ "${FORCE:-0}" = 1 ] || [ "${REANALYZE:-0}" = 1 ]; then ts=$(date +%s)   # stary wynik odkładany, nie kasowany
            [ -f "$MU_OUT/DONE" ] && mv -f "$MU_OUT/DONE" "$MU_OUT/DONE.old_$ts"; [ -f "$MU_OUT/res.tgz" ] && mv -f "$MU_OUT/res.tgz" "$MU_OUT/res.tgz.old_$ts"; fi
          cd "$MU" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_mu -maxjobs "$MAXJOBS" mu.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_mu|Total for klyczak'
          for s in $SAMPLES; do echo "parts $s: $(ls "$MU_OUT/parts/$s" 2>/dev/null | wc -l)"; done; ls -la "$MU_OUT" 2>/dev/null | head ;;
  *) echo "użycie: $0 check|submit|status"; exit 2 ;;
esac
