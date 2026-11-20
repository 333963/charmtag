#!/bin/bash
# submit_h5v20.sh check|submit|status   (na lxplus, w ~/private/charm_q/Q_lxplus/h5_v20)
#  check  — wejście ROOT v2.0 (DONE per próbka, liczba plików *dpt3to150*), format v2.0 (nSelBL, AFP), venv, kolejka (charm_h5v20),
#           kwota AFS < 85 %, ile części już na EOS ($H5V20/parts); kod wyjścia 0 = można wysłać (albo komplet)
#  submit — check + make_dag_h5v20.sh → condor_submit_dag -f (batch charm_h5v20) w CZYSTYM env (LCG psuje condor_dagman);
#           tylko brakujące części; komplet bez DONE → sam finish (DONE); DRYRUN=1 → DAG + komenda 1. procesu, bez wysyłki
#  status — kolejka, ogon dagman.out, brakujące części, meta/summary.json
# Zmienne: SAMPLES=mb_nd,fd_nd CHUNK=25 STEP=100MB MAX_EVENTS (test) FORCE=1 (nadpisz części) MAXJOBS=60 H5V20 RAW20 PTPROD
set -uo pipefail
HV=$(cd "$(dirname "$0")" && pwd); export HV
source "$HV/h5v20_common.sh"
cmd=${1:-check}
running() { condor_q -batch 2>/dev/null | grep -qE 'charm_h5v20|h5v20.dag'; }
check() {
  local ok=0 q f
  echo "== wejście $RAW20 (PTPROD=$PTPROD) → wyjście $PARTS"
  hv_lists || { echo "STOP: brak wejścia albo brak DONE pobierania (rucio_v20)"; ok=1; }
  grep -q dpt2to3 "$HV/lists/jobs.txt" && { echo "STOP: listy zawierają dpt2to3"; ok=1; }
  f=$(awk 'NR==1{print $2}' "$HV/lists/jobs.txt")
  if [ -n "$f" ]; then
    ( q_env >/dev/null 2>&1; python3 - "$f" <<'PY' ) || { echo "STOP: venv/format (uproot, h5py, gałęzie v2.0)"; ok=1; }
import sys, uproot, h5py, awkward
k = set(uproot.open(sys.argv[1])["CharmAnalysis"].keys())
need = ["EventInfo_nSelBL_NOSYS", "EventInfo_leadTrackPtBL_NOSYS", "EventInfo_xi_cal_A_NOSYS", "AfpProtons_xi_NOSYS", "AfpProtons_side_NOSYS", "EventInfo_lumiBlock"]
miss = [b for b in need if b not in k]
assert not miss, f"brak gałęzi v2.0 {miss}"
print("  format v2.0 OK (nSelBL, leadTrackPtBL, xi_cal_A, AfpProtons, lumiBlock); venv OK")
PY
  fi
  running && { echo "STOP: charm_h5v20 w kolejce"; ok=1; }
  q=$(fs lq ~ 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}'); q=${q:-0}
  echo "== kwota AFS: ${q}% (limit 85 %)"; [ "${q%.*}" -ge 85 ] && { echo "STOP: kwota AFS ≥ 85 %"; ok=1; }
  echo "== części: $(wc -l < "$HV/lists/jobs.txt") plików ROOT, brakujących $(wc -l < "$HV/lists/missing.txt"); DONE na EOS: $([ -s "$PARTS/DONE" ] && cat "$PARTS/DONE" || echo brak)"
  return $ok
}
case $cmd in
  check)  check && echo "OK — można wysłać (./submit_h5v20.sh submit)" ;;
  submit) check || exit 1
          [ -s "$PARTS/DONE" ] && [ "$(wc -l < "$HV/lists/missing.txt")" -eq 0 ] && { echo "== komplet + DONE — nic do wysłania"; exit 0; }
          chmod +x "$HV"/*.sh "$HV"/*.py
          "$HV/make_dag_h5v20.sh" || exit 1
          if [ "${DRYRUN:-0}" = 1 ]; then
            l=$(head -1 "$HV/lists/missing.txt"); [ -n "$l" ] && DRYRUN=1 "$HV/convert_v20.sh" $l; echo "DRYRUN: condor_submit_dag -f -batch-name charm_h5v20 -maxjobs $MAXJOBS h5v20.dag"; exit 0
          fi
          mkdir -p "$PARTS" "$Q/logs/h5_v20"
          cd "$HV" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_h5v20 -maxjobs "$MAXJOBS" h5v20.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_h5v20|Total for klyczak' || echo "(charm_h5v20: nic w kolejce)"
          tail -n 4 "$HV/h5v20.dag.dagman.out" 2>/dev/null
          hv_lists 2>/dev/null | tail -8
          ls "$HV"/h5v20.dag.rescue* 2>/dev/null && echo "UWAGA: plik rescue (nieudane węzły) — sprawdzić logs/h5_v20/*.err"
          [ -s "$H5V20/meta/summary.json" ] && head -c 1500 "$H5V20/meta/summary.json"; echo
          echo "DONE: $([ -s "$PARTS/DONE" ] && cat "$PARTS/DONE" || echo brak)" ;;
  *) echo "użycie: $0 check|submit|status"; exit 2 ;;
esac
