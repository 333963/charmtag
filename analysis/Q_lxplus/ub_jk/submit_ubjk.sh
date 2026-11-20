#!/bin/bash
# submit_ubjk.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/ub_jk)
#  check  — nic z charm_ubjk w kolejce, kwota AFS < 85 %, flaga $UB_OUT/UNBLIND_OK + hash cfg, $UB_OUT/DONE, wejścia odsłonięcia na EOS,
#           brak $UJ_OUT/DONE, kod ($PF/src z unblinded_jackknife) i ubjk_run.py
#  submit — check → make_dag_ubjk.sh → test importów (podpowłoka z env LCG) → condor_submit_dag w CZYSTYM env (batch charm_ubjk)
#  status / fetch — kolejka + pliki $UJ_OUT
set -uo pipefail
UJ=$(cd "$(dirname "$0")" && pwd); export UJ
source "$UJ/ubjk_common.sh"
cmd=${1:-check}
check() {
  local ok=0 q f
  condor_q -batch 2>/dev/null | grep -qE 'charm_ubjk|ubjk.dag' && { echo "STOP: UBJK w kolejce"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  uj_gate || ok=1
  for f in pair_inputs.json pair_inputs.npz pairs_full.npz fit2d_unblind.json; do
    [ -s "$UJ_PAIRS/$f" ] || { echo "STOP: brak $UJ_PAIRS/$f"; ok=1; }
  done
  if [ -s "$UJ_LEGS" ]; then echo "  nogi: $UJ_LEGS ($(du -h "$UJ_LEGS" | cut -f1)) → prep lekki"
  else echo "UWAGA: brak $UJ_LEGS → prep przebuduje pary w scratch (cmd_pairs --unblind, 16 GB, ~40 min; kontrola Σ == pairs_full.npz)"; fi
  [ -s "$UJ_FILES_JSON" ] && echo "  tabela plików (kontrola): $UJ_FILES_JSON" || echo "UWAGA: brak $UJ_FILES_JSON (tabela tylko z h5)"
  [ -s "$UJ_OUT/DONE" ] && { echo "STOP: UBJK już policzone ($(cat "$UJ_OUT/DONE")) — usuń $UJ_OUT/DONE tylko za zgodą"; ok=1; }
  [ -f "$UJ/ubjk_run.py" ] || { echo "STOP: brak ubjk_run.py"; ok=1; }
  grep -q 'def unblinded_jackknife' "$PF/src/charmtag/syst/files.py" 2>/dev/null || { echo "STOP: $PF/src bez syst/files.unblinded_jackknife"; ok=1; }
  grep -q '"track_id": "d0z0"' "$PF/src/charmtag/pairs/build.py" 2>/dev/null || { echo "STOP: $PF/src bez poprawki PF (d0z0)"; ok=1; }
  echo "kwota AFS: ${q:-?}%  GROUPS=$UJ_GROUPS UB_OUT=$UB_OUT UJ_OUT=$UJ_OUT NCHUNK=$UJ_NCHUNK MODES=$UJ_MODES REF=$UJ_REF NORM=$NORM"
  echo "  src pairs sha256 teraz: $(cat $(ls "$PF"/src/charmtag/pairs/*.py | sort) | sha256sum | cut -c1-64) (przy odsłonięciu: $(awk '/^pairs_src_sha256/{print $2}' "$UB_OUT/UNBLIND_OK" 2>/dev/null)) — różnica dopuszczalna, prep sprawdza odtworzenie μ"
  return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_ubjk.sh submit)" ;;
  submit) check || exit 1
          "$UJ/make_dag_ubjk.sh" || exit 1
          chmod +x "$UJ"/*.sh
          if [ "${DRYRUN:-0}" = 1 ]; then
            DRYRUN=1 "$UJ/ubjk_job.sh" prep; for ((i = 0; i < UJ_NCHUNK; i++)); do DRYRUN=1 "$UJ/ubjk_job.sh" fit "$i"; done
            DRYRUN=1 "$UJ/ubjk_job.sh" ref; DRYRUN=1 "$UJ/ubjk_job.sh" finish; exit 0
          fi
          ( uj_env && CHARM_UNBLIND=1 python3 "$UJ/ubjk_run.py" --help >/dev/null && echo "import OK" ) || exit 1   # podpowłoka: LCG nie przecieka do dagmana
          cd "$UJ" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_ubjk ubjk.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_ubjk|Total for klyczak'; tail -n 3 "$UJ/ubjk.dag.dagman.out" 2>/dev/null | cut -c1-140
          ls -la "$UJ_OUT" 2>/dev/null | head -20; [ -f "$UJ_OUT/ubjk.md" ] && head -40 "$UJ_OUT/ubjk.md" ;;
  fetch)  find "$UJ_OUT" -maxdepth 2 \( -name '*.md' -o -name '*.json' -o -name DONE -o -name '*.log' \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac
