#!/bin/bash
# submit_psf.sh check|submit|status  (lxplus, ~/private/charm_q/Q_lxplus/pairsf) — projekt A, pojedyncze D0←D*, bez par.
set -uo pipefail
PSF=$(cd "$(dirname "$0")" && pwd); export PSF
source "$PSF/psf_common.sh"
cd "$PSF"; mkdir -p "$Q/logs/pairsf"
status() { condor_q -batch 2>/dev/null | grep -E "charm_pairsf|Total for klyczak" || true; ls -la "$PSF_OUT" 2>/dev/null | tail -5 || true; return 0; }
check() {
  local ok=0 q s
  condor_q -batch 2>/dev/null | grep -qE 'charm_pairsf' && { echo "STOP: pairsf już w kolejce"; ok=1; }
  [ -s "$PSF_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$PSF_OUT/DONE")); FORCE=1 nadpisuje"; ok=1; }
  for f in psf_analyze.py btrig_common.py src/charmtag/dstar_refit.py psf.sub psf_job.sh; do [ -f "$PSF/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  for s in data_428770 data_427929 fd_nd; do
    n=$(ls $(psf_glob $s) 2>/dev/null | wc -l); echo "  wejście $s: $n plik(ów) [$(psf_glob $s)]"; [ "$n" -gt 0 ] || { echo "STOP: brak wejścia $s"; ok=1; }
  done
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  PSF_IN=$PSF_IN PSF_OUT=$PSF_OUT"; return $ok
}
case ${1:-check} in
  status) status ;;
  check)  check && echo "OK — gotowe (./submit_psf.sh submit)" ;;
  submit) check || exit 1
          chmod +x "$PSF"/*.sh
          env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit psf.sub
          sleep 3; status ;;
  *) sed -n 2p "$0" ;;
esac
