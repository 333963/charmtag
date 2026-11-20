#!/bin/bash
# submit_nw.sh check|submit|status — UB_NDSW (1 job). Submit w czystym env (bez LCG). Stary wynik odkładany, nie kasowany.
set -uo pipefail
NW=$(cd "$(dirname "$0")" && pwd); export NW
source "$NW/nw_common.sh"
mkdir -p "$Q/logs/ubndsw"
check() { local ok=0 n8 n7 nt
  n8=$(ls "$NW_H5_8"/data_428770_*_none.h5 2>/dev/null | wc -l); n7=$(ls "$NW_H5_7"/data_427929_*_mb.h5 2>/dev/null | wc -l); nt=$(ls "$NW_TRG"/trg_*.npz 2>/dev/null | wc -l)
  echo "h5 428770 v1.1: $n8; h5 427929: $n7; trg: $nt"; { [ "$n8" -gt 0 ] && [ "$n7" -gt 0 ] && [ "$nt" -gt 0 ]; } || { echo "STOP: brak wejść"; ok=1; }
  [ -s "$NW/ndstar_full_428770.json" ] || { echo "STOP: brak ndstar_full_428770.json"; ok=1; }
  [ -s "$Q/ub_trig/ubtrig_run.py" ] || { echo "STOP: brak ub_trig/ubtrig_run.py"; ok=1; }
  condor_q -batch 2>/dev/null | grep -q charm_ubndsw && { echo "STOP: charm_ubndsw w kolejce"; ok=1; }
  echo "AFS: $(fs lq ~ | awk 'NR==2{print $4}')"; [ $ok = 0 ] && echo "OK — można wysłać"; return $ok; }
case ${1:-check} in
  check) check ;;
  submit) check || exit 1
    if [ -e "$NW_OUT/DONE" ]; then old="$NW_OUT/old_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$old"; for f in DONE NDSW.md ndsw.json ndsw_hists.npz nw.log; do [ -e "$NW_OUT/$f" ] && mv "$NW_OUT/$f" "$old/"; done; fi
    cd "$NW" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit nw.sub ;;
  status) condor_q -batch 2>/dev/null | grep charm_ubndsw; ls -la "$NW_OUT" 2>/dev/null ;;
esac
