#!/bin/bash
# submit_uc.sh check | submit | status — DAG: jeden węzeł fit. Submit w czystym env (bez LCG), condor_submit_dag -f.
set -uo pipefail
UC=$(cd "$(dirname "$0")" && pwd); export UC
source "$UC/uc_common.sh"
mkdir -p "$Q/logs/ubcomb"
cmd=${1:-check}
check() {
  local ok=0 n
  uj_gate || ok=1
  [ -s "$UJ_LEGS" ] || { echo "STOP: brak $UJ_LEGS"; ok=1; }
  n=$(ls "$UC_TRG"/trg_*.npz 2>/dev/null | wc -l); echo "trg 428770: $n plików"; [ "$n" -ge 1 ] || { echo "STOP: brak $UC_TRG/trg_*.npz (UBTRIG)"; ok=1; }
  [ -s "$BTD/bt_data_428770.npz" ] && [ -s "$BTD/bt_fd_nd.npz" ] || { echo "STOP: brak bt_*.npz w $BTD"; ok=1; }
  for f in pairs/pair_inputs.json pairs/pair_inputs.npz pairs/pairs_full.npz ndstar_427929.json; do [ -s "$U4_OUT/$f" ] || { echo "STOP: brak $U4_OUT/$f (UB427)"; ok=1; }; done
  [ -s "$U4_OUT/jk_ev200/ubjk.json" ] && echo "  jackknife 427929: jest" || echo "  jackknife 427929: brak (f_jk = 1,18 założone)"
  echo "AFS: $(fs lq ~ | awk 'NR==2{print $4}')"
  condor_q -batch 2>/dev/null | grep -q charm_ubcomb && { echo "STOP: charm_ubcomb już w kolejce"; ok=1; }
  [ $ok = 0 ] && echo "OK — można wysłać"; return $ok
}
case $cmd in
  check) check ;;
  submit)
    check || exit 1
    if [ -e "$UC_OUT/DONE" ]; then                                 # poprzedni wynik nie kasowany, tylko odkładany
      old="$UC_OUT/old_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$old"
      for f in DONE UBCOMB.md ubcomb.json ubcomb.log; do [ -e "$UC_OUT/$f" ] && mv "$UC_OUT/$f" "$old/"; done; echo "stary wynik → $old"
    fi
    P="penv=\"UC_OUT=$UC_OUT OMP_NUM_THREADS=2\""
    { echo "JOB fit $UC/uc.sub"; echo "VARS fit step=\"fit\" cpus=\"2\" mem=\"12GB\" flavour=\"workday\" $P"; echo "RETRY fit 1"; } > "$UC/uc.dag"
    cd "$UC" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_ubcomb uc.dag ;;
  status) condor_q -batch 2>/dev/null | grep charm_ubcomb; ls -la "$UC_OUT" 2>/dev/null | head ;;
esac
