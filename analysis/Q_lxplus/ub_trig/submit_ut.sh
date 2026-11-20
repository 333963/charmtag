#!/bin/bash
# submit_ut.sh check | submit | status — DAG: trg × UT_NG → fit. Submit w czystym env (bez LCG), condor_submit_dag -f.
set -uo pipefail
UT=$(cd "$(dirname "$0")" && pwd); export UT
source "$UT/ut_common.sh"
mkdir -p "$Q/logs/ubtrig"
cmd=${1:-check}
check() {
  local ok=0 n
  uj_gate || ok=1
  n=$(ut_files | wc -l); echo "ROOT 428770 dpt3to150: $n plików"; [ "$n" -ge 90 ] || { echo "STOP: za mało ROOT 428770"; ok=1; }
  [ -s "$UJ_LEGS" ] || { echo "STOP: brak $UJ_LEGS"; ok=1; }
  [ -s "$BTD/bt_data_428770.npz" ] && [ -s "$BTD/bt_fd_nd.npz" ] || { echo "STOP: brak bt_*.npz w $BTD"; ok=1; }
  echo "AFS: $(fs lq ~ | awk 'NR==2{print $4}')"
  condor_q -batch 2>/dev/null | grep -c charm_ubtrig | sed 's/^/joby charm_ubtrig w kolejce: /'
  [ $ok = 0 ] && echo "OK — można wysłać"; return $ok
}
case $cmd in
  check) check ;;
  submit)
    check || exit 1
    if [ -e "$UT_OUT/DONE" ]; then                                 # poprzedni wynik nie kasowany, tylko odkładany
      old="$UT_OUT/old_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$old"
      for f in DONE UBTRIG.md ubtrig.json ubtrig.log; do [ -e "$UT_OUT/$f" ] && mv "$UT_OUT/$f" "$old/"; done; echo "stary wynik → $old"
    fi
    P="penv=\"UT_OUT=$UT_OUT UT_NG=$UT_NG OMP_NUM_THREADS=2\""
    {
      for ((g = 0; g < UT_NG; g++)); do [ -s "$UT_OUT/trg/trg_$g.npz" ] && continue; echo "JOB trg_$g $UT/ut.sub"; echo "VARS trg_$g step=\"trg\" g=\"$g\" cpus=\"1\" mem=\"4GB\" flavour=\"workday\" $P"; echo "RETRY trg_$g 1"; done
      echo "JOB fit $UT/ut.sub"; echo "VARS fit step=\"fit\" g=\"-\" cpus=\"2\" mem=\"8GB\" flavour=\"workday\" $P"; echo "RETRY fit 1"
      par=$(for ((g = 0; g < UT_NG; g++)); do [ -s "$UT_OUT/trg/trg_$g.npz" ] || printf 'trg_%d ' $g; done)
      [ -n "$par" ] && echo "PARENT ${par}CHILD fit"
    } > "$UT/ut.dag"
    cd "$UT" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_ubtrig ut.dag ;;
  nullsubmit)   # test zerowy: UT_NNULL ziaren (domyślnie 24) → nullagg; wymaga trg i prawdziwego ubtrig.json
    n=$(ls "$UT_OUT"/trg/trg_*.npz 2>/dev/null | wc -l); [ "$n" -eq "$UT_NG" ] || { echo "STOP: trg $n / $UT_NG"; exit 1; }
    [ -s "$UT_OUT/ubtrig.json" ] || { echo "STOP: brak $UT_OUT/ubtrig.json"; exit 1; }
    condor_q -batch 2>/dev/null | grep -q charm_ubtrnull && { echo "STOP: charm_ubtrnull już w kolejce"; exit 1; }
    if [ -e "$UT_OUT/null/NULL_DONE" ]; then old="$UT_OUT/null_old_$(date +%Y%m%d_%H%M%S)"; mv "$UT_OUT/null" "$old"; echo "stary test zerowy → $old"; fi
    NN=${UT_NNULL:-24}; P="penv=\"UT_OUT=$UT_OUT UT_NG=$UT_NG OMP_NUM_THREADS=2\""
    {
      for ((s = 1; s <= NN; s++)); do echo "JOB null_$s $UT/ut.sub"; echo "VARS null_$s step=\"null\" g=\"$s\" cpus=\"2\" mem=\"8GB\" flavour=\"workday\" bname=\"charm_ubtrnull\" $P"; echo "RETRY null_$s 1"; done
      echo "JOB nullagg $UT/ut.sub"; echo "VARS nullagg step=\"nullagg\" g=\"-\" cpus=\"1\" mem=\"2GB\" flavour=\"microcentury\" bname=\"charm_ubtrnull\" $P"
      echo "PARENT $(for ((s = 1; s <= NN; s++)); do printf 'null_%d ' $s; done)CHILD nullagg"
    } > "$UT/ut_null.dag"
    cd "$UT" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_ubtrnull ut_null.dag ;;
  status) condor_q -batch 2>/dev/null | grep charm_ubtrig; ls -la "$UT_OUT" 2>/dev/null | head; ls "$UT_OUT/trg" 2>/dev/null | wc -l ;;
esac
