#!/bin/bash
# submit_pz.sh check | submit | status — DAG PREZ-2: ext × 4 (D0/D+ × fd/mb) → train × 4 (gatunek × xgb/net) → report. Submit w czystym env (bez LCG).
# Części ext/train już na EOS są pomijane (FORCE=1 → wszystko od nowa). Stary wynik (DONE/report) odkładany do old_<data>, nic nie kasowane.
set -uo pipefail
PZ=$(cd "$(dirname "$0")" && pwd); export PZ
source "$PZ/pz_common.sh"
mkdir -p "$Q/logs/preznet"
cmd=${1:-check}
check() {
  local ok=0 sp f
  [ -s "$PZ_CFG" ] || { echo "STOP: brak configu nominału $PZ_CFG"; ok=1; }
  for sp in D0 Dplus; do
    for f in meta.json fold0.json fold4.json; do [ -s "$PZ_MODELS/${sp}_H/$f" ] || { echo "STOP: brak $PZ_MODELS/${sp}_H/$f"; ok=1; }; done
  done
  for f in fd_nd mb_nd; do n=$(ls /eos/user/k/klyczak/charm_h5_v11/parts/$f/${f}_*_none.h5 2>/dev/null | wc -l); echo "  h5 v1.1 $f: $n plików"; [ "$n" -gt 0 ] || ok=1; done
  [ -f "$PZ/src/charmtag/samples.py" ] || { echo "STOP: brak src/charmtag"; ok=1; }
  ( pz_env ) | tail -1 || ok=1
  echo "AFS: $(fs lq ~ | awk 'NR==2{print $4}')"
  condor_q -batch 2>/dev/null | grep -q charm_preznet && { echo "STOP: charm_preznet już w kolejce"; ok=1; }
  [ $ok = 0 ] && echo "OK — można wysłać"; return $ok
}
have() { [ "${FORCE:-0}" != 1 ] && [ -s "$1" ]; }
case $cmd in
  check) check ;;
  submit)
    check || exit 1
    if [ -e "$PZ_OUT/DONE" ]; then old="$PZ_OUT/old_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$old"; for f in DONE report; do [ -e "$PZ_OUT/$f" ] && mv "$PZ_OUT/$f" "$old/"; done; echo "stary wynik → $old"; fi
    P="penv=\"PZ_OUT=$PZ_OUT PZ_CPUS=8 PZ_EPOCHS=$PZ_EPOCHS OMP_NUM_THREADS=8\""
    {
      E=(); T=()
      for sp in D0 Dplus; do
        for x in fd mb; do
          have "$PZ_OUT/ext/${sp}_$x.npz" && continue
          echo "JOB ext_${sp}_$x $PZ/pz.sub"; echo "VARS ext_${sp}_$x step=\"ext\" sp=\"$sp\" x=\"$x\" cpus=\"1\" mem=\"16GB\" flavour=\"workday\" $P"; echo "RETRY ext_${sp}_$x 1"
          E+=("ext_${sp}_$x")
        done
        for x in xgb net; do
          have "$PZ_OUT/train_$sp/scores_$x.npz" && continue
          echo "JOB train_${sp}_$x $PZ/pz.sub"; fl=workday; [ "$x" = net ] && fl=tomorrow                       # sieć: 2 × 5 foldów × ≤ 30 epok na ~2,5 M kandydatów → do ~10 h CPU
          echo "VARS train_${sp}_$x step=\"train\" sp=\"$sp\" x=\"$x\" cpus=\"8\" mem=\"24GB\" flavour=\"$fl\" $P"; echo "RETRY train_${sp}_$x 1"
          T+=("train_${sp}_$x")
          for y in fd mb; do printf '%s\n' "${E[@]}" | grep -qx "ext_${sp}_$y" && echo "PARENT ext_${sp}_$y CHILD train_${sp}_$x"; done
        done
      done
      echo "JOB report $PZ/pz.sub"; echo "VARS report step=\"report\" sp=\"-\" x=\"-\" cpus=\"1\" mem=\"16GB\" flavour=\"longlunch\" $P"; echo "RETRY report 1"
      [ ${#T[@]} -gt 0 ] && echo "PARENT ${T[*]} CHILD report"
    } > "$PZ/pz.dag"
    cd "$PZ" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_preznet pz.dag ;;
  status) condor_q -batch 2>/dev/null | grep charm_preznet; ls -la "$PZ_OUT" "$PZ_OUT/ext" 2>/dev/null | head -20 ;;
esac
