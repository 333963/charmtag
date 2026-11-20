#!/bin/bash
# submit_u4.sh check | submit | status — DAG UB427:
#   filt × U4_NF ─────────────────────────────→ apply ─┐
#   gn mb_nd × U4_NMC, gn 428770 × U4_N8, gn 427929 × U4_N7 → cw ─┴→ pairs (pairs --unblind, fit2d --unblind, raport) → DONE
# Submit w czystym env (LCG psuje condor_dagman), condor_submit_dag -f. DRYRUN=1 → tylko DAG.
set -uo pipefail
U4=$(cd "$(dirname "$0")" && pwd); export U4
source "$U4/u4_common.sh"
mkdir -p "$Q/logs/ub427"
cmd=${1:-check}
check() {
  local ok=0 n
  n=$(u4_h5 | wc -l); echo "h5 v2.0 data_427929 *_none: $n"; [ "$n" -ge 400 ] || { echo "STOP: za mało h5 427929"; ok=1; }
  for s in mb_nd data_428770 data_427929; do n=$(u4_root $s | wc -l); echo "ROOT v2.0 $s (do gn_root): $n"; [ "$n" -gt 0 ] || { echo "STOP: brak ROOT $s"; ok=1; }; done
  [ -s "$U4_NOM_CFG" ] || { echo "STOP: brak $U4_NOM_CFG"; ok=1; }
  [ -s "$U4_NOM_MB" ] || { echo "STOP: brak $U4_NOM_MB"; ok=1; }
  for f in "$U4_NOM_PAIRS/pair_inputs.json" "$U4_NOM_PAIRS/pairs_full.npz" "$NOMINAL/scores/D0_H/signal_mc.h5" "$NOMINAL/scores/Dplus_H/incl_mc.h5" \
           "$NOMINAL/models/D0_H/meta.json" ${U4_ANCHORS//,/ }; do [ -s "$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -e "$U4_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: $U4_OUT/DONE istnieje (FORCE=1 nadpisuje)"; ok=1; }
  echo "AFS: $(fs lq ~ 2>/dev/null | awk 'NR==2{print $4}')"
  echo "joby charm_ub427 w kolejce: $(condor_q -batch 2>/dev/null | grep -c charm_ub427)"
  [ $ok = 0 ] && echo "OK — można wysłać"; return $ok
}
dag() {
  local P="penv=\"U4_OUT=$U4_OUT NORM=$NORM OMP_NUM_THREADS=2\"" g kids_f="" kids_g=""
  for ((g = 0; g < U4_NF; g++)); do echo "JOB filt_$g $U4/u4.sub"; echo "VARS filt_$g step=\"filt\" a1=\"$g\" a2=\"-\" cpus=\"1\" mem=\"4GB\" disk=\"8GB\" flavour=\"longlunch\" $P"; echo "RETRY filt_$g 1"; kids_f="$kids_f filt_$g"; done
  for spec in "mb_nd $U4_NMC 4GB workday" "data_428770 $U4_N8 3GB workday" "data_427929 $U4_N7 3GB workday"; do
    set -- $spec
    for ((g = 0; g < $2; g++)); do echo "JOB gn_${1}_$g $U4/u4.sub"; echo "VARS gn_${1}_$g step=\"gn\" a1=\"$1\" a2=\"$g\" cpus=\"1\" mem=\"$3\" disk=\"6GB\" flavour=\"$4\" $P"; echo "RETRY gn_${1}_$g 1"; kids_g="$kids_g gn_${1}_$g"; done
  done
  echo "JOB cw $U4/u4.sub";    echo "VARS cw step=\"cw\" a1=\"-\" a2=\"-\" cpus=\"1\" mem=\"16GB\" disk=\"10GB\" flavour=\"longlunch\" $P"; echo "RETRY cw 1"
  echo "JOB apply $U4/u4.sub"; echo "VARS apply step=\"apply\" a1=\"-\" a2=\"-\" cpus=\"4\" mem=\"8GB\" disk=\"10GB\" flavour=\"workday\" $P"; echo "RETRY apply 1"
  echo "JOB pairs $U4/u4.sub"; echo "VARS pairs step=\"pairs\" a1=\"-\" a2=\"-\" cpus=\"4\" mem=\"16GB\" disk=\"20GB\" flavour=\"workday\" $P"; echo "RETRY pairs 1"
  echo "PARENT$kids_f CHILD apply"; echo "PARENT$kids_g CHILD cw"; echo "PARENT apply cw CHILD pairs"
}
case $cmd in
  check) check ;;
  submit)
    check || exit 1
    dag > "$U4/u4.dag"; echo "u4.dag: $(grep -c '^JOB' "$U4/u4.dag") węzłów"
    [ "${DRYRUN:-0}" = 1 ] && exit 0
    cd "$U4" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_ub427 u4.dag ;;
  status) condor_q -batch 2>/dev/null | grep charm_ub427; ls "$U4_OUT" 2>/dev/null; for s in h5 gn/mc gn/d428 gn/d427; do echo "$s: $(ls "$U4_OUT/$s" 2>/dev/null | wc -l)"; done ;;
  dag) dag ;;
  *) echo "użycie: $0 check|submit|status|dag"; exit 2 ;;
esac
