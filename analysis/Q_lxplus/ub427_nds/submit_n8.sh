#!/bin/bash
# submit_n8.sh check|submit|status — UB427_NDS (1 job). Submit w czystym env (bez LCG). Stary wynik odkładany, nie kasowany.
set -uo pipefail
N8=$(cd "$(dirname "$0")" && pwd); export N8
source "$N8/n8_common.sh"
mkdir -p "$Q/logs/ub427nds"
check() { local ok=0 n
  n=$(ls "$U4_OUT"/h5/data_427929_*_mb.h5 2>/dev/null | wc -l); echo "h5 427929 (filtr UB427): $n"; [ "$n" -gt 0 ] || { echo "STOP: brak h5"; ok=1; }
  [ -s "$N8/ndstar_full_428770.json" ] || { echo "STOP: brak ndstar_full_428770.json"; ok=1; }
  condor_q -batch 2>/dev/null | grep -q charm_ub427nds && { echo "STOP: charm_ub427nds w kolejce"; ok=1; }
  echo "AFS: $(fs lq ~ | awk 'NR==2{print $4}')"; [ $ok = 0 ] && echo "OK — można wysłać"; return $ok; }
case ${1:-check} in
  check) check ;;
  submit) check || exit 1
    if [ -e "$N8_OUT/DONE" ]; then old="$N8_OUT/old_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$old"; for f in DONE NDS_427929.md ndstar_variants_427929.json n8.log; do [ -e "$N8_OUT/$f" ] && mv "$N8_OUT/$f" "$old/"; done; fi
    cd "$N8" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit n8.sub ;;
  status) condor_q -batch 2>/dev/null | grep charm_ub427nds; ls -la "$N8_OUT" 2>/dev/null ;;
esac
