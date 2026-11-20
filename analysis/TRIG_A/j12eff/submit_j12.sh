#!/bin/bash
# submit_j12.sh check|submit|status   (na lxplus, w ~/private/charm_q/Q_lxplus/j12eff)
set -uo pipefail
J12=$(cd "$(dirname "$0")" && pwd); export J12
source "$J12/j12_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s q r f
  for f in j12_extract.py j12_merge.py tj_kernel.py trg_root_extract.py btrig_common.py; do [ -f "$J12/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  for s in $J12_SAMPLES; do r=$(ls "$RAWV20/$s" 2>/dev/null | grep -c 'dpt3to150.*\.root$'); echo "  $s: ROOT $r, węzłów $(j12_np "$s")"; [ "$r" -gt 0 ] || { echo "STOP: brak ROOT $s"; ok=1; }; done
  condor_q -batch 2>/dev/null | grep -qE 'charm_j12eff' && { echo "STOP: J12EFF w kolejce"; ok=1; }
  [ -s "$J12_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone (FORCE=1 nadpisuje)"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  J12_OUT=$J12_OUT"; return $ok
}
mkdag() {
  local s k np nodes="" PENV
  PENV="J12_OUT=$J12_OUT RAWV20=$RAWV20 J12_SAMPLES=$(echo $J12_SAMPLES | tr ' ' ,) J12_NPARTS=$(echo $J12_NPARTS | tr ' ' ,)"
  { for s in $J12_SAMPLES; do np=$(j12_np "$s"); for ((k = 0; k < np; k++)); do
      echo "JOB ext_${s}_$k $J12/j12.sub"; echo "VARS ext_${s}_$k node=\"ext_${s}_$k\" args=\"ext $s $k\" mem=\"4GB\" flavour=\"workday\" penv=\"$PENV\""; echo "RETRY ext_${s}_$k 2"; nodes="$nodes ext_${s}_$k"
    done; done
    echo "JOB merge $J12/j12.sub"; echo "VARS merge node=\"merge\" args=\"merge\" mem=\"4GB\" flavour=\"longlunch\" penv=\"$PENV\""; echo "RETRY merge 1"
    echo "PARENT$nodes CHILD merge"; } > "$J12/j12.dag"
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_j12.sh submit)" ;;
  submit) check || exit 1
          mkdir -p "$Q/logs/j12eff"; chmod +x "$J12"/*.sh "$J12"/*.py
          ( j12_env && cd "$J12" && python3 tj_kernel.py >/dev/null 2>&1; python3 -c 'import sys; sys.path.insert(0, "'"$J12"'"); import j12_extract, j12_merge' ) || { echo "STOP: import"; exit 1; }
          mkdag; cd "$J12" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_j12eff j12.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_j12eff|Total for klyczak' || true; echo "części: $(ls "$J12_OUT/parts" 2>/dev/null | wc -l)"; return 0 2>/dev/null || exit 0 ;;
  *) echo "użycie: $0 check|submit|status"; exit 2 ;;
esac
