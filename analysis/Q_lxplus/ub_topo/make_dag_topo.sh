#!/bin/bash
# cfg/ub_topo.yaml (= pf_pairs/cfg/pf_unblind.yaml + out_dir=$UT_OUT + pairs.diff.{var,channels}) i topo.dag: pairs → fit2d → report. NIE wysyła.
set -euo pipefail
UT=$(cd "$(dirname "$0")" && pwd); export UT
source "$UT/topo_common.sh"
[ -s "$PF/cfg/pf_unblind.yaml" ] || { echo "brak $PF/cfg/pf_unblind.yaml (bieg odsłonięcia)"; exit 1; }
mkdir -p "$UT/cfg" "$Q/logs/ubtopo"
python3 "$PF/tools/merge_cfg.py" "$UT/cfg/ub_topo.yaml" "$PF/cfg/pf_unblind.yaml" -s "paths.out_dir=$UT_OUT" -s "pairs.diff.var=$UT_VAR" \
  -s "pairs.diff.channels=[$UT_CHANNELS]" >/dev/null
if ! diff <(grep -vE 'out_dir|^\s+diff:|^\s+(var|channels|charge_source|edges|expect):' "$PF/cfg/pf_unblind.yaml") \
          <(grep -vE 'out_dir|^\s+diff:|^\s+(var|channels|charge_source|edges|expect):|^\s+- (D0D0|DpD0|DpDp)$' "$UT/cfg/ub_topo.yaml") >/dev/null; then
  echo "STOP: cfg/ub_topo.yaml różni się od pf_unblind.yaml nie tylko out_dir/pairs.diff"; exit 1; fi
V="penv=\"UT_OUT=$UT_OUT UT_VAR=$UT_VAR UT_CHANNELS=$UT_CHANNELS UT_MU_REF=$UT_MU_REF\""
S=$UT/topo.sub
{
  echo "# $(date -u +%FT%TZ): UT_OUT=$UT_OUT VAR=$UT_VAR CHANNELS=$UT_CHANNELS NORM=$NORM"
  echo "JOB t_pairs $S";  echo "VARS t_pairs step=\"pairs\" cpus=\"4\" mem=\"16GB\" flavour=\"workday\" $V"; echo "RETRY t_pairs 1"
  echo "JOB t_fit2d $S";  echo "VARS t_fit2d step=\"fit2d\" cpus=\"2\" mem=\"8GB\" flavour=\"workday\" $V"
  echo "JOB t_report $S"; echo "VARS t_report step=\"report\" cpus=\"1\" mem=\"2GB\" flavour=\"espresso\" $V"
  echo "PARENT t_pairs CHILD t_fit2d"; echo "PARENT t_fit2d CHILD t_report"
} > "$UT/topo.dag"
echo "topo.dag: 3 węzły; cfg sha256 $(sha256sum "$UT/cfg/ub_topo.yaml" | cut -c1-16)"
