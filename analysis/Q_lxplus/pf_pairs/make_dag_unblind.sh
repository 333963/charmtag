#!/bin/bash
# Generuje cfg/pf_unblind.yaml (= cfg/pf.yaml z out_dir=$PF_OUT/unblind) i unblind.dag: pairs(--unblind) → fit2d(--unblind, toys) → report.
# NIE wysyła. Wymaga istniejącego cfg/pf.yaml (z biegu pf) — ten sam config, żadnych innych zmian.
set -euo pipefail
PF=$(cd "$(dirname "$0")" && pwd); export PF
source "$PF/pf_common.sh"
UB_OUT=${UB_OUT:-$PF_OUT/unblind}
[ -s "$PF/cfg/pf.yaml" ] || { echo "brak cfg/pf.yaml — najpierw bieg pf (make_dag_pf.sh)"; exit 1; }
python3 "$PF/tools/merge_cfg.py" "$PF/cfg/pf_unblind.yaml" "$PF/cfg/pf.yaml" -s "paths.out_dir=$UB_OUT" >/dev/null
# różnica względem cfg/pf.yaml może dotyczyć TYLKO paths.out_dir
if ! diff <(grep -v 'out_dir' "$PF/cfg/pf.yaml") <(grep -v 'out_dir' "$PF/cfg/pf_unblind.yaml") >/dev/null; then
  echo "STOP: cfg/pf_unblind.yaml różni się od cfg/pf.yaml nie tylko out_dir"; exit 1; fi
S=$PF/unblind.sub
{
  echo "# wygenerowane $(date -u +%FT%TZ): UB_OUT=$UB_OUT NORM=$NORM UB_NTOYS=${UB_NTOYS:-100}  (odsłonięcie — tylko po decyzji)"
  echo "JOB ub_pairs $S";  echo "VARS ub_pairs step=\"pairs\" cpus=\"4\" mem=\"16GB\" flavour=\"workday\""; echo "RETRY ub_pairs 1"
  echo "JOB ub_fit2d $S";  echo "VARS ub_fit2d step=\"fit2d\" cpus=\"2\" mem=\"8GB\" flavour=\"workday\""
  echo "JOB ub_report $S"; echo "VARS ub_report step=\"report\" cpus=\"1\" mem=\"2GB\" flavour=\"espresso\""
  echo "PARENT ub_pairs CHILD ub_fit2d"; echo "PARENT ub_fit2d CHILD ub_report"
} > "$PF/unblind.dag"
echo "unblind.dag: 3 węzły; cfg/pf_unblind.yaml (out_dir=$UB_OUT); sha256 cfg: $(sha256sum "$PF/cfg/pf_unblind.yaml" | cut -c1-64)"
