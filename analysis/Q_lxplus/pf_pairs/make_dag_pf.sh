#!/bin/bash
# Generuje cfg/pf.yaml, cfg/pf_mix.yaml i pf.dag — TYLKO węzły zależne od par, z poprawką PF (src = AA_merge + PF):
#   pairs → fit2d (ślepo) → toy1..NTOYJOBS (×NTOY) → toys_merge                      [T]
#   pairs → syst core_{D0D0,DpD0,DpDp}, files, wpb_<6> → wp ; norm (POST /bin/true) ; report   [U]
#   mix_pairs (pairs.mix_control) → mix_fit2d                                          [S5 nom_mix]
#   compare (po wszystkim) → $PF_OUT/compare_pf.{md,json}, DONE
# Bez treningu: modele/score/validation z $NOMINAL. Użycie: NTOYJOBS=50 NTOY=200 ./make_dag_pf.sh (NOSYST=1 / NOMIX=1 / NOTOYS=1)
set -euo pipefail
PF=$(cd "$(dirname "$0")" && pwd); export PF
source "$PF/pf_common.sh"
NTOYJOBS=${NTOYJOBS:-50}; NTOY=${NTOY:-200}
V=$PF/src/config/variants_bkg; mkdir -p "$PF/cfg" "$Q/logs/pf"
COMMON=(-s "paths.out_dir=$PF_OUT" -s "paths.models_dir=$NOMINAL/models" -s "paths.scores_dir=$NOMINAL/scores" -s "pairs.track_id=d0z0" -s "pairs.lumi_fraction=1.0")
python3 "$PF/tools/merge_cfg.py" "$PF/cfg/pf.yaml" "$BASE_CFG" "${COMMON[@]}" >/dev/null
python3 "$PF/tools/merge_cfg.py" "$PF/cfg/pf_mix.yaml" "$BASE_CFG" "$V/s5_common.yaml" "${COMMON[@]}" -s "paths.out_dir=$PF_OUT/mix" >/dev/null
S=$PF/pf.sub
node() { # node NAZWA step arg arg2 cpus mem flavour
  echo "JOB $1 $S"; echo "VARS $1 step=\"$2\" arg=\"$3\" arg2=\"$4\" cpus=\"$5\" mem=\"$6\" flavour=\"$7\""
}
{
  echo "# wygenerowane $(date -u +%FT%TZ): PF_OUT=$PF_OUT NOMINAL=$NOMINAL NORM=$NORM NTOYJOBS=$NTOYJOBS NTOY=$NTOY NOSYST=${NOSYST:-0} NOMIX=${NOMIX:-0} NOTOYS=${NOTOYS:-0}"
  all=""
  node pairs pairs - - 4 16GB workday; echo "RETRY pairs 1"
  node fit2d fit2d - - 2 8GB workday; echo "PARENT pairs CHILD fit2d"; all="$all fit2d"
  if [ "${NOTOYS:-0}" != 1 ]; then
    kids=""
    for i in $(seq 1 "$NTOYJOBS"); do node toy$i toy "$i" "$NTOY" 1 4GB tomorrow; echo "RETRY toy$i 2"; kids="$kids toy$i"; done
    echo "PARENT fit2d CHILD$kids"
    node toys_merge toys_merge - - 1 2GB espresso; echo "PARENT$kids CHILD toys_merge"; all="$all toys_merge"
  fi
  if [ "${NOSYST:-0}" != 1 ]; then
    for ch in D0D0 DpD0 DpDp; do node core_$ch syst core "$ch" 1 8GB workday; echo "PARENT pairs CHILD core_$ch"; echo "RETRY core_$ch 1"; all="$all core_$ch"; done
    node files syst files - 1 16GB workday; echo "PARENT pairs CHILD files"; echo "RETRY files 1"; all="$all files"
    wpk=""
    for v in wp_tight wp_loose arb85 arb95 tagdm07 tagdm15; do
      node wpb_$v syst wp-build "$v" 4 32GB workday; echo "PARENT pairs CHILD wpb_$v"; echo "RETRY wpb_$v 1"; wpk="$wpk wpb_$v"
    done
    node wp syst wp - 1 8GB workday; echo "PARENT$wpk CHILD wp"; all="$all wp"
    node norm syst norm - 1 2GB espresso; echo "PARENT pairs CHILD norm"; echo "SCRIPT POST norm /bin/true"   # brak truth A3/tune nie blokuje
    node report syst report - 1 2GB espresso; echo "PARENT$all norm CHILD report"; echo "SCRIPT POST report /bin/true"
    all="$all report"
  fi
  if [ "${NOMIX:-0}" != 1 ]; then
    node mix_pairs mix_pairs - - 4 16GB workday; echo "RETRY mix_pairs 1"
    node mix_fit2d mix_fit2d - - 2 8GB workday; echo "PARENT mix_pairs CHILD mix_fit2d"; all="$all mix_fit2d"
  fi
  node compare compare - - 1 2GB espresso; echo "PARENT$all CHILD compare"
  echo "PRIORITY pairs 10"
} > "$PF/pf.dag"
echo "zapisano $PF/pf.dag ($(grep -c '^JOB' "$PF/pf.dag") węzłów); cfg: $(ls "$PF/cfg" | tr '\n' ' ')"
