#!/bin/bash
# Jeden węzeł DAG PF2. Użycie: pf_job.sh <krok> [arg] [arg2]
#   pairs | fit2d | toy <seed> <ntoy> | toys_merge | syst <część> [arg] | mix_pairs | mix_fit2d | compare
# Zawsze ŚLEPO (bez --unblind). Modele/score nominalne z $NOMINAL (paths.models_dir/scores_dir), wyjścia $PF_OUT (EOS).
# KWOTA AFS: praca w mktemp w scratch, trap rm -rf; na AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
PF=$(cd "$(dirname "$0")" && pwd); export PF
source "$PF/pf_common.sh"
step=$1; a1=${2:-}; a2=${3:-}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/pf_${step}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
CFG=$PF/cfg/pf.yaml; MIXCFG=$PF/cfg/pf_mix.yaml
echo "[$(date -u +%FT%TZ)] pf $step $a1 $a2 host $(hostname) WD=$WD"
[ "${DRYRUN:-0}" = 1 ] || pf_env
case $step in
  pairs)      run charmtag pairs -c "$CFG" --norm "$NORM" ;;
  fit2d)      run charmtag fit2d -c "$CFG" --ntoys 0 ;;
  toy)        run python3 -u "$Q/charmtag_jobs/toys_job.py" -c "$CFG" --seed "$a1" --ntoy "$a2" ;;
  toys_merge) run python3 -u "$Q/charmtag_jobs/toys_merge.py" -c "$CFG" ;;
  syst)       case $a1 in
                core)     run charmtag syst -c "$CFG" core --channels "$a2" ;;
                files)    run charmtag syst -c "$CFG" files ;;
                wp-build) run charmtag syst -c "$CFG" wp-build --variant "$a2" --norm "$NORM" ;;
                wp)       run charmtag syst -c "$CFG" wp ;;
                norm)     run charmtag syst -c "$CFG" norm --truth-dir "$TRUTH" ;;
                report)   run charmtag syst -c "$CFG" report --truth-dir "$TRUTH" ;;
                *) echo "nieznana część syst $a1"; exit 2 ;;
              esac ;;
  mix_pairs)  run charmtag pairs -c "$MIXCFG" --norm "$NORM" ;;
  mix_fit2d)  run charmtag fit2d -c "$MIXCFG" --ntoys 0 ;;
  compare)    run python3 -u "$PF/tools/compare_pf.py" --pf "$PF_OUT" --nominal "$NOMINAL" --out "$WD/compare_pf"
              if [ "${DRYRUN:-0}" != 1 ]; then
                put_output "$WD/compare_pf.md" "$PF_OUT/compare_pf.md"; put_output "$WD/compare_pf.json" "$PF_OUT/compare_pf.json"
                rm -rf "$PF_OUT/cache" "$PF_OUT/mix/cache" 2>/dev/null || true
                date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$PF_OUT/DONE"
              fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK $step $a1 $a2"
