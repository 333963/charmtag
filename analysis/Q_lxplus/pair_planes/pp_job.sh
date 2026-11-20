#!/bin/bash
# Węzeł PAIR_PLANES. Użycie: pp_job.sh <próbka> <część> <nczęści> | pp_job.sh finish. Scratch mktemp + trap; wyjścia tylko EOS.
set -eo pipefail
PX=$(cd "$(dirname "$0")" && pwd); export PX
[ -n "${PP_SAMPLES:-}" ] && PP_SAMPLES=${PP_SAMPLES//,/ }
[ -n "${PP_NFILES:-}" ] && PP_NFILES=${PP_NFILES//,/ }
source "$PX/pp_common.sh"
step=$1; part=${2:-0}; np=${3:-1}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/pp_${step}_${part}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] PAIR_PLANES $step/$part/$np host $(hostname) PP_OUT=$PP_OUT"
if [ "${DRYRUN:-0}" != 1 ]; then pp_env; mkdir -p "$PP_OUT/parts" 2>/dev/null || true; fi
case $step in
  finish)
    run python3 -u "$PX/pp_merge.py" --parts "$PP_OUT/parts" --out "$WD" --samples "$PP_SAMPLES"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/*.npz "$WD"/*.meta.json "$WD"/PP_SUMMARY.json; do pp_put "$f" "$PP_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; pp_put "$WD/DONE" "$PP_OUT/DONE"
    fi ;;
  data_*|mb_*|fd_*)
    o=$step.p$part
    run python3 -u "$PX/pp_extract.py" --sample "$step" --h5dir "$H5V20/$step" --part "$part" --nparts "$np" --out "$WD/$o.npz" --wd "$WD"
    if [ "${DRYRUN:-0}" != 1 ]; then for x in npz meta.json; do pp_put "$WD/$o.$x" "$PP_OUT/parts/$o.$x"; done; fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK PAIR_PLANES $step/$part"
