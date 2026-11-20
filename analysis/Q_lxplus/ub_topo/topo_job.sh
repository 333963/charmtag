#!/bin/bash
# Węzeł UB_TOPO: topo_job.sh pairs|fit2d|report. Scratch mktemp + trap; wyjścia tylko EOS ($UT_OUT); AFS tylko logi. DRYRUN=1 → komendy.
set -eo pipefail
UT=$(cd "$(dirname "$0")" && pwd); export UT
source "$UT/topo_common.sh"
step=$1
set -u
CFG=$UT/cfg/ub_topo.yaml
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/ut_${step}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD CHARM_UNBLIND=1
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] UB_TOPO $step host $(hostname) UT_OUT=$UT_OUT VAR=$UT_VAR"
[ -f "$CFG" ] || { echo "brak $CFG — make_dag_topo.sh"; exit 2; }
if [ "${DRYRUN:-0}" != 1 ]; then
  [ -s "$UB_OUT/UNBLIND_OK" ] && [ -s "$UB_OUT/DONE" ] || { echo "STOP: brak $UB_OUT/UNBLIND_OK lub DONE"; exit 3; }
  pd_env; mkdir -p "$UT_OUT/pairs"
fi
case $step in
  pairs)  run charmtag pairs -c "$CFG" --unblind --norm "$NORM" --diff "$UT_VAR" ;;
  fit2d)  [ "${DRYRUN:-0}" = 1 ] || cp -f "$UB_OUT/pairs/fit2d_unblind.json" "$UT_OUT/pairs/fit2d_unblind.json" 2>/dev/null || true
          run charmtag fit2d -c "$CFG" --unblind --diff --ntoys 0 ;;
  report) run python3 -u "$UT/../pd_diff/tools/diff_fill.py" --diff "$UT_OUT" --unblind "$UB_OUT" --out "$WD/tfin_diff"
          run python3 -u "$UT/topo_report.py" --pairs "$UT_OUT/pairs" --gen "$UT/ref/topo_gen.json" --out "$WD/res" --mu-ref "$UT_MU_REF" \
            --forecast "σ_λ ≈ 0,85; czysty DPS vs Pythia ~1,2σ; nadwyżka-DPS (λ=0,33) vs Pythia ~0,5–0,75σ (ref/forecast.json)"
          if [ "${DRYRUN:-0}" != 1 ]; then
            for f in tfin_diff.md tfin_diff.json res/UBTOPO.md res/ubtopo.json; do put_output "$WD/$f" "$UT_OUT/$(basename "$f")"; done
            rm -rf "$UT_OUT/cache" 2>/dev/null || true
            date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$UT_OUT/DONE"
          fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK UB_TOPO $step"
