#!/bin/bash
# Węzeł UBCOMB: uc_job.sh fit → $UC_OUT/{ubcomb.json,UBCOMB.md,ubcomb.log,DONE}
set -eo pipefail
UC=$(cd "$(dirname "$0")" && pwd); export UC
source "$UC/uc_common.sh"
step=${1:-fit}
set -u
export CHARM_UNBLIND=1
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/uc_${step}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
echo "[$(date -u +%FT%TZ)] UBCOMB $step host $(hostname) UC_OUT=$UC_OUT"
uj_gate || exit 3
uj_env
export PYTHONPATH="$UC:$Q/ub_jk:${PYTHONPATH:-}"
python3 -c "import dstar_refit_aa, charmtag, os; print('dstar_refit:', dstar_refit_aa.__file__); print('charmtag:', charmtag.__file__)"
case $step in
  fit)
    cp "$UC_TRG"/trg_*.npz "$WD/"; cp "$BTD/bt_data_428770.npz" "$BTD/bt_fd_nd.npz" "$WD/"
    mkdir -p "$WD/p427"; cp "$U4_OUT"/pairs/{pair_inputs.json,pair_inputs.npz,pairs_full.npz} "$WD/p427/"
    JK=(); [ -s "$U4_OUT/jk_ev200/ubjk.json" ] && JK=(--jk427 "$U4_OUT/jk_ev200/ubjk.json")
    python3 -u "$UC/ubcomb_run.py" --p428 "$UJ_PAIRS" --legs428 "$UJ_LEGS" --trg428 "$WD/trg_*.npz" --bt428 "$WD/bt_data_428770.npz" \
      --bt-mc "$WD/bt_fd_nd.npz" --p427 "$WD/p427" --ndstar427 "$U4_OUT/ndstar_427929.json" "${JK[@]}" --variants "$UC_VARIANTS" \
      --out "$WD/res" 2>&1 | tee "$WD/ubcomb.log"
    for f in ubcomb.json UBCOMB.md; do uc_put "$WD/res/$f" "$UC_OUT/$f"; done
    uc_put "$WD/ubcomb.log" "$UC_OUT/ubcomb.log"
    date -u +%FT%TZ > "$WD/DONE"; uc_put "$WD/DONE" "$UC_OUT/DONE" ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK UBCOMB $step"
