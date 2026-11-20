#!/bin/bash
# Węzeł PREZ-2. pz_job.sh ext <gatunek> <fd|mb>   → $PZ_OUT/ext/<gatunek>_<fd|mb>.npz
#               pz_job.sh train <gatunek> <xgb|net> → $PZ_OUT/train_<gatunek>/{scores,info,imp}_<rodzina>.*
#               pz_job.sh report                  → $PZ_OUT/report/{*.png,*.pdf,TABLE.md,table.json}, $PZ_OUT/DONE
set -eo pipefail
PZ=$(cd "$(dirname "$0")" && pwd); export PZ
source "$PZ/pz_common.sh"
step=$1; sp=${2:--}; x=${3:--}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/pz_${step}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
echo "[$(date -u +%FT%TZ)] PREZ-2 $step $sp $x host $(hostname) PZ_OUT=$PZ_OUT"
pz_env
case $step in
  ext)
    cap=$PZ_CAP_SIG; [ "$x" = mb ] && { v=PZ_CAP_BKG_$sp; cap=${!v}; }
    python3 -u "$PZ/pz_ext.py" --cfg "$PZ_CFG" --models "$PZ_MODELS/${sp}_H" --species "$sp" --sample "$x" --cap "$cap" --workdir "$WD" \
      --out "$WD/${sp}_$x.npz"
    pz_put "$WD/${sp}_$x.npz" "$PZ_OUT/ext/${sp}_$x.npz" ;;
  train)
    cp "$PZ_OUT/ext/${sp}_fd.npz" "$PZ_OUT/ext/${sp}_mb.npz" "$WD/"
    python3 -u "$PZ/pz_train.py" --sig "$WD/${sp}_fd.npz" --bkg "$WD/${sp}_mb.npz" --family "$x" --out "$WD/res" --threads "$PZ_CPUS" --epochs "$PZ_EPOCHS" \
      2>&1 | tee "$WD/train_$x.log"
    for f in "$WD"/res/* "$WD/train_$x.log"; do pz_put "$f" "$PZ_OUT/train_$sp/$(basename "$f")"; done ;;
  report)
    for s in D0 Dplus; do mkdir -p "$WD/in/$s"; cp "$PZ_OUT/train_$s"/*.npz "$PZ_OUT/train_$s"/*.json "$WD/in/$s/"; done
    python3 -u "$PZ/pz_report.py" --in "$WD/in/D0" --in "$WD/in/Dplus" --out "$WD/report" 2>&1 | tee "$WD/report.log"
    cp "$PZ/model.md" "$WD/report/" 2>/dev/null || true
    for f in "$WD"/report/* "$WD/report.log"; do pz_put "$f" "$PZ_OUT/report/$(basename "$f")"; done
    date -u +%FT%TZ > "$WD/DONE"; pz_put "$WD/DONE" "$PZ_OUT/DONE" ;;
  *) echo "nieznany krok $step"; exit 1 ;;
esac
echo "[$(date -u +%FT%TZ)] OK PREZ-2 $step $sp $x"
