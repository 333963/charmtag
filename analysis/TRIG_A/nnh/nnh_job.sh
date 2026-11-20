#!/bin/bash
# Jedyny węzeł nnh: trening sieci i XGB (OOF) + metryki MC + ε_B/SF na danych (pojedyncze D0←D*) → $NNH_OUT/{NNH.md,nnh.json,nnh_scores_*.npz,DONE}.
set -eo pipefail
NNH=$(cd "$(dirname "$0")" && pwd); Q=$(cd "$NNH/.." && pwd)
source "$Q/common.sh"
: "${NNH_IN:=$EOS_OUT/charm_q/b_trig}" "${NNH_OUT:=$EOS_OUT/charm_q/nnh}" "${NNH_NPROC:=8}"
g() { if [ -s "$NNH_IN/bt_$1.npz" ]; then echo "$NNH_IN/bt_$1.npz"; else echo "$NNH_IN/parts/bt_$1_*.npz"; fi; }
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/nnh_XXXXXX"); trap 'cd /; rm -rf "$WD"' EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
echo "[$(date -u +%FT%TZ)] NNH host $(hostname) IN=$NNH_IN OUT=$NNH_OUT"
q_env
export PYTHONPATH="$NNH:$NNH/src:${PYTHONPATH:-}" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=$NNH_NPROC
python3 -c "import torch, xgboost, charmtag.dstar_refit as d, btrig_common; assert hasattr(d, 'fit_passfail'); print('torch', torch.__version__)"
python3 -u "$NNH/nnh_train.py" --fd "$(g fd_nd)" --mb "$(g mb_nd)" --data "$(g data_428770)@428770" --data "$(g data_427929)@427929" \
    --threads "$NNH_NPROC" --out "$WD/res"
mkdir -p "$NNH_OUT"
for f in "$WD"/res/*; do put_output "$f" "$NNH_OUT/$(basename "$f")"; done
date -u +%FT%TZ > DONE; put_output DONE "$NNH_OUT/DONE"
echo "[$(date -u +%FT%TZ)] OK NNH"
