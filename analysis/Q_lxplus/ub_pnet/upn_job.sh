#!/bin/bash
# Węzeł UBPNET.  upn_job.sh ext <g>  → $UPN_OUT/parts/pnf_data_428770_<g>.npz (pliki h5 v2.0 nr g·UPN_NG … ; okno sygnału otwarte)
#                upn_job.sh evset    → $UPN_OUT/evset.npz (event_number próbki nominalnej, h5 v1.1 data_none)
#                upn_job.sh fit      → $UPN_OUT/{ubpnet.json,UBPNET.md,ubpnet.log,DONE}
set -eo pipefail
UPN=$(cd "$(dirname "$0")" && pwd); export UPN
source "$UPN/upn_common.sh"
step=$1; g=${2:--}
set -u
export CHARM_UNBLIND=1
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/upn_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */upn_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
echo "[$(date -u +%FT%TZ)] UBPNET $step $g host $(hostname) UPN_OUT=$UPN_OUT"
uj_gate || exit 3
case $step in
  ext)
    # środowisko jak dext łańcucha v3 (pn_env = q_env, bez $PF/src)
    PN=$UPN_PN; export PN; source "$UPN_PN/pn_common.sh"; pn_env
    export PYTHONPATH="$UPN_PN:${PYTHONPATH:-}"
    mapfile -t ALL < <(upn_h5)
    FILES=("${ALL[@]:$((g * UPN_NG)):$UPN_NG}")
    [ ${#FILES[@]} -gt 0 ] || { echo "pusta grupa $g"; exit 2; }
    printf '%s\n' "${FILES[@]}" > "$WD/files_$g.txt"
    python3 -u "$UPN/upn_extract.py" --models "$PN_MODELS" --out "$WD/pnf_data_428770_$g.npz" "${FILES[@]}"
    upn_put "$WD/pnf_data_428770_$g.npz" "$UPN_OUT/parts/pnf_data_428770_$g.npz"
    upn_put "$WD/files_$g.txt" "$UPN_OUT/parts/files_$g.txt" ;;
  evset)
    q_env
    mapfile -t V11 < <(upn_v11)
    [ ${#V11[@]} -eq 39 ] || { echo "STOP: h5 v1.1 data_none: ${#V11[@]} plików (oczekiwane 39)"; exit 4; }
    python3 -u "$UPN/upn_evset.py" --out "$WD/evset.npz" "${V11[@]}"
    upn_put "$WD/evset.npz" "$UPN_OUT/evset.npz" ;;
  fit)
    uj_env                                                    # LCG + venv + $PF/src (charmtag snapshotu odsłonięcia — kontrola nominału)
    export OMP_NUM_THREADS=$UPN_CPUS MKL_NUM_THREADS=$UPN_CPUS OPENBLAS_NUM_THREADS=1
    export PYTHONPATH="$UPN:$UPN_PN:$Q/ub_jk:${PYTHONPATH:-}" UPN_PN
    python3 -c 'import torch, xgboost, iminuit, charmtag, os; print("torch", torch.__version__, "xgb", xgboost.__version__, "charmtag", charmtag.__file__)
assert charmtag.__file__.startswith(os.environ["PF"] + "/src"), charmtag.__file__'
    n=$(ls "$UPN_OUT"/parts/pnf_data_428770_*.npz 2>/dev/null | wc -l); ng=$(upn_ngroups)
    [ "$n" -eq "$ng" ] || { echo "STOP: części danych $n / $ng"; exit 4; }
    [ -s "$UPN_OUT/evset.npz" ] || { echo "STOP: brak evset.npz"; exit 4; }
    mkdir -p "$WD/mc" "$WD/d" "$WD/md"
    for f in "$P3"/parts/pn_*_nd_*.npz; do cp "$(fetch_input "$f" "$WD/mc")" "$WD/mc/" 2>/dev/null || true; done
    for f in "$UPN_OUT"/parts/pnf_data_428770_*.npz "$UPN_OUT/evset.npz"; do cp "$(fetch_input "$f" "$WD/d")" "$WD/d/" 2>/dev/null || true; done
    for f in pn_models_sbtrain.pt xgbd_fold0.json xgbd_fold1.json xgbd_fold2.json xgbd_fold3.json xgbd_fold4.json; do cp "$(fetch_input "$SB3/train/$f" "$WD/md")" "$WD/md/" 2>/dev/null || true; done
    cp "$HY3/res/hybsf.json" "$WD/hybsf.json"; cp "$HY3/ds/hybds.json" "$WD/hybds.json"
    ls -la "$WD/mc" "$WD/d" "$WD/md" | head -60
    python3 -u "$UPN/upn_fit.py" --mc "$WD/mc/pn_*_nd_*.npz" --data "$WD/d/pnf_data_428770_*.npz" --evset "$WD/d/evset.npz" \
      --modelsd "$WD/md/pn_models_sbtrain.pt" --xgbd "$WD/md" --hybsf "$WD/hybsf.json" --hybds "$WD/hybds.json" \
      --pairs-dir "$UJ_PAIRS" --legs "$UJ_LEGS" --K "$UPN_K" --nboot "$UPN_NBOOT" --threads "$UPN_CPUS" --out "$WD/res" 2>&1 | tee "$WD/ubpnet.log"
    [ "${PIPESTATUS[0]}" -eq 0 ] || { upn_put "$WD/ubpnet.log" "$UPN_OUT/ubpnet.log"; echo "BŁĄD upn_fit"; exit 5; }
    for f in ubpnet.json UBPNET.md; do upn_put "$WD/res/$f" "$UPN_OUT/$f"; done
    upn_put "$WD/ubpnet.log" "$UPN_OUT/ubpnet.log"
    date -u +%FT%TZ > "$WD/DONE"; upn_put "$WD/DONE" "$UPN_OUT/DONE" ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK UBPNET $step $g"
