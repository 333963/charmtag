#!/bin/bash
# Węzeł UBPNET427.  upn427_job.sh ext <g> → $UPN_OUT/parts/pnf_data_427929_<g>.npz (h5 427929 po filtrze UB427; okno sygnału otwarte)
#                   upn427_job.sh fit     → $UPN_OUT/{ubpnet.json,UBPNET427.md,ubpnet.log,DONE}
# Bez evset: h5 po filtrze zawierają tylko zdarzenia klasy MB (te same co pary UB427). Modele pnetd/xgbd z SBTRAIN v3 (tło 428770) bez retrenowania;
# test transferu: pary SB×SB 427929 vs części pasm 428770 ($SB3/parts). SF nóg: HYBSF/HYBDS dla data_427929; n_pair × f_ch (prawda mb_nd).
set -eo pipefail
UPN=$(cd "$(dirname "$0")" && pwd); export UPN
source "$UPN/upn427_common.sh"
step=$1; g=${2:--}
set -u
export CHARM_UNBLIND=1
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/upn4_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */upn4_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
echo "[$(date -u +%FT%TZ)] UBPNET427 $step $g host $(hostname) UPN_OUT=$UPN_OUT"
uj_gate || exit 3
case $step in
  ext)
    PN=$UPN_PN; export PN; source "$UPN_PN/pn_common.sh"; pn_env
    export PYTHONPATH="$UPN_PN:${PYTHONPATH:-}"
    mapfile -t ALL < <(upn_h5)
    FILES=("${ALL[@]:$((g * UPN_NG)):$UPN_NG}")
    [ ${#FILES[@]} -gt 0 ] || { echo "pusta grupa $g"; exit 2; }
    printf '%s\n' "${FILES[@]}" > "$WD/files_$g.txt"
    python3 -u "$UPN/upn_extract.py" --models "$PN_MODELS" --out "$WD/pnf_data_427929_$g.npz" "${FILES[@]}"
    upn_put "$WD/pnf_data_427929_$g.npz" "$UPN_OUT/parts/pnf_data_427929_$g.npz"
    upn_put "$WD/files_$g.txt" "$UPN_OUT/parts/files_$g.txt" ;;
  fit)
    uj_env
    export OMP_NUM_THREADS=$UPN_CPUS MKL_NUM_THREADS=$UPN_CPUS OPENBLAS_NUM_THREADS=1
    export PYTHONPATH="$UPN:$UPN_PN:$Q/ub_jk:${PYTHONPATH:-}" UPN_PN
    python3 -c 'import torch, xgboost, iminuit, charmtag, os; print("torch", torch.__version__, "xgb", xgboost.__version__, "charmtag", charmtag.__file__)
assert charmtag.__file__.startswith(os.environ["PF"] + "/src"), charmtag.__file__'
    n=$(ls "$UPN_OUT"/parts/pnf_data_427929_*.npz 2>/dev/null | wc -l); ng=$(upn_ngroups)
    [ "$n" -eq "$ng" ] || { echo "STOP: części danych $n / $ng"; exit 4; }
    mkdir -p "$WD/mc" "$WD/d" "$WD/md" "$WD/sb" "$WD/pairs"
    for f in "$P3"/parts/pn_*_nd_*.npz; do cp "$(fetch_input "$f" "$WD/mc")" "$WD/mc/" 2>/dev/null || true; done
    for f in "$UPN_OUT"/parts/pnf_data_427929_*.npz; do cp "$(fetch_input "$f" "$WD/d")" "$WD/d/" 2>/dev/null || true; done
    for f in "$SB3"/parts/pn_data_428770_*.npz; do cp "$(fetch_input "$f" "$WD/sb")" "$WD/sb/" 2>/dev/null || true; done
    for f in pn_models_sbtrain.pt xgbd_fold0.json xgbd_fold1.json xgbd_fold2.json xgbd_fold3.json xgbd_fold4.json; do cp "$(fetch_input "$SB3/train/$f" "$WD/md")" "$WD/md/" 2>/dev/null || true; done
    for f in pair_inputs.json pair_inputs.npz pairs_full.npz fit2d_unblind.json; do cp "$U427_PAIRS/$f" "$WD/pairs/"; done
    cp "$HY3/res/hybsf.json" "$WD/hybsf.json"; cp "$HY3/ds/hybds.json" "$WD/hybds.json"
    ls -la "$WD/mc" "$WD/d" "$WD/md" "$WD/pairs" | head -60; echo "SB 428770: $(ls "$WD/sb" | wc -l) części"
    python3 -u "$UPN/upn_fit.py" --mc "$WD/mc/pn_*_nd_*.npz" --data "$WD/d/pnf_data_427929_*.npz" \
      --modelsd "$WD/md/pn_models_sbtrain.pt" --xgbd "$WD/md" --hybsf "$WD/hybsf.json" --hybds "$WD/hybds.json" \
      --pairs-dir "$WD/pairs" --sample data_427929 --npair-scale "$UPN_FCH" --scale-now 1.0 --ref-sb "$WD/sb/pn_data_428770_*.npz" \
      --K "$UPN_K" --nboot "$UPN_NBOOT" --threads "$UPN_CPUS" --out "$WD/res" 2>&1 | tee "$WD/ubpnet.log"
    [ "${PIPESTATUS[0]}" -eq 0 ] || { upn_put "$WD/ubpnet.log" "$UPN_OUT/ubpnet.log"; echo "BŁĄD upn_fit"; exit 5; }
    cp "$WD/res/UBPNET.md" "$WD/res/UBPNET427.md"
    for f in ubpnet.json UBPNET427.md; do upn_put "$WD/res/$f" "$UPN_OUT/$f"; done
    upn_put "$WD/ubpnet.log" "$UPN_OUT/ubpnet.log"
    date -u +%FT%TZ > "$WD/DONE"; upn_put "$WD/DONE" "$UPN_OUT/DONE" ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK UBPNET427 $step $g"
