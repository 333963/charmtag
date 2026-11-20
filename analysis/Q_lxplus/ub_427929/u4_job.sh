#!/bin/bash
# Węzeł UB427. u4_job.sh <krok> [grupa]
#   filt <g>        h5 v2.0 427929 (grupa g z U4_NF) → u4_filter.py → $U4_OUT/h5/*_mb.h5 + fstats/f_<g>.json
#   gn <próbka> <g> ROOT v2.0 → gn_root.py (mc: mb_nd; data: data_428770 / data_427929) → $U4_OUT/gn/<mc|d428|d427>/*.npz
#   cw              u4_cw.py → $U4_OUT/cw/{cw.json, mb_rates_427929.json, w2d_*.npy}
#   apply           charmtag apply D0/Dplus H (modele nominalne) na h5 *_mb.h5 → $U4_OUT/scores/{D0_H,Dplus_H}/data_none.h5
#   pairs           cfg 427929 → charmtag pairs --unblind --norm direct → fit2d --unblind → u4_report.py → $U4_OUT/{UB427.md,ub427.json,pairs/,DONE}
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
U4=$(cd "$(dirname "$0")" && pwd); export U4
source "$U4/u4_common.sh"
step=$1; a1=${2:--}; a2=${3:--}
set -u
export CHARM_UNBLIND=1
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/u4_${step}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] UB427 $step $a1 $a2 host $(hostname) WD=$WD U4_OUT=$U4_OUT"
[ "${DRYRUN:-0}" = 1 ] || u4_env
put() { [ "${DRYRUN:-0}" = 1 ] || put_output "$1" "$2"; }
case $step in
  filt)
    mapfile -t F < <(u4_h5 | u4_group "$a1" "$U4_NF"); [ ${#F[@]} -gt 0 ] || { echo "pusta grupa $a1"; exit 2; }
    mkdir -p "$WD/in"; L=(); for f in "${F[@]}"; do L+=("$(fetch_input "$f" "$WD/in")"); done
    run python3 -u "$U4/u4_filter.py" "${L[@]}" --outdir "$WD/out" --stats "$WD/f_$a1.json"
    if [ "${DRYRUN:-0}" != 1 ]; then for f in "$WD"/out/*_mb.h5; do put "$f" "$U4_OUT/h5/$(basename "$f")"; done; put "$WD/f_$a1.json" "$U4_OUT/fstats/f_$a1.json"; fi ;;
  gn)
    smp=$a1; g=$a2
    case $smp in mb_nd) kind=mc; n=$U4_NMC; sub=mc ;; data_428770) kind=data; n=$U4_N8; sub=d428 ;; data_427929) kind=data; n=$U4_N7; sub=d427 ;;
      *) echo "próbka $smp?"; exit 2 ;; esac
    mapfile -t F < <(u4_root "$smp" | u4_group "$g" "$n"); [ ${#F[@]} -gt 0 ] || { echo "pusta grupa $smp $g"; exit 2; }
    mkdir -p "$WD/in"
    for f in "${F[@]}"; do
      tag=$(job_tag "$f"); out=$U4_OUT/gn/$sub/${smp}_$tag.npz
      if [ "${FORCE:-0}" != 1 ] && eos_exists "$out"; then echo "pomijam $out"; continue; fi
      in=$f; [ "${DRYRUN:-0}" = 1 ] || in=$(fetch_input "$f" "$WD/in")
      run python3 -u "$U4/gn_root.py" --kind "$kind" --tag "$tag" "$in" "$WD/${smp}_$tag.npz"
      put "$WD/${smp}_$tag.npz" "$out"; rm -f "$WD/in/$(basename "$f")" "$WD/${smp}_$tag.npz"
    done ;;
  cw)
    [ "${DRYRUN:-0}" = 1 ] || for s in mc d428 d427; do mkdir -p "$WD/gn/$s"; cp "$U4_OUT"/gn/$s/*.npz "$WD/gn/$s/"; done
    run python3 -u "$U4/u4_cw.py" --mc "$WD/gn/mc" --d428 "$WD/gn/d428" --d427 "$WD/gn/d427" --nominal-mb "$U4_NOM_MB" --anchors "$U4_ANCHORS" --out "$WD/cw"
    if [ "${DRYRUN:-0}" != 1 ]; then for f in "$WD"/cw/*; do put "$f" "$U4_OUT/cw/$(basename "$f")"; done; fi ;;
  apply|pairs)
    python3 "$Q/pf_pairs/tools/merge_cfg.py" "$WD/cfg427.yaml" "$U4_NOM_CFG" -s "paths.out_dir=$WD/run" -s "paths.scores_dir=$WD/scores" \
      -s "paths.samples.data_none=$U4_OUT/h5/data_427929_*_mb.h5" -s "paths.sample_info={}" -s "pairs.mb_rates_json=$WD/mb_rates_427929.json" >/dev/null
    if [ "$step" = apply ]; then
      for sp in D0 Dplus; do run charmtag apply -c "$WD/cfg427.yaml" --species "$sp" --variant H --samples data_none
        put "$WD/scores/${sp}_H/data_none.h5" "$U4_OUT/scores/${sp}_H/data_none.h5"; put "$WD/scores/${sp}_H/samples_meta.json" "$U4_OUT/scores/${sp}_H/samples_meta.json"; done
    else
      [ "${DRYRUN:-0}" = 1 ] || cp "$U4_OUT/cw/mb_rates_427929.json" "$WD/"
      [ "${DRYRUN:-0}" = 1 ] || for sp in D0 Dplus; do mkdir -p "$WD/scores/${sp}_H"
        cp "$NOMINAL/scores/${sp}_H/signal_mc.h5" "$NOMINAL/scores/${sp}_H/incl_mc.h5" "$WD/scores/${sp}_H/"     # MC: nominał (te same progi/szablony)
        cp "$U4_OUT/scores/${sp}_H/data_none.h5" "$WD/scores/${sp}_H/"; done
      run charmtag pairs -c "$WD/cfg427.yaml" --unblind --norm "$NORM"
      run charmtag fit2d -c "$WD/cfg427.yaml" --unblind --ntoys "${U4_NTOYS:-0}"
      run python3 -u "$U4/u4_ndstar.py" "$U4_OUT/h5/data_427929_*_mb.h5" --out "$WD/ndstar_427929.json"
      run python3 -u "$U4/u4_report.py" --pairs "$WD/run/pairs" --nominal "$U4_NOM_PAIRS" --cw "$U4_OUT/cw/cw.json" --fstats "$U4_OUT/fstats" \
          --ndstar "$WD/ndstar_427929.json" --out "$WD/rep"
      if [ "${DRYRUN:-0}" != 1 ]; then
        for f in "$WD"/run/pairs/*.json "$WD"/run/pairs/*.md "$WD"/run/pairs/*.npz "$WD"/run/pairs/charmtag.log "$WD/cfg427.yaml"; do [ -f "$f" ] && put "$f" "$U4_OUT/pairs/$(basename "$f")"; done
        put "$WD/rep/UB427.md" "$U4_OUT/UB427.md"; put "$WD/rep/ub427.json" "$U4_OUT/ub427.json"; put "$WD/ndstar_427929.json" "$U4_OUT/ndstar_427929.json"
        date -u +%FT%TZ > "$WD/DONE"; put "$WD/DONE" "$U4_OUT/DONE"
      fi
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK UB427 $step $a1 $a2"
