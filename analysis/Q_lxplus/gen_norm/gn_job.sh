#!/bin/bash
# Węzeł gen_norm. Użycie: gn_job.sh <step> <kind> <tag> <src>
#  root mc|data <tag> <plik.root>  → gn_root.py → $GN_OUT/parts/<kind>/<próbka>_<tag>.npz (istnieje → pomija, FORCE=1 liczy ponownie)
#  legs - incl_mc -                → gn_legs.py (charmtag legs_for_truth + mc_rates; scores nominalne) → $GN_OUT/legs/
#  combine - all -                 → gn_combine.py → $GN_OUT/{w2d.npy, acceptance.json, effmap_direct.npz, …}, run_info.json, DONE
# Scratch mktemp + trap (kasowany); wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
GN=$(cd "$(dirname "$0")" && pwd); export GN
source "$GN/gn_common.sh"
step=$1; kind=${2:--}; tag=${3:--}; src=${4:--}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/gn_${step}_${tag}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] gen_norm $step $kind $tag host $(hostname) WD=$WD GN_OUT=$GN_OUT"
[ "${DRYRUN:-0}" = 1 ] || gn_env
case $step in
  root)
    smp=mb_nd; [ "$kind" = data ] && smp=data_428770
    out=$GN_OUT/parts/$kind/${smp}_$tag.npz
    if [ "${FORCE:-0}" != 1 ] && eos_exists "$out"; then echo "pomijam (istnieje) $out"; exit 0; fi
    in=$src; [ "${DRYRUN:-0}" = 1 ] || in=$(fetch_input "$src" "$WD")
    run python3 -u "$GN/gn_root.py" --kind "$kind" --tag "$tag" ${MAX_EVENTS:+--max-events $MAX_EVENTS} "$in" "$WD/${smp}_$tag.npz"
    if [ "${DRYRUN:-0}" != 1 ]; then put_output "$WD/${smp}_$tag.npz" "$out"; fi ;;   # kopia wejścia w $WD → trap
  legs)
    run python3 -u "$GN/gn_legs.py" --cfg "$BASE_CFG" --out-dir "$NOMINAL" --out "$WD/legs" --ref-rates "$REF_RATES"
    if [ "${DRYRUN:-0}" != 1 ]; then for f in "$WD"/legs/*; do put_output "$f" "$GN_OUT/legs/$(basename "$f")"; done; fi ;;
  combine)
    run python3 -u "$GN/gn_combine.py" --mc "$GN_OUT/parts/mc" --data "$GN_OUT/parts/data" --legs "$GN_OUT/legs" --out "$WD/out" \
        --data-ref-tags "$DATA_REF_TAGS" ${REF_ANCHORS:+--ref-anchors "$REF_ANCHORS"} --fmt "$GN_FMT" ${GN_MB_EXCLUDE:+--mb-exclude "$GN_MB_EXCLUDE"}
    if [ "${DRYRUN:-0}" != 1 ]; then
      python3 - "$WD/out/run_info.json" <<'PY'
import json, os, sys, time
json.dump({"finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "gn_out": os.environ["GN_OUT"], "base_cfg": os.environ["BASE_CFG"],
           "nominal": os.environ["NOMINAL"], "data_ref_tags": os.environ["DATA_REF_TAGS"], "gn_fmt": os.environ.get("GN_FMT"),
           "ptprod": os.environ.get("PTPROD"), "raw": os.environ.get("GN_RAW"), "ref_anchors": os.environ.get("REF_ANCHORS"), "mb_exclude": os.environ.get("GN_MB_EXCLUDE"), "max_events": os.environ.get("MAX_EVENTS"),
           "src_snapshot": open(os.path.join(os.environ["GN"], "src", "SNAPSHOT.txt")).read()}, open(sys.argv[1], "w"), indent=1)
PY
      for f in "$WD"/out/*; do put_output "$f" "$GN_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$GN_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK gen_norm $step $kind $tag"
