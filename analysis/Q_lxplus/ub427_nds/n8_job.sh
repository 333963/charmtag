#!/bin/bash
# Węzeł UB427_NDS: n8_job.sh fit → $N8_OUT/{NDS_427929.md, ndstar_variants_427929.json, dstar_hist_427929.npz, DONE}
set -eo pipefail
N8=$(cd "$(dirname "$0")" && pwd); export N8
source "$N8/n8_common.sh"
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/n8_XXXXXX"); trap 'cd /; rm -rf "$WD"' EXIT; cd "$WD"
echo "[$(date -u +%FT%TZ)] UB427_NDS host $(hostname) N8_OUT=$N8_OUT"
n8_env
n=$(ls "$U4_OUT"/h5/data_427929_*_mb.h5 2>/dev/null | wc -l); [ "$n" -gt 0 ] || { echo "STOP: brak h5 427929 w $U4_OUT/h5"; exit 3; }
mkdir -p "$WD/h5"; cp "$U4_OUT"/h5/data_427929_*_mb.h5 "$WD/h5/"
python3 -u "$N8/n8_ndstar.py" "$WD/h5/*.h5" --ref "$N8/ndstar_full_428770.json" --out "$WD/res" 2>&1 | tee "$WD/n8.log"
mkdir -p "$N8_OUT"
for f in NDS_427929.md ndstar_variants_427929.json dstar_hist_427929.npz; do put_output "$WD/res/$f" "$N8_OUT/$f"; done
put_output "$WD/n8.log" "$N8_OUT/n8.log"
date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$N8_OUT/DONE"
echo "[$(date -u +%FT%TZ)] OK UB427_NDS"
