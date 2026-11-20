#!/bin/bash
# Węzeł UB_NDSW: h5 czytane wprost z EOS (fuse), wynik → $NW_OUT/{NDSW.md, ndsw.json, ndsw_hists.npz, nw.log, DONE}
set -eo pipefail
NW=$(cd "$(dirname "$0")" && pwd); export NW
source "$NW/nw_common.sh"
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/nw_XXXXXX"); trap 'cd /; rm -rf "$WD"' EXIT; cd "$WD"
echo "[$(date -u +%FT%TZ)] UB_NDSW host $(hostname) NW_OUT=$NW_OUT"
nw_env
python3 -u "$NW/ndsw.py" --h5-428770 "$NW_H5_8/data_428770_*_none.h5" --h5-427929 "$NW_H5_7/data_427929_*_mb.h5" --trg "$NW_TRG/trg_*.npz" \
  --ref "$NW/ndstar_full_428770.json" --out "$WD/res" 2>&1 | tee "$WD/nw.log"
mkdir -p "$NW_OUT"
for f in NDSW.md ndsw.json ndsw_hists.npz; do put_output "$WD/res/$f" "$NW_OUT/$f"; done
put_output "$WD/nw.log" "$NW_OUT/nw.log"
date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$NW_OUT/DONE"
echo "[$(date -u +%FT%TZ)] OK UB_NDSW"
