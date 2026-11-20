#!/bin/bash
# 1 proces = 1 plik ROOT v2.0 → h5 (schemat v1.2 + kolumny zdarzeniowe v2.0) → $H5V20/parts/<klucz>/<klucz>_<NNNNNN>_<presel>.h5
# Użycie: convert_v20.sh <klucz> <plik.root> <NNNNNN> <presele: none[,single_tag]>
# Zmienne: FORCE=1 (nadpisz), MAX_EVENTS=N (test), DRYRUN=1, OUTDIR (zamiast EOS). Scratch w mktemp + trap (AFS/węzeł czysty).
set -eo pipefail
HV=$(cd "$(dirname "$0")" && pwd)
source "$HV/h5v20_common.sh"
key=$1; src=$2; tag=$3; presels=${4:-none}; presels=${presels//+/,}   # w listach DAG presele łączone „+” (przecinek dzieli queue from)
set -u
OUTDIR=${OUTDIR:-$PARTS/$key}
WORK=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/h5v20_${key}_${tag}_XXXX")
trap 'rm -rf "$WORK"' EXIT
export PYTHONDONTWRITEBYTECODE=1
echo "[$(date -u +%FT%TZ)] H5v20 $key $tag presele $presels host $(hostname) src $src"
todo=()
for p in ${presels//,/ }; do
  out="$OUTDIR/${key}_${tag}_${p}.h5"
  if [ "${FORCE:-0}" != 1 ] && eos_exists "$out"; then echo "  pomijam (istnieje) $out"; else todo+=("$p"); fi
done
[ ${#todo[@]} -eq 0 ] && { echo "nic do zrobienia"; exit 0; }
if [ "${DRYRUN:-0}" = 1 ]; then
  for p in "${todo[@]}"; do echo "DRYRUN: python3 $HV/charmpp_h5_v20.py -o $WORK/${key}_${tag}_${p}.h5 --presel $p --step '$STEP' $src → $OUTDIR/"; done; exit 0
fi
q_env
python3 -c "import uproot, awkward, h5py" || { echo "BŁĄD: brak uproot/awkward/h5py w venv"; exit 2; }
in=$(fetch_input "$src" "$WORK"); echo "  wejście: $in"
mx=(); [ -n "${MAX_EVENTS:-}" ] && mx=(--max-events "$MAX_EVENTS")
for p in "${todo[@]}"; do
  loc="$WORK/${key}_${tag}_${p}.h5"
  python3 -u "$HV/charmpp_h5_v20.py" -o "$loc" --presel "$p" --step "$STEP" "${mx[@]}" "$in"
  python3 - "$loc" "$src" "$tag" "$key" <<'PY'
import sys, h5py
p, src, tag, key = sys.argv[1:]
with h5py.File(p, "a") as f:
    f.attrs["source_root"] = src; f.attrs["job"] = tag; f.attrs["sample"] = key; f.attrs["root_format"] = "v2.0 (CharmAnalysis charm-prod-v2.0, Oct28)"
    n = f["events"].shape[0] if "events" in f else 0
    print(f"  {p}: {n} wierszy, {len(f['events'].dtype.names) if 'events' in f else 0} kolumn, n_events_read={f.attrs.get('n_events_read')}, evtable={f.attrs.get('n_evtable')}")
PY
  put_output "$loc" "$OUTDIR/${key}_${tag}_${p}.h5"
  echo "  zapisano $OUTDIR/${key}_${tag}_${p}.h5 ($(du -h "$loc" | cut -f1))"
  rm -f "$loc"
done
[ "$in" != "$src" ] && rm -f "$in"
echo "[$(date -u +%FT%TZ)] OK"
