#!/bin/bash
# Węzeł J12EFF: j12_job.sh ext <próbka> <część> | merge. Scratch mktemp + trap; wyjścia tylko EOS.
set -eo pipefail
J12=$(cd "$(dirname "$0")" && pwd); export J12
[ -n "${J12_SAMPLES:-}" ] && J12_SAMPLES=${J12_SAMPLES//,/ }
[ -n "${J12_NPARTS:-}" ] && J12_NPARTS=${J12_NPARTS//,/ }
source "$J12/j12_common.sh"
step=$1; s=${2:-}; part=${3:-0}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/j12_${step}_${s}_${part}_XXXXXX")
trap 'cd /; rm -rf "$WD"' EXIT
cd "$WD"; export TMPDIR=$WD
echo "[$(date -u +%FT%TZ)] J12EFF $step $s $part host $(hostname)"
j12_env; mkdir -p "$J12_OUT/parts"
case $step in
  ext)
    np=$(j12_np "$s"); o=j12_${s}_p$part.npz
    ls "$RAWV20/$s"/*dpt3to150*.root | sort | awk -v p="$part" -v n="$np" '(NR-1)%n==p' > files.txt
    echo "plików: $(wc -l < files.txt)"; [ -s files.txt ]
    # shellcheck disable=SC2046
    python3 -u "$J12/j12_extract.py" --out "$WD/$o" $(cat files.txt)
    put_output "$WD/$o" "$J12_OUT/parts/$o" ;;
  merge)
    python3 -u "$J12/j12_merge.py" --parts "$J12_OUT/parts" --out "$WD/res"
    for f in "$WD"/res/*; do put_output "$f" "$J12_OUT/$(basename "$f")"; done
    date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$J12_OUT/DONE" ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK J12EFF $step $s $part"
