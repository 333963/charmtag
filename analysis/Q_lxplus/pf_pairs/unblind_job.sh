#!/bin/bash
# Węzeł ODSŁONIĘCIA (T-PREP). Użycie: unblind_job.sh <krok>   krok: pairs | fit2d | report
# Uruchamiany WYŁĄCZNIE przez unblind.dag po decyzji o odsłonięciu (plik-flaga UNBLIND_OK + wpis w unblinding_checklist.md).
# Wejścia: ten sam cfg co pf (cfg/pf.yaml), ale paths.out_dir = $PF_OUT/unblind (ślepe wyjścia pf nietknięte).
# Scratch mktemp + trap rm -rf; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
PF=$(cd "$(dirname "$0")" && pwd); export PF
source "$PF/pf_common.sh"
step=$1
set -u
UB_OUT=${UB_OUT:-$PF_OUT/unblind}
CFG=$PF/cfg/pf_unblind.yaml
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/ub_${step}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] UNBLIND $step host $(hostname) WD=$WD UB_OUT=$UB_OUT"
[ -f "$CFG" ] || { echo "brak $CFG — uruchom make_dag_unblind.sh"; exit 2; }
# twarda blokada: bez flagi decyzji na EOS nic nie liczymy (flaga = kopia UNBLIND_OK z Maca + hash configu)
if [ "${DRYRUN:-0}" != 1 ]; then
  [ -s "$UB_OUT/UNBLIND_OK" ] || { echo "STOP: brak $UB_OUT/UNBLIND_OK (decyzja o odsłonięciu) — nie odsłaniam"; exit 3; }
  want=$(awk '/^cfg_sha256/{print $2}' "$UB_OUT/UNBLIND_OK"); have=$(sha256sum "$CFG" | cut -c1-64)
  [ -z "$want" ] || [ "$want" = "$have" ] || { echo "STOP: cfg zmieniony po decyzji ($have ≠ $want)"; exit 3; }
  pf_env
fi
case $step in
  pairs)  run charmtag pairs -c "$CFG" --unblind --norm "$NORM" ;;
  fit2d)  run charmtag fit2d -c "$CFG" --unblind --ntoys "${UB_NTOYS:-100}" ;;
  report) run python3 -u "$PF/tools/tfin_fill.py" --pf "$PF_OUT" --unblind "$UB_OUT" --template "$PF/tools/tfin_template.md" --out "$WD/tfin"
          if [ "${DRYRUN:-0}" != 1 ]; then
            put_output "$WD/tfin.md" "$UB_OUT/tfin.md"; put_output "$WD/tfin.json" "$UB_OUT/tfin.json"
            date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$UB_OUT/DONE"
          fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK UNBLIND $step"
