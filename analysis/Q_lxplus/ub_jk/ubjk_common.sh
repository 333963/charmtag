# ubjk_common.sh — paczka UBJK (jackknife μ po plikach na ODSŁONIĘTEJ próbce par; 12.11; source, nie uruchamiać).
# Kod: $Q/pf_pairs/src (ten sam snapshot co odsłonięcie/PD — NIE kopiujemy src), cfg: $Q/pf_pairs/cfg/pf_unblind.yaml.
# Wymaga ../pf_pairs/pf_common.sh (→ ../common.sh: EOS_OUT, q_env, put_output, fetch_input).
UJ=${UJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$UJ/.." && pwd)}
PF=$Q/pf_pairs
# shellcheck disable=SC1091
source "$PF/pf_common.sh"                                 # PF_OUT, NORM, pf_env (PYTHONPATH=$PF/src + asercje)
: "${UB_OUT:=$PF_OUT/unblind}"                           # odsłonięcie 05.11 (pairs/pair_inputs.*, pairs_full.npz, fit2d_unblind.json, UNBLIND_OK, DONE)
: "${UJ_GROUPS:=files}"                                  # files (39 plików, nominalnie) | events:G (G grup zdarzeń z hasha — opcja, lepsza precyzja σ)
case $UJ_GROUPS in files) _ujd=jk ;; events:*) _ujd=jk_ev${UJ_GROUPS#events:} ;; *) _ujd=jk_x ;; esac
: "${UJ_OUT:=$UB_OUT/$_ujd}"                             # wyjścia (EOS): perfile.*, fit_<i>.json, ref_*.json, ubjk.{json,md}, DONE
: "${UJ_NCHUNK:=8}"                                      # porcje replik (40 replik = nominał + 39 plików; ~5 replik × 2 tryby na węzeł)
: "${UJ_MODES:=fixed,rederive}"                          # fixed = jak unblinded_jackknife; rederive = też dataleg/SB×SB bez pliku k
: "${UJ_REF:=1}"                                         # 1 → węzeł ref (wprost charmtag.syst.files.unblinded_jackknife; tylko groups=files)
[ "$UJ_GROUPS" = files ] || UJ_REF=0
UJ_CFG=$PF/cfg/pf_unblind.yaml
UJ_PAIRS=$UB_OUT/pairs
UJ_LEGS=$UJ_PAIRS/legs_data.parquet                       # zapisany przez `pairs --unblind` (gdy pyarrow) — inaczej przebudowa w prep
UJ_FILES_JSON=$PF_OUT/pairs/syst/files.json               # ślepy `syst files` (tabela 39 plików) — kontrola/fallback
export UJ Q PF UB_OUT UJ_GROUPS UJ_OUT UJ_NCHUNK UJ_MODES UJ_REF UJ_CFG UJ_PAIRS UJ_LEGS UJ_FILES_JSON
uj_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else mkdir -p "$(dirname "$2")" && cp -f "$1" "$2"; fi; }
# blokada odsłonięcia: flaga decyzji + hash configu (jak unblind_job.sh) + odsłonięcie zakończone
uj_gate() {
  [ -s "$UB_OUT/UNBLIND_OK" ] || { echo "STOP: brak $UB_OUT/UNBLIND_OK (decyzja o odsłonięciu)"; return 3; }
  [ -s "$UB_OUT/DONE" ] || { echo "STOP: brak $UB_OUT/DONE (odsłonięcie niezakończone)"; return 3; }
  [ -s "$UJ_CFG" ] || { echo "STOP: brak $UJ_CFG"; return 3; }
  local want have
  want=$(awk '/^cfg_sha256/{print $2}' "$UB_OUT/UNBLIND_OK"); have=$(sha256sum "$UJ_CFG" | cut -c1-64)
  [ -n "$want" ] && [ "$want" = "$have" ] || { echo "STOP: cfg ≠ cfg z decyzji odsłonięcia ($have vs ${want:-brak})"; return 3; }
  return 0
}
uj_env() {   # = pf_env (LCG + venv + $PF/src) + kontrola funkcji używanych przez ubjk_run.py
  pf_env || return 1
  export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
  python3 - <<'PY' || { echo "BŁĄD: brak funkcji charmtag potrzebnych UBJK"; return 1; }
import os, charmtag, pyarrow  # noqa: F401
from charmtag.syst.files import file_table, assign_files, per_file_hists, scaled_inputs, unblinded_jackknife  # noqa: F401
from charmtag.pairs.fit2d import fit_channel  # noqa: F401
from charmtag.pairs.frozen import sbsb_counts  # noqa: F401
assert charmtag.__file__.startswith(os.environ["PF"] + "/src"), charmtag.__file__
PY
}
