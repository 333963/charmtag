# mu_common.sh — paczka mult_charm (MULT_charm, 16.11.2026; source, nie uruchamiać). Wymaga ../common.sh (Q: EOS_OUT, q_env, put_output, fetch_input, eos_exists, job_tag).
# Na lxplus: ~/private/charm_q/Q_lxplus/mult_charm = zawartość charm_work/MULT_charm/lxplus + mult_common.py mult_extract.py mult_fit.py mult_analyze.py mult_gates.py
MU=${MU:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$MU/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"
: "${RAW20:=$EOS_OUT/charm_raw_v20}"                  # ROOT v2.0 (charm-prod-v2.0, Oct28): {data_428770,mb_nd,data_427929}/…/*.root
: "${MU_OUT:=$EOS_OUT/charm_q/mult_charm}"           # parts/<próbka>/<próbka>_<job>.npz, res/, DONE
: "${PTPROD:=dpt3to150}"                              # NIE mieszać d23 i d150 (osobne zadania grid, inne dzielenie plików)
: "${SAMPLES:=mb_nd data_427929 data_428770}"
: "${MAXJOBS:=40}"
export MU Q RAW20 MU_OUT PTPROD SAMPLES
mu_lists() {   # $MU/lists/<próbka>.txt (pliki *.root z nazwą zawierającą PTPROD)
  mkdir -p "$MU/lists"; local s
  for s in $SAMPLES; do
    if [ -d "$RAW20/$s" ]; then find "$RAW20/$s" -name "*${PTPROD}*.root" | sort > "$MU/lists/$s.txt"
    else : > "$MU/lists/$s.txt"; fi
  done
}
mu_env() {
  q_env
  export PYTHONDONTWRITEBYTECODE=1
  python3 -c "import uproot, awkward, scipy, numpy, pandas" || { echo "BŁĄD: brak uproot/awkward/scipy w LCG+venv"; return 1; }
}
