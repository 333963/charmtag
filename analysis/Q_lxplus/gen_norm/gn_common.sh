# gn_common.sh — paczka gen_norm (16.11.2026; v2.0: 17.11; source, nie uruchamiać). Wymaga ../common.sh (Q: EOS_OUT, RUNS, q_env, put_output, fetch_input).
# GN_FMT=v20 (domyślnie): ROOT v2.0 (CharmAnalysis Oct28) z $RAW_GN/{mb_nd,data_428770}/*${PTPROD}*.root — common.sh (RAW v1.1) NIE jest zmieniany,
#   nadpisanie tylko tutaj. GN_FMT=v11: dawne zachowanie (RAW v1.1 + config/samples.txt; ścieżka v1.1 na EOS 17.11 niedostępna).
GN=${GN:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$GN/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"
: "${GN_FMT:=v20}"
: "${PTPROD:=dpt3to150}"                             # v2.0: NIE mieszać dpt3to150 i dpt2to3 (te same zdarzenia, inne dzielenie plików)
: "${RAW_GN:=$EOS_OUT/charm_raw_v20}"                # v2.0: {mb_nd,data_428770}/…/*.root
if [ "$GN_FMT" = v20 ]; then
  : "${GN_OUT:=$EOS_OUT/charm_q/gen_norm_v20}"       # osobno od ewentualnego biegu v1.1
  : "${GN_MB_EXCLUDE:=}"                             # §3b (18.11): tagi mb_nd bez ROOT (np. 000037,000070) — wtedy N_MB_EXP = 80 − liczba tagów
  if [ -n "$GN_MB_EXCLUDE" ]; then : "${N_MB_EXP:=$((80 - $(echo "$GN_MB_EXCLUDE" | tr ',' '\n' | grep -c .)))}"; fi
  : "${N_MB_EXP:=80}"; : "${N_DATA_EXP:=97}"
  : "${DATA_REF_TAGS:=}"                             # v2.0: inne dzielenie plików → referencja po zdarzeniach (REF_ANCHORS), nie po tagach
  : "${REF_ANCHORS:=$EOS_OUT/charm_h5_v11/parts/data_428770/data_428770_000021_none.h5,$EOS_OUT/charm_h5_v11/parts/data_428770/data_428770_000030_none.h5}"
  MB_DIR=mb_nd; DATA_DIR=data_428770; GN_RAW=$RAW_GN
else
  : "${GN_OUT:=$EOS_OUT/charm_q/gen_norm}"
  : "${N_MB_EXP:=80}"; : "${N_DATA_EXP:=39}"
  : "${DATA_REF_TAGS:=000021,000030}"                # pliki danych, z których etap1/t4 wyznaczył w2 (441 691 zdarzeń)
  : "${REF_ANCHORS:=}"
  MB_DIR=$(awk '$1=="mb_nd"{print $2}' "$Q/config/samples.txt")
  DATA_DIR=$(awk '$1=="data_428770"{print $2}' "$Q/config/samples.txt"); GN_RAW=$RAW
fi
: "${BASE_CFG:=$Q/config/lxplus_v11.yaml}"           # config nominalny (h5 v1.1 parts, out_dir = nominalny bieg) — nogi ZAWSZE z h5 v1.1
: "${NOMINAL:=$RUNS}"                                # nominalne scores/<sp>_H/incl_mc.h5 (tylko odczyt)
: "${DATA_SCOPE:=all}"                               # v11: all = 39 plików, ref = tylko DATA_REF_TAGS; v20: zawsze all (ref wybierany w combine)
: "${REF_RATES:=1344,9534,8665,22,41,15}"            # n_D*, n_D0, n_D+, pary (PF2 pairs/pair_inputs.json meta.mc_rates_unweighted)
: "${MAXJOBS:=40}"
export GN_MB_EXCLUDE GN Q GN_FMT PTPROD RAW_GN GN_RAW GN_OUT BASE_CFG NOMINAL DATA_SCOPE DATA_REF_TAGS REF_ANCHORS REF_RATES MB_DIR DATA_DIR N_MB_EXP N_DATA_EXP
gn_lists() {   # listy ROOT → $GN/lists/{mb_nd,data_428770}.txt (własne; NIE nadpisuje Q/lists/*)
  mkdir -p "$GN/lists"
  local k d
  for k in mb_nd data_428770; do
    d=$MB_DIR; [ "$k" = data_428770 ] && d=$DATA_DIR
    if [ "$GN_FMT" = v20 ]; then   # find (podkatalogi) + filtr PTPROD
      if [ -d "$GN_RAW/$d" ]; then find "$GN_RAW/$d" -name "*${PTPROD}*.root" | sort > "$GN/lists/$k.txt"
      else echo "BŁĄD: nie widzę $GN_RAW/$d"; return 1; fi
    elif ls -d "$GN_RAW/$d" >/dev/null 2>&1; then ls -1 "$GN_RAW/$d" | grep '\.root$' | sort | sed "s|^|$GN_RAW/$d/|" > "$GN/lists/$k.txt"
    elif command -v xrdfs >/dev/null 2>&1; then xrdfs "${XRD#root://}" ls "$GN_RAW/$d" | grep '\.root$' | sort > "$GN/lists/$k.txt"
    else echo "BŁĄD: nie widzę $GN_RAW/$d"; return 1; fi
  done
}
gn_env() {     # LCG + venv (q_env) i NASZ src (snapshot charmtag 2.1.0) przed zainstalowanym charmtag
  q_env
  export PYTHONPATH="$GN/src:${PYTHONPATH:-}" PYTHONDONTWRITEBYTECODE=1
  python3 - <<'PY' || { echo "BŁĄD: charmtag nie z $GN/src albo brak uproot/awkward/pyarrow"; return 1; }
import os, charmtag, uproot, awkward, pyarrow  # noqa: F401
from charmtag.pairs import build
assert charmtag.__file__.startswith(os.environ["GN"] + "/src"), charmtag.__file__
assert hasattr(build, "check_event_key_unique"), "src bez check_event_key_unique (inny snapshot)"
PY
}
