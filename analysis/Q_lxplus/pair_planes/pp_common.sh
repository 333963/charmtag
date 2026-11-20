# pp_common.sh — paczka PAIR_PLANES (19.11, projekt A: płaszczyzny mas par D–D jak plots_doubleD z kodu bazowego charmpp, z wetem wspólnych torów,
# OS/LS, wariantem weta D* i |Δφ|/Δy; source, nie uruchamiać). Wymaga ../common.sh.
PX=${PX:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$PX/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"
: "${H5V20:=/eos/user/k/klyczak/charm_h5_v20/parts}"
: "${PP_OUT:=$EOS_OUT/charm_q/pair_planes}"                       # <próbka>.npz + .meta.json, PP_SUMMARY.json, DONE; części w parts/
: "${PP_SAMPLES:=data_427929:24 data_428770:8 mb_nd:6 fd_nd:3}"   # próbka:liczba węzłów
: "${PP_NFILES:=data_427929:458 data_428770:97 mb_nd:78 fd_nd:24}"
: "${PP_DISK:=6GB}"
export PX Q H5V20 PP_OUT PP_SAMPLES PP_NFILES PP_DISK
pp_kv() { local kv; for kv in $2; do [ "${kv%%:*}" = "$1" ] && { echo "${kv##*:}"; return; }; done; echo 0; }
pp_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else cp -f "$1" "$2" && rm -f "$1"; fi; }
pp_env() { q_env; export PYTHONDONTWRITEBYTECODE=1; python3 -c 'import numpy, h5py' || { echo "BŁĄD: brak numpy/h5py"; return 1; }; }
