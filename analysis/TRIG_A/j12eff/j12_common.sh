# j12_common.sh — paczka J12EFF (17.11, projekt A): wydajność L1_J12 z danych vs dżet torowy; tylko zdarzenia, bez par D–D. Wymaga ../common.sh.
J12=${J12:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$J12/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"
: "${RAWV20:=$EOS_OUT/charm_raw_v20}"            # wejście (odczyt): ROOT v2.0 z drzewem Tracks
: "${J12_OUT:=$EOS_OUT/charm_q/j12eff}"
: "${J12_SAMPLES:=data_428770 data_427929}"
: "${J12_NPARTS:=data_428770:12 data_427929:36}"
export J12 Q RAWV20 J12_OUT J12_SAMPLES J12_NPARTS
j12_np() { local kv; for kv in $J12_NPARTS; do [ "${kv%%:*}" = "$1" ] && { echo "${kv##*:}"; return; }; done; echo 4; }
j12_env() { q_env; export PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=${TMPDIR:-/tmp}
  python3 -c 'import numpy, uproot, awkward, numba, scipy' || { echo "BŁĄD: brak numpy/uproot/awkward/numba/scipy w q_env"; return 1; }; }
