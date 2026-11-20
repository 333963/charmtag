# uc_common.sh — paczka UBCOMB (wspólne μ z próbek bez triggera dżetowego; 19.11.2026; source, nie uruchamiać).
# Korzysta z ../ub_jk/ubjk_common.sh (UJ_PAIRS = pary 428770 05.11, UJ_LEGS, uj_gate = blokada UNBLIND_OK + hash cfg, uj_env = pf_env).
UC=${UC:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$UC/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/ub_jk/ubjk_common.sh"
: "${UC_OUT:=$EOS_OUT/charm_q/ub_comb}"                  # ubcomb.{json,md}, ubcomb.log, DONE
: "${UC_TRG:=$EOS_OUT/charm_q/ub_trig/trg}"               # trg_<g>.npz z UBTRIG (bity trigPassed 428770)
: "${BTD:=$EOS_OUT/charm_q/b_trig}"                       # bt_data_428770.npz, bt_fd_nd.npz
: "${U4_OUT:=$EOS_OUT/charm_q/ub_427929}"                 # pairs/ (427929 MB), ndstar_427929.json, jk_ev200/ubjk.json (opcjonalnie)
: "${UC_VARIANTS:=strict,nonj12,nonj12_noafp}"
export UC Q UC_OUT UC_TRG BTD U4_OUT UC_VARIANTS
uc_put() { uj_put "$1" "$2"; }
