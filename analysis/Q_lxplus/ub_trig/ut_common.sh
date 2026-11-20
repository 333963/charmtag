# ut_common.sh — paczka UBTRIG (μ per klasa triggera, 428770, odsłonięte; 18.11.2026; source, nie uruchamiać).
# Korzysta z ../ub_jk/ubjk_common.sh (UB_OUT, UJ_PAIRS, UJ_LEGS, UJ_CFG, uj_gate = blokada UNBLIND_OK + hash cfg, uj_env = pf_env).
UT=${UT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$UT/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/ub_jk/ubjk_common.sh"
: "${UT_OUT:=$EOS_OUT/charm_q/ub_trig}"                 # trg/trg_<g>.npz, ubtrig.{json,md}, DONE
: "${RAWV20:=$EOS_OUT/charm_raw_v20}"                    # ROOT v2.0 (dpt3to150) — tabela triggerów 428770
: "${BTD:=$EOS_OUT/charm_q/b_trig}"                      # bt_data_428770.npz, bt_fd_nd.npz (B_trig, tylko odczyt)
: "${UT_NG:=10}"                                          # grupy plików ROOT dla tabeli trg
export UT Q UT_OUT RAWV20 BTD UT_NG
ut_files() { ls "$RAWV20/data_428770"/*dpt3to150*.root 2>/dev/null; }
ut_put() { uj_put "$1" "$2"; }
