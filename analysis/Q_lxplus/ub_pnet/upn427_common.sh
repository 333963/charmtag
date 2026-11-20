# upn427_common.sh — UBPNET427 (19.11.2026, projekt A; zatwierdzone 18.11): sieć par PAIRNET v3 na czystym runie 427929 (klasa MB).
# Nakładka na upn_common.sh (bieg 428770 bez zmian): inne wejście (h5 v2.0 427929 po filtrze UB427: L_xy ≥ 0,25, okno D0, klasa MB),
# inny katalog wyjścia, kontrola nominału na parach UB427 (pairs_full.npz + fit2d_unblind.json + pair_inputs z N_D* = 2046 i c_w UB427).
# source, nie uruchamiać.
UPN=${UPN:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
# shellcheck disable=SC1091
source "$UPN/upn_common.sh"
UPN_OUT=${UPN427_OUT:-$EOS_OUT/charm_q/ub_pnet427}         # parts/pnf_data_427929_<g>.npz, ubpnet.{json,md}, UBPNET427.md, DONE
U4_DIR=${U4_DIR:-$EOS_OUT/charm_q/ub_427929}               # wyjście UB427: h5/ (458 × *_mb.h5, filtr), pairs/ (nominał 427929)
U427_PAIRS=$U4_DIR/pairs
UPN_NG=${UPN427_NG:-33}                                    # 458 plików (po filtrze ~12 MB) → 14 grup
UPN_FCH=${UPN_FCH:-1.032,1.266,1.060}                     # normalizacja z prawdy mb_nd f_ch (GEN_mcstat/PROCEDURA.md §6; 427929: f jak 428770)
export UPN_OUT U4_DIR U427_PAIRS UPN_NG UPN_FCH
upn_h5() { ls "$U4_DIR/h5/"data_427929_*_mb.h5 2>/dev/null | grep -v '\.part$' | sort; }
