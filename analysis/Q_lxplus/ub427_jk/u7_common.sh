# u7_common.sh — paczka UB427JK (jackknife μ po GRUPACH ZDARZEŃ dla runu 427929; 19.11.2026, projekt A; source, nie uruchamiać).
# Cel: czynnik σ_jk/σ_profil68 dla 427929 (w 428770 jk_ev200 dał ×1,18), żeby Δμ(427929 − 428770) z UB427 miało błąd z korelacjami zdarzeń.
# Kod: ../ub_jk/ubjk_run.py BEZ ZMIAN (prep --groups events:G z przebudową nóg; fit; finish), wejścia modelu z biegu UB427 ($U4_OUT/pairs):
#   pair_inputs.{json,npz}, pairs_full.npz, fit2d_unblind.json. legs_data.parquet 427929 nie był zapisany → prep przebudowuje pary
#   (cmd_pairs --unblind) z TYM SAMYM cfg co węzeł `pairs` UB427 (merge_cfg z cfg odsłonięcia + h5 *_mb + scores + mb_rates_427929.json);
#   kontrola ubjk: Σ grup == pairs_full.npz bitowo, fit nominalny == fit2d_unblind.json 427929 (|Δμ| < 2·10⁻³) → inaczej STOP.
# Pary w danych zatwierdzone (projekt A): 18.11.2026. Blokada: $U4_OUT/DONE (UB427 zakończony) i komplet wejść.
U7=${U7:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$U7/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/ub_427929/u4_common.sh"                        # U4_OUT, U4_NOM_CFG, pf_env, put_output, NOMINAL, NORM (przez pf_common)
: "${U7_G:=200}"                                         # grupy zdarzeń (jak jk_ev200 w 428770; precyzja σ_jk ±5 %)
: "${U7_NCHUNK:=25}"                                     # porcje replik (201 replik ≈ 8/porcję; tylko tryb fixed)
: "${U7_MODES:=fixed}"
: "${U7_OUT:=$U4_OUT/jk_ev$U7_G}"                        # EOS: perfile.*, fit_<i>.json, ubjk.{json,md}, legs_data_rebuilt.parquet, DONE
U7_PAIRS=$U4_OUT/pairs
UJ_RUN=$Q/ub_jk/ubjk_run.py
export U7 Q U7_G U7_NCHUNK U7_MODES U7_OUT U7_PAIRS UJ_RUN
u7_gate() {
  local f ok=0
  [ -s "$U4_OUT/DONE" ] || { echo "STOP: brak $U4_OUT/DONE (UB427 niezakończony)"; ok=3; }
  for f in pair_inputs.json pair_inputs.npz pairs_full.npz fit2d_unblind.json; do [ -s "$U7_PAIRS/$f" ] || { echo "STOP: brak $U7_PAIRS/$f"; ok=3; }; done
  [ -s "$U4_OUT/cw/mb_rates_427929.json" ] || { echo "STOP: brak $U4_OUT/cw/mb_rates_427929.json"; ok=3; }
  for f in D0_H Dplus_H; do [ -s "$U4_OUT/scores/$f/data_none.h5" ] || { echo "STOP: brak $U4_OUT/scores/$f/data_none.h5"; ok=3; }; done
  [ -s "$UJ_RUN" ] || { echo "STOP: brak $UJ_RUN (paczka ub_jk na lxplus)"; ok=3; }
  return $ok
}
u7_env() {   # = pf_env (LCG + venv + pf_pairs/src) + ubjk_run importowalny
  pf_env || return 1
  export PYTHONPATH="$Q/ub_jk:${PYTHONPATH:-}" OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  python3 -c "import charmtag, pyarrow, ubjk_run" || { echo "BŁĄD: importy UB427JK"; return 1; }
}
