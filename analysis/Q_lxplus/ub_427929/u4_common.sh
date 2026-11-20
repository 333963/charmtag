# u4_common.sh — paczka UB427 (μ par D–D w czystym runie 427929, klasa MB; 18.11.2026; source, nie uruchamiać).
# Ten sam charmtag (pf_pairs/src), modele, progi, szablony i model fit2d co odsłonięcie 428770 (05.11); zmieniają się dane (h5 v2.0 427929
# po cięciu L_xy/okna ≡ v1.1), N_D* (fit w 427929) i c_w (wagi w2 do rozkładu 427929). Pary w danych zatwierdzone: 18.11.2026.
U4=${U4:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$U4/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/pf_pairs/pf_common.sh"                         # PF_OUT, NOMINAL, NORM, pf_env, q_env, put_output, fetch_input
: "${U4_OUT:=$EOS_OUT/charm_q/ub_427929}"                 # h5/ (filtr), fstats/, gn/{mc,d428,d427}/, cw/, scores/, wynik UB427.{md,json}, DONE
: "${H5V20:=$EOS_OUT/charm_h5_v20/parts}"                 # data_427929/data_427929_NNNNNN_none.h5 (etap h5_v20)
: "${RAWV20:=$EOS_OUT/charm_raw_v20}"                     # ROOT v2.0 (dpt3to150): mb_nd (78/80), data_428770 (97), data_427929 (458)
: "${U4_NOM_PAIRS:=$PF_OUT/unblind/pairs}"                # odsłonięcie 428770 (porównanie, ten sam fit)
: "${U4_NOM_CFG:=$Q/pf_pairs/cfg/pf_unblind.yaml}"        # cfg odsłonięcia → baza cfg 427929 (merge_cfg z nadpisaniami)
: "${U4_NOM_MB:=$Q/inputs/mb_rates_2d.json}"              # zamrożone c_w nominału
: "${U4_ANCHORS:=$EOS_OUT/charm_h5_v11/parts/data_428770/data_428770_000021_none.h5,$EOS_OUT/charm_h5_v11/parts/data_428770/data_428770_000030_none.h5}"
: "${U4_NF:=12}"                                          # grupy plików h5 427929 (filtr)
: "${U4_NMC:=13}"                                         # grupy plików ROOT mb_nd (gn_root mc; 78 plików → 6/grupa)
: "${U4_N8:=10}"                                          # grupy ROOT 428770 (gn_root data, 97 plików)
: "${U4_N7:=12}"; : "${U4_STRIDE7:=4}"                    # ROOT 427929: co U4_STRIDE7-ty plik (≈115 z 458; mapa w2 — wystarczy podpróbka)
export U4 Q U4_OUT H5V20 RAWV20 U4_NOM_PAIRS U4_NOM_CFG U4_NOM_MB U4_ANCHORS U4_NF U4_NMC U4_N8 U4_N7 U4_STRIDE7
u4_h5() { ls "$H5V20/data_427929"/data_427929_*_none.h5 2>/dev/null; }
u4_root() {   # u4_root mb_nd|data_428770|data_427929
  case $1 in
    data_427929) ls "$RAWV20/$1"/*dpt3to150*.root 2>/dev/null | awk -v s="$U4_STRIDE7" '(NR - 1) % s == 0' ;;
    *) ls "$RAWV20/$1"/*dpt3to150*.root 2>/dev/null ;;
  esac
}
u4_group() { awk -v g="$1" -v n="$2" '(NR - 1) % n == g'; }   # stdin: lista → pliki grupy g z n
u4_env() {    # = pf_env (LCG + venv + pf_pairs/src) + uproot/awkward/h5py dla gn_root/filtra
  pf_env || return 1
  export PYTHONPATH="$U4:$Q/ub_jk:${PYTHONPATH:-}" OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
  python3 -c "import uproot, awkward, h5py, yaml, iminuit, charmtag, ubjk_run" || { echo "BŁĄD: importy UB427"; return 1; }
}
