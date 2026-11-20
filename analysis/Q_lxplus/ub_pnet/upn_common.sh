# upn_common.sh — paczka UBPNET (μ par D–D z kategoriami sieci par PAIRNET v3, 428770, odsłonięte; 19.11.2026; source, nie uruchamiać).
# Korzysta z ../ub_jk/ubjk_common.sh (UB_OUT, UJ_PAIRS, UJ_LEGS, uj_gate = blokada UNBLIND_OK + hash cfg, uj_env = pf_env: LCG + venv + $PF/src)
# oraz z paczki ../pairnet (pn_extract, pn_train, pn_fit, pn_hybsf, pn_sbeval, pn_dpcal, dn_feat, ext/, src/ — łańcuch v3; tylko odczyt).
UPN=${UPN:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$UPN/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/ub_jk/ubjk_common.sh"
: "${UPN_OUT:=$EOS_OUT/charm_q/ub_pnet}"                 # parts/pnf_data_428770_<g>.npz, evset.npz, ubpnet.{json,md}, ubpnet.log, DONE
: "${UPN_PN:=$Q/pairnet}"                                # kod PAIRNET (v3) na lxplus
: "${H5V20:=$EOS_OUT/charm_h5_v20/parts}"                # nogi par (jak łańcuch v3)
: "${H5V11:=$EOS_OUT/charm_h5_v11/parts}"                # próbka nominalna odsłonięcia (data_none, 39 plików) → zbiór zdarzeń
: "${PN_MODELS:=$RUNS/models}"                           # nominalne modele XGB H (D0_H, Dplus_H)
: "${P3:=$EOS_OUT/charm_q/pairnet_v3}"                   # części MC v3 (pn_*_nd_*.npz)
: "${SB3:=$EOS_OUT/charm_q/pairnet_sb3}"                 # train/pn_models_sbtrain.pt, xgbd_fold*.json; parts/pn_data_428770_*.npz (tylko pasma)
: "${HY3:=$EOS_OUT/charm_q/pairnet_hyb3}"                # res/hybsf.json, ds/hybds.json
: "${UPN_NG:=7}"                                          # plików h5 na węzeł ekstrakcji (jak dext v3: data_428770:7)
: "${UPN_K:=6}"; : "${UPN_NBOOT:=40}"; : "${UPN_CPUS:=8}"
export UPN Q UPN_OUT UPN_PN H5V20 H5V11 PN_MODELS P3 SB3 HY3 UPN_NG UPN_K UPN_NBOOT UPN_CPUS
upn_h5() { ls "$H5V20/data_428770/"*_none.h5 2>/dev/null | grep -v '\.part$' | sort; }
upn_v11() { ls "$H5V11/data_428770/"data_428770_*_none.h5 2>/dev/null | grep -v '\.part$' | sort; }
upn_ngroups() { local n; n=$(upn_h5 | wc -l); echo $(( (n + UPN_NG - 1) / UPN_NG )); }
upn_put() { uj_put "$1" "$2"; }
