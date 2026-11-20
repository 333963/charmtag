# nw_common.sh — paczka UB_NDSW (19.11.2026, projekt A): skala σ piku D* w klasach triggera 428770 vs 427929. source, nie uruchamiać.
NW=${NW:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$NW/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/pf_pairs/pf_common.sh"                         # EOS_OUT, q_env, put_output
: "${NW_H5_8:=$EOS_OUT/charm_h5_v11/parts/data_428770}"   # h5 v1.1 428770 (*_none.h5; próbka nominału N_D* 05.11) — tylko odczyt
: "${NW_H5_7:=$EOS_OUT/charm_q/ub_427929/h5}"             # h5 427929 po filtrze UB427 — tylko odczyt
: "${NW_TRG:=$EOS_OUT/charm_q/ub_trig/trg}"               # trg_*.npz (ub_trig) — tylko odczyt
: "${NW_OUT:=$EOS_OUT/charm_q/ub_ndsw}"                   # wynik: NDSW.md, ndsw.json, ndsw_hists.npz, nw.log, DONE
export NW Q NW_H5_8 NW_H5_7 NW_TRG NW_OUT
nw_env() { q_env >/dev/null 2>&1 || return 1; export PYTHONPATH="$NW:$Q/ub_trig" OMP_NUM_THREADS=2
  python3 -c "import h5py, iminuit, scipy, numpy, dstar_refit_nds, nds_var, ubtrig_run, btrig_common" || { echo "BŁĄD: importy UB_NDSW"; return 1; }; }
