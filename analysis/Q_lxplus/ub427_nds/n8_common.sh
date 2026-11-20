# n8_common.sh — paczka UB427_NDS (19.11.2026, projekt A, W2): systematyka N_D* w 427929 (skan wariantów NDS). source, nie uruchamiać.
N8=${N8:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$N8/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/pf_pairs/pf_common.sh"                         # EOS_OUT, q_env, put_output
: "${U4_OUT:=$EOS_OUT/charm_q/ub_427929}"                 # wejście (tylko odczyt): h5/data_427929_*_mb.h5 z filtra UB427
: "${N8_OUT:=$EOS_OUT/charm_q/ub427_nds}"                 # wynik: NDS_427929.md, ndstar_variants_427929.json, dstar_hist_427929.npz, DONE
export N8 Q U4_OUT N8_OUT
n8_env() { q_env >/dev/null 2>&1 || return 1; export PYTHONPATH="$N8" OMP_NUM_THREADS=2
  python3 -c "import h5py, iminuit, scipy, numpy, dstar_refit_nds, nds_var" || { echo "BŁĄD: importy UB427_NDS"; return 1; }; }
