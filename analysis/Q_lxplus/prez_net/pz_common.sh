# pz_common.sh — paczka PREZ-2 (19.11.2026, projekt A): klasyfikator pojedynczych D, BDT vs sieć × high/low-level. TYLKO MC. source, nie uruchamiać.
# Na lxplus: ~/private/charm_q/Q_lxplus/prez_net. src/ = migawka charmtag z AA_merge (src/SNAPSHOT.txt). Wymaga ../common.sh.
PZ=${PZ:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$PZ/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                          # EOS_OUT, RUNS, q_env, put_output
: "${PZ_OUT:=$EOS_OUT/charm_q/prez_net}"                       # wyjście: ext/, train_<gatunek>/, report/, DONE
: "${PZ_MODELS:=$RUNS/models}"                                 # nominalne XGB H (<gatunek>_H: fold*.json, meta.json, config_used.yaml) — tylko odczyt
: "${PZ_CFG:=$PZ_MODELS/D0_H/config_used.yaml}"                # config nominału (ta sama produkcja h5 v1.1, kolumny, prawda, weto D*)
: "${PZ_CAP_SIG:=1500000}"; : "${PZ_CAP_BKG_D0:=3000000}"; : "${PZ_CAP_BKG_Dplus:=5000000}"
: "${PZ_EPOCHS:=30}"; : "${PZ_CPUS:=8}"
export PZ Q PZ_OUT PZ_MODELS PZ_CFG PZ_CAP_SIG PZ_CAP_BKG_D0 PZ_CAP_BKG_Dplus PZ_EPOCHS PZ_CPUS
pz_env() {
  q_env
  export PYTHONPATH="$PZ:$PZ/src:${PYTHONPATH:-}" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=$PZ_CPUS MKL_NUM_THREADS=$PZ_CPUS
  python3 - <<'PY' || { echo "BŁĄD: środowisko pz (torch/xgboost/matplotlib/charmtag z \$PZ/src)"; return 1; }
import os, numpy, scipy, h5py, yaml, xgboost, torch, matplotlib, charmtag
from charmtag import samples, config, features
assert os.path.realpath(os.path.dirname(charmtag.__file__)).startswith(os.path.realpath(os.environ["PZ"] + "/src")), charmtag.__file__
print("pz_env: torch", torch.__version__, "xgboost", xgboost.__version__, "charmtag", os.path.dirname(charmtag.__file__))
PY
}
pz_put() { put_output "$1" "$2"; }
