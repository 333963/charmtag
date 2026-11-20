# pf_common.sh — paczka PF2 (powtórka par z poprawką weta wspólnych śladów; source, nie uruchamiać). Wymaga ../common.sh (Q).
PF=${PF:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$PF/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"
: "${PF_OUT:=$EOS_OUT/charm_q/pf_pairs}"        # wyjścia (EOS): pairs/, pairs_syst/, mix/, …
: "${NOMINAL:=$RUNS}"                            # nominalny bieg charmtag (modele/score/validation) — tylko odczyt
: "${BASE_CFG:=$Q/config/lxplus_v11.yaml}"
: "${NORM:=direct}"                              # direct (porównywalne z ½ lumi/wrześniem) | truth (gdy pairs.truth_norm_json w cfg)
export PF PF_OUT NOMINAL BASE_CFG NORM
pf_env() {   # LCG + venv (q_env) i NASZ src (AA_merge + PF) przed zainstalowanym charmtag (editable charmtag_src nietknięty)
  q_env
  export PYTHONPATH="$PF/src:${PYTHONPATH:-}" PYTHONDONTWRITEBYTECODE=1
  python3 - <<'PY' || { echo "BŁĄD: charmtag nie z $PF/src albo bez poprawki PF"; return 1; }
import os, charmtag
from charmtag.pairs import build
assert charmtag.__file__.startswith(os.environ["PF"] + "/src"), charmtag.__file__
assert build.pcfg({})["track_id"] == "d0z0", "brak poprawki PF (pairs.track_id)"
PY
}
