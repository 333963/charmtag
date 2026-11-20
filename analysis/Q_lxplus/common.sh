# common.sh — wspólne ścieżki i funkcje (source, nie uruchamiać). Etap Q, 2026-11-02.
# Każdą zmienną można nadpisać w środowisku przed `source`.
: "${Q:=/afs/cern.ch/user/k/klyczak/private/charm_q/Q_lxplus}"     # ten katalog (AFS)
: "${QBASE:=/afs/cern.ch/user/k/klyczak/private/charm_q}"
: "${VENV:=$QBASE/venv}"
: "${LCG_SETUP:=/cvmfs/sft.cern.ch/lcg/views/LCG_110_swan/x86_64-el9-gcc13-opt/setup.sh}"
: "${RAW:=/eos/user/m/malewick/tmp/charm_ntuples_v1.1}"           # surowe ntuple ROOT v1.1 (tylko odczyt)
: "${EOS_OUT:=/eos/user/k/klyczak}"
: "${H5_PARTS:=$EOS_OUT/charm_h5_v11/parts}"                      # h5: parts/<próbka>/<próbka>_<job>_<presel>.h5
: "${RUNS:=$EOS_OUT/charm_q/runs/lxplus_v11}"                     # charmtag out_dir
: "${TRUTH:=$EOS_OUT/charm_q/truth}"                              # Pythia gen + mapa wydajności
: "${XRD:=root://eosuser.cern.ch}"
export Q QBASE VENV LCG_SETUP RAW EOS_OUT H5_PARTS RUNS TRUTH XRD

q_env() {   # LCG_110_swan + venv (z --system-site-packages)
  set +u
  # shellcheck disable=SC1090
  source "$LCG_SETUP"
  if [ -f "$VENV/bin/activate" ]; then source "$VENV/bin/activate"; export PYTHONNOUSERSITE=1; fi   # bez ~/.local w venv
  export PYTHONPATH="$(python3 -c "import site;print(site.getsitepackages()[0])"):${PYTHONPATH:-}"  # VENV_FIRST (xgboost venv > LCG)
  set -u
  if [ -z "${PYTHIA8DATA:-}" ] && command -v pythia8-config >/dev/null 2>&1; then
    PYTHIA8DATA=$(pythia8-config --xmldoc 2>/dev/null || true); export PYTHIA8DATA
  fi
}

# Kopiuj wejście z EOS do katalogu roboczego (xrdcp); gdy się nie da → zwróć ścieżkę FUSE.
fetch_input() {  # fetch_input /eos/... katalog_docelowy  → echo ścieżka lokalna
  local src=$1 dst=$2/$(basename "$1")
  if command -v xrdcp >/dev/null 2>&1 && xrdcp -s -f "$XRD/$src" "$dst" 2>/dev/null; then echo "$dst"; else echo "$src"; fi
}

# Zapisz plik na EOS (xrdcp, fallback: cp przez FUSE). put_output lokalny /eos/.../nazwa
put_output() {
  local src=$1 dst=$2
  if command -v xrdcp >/dev/null 2>&1; then
    xrdfs "${XRD#root://}" mkdir -p "$(dirname "$dst")" 2>/dev/null || true
    xrdcp -s -f "$src" "$XRD/$dst" && return 0
  fi
  mkdir -p "$(dirname "$dst")" && cp -f "$src" "$dst"
}

# Czy plik istnieje na EOS i ma > 0 B
eos_exists() { [ -s "$1" ] || { command -v xrdfs >/dev/null 2>&1 && xrdfs "${XRD#root://}" stat "$1" >/dev/null 2>&1; }; }

job_tag() {   # user.malewick.51907090.EXT0._000011.CharmAnalysis....root → 000011
  basename "$1" | grep -o '\._[0-9]\{6\}\.' | head -1 | tr -d '._'
}
