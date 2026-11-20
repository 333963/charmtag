#!/bin/bash
# Jedyny węzeł pairsf: fity pass/fail (≈150) na PSF_NPROC rdzeniach → $PSF_OUT/{psf.json,PSF.md,DONE}. Scratch condora; AFS tylko logi.
set -eo pipefail
PSF=$(cd "$(dirname "$0")" && pwd); export PSF
source "$PSF/psf_common.sh"
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/psf_XXXXXX")
trap 'cd /; rm -rf "$WD"' EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
echo "[$(date -u +%FT%TZ)] PSF host $(hostname) WD=$WD IN=$PSF_IN OUT=$PSF_OUT nproc=$PSF_NPROC"
q_env
export PYTHONPATH="$PSF:$PSF/src:${PYTHONPATH:-}" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -c "import charmtag.dstar_refit as d, btrig_common; print('dstar_refit:', d.__file__); assert hasattr(d, 'fit_passfail')"
python3 -u "$PSF/psf_analyze.py" --data "$(psf_glob data_428770)@428770" --data "$(psf_glob data_427929)@427929" \
    --mc "$(psf_glob fd_nd)" --nproc "$PSF_NPROC" --out "$WD/res"
for f in "$WD"/res/*; do put_output "$f" "$PSF_OUT/$(basename "$f")"; done
date -u +%FT%TZ > DONE; put_output DONE "$PSF_OUT/DONE"
echo "[$(date -u +%FT%TZ)] OK PSF"
