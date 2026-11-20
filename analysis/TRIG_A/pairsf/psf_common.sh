# psf_common.sh — paczka TRIG_A/pairsf (projekt A: SF klasyfikatora w kategoriach score par per klasa triggera; pojedyncze D0←D*). source.
# Na lxplus: ~/private/charm_q/Q_lxplus/pairsf (= lokalnie charm_work/TRIG_A/pairsf + btrig_common.py z B_trig + src/charmtag z AA_merge).
PSF=${PSF:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$PSF/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"
: "${PSF_IN:=$EOS_OUT/charm_q/b_trig}"        # wyniki etapu b_trig (tylko odczyt): bt_<próbka>.npz, parts/bt_<próbka>_<g>.npz
: "${PSF_OUT:=$EOS_OUT/charm_q/pairsf}"       # psf.json, PSF.md, DONE
: "${PSF_NPROC:=8}"
export PSF Q PSF_IN PSF_OUT PSF_NPROC
psf_glob() {   # scalony plik, a gdy go nie ma — części
  if [ -s "$PSF_IN/bt_$1.npz" ]; then echo "$PSF_IN/bt_$1.npz"; else echo "$PSF_IN/parts/bt_$1_*.npz"; fi
}
