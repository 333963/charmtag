# topo_common.sh — paczka UB_TOPO (19.11.2026, projekt A; source). Mechanika = pd_diff (charmtag pairs/fit2d --unblind --diff), osobny cfg/out/batch.
UT=${UT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
# shellcheck disable=SC1091
source "$UT/../pd_diff/pd_common.sh"            # PF, Q, NORM, UB_OUT, pf_env/pd_env (src = pf_pairs/src z diff.py); PD_* nadpisane niżej
: "${UT_VAR:=dy}"                                # dy (jedyna zmienna z mocą, PROGNOZA.md); dphi tylko na życzenie (moc ≈ 0)
: "${UT_CHANNELS:=D0D0,DpD0,DpDp}"
: "${UT_OUT:=$EOS_OUT/charm_q/ub_topo/$UT_VAR}"
: "${UT_MU_REF:=1.50}"
VAR=$UT_VAR; CHANNELS=$UT_CHANNELS; PD_OUT=$UT_OUT; CHARGE_SRC=""
export UT UT_VAR UT_CHANNELS UT_OUT UT_MU_REF VAR CHANNELS PD_OUT CHARGE_SRC
