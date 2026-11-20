# PD: pary różniczkowe (like-sign/opposite-sign, Δφ) na odsłoniętej próbce — fit-only na nominalnych score, src = pf_pairs/src
#   (= AA_merge/charmtag: NDS z 5.11 + PD diff.py). Paczka Q_lxplus/pd_diff (etap PD, 05.11). Wymaga unblind DONE (sprawdza submit_diff.sh).
#   ENV: PD_VAR (sign|dphi|dy|sign_dphi; dom. sign), PD_CHARGE_SRC (h5 v1.2 → ładunek K córki 0 dla D0; dom. charm_h5_v12/parts), PD_OUT_SUFFIX.
source "$HERE/stages/_common.sh"
: "${PD_VAR:=sign}" "${PD_CHARGE_SRC:=/eos/user/k/klyczak/charm_h5_v12/parts}" "${PD_OUT_SUFFIX:=_v12}"   # 06.11: bieg 1 (diff/) bez ładunków D0 (bug charge_files podkatalog) → powtórka do diff_v12/
SRC=$QL/pd_diff; R=$RQ/pd_diff; PFR=$RQ/pf_pairs; EOSD=$EOSQ/pf_pairs/diff${PD_OUT_SUFFIX:-}; DST=$W/PD_pairs_diff/runs/lxplus${PD_OUT_SUFFIX:-}
PDENV="VAR=$PD_VAR CHARGE_SRC=$PD_CHARGE_SRC PD_OUT=$EOSD"
stage_available() { [ -f "$SRC/submit_diff.sh" ] && [ -f "$QL/pf_pairs/src/charmtag/pairs/diff.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/pd" && rpush --exclude 'diff.dag*' --exclude cfg/ "$SRC/" "$H:$R/" && rpush --delete --exclude __pycache__ --exclude '.pytest_cache' "$QL/pf_pairs/src/charmtag/" "$H:$PFR/src/charmtag/"; }   # src: tylko charmtag/ (config/ lxplus nietknięty)
stage_check()  { rssh "cd ~/$R && chmod +x *.sh tools/*.py && $PDENV ./submit_diff.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && $PDENV ./submit_diff.sh submit; $PDENV ./submit_diff.sh status | tail -5" | tail -12 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { rssh "condor_q -batch 2>/dev/null | grep -E 'charm_pd'; condor_q -totals 2>/dev/null | tail -1" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null; rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include 'pairs_full_diff.npz' --include 'pair_inputs_diff.json' "${mdjson[@]}" "$H:$EOSD/" "$DST/"; }
stage_result() { echo "PD ($PD_VAR): $DST/tfin_diff.md (EOS $EOSD/)"; }
