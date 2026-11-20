# PF2: powtórka par z poprawką weta wspólnych śladów (pairs → fit2d ślepo → toys → syst par → nom_mix) na nominalnych score z EOS.
# Paczkę i submit_pf.sh check|submit|status tworzy etap PF2 (charm_work/PF2_lxplus → Q_lxplus/pf_pairs). Brak paczki → etap pominięty.
source "$HERE/stages/_common.sh"
SRC=$QL/pf_pairs; R=$RQ/pf_pairs; EOSD=$EOSQ/pf_pairs; DST=$W/PF2_lxplus/runs/lxplus
stage_available() { [ -f "$SRC/submit_pf.sh" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/pf_pairs" && rpush --exclude cfg/ "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && ./submit_pf.sh check" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_pf.sh submit; ./submit_pf.sh status | tail -5" | tail -10 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull "${mdjson[@]}" --include '*.npz' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "wyniki → charm_work/PF2_lxplus/runs/lxplus (compare_pf.py)"; }
