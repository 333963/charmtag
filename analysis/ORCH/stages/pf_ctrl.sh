# PF3: kontrola izolująca efekt weta wspólnych śladów — dwa biegi fit-only (pairs.track_id pt_eta vs d0z0) na tym samym src 2.0.2 + compare_ctrl.
# Paczka Q_lxplus/pf_ctrl (etap PF3). Wymaga: pf_pairs done (nominalne score; wejście do bramki G6). Brak paczki → pominięte.
source "$HERE/stages/_common.sh"
SRC=$QL/pf_ctrl; R=$RQ/pf_ctrl; EOSD=$EOSQ/pf_ctrl; DST=$W/PF3_ctrl/runs/lxplus
stage_available() { [ -f "$SRC/submit_ctrl.sh" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/pf_ctrl && touch ~/$R/.orch_start" && rpush --exclude cfg/ --exclude 'ctrl.dag*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh tools/*.py && FORCE=${FORCE:-1} ./submit_ctrl.sh check" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && FORCE=${FORCE:-1} ./submit_ctrl.sh submit; ./submit_ctrl.sh status | tail -5" | tail -10 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
held_fix() {  # diag HELD (05.11): loguje powód; jednorazowy condor_release (licznik w state), potem już tylko log
  local n=0; [ -f "$STATE/pf_ctrl.release" ] && n=$(cat "$STATE/pf_ctrl.release")
  rssh "condor_q -hold -af ClusterId ProcId JobBatchName HoldReasonCode HoldReasonSubCode HoldReason 2>/dev/null | head -5" | sed 's/^/[held] /' | tee -a "$LOG" >/dev/null
  if [ "$n" -lt 1 ] && rssh "condor_q -hold -af JobBatchName 2>/dev/null | grep -q ."; then
    echo $((n+1)) > "$STATE/pf_ctrl.release"; log "[pf_ctrl] condor_release klyczak (jednorazowo)"
    rssh "condor_release klyczak 2>&1 | tail -2" | tee -a "$LOG" >/dev/null
  fi; }
stage_done()   { rssh "condor_q -batch 2>/dev/null | grep -E 'charm_pfctrl|charm_s5|charm_v|charm_s4|charm_ev'; condor_q -totals 2>/dev/null | tail -1" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null; held_fix; rssh "[ -f ~/$R/.orch_start ] && [ -f $EOSD/DONE ] && [ $EOSD/DONE -nt ~/$R/.orch_start ] && [ -n \"\$(find $EOSD/DONE -newermt '2026-11-05 03:45' 2>/dev/null)\" ]"; }   # 05.11: DONE sprzed poprawki ID śladu (5.11) nieważny
stage_fetch()  { mkdir -p "$DST" && rpull "${mdjson[@]}" "$H:$EOSD/" "$DST/"; }
stage_result() { echo "wyniki → charm_work/PF3_ctrl/runs/lxplus/compare_ctrl.md (bramka G6)"; }
