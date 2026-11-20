# UBTRIGNULL (19.11.2026, projekt A; zatwierdzone 18.11, wysyłka 19.11): test zerowy UBTRIG — klasy zdarzeń 428770 przetasowane losowo
#   w warstwach (nogi × D*), te same liczności co klasy triggera, 24 ziarna → czy małe klasy dają obciążone μ (model tła), czy efekt J12 jest
#   z triggera. Paczka Q_lxplus/ub_trig (kroki null/nullagg). Wynik → PF2_lxplus/runs/lxplus/unblind/ubtrig/null.
source "$HERE/stages/_common.sh"
SRC=$QL/ub_trig; R=$RQ/ub_trig; EOSD=$EOSQ/ub_trig/null; DST=$W/PF2_lxplus/runs/lxplus/unblind/ubtrig/null
stage_available() { [ -f "$SRC/ubtrig_null_agg.py" ] && [ -f "$SRC/submit_ut.sh" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/ubtrig" && rpush --exclude '__pycache__' --exclude '*.pyc' --exclude '*.dag*' --exclude '*.bak*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && source ut_common.sh && test \$(ls \$UT_OUT/trg/trg_*.npz | wc -l) -eq \$UT_NG && test -s \$UT_OUT/ubtrig.json && ! condor_q -batch 2>/dev/null | grep -q charm_ubtrnull && echo 'OK — test zerowy można wysłać'" | tail -4 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_ut.sh nullsubmit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_ubtrnull'; ls $EOSD 2>/dev/null | head -40" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/NULL_DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include 'UBTRIG_NULL.md' --include 'ubtrig_null.json' --include 'NULL_DONE' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "UBTRIGNULL: $DST/UBTRIG_NULL.md (test zerowy klas triggera, 24 ziarna)"; }
