# UB427JK (19.11.2026, projekt A; pary w danych zatwierdzone 18.11): jackknife μ po 200 grupach zdarzeń (hash(run, event)) dla runu 427929
#   na wejściach UB427 (charm_q/ub_427929/pairs) — czynnik σ_jk/σ_profil do Δμ(427929 − 428770) w TRIG_A/syst_trig.py (prerejestracja TRIG_A/SYST_TRIG.md).
#   Paczka Q_lxplus/ub427_jk (+ ub_jk/ubjk_run.py bez zmian, ub_427929/u4_common.sh). Prep przebudowuje nogi 427929 (cmd_pairs --unblind, 16 GB).
#   Wynik → PF2_lxplus/runs/lxplus/unblind/r427929/jk_ev200/ubjk.md. NIE dopisany do pipeline.yaml.
source "$HERE/stages/_common.sh"
SRC=$QL/ub427_jk; R=$RQ/ub427_jk; EOSD=$EOSQ/ub_427929/jk_ev200; DST=$W/PF2_lxplus/runs/lxplus/unblind/r427929/jk_ev200
stage_available() { [ -f "$SRC/submit_u7.sh" ] && [ -f "$SRC/u7_job.sh" ] && [ -f "$QL/ub_jk/ubjk_run.py" ] && [ -f "$QL/ub_427929/u4_common.sh" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/ub427jk" && rpush --exclude '__pycache__' --exclude '*.pyc' --exclude '*.dag*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && FORCE=${FORCE:-0} ./submit_u7.sh check" | tail -10 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && FORCE=${FORCE:-0} ./submit_u7.sh submit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_ub427jk'; echo '## dagman'; tail -n 8 ~/$R/u7.dag.dagman.out 2>/dev/null; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/ub427jk/*.err 2>/dev/null | head -3); do echo \"### \$f\"; tail -n 15 \$f; done; echo \"fit: \$(ls $EOSD 2>/dev/null | grep -c '^fit_')\"" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include 'ubjk.md' --include 'ubjk.json' --include 'perfile.json' --include 'DONE' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "UB427JK: $DST/ubjk.md (σ_jk/σ dla 427929; dalej TRIG_A/syst_trig.py --ub427jk $DST/ubjk.json)"; }
