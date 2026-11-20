# UBTRIG (18.11.2026, projekt A; pary w danych zatwierdzone 18.11): μ par D–D per rozłączna klasa triggera w 428770
#   (J12only / MBTS_j20 / MBTS_noj20 / other) na odsłoniętej próbce 05.11 — ten sam model/fit, N_D* per klasa z fitu Δm.
#   Paczka Q_lxplus/ub_trig (+ ub_jk: blokada UNBLIND_OK i env). Wynik → PF2_lxplus/runs/lxplus/unblind/ubtrig.
source "$HERE/stages/_common.sh"
SRC=$QL/ub_trig; R=$RQ/ub_trig; EOSD=$EOSQ/ub_trig; DST=$W/PF2_lxplus/runs/lxplus/unblind/ubtrig
stage_available() { [ -f "$SRC/submit_ut.sh" ] && [ -f "$SRC/ubtrig_run.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/ubtrig" && rpush --exclude '__pycache__' --exclude '*.pyc' --exclude '*.dag*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && ./submit_ut.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_ut.sh submit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_ubtrig'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/ubtrig/*.err ~/$RQ/logs/ubtrig/*.out 2>/dev/null | head -4); do echo \"### \$f\"; tail -n 15 \$f; done; ls $EOSD/trg 2>/dev/null | wc -l" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.md' --include 'ubtrig.json' --include 'ubtrig.log' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "UBTRIG: $DST/UBTRIG.md (μ per klasa triggera, 428770)"; }
