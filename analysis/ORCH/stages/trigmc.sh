# TRIGMC (17.11.2026, projekt A): zamknięcie wag w2(n_trk, lead_pt) wobec emulacji triggera L1_J12 na poziomie generatora.
#   Pythia 8 A3 ND, 200 × 1 M zdarzeń, z oknem kalorymetrycznym 0,8×0,8 i dżetami naładowanymi; combine → κ = [pary/D*]_wyzwolone / [pary/D*]_w2.
#   Bez danych (tylko generator). Paczka charm_work/TRIG_A/trigmc → ~/private/charm_q/Q_lxplus/trigmc; wynik → TRIG_A/trigmc/runs/lxplus.
source "$HERE/stages/_common.sh"
SRC=$W/TRIG_A/trigmc; R=$RQ/trigmc; EOSD=$EOSQ/trigmc; DST=$W/TRIG_A/trigmc/runs/lxplus
stage_available() { [ -f "$SRC/submit_tm.sh" ] && [ -f "$SRC/pythia_dd_jet.cc" ] && [ -f "$SRC/tm_combine.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/trigmc && touch ~/$R/.orch_start" \
                 && rpush --exclude 'bin/' --exclude 'runs/' --exclude 'tm.dag*' --exclude '*.bak*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && ./submit_tm.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_tm.sh submit" | tail -12 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "cd ~/$R 2>/dev/null && ./submit_tm.sh status 2>/dev/null | head -4; echo '## dagman.out'; tail -n 8 tm.dag.dagman.out 2>/dev/null; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/trigmc/*.err 2>/dev/null | head -2); do echo \"### \$f\"; tail -n 12 \$f; done" > "$DST/diag.txt" 2>&1
                 head -2 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rssh "[ -f ~/$R/.orch_start ] && [ -f $EOSD/DONE ] && [ $EOSD/DONE -nt ~/$R/.orch_start ]"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.md' --include '*.json' --exclude '*' "$H:$EOSD/res/" "$DST/"; }
stage_result() { echo "TRIGMC: $DST/TRIGMC.md + trigmc.json (κ zamknięcia w2 per emulacja triggera; generator, bez danych)"; }
