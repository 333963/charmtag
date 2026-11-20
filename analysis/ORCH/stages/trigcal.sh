# TRIGCAL (17.11.2026, projekt A): emulacja L1_J12 skalibrowana na danych — ponowny combine próbki trigmc (Pythia A3 ND, 200 M zdarzeń, już na EOS)
#   z siatką drobną progów; χ² przewidywanej ε(J12 | zdarzenie z D*, pT dżetu naładowanego) względem pomiaru j12eff (427929, 5–20 GeV) → κ dla emulacji
#   zgodnych z danymi. Tylko generator + gotowa tabela wydajności zdarzeń (bez par). Wynik: EOS charm_q/trigmc/res_cal → TRIG_A/trigmc/runs/lxplus_cal/TRIGCAL.md
source "$HERE/stages/_common.sh"
SRC=$W/TRIG_A/trigmc; R=$RQ/trigmc; EOSD=$EOSQ/trigmc/res_cal; DST=$W/TRIG_A/trigmc/runs/lxplus_cal
stage_available() { [ -f "$SRC/submit_tm.sh" ] && [ -f "$SRC/tm_combine.py" ] && [ -s "$W/TRIG_A/j12eff/runs/lxplus/j12eff.json" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/trigmc" && rpush --exclude 'runs/' --exclude 'bin/' --exclude '*.bak*' --exclude 'tm.dag*' "$SRC/" "$H:$R/" \
                 && rpush --include 'j12eff.json' --exclude '*' "$W/TRIG_A/j12eff/runs/lxplus/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && ls -la j12eff.json tm_combine.py | tail -2" | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && chmod +x *.sh && FORCE=${FORCE:-0} ./submit_tm.sh calib" | tail -5 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_trigcal'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/trigmc/calib.*.err ~/$RQ/logs/trigmc/calib.*.out 2>/dev/null | head -2); do echo \"### \$f\"; tail -n 12 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.md' --include '*.json' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "TRIGCAL: $DST/TRIGCAL.md (κ dla emulacji J12 zgodnych z pomiarem j12eff)"; }
