# PAIRFIT (17.11.2026, projekt A): oczekiwany błąd μ (fit Asimova 2D mas, informacja Fishera) dla kategorii score pary vs nominalne kategorie nóg.
#   TYLKO MC (nogi i score par z etapu pairnet); normalizacja do istniejących liczb nominalnych (n_pair, Y, Σ sbsb z pair_inputs.json) — bez nowych obliczeń na danych.
#   Paczka charm_work/PAIRNET (pn_fit.py, submit_pnfit.sh) → Q_lxplus/pairnet; wynik EOS charm_q/pairnet/fit → PAIRNET/runs/lxplus/fit/PAIRFIT_K*.md
source "$HERE/stages/_common.sh"
SRC=$W/PAIRNET; R=$RQ/pairnet; EOSD=$EOSQ/pairnet/fit_v2; DST=$W/PAIRNET/runs/lxplus/fit_v2
stage_available() { [ -f "$SRC/submit_pnfit.sh" ] && [ -f "$SRC/pn_fit.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/pairnet" && rpush --exclude 'runs/' --exclude '__pycache__' --exclude '*.pyc' --exclude 'pn.dag*' --exclude 'lists/' --exclude '*.bak*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && FORCE=${FORCE:-0} ./submit_pnfit.sh check" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && FORCE=${FORCE:-0} ./submit_pnfit.sh submit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_pairfit'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/pairnet/pn_fit.*.err ~/$RQ/logs/pairnet/pn_fit.*.out 2>/dev/null | head -2); do echo \"### \$f\"; tail -n 15 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.md' --include '*.json' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "PAIRFIT: $DST/PAIRFIT_K6.md, PAIRFIT_K12.md (σ_μ Asimova: kategorie score pary vs kategorie nóg; tylko MC)"; }
