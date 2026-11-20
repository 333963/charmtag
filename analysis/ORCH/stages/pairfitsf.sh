# PAIRFITSF (17.11.2026, projekt A): fit Asimova σ_μ na parach z MC (pairnet_v2) z sygnałem przeważonym SF kategorii score nogi D+ z dpcal2
#   (pojedyncze D+ w danych vs MC) — ile zysku kategorii (noga D+ z cechami Dalitza, score pary) przeżywa migrację zmierzoną w danych.
#   ŻADNYCH PAR Z DANYCH (pary tylko MC; z danych tylko gotowe SF pojedynczych D+). Paczka charm_work/PAIRNET → Q_lxplus/pairnet;
#   EOS charm_q/pairnet_v2/fit_sf → PAIRNET/runs/lxplus_v2/fit_sf/.
source "$HERE/stages/_common.sh"
SRC=$W/PAIRNET; R=$RQ/pairnet; EOSD=$EOSQ/pairnet_v2/fit_sf; DST=$W/PAIRNET/runs/lxplus_v2/fit_sf
stage_available() { [ -f "$SRC/submit_pnfitsf.sh" ] && [ -f "$SRC/pn_fit.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/pairnet" && rpush --exclude 'runs/' --exclude '__pycache__' --exclude '*.pyc' --exclude 'pn.dag*' --exclude 'dpcal.dag*' --exclude 'lists/' --exclude '*.bak*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && FORCE=${FORCE:-0} ./submit_pnfitsf.sh check" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && FORCE=${FORCE:-0} ./submit_pnfitsf.sh submit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_pairfitsf'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/pairnet/pn_fitsf.*.err ~/$RQ/logs/pairnet/pn_fitsf.*.out 2>/dev/null | head -2); do echo \"### \$f\"; tail -n 15 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.md' --include '*.json' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "PAIRFITSF: $DST/PAIRFIT_SF_F3.md, PAIRFIT_SF_F0.md (zysk kategorii po SF nogi D+ z danych; pary tylko MC)"; }
