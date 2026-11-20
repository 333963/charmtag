# UBNDSW (19.11.2026, projekt A; kontrole luk zatwierdzone 19.11): skala σ piku D* w klasach triggera 428770 (MBTS / J12only / other)
#   vs 427929 — czy szeroki pik D* w 427929 (1,40 × MC) jest prawdziwy. Kryterium: Q_lxplus/ub_ndsw/README.md (zapisane przed wynikiem).
#   Paczka Q_lxplus/ub_ndsw (+ ub_trig: load_trg, lookup_cls; wejścia tylko do odczytu). Wynik → PF2_lxplus/runs/lxplus/unblind/r427929/ndsw.
source "$HERE/stages/_common.sh"
SRC=$QL/ub_ndsw; R=$RQ/ub_ndsw; EOSD=$EOSQ/ub_ndsw; DST=$W/PF2_lxplus/runs/lxplus/unblind/r427929/ndsw
stage_available() { [ -f "$SRC/submit_nw.sh" ] && [ -f "$SRC/ndsw.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/ubndsw" && rpush --exclude '__pycache__' --exclude '*.pyc' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && ./submit_nw.sh check" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_nw.sh submit" | tail -4 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_ubndsw'; for f in \$(ls -t ~/$RQ/logs/ubndsw/*.err ~/$RQ/logs/ubndsw/*.out 2>/dev/null | head -2); do echo \"### \$f\"; tail -n 20 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include 'NDSW.md' --include 'ndsw.json' --include 'ndsw_hists.npz' --include 'nw.log' --include 'DONE' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "UBNDSW: $DST/NDSW.md (skala σ D* w klasach triggera 428770 vs 427929)"; }
