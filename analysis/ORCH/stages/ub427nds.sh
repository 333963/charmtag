# UB427NDS (19.11.2026, projekt A, studium W2 „427929 jako pomiar główny”): systematyka N_D* w 427929 — skan wariantów NDS jak dla 428770
#   (tło, zakres, binowanie, sygnał G+CB / skala σ stała, okno D0) + fity per bin pT ze stałą skalą σ (w 427929 skala wolna = 1,40 vs 1,02 w 428770).
#   Paczka Q_lxplus/ub427_nds (wejście tylko do odczytu: charm_q/ub_427929/h5). Wynik → PF2_lxplus/runs/lxplus/unblind/r427929/nds.
source "$HERE/stages/_common.sh"
SRC=$QL/ub427_nds; R=$RQ/ub427_nds; EOSD=$EOSQ/ub427_nds; DST=$W/PF2_lxplus/runs/lxplus/unblind/r427929/nds
stage_available() { [ -f "$SRC/submit_n8.sh" ] && [ -f "$SRC/n8_ndstar.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/ub427nds" && rpush --exclude '__pycache__' --exclude '*.pyc' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && ./submit_n8.sh check" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_n8.sh submit" | tail -4 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_ub427nds'; for f in \$(ls -t ~/$RQ/logs/ub427nds/*.err ~/$RQ/logs/ub427nds/*.out 2>/dev/null | head -2); do echo \"### \$f\"; tail -n 20 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include 'NDS_427929.md' --include 'ndstar_variants_427929.json' --include 'dstar_hist_427929.npz' --include 'n8.log' --include 'DONE' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "UB427NDS: $DST/NDS_427929.md (systematyka N_D* 427929)"; }
