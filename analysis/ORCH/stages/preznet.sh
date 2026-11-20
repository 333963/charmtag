# PREZNET (19.11.2026, projekt A; prezentacja wg kartki, pkt 3): klasyfikator pojedynczych D — BDT (XGBoost) vs sieć (MLP / Deep Sets)
#   × wejścia high-level (18 cech H nominału) vs low-level (ślady córek + L_xy). TYLKO MC (fd_nd sygnał, mb_nd tło kombinatoryczne), h5 v1.1, foldy nominału.
#   Paczka Q_lxplus/prez_net. Wynik → PREZ_2026-11-19/p3/runs (ROC, AUC, rozkłady score, ważność, rzeźbienie, TABLE.md).
source "$HERE/stages/_common.sh"
SRC=$QL/prez_net; R=$RQ/prez_net; EOSD=$EOSQ/prez_net; DST=$W/PREZ_2026-11-19/p3/runs
stage_available() { [ -f "$SRC/submit_pz.sh" ] && [ -f "$SRC/pz_train.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/preznet" && rpush --exclude '__pycache__' --exclude '*.pyc' --exclude '*.dag*' --exclude '*.bak*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && ./submit_pz.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_pz.sh submit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_preznet'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/preznet/*.err ~/$RQ/logs/preznet/*.out 2>/dev/null | head -4); do echo \"### \$f\"; tail -n 15 \$f; done; ls $EOSD/ext $EOSD/train_D0 $EOSD/train_Dplus 2>/dev/null" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.png' --include '*.pdf' --include '*.md' --include '*.json' --include '*.log' --include 'DONE' --exclude '*' "$H:$EOSD/report/" "$DST/" && rpull --include 'DONE' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "PREZNET: $DST/TABLE.md (BDT vs NN × high/low-level, pojedyncze D, MC)"; }
