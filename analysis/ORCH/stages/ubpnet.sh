# UBPNET (19.11.2026, projekt A; pary w danych zatwierdzone 18.11): μ par D–D w 428770 (odsłonięta próbka 05.11) z kategoriami
#   sieci par PAIRNET v3 (pnetd = tło SB×SB danych; kontrola xgbd/base_sum) vs nomcat w tym samym modelu, SF_rel z HYBSF/HYBDS v3.
#   Paczka Q_lxplus/ub_pnet (+ ub_jk: blokada UNBLIND_OK i env; + pairnet: kod v3). Wynik → PF2_lxplus/runs/lxplus/unblind/ubpnet.
#   Wymaga zsynchronizowanej paczki pairnet (etap pairnet/hybds) na lxplus — check to sprawdza.
source "$HERE/stages/_common.sh"
SRC=$QL/ub_pnet; R=$RQ/ub_pnet; EOSD=$EOSQ/ub_pnet; DST=$W/PF2_lxplus/runs/lxplus/unblind/ubpnet
stage_available() { [ -f "$SRC/submit_upn.sh" ] && [ -f "$SRC/upn_fit.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/ubpnet" && rpush --exclude '__pycache__' --exclude '*.pyc' --exclude '*.dag*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && ./submit_upn.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_upn.sh submit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_ubpnet'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/ubpnet/*.err ~/$RQ/logs/ubpnet/*.out 2>/dev/null | head -4); do echo \"### \$f\"; tail -n 15 \$f; done; ls $EOSD/parts 2>/dev/null | grep -c npz" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.md' --include 'ubpnet.json' --include 'ubpnet.log' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "UBPNET: $DST/UBPNET.md (μ z kategoriami sieci par vs nomcat, 428770)"; }
