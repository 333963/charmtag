# UBPNET427 (19.11.2026, projekt A; pary w danych zatwierdzone 18.11): sieć par PAIRNET v3 (pnetd, bez retrenowania; tło SBTRAIN z 428770)
#   na czystym runie 427929 (klasa MB, filtr UB427) vs nomcat w tym samym modelu; SF nóg HYBSF/HYBDS data_427929; n_pair × f_ch (prawda mb_nd);
#   test transferu tła (pary SB×SB 427929 vs 428770). Paczka Q_lxplus/ub_pnet (pliki upn427_*; bieg 428770 bez zmian).
#   Wynik → PF2_lxplus/runs/lxplus/unblind/r427929/pnet. Wymaga etapów ub427 (h5 + pary) i pairnet/hybds na EOS — check sprawdza.
source "$HERE/stages/_common.sh"
SRC=$QL/ub_pnet; R=$RQ/ub_pnet; EOSD=$EOSQ/ub_pnet427; DST=$W/PF2_lxplus/runs/lxplus/unblind/r427929/pnet
stage_available() { [ -f "$SRC/submit_upn427.sh" ] && [ -f "$SRC/upn427_job.sh" ] && [ -f "$SRC/upn_fit.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/ubpnet427" && rpush --exclude '__pycache__' --exclude '*.pyc' --exclude '*.dag*' --exclude '*.bak*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && ./submit_upn427.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_upn427.sh submit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_ubpnet427'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/ubpnet427/*.err ~/$RQ/logs/ubpnet427/*.out 2>/dev/null | head -4); do echo \"### \$f\"; tail -n 15 \$f; done; ls $EOSD/parts 2>/dev/null | grep -c npz" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include 'UBPNET427.md' --include 'ubpnet.json' --include 'ubpnet.log' --include 'DONE' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "UBPNET427: $DST/UBPNET427.md (sieć par na 427929: μ pnetd vs nomcat, transfer tła, SF syst)"; }
