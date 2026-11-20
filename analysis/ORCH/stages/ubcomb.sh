# UBCOMB (19.11.2026, projekt A; zatwierdzone 18.11): wspólne μ par D–D z próbek bez triggera dżetowego — 427929 (MB) + podzbiory 428770
#   (strict = łańcuchy bez obciążenia; nonj12 = klasy UBTRIG MBTS_noj20 ⊕ other; nonj12_noafp), osobne N_D*/c_w/tło, wspólne μ (suma profili).
#   Paczka Q_lxplus/ub_comb (+ ub_jk: blokada UNBLIND_OK i env; wejścia z UBTRIG trg, B_trig bt, UB427 pairs/ndstar, opcjonalnie UB427JK).
#   Wynik → PF2_lxplus/runs/lxplus/unblind/ubcomb/UBCOMB.md
source "$HERE/stages/_common.sh"
SRC=$QL/ub_comb; R=$RQ/ub_comb; EOSD=$EOSQ/ub_comb; DST=$W/PF2_lxplus/runs/lxplus/unblind/ubcomb
stage_available() { [ -f "$SRC/submit_uc.sh" ] && [ -f "$SRC/ubcomb_run.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/ubcomb" && rpush --exclude '__pycache__' --exclude '*.pyc' --exclude '*.dag*' --exclude '*.bak*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && ./submit_uc.sh check" | tail -10 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_uc.sh submit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_ubcomb'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/ubcomb/*.err ~/$RQ/logs/ubcomb/*.out 2>/dev/null | head -2); do echo \"### \$f\"; tail -n 20 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.md' --include 'ubcomb.json' --include 'ubcomb.log' --include 'DONE' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "UBCOMB: $DST/UBCOMB.md (wspólne μ: 427929 + 428770 bez J12)"; }
