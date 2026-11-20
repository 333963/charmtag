# UB427 (18.11.2026, projekt A; pary w danych zatwierdzone 18.11): μ par D–D w czystym runie minimum-bias 427929 (klasa MB: GRL ∧
#   HLT_mb_mbts_L1MBTS_1) — kontrola nieobciążona triggerem dla μ = 1,50 z 428770. h5 v2.0 po cięciu L_xy(D0, D*) ≥ 0,25 mm i oknie D0 (≡ v1.1),
#   te same modele/progi/szablony/fit; N_D* z 427929; c_w = c_w(nom) × c(427929)/c(428770) z wag w2. Paczka Q_lxplus/ub_427929 (+ pf_pairs, ub_jk).
#   Runów nie łączy. Wynik → PF2_lxplus/runs/lxplus/unblind/r427929/UB427.md.
source "$HERE/stages/_common.sh"
SRC=$QL/ub_427929; R=$RQ/ub_427929; EOSD=$EOSQ/ub_427929; DST=$W/PF2_lxplus/runs/lxplus/unblind/r427929
stage_available() { [ -f "$SRC/submit_u4.sh" ] && [ -f "$SRC/u4_report.py" ] && [ -d "$QL/ub_jk" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/ub427" && rpush --exclude '__pycache__' --exclude '*.pyc' --exclude '*.dag*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && FORCE=${FORCE:-0} ./submit_u4.sh check" | tail -14 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && FORCE=${FORCE:-0} ./submit_u4.sh submit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_ub427'; echo '## dagman'; tail -n 8 ~/$R/u4.dag.dagman.out 2>/dev/null; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/ub427/*.err 2>/dev/null | head -3); do echo \"### \$f\"; tail -n 15 \$f; done; for s in h5 gn/mc gn/d428 gn/d427 cw scores; do echo \"\$s: \$(ls $EOSD/\$s 2>/dev/null | wc -l)\"; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include 'UB427.md' --include 'ub427.json' --include 'ndstar_427929.json' --include 'cw/' --include 'cw/cw.json' --include 'pairs/' --include 'pairs/*.json' --include 'pairs/*.md' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "UB427: $DST/UB427.md (μ par w 427929 vs 428770)"; }
