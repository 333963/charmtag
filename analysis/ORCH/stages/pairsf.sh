# PAIRSF (17.11.2026, projekt A): SF klasyfikatora H dla D0 z tagu D* w kategoriach score używanych przez pary (c0/c1/c2/c4, ≥WP65, ≥WP90),
#   per klasa triggera (428770: all / J12only / MBTSbased / other) i run 427929 (MBTS). Pojedyncze D — ŻADNYCH par D–D.
#   Wejście: wyniki b_trig na EOS (tylko odczyt). Paczka charm_work/TRIG_A/pairsf (+ btrig_common.py z B_trig, src/charmtag z AA_merge)
#   → ~/private/charm_q/Q_lxplus/pairsf; wynik → TRIG_A/pairsf/runs/lxplus/{PSF.md,psf.json}.
source "$HERE/stages/_common.sh"
SRC=$W/TRIG_A/pairsf; R=$RQ/pairsf; EOSD=$EOSQ/pairsf; DST=$W/TRIG_A/pairsf/runs/lxplus
stage_available() { [ -f "$SRC/submit_psf.sh" ] && [ -f "$SRC/psf_analyze.py" ] && [ -f "$W/B_trig/btrig_common.py" ] && [ -f "$W/AA_merge/src/charmtag/dstar_refit.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R/src ~/$RQ/logs/pairsf && touch ~/$R/.orch_start" \
                 && rpush --exclude 'runs/' --exclude '*.bak*' "$SRC/" "$H:$R/" \
                 && rpush --include 'btrig_common.py' --exclude '*' "$W/B_trig/" "$H:$R/" \
                 && rpush --exclude '*.bak*' --exclude '*.pyc' "$W/AA_merge/src/charmtag" "$H:$R/src/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && ./submit_psf.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_psf.sh submit" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_pairsf'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/pairsf/*.err ~/$RQ/logs/pairsf/*.out 2>/dev/null | head -2); do echo \"### \$f\"; tail -n 15 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rssh "[ -f ~/$R/.orch_start ] && [ -f $EOSD/DONE ] && [ $EOSD/DONE -nt ~/$R/.orch_start ]"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.md' --include '*.json' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "PAIRSF: $DST/PSF.md + psf.json (SF w kategoriach score par per klasa triggera; pojedyncze D0←D*)"; }
