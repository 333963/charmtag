# NNH (17.11.2026, projekt A): sieć neuronowa vs drzewa — MLP (3×64 wg CharmNet, 4×256) i XGB na tych samych 18 cechach H, kandydatach i foldach
#   (D0 z par D0←D*, wyniki b_trig); metryki MC, ε_B z pasma bocznego danych, SF z tagu D*. Pojedyncze D — ŻADNYCH par D–D.
#   Paczka charm_work/TRIG_A/nnh (+ btrig_common.py z B_trig, src/charmtag z AA_merge) → Q_lxplus/nnh; wynik → TRIG_A/nnh/runs/lxplus.
source "$HERE/stages/_common.sh"
SRC=$W/TRIG_A/nnh; R=$RQ/nnh; EOSD=$EOSQ/nnh; DST=$W/TRIG_A/nnh/runs/lxplus
stage_available() { [ -f "$SRC/submit_nnh.sh" ] && [ -f "$SRC/nnh_train.py" ] && [ -f "$W/B_trig/btrig_common.py" ] && [ -f "$W/AA_merge/src/charmtag/dstar_refit.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R/src ~/$RQ/logs/nnh && touch ~/$R/.orch_start" \
                 && rpush --exclude 'runs/' --exclude '*.bak*' "$SRC/" "$H:$R/" \
                 && rpush --include 'btrig_common.py' --exclude '*' "$W/B_trig/" "$H:$R/" \
                 && rpush --exclude '*.bak*' --exclude '*.pyc' "$W/AA_merge/src/charmtag" "$H:$R/src/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && ./submit_nnh.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_nnh.sh submit" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_nnh'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/nnh/*.err ~/$RQ/logs/nnh/*.out 2>/dev/null | head -2); do echo \"### \$f\"; tail -n 15 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rssh "[ -f ~/$R/.orch_start ] && [ -f $EOSD/DONE ] && [ $EOSD/DONE -nt ~/$R/.orch_start ]"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.md' --include '*.json' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "NNH: $DST/NNH.md + nnh.json (sieć vs XGB na cechach H; score na EOS charm_q/nnh)"; }
