# J12EFF (17.11.2026, projekt A): wydajność L1_J12 z DANYCH (oba runy) w funkcji wiodącego dżetu torowego / okna torowego / toru wiodącego,
#   względem zdarzeń MBTS; histogramy zdarzeń, ŻADNYCH par D–D. Kalibruje emulację triggera w etapie trigmc.
#   Paczka charm_work/TRIG_A/j12eff (+ tj_kernel.py z Q_lxplus/trkjet, trg_root_extract.py i btrig_common.py z B_trig) → Q_lxplus/j12eff.
source "$HERE/stages/_common.sh"
SRC=$W/TRIG_A/j12eff; R=$RQ/j12eff; EOSD=$EOSQ/j12eff; DST=$W/TRIG_A/j12eff/runs/lxplus
stage_available() { [ -f "$SRC/submit_j12.sh" ] && [ -f "$SRC/j12_extract.py" ] && [ -f "$QL/trkjet/tj_kernel.py" ] && [ -f "$W/B_trig/trg_root_extract.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/j12eff && touch ~/$R/.orch_start" \
                 && rpush --exclude 'runs/' --exclude 'j12.dag*' --exclude '*.bak*' "$SRC/" "$H:$R/" \
                 && rpush --include 'tj_kernel.py' --exclude '*' "$QL/trkjet/" "$H:$R/" \
                 && rpush --include 'btrig_common.py' --include 'trg_root_extract.py' --exclude '*' "$W/B_trig/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && ./submit_j12.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_j12.sh submit" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "cd ~/$R 2>/dev/null && ./submit_j12.sh status 2>/dev/null | head -3; echo '## dagman.out'; tail -n 6 j12.dag.dagman.out 2>/dev/null; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/j12eff/*.err 2>/dev/null | head -2); do echo \"### \$f\"; tail -n 12 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rssh "[ -f ~/$R/.orch_start ] && [ -f $EOSD/DONE ] && [ $EOSD/DONE -nt ~/$R/.orch_start ]"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.md' --include '*.json' --include 'j12eff_*.npz' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "J12EFF: $DST/J12EFF.md + j12eff.json (ε L1_J12 vs dżet torowy, dane obu runów; bez par)"; }
