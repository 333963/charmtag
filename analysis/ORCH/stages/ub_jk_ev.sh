# UBJK-EV (12.11): jackknife μ po 200 grupach zdarzeń (hash) na ODSŁONIĘTEJ próbce — ta sama paczka Q_lxplus/ub_jk co etap ub_jk,
#   z UJ_GROUPS=events:200 UJ_MODES=fixed UJ_NCHUNK=25 (precyzja czynnika ~±5 %; ~7–13 h CPU na 25 węzłach). Po ub_jk (wspólny ubjk.dag
#   i batch charm_ubjk → nie równolegle). Wyniki: EOS pf_pairs/unblind/jk_ev200/ → lokalnie unblind/jk_ev200. Zamki jak ub_jk.
source "$HERE/stages/_common.sh"
SRC=$QL/ub_jk; R=$RQ/ub_jk; EOSD=$EOSQ/pf_pairs/unblind/jk_ev200; DST=$W/PF2_lxplus/runs/lxplus/unblind/jk_ev200; FLAG=$HERE/UNBLIND_OK; UJE="UJ_GROUPS=events:200 UJ_MODES=fixed UJ_NCHUNK=25"
stage_available() { [ -f "$SRC/submit_ubjk.sh" ] && [ -f "$SRC/ubjk_run.py" ] && [ -s "$FLAG" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/ubjk" && rpush --delete --exclude '__pycache__' --exclude '*.pyc' --exclude 'ubjk.dag*' --exclude '*.bak*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && $UJE ./submit_ubjk.sh check" | tail -10 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && $UJE ./submit_ubjk.sh submit; $UJE ./submit_ubjk.sh status | tail -4" | tail -12 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { rssh "condor_q -batch 2>/dev/null | grep -E 'charm_ubjk'; condor_q -totals 2>/dev/null | tail -1" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
  mkdir -p "$DST"; rssh "cd ~/$R; echo '## dagman.out'; tail -n 20 ubjk.dag.dagman.out 2>/dev/null; echo '## logs/ubjk'; ls -lt ~/$RQ/logs/ubjk | head -10; for f in \$(ls -t ~/$RQ/logs/ubjk/*.err 2>/dev/null | head -3); do echo \"### \$f\"; tail -n 20 \$f; done; echo '## EOS'; ls -la $EOSD 2>/dev/null | head -20" > "$DST/diag.txt" 2>&1
  rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include 'logs/' --include 'logs/*.log' "${mdjson[@]}" "$H:$EOSD/" "$DST/"; }   # md/json/DONE + logi; bez perfile.npz/parquet (zostają na EOS)
stage_result() { echo "UBJK-EV: $DST/ubjk.md (200 grup zdarzeń; σ_jk/σ_Hesse per kanał + wspólne μ; EOS $EOSD/)"; }
