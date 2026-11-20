# MULT (17.11.2026, etap MULT_charm): uzysk D (pojedyncze, dane vs MC) i par (TYLKO MC) vs krotność z ROOT v2.0
#   (/eos/user/k/klyczak/charm_raw_v20/{mb_nd,data_427929,data_428770}, PTPROD=dpt3to150). Paczka na lxplus = MULT_charm/lxplus/* +
#   mult_{common,extract,fit,analyze,gates}.py → ~/private/charm_q/Q_lxplus/mult_charm (batch charm_mu) → EOS charm_q/mult_charm
#   → lokalnie MULT_charm/runs/lxplus (res.tgz rozpakowany do res/).
#   BLOKADA: bez pliku ORCH/MULT_OK etap NIC nie robi na lxplus (sync/check/done to no-op) i zostaje „waiting”
#   — po `touch ORCH/MULT_OK` rusza sam w następnym obiegu (bez resetu state). Pary na DANYCH (bramka G5) nigdy w tym etapie:
#   mu_job.sh nie przekazuje --pairs-data, a stage_check odmawia, gdy paczka zawiera APPROVAL_PAIRS_DATA.txt lub wywołanie --pairs-data.
source "$HERE/stages/_common.sh"
MUL=$W/MULT_charm; SRC=$MUL/lxplus; R=$RQ/mult_charm; EOSD=$EOSQ/mult_charm; DST=$MUL/runs/lxplus
MU_PY=(mult_common.py mult_extract.py mult_fit.py mult_analyze.py mult_gates.py)
mult_ok() { [ -f "$HERE/MULT_OK" ]; }
# 19.11: REANALYZE=1 — tylko analyze na istniejących parts (poprawka fitów masy: granice kształtu z prawdy MC, d0_hash_keep)
stage_available() { [ -f "$SRC/submit_mu.sh" ] && [ -f "$SRC/mu_job.sh" ] && for f in "${MU_PY[@]}"; do [ -f "$MUL/$f" ] || return 1; done; }
stage_sync()   { mult_ok || return 0
  rssh "mkdir -p ~/$R ~/$RQ/logs/mult_charm" && rpush --delete --exclude '__pycache__' --exclude '*.pyc' --exclude 'mu.dag*' --exclude 'lists/' --exclude '*.bak*' \
    "$SRC/" "${MU_PY[@]/#/$MUL/}" "$H:$R/"; }
stage_check()  { if ! mult_ok; then log "[mult] zablokowane: brak ORCH/MULT_OK (nic nie wysyłam)"; return 1; fi
  if [ -e "$SRC/APPROVAL_PAIRS_DATA.txt" ] || grep -q -- '--pairs-data' "$SRC/mu_job.sh"; then log "[mult] STOP G5: pary na danych nie w tym etapie"; return 1; fi
  rssh "cd ~/$R && test ! -e APPROVAL_PAIRS_DATA.txt && chmod +x *.sh *.py && REANALYZE=1 ./submit_mu.sh check" | tail -10 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { mult_ok || return 1
  rssh "cd ~/$R && REANALYZE=1 ./submit_mu.sh submit; ./submit_mu.sh status | tail -5" | tail -12 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mult_ok || return 1
  rssh "condor_q -batch 2>/dev/null | grep -E 'charm_mu'; condor_q -totals 2>/dev/null | tail -1" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
  mkdir -p "$DST"; rssh "cd ~/$R; echo '## dagman.out'; tail -n 30 mu.dag.dagman.out 2>/dev/null; echo '## mu.dag.rescue'; ls mu.dag.rescue* 2>/dev/null; echo '## logs/mult_charm'; ls -lt ~/$RQ/logs/mult_charm | head -12; for f in \$(ls -t ~/$RQ/logs/mult_charm/*.err 2>/dev/null | head -3); do echo \"### \$f\"; tail -n 25 \$f; done; echo '## EOS'; ls -la $EOSD 2>/dev/null | head -20; for s in mb_nd data_427929 data_428770 fd_nd mb_sd mb_dd mb_cd; do echo \"parts \$s \$(ls $EOSD/parts/\$s 2>/dev/null | wc -l)\"; done" > "$DST/diag.txt" 2>&1
  rssh "test -f $EOSD/DONE && test -f ~/$R/mu.dag && test $EOSD/DONE -nt ~/$R/mu.dag"; }   # stare DONE (sprzed mu.dag) ignorowane
stage_fetch()  { mkdir -p "$DST" && rpull --exclude 'parts/' --include 'res.tgz' "${mdjson[@]}" "$H:$EOSD/" "$DST/" && { [ ! -f "$DST/res.tgz" ] || tar xzf "$DST/res.tgz" -C "$DST"; }; }
stage_result() { echo "MULT: wyniki → $DST/res (gates.json, r<run>_<schemat>_<estymator>_w2<w>.*); pary tylko MC — pary na danych wyłącznie po G5 (PROCEDURA MULT_charm)"; }
