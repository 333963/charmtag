# TUNE (10.11): warianty tune Pythii (MPIup/down, CRup/down, CRqcd, Monash; 6 × 100 × 1 M zdarzeń) → combined_syst_<t>/truth_norm.json
#   → systematyka normalizacji μ (wiersze „CZEKA” w syst_blind.md). Paczka Q_lxplus/syst (make_dag_tune.sh, submit_tune.sh); tylko gen+combine.
#   Lokalnie po fetch: ../.venv_h5/bin/python TUNE_pythia/tune_syst.py
source "$HERE/stages/_common.sh"
SRC=$QL/syst; R=$RQ/syst; EOSD=$EOSQ/truth; DST=$CH/runs/lxplus/truth
stage_available() { [ -f "$SRC/submit_tune.sh" ] && [ -f "$SRC/make_dag_tune.sh" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/syst" && rpush --exclude 'bin/' --exclude '*.dag*' --exclude '*.bak*' --exclude '__pycache__' --exclude '*.pyc' "$SRC/" "$H:$R/"; }
# Jednorazowo (zatwierdzone 10.11): usuń martwy syst.dag z 03.11 (klaster 12143398; PRE-skrypty w pętli DEFER od tygodnia, 0 jobów).
# Tylko gdy to nadal condor_dagman z syst.dag i NIE ma żadnych jobów potomnych; flaga state/tune.zombie_rm → nigdy drugi raz.
zombie_rm() {
  [ -f "$STATE/tune.zombie_rm" ] && return 0
  rssh 'c=12143398; cmd=$(condor_q $c -af Cmd JobBatchName 2>/dev/null | head -1); kids=$(condor_q -constraint "DAGManJobId == $c" -af ClusterId 2>/dev/null | wc -l)
        if [ -z "$cmd" ]; then echo "zombie_rm: klastra $c nie ma w kolejce — nic do zrobienia"; exit 0; fi
        case "$cmd" in *condor_dagman*" syst.dag+12143398") ;; *) echo "zombie_rm: $c to nie syst.dag ($cmd) — NIE usuwam"; exit 3 ;; esac
        if [ "$kids" -gt 0 ]; then echo "zombie_rm: $c ma $kids jobów potomnych — NIE usuwam"; exit 3; fi
        condor_rm $c 2>&1 | tail -1; sleep 8; condor_q -batch 2>/dev/null | grep -c "syst.dag" | sed "s/^/zombie_rm: syst.dag w kolejce po rm: /"' 2>&1 | sed 's/^/[tune] /' | tee -a "$LOG" >/dev/null
  local rc=${PIPESTATUS[0]}
  [ "$rc" = 0 ] && date +%s > "$STATE/tune.zombie_rm"
  return 0
}
stage_check()  { zombie_rm; rssh "cd ~/$RQ && chmod +x syst/*.sh && ./syst/submit_tune.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$RQ && ./syst/submit_tune.sh submit" | tail -14 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
# stage_done (12.11): jedno ssh = log kolejki + decyzja. Gotowe (rc 0) gdy tune_DONE na EOS, ALBO brak charm_tune w kolejce
#   i ≥ 5 combined_syst_*/truth_norm.json (CRqcd może nie dojść). condor_q z błędem → rc 2 (nie rozstrzygamy, nie pobieramy).
stage_done()   { rssh "q=\$(condor_q -batch 2>/dev/null) || { echo 'tune: condor_q błąd — nie rozstrzygam'; exit 2; }
    echo \"\$q\" | grep -E 'charm_tune'; condor_q -totals 2>/dev/null | tail -1
    (cd ~/$RQ && ./syst/submit_tune.sh status 2>/dev/null | grep -E 'gen |DAG status' | tail -7)
    n=\$(ls $EOSD/combined_syst_*/truth_norm.json 2>/dev/null | wc -l)
    if [ -f $EOSD/tune_DONE ]; then echo \"tune: tune_DONE jest, combined \$n/6\"; exit 0; fi
    if echo \"\$q\" | grep -q 'charm_tune'; then echo \"tune: charm_tune w kolejce, combined \$n/6\"; exit 1; fi
    if [ \"\$n\" -ge 5 ]; then echo \"tune: kolejka pusta, combined \$n/6 (≥5) — pobieram częściowe\"; exit 0; fi
    echo \"tune: kolejka pusta, combined tylko \$n/6 — nie gotowe\"; exit 1" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
  return "${PIPESTATUS[0]}"; }
stage_fetch()  { mkdir -p "$DST" && rpull --prune-empty-dirs --include 'combined_syst_*/' --include 'combined_syst_*/truth_norm.json' --include 'combined_syst_*/truth_norm.md' --include 'tune_DONE' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { local n; n=$(ls -d "$DST"/combined_syst_*/ 2>/dev/null | wc -l | tr -d ' ')
  if [ "$n" -lt 6 ]; then echo "TUNE: UWAGA pobrano tylko $n/6 wariantów (combined_syst_*/truth_norm.json) → $DST (analiza: charm_work/TUNE_pythia/tune_syst.py)"
  else echo "TUNE: pobrano $n/6 wariantów combined_syst_*/truth_norm.json → $DST (analiza: charm_work/TUNE_pythia/tune_syst.py)"; fi; }
