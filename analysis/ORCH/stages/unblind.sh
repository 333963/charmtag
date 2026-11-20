# UNBLIND (T-PREP): odsłonięcie par — TYLKO po decyzji o odsłonięciu. Dwa zamki: (1) plik-flaga ORCH/UNBLIND_OK na Macu (tworzony ręcznie:
#   echo "decyzja ODSŁOŃ $(date) <imię>" > charm_work/ORCH/UNBLIND_OK), (2) submit_unblind.sh na lxplus sprawdza pf DONE, brak unblind/DONE,
#   hash configu; unblind_job.sh odmawia bez $UB_OUT/UNBLIND_OK. Wyniki do $PF_OUT/unblind (ślepe pliki pf nietknięte).
# UWAGA: pipeline.sh parsuje pipeline.yaml raz przy starcie → ten etap działa dopiero po ponownym uruchomieniu pętli
#   (albo ręcznie: dwuklik START_lxplus_unblind.command). Nie restartować pętli, gdy pf_pairs biegnie.
source "$HERE/stages/_common.sh"
SRC=$QL/pf_pairs; R=$RQ/pf_pairs; EOSD=$EOSQ/pf_pairs/unblind; DST=$W/PF2_lxplus/runs/lxplus/unblind; FLAG=$HERE/UNBLIND_OK
stage_available() { [ -f "$SRC/submit_unblind.sh" ] && [ -s "$FLAG" ]; }     # bez flagi etap jest pomijany (skipped), nie failed
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/pf" && cp "$FLAG" "$SRC/UNBLIND_OK" && rpush --exclude src "$SRC/" "$H:$R/" && rpush --delete --exclude __pycache__ "$SRC/src/" "$H:$R/src/"; }   # 05.11: src = 2.0.2 (dstar_refit N_D*, ślepy zapis tylko SB) — ten sam kod co pf_ctrl (gałąź d0z0 = ślepe odniesienie)
stage_check()  { rssh "cd ~/$R && chmod +x *.sh tools/*.py && ./submit_unblind.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_unblind.sh submit; ./submit_unblind.sh status | tail -5" | tail -12 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { rssh "condor_q -batch 2>/dev/null | grep -E 'charm_unblind|unblind'; condor_q -totals 2>/dev/null | tail -1; tail -2 ~/$R/unblind.dag.dagman.out 2>/dev/null | cut -c1-120; ls $EOSD 2>/dev/null | tr '\n' ' '" | sed 's/^/[uq] /' | tee -a "$LOG" >/dev/null; rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include 'UNBLIND_OK' --include 'pairs_full.npz' --include 'pair_inputs.npz' --include '*.png' "${mdjson[@]}" "$H:$EOSD/" "$DST/"; }   # include przed mdjson (rsync: pierwsza pasująca reguła; mdjson kończy się --exclude '*')
stage_result() { echo "odsłonięte: $DST/tfin.md (kopia EOS $EOSD/)"; }
