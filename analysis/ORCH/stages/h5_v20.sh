# H5_V20 (17.11.2026): konwersja ROOT v2.0 (CharmAnalysis Oct28, dpt3to150; wejście rucio_v20 → /eos/user/k/klyczak/charm_raw_v20) → h5 w schemacie v1.2
#   + kolumny zdarzeniowe v2.0 (AFP afp_*_A/C + tag minbias, xi_cal_A/C, nSelBL, leadTrackPtBL[GeV], nPileupVertices, passGRL, lumiBlock, trigPassed/trigPS/trigTBP)
#   + dataset evtable → /eos/user/k/klyczak/charm_h5_v20/parts/<próbka>/. Paczka Q_lxplus/h5_v20 (README.md). DAG: węzeł per (próbka, kawałek ≤25 plików) + FINAL finish.
#   Gotowe = parts/DONE na EOS (pisze finish tylko przy komplecie). fetch: tylko meta/*.json (h5 zostają na EOS).
source "$HERE/stages/_common.sh"
SRC=$QL/h5_v20; R=$RQ/h5_v20; H5O=/eos/user/k/klyczak/charm_h5_v20; DONEF=$H5O/parts/DONE; DST=$SRC/status
stage_available() { [ -f "$SRC/submit_h5v20.sh" ] && [ -f "$SRC/charmpp_h5_v20.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/h5_v20" && rpush --exclude status --exclude 'lists/' --exclude 'h5v20.dag*' --exclude '*.pyc' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && ./submit_h5v20.sh check" | tail -12 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_h5v20.sh submit" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   {  # DONE na EOS → gotowe; DAG w kolejce → czekaj; pusto bez DONE → dosyłka brakujących (maks. 2 razy), potem FAIL przez max_hours
  rtest "$DONEF" && return 0
  queued charm_h5v20 && return 1
  [ "$(get h5_v20)" = running ] || return 1          # przed submitem: bez skutków ubocznych
  local n=0; [ -f "$STATE/h5_v20.retry" ] && n=$(cat "$STATE/h5_v20.retry")
  if [ "$n" -lt 2 ]; then echo $((n+1)) > "$STATE/h5_v20.retry"; log "[h5_v20] kolejka pusta, brak DONE → ./submit_h5v20.sh submit (dosyłka brakujących, próba $((n+1)))"
    rssh "cd ~/$R && ./submit_h5v20.sh status | tail -12; ./submit_h5v20.sh submit" | tail -16 | tee -a "$LOG" >/dev/null; fi
  rtest "$DONEF"
}
stage_fetch()  { mkdir -p "$DST" && rpull --include '*/' --include '*.json' --exclude '*' "$H:$H5O/meta/" "$DST/meta/"
  rssh "cd ~/$R && ./submit_h5v20.sh status" > "$DST/status.txt" 2>&1; tail -4 "$DST/status.txt" | tee -a "$LOG" >/dev/null; return 0; }
stage_result() { echo "H5_V20: h5 (v1.2 + kolumny v2.0) na EOS $H5O/parts/{data_428770,data_427929,mb_nd,fd_nd,mb_sd,mb_dd,mb_cd}; meta → Q_lxplus/h5_v20/status/meta/summary.json"; }
