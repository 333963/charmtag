# RUCIO_V20 (16.11.2026): v1.1 na gridzie wygasło (kontenery puste) → pobranie ROOT v2.0 (CharmAnalysis Oct28, produkcja ntupli v2.0) z SCRATCHDISK
#   (też wygaśnie) na EOS: data 428770 (32 GB), mb ND (27 GB), data 427929 (41 GB) → /eos/user/k/klyczak/charm_raw_v20/<nazwa>.
#   Proxy VOMS ~/private/x509up_charm (96 h od 16.11 21:2x). Wejście dla gen_norm (RAW → charm_raw_v20, format v2.0 do sprawdzenia).
source "$HERE/stages/_common.sh"
SRC=$QL/rucio_v20; R=$RQ/rucio_v20; EOSD=/eos/user/k/klyczak/charm_raw_v20; DST=$W/Q_lxplus/rucio_v20/status
stage_available() { [ -f "$SRC/submit_rv.sh" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/rucio_v20" && rpush --delete "$SRC/" --exclude status "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh && ./submit_rv.sh check" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_rv.sh submit && ./submit_rv.sh status" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "cd ~/$R && ./submit_rv.sh status; for f in ~/$RQ/logs/rucio_v20/*.out; do echo \"### \$f\"; tail -n 5 \$f; done" > "$DST/status.txt" 2>&1
  sed 's/^/[rv] /' "$DST/status.txt" | grep -E 'plików|charm_rv' | tee -a "$LOG" >/dev/null
  rssh "for n in data_428770 mb_nd data_427929 fd_nd mb_sd mb_dd mb_cd; do test -f $EOSD/\$n/DONE || exit 1; done"; }
stage_fetch()  { :; }
stage_result() { echo "RUCIO_V20: ROOT v2.0 na EOS $EOSD/{data_428770,mb_nd,data_427929} (status: Q_lxplus/rucio_v20/status/status.txt). Dalej: gen_norm RAW → v2.0 (format do sprawdzenia)."; }
