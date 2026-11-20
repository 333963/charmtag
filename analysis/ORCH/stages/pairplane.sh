# PAIRPLANE (19.11, projekt A, prezentacja pkt 12 „corr. excess dla danych”): płaszczyzny mas par D–D jak plots_doubleD z kodu bazowego charmpp,
#   ale z wetem wspólnych torów (trackId), OS/LS, wariantem weta D* dla D⁺ (Δm > 160 MeV), |Δφ|/Δy pik×pik vs pasma, poziom matched w MC.
#   Wejście h5 v2.0 (data_427929, data_428770, mb_nd, fd_nd); paczka Q_lxplus/pair_planes (batch charm_pp; 41 węzłów + finish).
#   Wykresy lokalnie: .venv_h5/bin/python charm_work/PLOTS_doubleD/plot_planes.py --in $DST --out $DST/out
source "$HERE/stages/_common.sh"
SRC=$QL/pair_planes; R=$RQ/pair_planes; EOSD=$EOSQ/pair_planes; DST=$W/PLOTS_doubleD/runs/lxplus
stage_available() { [ -f "$SRC/submit_pp.sh" ] && [ -f "$SRC/pp_extract.py" ] && [ -f "$SRC/pp_merge.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/pair_planes" && rpush --delete --exclude '__pycache__' --exclude '*.pyc' --exclude 'pp.dag*' --exclude '*.bak*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && ./submit_pp.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_pp.sh submit; ./submit_pp.sh status | tail -4" | tail -10 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { rssh "condor_q -batch 2>/dev/null | grep -E 'charm_pp'; condor_q -totals 2>/dev/null | tail -1" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
  mkdir -p "$DST"; rssh "cd ~/$R; echo '## dagman.out'; tail -n 20 pp.dag.dagman.out 2>/dev/null; echo '## logs'; ls -lt ~/$RQ/logs/pair_planes | head -8; for f in \$(ls -t ~/$RQ/logs/pair_planes/*.err 2>/dev/null | head -3); do echo \"### \$f\"; tail -n 20 \$f; done; echo '## EOS'; ls -la $EOSD 2>/dev/null | head -30; ls $EOSD/parts 2>/dev/null | wc -l" > "$DST/diag.txt" 2>&1
  rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --exclude 'parts/' --include '*.npz' --include 'DONE' "${mdjson[@]}" "$H:$EOSD/" "$DST/"; }
stage_result() { echo "PAIRPLANE: <próbka>.npz + meta → $DST (wykresy: .venv_h5/bin/python charm_work/PLOTS_doubleD/plot_planes.py --in $DST --out $DST/out)"; }
