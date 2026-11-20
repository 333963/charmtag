# PAIRNET2 (17.11.2026, projekt A): sieć par D–D v2 — noga D+ z cechami Dalitza (H + DLZ8; XGB, sygnał fd, tło = dalekie SB POJEDYNCZYCH D+ w danych
#   428770 z tabel DNET — jak nominalny H; żadnych par z danych) → pary z MC (fd_nd, mb_nd) → sieć par / XGB z rozszerzonymi wejściami nóg → fit Asimova
#   (σ_μ: kategorie score pary i kategorie nóg z nowym D+ vs nominalne kategorie nóg). Ta sama paczka co pairnet (charm_work/PAIRNET → Q_lxplus/pairnet),
#   PN_V2=1, wyjście EOS charm_q/pairnet_v2 (batch charm_pairnet2; 12 ext + legdp + train + fit). Wynik → PAIRNET/runs/lxplus_v2/.
source "$HERE/stages/_common.sh"
SRC=$W/PAIRNET; DN1=$QL/dnet1; R=$RQ/pairnet; EOSD=$EOSQ/pairnet_v2; DST=$W/PAIRNET/runs/lxplus_v2
PNV="PN_V2=1 PN_OUT=$EOSD PN_BATCH=charm_pairnet2"
stage_available() { [ -f "$SRC/submit_pn.sh" ] && [ -f "$SRC/pn_legdp.py" ] && [ -f "$SRC/pn_fit.py" ] && [ -f "$DN1/dn_feat.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/pairnet" \
                 && rpush --exclude 'runs/' --exclude '__pycache__' --exclude '*.pyc' --exclude 'pn.dag*' --exclude 'lists/' --exclude '*.bak*' "$SRC/" "$H:$R/" \
                 && rpush --exclude '__pycache__' --exclude '*.pyc' --include 'dn_feat.py' --include 'ext/' --include 'ext/*' --include 'src/' --include 'src/**' --exclude '*' "$DN1/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && $PNV FORCE=${FORCE:-0} ./submit_pn.sh check" | tail -9 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && $PNV FORCE=${FORCE:-0} ./submit_pn.sh submit" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "cd ~/$R; condor_q -batch 2>/dev/null | grep -E 'charm_pairnet2'; echo '## dagman.out'; tail -n 12 pn.dag.dagman.out 2>/dev/null; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/pairnet/*.err ~/$RQ/logs/pairnet/pn_train.*.out ~/$RQ/logs/pairnet/pn_legdp.*.out 2>/dev/null | head -4); do echo \"### \$f\"; tail -n 12 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/fit_v2/DONE"; }
stage_fetch()  { mkdir -p "$DST/fit" && rpull --include '*.md' --include '*.json' --exclude '*' "$H:$EOSD/res/" "$DST/" && rpull --include '*.md' --include '*.json' --exclude '*' "$H:$EOSD/fit_v2/" "$DST/fit/" && rpull --include 'meta.json' --exclude '*' "$H:$EOSD/legdp/" "$DST/"; }
stage_result() { echo "PAIRNET2: $DST/PAIRNET.md, fit/PAIRFIT_K*_F*.md, meta.json (noga D+ z cechami Dalitza)"; }
