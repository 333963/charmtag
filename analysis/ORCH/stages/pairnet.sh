# PAIRNET (17.11.2026, projekt A): sieć stricte pod pary D–D (double charm). TYLKO MC (fd_nd, mb_nd, h5 v2.0) — żadnych par z danych.
#   Nogi po preselekcji nominalnego H ≥ WP90 → pary w zdarzeniu → sieć par (Deep Sets na 18 cechach H obu nóg; wariant z kinematyką pary),
#   kontrola XGB, odniesienie = złożenie score nóg H i punkt nominalny (obie nogi ≥ WP65). Metryki: tło par przy tej samej wydajności sygnału,
#   rzeźbienie, zależność ε_S od Δφ/ΔR/LS-OS. Paczka charm_work/PAIRNET (+ dn_feat.py, ext/, src/ z Q_lxplus/dnet1) → Q_lxplus/pairnet
#   (batch charm_pairnet; 12 ext + train). Wynik: EOS charm_q/pairnet/res → PAIRNET/runs/lxplus/{PAIRNET.md,pairnet.json}.
source "$HERE/stages/_common.sh"
SRC=$W/PAIRNET; DN1=$QL/dnet1; R=$RQ/pairnet; EOSD=$EOSQ/pairnet; DST=$W/PAIRNET/runs/lxplus
stage_available() { [ -f "$SRC/submit_pn.sh" ] && [ -f "$SRC/pn_extract.py" ] && [ -f "$SRC/pn_train.py" ] && [ -f "$DN1/dn_feat.py" ] && [ -d "$DN1/src/charmbkg" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/pairnet" \
                 && rpush --exclude 'runs/' --exclude '__pycache__' --exclude '*.pyc' --exclude 'pn.dag*' --exclude 'lists/' --exclude '*.bak*' "$SRC/" "$H:$R/" \
                 && rpush --exclude '__pycache__' --exclude '*.pyc' --include 'dn_feat.py' --include 'ext/' --include 'ext/*' --include 'src/' --include 'src/**' --exclude '*' "$DN1/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && FORCE=${FORCE:-0} ./submit_pn.sh check" | tail -9 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && FORCE=${FORCE:-0} ./submit_pn.sh submit; ./submit_pn.sh status | tail -4" | tail -10 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "cd ~/$R; condor_q -batch 2>/dev/null | grep -E 'charm_pairnet'; echo '## dagman.out'; tail -n 12 pn.dag.dagman.out 2>/dev/null; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/pairnet/*.err ~/$RQ/logs/pairnet/pn_train.*.out 2>/dev/null | head -3); do echo \"### \$f\"; tail -n 15 \$f; done; echo parts \$(ls $EOSD/parts 2>/dev/null | wc -l)" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include '*.md' --include '*.json' --exclude '*' "$H:$EOSD/res/" "$DST/"; }
stage_result() { echo "PAIRNET: $DST/PAIRNET.md + pairnet.json (sieć par D–D na MC; score i modele na EOS charm_q/pairnet/res)"; }
