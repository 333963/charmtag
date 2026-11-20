# PAIRNET3 (18.11.2026, projekt A): łańcuch sieci par z arbitrażem hipotez D0 (jedna hipoteza Kπ/πK na parę torów, jak bieg nominalny):
#   ekstrakcja MC + nogi z pasm bocznych danych → trening z tłem z par SB×SB 428770 (zatwierdzone 18.11) → pomiar wydajności metodą pary hybrydowej
#   (pojedyncze D). Okno sygnału par w danych nieczytane. Paczka charm_work/PAIRNET → Q_lxplus/pairnet; wyniki → PAIRNET/runs/v3/{sbtrain,hybsf}.
source "$HERE/stages/_common.sh"
SRC=$W/PAIRNET; R=$RQ/pairnet; EOSD=$EOSQ/pairnet_hyb3; EOST=$EOSQ/pairnet_sb3/train; DST=$W/PAIRNET/runs/v3
stage_available() { [ -f "$SRC/submit_pn3.sh" ] && [ -f "$SRC/pn_hybsf.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/pairnet" && rpush --exclude 'runs/' --exclude '__pycache__' --exclude '*.pyc' --exclude '*.dag*' --exclude 'lists/' --exclude '*.bak*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && FORCE=${FORCE:-0} ./submit_pn3.sh check" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && FORCE=${FORCE:-0} ./submit_pn3.sh submit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_pn3'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/pairnet/pn_sbtrain.*.out ~/$RQ/logs/pairnet/pn_hybsf.*.out ~/$RQ/logs/pairnet/pn_hybsf.*.err 2>/dev/null | head -3); do echo \"### \$f\"; tail -n 12 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST/sbtrain" "$DST/hybsf" && rpull --include '*.md' --include 'sbtrain.json' --exclude '*' "$H:$EOST/" "$DST/sbtrain/" && rpull --include '*.md' --include 'hybsf.json' --exclude '*' "$H:$EOSD/res/" "$DST/hybsf/"; }
stage_result() { echo "PAIRNET3: $DST/sbtrain/SBTRAIN.md, $DST/hybsf/HYBSF.md (sieć par z arbitrażem hipotez D0)"; }
