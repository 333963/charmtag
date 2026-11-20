# UB_TOPO (19.11.2026, projekt A; pary w danych zatwierdzone 18.11): kształt |Δη| sygnału par D–D (428770, odsłonięta próbka 05.11) vs Pythia/DPS
#   z generatora (Q_lxplus/ub_topo/ref/topo_gen.json). Mechanika pd_diff (pairs/fit2d --unblind --diff dy), raport topo_report.py (skan λ).
#   Prognoza mocy (ref/forecast.json, PROGNOZA.md): czysty DPS vs Pythia ~1,2σ; klasy triggera i 427929 — brak mocy (nie liczone).
#   Paczka Q_lxplus/ub_topo (używa ../pd_diff/{pd_common.sh,tools/diff_fill.py} i ../pf_pairs/src na lxplus — tylko odczyt). Wynik → PF2_lxplus/runs/lxplus/unblind/ubtopo.
source "$HERE/stages/_common.sh"
SRC=$QL/ub_topo; R=$RQ/ub_topo; EOSD=$EOSQ/ub_topo/dy; DST=$W/PF2_lxplus/runs/lxplus/unblind/ubtopo
stage_available() { [ -f "$SRC/submit_topo.sh" ] && [ -f "$SRC/topo_report.py" ] && [ -s "$SRC/ref/topo_gen.json" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/ubtopo" && rpush --exclude '__pycache__' --exclude '*.pyc' --exclude '*.dag*' --exclude 'cfg/' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && ./submit_topo.sh check" | tail -8 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && ./submit_topo.sh submit" | tail -6 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_ubtopo'; echo '## logs'; for f in \$(ls -t ~/$RQ/logs/ubtopo/*.err ~/$RQ/logs/ubtopo/*.out 2>/dev/null | head -4); do echo \"### \$f\"; tail -n 15 \$f; done" > "$DST/diag.txt" 2>&1
                 head -1 "$DST/diag.txt" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
                 rtest "$EOSD/DONE"; }
stage_fetch()  { mkdir -p "$DST" && rpull --include 'UBTOPO.md' --include 'ubtopo.json' --include 'tfin_diff.md' --include 'tfin_diff.json' --include 'DONE' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "UB_TOPO: $DST/UBTOPO.md (kształt |Δη| par vs Pythia/DPS, 428770)"; }
