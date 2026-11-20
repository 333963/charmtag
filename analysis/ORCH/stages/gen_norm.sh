# GEN_NORM (16.11.2026, etap GEN_mcstat; 17.11 → ROOT v2.0 dpt3to150: mb_nd 80 + 428770 97, okno D0 1,75–2,00, B1 v2.0 = B1a/B1b/B1c w PROCEDURA): odczyt ROOT mb_nd + 428770 + nogi direct z h5/scores nominalnych →
#   (a) mapa wag w2d.npy (+ test akceptacji c_w = mb_rates_2d.json), (b) mapy ε selekcji direct z mb_nd (z flagą okna D*) + tabele prawda↔nogi,
#   (c) n_ch→n_trk tylko gdy TruthParticles (v1.1: brak → pominięte). Paczka Q_lxplus/gen_norm (src = snapshot charmtag 2.1.0).
#   Wynik → EOS charm_q/gen_norm → charm_work/GEN_mcstat/runs/lxplus/. Dalej lokalnie: GEN_mcstat/PROCEDURA.md (prerejestracja).
#   BLOKADA: etap startuje tylko gdy JEDNOCZEŚNIE ORCH/GEN_NORM_OK istnieje i state/gen_norm ≠ hold.
#   Odblokowanie: touch ORCH/GEN_NORM_OK && echo "pending $(date +%s) 0" > ORCH/state/gen_norm   (bez flagi → „skipped”, nic nie wysyła)
source "$HERE/stages/_common.sh"
SRC=$QL/gen_norm; R=$RQ/gen_norm; EOSD=$EOSQ/gen_norm_v20; DST=$W/GEN_mcstat/runs/lxplus_v20   # 17.11: ROOT v2.0 (GN_FMT=v20, GN_OUT=charm_q/gen_norm_v20); v1.1: .bak_v11
# 18.11 (zatwierdzone działanie na danych + §3a/§3b PROCEDURA): wariant 78 plików mb_nd — GN_MB_EXCLUDE=000037,000070 (pliki skasowane z gridu)
# 19.11: RECOMBINE=1 — ponowny combine po poprawce przetagowania nóg po eventNumber (tagi h5 000003–000010 ≠ ROOT v2.0); części ROOT na EOS zostają
stage_available() { [ -f "$SRC/submit_gn.sh" ] && [ -f "$SRC/src/charmtag/pairs/build.py" ] && [ -f "$HERE/GEN_NORM_OK" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/gen_norm" && rpush --delete --exclude '__pycache__' --exclude '*.pyc' --exclude 'gn.dag*' --exclude 'lists/' --exclude '.pytest_cache' --exclude '*.bak*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && GN_MB_EXCLUDE=000037,000070 RECOMBINE=1 ./submit_gn.sh check" | tail -10 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && GN_MB_EXCLUDE=000037,000070 RECOMBINE=1 ./submit_gn.sh submit; GN_MB_EXCLUDE=000037,000070 ./submit_gn.sh status | tail -4" | tail -12 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { rssh "condor_q -batch 2>/dev/null | grep -E 'charm_gn'; condor_q -totals 2>/dev/null | tail -1" | sed 's/^/[q] /' | tee -a "$LOG" >/dev/null
  mkdir -p "$DST"; rssh "cd ~/$R; echo '## dagman.out'; tail -n 30 gn.dag.dagman.out 2>/dev/null; echo '## gn.dag.rescue'; ls gn.dag.rescue* 2>/dev/null; echo '## logs/gen_norm'; ls -lt ~/$RQ/logs/gen_norm | head -12; for f in \$(ls -t ~/$RQ/logs/gen_norm/*.err 2>/dev/null | head -3); do echo \"### \$f\"; tail -n 25 \$f; done; echo '## EOS'; ls -la $EOSD $EOSD/legs 2>/dev/null | head -40; echo \"parts mc \$(ls $EOSD/parts/mc 2>/dev/null | wc -l) data \$(ls $EOSD/parts/data 2>/dev/null | wc -l)\"" > "$DST/diag.txt" 2>&1
  rssh "test -f $EOSD/DONE && test -f ~/$R/gn.dag && test $EOSD/DONE -nt ~/$R/gn.dag"; }   # stare DONE (sprzed gn.dag) ignorowane
stage_fetch()  { mkdir -p "$DST" && rpull --exclude 'parts/' --include '*.parquet' --include '*.npy' --include '*.npz' --include 'DONE' "${mdjson[@]}" "$H:$EOSD/" "$DST/"; }
stage_result() { local v; v=$(grep -o '"verdict": "[A-Z]*"' "$DST/acceptance.json" 2>/dev/null | cut -d'"' -f4)
  echo "GEN_NORM: wyniki → $DST (summary.md; akceptacja w2: ${v:-?}); dalej GEN_mcstat/PROCEDURA.md — μ na danych dopiero po decyzji"; }
