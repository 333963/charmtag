#!/bin/bash
# submit_upn427.sh check | submit | status — DAG UBPNET427: ext × ⌈458/UPN_NG⌉ → fit. Submit w czystym env (bez LCG), condor_submit_dag -f.
# Węzły ext z istniejącą częścią na EOS pomijane; FORCE=1 → stary wynik odkładany do old_<data> (nie kasowany).
set -uo pipefail
UPN=$(cd "$(dirname "$0")" && pwd); export UPN
source "$UPN/upn427_common.sh"
mkdir -p "$Q/logs/ubpnet427"
cmd=${1:-check}
check() {
  local ok=0 n f q T
  uj_gate || ok=1
  n=$(upn_h5 | wc -l); echo "h5 427929 (filtr UB427, klasa MB): $n plików → $(upn_ngroups) grup po $UPN_NG"; [ "$n" -eq 458 ] || { echo "STOP: h5 427929 ≠ 458"; ok=1; }
  for f in pn_extract.py pn_train.py pn_fit.py pn_hybsf.py dn_feat.py pn_common.sh; do [ -f "$UPN_PN/$f" ] || { echo "STOP: brak $UPN_PN/$f"; ok=1; }; done
  grep -q "keep &= sb_mask(T\[\"m\"\], isp)" "$UPN_PN/pn_extract.py" || { echo "STOP: pn_extract.py bez maski sb_mask w znanym miejscu"; ok=1; }
  for f in "$SB3/train/pn_models_sbtrain.pt" "$SB3/train/xgbd_fold4.json" "$HY3/res/hybsf.json" "$HY3/ds/hybds.json" "$U427_PAIRS/pair_inputs.json" \
           "$U427_PAIRS/pair_inputs.npz" "$U427_PAIRS/pairs_full.npz" "$U427_PAIRS/fit2d_unblind.json" "$PN_MODELS/D0_H/meta.json" "$PN_MODELS/Dplus_H/meta.json"; do
    [ -s "$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  grep -q '"data_427929|DpD0|D+|pnetd"' "$HY3/res/hybsf.json" || { echo "STOP: hybsf.json bez kluczy data_427929"; ok=1; }
  grep -q '"data_427929"' "$HY3/ds/hybds.json" || { echo "STOP: hybds.json bez wierszy data_427929"; ok=1; }
  [ "$(ls "$P3"/parts/pn_*_nd_*.npz 2>/dev/null | wc -l)" -ge 2 ] || { echo "STOP: brak części MC v3 ($P3/parts)"; ok=1; }
  [ "$(ls "$SB3"/parts/pn_data_428770_*.npz 2>/dev/null | wc -l)" -ge 1 ] || { echo "STOP: brak części pasm 428770 ($SB3/parts) — test transferu"; ok=1; }
  T=$(mktemp -d /tmp/upn4_self_XXXXXX)
  ( uj_env >/dev/null && cd "$UPN" && CHARM_UNBLIND=1 UPN_PN=$UPN_PN PYTHONPATH="$UPN:$UPN_PN:$Q/ub_jk:${PYTHONPATH:-}" TMPDIR=$T \
      nice python3 upn_fit.py --selftest --no-torch --schemes base_sum --nboot 4 --threads 4 --sample data_427929 --npair-scale "$UPN_FCH" \
      --scale-now 1.0 --ref-sb self --out "$T/res" 2>&1 | tail -2 && grep -q "SELFTEST: OK" "$T/res/UBPNET.md" ) \
    || { echo "STOP: selftest upn_fit (tryb 427929)"; ok=1; }
  rm -rf "${T:?}"
  condor_q -batch 2>/dev/null | grep -q charm_ubpnet427 && { echo "STOP: charm_ubpnet427 w kolejce"; ok=1; }
  [ -s "$UPN_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($UPN_OUT/DONE); FORCE=1 → nowy fit"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "AFS: ${q:-?}%  części danych już na EOS: $(ls "$UPN_OUT"/parts/pnf_data_427929_*.npz 2>/dev/null | wc -l)"
  [ $ok = 0 ] && echo "OK — można wysłać"; return $ok
}
case $cmd in
  check) check ;;
  submit)
    check || exit 1
    if [ -e "$UPN_OUT/DONE" ]; then
      old="$UPN_OUT/old_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$old"
      for f in DONE UBPNET427.md ubpnet.json ubpnet.log; do [ -e "$UPN_OUT/$f" ] && mv "$UPN_OUT/$f" "$old/"; done; echo "stary wynik → $old"
    fi
    P="penv=\"UPN427_OUT=$UPN_OUT UPN427_NG=$UPN_NG UPN_K=$UPN_K UPN_NBOOT=$UPN_NBOOT UPN_CPUS=$UPN_CPUS UPN_FCH=$UPN_FCH OMP_NUM_THREADS=4\""
    NG=$(upn_ngroups); par=""
    {
      for ((g = 0; g < NG; g++)); do
        [ -s "$UPN_OUT/parts/pnf_data_427929_$g.npz" ] && continue
        echo "JOB ext_$g $UPN/upn427.sub"; echo "VARS ext_$g step=\"ext\" g=\"$g\" cpus=\"4\" mem=\"12GB\" disk=\"6GB\" flavour=\"workday\" $P"; echo "RETRY ext_$g 2"
        par+="ext_$g "
      done
      echo "JOB fit $UPN/upn427.sub"; echo "VARS fit step=\"fit\" g=\"-\" cpus=\"$UPN_CPUS\" mem=\"32GB\" disk=\"16GB\" flavour=\"tomorrow\" $P"; echo "RETRY fit 1"
      [ -n "$par" ] && echo "PARENT ${par}CHILD fit"
    } > "$UPN/upn427.dag"
    echo "upn427.dag: ext $(grep -c '^JOB ext_' "$UPN/upn427.dag") / $NG, fit 1"
    cd "$UPN" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_ubpnet427 upn427.dag ;;
  status) condor_q -batch 2>/dev/null | grep charm_ubpnet427; ls -la "$UPN_OUT" 2>/dev/null | head; ls "$UPN_OUT/parts" 2>/dev/null | grep -c npz ;;
  *) echo "użycie: $0 check|submit|status"; exit 2 ;;
esac
