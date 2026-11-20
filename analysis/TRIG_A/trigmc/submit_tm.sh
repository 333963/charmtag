#!/bin/bash
# submit_tm.sh check|submit|status   (lxplus, ~/private/charm_q/Q_lxplus/trigmc) — projekt A, tylko generator (bez danych).
#  check  — condor_q (brak charm_trigmc), kwota AFS < 85 %, brak DONE (FORCE=1 nadpisuje)
#  submit — check → build pythia_dd_jet w podpowłoce LCG + test 300 zdarzeń → tm.dag → condor_submit_dag -f w czystym env
set -uo pipefail
TM=$(cd "$(dirname "$0")" && pwd); export TM
source "$TM/tm_common.sh"
cd "$TM"; mkdir -p "$Q/logs/trigmc"
status() {
  condor_q -batch 2>/dev/null | grep -E 'charm_trigmc|Total for klyczak' || true
  echo "gen: $(ls "$TM_OUT/gen" 2>/dev/null | grep -c 'tsv\.gz\.json$') / $TM_NJ plików; res: $(ls "$TM_OUT/res" 2>/dev/null | tr '\n' ' '); DONE: $(cat "$TM_OUT/DONE" 2>/dev/null || echo nie)"
  [ -f tm.dag.dagman.out ] && tail -3 tm.dag.dagman.out || true
}
check() {
  local ok=0 q
  condor_q -batch 2>/dev/null | grep -qE 'charm_trigmc' && { echo "STOP: trigmc już w kolejce"; ok=1; }
  [ -s "$TM_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$TM_OUT/DONE")); FORCE=1 nadpisuje"; ok=1; }
  for f in pythia_dd_jet.cc tm_combine.py tm_job.sh tm_gen.sub tm.sub; do [ -f "$TM/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  TM_OUT=$TM_OUT  NJ=$TM_NJ NEV=$TM_NEV seed0=$TM_SEED0  binarka: $([ -x bin/pythia_dd_jet ] && echo jest || echo 'do zbudowania')"
  return $ok
}
build() {   # podpowłoka: LCG nie może zostać w środowisku submitu
  ( set -e; q_env; mkdir -p "$TM/bin"
    g++ -O2 -std=c++17 $(pythia8-config --cxxflags) "$TM/pythia_dd_jet.cc" -o "$TM/bin/pythia_dd_jet" $(pythia8-config --libs) -ldl
    T=$(mktemp -d); cd "$T"
    "$TM/bin/pythia_dd_jet" 300 1 t.tsv A3 13600 nd > t.log 2>&1 || { tail -20 t.log; exit 1; }
    ne=$(grep -c '^E' t.tsv || true); nh=$(wc -l < t.tsv.h3em)
    awk -F'\t' '$1=="E" && NF!=11 {bad=1} $1=="D" && NF!=9 {bad=1} END{exit bad}' t.tsv || { echo "zły format TSV"; exit 1; }
    [ "$nh" -gt 10 ] || { echo "pusty histogram h3em"; exit 1; }
    gzip t.tsv; mkdir g; mv t.tsv.gz g/gen_1.tsv.gz; mv t.tsv.json g/gen_1.tsv.gz.json; mv t.tsv.h3em g/gen_1.tsv.h3em; mv t.tsv.h3all g/gen_1.tsv.h3all
    echo "pythia_dd_jet test OK: $ne zdarzeń z D / 300, komórek h3em: $nh"; cd /; rm -rf "$T" )
}
case ${1:-check} in
  status) status ;;
  check)  check && echo "OK — gotowe (./submit_tm.sh submit)" ;;
  submit) check || exit 1
          build || { echo "STOP: build/test pythia_dd_jet nieudany"; exit 1; }
          chmod +x "$TM"/*.sh
          { echo "# wygenerowane $(date -u +%FT%TZ) NJ=$TM_NJ NEV=$TM_NEV seed0=$TM_SEED0"
            echo "JOB gen $TM/tm_gen.sub"; echo "VARS gen nj=\"$TM_NJ\" nev=\"$TM_NEV\" seed0=\"$TM_SEED0\" tmout=\"$TM_OUT\""
            echo "SCRIPT POST gen /bin/true"          # pojedyncze nieudane joby generacji nie blokują combine
            echo "JOB combine $TM/tm.sub"; echo "VARS combine node=\"combine\" args=\"combine\" tmout=\"$TM_OUT\""; echo "RETRY combine 1"
            echo "PARENT gen CHILD combine"; } > tm.dag
          env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_trigmc tm.dag
          sleep 3; status ;;
  calib)  # trigcal: ponowny combine z siatką drobną i kalibracją na j12eff.json (bez generacji)
          [ -s "$TM/j12eff.json" ] || { echo "STOP: brak j12eff.json"; exit 1; }
          [ "$(ls "$TM_OUT/gen" 2>/dev/null | grep -c 'tsv\.gz\.json$')" -ge 150 ] || { echo "STOP: za mało plików gen"; exit 1; }
          condor_q -batch 2>/dev/null | grep -qE 'charm_trigcal' && { echo "STOP: trigcal w kolejce"; exit 1; }
          [ -s "$TM_OUT/res_cal/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone"; exit 1; }
          chmod +x "$TM"/*.sh
          env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit tm.sub node=calib args=calib tmout="$TM_OUT" -batch-name charm_trigcal -append '+JobFlavour="workday"' ;;
  *) sed -n 2,4p "$0" ;;
esac
