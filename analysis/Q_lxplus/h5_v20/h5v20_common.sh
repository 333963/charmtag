# h5v20_common.sh — paczka h5_v20 (17.11.2026; source, nie uruchamiać). Wymaga ../common.sh (Q: EOS_OUT, q_env, fetch_input, put_output, eos_exists, job_tag).
# common.sh NIE jest zmieniany — ścieżki v2.0 nadpisywane tylko tutaj.
HV=${HV:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$HV/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"
: "${PTPROD:=dpt3to150}"                              # NIE mieszać z dpt2to3 (te same zdarzenia, inne dzielenie plików)
: "${RAW20:=$EOS_OUT/charm_raw_v20}"                  # wejście (rucio_v20), tylko odczyt
: "${H5V20:=$EOS_OUT/charm_h5_v20}"                   # wyjście: $H5V20/parts/<próbka>/<próbka>_<NNNNNN>_<presel>.h5, $H5V20/meta/*.json, parts/DONE
: "${CHUNK:=25}"                                      # plików ROOT na węzeł DAG
: "${STEP:=100MB}"                                    # porcja uproot (bez spacji: trafia do environment w DAG) (dane 428770: RSS ~4 GB przy 200 MB)
: "${SAMPLES:=}"                                      # pusty = wszystkie z samples_v20.txt; np. SAMPLES=mb_nd,data_428770 (przecinki)
: "${MAXJOBS:=60}"
export HV Q PTPROD RAW20 H5V20 CHUNK STEP SAMPLES MAXJOBS
PARTS=$H5V20/parts
hv_samples() { awk '!/^#/ && NF>=4 {print $1, $2, $3, $4}' "$HV/samples_v20.txt" | while read -r k d p n; do
  if [ -z "$SAMPLES" ] || [[ ",$SAMPLES," == *",$k,"* ]]; then echo "$k $d $p $n"; fi; done; }
hv_lists() {  # lists/<klucz>.txt (pliki ROOT) + lists/jobs.txt (klucz plik tag presele) + lists/missing.txt (brakujące części)
  mkdir -p "$HV/lists"; : > "$HV/lists/jobs.txt"; : > "$HV/lists/missing.txt"; local rc=0
  while read -r k d p n; do
    if [ ! -d "$RAW20/$d" ]; then echo "  $k: BRAK $RAW20/$d"; rc=1; continue; fi
    find "$RAW20/$d" -name "*${PTPROD}*.root" | sort > "$HV/lists/$k.txt"
    local m; m=$(wc -l < "$HV/lists/$k.txt"); local dn="brak DONE"; [ -f "$RAW20/$d/DONE" ] && dn="DONE"
    [ "$m" -eq "$n" ] || { echo "  $k: UWAGA $m plików (oczekiwane $n)"; }
    [ -f "$RAW20/$d/DONE" ] || rc=1
    local miss=0
    while read -r f; do
      local t; t=$(job_tag "$f"); [ -n "$t" ] || { echo "  STOP: brak numeru pliku w $f"; rc=1; continue; }
      echo "$k $f $t $p" >> "$HV/lists/jobs.txt"
      local need=0 x
      for x in ${p//,/ }; do [ -s "$PARTS/$k/${k}_${t}_${x}.h5" ] || need=1; done
      [ $need = 1 ] && { echo "$k $f $t $p" >> "$HV/lists/missing.txt"; miss=$((miss+1)); }
    done < "$HV/lists/$k.txt"
    printf "  %-12s ROOT %4d (oczek. %4d, %s)  brakujących części: %d\n" "$k" "$m" "$n" "$dn" "$miss"
  done < <(hv_samples)
  return $rc
}
