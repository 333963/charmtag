#!/bin/bash
# Węzeł mult_charm (próbki: SAMPLES w mu_common.sh; zalecane dodać fd_nd mb_sd mb_dd mb_cd — kind mc). Użycie: mu_job.sh extract <mc|data> <próbka> <tag> <plik.root>   |   mu_job.sh analyze - - all -
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → komendy. Pary danych: NIGDY (bramka G5 w mult_analyze.py).
set -eo pipefail
MU=$(cd "$(dirname "$0")" && pwd); export MU
source "$MU/mu_common.sh"
step=$1; kind=${2:--}; smp=${3:--}; tag=${4:--}; src=${5:--}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/mu_${step}_${tag}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] mult_charm $step $kind $smp $tag host $(hostname) WD=$WD MU_OUT=$MU_OUT"
[ "${DRYRUN:-0}" = 1 ] || mu_env
case $step in
  extract)
    out=$MU_OUT/parts/$smp/${smp}_$tag.npz
    if [ "${FORCE:-0}" != 1 ] && eos_exists "$out"; then echo "pomijam (istnieje) $out"; exit 0; fi
    in=$src; [ "${DRYRUN:-0}" = 1 ] || in=$(fetch_input "$src" "$WD")
    run python3 -u "$MU/mult_extract.py" --kind "$kind" --tag "$tag" ${MAX_EVENTS:+--max-events $MAX_EVENTS} "$in" "$WD/${smp}_$tag.npz"
    [ "${DRYRUN:-0}" = 1 ] || put_output "$WD/${smp}_$tag.npz" "$out" ;;
  analyze)
    # macierz prerejestrowana (PROCEDURA §3): run × schemat × estymator × w2; tylko pojedyncze D na danych + pary MC
    mkdir -p "$WD/res"
    # MC: mieszanka σ/N_gen (PROCEDURA §2), gdy są części mb_sd/dd/cd; ε_rel z fd_nd, gdy są części fd_nd
    MCSPEC="$MU_OUT/parts/mb_nd:2.781e-6"; EPS=()
    for x in "mb_sd:2.143e-6" "mb_dd:1.770e-6" "mb_cd:1.28e-7"; do
      ls "$MU_OUT/parts/${x%%:*}"/*.npz >/dev/null 2>&1 && MCSPEC+=",$MU_OUT/parts/${x%%:*}:${x#*:}"
    done
    ls "$MU_OUT/parts/fd_nd"/*.npz >/dev/null 2>&1 && EPS=(--eps-mc "$MU_OUT/parts/fd_nd")
    echo "MC: $MCSPEC  ε: ${EPS[*]:-mb_nd}"
    G=(); for r in 427929 428770; do ls "$MU_OUT/parts/data_$r"/*.npz >/dev/null 2>&1 && G+=(--data$r "$MU_OUT/parts/data_$r"); done
    run python3 -u "$MU/mult_gates.py" --mc "$MU_OUT/parts/mb_nd" "${G[@]}" --out "$WD/res/gates.json" > "$WD/res/gates.log" 2>&1 || echo "BŁĄD gates"
    for run_ in 427929 428770; do
      d=$MU_OUT/parts/data_$run_; [ -d "$d" ] && ls "$d"/*.npz >/dev/null 2>&1 || { echo "brak części data_$run_ — pomijam"; continue; }
      case $run_ in 427929) schemes="mbts mbts2 all" ;; 428770) schemes="sptrk4 sptrk2 mbts j12 all" ;; esac
      for s in $schemes; do for e in opp blx bl trk500; do for w in none self; do
        run python3 -u "$MU/mult_analyze.py" --data "$d" --mc "$MCSPEC" "${EPS[@]}" --scheme "$s" --est "$e" --w2 "$w" \
            --label "r${run_}_${s}_${e}_w2${w}" --out "$WD/res" > "$WD/res/r${run_}_${s}_${e}_w2${w}.log" 2>&1 || echo "BŁĄD r$run_ $s $e $w"
      done; done; done
    done
    # 19.11 D0-tight (pT(D0) > 6 GeV, prerejestrowane z MC): tylko konfiguracje główne; etykieta _d0t
    for cfg in "427929 mbts" "428770 sptrk4"; do set -- $cfg; d=$MU_OUT/parts/data_$1; ls "$d"/*.npz >/dev/null 2>&1 || continue
      for w in none self; do
        run python3 -u "$MU/mult_analyze.py" --data "$d" --mc "$MCSPEC" "${EPS[@]}" --scheme "$2" --est opp --w2 "$w" --d0-tight \
            --label "r$1_$2_opp_w2${w}_d0t" --out "$WD/res" > "$WD/res/r$1_$2_opp_w2${w}_d0t.log" 2>&1 || echo "BŁĄD d0t r$1 $2 $w"
      done
    done
    if [ "${DRYRUN:-0}" != 1 ]; then
      ( cd "$WD" && tar czf res.tgz res ) && put_output "$WD/res.tgz" "$MU_OUT/res.tgz"
      date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$MU_OUT/DONE"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK mult_charm $step $smp $tag"
