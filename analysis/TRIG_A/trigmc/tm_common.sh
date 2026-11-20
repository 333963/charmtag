# tm_common.sh — paczka TRIG_A/trigmc (projekt A: zamknięcie wag w2 wobec emulacji L1_J12 na generatorze). source; wymaga ../common.sh.
# Na lxplus: ~/private/charm_q/Q_lxplus/trigmc (= lokalnie charm_work/TRIG_A/trigmc, kopiowane przez etap ORCH trigmc).
TM=${TM:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$TM/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"
: "${TM_OUT:=$EOS_OUT/charm_q/trigmc}"      # gen/gen_<seed>.tsv.gz(+.json,.h3em,.h3all), res/{trigmc.json,TRIGMC.md}, DONE
: "${TM_NJ:=200}"; : "${TM_NEV:=1000000}"    # 200 × 1 M zdarzeń ND
: "${TM_SEED0:=70000}"                       # ziarna 70001… (rozłączne z A3: 1–500 i tune)
export TM Q TM_OUT TM_NJ TM_NEV TM_SEED0
