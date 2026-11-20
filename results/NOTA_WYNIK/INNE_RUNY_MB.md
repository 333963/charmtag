# Inne runy physics_MinBias 2022 (13,6 TeV) — przegląd 19.11.2026 (tylko odczyt: EOS + rucio list-*)

## Co mamy
EOS `charm_raw_v20/`: dane tylko **data_427929** (459 pl., 39 GB) i **data_428770** (98 pl., 30 GB, strumień physics_Main) + MC (fd_nd, mb_nd/sd/dd/cd).
`charm_raw_v20_d23` — te same dwa runy. Produkcja ntupli v2.0 obejmuje tylko 428770 i 427929 (mirror charmpp, `samples_v2.0.txt`).
GitLab mirror: inne runy MinBias Run 3 nigdzie niewymienione.

## Kontenery rucio `data22_13p6TeV.*.physics_MinBias.merge.AOD*` z niepustą zawartością
| run | tag | pliki | zdarzenia AOD | replika | uwagi |
|---|---|---|---|---|---|
| 427929 | r14956_p5604 (też r15138_p5949 = v2.0) | 450 | 57,9 M | (v2.0 już na EOS) | obecny wynik główny |
| **428700** | r14956_p5604 | 554 | **44,1 M** | tylko taśma INFN-T1 | ten sam okres co 428770; jest też physics_Main i ZeroBias |
| **435333** | r14956_p5604 | 411 | **63,7 M** | tylko taśma INFN-T1 | |
| **428770** | r16546_p6832 / r14956_p5604 | 141 | **16,5 M** | — | strumień **MinBias** tego samego runu (my używamy physics_Main z J12) |
| 428747 | r16487_p5949 | 46 | 4,6 M | taśma INFN-T1 | mały |
| 434876 | r14956_p5604 | 1 | 0,3 M | taśma FZK | pomijalny |
| 435229 | r14956_p5604 (+ part08, r14245, r14470) | 10157 | **1574 M** | taśmy (FZK, IN2P3, PIC, RAL, TRIUMF) | bardzo duży — najpewniej run specjalny o innym μ/menu; wymaga sprawdzenia w runquery/COMA |
Puste (0 plików pod r14956_p5604): 425504, 426646, 426735, 427884, 428580, 428759, 434977, 435075–435180, 435678.

μ (pile-up), menu L1 (czy MBTS bez preskali jak w 427929) i GRL **nie zostały sprawdzone** — tylko rucio, bez AMI/COMA.

## Czy warto
Stat. μ skaluje się ~1/√N (przy podobnym uzysku D na zdarzenie jak 427929, 0,14–0,30× uzysku 428770 physics_Main):
- +428700 +435333 +428770-MinBias ≈ +124 M → 182 M (×3,1) → σ_stat(μ) z ±0,39 do **≈ ±0,22**; z N_D* (±0,21, też stat.) podobnie spada → całość ≈ ±0,35 zamiast ±0,52.
- 435229 (1,57 G zd., ×27) dałby ±0,08 stat., ale tylko jeśli to niski μ i MBTS — bez tego nieużyteczny (przy wysokim pile-up pary z różnych wierzchołków).
Wniosek: **warto**, to jedyna duża dźwignia na wynik główny (dziś zdominowany stat.); najtańsze: 428770-MinBias (ten sam run, ta sama kalibracja, nowy jest tylko strumień) i 428700.

## Co trzeba by uruchomić (nic nie wysłano)
1. Sprawdzić runy w runquery/COMA (μ, menu L1_MBTS_1 PS, GRL, lumi) — strona WWW, tylko odczyt.
2. Staging z taśmy (rucio rule na SCRATCHDISK/LOCALGROUPDISK) — akcja na gridzie, **do zatwierdzenia** (kwota, czas dni).
3. Produkcja CharmAnalysis `charm-prod-v2.0` (config gałęzi minimum-bias) na nowych kontenerach — w ramach produkcji v2.0 albo własne zadania grid; potem rucio → EOS `charm_raw_v20/data_<run>`.
4. h5 (etap jak dla 427929) → ub427nds/ubpnet-podobne etapy ORCH: N_D* i fit par per run, kombinacja μ (wspólny Pythia, stat. niezależne).
Szacunek: tygodnie (staging + grid), lokalnie ~1 dzień po dostarczeniu ntupli.
