# H5_v20 — ROOT v2.0 → h5 (schemat v1.2 + kolumny zdarzeniowe v2.0), 17.11.2026

Projekt B (charmbkg). Wejście: ROOT v2.0 (CharmAnalysis `charm-prod-v2.0` = 6007843, Oct28),
produkcja `dpt3to150`, `/eos/user/k/klyczak/charm_raw_v20/<próbka>/` (etap `rucio_v20`). Wyjście:
`/eos/user/k/klyczak/charm_h5_v20/parts/<próbka>/<próbka>_<NNNNNN>_<presel>.h5` (jak v1.2: `charm_h5_v12/parts/...`),
`.../charm_h5_v20/meta/{<próbka>,summary}.json`, `.../parts/DONE` (tylko gdy komplet). NIC nie wysłano — wysyła ORCH (`stages/h5_v20.sh`).

Konwerter `charmpp_h5_v20.py` = snapshot `Q_lxplus/h2_pdg/charmpp_h5_v12.py` (repo A / charmtag nietknięte) + dopiski.
Wszystkie 79 kolumn v1.2 liczone identycznie (test: 192 469 dopasowanych wierszy mb_nd _000002, 0 różnic w 76 kolumnach liczbowych).

## Próbki (`samples_v20.txt`)
| klucz | pliki | presele | uwaga |
|---|---|---|---|
| data_428770 | 97 (v1.1: 39) | none, single_tag | physics_Main, zdominowany L1_J12 (V20_FORMAT §5) |
| data_427929 | 458 | none, single_tag | physics_MinBias; LB 94–118 poza GRL (kolumna, bez filtra) |
| mb_nd | 78 (v1.1: 80) | none | |
| fd_nd | 24 | none | eventNumber od 1 w każdym AOD → klucz `event_uid` |
| mb_sd / mb_dd / mb_cd | 10 / 10 / 20 | none | |
fd_sd/dd/cd nie ma w v2.0 na EOS. Numer części `NNNNNN` = numer pliku grid v2.0 (**inny podział niż v1.1** — ten sam numer ≠ te same zdarzenia, poza mb_nd).

## Nowe kolumny (per wiersz kandydata; te same w datasecie `evtable`)
| kolumna | typ | źródło / definicja |
|---|---|---|
| `event_entry` | int64 | wpis drzewa `CharmAnalysis` (= wpis `Tracks`) |
| `event_uid` | uint64 | `(task<<36) \| (nr pliku<<22) \| wpis` (jak charmpp v2.x) — klucz zdarzenia MC/fd |
| `lumiBlock` | uint32 | `EventInfo_lumiBlock` (MC = 1) |
| `passGRL_MinBiasAFP` | int8 | dane: dekoracja GRL; MC = 1 (jak charmpp: brak gałęzi → pass) |
| `nPrimaryVertices`, `nPileupVertices` | int32 | weto pile-up: `nPileupVertices == 0` |
| `nSelBL`, `nTracksPt500Eta25` | int32 | MinBias_BL / pT>0,5 |η|<2,5 |
| `leadTrackPtBL` | float32 | **GeV** (ROOT: MeV; jak charmpp v2.1) |
| `pv_z`, `xi_cal_A`, `xi_cal_C` | float32 | mm; ξ kalorymetryczne (definicja minbias, strony poprawione H8) |
| `afp_nProtons_A/C` | int32 | liczba protonów AFP (`AfpProtons_side` 0 = A, 1 = C) |
| `afp_nInWindow_A/C` | int32 | protony z 0,030 ≤ ξ ≤ 0,100 |
| `afp_xi1_A/C`, `afp_xi2_A/C` | float32 | najmniejsze / drugie ξ na stronie, −1 brak |
| `afp_tag_A/C` | uint8 | **tag minbias z charmpp**: dokładnie 1 proton na stronie i w oknie (charmpp `Structs.h` `tag_reco`, `Preselection::afpPreselection`) |
| `trigPassed_<chain>` | int8 | 11 łańcuchów `mult_common.TRIG_STORE`; MC / brak gałęzi = −1 |
| `trigPS_<chain>` | float32 | prescale HLT×L1 w LB (0 = łańcuch wyłączony); MC = −1 |
| `trigTBP_<L1>` | int8 | `L1_{MBTS_1,MBTS_2,J12,AFP_A_OR_C,RD0_FILLED}` przed prescale; MC = −1 |
`actual_mu`/`average_mu` były już w v1.2 (w MC = 0). Łańcuchy: HLT_j20_L1MBTS_2, HLT_mb_mbts_L1MBTS_1, HLT_mb_sptrk_L1MBTS_1,
HLT_mb_sptrk_pt2_L1MBTS_2, HLT_mb_sptrk_pt4_L1AFP_A_OR_C, HLT_mb_sptrk_pt4_L1MBTS_2, HLT_noalg_L1AFP_A_OR_C, HLT_noalg_L1J12,
HLT_noalg_L1RD0_FILLED, HLT_noalg_mb_L1MBTS_1, HLT_noalg_mb_L1MBTS_2 (wszystkie obecne w 428770 i 427929).

**`evtable`** (tylko `*_none.h5`): 1 wiersz / wpis drzewa (także bez kandydata) — `event_number, run_number, actual_mu, average_mu` + kolumny powyżej;
mianowniki (MULT, w2, tag AFP). `charmpp_h5_v20.py --merge` scala też `evtable` (ten sam offset `event_number += job<<32` w MC).
Atrybuty pliku: `schema = v1.2+v2.0ev`, `n_events_read`, `n_evtable`, `afp_xi_window`, `afp_tag_def`, `d0_window`, `trig_missing`, `source_root`, `job`, `sample`.

## Zmiany zachowania względem h5 v1.2 (świadome)
1. **`single_tag` = `afp_tag_A | afp_tag_C`** (definicja z charmpp, 28.10: dokładnie 1 proton, 0,030–0,100). v1.2: „≥ 1 proton w (0,03; 0,11)”.
   Również `--presel A|C|double_tag`. (ξ w danych v1.1 i tak było złe — H10.)
2. **Kandydaci jak w ROOT v2.0, bez cięć dodatkowych**: D0 1,70–2,05 GeV i L_xy ≥ 0,10 mm (patrz niżej). Okno v1.1 dla samodzielnych D0:
   `--d0-window 1.75,2.00` (nie używane w DAG). Odtworzenie v1.2 bit w bit: `Lxy ≥ 0.25` + D0 `1.75 < m < 2.00`.

## Różnica kandydatów D0/D* v2.0 vs v1.1/v1.2 — przyczyna (sprawdzone, dane NIE poprawiane)
Plik mb_nd _000002 (te same 246 498 zdarzeń w v1.2 i v2.0), klucz `(event_number, decayType, trackId0..2, pdgId córki 0)`:
| | v1.2 | v2.0 | v2.0 tylko |
|---|---|---|---|
| D⁺ wiersze / prawdziwe | 139 078 / 129 | 139 078 / 129 | 0 (identyczne) |
| D0 w 1,75–2,00: wiersze / prawdziwe poprawne (truth 0/1, hypo_swap 0) | 39 725 / 97 | 84 349 / 155 | 44 624 / 58 |
| D* wiersze / prawdziwe | 13 666 / 14 | 28 684 / 23 | 15 018 / 9 |
- **Każdy** kandydat v1.2 jest w v2.0 z identycznymi wartościami wszystkich kolumn; nadmiarowe v2.0 mają **wszystkie** L_xy ∈ [0,10; 0,25) mm.
- Przyczyna: `CharmAnalysis/python/DSelectionConfig.py`, commit 6007843 (charm-prod-v2.0) vs f6df865/cb5ea3d (v1.1, Aug7):
  `MinLxyD0 = 0.25 mm → 0.10 mm` (i `MinLxyD0InDstar = MinLxyD0` → dotyczy też D*), plus okno samodzielnych D0 `1.75–2.00 → 1.70–2.05`
  (D0-w-D* zostaje 1,75–2,00). Kod C++ `DSelectionAlg.cxx` bez zmian; Λc→Λπ wyłączone (nie wpływa na D).
  Opis commita wymienia tylko okno (M9) — obniżenia L_xy w nim brak; `MULT_charm/V20_FORMAT.md` §3 („D* identyczne”, „podłoga 0,10 jak w v1.1”) jest w tym punkcie nieścisły.
- Po cięciu `Lxy ≥ 0.25` i D0 1,75–2,00 v2.0 ≡ v1.2 (D0 39 725, D* 13 666, poprawne 97/14).
- Liczby z notatki zlecenia (D0 155 vs 163, D* 20 vs 26) — 155 = poprawne D0 v2.0 w oknie 1,75–2,00 (zgodne); 163/20/26 nie odtworzyłem
  (v1.2: 97/14, v1.1 bez hypo_swap: 148/14). Kierunek jest odwrotny: v2.0 ma **więcej** kandydatów.
- fd_nd: tego samego pliku nie da się porównać po numerze (inny podział AOD→pliki: _000002 v1.2 99 994 vs v2.0 99 990 zdarzeń;
  eventNumber od 1 w każdym AOD) — porównanie tylko po zdarzeniach mb_nd lub (run, event) w danych.
- Konsekwencja: wydajności / N_D0, N_D* liczone na v1.2 nie przenoszą się 1:1; dla zgodności z nominałem ciąć `Dmeson_Lxy ≥ 0.25`
  (D0, D*) i D0 1,75–2,00 — decyzja analizy, nie konwersji.

## Uruchomienie (lxplus, `~/private/charm_q/Q_lxplus/h5_v20`; normalnie przez ORCH)
`./submit_h5v20.sh check` → wejście (DONE z rucio_v20, liczby plików, tylko `*dpt3to150*`), format v2.0, venv, kolejka, kwota AFS < 85 %, brakujące części.
`./submit_h5v20.sh submit` → `make_dag_h5v20.sh`: `h5v20.dag` = węzeł per (próbka, kawałek ≤ `CHUNK`=25 plików, tylko brakujące części;
`h5v20.sub`, `queue … from lists/chunks/<próbka>_cNNN.txt`) + `FINAL finish` (`finish_h5v20.sh` → `meta_h5v20.py` → meta json, DONE gdy komplet);
`condor_submit_dag -f -batch-name charm_h5v20` w czystym env. `DRYRUN=1` = DAG + komenda 1. procesu. `SAMPLES=mb_nd,fd_nd` zawęża (bez DONE).
`./submit_h5v20.sh status`. Job: `convert_v20.sh <klucz> <root> <NNNNNN> <none[+single_tag]>` (mktemp + trap, xrdcp na EOS, pomija istniejące; `FORCE=1`).
Logi: `Q_lxplus/logs/h5_v20/`. Zasoby: 6 GB RAM, `longlunch`, `STEP=100MB`.

## Test (lxplus, /tmp, skasowane)
- mb_nd _000002: 246 498 zdarzeń → 285 737 wierszy, 128 kolumn, 84 MB, 29 s, RSS 1,9 GB; v1.2 ⊂ v2.0 bit w bit (wyżej).
- fd_nd _000002: 537 492 wierszy, 144 MB, 31 s, RSS 2,3 GB.
- data_428770 _000003: 83 847 zdarzeń → 1,14 M wierszy, 262 MB, 44 s (RSS 4,3 GB przy step 200 MB → w DAG 100 MB); single_tag 12 214 wierszy = kandydaci z `afp_tag`.
  tag A/C 7,8/6,0 % zdarzeń; PS: sptrk_pt4_MBTS_2 25, AFP_A_OR_C 120, mbts_L1MBTS_1 5·10⁴, L1J12 1 (zgodne z V20_FORMAT).
- data_427929 _000003: 45 727 zdarzeń → 55 505 wierszy, 15 MB, 11 s; PS ≈ 1, HLT_mb_mbts_L1MBTS_1 99,9 %.
- `check`, `DRYRUN=1 submit` (31 węzłów, 697 procesów), job MAX_EVENTS=5000, `finish` (meta json), `--merge` z `evtable` — OK w kopii /tmp.
Szacunek pełnego biegu: 697 procesów × 10–45 s (CPU ≈ 4–5 h; przy 60 równolegle ~0,5–1 h + kolejka), wyjście ≈ 43 GB na EOS
(data_428770 ≈ 25, data_427929 ≈ 7, mb_nd ≈ 6,5, fd_nd ≈ 3,5, mb_sd/dd/cd ≈ 1). Kwota EOS: 280 GB / 1 TB (17.11).

## Nieobjęte
TruthProtons (MC, `truthProton_*`/`truth_tag_*` z charmpp: wymaga cięcia E > 0,66 E_beam) — nie dodane. `AfpProtonsAOD_*`, `AfpTracks_*` — nie.
Kandydaci D z `dpt2to3` — nie (osobna produkcja, nie mieszać).
