# Ntuple ROOT v2.0 (CharmAnalysis `charm-prod-v2.0`, Oct28) — format i co zmienić w `gen_norm` (16.11.2026)

Źródła (tylko odczyt, lokalny mirror `~/Projects/charm/gitlab_mirror`): CharmAnalysis tag `charm-prod-v2.0` = commit `6007843`
(gałąź minimum-bias, poprzednik `cb5ea3d` = stan v1.1): `CharmAnalysis/python/{Branches,OutputConfig,TrackDedxConfig}.py`,
`data/charm_minbias_dpt{3to150,2to3}.yaml`, `data/cp_config_minbias_base.yaml`, `src/{EventInfoAlg,DecorateDMesons}.cxx`;
charmpp, gałąź minimum-bias: `docs/CHARMPP_v2.0_PRODUCTION.md`, `docs/CHARMPP_MINBIAS.md`, `scripts/minbias_v2/samples_v2.0.txt`.
Porównanie z rzeczywistym plikiem v1.1 (`~/Documents/charm/mc23_minbiasND_Dpt3to150GeV_v9/*.root`, 126 gałęzi) — odczyt uproot.
**Plików v2.0 jeszcze nie widziałem** (`/eos/user/k/klyczak/charm_raw_v20/` 16.11 22:50: nie istnieje) → pierwszy krok na lxplus:
`submit_mu.sh check` (sprawdza nSelBL/leadTrackPtBL w 1. pliku każdej próbki).

## 1. Pliki i próbki
- Kontenery grid: `user.malewick.<data.minbias.428770|mc23.minbiasND|…>.CharmAnalysis.Dpt{3to150,2to3}GeV.v2.0.Oct28_EXT0`.
  Nazwa pliku: `user.malewick.<task>.EXT0._NNNNNN.CharmAnalysis.{data.2022|mc.MC23a}.charm_minbias_dpt{3to150|2to3}.root`.
- **Dwie produkcje pT D (nowość):** `dpt3to150` i `dpt2to3` — osobne zadania grid, inne dzielenie plików → te same zdarzenia w obu,
  ale numery plików/wpisów różne. **Nie mieszać** (zdarzenia liczone 2×). Wszystko poniżej: `dpt3to150` (jak v1.1), zmienna `PTPROD`.
- Liczba plików (3–150): data 428770 **97** (v1.1: 39), mb_nd 80 (=), fd_nd 24, data 427929 ~583 (zadanie mogło nie być skończone 1.11).
- Wpisy drzewa: data 428770 7 989 496 (z 11 503 705 w AOD; filtr nie zależy od D — ta sama liczba w d23 i d150), mb_nd **19 714 175 = v1.1**
  (te same zdarzenia, test akceptacji gen_norm liczy 19 714 175 ✓).
- Drzewa: `CharmAnalysis` + `Tracks` (wpis i ↔ wpis i; `storeAllTracksTree: True`, widok `SavedTracks` pT > 200 MeV — jak v1.1) +
  histogramy `metadata`, `DSelection_Vertices` (oczekiwane jak v1.1 — do potwierdzenia `check`).

## 2. Gałęzie zdarzeniowe (sufiks `_NOSYS` przy wszystkich `%SYS%`, jak w v1.1)
| gałąź | v1.1 | v2.0 | uwagi |
|---|---|---|---|
| `EventInfo_eventNumber`, `_runNumber`, `_lumiBlock`, `_bcid` | ✓ | ✓ | MC: runNumber 420000/RandomRunNumber 0; eventNumber w fd powtarza się (klucz = plik×wpis) |
| `EventInfo_nPrimaryVertices_NOSYS`, `_nPileupVertices_NOSYS`, `_pv_{x,y,z}_NOSYS` | ✓ | ✓ | npv = 1 zawsze; weto: `nPileupVertices == 0` (dane 5,9 % >0, MC 0,65 %) |
| `EventInfo_nTracksPt500Eta25_NOSYS` | ✓ | ✓ | ślady pT>0,5, |η|<2,5 bez wymogów jakości/IP |
| **`EventInfo_nSelBL_NOSYS`** (int) | – | ✓ | n_sel^BL analizy minbias: pT>0,5, |η|<2,5, ≥1 pix, ≥6 SCT, hit wewn., |d0_BL|<1,5 mm, χ²-prob>0,01 dla pT>10 |
| **`EventInfo_leadTrackPtBL_NOSYS`** (float) | – | ✓ | wiodący pT śladów BL, **MeV** w ROOT (charmpp v2.1: GeV) |
| `EventInfo_xi_cal_{A,C}_NOSYS` (+ `_t15`, `_t15_E200`, `_t15_E500`) | `xi_cal_{plus,minus}[_E200/_E500]` | zmienione | strony v1.1 były zamienione (H8); w MC ≈ 0 |
| `EventInfo_passGRL_MinBiasAFP` (char) | – | dane | dekoracja, bez filtrowania (428770: wszystko =1; 427929: LB 94–118 poza GRL) |
| `EventInfo_trigPassed_<chain>_NOSYS` | 28 łańcuchów | **106** (dane) | lista z `cp_config_minbias_base.yaml`; m.in. `HLT_mb_mbts_L1MBTS_1`, `HLT_mb_sptrk_pt{2,4,6,8}_L1MBTS_2`, `HLT_noalg_L1J12`, `HLT_noalg_L1AFP_A_OR_C` |
| **`EventInfo_trigPS_<chain>`** (float) | – | 106 (dane) | całkowity prescale HLT×L1 w LB zdarzenia (np. AFP_A_OR_C 120/60, sptrk_pt4_MBTS_2 25/50, mbts_L1MBTS_1 ≈ 5·10⁴ w 428770) |
| **`EventInfo_trigTBP_<L1>`** (char) | – | 31 (dane) | L1 przed prescale |
| `AfpProtons_*_NOSYS` | ✓ | ✓ (optyka 2022) | ξ w danych poprawne dopiero w v2.0 (H10) |
| `AfpTracks_*`, `AfpProtonsAOD_*` (dane), `TruthProtons_*` (MC) | – | ✓ | nowe |
| `EventInfo_mcEventWeights`, `_weight_mc_NOSYS`, `_mcChannelNumber`, `_pdfInfo_*` | ✓ | ✓ | wagi = 1 |

## 3. Kandydaci D (`DMesons_*_NOSYS`) i prawda
- Zestaw gałęzi **identyczny z v1.1** (te same nazwy z `_NOSYS`: `decayType`, `pdgId`, `pt`, `eta`, `phi`, `m`, `DeltaMass`, `D0Index`,
  `mKpi1/2`, `mPhi1/2`, `fitOutput__{Lxy,LxyErr,Chi2,Impact,ImpactSignificance,…}`, `cosThetaP`, `costhetastar`, `ptcone40…100`,
  `daughterInfo__{pdgId,charge,pt,eta,phi,d0,trackId,z0SinTheta,…,dEdx,npts_dEdx}`, MC: `TruthMatchCategory`, `truthPdgId`,
  `daughterInfo__truthPdgId`, `daughter_truthMatchProb`, `truthAncestor*`). Jednostki MeV/mm.
- **Zmiany zawartości:** okno m(Kπ) D0 **1,70–2,05 GeV** (v1.1: 1,75–2,00) → D0 ×1,30–1,40; D⁺ i D* **identyczne** z v1.2
  (charmpp `compare_v12.py`, każda próbka i bin pT); Λc→Λπ wyłączony (`doLambdaCLambdaPi: False`); podłoga L_xy D0/D* 0,10 mm jak w v1.1.
- `daughterInfo__trackId` = indeks w kontenerze `OutTracks` (pT > 100 MeV), **nie** w drzewie `Tracks` (pT > 200 MeV) → nie da się
  wprost wyciąć córek z `Tracks`; MULT odejmuje córki po (pT, η) z `daughterInfo`.
- `TruthDmesons_*` (MC): te same 22 gałęzie co v1.1 (`pdgId, pt, eta, phi, m, decayType, fromBottom, …, recoMatchIndices, …`).
- **`TruthParticles_*` nadal brak** (`storeTruthParticles: False`; nawet włączone byłyby tylko HF `TruthHFWithDecayParticles`) →
  prawdziwa krotność naładowanych n_ch **nie jest** dostępna w ntuplach MC v1.1 ani v2.0. Jedyne źródło n_ch: generator A3
  (`AFP_ddpred/gen`, kolumna n_ch + `nch_hist`).
- Brak `event_uid` (to kolumna charmpp, nie ROOT); w danych (run, event) jest kluczem, w MC — (plik, wpis).

## 4. Co zmienić w `Q_lxplus/gen_norm`, żeby działał na v2.0
Kod czytający (`gn_root.py`, `effmap_gn.py`) **w większości działa bez zmian** (res(): `_NOSYS`/bez; Tracks, TruthDmesons, DMesons — te same nazwy).
Do zmiany/decyzji:
1. **Ścieżki i listy:** `RAW` → `/eos/user/k/klyczak/charm_raw_v20`; katalogi próbek `{mb_nd,data_428770,data_427929}` zamiast nazw
   kontenerów v1.1 z `config/samples.txt` → nowy `config/samples_v20.txt` albo `MB_DIR/DATA_DIR` z env. `gn_lists` musi filtrować
   `*dpt3to150*.root` (inaczej wciągnie `dpt2to3` = podwójne zdarzenia) i schodzić do podkatalogów (`find` zamiast `ls`).
2. **Liczba plików danych:** check oczekuje 39 (428770) → 97. `DATA_REF_TAGS=000021,000030` (referencja w2 = 441 691 zdarzeń etap1/t4)
   **nie istnieje w v2.0** (inne dzielenie plików) → test akceptacji B1 (c_w 0,935/0,793/1,374/1,201/0,981 bit w bit) jest możliwy
   **tylko na v1.1**. Na v2.0: referencja = lista (run, eventNumber) z plików v1.1 _000021/_000030 (dopisać `--ref-events`
   w `gn_combine`) albo `DATA_SCOPE=all` (= istniejący wariant `w2d_alldata`).
3. **D0 (n_corr_d0, c_w, mb_rates):** szersze okno → liczniki D0 inne niż w referencji (9448 → ~×1,4). Do odtworzenia nominału dodać w `run_mc`
   `& (m > 1750) & (m < 2000)` dla `dt == 10` (w `sel` etap1); D⁺/D* bez zmian. Bez tego B1 dla D0 nie przejdzie.
4. **`legs` (h5 + scores):** nogi idą z h5 v1.1 + nominalnych `scores/*/incl_mc.h5` (charmtag) — **nie ma h5 v2.0 na EOS**
   (mini-ntuple v2.0/v2.1 są na NERSC). Węzeł `legs` zostaje na v1.1; łączenie nóg v1.1 z prawdą v2.0 po (job, eventNumber) jest
   **niebezpieczne** (inne numery jobów) → mapa ε trybu direct liczona nadal z v1.1 (te same zdarzenia mb_nd, więc bez straty).
5. **Nowe zmienne (opcjonalnie, zalecane):** do `ev_*` dopisać `nSelBL`, `leadTrackPtBL/1000`, `nPileupVertices`, `passGRL`, oraz w danych
   `trigPassed/trigPS` wybranych łańcuchów (`mult_common.TRIG_STORE`) — pozwala zrobić wariant w2(nSelBL, leadTrackPtBL) i w2 per trigger
   (patrz §5). Wzór: `mult_extract.py` (ta sama pętla porcji, `resolve_event_branches`).
6. **(c) n_ch → n_trk:** w v2.0 nadal `skipped` (brak TruthParticles) — bez zmian; proxy ε_trk·n_ch zostaje (GEN_mcstat B6).
7. `hist_values`: `metadata`/`DSelection_Vertices` — zostawić try/except (nie sprawdzone na v2.0).
8. Klucz zdarzeń MC: mb_nd v2.0 — `check_event_key_unique` jak dotąd; fd (`eventNumber` od 1 w każdym AOD) — tylko (plik, wpis).

## 5. Uwaga fizyczna wynikła z formatu (ważne dla w2 i ISO_sf)
Strumień **physics_Main** runu 428770 jest zdominowany triggerami nie-minbias: w zdarzeniach z kandydatem (3–150) `HLT_noalg_L1J12`
odpalił w 4,14 M z 5,92 M, `HLT_j20_L1MBTS_2` 1,01 M, `sptrk_pt4_L1MBTS_2` 0,76 M (charmpp v2.0 §5.4). v1.1 nie miał `trigPS` ani
łańcuchów J12 → **dotychczasowe porównania dane/MC (w2, ISO_sf: „n_trkid ≥ 15: 0,71 vs 0,23”) mieszają efekt fizyki z obciążeniem
triggerem jetowym L1_J12** (wybiera aktywne zdarzenia). w2(n_trk, lead_pt) częściowo to wchłania (wymiar lead_pt), ale nie w pełni.
Stąd w MULT schematy triggera (`mult_common.TRIG_SCHEMES`) i run 427929 (physics_MinBias, MBTS bez prescale) jako nominał.
