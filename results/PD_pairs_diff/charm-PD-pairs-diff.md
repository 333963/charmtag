# charm-PD — pary różniczkowe: like-sign/opposite-sign, Δφ, Δy (DPS vs SPS)

Etap PD, 2026-11-05 (bez dostępu do lxplus/h5). Baza: snapshot `AA_merge/src` (charmtag 2.0.3, po poprawce ID śladu i NDS z 5.11).
Wyjścia w `out/`: `PD_pairs_diff.diff` (patch na src, tylko addytywne zmiany + nowe pliki), `pkg/pd_diff/` (paczka lxplus),
`PD_docs.diff` (sekcja 2.11 w `docs/double_charm.md`), `tfin_diff_SYNTETYK_przyklad.md` (raport z syntetyku — NIE fizyka), `pytest_final.txt`.

## 1. Cel i werdykt

Cel: na odsłoniętej próbce par (wspólne μ = 1,25 [0,90; 1,51], Z = 4,06σ; S oczek. 58/90/26 par w D0D0/D⁺D0/D⁺D⁺) dołożyć
pomiar różniczkowy rozróżniający DPS (dwa niezależne procesy twarde) od SPS (jeden proces, g → cc̄), **bez zmiany nominału**.

Werdykt (statystyka, p. 5): przy ~220 parach sygnału
- **LS/OS jest jedyną rozstrzygającą obserwablą** (SPS → tylko OS; DPS → 1:1; LHCb LS/OS ≈ 0,10–0,13, czyli f_LS ≈ 0,09–0,12) — zgodnie z INT.
  Z **v1.2** (ładunek K córki 0 dla D0 → wszystkie 3 kanały) σ(f_LS) ≈ 0,13 → separacja „czysty DPS (0,5)” od „SPS-like (0,1)” ≈ **3σ**
  (INT: 3–4σ). Z samym **v1.1** (tylko D⁺D⁺, S ≈ 36 par przy μ = 1,39) σ(f_LS) ≈ 0,4 → ≈ 1σ, niekonkluzywne.
- **Δφ w 2 binach** (< π/2, ≥ π/2): σ(f_away) ≈ 0,13; Pythia (f_away ~ 0,6) vs czysty DPS (0,5) < 1σ; dopiero Pythia-SPS (0,7) ≈ 1,5σ. Δy podobnie.
  Δφ ma sens jako drugi wariant i wewnątrz OS (`sign_dphi`: g → cc̄ przy Δφ < π/2 vs 2→2/DPS), nie jako główny pomiar.
- **Rekomendowane domyślne**: `pairs.diff.var: sign`, wszystkie kanały, `expect: mc`, `d0_flavour: closest_mass`, `charge_source` = h5 v1.2
  (bez niego nogi D0 „nieokreślone” → tylko D⁺D⁺). Δφ: 2 biny, kanały D0D0 + D⁺D0 (D⁺D⁺ tylko łącznie; per bin σ_μ ≈ 1,3).

## 2. Projekt (co liczy kod)

Zmienne pary (symetryczne — kolejność nóg bez znaczenia): Δφ = |φ1 − φ2| zawinięte do [0, π]; Δy ≈ |η1 − η2|; znak zapachu c nogi:
D⁺ → `Dmeson_charge` (±1, jest w v1.1), D0 → −(ładunek K = córki 0) (v1.2 `Dmeson_daughter0_charge`; v1.1 → NaN → „nieokreślony”).
`like_sign` = oba znaki równe (c c lub c̄ c̄ → tylko DPS); D⁺D0 z K⁻ to c c = like.

Biny (`pairs.diff`): `sign` → [like, opposite]; `dphi`/`dy` → z `edges` (domyślnie 2 biny: [0, π/2, π], [0, 1, 5]);
`sign_dphi` → [like, opp:[0,π/2), opp:[π/2,π]]. Maska binu nakładana na pary **po** wecie wspólnych śladów i po ślepocie — ta sama pętla
co histogramy nominalne (`hists.build_pair_hists(diff=)`), więc Σ_bin H_bin + nieokreślone = H nominalne (test).

Oczekiwanie per bin: SS_i = μ_i · N_par · f_i, f_i = udział binu w parach MC minbias Pythii (`mc_rates(diff=)`, ta sama klasyfikacja
co w danych; przy N_MC < 10 albo pustym binie → płasko z notką), SB/BS: n_real × g_i (g_i = udział SB×SB w binie, dolna granica 1/(2 n_bin);
κ wolne per kategoria absorbuje resztę), BB: GLM na SB×SB per bin. Szablony sygnału, Y, kształty fałszywych nóg, N_D* — wspólne (nominalne).

Fit (`diff.run_diff`): per (kanał × bin) ten sam `PairFit` (plan B, pełna płaszczyzna): μ̂_i, Hesse, Z_i (μ=0), 68 % z profilu na siatce
[−1,5; 10] co 0,25. Łączenie po kanałach per bin = suma profili (jak nominalne `combine`). Testy kształtu bez ponownych fitów — z profili:
dla hipotezy udziałów h (Pythia = f_MC; DPS = 1:1 / płasko w Δφ) Δχ² = 2·min_μ Σ ΔNLL_i(μ·h_i/f_i) − Σ_i 2·ΔNLL_i(min) ~ χ²(n_bin − 1);
łącznie w kanałach: model „μ_k per bin wspólne w kanałach z kształtem h” (n_bin par.) vs „jedno μ globalne z h” (1 par.) → Δχ², p.
Blok like-sign: N_LS = μ_LS·N_par,LS, N_OS analog., f_LS = N_LS/(N_LS+N_OS) ± (Hesse, propagacja), f_LS^Pythia ± binomialnie, LHCb 0,09–0,12,
Z odrzucenia Pythii/DPS = √χ²(1) z p testu. Kontrola: μ wspólne we wszystkich (kanał × bin) ≈ nominał.

### Swapy K↔π i znak D0 (ważne)
Udział swapów w szablonach D0 wynosi 27–39 % (kat. 0–2); swap odwraca ładunek „K” → bez korekty f_LS w D0D0 byłoby rozmyte niemal do 0,5
(p_flip ≈ 0,42), w D⁺D0 ×0,4 — a MC (truth-correct) nie. Dlatego `d0_flavour: closest_mass` (domyślnie): ładunek K brany z hipotezy
**bliższej m0 spośród rzędów tej samej pary torów** (obie hipotezy są w h5 jako osobne wiersze), zastosowany identycznie w danych
(`make_legs`, przed arbitrażem) i w MC minbias (`legs_for_truth`: ładowane też truth 100/101, po wyborze odrzucane). Zmienia się
**tylko `dq0`** (klasyfikacja LS/OS); arbitraż (hash, ślepy na masę), masa, kategoria, szablony — bez zmian. Resztkowa dylucja
(swap bliżej m0 niż poprawna hipoteza; szerokość swapu ~70 MeV vs ~15 MeV) wchodzi do f_MC przez tę samą procedurę. `hypothesis` = stare
(ładunek z rzędu po arbitrażu) do porównania.

### charge_source (v1.1 score ↔ v1.2 ładunki)
`docs/data_format.md`: v1.2 = „te same 60 kolumn bit w bit + 19 nowych”, czyli **te same wiersze w tej samej kolejności** → join po `row`
(globalny wiersz z `open_events`). `pairs.diff.charge_source`: katalog z plikami o nazwach jak `paths.samples` **albo** słownik
`{data_none: plik|glob|lista, incl_mc: …, signal_mc: …}`. `_ChargeSource` czyta te same porcje wierszy i **sprawdza** run/event/decayType/masę
w wybranych wierszach (różne pliki/kolejność → `SystemExit` z liczbą niezgodnych; test z przetasowanym plikiem). Brak klucza → NaN (nieokreślone).

## 3. Co zmieniono (pliki / linie; wszystko addytywne, domyślnie wyłączone)

| plik | zmiana |
|---|---|
| `charmtag/pairs/diff.py` (nowy, 491 l.) | spec (`spec_from_cfg`, `labels`, `nbins`, `flat_fractions`), zmienne (`wrap_dphi`, `charm_sign`, `pair_vars`, `diff_bin`), oczekiwania (`fractions_mc`, `fit_fractions`, `diff_inputs_json`, `inputs_for_bin`), fit (`profile_at`, `combine_bins`, `shape_tests`, `run_diff`, `data_fractions`, `ratio_of_bins`, `ls_summary`), raport (`write_report`), `load_diff` |
| `charmtag/pairs/build.py` | `DEFAULTS["diff"]=None`; `scan_source(kin=, charge_src=)` (+ `_ChargeSource`, `charge_files`); `d0_flavour_closest_mass`; `make_legs` (D0: flavour przed arbitrażem, tylko gdy `dq0` w nogach); `mc_rates(diff=)`; `cmd_pairs(diff=)` (wymaga `unblind`; zapis `pairs_{region}_diff.npz`, `pair_inputs_diff.json`, `meta.diff`); `legs_for_export` (+ kolumny kin, gdy są); `legs_for_truth(kin=)` |
| `charmtag/pairs/hists.py` | `build_pair_hists(diff=None)` → przy `diff` trzeci element `{ch__k: hist}` + `counts.diff_kept/diff_undefined`; bez `diff` kod i wynik identyczne |
| `charmtag/pairs/frozen.py` | `sbsb_counts(diff=None)` → per bin; wspólna `_sbsb_from_pairs` |
| `charmtag/pairs/model.py` | `Expectation.expected`: `GLM["b"].get(…, 0.0)` dla pustej kategorii (rzadkie/syntetyczne wejścia; wcześniej KeyError także w nominale) |
| `charmtag/pairs/cli.py` | `--diff [VAR]` (pairs: `cmd_pairs(diff=)`; fit2d: `run_diff` + `write_report`); `--diff` bez `--unblind` → `SystemExit` |
| `tests/test_pd_diff.py` (nowy, 11 testów) | p. 7 |
| `CHANGELOG.md` | wpis „Niewydane — PD” |
| `docs/double_charm.md` (`PD_docs.diff`) | sekcja 2.11 |

Nie dotykane: `columns:` w `default.yaml` (kandydaci kolumny ładunku córki w `diff.D_CHARGE_CANDS`, nadpisywalne `columns.d_charge`) —
bo `samples.load_species` haszuje `cfg["columns"]` do klucza cache i dodanie klucza unieważniłoby cache nominalnego biegu.

## 4. Pliki wyjściowe i raport

`pairs --unblind --diff VAR` → `pairs_full.npz` (bez zmian), `pairs_full_diff.npz` (klucze `D0D0__0`, …), `pair_inputs.json` (+`meta.diff`),
`pair_inputs_diff.json` (spec, etykiety, N par MC per bin i nieokreślone, per bin: `n_pair`, `n_real`, `g_sbsb`, `sbsb`, `f_fit/f_mode/f_mc/f_mc_err/f_flat`, `pairs_kept`),
`pairs_counts.json` (+`diff_kept`, `diff_undefined`, `mc_undefined`), `legs_data.parquet` (+ `eta, phi, charge, dq0`).
`fit2d --unblind --diff` → `fit2d_diff.json` (per kanał: biny {μ, 68 %, Z, dev/ndf, κ, profil}, `tests`, `data_fractions`, `ratio_0_over_1`;
`bins` (łączenie po kanałach); `tests` łączne; `like_sign`; `combined_all_bins`) i `tfin_diff.md` (tabela kanał × bin, testy, blok like-sign).
Paczka dokłada `$PD_OUT/tfin_diff.md` z nagłówkiem spójności z nominałem (Σ μ_i f_i vs μ z `fit2d_unblind.json`).

## 5. Statystyka: domyślne biny i oczekiwana czułość

Z `fit2d_unblind.json`: S = 58,1 / 90,4 / 25,8; σ_μ = 0,55 / 0,42 / 0,93; μ̂ = 1,16 / 1,16 / 1,39 (wspólne 1,25). Efektywne tło w jednostkach
par: B_eff = σ_μ²S² − μS ≈ 950 / 1330 / 540. Dla binu o udziale sygnału f i tła g: σ(N_i) ≈ √(μ S f + B_eff g).

| pomiar | kanały | σ per bin (f = ½) | wielkość | σ | hipotezy | separacja |
|---|---|---|---|---|---|---|
| LS/OS, v1.1 | D⁺D⁺ | σ_μ ≈ 1,3 | f_LS | 0,40 | DPS 0,5 vs SPS-like 0,1 | ~1σ |
| LS/OS, v1.2 | D0D0 + D⁺D0 + D⁺D⁺ | 0,78 / 0,59 / 1,3 | f_LS | 0,13 | DPS 0,5 vs 0,1 (0,15; 0,3) | **3,0σ** (2,6σ; 1,5σ) |
| Δφ 2 biny | D0D0 + D⁺D0 | 0,78 / 0,59 | f_away | 0,13 | Pythia 0,6 vs DPS 0,5 | 0,75σ (0,7 vs 0,5: 1,5σ) |
| Δφ 3 biny | — | σ_μ(f=0,15) ≈ 1,9–2,5 | — | — | — | bez sensu |

Udziały Pythii f_MC per bin będą znane dopiero z biegu (MC minbias ND ma tylko 22/41/15 par po wecie → błąd f_MC ±0,08–0,13; w raporcie
osobno, nie wliczany w μ_i). Flagi pochodzenia DPS/SPS nie ma ani w v1.1, ani v1.2 (`data_format.md`: `truthParentPdgId` córek, nie partonów)
— stąd test „Pythia” = kształt z MC, test „DPS” = 1:1 / płasko. Pomiar LS/OS jest zatem **wykonalny statystycznie tylko z v1.2**;
z v1.1 wynik będzie jakościowy (D⁺D⁺). Jeśli f_LS^dane wyjdzie ≈ 0,1 z σ 0,13 — odrzucenie czystego DPS na ~3σ; ≈ 0,5 — odrzucenie SPS-like.

## 6. Jak uruchomić

Lokalnie: `charmtag pairs -c CFG --unblind --diff sign` (lub `pairs.diff:` w YAML + `--diff`), potem `charmtag fit2d -c CFG --unblind --diff`.
Config przykładowy:
```yaml
pairs:
  diff:
    var: sign                 # sign | dphi | dy | sign_dphi
    expect: mc                # mc | flat
    d0_flavour: closest_mass
    charge_source: /eos/.../h5_v12        # albo {data_none: ..., incl_mc: ..., signal_mc: ...}; brak → tylko D⁺D⁺
    # channels: [D0D0, DpD0]  # domyślnie: sign → wszystkie; dphi/dy → D0D0, DpD0
    # edges: [0, 1.5708, 3.1416]
```
lxplus (`pkg/pd_diff/`, wzór `pf_pairs` + `unblind`): patch PD na `$PF/src`, potem `./submit_diff.sh check`, `VAR=sign CHARGE_SRC=… ./submit_diff.sh submit`
→ `diff.dag`: `pd_pairs → pd_fit2d → pd_report`, wyjścia `$PF_OUT/diff/`; `DRYRUN=1` tylko generuje cfg/DAG i wypisuje komendy. Blokada:
wymaga `$UB_OUT/UNBLIND_OK` + `DONE`. Config = `pf_unblind.yaml` + `out_dir` + `pairs.diff` (skrypt odmawia przy innej różnicy).
**Nic nie uruchomiono** — paczka przetestowana w DRYRUN na atrapie `Q/` i pełnym biegu na syntetyku (`charmtag pairs/fit2d` z wygenerowanym `cfg/pd_diff.yaml`).
Czas na lxplus: pairs 1–2 h (+ skan v1.2), fit2d ≈ 0,5–3 min na (kanał × bin) → 10–20 min.

## 7. Testy

`tests/test_pd_diff.py` (11, wszystkie bez `slow`, ~30 s): `wrap_dphi`/`charm_sign`; `diff_bin` dla dphi/dy/sign/sign_dphi (π domknięte,
nieokreślone); spec (domyślne `sign`, błędy); **bitowa identyczność** histogramów bez `diff` z hashem z kodu sprzed PD
(blind `b6d566d9…`, full `f5234dcd…`) także z dodatkowymi kolumnami nóg; podział binów = całość (Δφ zawsze określone; v1.1 → D0 nieokreślone,
D⁺D⁺ kompletne; ślepo → centrum puste w binach); SB×SB i `mc_rates` per bin sumują się; **odzyskanie μ_i** (D⁺D⁺, 2 biny, μ = 2,5 / 0,6,
toy Poisson z modelu per bin → |μ̂_i − μ_i| < 2,5σ, f_LS ≈ 0,64 oczek., testy i blok like-sign obecne, raport); `combine_bins` (parabola:
μ, 68 %, skalowanie h/f, suma); end-to-end na syntetyku: `--diff` bez `--unblind` → błąd, v1.1 (D⁺D⁺ określone, D0 nie; `legs_data.parquet`
z kinematyką; bez `--diff` histogramy identyczne i brak `meta.diff`), „v1.2” z `charge_source` (D0 określone, 0 niezgodnych) i plik
przetasowany → `SystemExit`; CLI fit2d (`--diff` bez `--unblind`, brak wejść). `d0_flavour_closest_mass` (grupa hipotez, `make_legs`).
Pełny zestaw: **145 passed, 4 skipped** (było 134 + 11), `ruff check --select F,E9` czysto (`out/pytest_final.txt`).
Syntetyk end-to-end przez CLI (`fit2d --diff`, 2 kanały × 2 biny, 4 min): raport `out/tfin_diff_SYNTETYK_przyklad.md` (ładunki losowe → f_LS ≈ 0,5, jak powinno).

## 8. Ograniczenia i otwarte punkty

- v1.1: brak ładunków córek → LS/OS tylko D⁺D⁺; Δφ/Δy pełne. v1.2 wymagane dla D0 (konwersja `charmpp_h5_v12.py` na lxplus, te same pliki).
- Udziały Pythii per bin z 22/41/15 par MC po wecie — duży błąd (nie w μ_i); warto policzyć f_MC z większej próbki MC (np. `truth_norm_json`
  z etapu T, jeśli ma Δφ/ładunki) — hak: `pair_inputs_diff.json.bins[k].f_mc` można nadpisać ręcznie przed `fit2d --diff`.
  **Zrobione 12.11** dla f_LS (prawda Pythii 5·10⁸ zd., `PD_fls_truth/fls_truth.md`) — §9, „Kontrola post hoc 12.11”; ponowny bieg niepotrzebny.
- Swapy: `closest_mass` usuwa dominującą dylucję, resztkowa (swap bliżej m0) w f_MC; nie liczono jej osobno — zobaczyć `n_mc_pairs` vs oczekiwane.
- D⁺D0 „like” (c c) zawiera też pary ze złym znakiem D0 — ta sama dylucja w MC; D⁺D⁺ czyste (ładunek z fitu).
- Nuisance per bin niezależne (osobne fity) — μ_i nieco szersze niż gdyby κ/kształty były wspólne; `combined_all_bins` vs nominał pokaże skalę.
- Δχ² z profili na siatce 0,25 (interpolacja liniowa) — dokładność p ~0,01; dla ndf = 1 Z = √Δχ².
- `model.py`: drobna zmiana robustności (pusta kategoria) — nominał na prawdziwych danych nietknięty (wszystkie kategorie zapełnione).

## 9. Wynik na danych (bieg lxplus v1.2, 2026-11-09)

Bieg: `sign`, h5 v1.2 (`charge_source` = `charm_h5_v12/parts`), fit-only na nominalnych score; EOS `charm_q/pf_pairs/diff_v12/`,
lokalnie `PD_pairs_diff/runs/lxplus_v12/` (`tfin_diff.md/json`, `pairs/`). Ładunek córki 0 sprawdzony: 0 niezgodnych (fd 829 070, dane 2 563 101,
mb 22 848 wierszy). **Nogi D0 określone** — 0 par nieokreślonych (bieg 1 z 05.11, `diff/`: 38 128 nieokreślonych, tylko D⁺D⁺).

| kanał | N_LS sig | N_OS sig | f_LS dane | f_LS Pythia (MC) | Z odrz. Pythii | Z odrz. DPS (0,5) | μ like | μ opposite |
|---|---|---|---|---|---|---|---|---|
| D0D0 | 31,5 | 112,5 | 0,22 ± 0,27 | 0,36 ± 0,10 | 0,59 | 1,16 | 0,96 ± 1,42 | 1,95 ± 0,91 |
| D⁺D0 | 83,2 | 158,2 | 0,34 ± 0,18 | 0,32 ± 0,07 | 0,26 | 0,89 | 1,78 ± 1,25 | 1,57 ± 0,62 |
| D⁺D⁺ | 39,4 | 40,6 | 0,49 ± 0,35 | 0,20 ± 0,10 | 0,75 | 0,11 | 4,48 ± 4,47 | 1,15 ± 1,15 |

- Łącznie (fit kanały × biny): Z_odrz. Pythii = 0,00σ, Z_odrz. czystego DPS = 1,32σ. Test `pythia_mc`: μ = 1,50 [1,20; 2,04], Z = 4,21σ,
  Δχ²(per bin vs kształt) = 0,00/1. Test `dps_flat`: μ = 1,64 [1,22; 2,09], Z = 4,07σ, Δχ² = 1,74/1 (p = 0,19).
- Biny łącznie: like μ = 1,50 [0,69; 2,53] (Z 1,83), opposite μ = 1,50 [1,16; 2,09] (Z 3,79).
- Średnia ważona f_LS z 3 kanałów (odwrotne wariancje, bez korelacji — rachunek poza raportem): **f_LS = 0,33 ± 0,14**
  (Pythia ważona ≈ 0,31). Od czystego DPS (0,5) ≈ 1,2σ, od LHCb (0,09–0,12) ≈ 1,6σ. σ ≈ 0,14 zgodne z prognozą z p. 5 (0,13).
- **Werdykt:** dane zgodne z Pythią (SPS+MPI) w kształcie LS/OS; ani czysty DPS, ani SPS-like (LHCb) nie są odrzucone (< 2σ).
  Pomiar LS/OS wykonalny, ale przy tej statystyce niekonkluzywny — nie rozstrzyga σ_eff.
- **Normalizacja:** ten bieg ma N_D* = 4022 (nominał 5584) → μ_i zawyżone ≈ ×1,39. Potwierdza to kontrola spójności: Σ μ_i f_i / μ_nom =
  1,59/1,16 = 1,37 (D0D0), 1,64/1,16 = 1,41 (D⁺D0), 1,82/1,39 = 1,31 (D⁺D⁺). **f_LS i Z testów kształtu są stosunkami — korekta ich nie zmienia**;
  μ_i i μ globalne nie są porównywalne z nominałem bez przeskalowania (nie robić tego ręcznie bez ustalenia wspólnego N_D*, por. NDS 4443).
- **Normalizacja — rozstrzygnięte (etap ndstar, 09.11, `ndstar/ndstar.md`):** wspólna baza N_D* = 4443 ± 612 ± 681 (NDS + syst źródła kształtu MC ±339).
  4022 w tym biegu = ten sam `fit_single`, ale kształt Δm z MC potoku par (k 2,11 vs 2,73) → −8 %. μ ∝ 1/N_D* dokładnie (κ wolne pochłaniają skalę SB/BS),
  więc k_PD = 4022/4443 = 0,905. Na bazie 4443: **μ_PD = 1,45 [1,08; 1,84] vs μ_nom = 1,49 [1,10; 1,90]** (minima z gładkiego profilu — wartości
  1,50/1,25 z JSON są skwantowane do siatki 0,25); spójność Σμ_i f_i/μ_nom = 0,99/1,01/0,95. μ per (kanał × bin) na bazie 4443:
  D0D0 0,87/1,77, D⁺D0 1,61/1,42, D⁺D⁺ 4,05/1,04 (like/opp); like łącznie 1,43 [0,62; 2,29], opposite 1,46 [1,05; 1,89]. f_LS, Z, testy kształtu bez zmian.
- Uwagi: f_MC z 8/14, 13/28, 3/12 par MC minbias (like/opp) — błąd f_MC ±0,07–0,10 nie jest wliczony w μ_i. D⁺D⁺ like praktycznie niemierzalne (μ ± 4,5).
- Incydent biegu: węzły condora (`getenv=False`) nie dostawały `PD_OUT` → krok report czytał stary `diff/` i tam zapisał tfin/DONE (nadpisane
  `diff/tfin_diff.*`, `diff/DONE`, `diff/pairs/fit2d_unblind.json`; nie sprzątane). Naprawa w paczce: `diff.sub` `environment="$(penv)"`,
  `make_dag_diff.sh` przekazuje PD_OUT/UB_OUT/VAR/CHARGE_SRC/PD_NTOYS, `submit_diff.sh report` (DAG tylko z krokiem report). Report powtórzony 05:06 CET.
- **Kontrola post hoc 12.11 — f_LS z prawdy Pythii** (`PD_fls_truth/fls_truth.md`; lokalnie, bez nowego biegu; nominał i wyniki PD bez zmian):
  - `LS_fraction` w `truth_norm.json` to czysta prawda generatora (pT > 3 GeV, |η| < 2,5, z b → c, bez ε); wersja ważona ε (473 pliki A3, ≈ 4,7·10⁸ zd.,
    mapa `combined_A3/eff_map.npz`): **0,187 / 0,189 / 0,184** (D0D0 / D⁺D0 / D⁺D⁺), ± 0,018 / 0,017 / 0,020 (stat ⊕ ε-tag ⊕ tune; dominują MPI i CRqcd).
  - Dylucja znaku D0 przez `closest_mass` (zła hipoteza K↔π, w = 0,13, zakres 0,11–0,16) → oczekiwane f_LS obserwowalne
    **0,329 ± 0,031 / 0,270 ± 0,025 / 0,184 ± 0,020**, śr. ważona 0,273 ± 0,026. f_MC biegu (0,36 / 0,32 / 0,20) zgodne: pulle +0,33 / +0,62 / +0,15σ
    — duża różnica D0D0 vs D⁺D⁺ w MC to dylucja, nie fizyka.
  - Dane (0,22 ± 0,27 / 0,34 ± 0,18 / 0,49 ± 0,35; śr. 0,33 ± 0,14) vs prawda + dylucja: pulle −0,41 / +0,41 / +0,88, łącznie **+0,43σ** (χ²(3) = 1,10);
    z zapisanych profili fit2d 0,44 / 0,44 / 0,78, łącznie Δχ² = 0,04 → **0,21σ** (kontrola post hoc, nie nowy wynik).
  - **Werdykt §9 bez zmian:** zgodność z Pythią; czysty DPS 1,2σ, LHCb-podobne (0,09–0,12) 1,6σ. Ponowny `fit2d --diff` z nadpisanym `f_mc`
    niepotrzebny (f_LS danych niezmiennicze; Z „Pythia” przeliczone z profili). Oczekiwanie Pythii podawać odtąd z prawdy (σ ≈ 0,02–0,03 zamiast 0,07–0,10).
  - Ograniczenia: prawda × mapa ε bez ε klasyfikatora (WP) i SF (≲ 0,01); udział b → c ± 10 pp → ± 0,015–0,018 (niewliczone); w z MC fd przy WP 0,8/0,7/0,5, PD używa 0,9.
