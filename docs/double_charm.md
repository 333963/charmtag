# Double charm — pomiar par D–D w zderzeniach minimum-bias (ATLAS, dane 2022, ntuple v1.1)

> **Uwaga (20.11.2026).** Ten dokument opisuje metodę i stan z 15.11.2026 (run 428770). Od 19.11.2026 wynikiem głównym jest
> run 427929 (czysty minimum bias): μ = 0,50, Z = 1,60 — zob. `README.md` w katalogu głównym, `results/README.md`
> i `results/NOTA_WYNIK/nota.md`. Liczby poniżej dla 428770 są zapisem historycznym (obserwacja obciążona triggerem L1_J12).

Dokument dla osoby, która wchodzi do projektu od zera: co mierzymy, jak, jakim kodem, co jest aktualne i co jeszcze otwarte. Liczby w tekście pochodzą z notatek roboczych projektu (nazwy w nawiasach, np. *[T]*); gdzie to możliwe podano też nazwę modułu w pakiecie `charmtag`.

Stan na 15 listopada 2026 (wynik: zapis z 13.11, §4). Pary D–D **odsłonięte 5.11 20:03 CET** (§4); liczby z wcześniejszych biegów (§2.x) są ślepe — oczekiwane (Asimov), nie obserwowane.

---

## 1. Problem fizyczny

**Co mierzymy.** Produkcję par mezonów powabnych w jednym zderzeniu pp (double open charm): kanały **D⁰D⁰**, **D⁺D⁰**, **D⁺D⁺** (łącznie z sprzężeniami ładunkowymi), w danych minimum-bias 2022 o niskim pile-upie z tagiem protonu AFP (v1.1 ntupli), z pT(D) w 3–30 GeV. Obserwabla: liczba par znormalizowana do liczby D* w tej samej próbce, porównana z Pythią (MC minbias ND/SD/DD/CD oraz forced-decay):

  μ = N_par(dane) / N_par oczekiwane,  N_par oczekiwane = N_D*(dane) × [pary/D*]_MC × c_w × SF.

μ ≈ 1 oznacza zgodność z Pythią (MPI, podwójna produkcja cc̄); odchylenie mierzy nadwyżkę/niedobór podwójnego powabu.

**Dlaczego to trudne.**
- **S/B par ≈ 0,05 %** przed cięciami: D-mezon w minbias ma S/B rzędu 10⁻³, para — iloczyn. Potrzebny silny klasyfikator pojedynczego D i dwuwymiarowy fit mas (SS/SB/BS/BB).
- **Swapy K↔π w D⁰** (hipoteza mas zamieniona): ta sama para torów daje dwóch kandydatów, prawda MC w v1.1 oznacza oba jako poprawne w 65 % przypadków *[DD]* → potrzebna heurystyka (`closest_mass`) i szablon refleksji w ficie.
- **Wspólne ślady**: dwa kandydaty w tym samym zdarzeniu często dzielą tor (np. jeden pion w dwóch D⁺). Takie „pary” to nie są pary cząstek, a wprowadzają korelację mas m₁–m₂. W v1.1 pędy córek są zapisane po refit wierzchołka kandydata, więc identyfikacja śladu po (pT, η) **nie działa** — to był błąd w kodzie do 4 listopada (sekcja 2.6, *[PF]*).
- **Normalizacja**: liczba par w MC inkluzywnym jest mała (dziesiątki) → osobna produkcja Pythii na poziomie truth (5×10⁸ zdarzeń) *[T]*.
- **Kontrola wydajności**: jedyny czysty sygnał w danych to D⁰ z tagu D* (Δm); fit Δm przy niskim pT jest kruchy *[V]*.

---

## 2. Metoda krok po kroku

### 2.1 Dane i format

Mini-ntuple ROOT z CharmAnalysis/charmpp (gałąź minbias) konwertowane do HDF5 (`events`, 60 kolumn): kandydaci D⁺ (decayType 4, Kππ), D* (1, Δm), D⁰ (10, Kπ), Λc (11); masa i pT w GeV, pędy córek i m(Kπ) w MeV, `Dmeson_truth_category` (0 c, 1 b, ≥8 tło), parametry torów córek (d0, z0·sinθ względem PV, dE/dx). Próbki: dane (`data`), MC inkluzywne minbias (`mb_*`), MC forced-decay (`fd_*`, sygnał). Brak w v1.1: ładunków i pdgId córek, indeksów śladów (dodaje konwersja v1.2 *[H2]*). Audyt: 0 duplikatów kandydatów w v1.1 *[DD]* — klucz audytu (bity pT/η/φ córek) jest wrażliwy na refit tak jak stary identyfikator śladu, więc audyt jest **powtarzany kluczem (d0, z0)** (etap DD2); do tego czasu wynik „0” traktować jako wstępny. Opis kolumn: `docs/data_format.md`.

Kandydatów wrong-sign w ntuplach nie ma — findery odrzucają złe znaki przy budowie *[Z]*; tło kontrolne bierze się z sidebandów.

### 2.1a Okna mas i progi — jedna tabela definicji

| nazwa | D⁰ | D⁺ | gdzie używane |
|---|---|---|---|
| zakres histogramu / fitu masy | 1,75–2,00 GeV | 1,70–2,20 GeV | `species.mass_range`; fit masy, pary (25 binów) |
| **SR** (region sygnału, zliczanie i kontrole) | \|m−m₀\| < 55 MeV | < 45 MeV | `species.signal_window`; ε_B pod pikiem, closure generatorów (projekt B) |
| **SB** (tło treningowe, ślepota par) | \|m−m₀\| > 75 MeV | > 70 MeV | `species.sideband_min`; sidebandy do treningu, `pairs`: centrum = obie nogi w oknie 75/70 (zaślepione) |
| pasmo walidacyjne (nieużywane w treningu) | 40–75 MeV | 45–70 MeV | `species.inner_band`; test rzeźbienia `R_inner_sym`, C2ST dane vs model w projekcie B (55–75 MeV dla D⁰) |
| okno D⁰ w normalizacji D* | \|m(D⁰)−m₀\| < 40 MeV | — | `pairs.dstar_norm.d0_window` |
| tag D* | \|Δm − 145,43 MeV\| < 1 MeV | — | kategorie nóg D⁰ 3/4 |
| WP80 / WP70 / WP50 | próg score o ε_S = 80/70/50 % per bin pT (ważone MC) | | fit masy, ε_B, SF |
| WP90 / WP65 | ε_S = 90 % (D⁰: wejście do arbitrażu hipotez i nóg) / 65 % (D⁺: próg nogi; D⁰: kategoria 3 z tagiem) — progi **globalne w pT**, nieważone (jak w fit2d M7) | | `pairs/build.cat_of` |
| kategorie nóg | przedziały ε_S [0,30), [30,45), [45,65) % + D⁰: tag D* (3: tag ∧ ≥WP65, 4: tag ∧ WP90–65) | | fit 2D |
| biny pT | 3–5, 5–8, 8–30 GeV | | WP, SF, walidacja |

### 2.2 Klasyfikator pojedynczego D — `charmtag` (moduły `features`, `samples`, `model`, `pipeline`) *[L]*

- **Sygnał**: `fd` z truth 0/1. **Tło**: sidebandy danych (|m−m₀| > 75 MeV D⁰ / 70 MeV D⁺), ważone do (pT, |η|) sygnału; dla D⁺ weto D* (min(m − m_Kπ) < 160 MeV, m_Kπ z `Dmeson_mKpi1/2`; usuwa kandydatów D⁺, których podpara Kπ pochodzi z D⁰ z D*). Izolacja (ptcone) jest zakazana jako cecha.
- **Cechy H** (18): log pT, |η|, χ² wierzchołka, Lxy, istotność parametru zderzenia, |z0 sinθ|, max|d0| córek, pointing, |cos θ*|, |IPsig|, min|d0|, max|z0PV|, max|η| córek, dE/dx kaonu (+log p, liczba punktów), max dE/dx pionów, log min p pionów. Wariant F7 (7 cech wierzchołkowych) jako odniesienie.
- **Model**: XGBoost; **k-fold = 5** po hashu (run, event); każdy kandydat punktowany modelem foldu, w którym nie był trenowany (OOF wszędzie — także MC). Progi WP 80/70/50 % wydajności sygnału **per bin pT** (3–5, 5–8, 8–30 GeV).
- **Wynik** (OOF, podpróbki v1.1): ε_B@WP80 ≈ 0,02–0,03 (D⁺) i 0,06–0,08 (D⁰) wobec 0,32–0,36 dla F7 → zysk **×5–6 przy WP80, ×12–24 przy WP50** *[L]*. (Te same liczby w `history.md`.)
- **Kontrola rzeźbienia masy** (`sculpt`, *[X, VAL]*): klasyfikator nie może kształtować piku z tła. Stara metryka R_win (okno/asymetryczne SB) dawała 1,24 dla D⁺ i wyglądała na rzeźbienie — to artefakt: ε_B tła maleje gładko z masą (×2,5 na 0,5 GeV), a SB są asymetryczne. Poprawne miary: `R_win_sym` ≈ 1,01, `pol2_closure` (SB → okno po cięciu) 0,98–1,04, `slope_eff`. **Cechy H zostają bez zmian.** Dodatkowo niezależna kontrola: ręcznie zadany dyskryminant bez cech skorelowanych z masą (`hand`, *[W]*).

### 2.3 Fit masy pojedynczego D (`massfit`)

Binned extended NLL (iminuit): sygnał 2 Gaussy z MC (σ z więzem 10 % do MC), tło exp/pol2, dla D⁰ dodatkowo szablon swapów K↔π (z MC; wariant hybrydowy z danych dla ogonów *[S5]*). Fit bez cięcia i po WP → N_sig, ε z fitu. Model tła exp vs pol2 = systematyka (±3–4 % @WP70 dla D⁺ *[Y]*); szablon pseudo-wrong-sign z zamiany K↔π **nie** nadaje się do fitu (χ²/ndf 3–6) *[Y]*.

### 2.4 Wydajność w danych — kontrola D⁰ z D* i SF (`dstar_refit`; wyznaczanie SF per WP — notatka N, moduł w przygotowaniu) *[V, N]*

D⁰ z tagu D* (|Δm−145,43| < 1 MeV) to czysta próbka do pomiaru ε(WP) w danych. Pułapka: fit Δm jednym Gaussem z zamrożonym kształtem przy pT 3–5 GeV ucieka do granicy σ (2,5 MeV) i zawyża mianownik — pozorna strata 0,43/0,80 dane/MC była artefaktem *[V]*. Poprawny fit (`dstar_refit`): 2 Gaussy z kształtem z MC × skala, tło progowe, **jednoczesny fit pass/fail** z ε jako parametrem. Z niego: SF_D⁰(WP) per bin lub zintegrowane; SF_D⁺ metodą względną ε(WP)/ε(WP90) z domknięciem na D⁰ (referencja wrześniowa: SF_D⁰ 0,88±0,05±0,04, SF_D⁺ 0,95±0,07±0,06). Osobnego modułu SF jeszcze nie ma w repo — liczby SF wchodzą do par przez `pairs.sf_json` (format jak w notatce N) *[N]*.

**Stan SF (4 listopada, N2/B2).** Na statystyce dostępnej lokalnie (ok. 1/12 lumi D⁰) SF jest mierzalne tylko w 8–30 GeV: D⁰ 0,97/0,97/1,14/1,01/1,03 ±0,16–0,21 @WP80/70/50/65/90, D⁺ względnie 1,01/0,98/0,91/0,97 ±0,05 (zgodne z wrześniem). Pełna lumi (histogramy pass/fail D*, B2) daje D⁰ 8–30: **SF_rel(÷WP95) 0,97±0,05 (WP70–90), 0,94±0,05 (WP50)**, trzy metody × dwa okna zgodne w ±0,03; 5–8: 1,0±0,15; 3–5: niemierzalne (przyjąć jak 5–8). Artefakt 0,43/0,80 @3–5 (stary fit) nie potwierdza się nowym fitem (0,57±0,27), ale bin pozostaje niemierzalny. Do par wchodzi `sf_for_pairs_v3.json` = wartości wrześniowe + niepewność 0,07–0,10 per WP (±8–13 % na Z par jako systematyka, bez przesunięcia). **Czerwona flaga:** 5–8 GeV, WP95: ε_dane 0,74±0,07 vs MC 0,95 w oknie ±55 MeV, ≈1 w ±36 MeV, stosunek płaski w WP → ~20 % piku Δm w pierścieniu 36–55 MeV od m_D⁰ H odrzuca jak tło; w MC to asymetryczne rozpady z ogonem masy (|cosθ*| 0,84–0,91), nie swapy. Hipoteza: ogony m(D⁰) w danych szersze niż w MC przy pT < 8 (systematyka okna, nie cech) — rozstrzyga analiza zdarzeniowa par D⁰←D* na pełnej lumi (`ev_dstar`). Swapy K↔π (8 % piku Δm) mają ε_H(WP50) 0,69 vs 0,49 dla czystych D⁰ → SF czystego D⁰ niższy o 0,02–0,03 przy WP50.

### 2.5 Pary — budowa (`pairs/build`) *[O]*

1. Nogi: D⁰ ≥ WP90, D⁺ ≥ WP65 (po wecie D*), masa w zakresie histogramu (D⁰ 1,75–2,00, D⁺ 1,70–2,20 GeV, 25 binów).
2. **Arbitraż hipotez D⁰**: jedna hipoteza K/π na (zdarzenie, para śladów), wybór hashem (ślepo na masę).
3. **Weto wspólnych śladów** (sekcja 2.6).
4. Kategorie score (przedziały ε_S z MC: [0,30), [30,45), [45,65) %), tag D* jako osobna kategoria D⁰.
5. Histogram 2D (m₁, m₂) per kanał × para kategorii. **Ślepota**: komórki z obiema masami w oknie (|m−m₀| < 75/70 MeV) są wycinane przy budowie; w trybie ślepym fituje się nuisance’y przy μ = 1 i raportuje tylko dev/ndf oraz Z_Asimov.

### 2.6 Weto wspólnych śladów — poprawka PF (4 listopada) *[PF, PF2]*

Do 4 listopada identyfikator śladu córki był liczony z bitów float32 (pT, η). Pędy w h5 są **po refit wierzchołka kandydata**, więc ten sam tor w dwóch kandydatach ma inne bity → weto usuwało tylko dokładne kopie (swapy D⁰), **0 % w D⁺D⁺ i D⁺D⁰**. Poprawny klucz: bity **(d0, z0·sinθ_PV)** córki — parametry toru względem PV, niezależne od refitu (`pairs.track_id: d0z0`, domyślnie).

Skutki (½ lumi ∩ podpróbki, ślepo):

| kanał | par tła usuniętych (dane) | corr(m₁,m₂) w SB×SB przed → po |
|---|---|---|
| D⁰D⁰ | +23 % ponad swapy | 0,055 → −0,01 |
| D⁺D⁰ | 33 % | ≈ 0 |
| D⁺D⁺ | 82 % (MC: 91 %) | MC: 0,24 → 0,01; dane: 0,26 → 0,19 (218 par, ±0,07) |

Prawdziwe pary DD (obie nogi truth) flagowane: 0 na 467 k / 276 k, 3 na 40 k → weto nie kosztuje sygnału. Tag D* i normalizacja truth bez zmian. **Korelacja ρ w tle D⁺D⁺ (model tilt+rho z [V]) pochodziła właśnie ze wspólnych śladów.** Wszystkie wcześniejsze wyniki par (sekcja 4) były liczone ze starym wetem i są powtarzane (`pf.dag`).

### 2.7 Fit 2D par (`pairs/model`, `pairs/fit2d`) *[O, V, U]*

Model na histogramie (m₁, m₂) per kanał i parę kategorii: SS = μ·E_SS (szablon sygnał×sygnał z MC i normalizacji); SB/BS = κ·d·E (prawdziwy D × fałszywa noga; d więzione ±10 %); BB wolne per komórka (start z GLM sideband×sideband); fałszywa noga exp(poly2) per gatunek × kategoria; szablon D⁰ = poprawne + swap. κ: „plan B” (per kategoria prawdziwej nogi, nominał). Rozszerzenie *[V]*: **tilt** (nachylenie fałszywej nogi D⁺ zależne od kategorii partnera) + **ρ** (korelacja mas w BB D⁺D⁺) → dev/ndf 0,79/1,27/1,28 zamiast 0,79/1,77/3,06; po PF ρ(m₁,m₂) w SB×SB D⁺D⁺ = 0,003 (mixing, pełna lumi) → tilt+ρ zbędne jako nominał, zostaje cross-check. Wynik: profil μ per kanał, wspólne μ (suma profili), Asimov, toys (bias, pull, pokrycie), przesunięcia Asimov.

### 2.8 Normalizacja (`pairs/truthnorm`) *[T]*

N_par oczekiwane = N_D*(dane, fit Δm w oknie D⁰ ±40 MeV) × [pary/D*]_MC × c_w, gdzie c_w = stosunek (X/D*) ważonego wagami krotności zdarzenia (n_trk, pT wiodącego śladu; MC 09-26) do nieważonego — h5 v1.1 nie ma tych zmiennych, więc efekt wag przenosimy jako czynnik. Stosunek [pary/D*] z produkcji Pythii na poziomie truth (5×10⁸ zdarzeń, A3): D⁰D⁰ 0,01257, D⁺D⁰ 0,0225, D⁺D⁺ 0,01049 na reco D*, stat. MC < 0,2 % (zamiast 27–37 % z MC inkluzywnego); warianty mapy ε ≤ 2 %. c_w = poprawka wag krotności (ważone/nieważone). N_D* liczy `dstar_refit.fit_single` (kształt z MC, σ swobodna; tło próg × pol2 w t=u/20 MeV, całki GL — NDS); stary `fit_dm` (limit σ ≤ 2,5 MeV + sztywne tło, artefakt *[V]*) jest zapisywany obok w `meta.dstar_fit_data.old` tylko do porównania. Na pełnej lumi 2022: N_D* = 4443 ± 612 (stat) ± 591 (syst, warianty tła/zakresu/binów/sygnału), fit_dm dawał 5588 (+26 %) *[NDS]*; biegi `pf_pairs`/`unblind` na src ≤ 2.0.2 z poprawką ID śladu z 5.11 wpadały w fallback fit_dm (fit nie zbiegał) → μ korygowane post hoc ×N_użyte/N_NDS (`fill_unblind.py --ndstar`), Z nie zależy od N_D*.

### 2.9 Systematyki (`syst`) *[U]* — ślepo, Δμ przy μ = 1, pełna lumi, model tilt+rho, stare weto

| źródło | D⁰D⁰ | D⁺D⁰ | D⁺D⁺ | wspólne |
|---|---|---|---|---|
| szerokość/skala masy sygnału | 0,15 | 0,06 | 0,11 | 0,08 |
| model tła (poly3 ↔ tilt+rho) | 0,04 | **0,82** | 0,44 | 0,46 |
| κ plan B ↔ kanał | 0,36 | 0,14 | 1,43 | 0,11 |
| SF D⁰/D⁺/tag | 0,18 | 0,14 | 0,20 | 0,14 |
| normalizacja (MC direct, N_D* 10 %, n_trk, ε) | 0,24 | 0,19 | 0,28 | 0,21 |
| tune Pythii (MPI±, CR±, Monash) | w toku | | | |
| **suma (bez tune)** | 0,49 | 0,87 | 1,54 | 0,55 |

Model tła był źródłem nr 1 w D⁺D⁰ — spodziewana redukcja po PF (korelacja tła znika). Dodatkowe źródła z projektu B: kształt tła pod pikiem z SALAD ±5 % @WP80 / ±10 % @WP70 *[S8]*, exp vs pol2 ±3–4 % *[Y]*. Klastrowanie par D⁺ (wiele kandydatów w zdarzeniu) inflowało wariancję ×1,54 — w dużej części to wspólne ślady; po PF zmierzyć ponownie (bootstrap zdarzeń).

### 2.10 Walidacja i strażnicy

- `validate`: ε_B OOF, AUC, `R_win_sym`, `pol2_closure`, `slope_eff`, ε_B pod pikiem (symetryczne SB / SALAD), kontrola D⁰ z D* (`dstar_refit`), opcjonalnie ręczny dyskryminant *[VAL]*.
- `dedup` (`io.dedup: report|true`): klucz (run, event, decayType, bity **d0/z0 córek w kolejności hipotezy**) — ten sam identyfikator śladu co weto PF, odporny na refit; v1.1: **0 duplikatów / 4,75 M wierszy** obiema metodami (DD, DD2); obie hipotezy K↔π tej samej pary (10 % D⁰ w danych, 24 % w fd) i D⁰ wewnątrz D* (osobny `decayType`) to nie duplikaty; włączyć `report` przy innych produkcjach (v9 miało 10,8 % kopii D⁰ z finderu D*) *[DD, DD2]*.
- `pairs_counts.json`: liczba par przed wetem / ze wspólnym śladem / zachowanych — kontrola, że weto działa.

---

### 2.11 Pary różniczkowe — LS/OS, Δφ, Δy (`pairs/diff`) *[PD]* — rozszerzenie, domyślnie wyłączone

Cel: rozróżnić DPS (dwa niezależne twarde procesy) od SPS (jeden proces, np. g → cc̄) na odsłoniętej próbce (~220 par sygnału).
Jedyna rozstrzygająca obserwabla przy tej statystyce to **podział like-sign / opposite-sign** (znak kwarku c obu nóg: SPS daje tylko OS,
DPS 1:1; LHCb LS/OS ≈ 0,10–0,13) *[INT]*; Δφ w 2 binach (< π/2, ≥ π/2) daje < 1σ i ma sens tylko wewnątrz OS (`sign_dphi`).
- Znak nogi: D⁺ → `Dmeson_charge` (jest w v1.1); D0 → −ładunek K (córka 0) — kolumna `Dmeson_daughter0_charge` tylko w v1.2:
  `pairs.diff.charge_source` (h5 v1.2 — te same 60 kolumn i wiersze; join po `row` z kontrolą run/event/decayType/masy). Bez v1.2
  nogi D0 są „nieokreślone" → LS/OS tylko dla D⁺D⁺. Swap K↔π odwraca znak, więc ładunek K bierzemy z hipotezy bliższej m0 spośród
  rzędów tej samej pary torów (`d0_flavour: closest_mass`; tylko do klasyfikacji — arbitraż, masa i kategoria bez zmian; tak samo w MC).
- Budowa: te same nogi i pary, maska binu PO wecie śladów; histogramy `(n_cat, 25, 25)` per kanał × bin (`pairs_full_diff.npz`),
  SB×SB per bin, oczekiwanie per bin `N_par · f_i` (f_i = udział binu w parach MC minbias Pythii; fallback płasko), `pair_inputs_diff.json`.
- Fit: ten sam model 2D per bin (μ_i, Z_i, 68 %), łączenie po kanałach per bin (suma profili), testy kształtu z profili:
  Δχ² = 2·min_μ Σ ΔNLL_i(μ·h_i/f_i) − Σ 2·ΔNLL_i(min) ~ χ²(n_bin−1) dla h = Pythia (`pythia_mc`) i czysty DPS (`dps_flat`: 1:1 / płasko w Δφ);
  blok like-sign: f_LS = N_LS/(N_LS+N_OS) ± vs Pythia (± binomialnie z N par MC) i LHCb → `fit2d_diff.json`, `tfin_diff.md`.
- Uruchomienie: `charmtag pairs -c CFG --unblind --diff [sign|dphi|dy|sign_dphi]` → `charmtag fit2d -c CFG --unblind --diff`
  (`--diff` wymaga `--unblind`; nominalne wyjścia bez `--diff` bitowo bez zmian — test z hashem). lxplus: `pkg/pd_diff` (DAG 3 węzły).
- Oczekiwana czułość (z σ_μ odsłoniętego fitu): v1.1 (tylko D⁺D⁺) σ(f_LS) ≈ 0,4 → ~1σ; v1.2 (3 kanały) σ(f_LS) ≈ 0,13 → DPS (0,5) vs
  SPS-like (0,1) ≈ 3σ; Δφ 2 biny: σ(f_away) ≈ 0,13 → Pythia 0,6 vs DPS 0,5 < 1σ. Notatka: `PD_pairs_diff/charm-PD-pairs-diff.md`.

---


## 3. Jak uruchomić

Instalacja: `pip install -e .[xgb,fit,dev]` (Python ≥ 3.10; `[flow]` tylko dla generatorów z projektu B). Konfiguracja: skopiować `configs/examples/*.yaml`, uzupełnić ścieżki próbek h5 (`paths.samples`) i `paths.out_dir`; pliki referencyjne normalizacji są prywatne: `pairs.mb_rates_json` (czynniki c_w) jest **wymagany** zawsze (przykład formatu: `configs/examples/mb_rates_2d.example.json`), `truth_norm_json` tylko dla `--norm truth`; bez niego pary liczą się z normalizacją `direct` (stosunki X/D* z MC inkluzywnego h5).

Kolejność (nazwy podkomend wg CLI pakietu; po reorganizacji wszystkie pod jednym wejściem `charmtag`):

```
charmtag inspect  -c cfg.yaml                         # kolumny, próbki, zakresy
charmtag train    -c cfg.yaml --species D0|Dplus      # XGB, 5-fold OOF
charmtag apply    -c cfg.yaml --species D0|Dplus      # score dla data/mb/fd (OOF)
charmtag fit      -c cfg.yaml --species D0|Dplus      # fit masy bez cięcia i po WP
charmtag validate -c cfg.yaml --species D0|Dplus      # ε_B, sculpt, D* control
charmtag report   -c cfg.yaml                         # markdown
charmtag pairs    -c cfg.yaml [--norm direct|truth]   # nogi, weto d0z0, histogramy 2D (ślepo)
charmtag fit2d    -c cfg.yaml [--ntoys N]             # fit ślepy, Asimov, toys
charmtag syst     -c cfg.yaml                         # systematyki ślepe
```

Odsłonięcie: `pairs --unblind` i `fit2d --unblind` — tylko po decyzji analizy; tryb ślepy odmawia bez pełnego histogramu.

Na klastrze (HTCondor): DAG-i i skrypty submitu (trening + apply + fit + pairs + toys; osobny DAG powtórki par bez treningu) nie są częścią repo — zawierają ścieżki prywatne; ich budowa jest opisana ogólnie w `docs/running_on_lxplus.md`. Zasady: wejścia i wyjścia na EOS, scratch kasowany, AFS tylko na logi (kwota 2 GB). Szczegóły: `docs/running_on_lxplus.md`. Podpróbki (≤ 1,5 GB) wystarczą do rozwoju na laptopie; pełna statystyka (~40 M kandydatów) tylko na klastrze.

---

## 4. Stan wyników — co jest aktualne

| data | co | wynik | status |
|---|---|---|---|
| sierpień 2026 | poprzednia iteracja (`double_charm_v11`, modele 7-cechowe HGB/DDT) | NO-GO: ESS MC sygnału 41 przy wymaganym 500 | historyczne |
| wrzesień 2026 | odsłonięcie (klasyfikator D⁺ wariant F, D⁰ H; fit2d plan B; normalizacja z MC inkluzywnego) | **3,6σ** = obserwowana istotność łączna (suma profili trzech kanałów; 2,8–3,8σ z systematykami/jackknife; D⁺D⁰ 2,9σ, D⁰D⁰ 1,9σ, D⁺D⁺ 1,0σ); **2,63σ** = oczekiwana istotność Asimova wspólnego μ (inna wielkość, nie sprzeczność); μ 1,50 [1,10; 1,98]; dane/MC 0,88 −0,36 +0,40 *[M, O]* | **stare weto śladów** → do powtórki |
| 2–3 listopada | port do `charmtag`, ½ lumi (23/47 plików) ślepo | Z_Asimov 1,65σ *[T]* | stare weto |
| 3 listopada | pełna lumi ślepo, model tilt+rho | Z_Asimov 2,27σ (2,07σ z więzami wycentrowanymi), z systematykami ≈1,5σ *[V, U]* | stare weto |
| 4 listopada | poprawka PF; `pf.dag` na pełnej lumi (pairs → fit2d → 50×200 toys → syst U → mixing), **ŚLEPO** | weto `d0z0` usuwa 76,5 % par; SB×SB: D⁰D⁰ 2382→1385, D⁺D⁰ 17744→7263, D⁺D⁺ 92077→14005; **Z_Asimov 1,90 / 2,71 / 1,10 → łącznie 3,49σ (z systematykami kształtu/κ 3,30σ)**, przed PF 1,38/1,60/0,47 → 2,17σ; σ_μ 0,57/0,41/0,90; dev/ndf 0,65/0,94/1,07; toys: σ(pull) 1,01–1,03, bias ≤0,08σ, pokrycie 0,68–0,69; syst/stat 0,56/0,67/0,41; mixing χ²/ndf 0,79/0,94/1,20; korelacja m₁–m₂ w D⁺D⁺ 0,26 → 0,003 (tilt+ρ zbędne) | **ślepo**; bramki: G1–G5, G7, G8, G10, G11 ✔, G6 D⁺D⁰ do wyjaśnienia (usunięto 59 % vs oczekiwane 25–40 % — nominał z głównego DAG różni się nie tylko wetem; test `track_id: pt_eta` fit-only), G9/G12 częściowo; N_D* w tym biegu jeszcze z `fit_dm` (skala μ) |
| 5 listopada, 03–17 | kontrola `pf_ctrl` fit-only `track_id: pt_eta` vs `d0z0` na tym samym kodzie (po poprawce ID śladu z 5.11: ID ≥ 0) | d0z0 Z_A 1,90/2,73/1,10 → 3,50σ ≈ pf_pairs; usunięte SB×SB 42/60/85 %; D⁺D⁰ 59 % = 29,5 % 1-śladowych (pas 21–50 %) + 30,2 % 2-śladowych (D⁰ z pary Kπ kandydata D⁺, legalne) → **G6 ✔**; G9 ✔ (N_D* NDS jako systematyka 19 %), G12 ✔ (hash kodu) | 12/12 bramek → rekomendacja GO; decyzja analizy 17:21 CET |
| **5 listopada, 20:03** | **odsłonięcie** — węzeł `unblind` (`pairs --unblind` → `fit2d --unblind --ntoys 100` → raport) na zamrożonym cfg (identyczny ze ślepym poza `out_dir`; sha256 cfg/src/compare w pliku flagi), bez zmian kodu/konfiguracji po decyzji | **Z_obs = 2,25 / 3,01 / 1,57σ (D⁰D⁰ / D⁺D⁰ / D⁺D⁺); wspólne μ = 1,25 [0,90; 1,51], Z = 4,06σ (p = 2·10⁻⁵)**; oczekiwane 3,49σ; z systematykami kształtu profilowanymi μ 1,25, Z 4,11σ; wspólna systematyka μ (U) ±0,28; więz κ partnera σ=0,3/0,5 → Z 3,42/3,13σ (dolne ograniczenie). μ per kanał 1,16 / 1,16 / 1,39 (σ_stat 0,55/0,42/0,93; dev/ndf 0,69/0,96/1,07). Korekta normalizacji N_D* post hoc (fit2d użył fallbacku `fit_dm` 5584; NDS 4443 ± 612 ± 591): ×1,257 → **μ = 1,57 ± 0,38 (stat) ± 0,30 (N_D*)**, Z bez zmian. Kontrole: bias toys 0,01/0,09/−0,11, σ(pull) 0,81/0,99/1,04, pokrycie 68 % 0,80/0,66/0,64; bez kategorii z tagiem D* μ 1,32/0,94/1,39; κ per kanał Z 2,28/2,79/1,55 | odsłonięte; zapis μ zastąpiony 9–13.11 (wiersz niżej); otwarte wtedy: systematyka tune'u Pythii, pokrycie toys (1000), N_D* bezpośrednio z NDS w następnej wersji `pairs` |
| **9–13 listopada** | kontrole post hoc po odsłonięciu (bez zmian modelu/cfg): N_D* (NDS, poprawka NLL `fit_single`, 12.11), tune Pythii (6 wariantów), jackknife po 200 grupach zdarzeń (13.11) | N_D* = 4436 (±21,9 %); μ = 1,50 +0,41/−0,39 (stat) ± 0,42 (syst, w tym tune ±0,17) ± 0,33 (N_D*); 13.11 stat × 1,18 (jackknife) → **μ = 1,50 +0,48/−0,46 ± 0,42 ± 0,33**; Z = 4,06σ bez zmian, z korelacjami zdarzeń ≈ 3,4σ; pary różniczkowe f_LS = 0,33 ± 0,14 (zgodne z Pythią) | **aktualny zapis wyniku** (źródło: `UNB_WYNIK.md` §1b, §7) |

Dlaczego wrześniowe 3,6σ nie jest ostateczne: (1) weto wspólnych śladów nie działało → tło D⁺D⁺ ×5–10 za duże i skorelowane, D⁺D⁰ ×1,5; (2) normalizacja z MC inkluzywnego miała 27–37 % niepewności statystycznej (teraz < 0,2 %); (3) SF i N_D* z kruchego fitu Δm. Po PF (pełna lumi, ślepo) oczekiwana istotność Asimova wzrosła z 2,17σ do 3,49σ; obserwowana po odsłonięciu 5.11: 4,06σ (μ 1,25 → 1,57 po korekcie N_D*). Wrzesień dawał μ 1,50 przy tle D⁺ ×1,5–10 za dużym; nowe μ bliżej Pythii, z istotnością wyższą dzięki czystszemu tłu.

Dwa „problemy klasyfikatora”, które wyglądały na realne, okazały się artefaktami pomiaru: strata sygnału przy pT 3–5 GeV (fit Δm, *[V]*) i rzeźbienie masy D⁺ (definicja R_win, *[X]*). Cechy H zostają.

---

## 5. Otwarte punkty i ryzyka

1. ~~Wynik PF na pełnej lumi~~ → odsłonięte (§4). Do noty: toys 1000 (pokrycie 68 % 0,64–0,80 przy 100), systematyka tune'u Pythii w normalizacji, N_D* z `dstar_refit` (NDS) bezpośrednio w `pairs` (dotąd korekta post hoc).
2. **SF z danych** (*[N, N2, B2]*): pełna lumi 8–30 zmierzone (0,97±0,05), 5–8 ±0,15, 3–5 niemierzalne; `sf_for_pairs_v3.json` (wrzesień + niepewność) jako systematyka ±8–13 % Z par; anomalia 5–8/WP95 do rozstrzygnięcia zdarzeniowo (`ev_dstar`).
3. **N_D*** nowym fitem — wpływa liniowo na S oczekiwane i μ.
4. **Swapy D⁰**: prawda z `truthPdgId` po konwersji v1.2 *[H2]* zamiast `closest_mass` (65 % par z obiema hipotezami oznaczonymi jako prawda).
5. **Model tła 2D**: po PF dev/ndf 0,65/0,94/1,07 — plan κ wystarcza; tilt+ρ tylko jako cross-check; dyspersja plików D⁺D⁺ (σ_jk/σ≈3) → `mu_err` = max(Hesse·infl., profil).
6. ~~Odsłonięcie~~ wykonane 5.11 20:03 CET po pisemnej decyzji analizy (12/12 bramek). Zasada po odsłonięciu: żadnych zmian modelu/kategorii/WP; dozwolone cross-checki i wykresy.
7. Zamknięte, nie wracać: wrong-sign z v1.1 (tylko Athena), pseudo-WS w ficie, kanał słabo nadzorowany dla par (S/B ≈ 0,05 %), ręczny dyskryminant jako klasyfikator roboczy.

---

## Skróty

| skrót | znaczenie |
|---|---|
| SS / SB / BS / BB | obie nogi sygnał / sygnał×tło / tło×sygnał / obie tło (komórki fitu 2D) |
| WP80/70/50 | próg score o wydajności sygnału 80/70/50 % (per bin pT) |
| OOF | out-of-fold: score z modelu, który nie widział danego kandydata |
| SR / SB | region sygnału (okno masy) / sidebandy |
| ε_S, ε_B | wydajność sygnału / tła po cięciu |
| R_win, R_win_sym | stosunek ε po cięciu okno/SB (asymetryczne) / (symetryczne) — kontrola rzeźbienia |
| pol2_closure | tło w oknie przewidziane fitem SB vs prawda, po cięciu |
| SF | współczynnik dane/MC wydajności klasyfikatora |
| κ, d, μ | parametry fitu 2D: stosunek fałszywych nóg, więz normalizacji, siła sygnału |
| Z_Asimov | oczekiwana istotność z danych Asimova (ślepo) |
| ESS | effective sample size MC (sierpniowe kryterium) |
| PF | poprawka weta wspólnych śladów (klucz d0/z0) |
| AFP | ATLAS Forward Proton (tag protonu) |

Tagi w nawiasach [L], [O], [PF]… odsyłają do notatek roboczych spoza repo; ich chronologię, pytania i wyniki streszcza `docs/history.md`. Źródła: L (klasyfikator), O (pary), T (normalizacja, gotowość), U (systematyki), V (D*, fit2d), N (SF), PF/PF2 (weto śladów, powtórka), DD (duplikaty), VAL (walidacja), X/W (rzeźbienie), Y (pseudo-WS), Z (wrong-sign), H2 (v1.2), M (wrzesień), R (pełna lumi), AA/REVIEW (kod).
