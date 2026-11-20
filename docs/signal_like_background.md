> **Przeniesione 15.11.2026:** projekt B ma własne repo **charmbkg** (pakiet `charmbkg`) — tam jest aktualna wersja tego dokumentu,
> HYP (`docs/hyp.md`) i stan otwartych punktów. Ten plik to kopia ze stanu 4.11 (charmtag 2.0.3), zostawiona dla odnośników z historii.

# Projekt B — obiekty „podobne do sygnału”: tło pod pikiem bez zaglądania do sygnału

Stan: 2026-11-04 (po B1/B2). Kod: pakiet `charmtag.bkggen` (CLI `charmtag-bkggen`). Dokumenty źródłowe w Projekcie: C, S, S2, S4, S5, S6, S7, S8, S9, W, X, Y, Z, F.

## 1. Problem

Mierzymy mezony D (D⁰→Kπ, D⁺→Kππ, D*⁺, Λc⁺) w danych minimum-bias ATLAS 2022 (ntuple v1.1), a docelowo pary D–D. Pod pikiem masy (SR) sygnał to ułamek procenta kandydatów; resztę stanowi tło kombinatoryczne i refleksje (błędna hipoteza K↔π). Klasyfikator i fit potrzebują wiarygodnego opisu tego tła **w SR**, a SR jest zaślepiony (dane) albo statystycznie ubogi (MC). Pomysł: wytworzyć obiekty, które od prawdziwego tła pod pikiem różnią się **jednym, kontrolowanym składnikiem** — masą (sidebandy), hipotezą (swap), pochodzeniem (MC przeważone), jednym torem (rotacja/mixing) — i sprawdzić, czy w pozostałych cechach są nieodróżnialne od prawdy.

Trzy cele, w kolejności trudności:
1. **Trening** — tło z danych zamiast MC (unikamy mismodellingu).
2. **Fit / odjęcie** — kształt i poziom tła pod pikiem po cięciu klasyfikatora.
3. **Zysk Z** — czy lepsze tło przekłada się na istotność pojedynczych D i par.

## 2. Katalog metod

Wspólne metryki: **C2ST** (AUC klasyfikatora „kontrola vs prawdziwe tło w SR” na MC; 0,50 = nieodróżnialne, cel < 0,52), **ε_B ratio** = ε_B(kontrola)/ε_B(prawda) po cięciu na WP 80/70/50 % ε_S, χ²/ndf kształtu masy, **Zproxy** = ε_S/√ε_B. „Prawda” = MC minbias z `truth_category ≥ 8` w |m−m₀| < 55 (D⁰) / 45 MeV (D⁺).

| Metoda | Czym różni się od tła w SR | Założenie | Wynik (D⁰, o ile nie zaznaczono) | Status |
|---|---|---|---|---|
| **Sidebandy (SB)** | masą | cechy H słabo zależą od m | C2ST 0,55; ε_B ratio oba SB 1,02/1,07/1,32 @80/70/50; bliski SB zanieczyszczony swapami i truth 2–7 | baza |
| **SB daleki + weto refleksji** \|m_swap−m₀\|>25 MeV | masą, bez refleksji | jw. | 0,91/0,93/1,06; **najlepsze tło treningowe** (S, S7 odtworzone) | **nominał treningu** |
| **Szablon refleksji z danych** (hybryda S5: ogony z danych, centrum z MC skalowane) | hipotezą K↔π | refleksja = ta sama para torów z zamienionymi masami | wynik S „N 1389 vs 1443, χ²/ndf 1,2” **obalony w S7** (z błędem Poissona szablon z MC na ½ statystyki to 138 ± 184, χ²/ndf 4,2 — szablonu z danych nie ma); zostaje hybryda S5: domknięcie normalizacji 0,98, Z ±1 %; właściwy szablon z `hypo_swap`/`truthPdgId` po v1.2 (H2) | systematyka fitu D⁰ (hybryda S5); szablon po v1.2 |
| **Mixing nóg** (pary z różnych zdarzeń) | korelacjami same-event | tło par = iloczyn nóg | closure ABCD 1,002/1,003, corr(m1,m2)=0,001; ale Δφ i korelacje realne złamane: C2ST vs prawdziwe pary 0,63–0,91 (S9) | kontrola BB, nie tło odniesienia |
| **SALAD** — reweighting MC→dane klasyfikatorem warunkowanym m, interpolacja do SR (2001.05001) | pochodzeniem (MC) | mismodelling gładki w m **i gładki skład tła w m** | **closure C2ST 0,507 (D⁰, S4) / 0,510 (D⁺, S8)**, podłoga 0,50; ε_B(SR) pred/true 0,98/0,94/0,88 (D⁰), 0,98/0,93/1,04 (D⁺); na danych C2ST w paśmie 0,510–0,515. **B1 (closure adwersarialny)**: skład kategorii zniekształcony lokalnie pod pikiem (+4,7 pp MisMatch*) jest niewidoczny dla C2ST (0,507 = podłoga), a ε_B pred/true 0,85/0,81/0,70 @80/70/50 | **v1.0 zamrożony: estymator ε_B + systematyka** (§3a) |
| **FETA** — flow transportujący MC→dane w SB, zastosowany do MC w SR (2212.11285) | pochodzeniem, transport ciągły | flow bazowy p_MC(z\|m) < 0,52 | 0,63 → 0,60 (CPU); limit = flow bazowy (0,618) | po treningu GPU |
| **Flow p(x\|m) + CATHODE** (2109.00546) | generowane z SB | flow odtwarza więzy (derive_H) | closure 0,61 (CPU), nie bije SB; CATHODE na danych przy S/B ≈ 0,2 % ślepy; po preselekcji (S/B ≈ 3,8 %) z tłem odniesienia SALAD widzi pik D⁰: Z 5,0 vs 3,7 sama preselekcja, ale < nadzorowany XGB | tylko słabo nadzorowany kanał |
| **candrot** — rotacja φ, odbicie η, swap K↔π, mixing toru (S6) | jednym torem | wierzchołek pasuje do nowej kinematyki | NIE: zachowany wierzchołek nie pasuje → C2ST 0,7–0,97, ε_B ×0,2–0,7, Z −7…−32 %; D⁰ swap = duplikaty hipotez (99,99 %) | odrzucone jako tło treningowe |
| **pseudo-WS D⁺** (swap K↔π_j) jako kształt tła w fit | hipotezą | kształt m jak tło | S6: N_pred/N_true 0,99/0,97/0,97 (lepiej od pol2 1,00/1,09/1,12); **Y: w pełnym ficie χ²/ndf 3–6, N_sig ×1,4–2,5, `mix` wybiera f=1 (czysty pol2)** — metryka S6 nie testowała kształtu | obalone |
| **Wrong-sign z ntupli** | ładunkiem | WS zapisany w ntupli | findery CharmAnalysis odrzucają złe znaki przed zapisem (`D0Finder.cxx:28`, `DplusFinder.cxx:30/54`, `DstarFinder.cxx:31/170`) | wymaga produkcji Athena |
| **Ręczny dyskryminant** (wagi pisane z fizyki, 0 etykiet; W) | — | fizyka wystarczy | ε_B@80 D⁺ 0,108 vs XGB 0,021 (×5 gorzej); R_win 1,00±0,02 → kontrola rzeźbienia i prior do fine-tuningu (D⁺ 0,17→0,108 po 2 epokach) | narzędzie kontrolne |
| **Dekorelacja masy** (X) | — | klasyfikator rzeźbi pik | R_win 1,24–1,30 = artefakt definicji (asymetryczne SB × gładkie nachylenie ε_B(m) ×2,5/0,5 GeV); symetryczne R_win 1,01/0,99/1,04, closure pol2 0,98/1,00/0,94 → **H bez zmian**; CDF-flattening kalibrowany na danych SB tworzy fałszywą strukturę | zamknięte, metryki w `validate` |

## 3. Co z tego wynika

**Hipoteza „lepsze tło → większe Z” dla klasyfikatora XGB H nie potwierdziła się.** Zproxy @WP70 (D⁰, MC closure) z tłem treningowym:

| tło treningowe | SB (nominał) | SB daleki+weto | SALAD | FETA | flow S2 | kopuła | idealne |
|---|---|---|---|---|---|---|---|
| Zproxy | 4,93 (S4) / 3,45 (S) | 3,43 vs 3,46 (S5) | 4,96 | 4,83 | 3,20 (−7 %) | 2,56 (−26 %) | 4,95 |

Klasyfikator jest nieczuły na mismodelling tej skali: SB = SALAD = ideał w granicach ±1 %. Na danych (fit po cięciu, S4): SB 1,6/2,7/3,5σ, FETA 1,6/2,4/3,4, SALAD 1,1/2,3/3,1 — brak zysku.

Właściwe zastosowania generatora:
- **Systematyka kształtu tła w fit**: SALAD vs SB ±5 % @WP80, ±10 % @WP70 (S8), ±25 % @WP50 (B1; składnik zlokalizowany pod pikiem przechodzi do ε_B bez korekty: −3,3/−4,1/−6,4 % na +1 pp MisMatch*); exp vs pol2 ±3–4 % (Y); refleksja hybryda S5 vs MC ±1 % Z.
- **Estymator ε_B pod pikiem**: proxy z asymetrycznych SB zaniżony o ~20 % (ε_B tła spada z m ×2,7); symetryczne SB lub SALAD dają 1,0–1,1 prawdy (S8, VAL).
- **Tło odniesienia dla CATHODE** po preselekcji podnoszącej S/B ≳ 3 % — SALAD ≈ idealne (SIC 7,2 vs 6,4 @S/B 3 %); flow od zera wymaga ≳ 3 % i AUC closure < 0,52.
- **Mixing** jako kontrola par (ρ(m1,m2), zamknięcie BB), nie jako tło odniesienia — S9: C2ST 0,63–0,91 vs pary same-event.
- **Anomaly detection dla par D–D — zamknięte** (S9): S/B par w danych ≈ 0,05 % (ciasne okno 0,6 %), 10–20× poniżej progu SIC > 1; fit 2D wyczerpuje informację (score + kinematyka: +1–5 %).

Efekt uboczny prac: dwa „problemy klasyfikatora” okazały się artefaktami pomiaru (R_win — X; kontrola D* przy 3–5 GeV — `fit_dm`, V), a S9 ujawnił niedziałające weto wspólnych śladów w parach (PF, projekt A).

## 3a. Sufit informacyjny tła i strona sygnałowa (B1/B2, 4 listopada)

**Sufit po stronie tła.** Zproxy @WP70 z tłem treningowym SB 4,93, SALAD 4,96, **idealne (prawdziwe) tło 4,95** — klasyfikator H trenowany na sidebandach danych jest już na suficie informacji dostępnej z tła. Kolejne generatory tła (FETA na GPU, flow, CATHODE, kopuła) nie mogą dać Z — zamknięte jako dźwignia.

**SALAD v1.0 — zamrożony artefakt** (`config/bkggen/salad_v1.yaml`, `charmtag/bkggen/salad_v1.py`, `make_generator("salad_v1")`): kind `hgb` (produkcyjny; `sk_mlp` tylko cross-check — rozrzut ε_B między podpróbkami ~10 %), 18 cech H, t = (m−m₀)/0,1, clip wag 0,05–20, seed 0, truth ≥ 8; wyjście `salad_w_<sp>_<sample>_v1.npz` z sha256 wejść, n_eff/N, AUC_SB, C2ST pasma i flagą `ok`. Na danych D⁺ (podpróbka): AUC_SB 0,643, n_eff/N 0,67, C2ST pasma 0,507.

**Closure adwersarialny (B1, D⁺).** Pseudo-dane = MC ze zniekształceniem *składu* kategorii tła, nie gładkiej funkcji m:

| zniekształcenie | C2ST SALAD SR | ε_B pred/true @80 / 70 / 50 |
|---|---|---|
| gładkie w m (jak S8) | 0,520 | 1,01 / 1,03 / 1,10 |
| skład globalny (MisMatch* ×1,5 wszędzie) | 0,501 | 0,91 / 0,86 / 0,80 |
| **skład lokalny pod pikiem** (+4,7 pp MisMatch*, σ = 35 MeV) | **0,507** | **0,85 / 0,81 / 0,70** |
| próg w m (lewy SB ×1,5) | 0,503 | 0,94 / 0,93 / 0,85 |
| cięcia fizyczne lokalnie (Lxy, IPsig) | 0,536 | 0,92 / 0,95 / 0,94 |

Wnioski: (1) **C2ST jest ślepy na skład tła** — jest warunkiem wstępnym, nie testem; (2) czułość ε_B: −3,3 / −4,1 / −6,4 % na +1 pp MisMatch* pod pikiem; (3) jedyny test, który skład wykrywa, to ε_B(SR) z fitu pol2 danych po cięciu — na danych SALAD/dane = 1,02±0,03 / 1,09±0,04 / 1,25±0,13 (hgb), co ogranicza nieznany skład do |Δf_MM| ≲ 1,5 pp @80; (4) SALAD uczy się sygnału z krawędzi SB (2 % domieszki w 70–90 MeV → +2,4/+3,8/+6,5 % ε_B, C2ST nie reaguje) → weto refleksji/D* przed fitem obowiązkowe, SB_min 70 MeV. **Systematyka kształtu tła: ±5 % @WP80, ±10 % @WP70, ±25 % @WP50.** Warunek użycia v1.0: zgodność ε_B z fitem danych w ±5/10 %; C2ST ≤ 0,53 tylko wstępnie.

**Strona sygnałowa (B2): D⁰ z tagu D* jako sygnał z danych.** Jedyna próbka sygnału z danych to D⁰ z Δm (sWeights/COWs per bin pT). Pełna lumi, 8–30 GeV: wszystkie 18 cech H |Δμ| ≤ 0,2 σ_MC (abs_cos_theta_star −0,35 σ, zależne od poziomu tła α); SF_rel(÷WP95) = 0,97±0,05 (WP70–90), 0,94±0,05 (WP50), trzy metody × dwa okna zgodne w ±0,03. Bin 5–8 niemierzalny w cechach (±0,15 w SF), 3–5 niemierzalny. Górna granica zysku z reweightingu MC→dane-sW i retreningu H: ≤3–5 % @WP≥70, ≤5–8 % @WP50 — poniżej błędu SF → **NO-GO dla reweightingu sygnału (B3)**; SF_rel wchodzi jako systematyka ε_S. Zamknięte także: augmentacja generatywna sygnału (2 M sygnału MC; zysk GANplification ≈ 0).

**Otwarte (EV):** w 5–8 GeV ε_dane(WP95) = 0,74±0,07 vs MC 0,95 w oknie ±55 MeV, ale ≈1 w ±36 MeV, stosunek płaski w WP → ~20 % piku Δm w pierścieniu 36–55 MeV od m_D⁰ to kandydaci, których H odrzuca jak tło. W MC są to asymetryczne rozpady (|cosθ*| 0,84–0,91 vs 0,52) z ogonem masy, nie swapy. Hipoteza: ogony m(D⁰) w danych przy pT < 8 szersze niż w MC → systematyka okna masy, nie cech; rozstrzyga analiza zdarzeniowa par D⁰←D* na pełnej lumi (paczka `ev_dstar`, DAG < 1 h).

## 4. Jak uruchomić

Po reorganizacji (REVIEW) wszystkie generatory mają jedno API:

```python
from charmtag.bkggen.api import make_generator, available   # available() → ['candrot', 'copula', 'flow', 'salad', 's_copula', 's_flow']
from charmtag.bkggen.c2st import c2st
gen = make_generator("salad", species="Dplus", kind="sk_mlp")   # sklearn, bez torcha; kind="mlp" = torch (S4)
gen.fit(SB_df, mc=MC_df)          # SB_df: dane w sidebandach (cechy H + "m"); MC_df: MC tła w całym zakresie mas
syn = gen.sample(m_sr)            # SALAD: prawdziwe zdarzenia MC w SR z wagą kolumny "w"; copula/flow: cechy generowane
auc = c2st(syn, truth_sr, feats=gen.feats_out, wa=syn["w"].to_numpy())   # closure (MC); na danych: pasmo walidacyjne
# FETA: charmtag.bkggen.feta.FetaTransport (wymaga [flow]); próbki kontrolne SB/refleksje/mixing: charmtag.bkggen.control
```

CLI (sidebandy, refleksje, D*, mixing, generatory, metryki):

```
charmtag-bkggen --method sideband --which far -c config.yaml --sample data_none     # próbki kontrolne SB (+ weto refleksji wg configu)
charmtag-bkggen --method generate --model copula|flow -c config.yaml                 # generatory S2 (SALAD: przez API lub `charmtag s8salad`)
charmtag-bkggen --method candrot --species Dplus --npz kandydaci.npz -o kontrole.npz  # obiekty kontrolne z kandydatów (S6)
charmtag-bkggen --method validate -c config.yaml                                     # C2ST, ε_B ratio, χ² m
charmtag s8salad -c config.yaml --species Dplus --kind sk_mlp                        # SALAD → syst/bkg_shape_salad_<sp>.json (±ε_B pod pikiem)
charmtag-bkggen --help   # pełna lista metod i opcji
```

Użycie w charmtag: `training.bkg_source: sideband|sideband_far_veto|mc|synthetic` + `training.bkg_options.synthetic_model: copula|flow|salad` (SALAD = MC tła ważone, bez torcha); systematyka kształtu tła po cięciu: `charmtag s8salad …` albo `validate.salad: true` w configu (`pipeline.cmd_validate` → `s8salad.shape_syst`) → `syst/bkg_shape_salad_<sp>.json`.

Wejścia: npz/h5 kandydatów z cechami H (`charmtag.features.compute`), próbka danych z sidebandami i MC minbias (truth ≥ 8). Flow/FETA wymagają extras `[flow]` (torch, zuko) i GPU dla pełnej statystyki; SALAD i kopuła działają na CPU w minuty.

## 5. Ryzyka metodyczne (S7 — przegląd niezależny)

1. **Wyciek sygnału do SB**: bliski SB zawiera swapy i prawdziwy charm (truth 2–7) → ε_B ratio 1,15 zamiast 1,07; weto refleksji obowiązkowe, gęstość masy fitować po wecie.
2. **Rzeźbienie przez sam generator/klasyfikator AD**: fit masy po CATHODE z SR ocenianym OOF, a SB ensemblem, dał fałszywy pik „16σ” przy 124 prawdziwych D⁰; CDF-flattening kalibrowany na danych SB — to samo. Zasada: SR i SB tym samym klasyfikatorem, test na MC bez sygnału przed każdym fitem.
3. **Heurystyka `closest_mass`** dla swapów D⁰ jest definicyjna (65 % par torów ma obie hipotezy truth-matched); szablon z MC na ½ statystyki jest pusty (138 ± 184). Rozwiązanie: `truthPdgId` córek w h5 v1.2 (H2).
4. **Statystyka**: SIC/SB przy top 1 % z B rzędu kilku zdarzeń jest bez sensu — podawać przedziały (bootstrap) i tylko przy B ≥ 20; ε_B ratio @WP50 ma ±3–4 % stat.
5. **Więzy kinematyczne**: generowanie cech H wprost łamie więzy (C2ST 0,96–0,999); generować zmienne bez więzów i wyprowadzać H (`derive_H`, `kin3` dla D⁺).
6. **Kontrole candrot zachowują wierzchołek** → dobre tylko do testów masy/kinematyki, nigdy do cech wierzchołkowych.

## 6. Dalsze kroki

- SALAD w pełnej luminozji jako węzeł `salad_{D0,Dplus}` w DAG S5 (`S5_SALAD=1`) → liczba systematyczna do raportu.
- SALAD v1.0 D⁰: closure adwersarialny (kod gotowy); bin 8–30 na pełnym MC (lxplus).
- Anomalia 5–8/WP95 zdarzeniowo na pełnej lumi (`ev_dstar`); SF_rel (B2) do N/U jako systematyka; po v1.2 f_MM w danych wprost (`hypo_swap`) → czułość B1 zamienia się w systematykę liczoną.
- FETA (s4) i flow S2 na pełnej statystyce — wyłącznie cross-check SALAD w SR; nie jako nominał.
- `validate`: symetryczny ε_B / SALAD jako nominalny estymator ε_B pod pikiem (wdrożone w VAL).
- Po v1.2: szablon refleksji z prawdy (`hypo_swap`) zamiast heurystyki.
- Nie wracać do: AD dla par, pseudo-WS w fit, rotacje/mixing toru jako tło treningowe, ręczny dyskryminant jako klasyfikator, WS z v1.1, **nowe generatory tła jako dźwignia Z, augmentacja generatywna sygnału, retrening H na D*-sW (B3) bez nowej statystyki**.

## Literatura

- CATHODE: 2109.00546 · CURTAINs: 2203.09470 · FETA: 2212.11285 · SALAD: 2001.05001 · F4F: 2305.04646
- DisCo: 2001.05310 · ABCDisCo · event mixing w pp (STAR 1204.4244) · ATLAS dE/dx: 1005.5254, 2205.06013
- LHCb double charm 1205.0975 · dyfrakcyjny charm z AFP: 2402.09102
