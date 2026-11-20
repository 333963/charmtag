# Pary mezonów powabnych D–D w zderzeniach minimum bias pp, √s = 13,6 TeV (ATLAS 2022) — szkic noty wynikowej

Każda liczba ma odnośnik do pliku źródłowego (ścieżki względem `charm_work/`). Wszystkie pola wypełnione (stan §9).

---

## 1. Wstęp

Pary hadronów powabnych w jednym zderzeniu pp powstają przez:
- **SPS** (single parton scattering), czyli jedno twarde rozproszenie dające c c̄: kreacja zapachu (c i c̄ w przeciwne strony, |Δφ| ≈ π), rozszczepienie gluonu g → c c̄ (pary kolinearne, małe Δφ, Δη) i ekscytację zapachu;
- **DPS** (double parton scattering), czyli dwa niezależne rozproszenia partonów w tym samym zderzeniu protonów. W przybliżeniu „pocket formula” σ_DPS = σ_cc̄² / (2 σ_eff), a pary są nieskorelowane kinematycznie i często mają ten sam znak (LS: cc lub c̄c̄).

LHCb mierzył pary DD do przodu (2,0 < y < 4,5); udział LS 0,09–0,12, zob. `PD_pairs_diff/charm-PD-pairs-diff.md`. Pomiar w centralnym obszarze (|η| < 2,5), w zderzeniach minimum bias i z niskim progiem pT (3 GeV) uzupełnia ten obraz.

Mierzona wielkość to siła sygnału względem Pythii 8 A3:

μ = [N_par / N_D*]_dane / [N_par / N_D*]_MC

Normalizacja przez D* z tych samych danych skraca świetlność i większość wydajności rekonstrukcji. μ = 1 oznacza zgodność z Pythią A3, w której σ_eff ≈ 28 mb (notatki projektu, „Wyniki dotąd”).

## 2. Dane i symulacja

**Dane** (ROOT CharmAnalysis v2.0, h5 v1.1/v1.2; pamięć `charm-v21-production`):

| run | typ | trigger dominujący | uwagi |
|---|---|---|---|
| 428770 | physics_Main, pełna świetlność 2022 dla tej analizy | **HLT_noalg_L1J12** (dżet L1 > 12 GeV): 77 % N_D* | próbka nominalna (odsłonięcie 05.11). Klasa MBTS to w ~93 % HLT_j20, czyli też dżetowa (`PF2_lxplus/runs/lxplus/unblind/ubtrig/UBTRIG.md`) |
| 427929 | physics_MinBias | **HLT_mb_mbts_L1MBTS_1** (PS ≈ 1) | czysty minimum bias, uzysk D* 0,14–0,30 × 428770; klasa MB = GRL ∧ MBTS (`PF2_lxplus/runs/lxplus/unblind/r427929/UB427.md`) |

Rozkład klas triggera w 428770 (`ubtrig/UBTRIG.md`), wg udziału N_D* (r_c):

| klasa | r_c |
|---|---|
| J12only | 0,769 ± 0,019 |
| MBTS_j20 | 0,135 ± 0,015 |
| MBTS_noj20 | 0,010 ± 0,004 |
| other | 0,084 ± 0,012 |

**Symulacja (Pythia 8, tune A3):**
- **mb_nd**: inkluzywny non-diffractive. 80 plików, z tego 78 w ROOT v2.0, 19,2 M zdarzeń (`GEN_mcstat/runs/lxplus_v20/summary.md`). Służy do oczekiwania [pary/D*] i do map wag zdarzeń.
- **fd_nd**: wymuszone rozpady D. Daje szablony sygnału, wydajności kategorii Y i trening XGB.
- **Generator A3** (500 plików, prawda) daje normalizację z prawdy (`GEN_mcstat`) i topologie (`PREZ_2026-11-19/p6`).

**Wagi zdarzeń:** mapa 2D w2(n_trk, lead_pt) dane/MC. Czynniki c_w to ważone/nieważone [X/D*]_MC (`VER_budget/RAPORT.md`):

| D⁰ | D⁺ | D⁰D⁰ | D⁺D⁰ | D⁺D⁺ |
|---|---|---|---|---|
| 0,935 | 0,793 | 1,374 | 1,201 | 0,981 |

## 3. Selekcja pojedynczych D

**Kandydaci:**
- D⁰ → K⁻π⁺;
- D⁺ → K⁻π⁺π⁺, z wetem D*;
- D*⁺ → D⁰π⁺ przez Δm.

pT 3–30 GeV, |η| < 2,5.

**Definicje prawdy** (`PREZ_2026-11-19/p1/defs.md`):

| grupa | kody prawdy |
|---|---|
| sygnał | 0, 1 |
| refleksja K↔π | 100, 101 |
| inny charm / inny kanał | 2, 3, 6, 7 |
| częściowo zrekonstruowane | 8–11 |
| kombinatoryka | 4, 5, 12–15 |

**Klasyfikator XGB H:**
- 18 cech high-level: log pT, |η|, χ² wierzchołka, L_xy, istotność parametru zderzenia, pointing, |cos θ*|, d0/z0 córek, dE/dx i p kaonu, …
- Sygnał to fd_nd z prawdą, tło to pasma boczne danych, trening k-fold.
- Punkty pracy liczone per bin pT:
  - pojedyncze D: WP70;
  - nogi par: kategorie score;
  - D⁰ z tagiem D*: WP90.

**Arbitraż hipotez D⁰ (K↔π):**
- pary: hash;
- pojedyncze D⁰: keep_both.

Arbiter v2 na danych z tagiem D* daje SF = 0,977 ± 0,041 i po poprawce r(m) nie przewyższa hasha (notatki projektu, `HYP_d0`).

Liczby z wykresów `PREZ_2026-11-19/p1/` (`p1/cutflow.md`, `p1/numbers.json`):
- **S/B w oknie sygnału, mb_nd** (`p2_mass_*`, `cutflow.md`):

  | | bez cięcia score | WP70 |
  |---|---|---|
  | D⁰ | 0,0063 | 0,111 |
  | D⁺ | 0,0041 | 0,26 |

- **Odrzucenie tła przy ε_S = 0,70** (`p4_score_*`): D⁰ ε_B = 0,040, D⁺ ε_B = 0,0094.
- **Dane 428770** (`p5_data_mass_goodbad_*`): pik widoczny tylko w GOOD; stosunek okno/pasma 1,37 (D⁰) i 1,46 (D⁺).
- Porównanie drzewo vs sieć (low-level/high-level): `PREZ_2026-11-19/p3/` → MC (okno sygnału, OOF): zysk odrzucenia tła vs nominał BDT H przy ε_S 90 % — D⁰: MLP (te same 18 cech) 1,06, BDT na śladach 1,12, Deep Sets 1,16; D⁺: 1,09 / 1,57 / 1,80; przy tych samych wejściach sieć ≈ drzewo, zysk pochodzi z wejść niskopoziomowych (otoczenie śladów), które na danych tracą kalibrację (projekt B: 1,1–1,4× po wymuszeniu SF); rzeźbienie masy bez zmian (PREZ_2026-11-19/p3/runs/TABLE.md) (etap `preznet` w toku).

## 4. Pary i fit 2D

Kanały D⁰D⁰, D⁺D⁰, D⁺D⁺. Dla każdej pary fit 2D mas (m₁, m₂) w kategoriach score obu nóg (model planB, `charmtag.pairs`):
- **SS** (obie nogi prawdziwe): SS = n_pair · kk · Y_a,i · Y_b,j, przy kk = 1 nominalnie (`VER_budget/RAPORT.md`, `pairs/model.py:206`).
- **SB/BS**: jedna noga prawdziwa, druga z tła. SB = n_real · Y · NF.
- **BB**: tło × tło, zakotwiczone w SB×SB z danych przez GLM.

**Normalizacja (tryb direct):**

n_pair = N_D* · (n_pairs/n_D*)_MC · c_w,pair

- Pary MC pochodzą z poprawnych D w mb_nd bez cięcia score: 21 / 41 / 15 par, n_eff 8,2 / 13,5 / 7,5 (`VER_budget/RAPORT.md`).
- N_D* w danych = 4436, ±21,9 % (NDS po poprawce NLL; `NDS_ndstar/rerun_2026-11-12/rerun.md`).
- Fit był wykonany przy N_D* = 5584, a μ przeskalowano: × 5584/4436 = × 1,26 (`UNB_WYNIK.md` §1b, `ubtrig/UBTRIG.md`).

Weto śladów d0/z0 usuwa pary dzielące ślady. μ wspólne to minimum sumy profili kanałów; gładkie minimum to 1,19 [0,88; 1,51] przy bazie 5584 (`PD_pairs_diff/ndstar/ndstar.md`).

**Normalizacja z prawdy (nowa, `GEN_mcstat`):**
- [pary/D*] z generatora A3 (5·10⁸) z mapą ε direct.
- Akceptacja B1 v2.0 (wariant 78 plików, §3b): **PASS** (`GEN_mcstat/runs/lxplus_v20/acceptance.json`, 19.11).
- Domknięcie split-half truth/direct (`GEN_mcstat/runs/lxplus_v20/summary.md`):

  | D⁰D⁰ | D⁺D⁰ | D⁺D⁺ |
  |---|---|---|
  | 1,04 | 1,32 | 1,19 |

  Pulle 0,1 / 1,1 / 0,5.
- Przeliczenie μ i syst. normalizacji: 428770: μ = 1,27 [0,94; 1,62]; 427929: μ = 0,50 [0,17; 0,86] (normalizacja z prawdy **mb_nd**, decyzja z 19.11; A3 odrzucone — B7 FAIL, `VER_budget/RAPORT_gen_norm.md`; fit: `PF2_lxplus/runs/lxplus/unblind/truthnorm/TRUTHNORM.md`), syst ±0,61 → ±0,42 przy μ = 1,50 (±0,36 przy μ = 1,27): stat. MC par ±0,18 (n_eff ×5), model wag ±0,20 (dalej wg `GEN_mcstat/PROCEDURA.md`).

## 5. Systematyki — budżet wyniku głównego (427929)

Budżet μ(427929) = 0,50, normalizacja z prawdy mb_nd (`PF2_lxplus/runs/lxplus/unblind/r427929/BUDZET_427929.md`, `budzet_427929.py`; względne syst przeniesione do zapisu w `truthnorm/TRUTHNORM.md`). Δμ/μ wspólne (wagi kanałów 0,334 / 0,478 / 0,188), malejąco; Δμ bezwzględne przy μ = 0,50.

| źródło | Δμ/μ | Δμ | uwaga |
|---|---|---|---|
| **statystyka danych (profil × jackknife 1,13)** | **+0,83 / −0,73** | **+0,41 / −0,37** | dominuje |
| **N_D* (NDS 427929: stat 23,7 % ⊕ 14 wariantów 33,9 %)** | **±0,414** | **±0,21** | `r427929/nds/NDS_427929.md`; skala σ piku 1,40 × MC (`ndsw/NDSW.md`) |
| statystyka MC par (n_eff, baza nominału) | 0,316 | 0,16 | direct; propozycja z prawdy mb_nd: 0,162 — przypis ¹ |
| tune Pythii — strona górna (MPI ⊕ fragmentacja ALICE) | 0,409 | 0,20 | tylko +; k = 0,65 / 0,73 / 0,81 (`unblind/frag/FRAG.md`) |
| SF D⁰ (score) ±1σ | 0,086 | 0,04 |  |
| SF D⁺ (score) ±1σ | 0,085 | 0,04 |  |
| κ: plan B ↔ jedno κ na kanał | 0,064 | 0,03 |  |
| wydajności kategorii score ±5 % | 0,061 | 0,03 |  |
| tune Pythii — strona dolna (MPI ⊕ CR↓ ⊕ Monash) | 0,059 | 0,03 | tylko − |
| SF tagu D* ±14 % | 0,052 | 0,03 |  |
| przeniesienie SF 428770 → 427929 | 0,049 | 0,02 | HYBSF v3 |
| model wag (w2 ↔ w1) | 0,046 | 0,02 | w 427929 wagi ≈ 1; cały |c_w − 1| |
| κ_cls (korelacja ε nóg) | 0,045 | 0,02 | jednostronnie w dół |
| tło fałszywych nóg exp(poly2) → exp(poly3) | 0,043 | 0,02 |  |
| szerokość piku D⁰ ±6 % | 0,041 | 0,02 |  |
| pozostałe (swapy K↔π, szerokość/skala D⁺, refleksje, ε(n_trk), mapa ε, trigger MBTS, pile-up, dyfrakcja) | 0,040 | 0,02 | każdy ≤ 0,022 |
| *suma grup: kształty 0,070 · SF 0,141 · normalizacja 0,328 (dół) · κ/kategorie 0,089* | | | |
| **syst razem (bez N_D*)** | **+0,551 / −0,374** | **+0,28 / −0,19** | góra z fragmentacją ALICE (decyzja z 19.11) |

Wiersze kształtów, SF i κ przeniesione zachowawczo z 428770 (Asimov); w 427929 S/B par ≈ 4× lepsze, więc realnie mniejsze. Wiersz triggera nie występuje (MBTS, PS = 1, D* i pary ten sam trigger → 0,004).

¹ **Propozycja, czeka na decyzję (nominał §7 bez zmian):** wiersz „statystyka MC par” liczony nadal metodą bezpośrednią (n_eff bazy nominału, 0,316), choć μ pochodzi z normalizacji z prawdy mb_nd. Statystyka prawdy mb_nd (jackknife 26 grup ⊕ bootstrap map, `GEN_mcstat/out/truth_mbnd.json`, wagi w2 dziedziczone z bazy) daje 0,167 / 0,158 / 0,165 → wspólne **0,162** → syst **+0,24 / −0,13** (zamiast +0,28 / −0,19). Dolne oszacowanie przy wagach 427929 ≈ 1 (Poisson prawdy ⊕ bootstrap map): 0,083 → +0,23 / −0,11. Skrypt `r427929/mcstat_truth_prop.py`.

Budżet 428770 (stary §5) — Aneks A.

## 6. Kontrole

| kontrola | wynik | źródło |
|---|---|---|
| jackknife po zdarzeniach (200 grup), 428770 | σ_jk/σ = 1,18 ± 0,06 (kanały 1,12 / 1,10 / 1,21); Z z korelacjami ≈ 3,4σ | `PF2_lxplus/runs/lxplus/unblind/jk_ev200/ubjk.md` |
| toy: bias / pokrycie | bias 0,012 / 0,094 / −0,113; pokrycie 68 %: 0,80 / 0,66 / 0,64 | `UNB_WYNIK.md` §2 |
| bez kategorii z tagiem D* | 1,32 ± 0,76 / 0,94 ± 0,53 / 1,39 ± 0,93 | `UNB_WYNIK.md` §2 |
| więz κ partnera σ = 0,3 / 0,5 | μ = 1,00, Z = 3,42σ / 3,13σ (baza 5584) | `UNB_WYNIK.md` §1 |
| LS/OS (PD) | μ = 1,46 [1,08; 1,84]; f_LS = 0,33 ± 0,14, zgodne z Pythią | `PD_pairs_diff/charm-PD-pairs-diff.md` |
| emulacja triggera J12 (TRIGCAL) | κ = 1,146–1,388 → κ_emu = 1,24 ± 0,10 | `TRIG_A/SYST_TRIG.md` |
| **μ per klasa triggera, 428770 (UBTRIG)** | J12only 1,63 [1,25; 2,02]; MBTS_j20 −0,72 [−1,24; −0,18]; MBTS_noj20 1,49 [0,47; 3,32]; other −0,66 [−1,49; 0,14]; χ²/ndf = 15,66/3 (p = 0,001); J12only − reszta = +3,71σ (≈ +3,1σ z czynnikiem 1,18). Skala próbki 05.11, μ_nom = 1,19 | `PF2_lxplus/runs/lxplus/unblind/ubtrig/UBTRIG.md` |
| **czysty MB, 427929 (UB427)** | μ = 0,53 [0,17; 0,92] (Z = 1,57; normalizacja direct — wynik główny z prawdy mb_nd 0,50, §7) vs 428770: 1,50 [1,10; 1,90]; Δ = −0,97 (−1,60σ; −1,47σ z N_D*). σ 427929 bez czynnika korelacji zdarzeń | `PF2_lxplus/runs/lxplus/unblind/r427929/UB427.md` |
| jackknife 427929 | σ_jk/σ = 1,13 (kanały 1,28 / 1,16 / …); μ = 0,567 [0,19; 0,99] | etap `ub427jk` → `r427929/jk_ev200/` |
| **test zerowy UBTRIG** (klasy losowe tej samej liczności, 24 ziarna, lxplus) | klasa wielkości J12: 1,06 ± 0,18 (prawdziwe 1,63, Z = +3,25; 24/24 ziaren poniżej); J12 − reszta −0,43 ± 0,50σ (prawdziwe +3,71σ, Z ≈ +8); μ reszty w teście zerowym 1,69 ± 0,64 (prawdziwe −0,53) → estymator w małych klasach obciążony W GÓRĘ (+0,5 względem μ_nom 1,19), więc nie tłumaczy ujemnego μ reszty — efekt J12 nie jest artefaktem metody; małe klasy szumią (rozrzut 1,7–2,7). 6 ziaren lokalnie dawało Z = +5,3 przy zaniżonym rozrzucie — wartość z 24 ziaren zastępuje | `ubtrig/null/UBTRIG_NULL.md` |
| prerejestrowana miara triggera (syst_trig) | d = +1,04 ± 0,32 (Z = 3,24): d₁ (UB427) = +0,65 ± 0,47, d₂ (UBTRIG) = +1,39 ± 0,44; zgodność z emulacją +2,75σ. Status: „efekt” | `TRIG_A/syst_trig_out/` (19.11) |
| topologia \|Δη\| (λ = udział DPS ponad Pythię) | λ̂ = −1,11 [−1,50; −0,51]: pary bardziej kolinearne niż w Pythii; czysty DPS odrzucony 2,8σ, „nadwyżka = DPS” 2,1σ (`unblind/ubtopo/UBTOPO.md`) (prognoza σ_λ ≈ 0,85) | `Q_lxplus/ub_topo/PROGNOZA.md`, etap `ubtopo` |
| μ z siecią par (PAIRNET v3, pnetd) | 428770, te same pary: sieć par μ = 1,43 [1,04; 1,83] vs kategorie nominalne 0,78 [0,37; 1,19]; σ_μ 0,97 × (nie 0,79); Δμ = 0,65 ± 0,38 (1,7σ); lepszy opis D⁺D⁰ (dewiancja/ndf 1,01 vs 1,11) → kontrola, nie nominał. 427929: sieć par μ = 0,10 [−0,11; 0,40] (σ_68 0,53× kategorii nominalnych, bootstrap 0,24 vs 0,70), BDT par 0,10, kategorie nominalne bez kategorii tagu D* −0,25 [−0,71; 0,24], suma score nóg 0,72; nominał charmtag (z tagiem D*) 0,57; różnice względem nomcat w granicach bootstrapu (pnetd Δ = 0,34 ± 0,69); wszystkie schematy μ ≤ 1 — brak nadwyżki; rozrzut schematów ~0,8 przy małej statystyce; transfer tła D⁺D⁺ do 428770 słaby (χ² 37/6) — kontrola, nie nominał (PF2_lxplus/runs/lxplus/unblind/r427929/pnet/UBPNET427.md) (oczekiwane σ_μ ≈ 0,79 × nominał) | etap `ubpnet` → `unblind/ubpnet/UBPNET.md` |
| normalizacja z prawdy Pythia A3 z korektą decayType (kontrola, 20.11) | B7 po korekcie PASS (D*/zd. mb_nd/A3 = 0,99; było 1,36 — artefakt definicji: A3 ekskluzywne, mc23 topologiczne); D⁺/D* A3/mc23 ×1,09–1,11; f_ch A3 kor. 1,22 / 1,36 / 1,94 vs mb_nd 1,03 / 1,27 / 1,06 → pełny fit: μ(427929) = 0,375 [0,11; 0,67] (nominał mb_nd 0,50 bez zmian), 428770 = 1,07 (skala zapisu). Różnica ≈ 0,12 — rzędu wiersza tune, < stat. | `GEN_mcstat/A3_KOREKTA_2026-11-20.md` |

Zastrzeżenie do UBTRIG: μ < 0 w klasach bez J12 jest niefizyczne (oczekiwanie Pythii to 1). Test zerowy (24 ziarna) pokazał, że estymator w małych klasach jest obciążony W GÓRĘ (+0,5), więc obciążenie metody nie tłumaczy ujemnego μ reszty; ujemne wartości to szum małych klas (`ubtrig/null/UBTRIG_NULL.md`).

## 7. Wynik — wariant (B): wynik główny 427929 (decyzja z 19.11)

> **Status:** decyzja z 19.11 (~19:15): (1) wynik główny = 427929, 428770 = obserwacja obciążona triggerem J12; (2) wiersz tune Pythii poszerzony asymetrycznie w górę do δ ≈ −29 % par/D* (fragmentacja ALICE, `unblind/frag/FRAG.md`; zastępuje górną stronę CRqcd); (3) nominał μ = 0,50 (kryterium `ubndsw` nierozstrzygnięte), z uwagą o wysokim pT (≈ 0,72). Zapis w `UNB_WYNIK.md` (sekcja „Zapis od 19.11”).

### (B) 427929 — wynik główny (decyzja z 19.11)
Efekt triggera potwierdzony: syst_trig Z = 3,24 („efekt”), test zerowy 24 ziarna (klasa wielkości J12 losowo 1,06 ± 0,18 vs 1,63, Z = +3,25; J12 − reszta Z ≈ +8), topologia kolinearna (λ̂ = −1,11).

> W próbce czystego minimum bias (run 427929, HLT_mb_mbts_L1MBTS_1), z normalizacją z prawdy mb_nd: μ = **0,50 +0,41/−0,37 (stat, × 1,13 jackknife) +0,28/−0,19 (syst) ± 0,21 (N_D*)**, zgodne z Pythią A3 (−0,6σ przy skalach syst. w μ = 1, zachowawczo; −0,9σ przy skalach w μ̂). Syst asymetryczny: strona górna zawiera fragmentację ALICE (tune ↑ 41 % zamiast 12 %). Istotność względem braku par Z = 1,60 (≈ 1,4 z korelacjami zdarzeń).

- Przypis: wiersz statystyki MC liczony metodą direct (zachowawczo); z prawdy mb_nd syst byłby +0,24/−0,13 — **propozycja, czeka na decyzję** (§5, przypis ¹).
- Normalizacja direct (kontrola): μ = 0,53 +0,44/−0,40 +0,29/−0,20 ± 0,22 (`r427929/BUDZET_427929.md`).
- **Warianty interpretacji (nie nominał):** (i) szerokość piku D* z MC (wskazanie danych przy wysokim pT, `ndsw/NDSW.md`): μ ≈ 0,72; (ii) fragmentacja ALICE (μ/k, k = 0,71): μ ≈ 0,70; oba razem ≈ 1,01. Wszystkie warianty zgodne z Pythią — brak nadwyżki par w czystym minimum bias.
- N_D* = 2046 ± 485 (stat) ± 677 (syst, 14 wariantów NDS) → 41 % (`r427929/nds/NDS_427929.md`; 428770: 22 %).
- **Systematyka — szerokość piku D\*:** fit z wolną skalą σ daje 1,40 × MC (428770: 1,02, ten sam detektor). Przy skali stałej z MC N_D* = 1417 → μ ≈ 0,72 (direct: ≈ 0,77). Rozstrzyga etap `ubndsw` (szerokość piku D* w klasach MBTS 428770; kryterium zapisane przed wynikiem): wynik wg kryterium zapisanego przed wynikiem: NIEROZSTRZYGNIĘTE — skala σ w zdarzeniach MBTS 428770 = 0,96 ± 0,24 (zgodna z 1,02 i z 1,40 w 2σ; 1,40 odległe o 1,8σ); kontrole OK (428770 całość 1,02, 427929 1,41); pomocniczo: przy wysokim pT (8–30 GeV, małe tło) 427929 ma 1,12 [0,99; 1,26], MBTS 428770 0,87 [0,76; 1,01] (1,40 wykluczone 3,3σ) — szeroki pik tylko w binach z dużym tłem (3–8 GeV), co wskazuje na tło wciągnięte w sygnał, ale formalnie bez rozstrzygnięcia (PF2_lxplus/runs/lxplus/unblind/r427929/ndsw/NDSW.md). Wpływ na μ: syst N_D* pozostaje ±0,21–0,22 (obejmuje oba warianty); nominał 0,50 bez zmian, wariant skali z MC μ ≈ 0,72.
- Sieć par na 427929 (kontrola schematu): sieć par μ = 0,10 [−0,11; 0,40] (σ_68 0,53× kategorii nominalnych, bootstrap 0,24 vs 0,70), BDT par 0,10, kategorie nominalne bez kategorii tagu D* −0,25 [−0,71; 0,24], suma score nóg 0,72; nominał charmtag (z tagiem D*) 0,57; różnice względem nomcat w granicach bootstrapu (pnetd Δ = 0,34 ± 0,69); wszystkie schematy μ ≤ 1 — brak nadwyżki; rozrzut schematów ~0,8 przy małej statystyce; transfer tła D⁺D⁺ do 428770 słaby (χ² 37/6) — kontrola, nie nominał (PF2_lxplus/runs/lxplus/unblind/r427929/pnet/UBPNET427.md).

Źródła: `truthnorm/TRUTHNORM.md`, `r427929/BUDZET_427929.md`, `r427929/UB427.md`, `r427929/jk_ev200/`, `r427929/nds/`.

### (A) 428770 — obserwacja obciążona triggerem J12 (pomiar pomocniczy)

> μ(428770) = **1,27 +0,41/−0,40 (stat, × 1,18) ± 0,36 (syst, prawda mb_nd) ± 0,28 (N_D*)**; Z = 4,05 (≈ 3,4 z korelacjami zdarzeń); względem Pythii +0,4σ (bez wiersza triggera). Zapis bieżący w `UNB_WYNIK.md` (direct, 16.11): 1,50 +0,48/−0,46 ± 0,61 ± 0,33.

- Nadwyżka siedzi w zdarzeniach J12only: μ = 1,63 [1,25; 2,02] (prawda mb_nd: 1,37 [1,05; 1,71]); J12only − reszta +3,71σ; χ²/ndf klas 15,66/3.
- Miara triggera d = 1,04 ± 0,32 (Z = 3,24) → wg prerejestracji status „efekt”; wariant μ_corr = −0,05 ± 0,48. Emulacja J12 przewidywała tylko d_emu = 0,16 → dane 2,8σ ponad emulacją.
- **Klasa MBTS_noj20** (zdarzenia MBTS bez HLT_j20, r_c = 1 %): μ = 1,49 [0,47; 3,32] (prawda mb_nd 1,28 [0,40; 2,87]) — nominalnie nie-dżetowa, a wysoka. Test zerowy (24 ziarna): losowe klasy tej wielkości dają 0,74 ± 2,69 (13 określonych ziaren), więc to fluktuacja małej próbki, nie przeczy obrazowi dżetowemu — ale też go nie wspiera. Ujemne μ MBTS_j20 (−0,72) i other (−0,66) leżą nisko w szumie małych klas (test zerowy 1,92 ± 1,67 i 2,40 ± 2,44; obciążenie w górę).
- Topologia |Δη|: λ̂ = −1,11 [−1,50; −0,51] — pary bardziej kolinearne niż w Pythii, zgodne z rozszczepieniem gluonów w dżetach.

### (C) Wspólny fit próbek bez triggera dżetowego — nie wykonany
427929 (MB) ⊕ 428770 bez klas J12only i HLT_j20. Analiza wykonalności (`ub_comb`, 19.11): czysty podzbiór 428770 to 0,84 % zdarzeń (N_D* ≈ 19); kombinacja z klasami nie-J12 daje ≈ 0,40, czyli praktycznie sam 427929. Bieg nie wysłany — (C) ≈ (B).

## 8. Interpretacja (σ_eff)

Przy μ zgodnym z Pythią A3 dane są zgodne z σ_eff ≈ 28 mb (notatki projektu). Nadwyżka μ ≈ 1,5 odpowiadałaby σ_eff ≈ 6–16 mb (notatki projektu, „Interpretacja”).

Za zgodnością z Pythią przemawiają też:
- udział LS w danych 428770 f_LS = 0,33 ± 0,14 (kanały 0,22 / 0,34 / 0,49, χ² = 0,38/2), Pythia po poprawce znaku D⁰ 0,27 ± 0,03 (+0,4σ); f_LS nie wchodzi do fitu μ (`unblind/fls/FLS.md`);
- topologia |Δη| λ̂ = −1,11 [−1,50; −0,51]: pary bardziej kolinearne niż w Pythii; czysty DPS odrzucony 2,8σ, „nadwyżka = DPS” 2,1σ (`unblind/ubtopo/UBTOPO.md`).

**Wnioski zależne od wariantu:**
- **(A):** nie interpretujemy w σ_eff — nadwyżka 428770 przypisana selekcji triggerem J12, nie DPS.
- **(B):** dolna granica **σ_eff > ~6 mb (95 % CL)** przy f_DPS = 0,40 (syst końcowy z asymetrią fragmentacji; rozpiętość wg metody toya i f_DPS 4–11 mb; jawny toy 9,4 mb), przy założeniu, że całe odchylenie μ od Pythii to DPS (μ = (1 − f) + f·28 mb/σ_eff). Warianty interpretacji (pik D* z MC μ ≈ 0,72; fragmentacja ALICE μ ≈ 0,70) > ~5 mb. Granica słaba: wyklucza tylko σ_eff ≲ 5–6 mb (`unblind/r427929/SIGMA_EFF.md`, sekcja „Aktualizacja 19.11 wieczór”).
  - Założenie fragmentacji: granica przyjmuje udział mezonów D na kwark c jak w Pythii. W pp obserwowany jest nadmiar barionów charmowych względem e⁺e⁻ (ALICE: Λc 20 % vs 6 %, D⁰ 40 % vs 61 %); jeśli Pythia przecenia f(c→D), stosunek par/D* w danych spada (∝ f) i μ < 1 jest oczekiwane bez DPS — część deficytu μ może być fragmentacją; od 19.11 górna strona wiersza tune obejmuje fragmentację ALICE (δ ≈ −29 % par/D*, μ +41 % względnie; `unblind/frag/FRAG.md`), więc efekt jest w syst, a μ/k = 0,70 podajemy jako wariant interpretacji.

Frakcje fragmentacji Pythii (D⁺/D* 1,33 wobec 0,99 w e⁺e⁻) należą do odniesienia. Przy frakcjach e⁺e⁻ μ_wsp. ≈ 1,9 (`VER_budget/RAPORT.md` (C)); to informacja, nie systematyka.

### Co odrzucono i dlaczego
Pełna lista: `~/Projects/charm/charm-a-nieudane/00_INDEX.md`.

- Pary z tagiem AFP: NO-GO.
- Sieć zamiast drzewa dla pojedynczego D⁰ (nnh, 17.11): odrzucona, brak zysku.
- Sieć par z kinematyką pary (pnet_kin): odrzucona.
- Noga D⁺ z cechami Dalitza: zysk zawyżony przez migrację.
- Sieć par na tle z MC redukuje tło w danych o połowę słabiej niż w MC. Naprawione częściowo treningiem na pasmach danych (v3).
- Transfer D⁺ przez Pythię jako osobny wiersz syst.: odrzucony, pokryty przez N_D* i model wag (`VER_budget/RAPORT.md` (C)).
- Topologia par per klasa triggera i w 427929: odrzucona przed liczeniem, brak mocy (σ_λ ≈ 1,9).
- Interpretacja μ = 1,50 jako fizyki double charm w 428770: podważona testem triggera (19.11), potwierdzone testem zerowym.
- Sieć par (PAIRNET v3) jako schemat główny: na danych σ_μ 0,97 × zamiast 0,79 ×; zostaje kontrolą (`unblind/ubpnet/UBPNET.md`).
- Normalizacja z prawdy A3: B7 FAIL (D⁺/D* × 1,42 vs mc23) → prawda mb_nd. *Uwaga 20.11: część tej różnicy to definicja `decayType` (topologiczny; D*/zd. ×1,36 = 1/0,734, D⁺/D* po korekcie ≈ ×1,11) — nominał mb_nd nietknięty (ta sama definicja w liczniku i mianowniku), ocena B7 dla A3 do powtórzenia (`GEN_mcstat/PROCEDURA.md`, uwaga 20.11).*
- Wariant (C): brak czystego podzbioru w 428770 (N_D* ≈ 19).
- Błędy techniczne z 18–19.11 (ubtrig NaN, dopasowanie tagów gen_norm, fity MULT): naprawione, bez wpływu na μ.

## 9. Stan punktów otwartych (19.11, wieczór)

Wypełnione 19.11: test zerowy (6 ziaren lokalnie), jackknife 427929, normalizacja z prawdy mb_nd, wiersz triggera, syst. i N_D* 427929 (ub427nds), sieć par 428770, topologia, f_LS.

1. **Decyzja (19.11, ZROBIONE):** wariant (B) wynikiem głównym, syst fragmentacji, nominał 0,50; zapis w `UNB_WYNIK.md` (sekcja „Zapis od 19.11”).
2. **Szerokość piku D\* (`ubndsw`)** → wynik wg kryterium zapisanego przed wynikiem: NIEROZSTRZYGNIĘTE — skala σ w zdarzeniach MBTS 428770 = 0,96 ± 0,24 (zgodna z 1,02 i z 1,40 w 2σ; 1,40 odległe o 1,8σ); kontrole OK (428770 całość 1,02, 427929 1,41); pomocniczo: przy wysokim pT (8–30 GeV, małe tło) 427929 ma 1,12 [0,99; 1,26], MBTS 428770 0,87 [0,76; 1,01] (1,40 wykluczone 3,3σ) — szeroki pik tylko w binach z dużym tłem (3–8 GeV), co wskazuje na tło wciągnięte w sygnał, ale formalnie bez rozstrzygnięcia (PF2_lxplus/runs/lxplus/unblind/r427929/ndsw/NDSW.md), syst N_D* pozostaje ±0,21–0,22 (obejmuje oba warianty); nominał 0,50 bez zmian, wariant skali z MC μ ≈ 0,72. Może przesunąć μ(427929) z 0,50 do ≈ 0,72.
3. **Test zerowy, 24 ziarna (`ubtrignull`)** → J12: 1,06 ± 0,18 vs 1,63 (Z = +3,25); J12 − reszta Z ≈ +8; małe klasy obciążone w górę (+0,5), ujemnego μ reszty nie tłumaczą (potwierdzenie 6 ziaren lokalnych).
4. **Sieć par na 427929 (`ubpnet427`)** → sieć par μ = 0,10 [−0,11; 0,40] (σ_68 0,53× kategorii nominalnych, bootstrap 0,24 vs 0,70), BDT par 0,10, kategorie nominalne bez kategorii tagu D* −0,25 [−0,71; 0,24], suma score nóg 0,72; nominał charmtag (z tagiem D*) 0,57; różnice względem nomcat w granicach bootstrapu (pnetd Δ = 0,34 ± 0,69); wszystkie schematy μ ≤ 1 — brak nadwyżki; rozrzut schematów ~0,8 przy małej statystyce; transfer tła D⁺D⁺ do 428770 słaby (χ² 37/6) — kontrola, nie nominał (PF2_lxplus/runs/lxplus/unblind/r427929/pnet/UBPNET427.md).
5. **Prezentacja, pkt 3 (`preznet`)** → MC (okno sygnału, OOF): zysk odrzucenia tła vs nominał BDT H przy ε_S 90 % — D⁰: MLP (te same 18 cech) 1,06, BDT na śladach 1,12, Deep Sets 1,16; D⁺: 1,09 / 1,57 / 1,80; przy tych samych wejściach sieć ≈ drzewo, zysk pochodzi z wejść niskopoziomowych (otoczenie śladów), które na danych tracą kalibrację (projekt B: 1,1–1,4× po wymuszeniu SF); rzeźbienie masy bez zmian (PREZ_2026-11-19/p3/runs/TABLE.md).
6. **Interpretacja** → σ_eff > ~6 mb (95 % CL, f_DPS = 0,40; syst końcowy; rozpiętość 4–11 mb), warianty (pik z MC, fragmentacja ALICE) > ~5 mb (`r427929/SIGMA_EFF.md`, aktualizacja 19.11 wieczór).

## Aneks A. Budżet 428770 (dawny §5) — obserwacja obciążona triggerem J12

Pomiar pomocniczy, obserwacja obciążona triggerem J12 (§7 (A)); nie dotyczy wyniku głównego. Budżet przejściowy z 16.11, wspólne μ = 1,50 (`VER_budget/RAPORT.md`, `UNB_WYNIK.md` aktualizacja 16.11):

| źródło | Δμ (wspólne) | źródło liczby | uwaga |
|---|---|---|---|
| statystyka MC par (n_eff) | ±0,46 | VER_budget (A) | ma zastąpić go normalizacja z prawdy |
| model wag zdarzeń (w2 ↔ w1) | ±0,24 | VER_budget (A) | proxy z pojedynczych D |
| tune Pythii (MPI ⊕ CR ⊕ Monash) | ±0,17 | `TUNE_pythia/budget_with_tune.md` | asymetria +0,43/−0,40 łącznie |
| κ_cls (korelacja ε między nogami) | −0,06 (jednostronnie) | VER_budget (B) | MC 1,03–1,06, mechanizm niewyjaśniony |
| kształty, SF, kategorie, κ partnera | w sumie U | `UNB_WYNIK.md` §4 | względnie D⁰D⁰ 0,325 / D⁺D⁰ 0,279 / D⁺D⁺ 0,366 |
| **suma syst. (bez N_D*)** | **±0,61** | VER_budget | |
| N_D* (NDS) | ±0,33 | `UNB_WYNIK.md` | 21,9 % |
| **normalizacja z prawdy (zastępuje mc_stat + model wag)** | **±0,27 (stat. MC par ±0,18 ⊕ model wag ±0,20; razem syst ±0,42 przy μ = 1,50)** | `GEN_mcstat` | oczekiwane ~±0,15–0,27 (`UNB_WYNIK.md`, 16.11) |
| **trigger / selekcja zdarzeń** | **nie jako wiersz syst.: status „efekt” (d = 1,04 ± 0,32, Z = 3,24 > 3) → 428770 traktowany jako obciążony triggerem; wariant μ_corr = −0,05 ± 0,48** | `TRIG_A/syst_trig_out/` | reguły z prerejestracji (`TRIG_A/SYST_TRIG.md`, 19.11) |

Reguła dla wiersza triggera (`TRIG_A/SYST_TRIG.md`):
- |Z_trig| < 2 i zgodność z emulacją → −μ·d_emu ≈ −0,23 (κ_emu = 1,24 ± 0,10, d_emu = 0,156).
- 2 ≤ Z < 3 → wiersz łączny „selekcja zdarzeń”, który zastępuje model wag.
- Z ≥ 3 → dodatkowo wariant μ_corr = μ(1 − d). Zmiana nominału tylko osobną decyzją.
