# VER_budget: niezależna weryfikacja trzech kandydatów na braki w budżecie μ (16.11.2026)

Weryfikacja lokalna, tylko do odczytu; wyniki nominalne bez zmian. Nic nie poszło na lxplus.
Nowe obliczenie jest jedno: `kcls_check.py` → `kcls_check.json` (κ_cls odtworzone, test mechanizmu pT, bootstrap po zdarzeniach).
Konwencja „wspólne” jest taka jak w `charmtag/syst/report.py`: Σ w·δ, w ∝ 1/σ²_stat(Asimov ślepy) = 0,298 / 0,580 / 0,122. Skala × 1,495.

Nominał (bez zmian): **μ = 1,50 +0,48/−0,46 (stat) ± 0,42 (syst) ± 0,33 (N_D*)**, syst wspólne względne 0,284, w tym „statystyka MC par” 0,188.

## Jak oczekiwanie wchodzi do fitu (kod, sprawdzone)
- `pairs/model.py:206`: `SS = n_pair[ch] · kk · Y[(a,i)] · Y[(b,j)]`, przy czym kk = 1 nominalnie.
- `build.py` `cmd_pairs` (tryb `direct`): `n_pair = N_D* · n_pairs/n_dstar · c_w_pair`. Pary są liczone z `legs_for_truth`, czyli z poprawnych D w mb_nd **bez cięcia score**.
  Czynnik `c_w_pair` (`weight_factors`) to stosunek ważony/nieważony z `mb_rates_2d.json` (wagi 2D w2) na **tych samych** 21/41/15 parach.
- `Y` = N(kat)/N(poprawnych) z pojedynczych D w signal_mc (fd_nd) × SF (`signal_templates`, `sf_cat`).
- `syst/norm.py`: źródło `mc_stat` = `rel_stat` = 1/√n_pairs (nieważone) ⊕ 1/√n_D*. Brak źródła „model wag zdarzeń”.
  `eff_ntrk` to wariant mapy ε w trybie truth, a nie wag. Tune to warianty generatora bez zmiany wag.

## (A) Model mieszanki zdarzeń (w1↔w2), deklarowane ~±0,22: **częściowo potwierdzony**
- **Brak wiersza się potwierdza.** Żaden wiersz budżetu (tune, eff_ntrk, dyfrakcja, mc_stat) nie zmienia modelu wag. Tune nie zmienia też D⁺/D* (niżej, C).
- **Rozmiar ±0,22 jest w większości szumem.** Wynika z różnic pair/D* równych 1,138 / 0,922 / 0,474 na tych samych 21/41/15 parach.
  Przy w1 n_eff wynosi zaledwie 4,0 / 18,9 / 4,7, a przy w2 8,2 / 13,5 / 7,5 (`kontrole_t1_t3_t4.../t4_rates.json`).
  Oczekiwany szum ln(w2/w1) to ~0,3–0,6 dla D⁰D⁰ i D⁺D⁺ oraz ~0,25 dla D⁺D⁰. Zmierzone przesunięcia mieszczą się w ≲0,5σ (D⁰D⁰, D⁺D⁰) i ~1,5–2σ (D⁺D⁺).
- **Jest silnie skorelowany z wierszem MC par**, bo ciężkie wagi tych samych kilku zdarzeń dają oba efekty. Dodanie go w kwadraturze liczy szum dwa razy.
- **Ważniejsza obserwacja: wiersz „statystyka MC par” jest dziś zaniżony.** W trybie direct n_pair ∝ Σw par, więc właściwa niepewność to 1/√n_eff = 0,35 / 0,27 / 0,37, a nie 1/√n = 0,215 / 0,159 / 0,260.
  Te same 35/27/37 % podawał `fizyka_blok1` (09-28), a `norm.py` używa ich tylko w trybie frozen. Wspólne: 0,188 → 0,308.
- **Część fizyczna jest realna i da się ją oszacować bez szumu par.** Na pojedynczych D (9,4k / 8,7k / 1,3k) w2/w1 daje D⁰/D* 0,974 i D⁺/D* 0,833.
  Proxy faktoryzowane dla par: δ = 0,052 / 0,189 / 0,307, wspólne 0,163.
  Jest to oszacowanie zachowawcze z dwóch powodów: w1 gorzej domyka D* (R = 3,06 wobec 1,70), a proxy pomija zależność krotności par.
- **Rekomendacja wg reguły:** bez poprawki centralnej (w2 jest modelem preferowanym przez dane, alternatywa nie jest istotna).
  Zamiast „±0,22 obok 1/√n” proponuję dwie zmiany:
  1. mc_stat → 1/√n_eff (0,308 wspólne);
  2. nowy wiersz „model wag (w2↔w1), proxy pojedynczych D”, 0,163 wspólne, w pełni skorelowany między kanałami i niezależny od tune i N_D*.
  Uwaga: konwencja liniowa traktuje statystykę MC jako skorelowaną między kanałami, choć pary są rozłączne. Przy kanałach niezależnych n_eff daje 0,195 zamiast 0,308.

## (B) κ_cls (korelacja ε klasyfikatora między nogami): **potwierdzony jako efekt w MC, mechanizm niezrozumiany**
- **Nie ma podwójnego liczenia.** n_pair jest liczone przed klasyfikatorem, a SS = n_pair·Y·Y przy kk = 1, więc korelacja nie jest wbudowana. Parametry κ planu B mnożą SB/BS, nie SS.
- **Odtworzyłem κ_cls przy WP65** (fd_nd, v11 H): 1,061 / 1,029 / 1,058, z bootstrapem po zdarzeniach ±0,010 / 0,009 / 0,021. Istotność 6σ / 3σ / 2,8σ statystycznie.
- **Test mechanizmu.** Przewidywanie z ε(pT) faktoryzowanego (20 binów pT, pary fd_nd) daje **1,007 / 1,002 / 0,997**.
  Twardsze widmo nóg par nie tłumaczy więc efektu. Zostaje korelacja na poziomie zdarzenia (zmienne zdarzenia lub otoczenia w klasyfikatorze), niezbadana.
- **Kategorie z tagiem D* (D⁰ do WP90)** mają κ = 1,018 / 1,008 / 1,010.
- **Przeniesienie z fd_nd na skim danych** (twardsze zdarzenia) nie jest sprawdzone. SF są mierzone na pojedynczych nogach, więc tego nie łapią.
- **Rekomendacja:** warunek (iii) nie jest spełniony, więc bez poprawki centralnej. Proponuję nowy wiersz ±|κ−1| (pełny, nie 50 %): 0,061 / 0,029 / 0,058, wspólne 0,042 → **±0,06**.
  Kierunek: jeśli efekt się przenosi, μ maleje. Wiersz jest niezależny od pozostałych; częściowe pokrycie z `cat_Y` (±5 % niezależnie per kategoria) można pominąć.

## (C) Transfer D⁺ przez Pythię (DR⁺ = 0,76), deklarowane ~±20 %: **odrzucony jako osobny wiersz, mechanizm przejmuje (A)**
- **Część wspólna** (~0,8–0,9) to N_D*. Przyznaje to sam raport FR i jest ona już w ±0,33.
- **Część zależna od gatunku**, DR⁺/DR⁰ = 0,89 ± 0,21, różni się od 1 o 0,5σ, więc warunek (ii) nie jest spełniony.
- **FR porównuje z MC nieważonym, a nominał używa MC ważonego.** c_w_single (`pair_inputs.meta`) wynosi D⁺ 0,793 i D⁰ 0,935.
  Po ich zastosowaniu (całkowanie po pT) DR⁺ ≈ 0,96, DR⁰ ≈ 0,95, a DR⁺/DR⁰ ≈ 1,01. Domniemany deficyt D⁺ ma rozmiar korekty wag, którą nominał już stosuje.
  Per bin pT (8–30 GeV) nie jest to sprawdzone, bo lokalnie brak wag per kandydat.
- **Tune nie pokrywa frakcji fragmentacji.** Sprawdziłem w `runs/lxplus/truth/combined_*`: D⁺/D* wynosi 1,321–1,354 we wszystkich wariantach (A3 1,326), a D⁰/D* 2,52–2,56. ±0,17 nie obejmuje więc tej kwestii.
  Główną zidentyfikowaną zależność MC dla D⁺/D* (−17 % w2↔w1) obejmuje natomiast nowy wiersz modelu wag z (A). Osobny wiersz C liczyłby ją drugi raz.
- **SF:** test FR stosuje SF per gatunek (D⁺ 1,01 przy WP70, D⁰ 0,93 przy WP80), zbliżone do nóg par (WP65 / WP90). Niepewność SF (9 %) jest już w wierszu SF budżetu, więc nie wolno jej dodawać ponownie przez ±0,21.
- **Domknięcie mapy fd→mb** (D⁺ 0,949 wobec D⁰ i D* 0,72, czyli ×1,32) liczbowo równa się 1/DR⁺. Nominał i FR biorą jednak D⁺/D* z tej samej rekonstrukcji mb, więc to nie tłumaczy DR⁺ wobec nominału.
  Gdyby to był artefakt rekonstrukcji lub dopasowania prawdy w mb, oznaczałoby to bias. Warto to zbadać, ale dziś nie jest to wiersz syst.
- **Interpretacja:** μ jest zdefiniowane jako dane/Pythia A3, więc frakcje fragmentacji Pythii (D⁺/D* 1,33 wobec 0,99 w e⁺e⁻) należą do odniesienia.
  Przy frakcjach e⁺e⁻ μ_wsp. wyniosłoby ≈ 1,9. Można to podać informacyjnie, nie jako syst.
- **Rekomendacja:** brak zmian centralnych i brak osobnego wiersza. Kontrola rozstrzygająca to DR⁺, DR⁰ w 8–30 GeV z MC ważonym tymi samymi wagami 2D per bin pT.

## Zapis μ po rekomendacjach (propozycja, nominał niezmieniony)
Syst. wspólne bez N_D*: 0,284² − 0,188² (stara mc_stat) + 0,308² (n_eff) + 0,163² (model wag) + 0,042² (κ_cls) → **0,411 → ±0,61**.

> **μ = 1,50 +0,48/−0,46 (stat) ± 0,61 (syst) ± 0,33 (N_D*)**, z rozbiciem: stat. MC par (n_eff) ±0,46, model wag ±0,24, κ_cls ±0,06 (jednostronny w dół), tune ±0,17. Z bez zmian.

Warianty dla porównania:

| wariant | syst |
|---|---|
| tylko mc_stat → n_eff (bez wiersza modelu) | ±0,56 |
| kanały statystyki MC niezależne | ±0,50 |
| propozycja SN (±0,22 ⊕ 50 % κ obok 1/√n), podwójne liczenie szumu | ±0,48 |

Dominująca pozycja budżetu to dziś statystyka MC par. Wiersze (A), oba mc_stat i model wag, zmniejszyłaby prawda 5·10⁸ ważona wagami zdarzeń (SN §6 p. 4).
