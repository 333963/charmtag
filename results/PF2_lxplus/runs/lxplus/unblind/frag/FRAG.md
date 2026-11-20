# FRAG — frakcje fragmentacji c→hadron a μ (19.11.2026, projekt A)

Pytanie: μ = (pary/D*)_dane / (pary/D*)_Pythia. Jeśli Pythia (fragmentacja strojona do e⁺e⁻) przecenia udział mezonów D na kwark c
(w pp nadmiar barionów charmowych), to (pary/D*) w danych jest niższe bez żadnej zmiany DPS. Skrypt: `frag.py`.

## Wejścia
| | f(D⁰) | f(D⁺) | f(D*⁺) | źródło |
|---|---|---|---|---|
| ALICE pp 5,02 TeV, |y|<0,5, pT-całkowane | 39,1 ± 1,7 +2,5/−3,7 % | 17,3 ± 1,8 +1,7/−2,1 % | 15,5 ± 1,2 +4,1/−1,9 % | arXiv 2105.06335 Tab. 1 (Λc 20,4 %, Ξc⁰ 8,0 %) |
| Pythia (= e⁺e⁻, Lisovyi 2016) | 60,9 % | 24,0 % | 24,3 % | sprawdzone na prawdzie mc23 mb_nd: f(D⁰)/f(D*) = 2,55 (e⁺e⁻ 2,51), f(D⁺)/f(D*) = 1,03 (0,99) — z TRUTH_NPAIR (na D*, po BR; bez poprawki akceptancji) |

Absolutnej frakcji bariony w mb_nd nie liczyłem (brak truth Λc w h5) — przyjęte jak e⁺e⁻; stosunki mezonów zgadzają się z e⁺e⁻ do 2–4 %.

## Czynnik k = (f_a f_b / f_D*)_ALICE / (…)_Pythia (c i c̄ fragmentują niezależnie; D⁰ inkluzywnie z D*)
| kanał | k (centralne) | k (toy, błędy ALICE niezależne — zachowawczo) |
|---|---|---|
| D⁰D⁰ | 0,647 | 0,61 [0,46; 0,79] |
| D⁺D⁰ | 0,725 | 0,68 [0,51; 0,88] |
| D⁺D⁺ | 0,812 | 0,77 [0,52; 1,08] |
| wspólne (wagi 0,334/0,478/0,188, 1/k ważone) | **0,71** | 0,66 [0,50; 0,85] |

Prawie cały efekt to absolutne f(D⁰) (bariony + Ds: ALICE ~36 % c w barionach vs ~8 % w Pythii); stosunki D*/D⁰, D⁺/D⁰ są bliskie Pythii.

## μ po przeważeniu Pythii do frakcji ALICE (μ_frag = μ / k)
| wynik | μ | μ_frag (ALICE, k = 0,71) | μ (Pythia CRqcd, δ = −11 %, k = 0,89) |
|---|---|---|---|
| 427929 (prawda mb_nd) | 0,50 | **0,70** | 0,56 |
| 427929 (direct) | 0,53 | 0,75 | 0,60 |
| 427929, szerokość piku z MC | 0,72 | 1,01 | 0,81 |
| 428770 (prawda mb_nd) | 1,27 | 1,79 | 1,43 |
| 428770 (zapis 16.11) | 1,50 | 2,11 | 1,69 |

σ_eff (95 % CL, f_DPS 0,40, model z SIGMA_EFF.md): μ_up 2,12 → 2,99 po przeważeniu → σ_eff > ~4,7 mb (zamiast 7,4 mb).

## Podwójne liczenie z wierszem tune CR
CRqcd = `ColourReconnection:mode = 1; BeamRemnants:remnantMode = 1` (Q_lxplus/syst/variants.sh) — CR „QCD-based”, który w Pythii
daje wzmocnienie barionów (ALICE: CR Mode 2 opisuje Λc/D⁰). Wiersz tune_CR (δ[pary/D*]_reco = −11 %, μ +0,12/−0,04 względnie) mierzy
więc TEN SAM mechanizm w naszej akceptancji, ale z domyślnymi parametrami trybu 1 — słabiej niż pomiar ALICE (−11 % vs −29 %).
**Nie dodawać obu.** Fragmentacja ALICE zastępuje (rozszerza) górną stronę wiersza tune_CR.

## Zastrzeżenia
- ALICE: 5,02 TeV, |y|<0,5, pT-całkowane; my: 13,6 TeV, |η|<2,5, pT > 3 GeV. Λc/D⁰ maleje z pT (ok. 0,5 przy 2–4 GeV → ok. 0,3 przy 8–12 GeV)
  i jest niższe w przód (LHCb) → efekt w naszej akceptancji raczej mniejszy niż 29 %; realny zakres k ≈ 0,71 (ALICE) … 0,89 (CRqcd).
- Błędy ALICE traktowane jako niezależne (w rzeczywistości frakcje sumują się do 1 i dzielą BR/lumi → korelacje) — przedział zachowawczy.
- Frakcje „Pythia” = e⁺e⁻ zweryfikowane tylko stosunkami mezonów; absolutny udział barionów w mc23 mb_nd niesprawdzony.

## Rekomendacja
1. **Nie poprawiać nominału** (μ jest zdefiniowane względem Pythii; zmiana odniesienia = osobna decyzja).
2. **Syst:** wiersz tune_CR rozszerzyć asymetrycznie w górę do δ ≈ −29 % (μ +0,41 względnie, tj. +0,20 przy μ = 0,50), w dół bez zmian; nie dodawać osobnego wiersza.
3. **Interpretacja (nota):** podać μ_frag obok μ i opisać, że deficyt μ < 1 w 427929 jest w pełni zgodny z fragmentacją z pp (μ_frag ≈ 0,7–1,0), a nadwyżka w 428770 po przeważeniu rośnie (1,8–2,1), co wzmacnia wniosek o obciążeniu triggerem.

Proponowany tekst do noty (§8, po założeniu fragmentacji):
> Przeważenie Pythii do frakcji fragmentacji zmierzonych w pp przez ALICE (f(D⁰) 39 %, f(D⁺) 17 %, f(D*⁺) 16 %, Λc 20 %) obniża oczekiwane pary/D*
> o czynnik k = 0,71 (0,65 / 0,73 / 0,81 dla D⁰D⁰ / D⁺D⁰ / D⁺D⁺; przedział 0,50–0,85 z błędów ALICE). Wtedy μ(427929) = 0,70 (1,01 przy szerokości
> piku z MC) — zgodne z brakiem dodatkowego DPS — a μ(428770) = 1,8. Wariant CRqcd w naszej akceptancji daje słabszy efekt (k = 0,89), więc
> realny czynnik leży między 0,71 a 0,89. Granica σ_eff spada do ok. 4,7 mb przy frakcjach ALICE.
