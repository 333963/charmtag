# UBPNET — μ par D–D w 428770 z kategoriami sieci par (PAIRNET v3) — odsłonięte (zatwierdzone 18.11)

Ten sam model fitu 2D (sygnał MC × μ + 3 składowe tła z wolnymi normalizacjami i kształtami exp(poly2) per kategoria) i te same pary danych (h5 v2.0, zdarzenia próbki nominalnej v1.1) dla wszystkich schematów; nomcat_* = nominalne kategorie nóg (bez kategorii z tagiem D*), pnetd = sieć par z tłem z SB×SB danych (schemat główny), xgbd / base_sum = kontrola. Normalizacja n_pair z pair_inputs odsłonięcia; schemat sieci × SF_rel (pary hybrydowe: D⁺ HYBSF v3, D⁰ HYBDS v3).

**Kontrola nominału (charmtag, pairs_full.npz):** D0D0 1.1583 (ref 1.1583), DpD0 1.1648 (ref 1.1648), DpDp 1.3856 (ref 1.3856); wspólne μ_nom = 1.1889; zgodność: TAK. Skala do zapisu bieżącego ×1.262.

Dane: nogi 3406212 → w zdarzeniach próbki nominalnej 3406212 (zdarzeń z nogami 970501 z 970501); pary w zakresie mas {'D0D0': 558687, 'DpD0': 1717434, 'DpDp': 1557887}; pokrycie zdarzeń nóg nominału przez nogi v2.0: 0.740.

| schemat | μ D0D0 | μ DpD0 | μ DpDp | μ wspólne [68 %] | Z | σ_68 / σ(nomcat_mc) | σ Asimov / σ Asimov(nomcat_mc) | σ bootstrap | μ × skala (zapis bieżący) |
|---|---|---|---|---|---|---|---|---|---|
| nomcat_fact | 0.06 ± 0.75 | 0.67 ± 0.60 | 1.99 ± 0.92 | **0.797** [0.389; 1.207] | 1.98 | 1.001 | 1.001 | 0.360 | 1.005 |
| nomcat_mc | 0.06 ± 0.73 | 0.66 ± 0.60 | 1.99 ± 0.92 | **0.778** [0.373; 1.191] | 1.94 | 1.000 | 1.000 | 0.358 | 0.982 |
| pnetd | 0.57 ± 0.79 | 1.12 ± 0.56 | 3.14 ± 0.93 | **1.432** [1.037; 1.833] | 3.71 | 0.973 | 0.969 | 0.347 | 1.807 |
| xgbd | 0.36 ± 0.76 | 0.94 ± 0.57 | 2.88 ± 0.91 | **1.236** [0.836; 1.643] | 3.15 | 0.987 | 0.950 | 0.345 | 1.560 |
| base_sum | -0.32 ± 0.62 | 0.98 ± 0.55 | 1.93 ± 0.90 | **0.781** [0.393; 1.179] | 2.04 | 0.962 | 0.938 | 0.366 | 0.986 |

## Różnica względem nomcat_mc (te same pary, ten sam model) i systematyki schematu

| schemat | Δμ = μ − μ(nomcat_mc) | σ_Δ (bootstrap) | Δ/σ_Δ | korelacja (bootstrap) | syst SF norm. (δμ) | syst nachylenie SF (δμ) | syst tło poly3 (δμ) | Asimov: μ0 → μ̂ |
|---|---|---|---|---|---|---|---|---|
| nomcat_fact | 0.018 | 0.012 | 1.58 | 1.00 | – | – | 0.002 | 0.797 → 0.796 |
| pnetd | 0.654 | 0.378 | 1.73 | 0.43 | 0.122 | -0.016 | -0.017 | 1.432 → 1.432 |
| xgbd | 0.458 | 0.322 | 1.42 | 0.58 | 0.102 | -0.014 | 0.031 | 1.236 → 1.235 |
| base_sum | 0.003 | 0.265 | 0.01 | 0.73 | 0.063 | -0.042 | -0.017 | 0.781 → 0.781 |

## SF_rel schematów (SF schematu / SF nomcat; para = iloczyn nóg)

| kanał | schemat | SF_rel | ± | nachylenie góra / dół | ρ(schemat, nomcat) MC |
|---|---|---|---|---|---|
| D0D0 | pnetd | 1.030 | 0.086 | 1.091 / 0.913 | 0.80 |
| D0D0 | xgbd | 1.002 | 0.085 | 1.114 / 0.892 | 0.81 |
| D0D0 | base_sum | 1.047 | 0.089 | 1.032 / 0.968 | 0.81 |
| DpD0 | pnetd | 0.989 | 0.070 | 0.952 / 1.048 | 0.76 |
| DpD0 | xgbd | 1.002 | 0.067 | 0.969 / 1.033 | 0.79 |
| DpD0 | base_sum | 1.023 | 0.068 | 1.017 / 0.978 | 0.79 |
| DpDp | pnetd | 0.896 | 0.089 | 1.002 / 0.995 | 0.81 |
| DpDp | xgbd | 0.919 | 0.086 | 0.952 / 1.044 | 0.83 |
| DpDp | base_sum | 1.021 | 0.092 | 0.884 / 1.131 | 0.85 |

## Dobroć fitów (dewiancja / ndf, kanały)

- nomcat_fact: D0D0 3442/3707, DpD0 6171/5561, DpDp 4030/3707
- nomcat_mc: D0D0 3442/3707, DpD0 6171/5561, DpDp 4030/3707
- pnetd: D0D0 3277/3707, DpD0 3761/3707, DpDp 3897/3707
- xgbd: D0D0 3259/3707, DpD0 3824/3707, DpDp 3856/3707
- base_sum: D0D0 3408/3707, DpD0 3923/3707, DpDp 4042/3707

Uwagi: σ z profilu bez korekty korelacji zdarzeń (bootstrap Poissona po zdarzeniach ją zawiera); normalizacja i tło D* tag poza modelem (nomcat_* bez kategorii z tagiem D* → nie jest identyczny z nominałem); progi sieci z MC, więc E_k = 0,65²/K × SF_rel.
