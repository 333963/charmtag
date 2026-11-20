# UBPNET — μ par D–D w 427929 z kategoriami sieci par (PAIRNET v3) — odsłonięte (zatwierdzone 18.11)

Próbka data_427929: SF z HYBSF/HYBDS dla tej próbki; mnożniki n_pair (normalizacja z prawdy mb_nd f_ch, PROCEDURA §6): 1.032,1.266,1.060; nomcat_mc_f1 = nomcat_mc bez mnożników (kontrola zgodności z nominałem charmtag).

Ten sam model fitu 2D (sygnał MC × μ + 3 składowe tła z wolnymi normalizacjami i kształtami exp(poly2) per kategoria) i te same pary danych (h5 v2.0 427929 po filtrze UB427: L_xy ≥ 0,25, okno D0, klasa MB) dla wszystkich schematów; nomcat_* = nominalne kategorie nóg (bez kategorii z tagiem D*), pnetd = sieć par z tłem z SB×SB danych (schemat główny), xgbd / base_sum = kontrola. Normalizacja n_pair z pair_inputs odsłonięcia; schemat sieci × SF_rel (pary hybrydowe: D⁺ HYBSF v3, D⁰ HYBDS v3).

**Kontrola nominału (charmtag, pairs_full.npz):** D0D0 0.5458 (ref 0.5458), DpD0 0.7667 (ref 0.7667), DpDp -0.0712 (ref -0.0712); wspólne μ_nom = 0.5674; zgodność: TAK. Skala do zapisu bieżącego ×1.000.

Dane: nogi 685685 → w zdarzeniach próbki nominalnej 685685 (zdarzeń z nogami 234894 z 234894); pary w zakresie mas {'D0D0': 62822, 'DpD0': 190482, 'DpDp': 170195}; pokrycie zdarzeń nóg nominału przez nogi v2.0: –.

| schemat | μ D0D0 | μ DpD0 | μ DpDp | μ wspólne [68 %] | Z | σ_68 / σ(nomcat_mc) | σ Asimov / σ Asimov(nomcat_mc) | σ bootstrap | μ × skala (zapis bieżący) |
|---|---|---|---|---|---|---|---|---|---|
| nomcat_fact | 0.09 ± 1.16 | -0.11 ± 0.58 | -1.00 ± 1.01 | **-0.244** [-0.703; 0.239] | 0.00 | 0.991 | 0.994 | 0.664 | -0.244 |
| nomcat_mc | 0.05 ± 1.21 | -0.11 ± 0.58 | -0.97 ± 1.02 | **-0.246** [-0.710; 0.241] | 0.00 | 1.000 | 1.000 | 0.704 | -0.246 |
| nomcat_mc_f1 | 0.05 ± 1.25 | -0.13 ± 0.74 | -1.03 ± 1.08 | **-0.309** [-0.861; 0.265] | 0.00 | 1.183 | 1.129 | 0.797 | -0.309 |
| pnetd | -0.02 ± 0.28 | 0.22 ± 0.34 | 0.51 ± 0.81 | **0.098** [-0.108; 0.395] | 0.37 | 0.529 | 0.476 | 0.235 | 0.098 |
| xgbd | -0.00 ± 0.07 | 0.18 ± 0.27 | -0.19 ± 0.86 | **0.100** [-0.062; 0.393] | 0.53 | 0.478 | 0.414 | 0.518 | 0.100 |
| base_sum | 1.54 ± 0.61 | 0.39 ± 0.55 | 0.41 ± 1.23 | **0.715** [0.319; 1.262] | 1.96 | 0.991 | 0.827 | 0.714 | 0.715 |

## Różnica względem nomcat_mc (te same pary, ten sam model) i systematyki schematu

| schemat | Δμ = μ − μ(nomcat_mc) | σ_Δ (bootstrap) | Δ/σ_Δ | korelacja (bootstrap) | syst SF norm. (δμ) | syst nachylenie SF (δμ) | syst tło poly3 (δμ) | Asimov: μ0 → μ̂ |
|---|---|---|---|---|---|---|---|---|
| nomcat_fact | 0.002 | 0.253 | 0.01 | 0.93 | – | – | -0.185 | -0.244 → -0.244 |
| nomcat_mc_f1 | -0.063 | 0.351 | -0.18 | 0.90 | – | – | -0.205 | -0.309 → -0.309 |
| pnetd | 0.344 | 0.688 | 0.50 | 0.24 | 0.007 | -0.013 | -0.099 | 0.098 → 0.098 |
| xgbd | 0.346 | 0.592 | 0.58 | 0.57 | 0.007 | 0.006 | -0.027 | 0.100 → 0.100 |
| base_sum | 0.961 | 0.795 | 1.21 | 0.37 | 0.049 | 0.051 | -0.587 | 0.715 → 0.716 |

## SF_rel schematów (SF schematu / SF nomcat; para = iloczyn nóg)

| kanał | schemat | SF_rel | ± | nachylenie góra / dół | ρ(schemat, nomcat) MC |
|---|---|---|---|---|---|
| D0D0 | pnetd | 1.079 | 0.083 | 1.021 / 0.979 | 0.80 |
| D0D0 | xgbd | 1.085 | 0.082 | 0.950 / 1.052 | 0.81 |
| D0D0 | base_sum | 1.019 | 0.077 | 0.987 / 1.013 | 0.81 |
| DpD0 | pnetd | 1.009 | 0.061 | 0.950 / 1.050 | 0.76 |
| DpD0 | xgbd | 0.998 | 0.057 | 0.900 / 1.103 | 0.79 |
| DpD0 | base_sum | 1.011 | 0.058 | 0.941 / 1.063 | 0.79 |
| DpDp | pnetd | 0.922 | 0.072 | 1.083 / 0.920 | 0.81 |
| DpDp | xgbd | 0.928 | 0.069 | 0.996 / 1.005 | 0.83 |
| DpDp | base_sum | 0.899 | 0.063 | 1.027 / 0.962 | 0.85 |

## Transfer tła SBTRAIN (pary SB×SB: ta próbka vs referencja 428770, ten sam model sieci)

Tło w ficie jest wolne per kategoria — różnica rozkładów nie obciąża μ, mierzy tylko, czy sieć rozdziela tło tak jak na próbce treningowej.

| kanał | schemat | n SB×SB (ref) | χ²/ndf | w schemacie: próbka / ref | górna połowa kategorii: próbka / ref |
|---|---|---|---|---|---|
| D0D0 | nomcat_mc | 9594 (85057) | 3.4/6 | 0.014 / 0.014 | 0.135 / 0.162 |
| D0D0 | pnetd | 9594 (85057) | 5.2/6 | 0.003 / 0.005 | 0.094 / 0.126 |
| D0D0 | xgbd | 9594 (85057) | 6.6/6 | 0.003 / 0.005 | 0.100 / 0.113 |
| D0D0 | base_sum | 9594 (85057) | 1.9/6 | 0.006 / 0.006 | 0.123 / 0.123 |
| DpD0 | nomcat_mc | 50494 (460683) | 10.9/9 | 0.014 / 0.015 | 0.108 / 0.131 |
| DpD0 | pnetd | 50494 (460683) | 10.0/6 | 0.004 / 0.005 | 0.113 / 0.134 |
| DpD0 | xgbd | 50494 (460683) | 9.5/6 | 0.004 / 0.005 | 0.093 / 0.133 |
| DpD0 | base_sum | 50494 (460683) | 2.4/6 | 0.008 / 0.008 | 0.130 / 0.134 |
| DpDp | nomcat_mc | 78714 (743467) | 15.9/6 | 0.017 / 0.019 | 0.269 / 0.241 |
| DpDp | pnetd | 78714 (743467) | 37.1/6 | 0.005 / 0.007 | 0.167 / 0.146 |
| DpDp | xgbd | 78714 (743467) | 43.4/6 | 0.006 / 0.008 | 0.153 / 0.143 |
| DpDp | base_sum | 78714 (743467) | 5.3/6 | 0.010 / 0.010 | 0.181 / 0.174 |

## Dobroć fitów (dewiancja / ndf, kanały)

- nomcat_fact: D0D0 2019/3707, DpD0 4020/5561, DpDp 3343/3707
- nomcat_mc: D0D0 2019/3707, DpD0 4020/5561, DpDp 3343/3707
- nomcat_mc_f1: D0D0 2019/3707, DpD0 4020/5561, DpDp 3343/3707
- pnetd: D0D0 1129/3707, DpD0 2068/3707, DpDp 2241/3707
- xgbd: D0D0 1118/3707, DpD0 2173/3707, DpDp 2435/3707
- base_sum: D0D0 1501/3707, DpD0 2731/3707, DpDp 3112/3707

Uwagi: σ z profilu bez korekty korelacji zdarzeń (bootstrap Poissona po zdarzeniach ją zawiera); normalizacja i tło D* tag poza modelem (nomcat_* bez kategorii z tagiem D* → nie jest identyczny z nominałem); progi sieci z MC, więc E_k = 0,65²/K × SF_rel.
