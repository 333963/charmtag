# TRUTHNORM — μ par D–D z normalizacją z prawdy mb_nd (decyzja z 19.11; GEN_mcstat/PROCEDURA.md §6)

n_pair(ch) × f_ch, f = {'D0D0': 1.031693445545306, 'DpD0': 1.2659039287655927, 'DpDp': 1.059714967709065} (`GEN_mcstat/out/truth_mbnd.json`); tło bez zmian; pełny ponowny fit (`truthnorm_run.py`, ub_jk.replica_fit, charmtag PF_pairs/src). Kontrola f ≡ 1 = odtworzenie nominału.

## (a) 428770 i (c) 427929

| bieg | μ D0D0 | μ DpD0 | μ DpDp | μ wspólne [68 %, profil] | Z |
|---|---|---|---|---|---|
| 428770 kontrola f ≡ 1 (próbka 05.11) | 1.16 ± 0.55 | 1.16 ± 0.42 | 1.39 ± 0.93 | **1.19** [0.88; 1.51] | 4.07 |
| 428770 prawda mb_nd (próbka 05.11, N_D* 5584) | 1.12 ± 0.54 | 0.92 ± 0.33 | 1.31 ± 0.87 | **1.01** [0.74; 1.29] | 4.05 |
| 427929 kontrola f ≡ 1 | 0.55 ± 0.70 | 0.77 ± 0.59 | -0.07 ± 0.94 | **0.57** [0.19; 0.99] | 1.57 |
| 427929 prawda mb_nd | 0.53 ± 0.68 | 0.61 ± 0.47 | -0.07 ± 0.89 | **0.50** [0.17; 0.86] | 1.60 |
| 428770 prawda mb_nd — skala zapisu (N_D* 4436) | 1.41 ± 0.67 | 1.16 ± 0.42 | 1.65 ± 1.10 | **1.27** [0.94; 1.62] | 4.05 |

Kontrole: 428770 f ≡ 1 → 1.1889 (ref. 1,1889: OK); 427929 f ≡ 1 → 0.5674 (ref. jackknife 0,567, nogi przebudowane; UB427 0,53: OK).

## Zapis wyników (stat × jackknife zdarzeń; syst bez N_D* i N_D* względnie z budżetów)

- **428770:** μ = 1.27 +0.41/−0.40 (stat, ×1.18) ± 0.36 (syst, prawda mb_nd) ± 0.28 (N_D*); Z = 4.05 (profil; ≈ 3.43 z korelacjami zdarzeń). Trigger: osobno (TRIG_A/SYST_TRIG.md).
- **427929 (wynik główny, decyzja z 19.11):** μ = 0.50 +0.41/−0.37 (stat, ×1.13) +0.28/−0.19 (syst) ± 0.21 (N_D*); Z = 1.60 (≈ 1.42); vs Pythia (μ = 1): -0.93σ (σ_tot górne ≈ 0.54, syst względne z BUDZET_427929 przy μ jak w biegu).

## (b) 428770 per klasa triggera (prawda mb_nd; r_c, nogi i tło jak UBTRIG)

| klasa | r_c | pary D0D0/DpD0/DpDp (= UBTRIG?) | μ D0D0 | μ DpD0 | μ DpDp | μ wspólne [68 %] | Z |
|---|---|---|---|---|---|---|---|
| J12only | 0.769 ± 0.019 | 8163/24594/26984 (tak) | 1.91 ± 0.68 | 1.14 ± 0.40 | 1.47 ± 1.03 | **1.37** [1.05; 1.71] | 4.53 |
| MBTS_j20 | 0.135 ± 0.015 | 838/2650/3128 (tak) | -4.73 ± 1.58 | 0.17 ± 0.75 | -1.25 ± 1.60 | **-0.51** [-0.95; -0.05] | 0.00 |
| MBTS_noj20 | 0.010 ± 0.004 | 47/117/167 (tak) | -20.00 ± 39.11 | 2.01 ± 1.90 | 4.66 ± 7.46 | **1.28** [0.40; 2.87] | 1.57 |
| other | 0.084 ± 0.012 | 379/1340/1527 (tak) | -1.93 ± 2.44 | -0.78 ± 1.07 | 4.22 ± 3.26 | **-0.59** [-1.31; 0.11] | 0.00 |

J12only vs reszta (próbka 05.11): 1.37 ± 0.33 vs -0.40 ± 0.37 → +3.56σ (bez czynnika jackknife; ≈ +3.02σ z nim). Na skali zapisu (×1.259): J12only 1.73, reszta -0.50.

Uwagi: (1) 427929 — f_ch z 428770 (metoda stosunku R_X; przybliżenie: R_X liczony na kandydatach reko, przyjęty równy dla prawdy). (2) Syst względne przeniesione z budżetów (RAPORT_gen_norm: 0,42/1,50; BUDZET_427929: 0,387 bez N_D*, N_D* 0,284 — ten ostatni czeka na ub427nds). (3) Nominał w UNB_WYNIK.md NIE zmieniony.
