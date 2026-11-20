# Pary D–D różniczkowo — znak (like/opposite) (ODSŁONIĘTE, etap PD)

- Zmienna: `sign`, biny: like, opposite; normalizacja μ_i: `mc` (SS_i = μ_i · N_par · f_i); kanały fitowane: D0D0, DpD0, DpDp.
- Pary MC minbias (Pythia) per bin: D0D0 [0, 0]; DpD0 [0, 0]; DpDp [3, 12]; nieokreślone w danych (poza binami): 38128.
- Nominał (łącznie, bez binów): wspólne μ = 1.25 [0.90, 1.51], Z = 4.06σ.

| kanał | bin | N par | S oczek. | f_fit (tryb) | f_MC ± | μ_i | 68 % | Z_i | dev/ndf | f_i dane ± |
|---|---|---|---|---|---|---|---|---|---|---|
| D0D0 | like | — | — | — | — | za mało par (0) albo S_oczek = 0 | | | | |
| D0D0 | opposite | — | — | — | — | za mało par (0) albo S_oczek = 0 | | | | |
| DpD0 | like | — | — | — | — | za mało par (0) albo S_oczek = 0 | | | | |
| DpD0 | opposite | — | — | — | — | za mało par (0) albo S_oczek = 0 | | | | |
| DpDp | like | 15160 | 3.7 | 0.20 (mc) | 0.20 ± 0.10 | 4.47 ± 4.47 | 0.20–5.00 | 1.05 | 1.03 | 0.49 ± 0.35 |
| DpDp | opposite | 16646 | 14.9 | 0.80 (mc) | 0.80 ± 0.10 | 1.15 ± 1.15 | 0.06–2.35 | 1.05 | 1.05 | 0.51 ± 0.35 |
| DpDp | μ(like)/μ(opposite) | | | | | 3.88 ± 5.49 | | | | |
| DpDp | test `dps_flat` h=[0.5, 0.5] | | | | | μ=1.80 | 0.58–3.15 | | χ²/ndf=0.01/1 | p=0.91 |
| DpDp | test `pythia_mc` h=[0.2, 0.8] | | | | | μ=1.50 | 0.30–2.56 | | χ²/ndf=0.56/1 | p=0.46 |

Łączenie po kanałach per bin (suma profili, jak nominalne `combine`):

| bin | kanały | μ | 68 % | Z |
|---|---|---|---|---|
| like | DpDp | 4.50 | 0.20–9.17 | 1.05 |
| opposite | DpDp | 1.25 | 0.06–2.35 | 1.05 |

- Test łączny `pythia_mc`: μ globalne = 1.50 [0.30, 2.56] (Z = 1.29σ); Δχ²(μ per bin wspólne w kanałach vs jedno μ z kształtem) = 0.55 / ndf 1 → p = 0.46.
- Test łączny `dps_flat`: μ globalne = 1.80 [0.58, 3.15] (Z = 1.48σ); Δχ²(μ per bin wspólne w kanałach vs jedno μ z kształtem) = 0.01 / ndf 1 → p = 0.94.

Like-sign (INT): f_LS = N_LS/(N_LS+N_OS) w sygnale; LHCb LS/OS 0.10–0.13 → f_LS 0.091–0.115; DPS czyste → 0,5; SPS czyste → 0.

| kanał | N_LS sig | N_OS sig | f_LS dane ± | f_LS Pythia ± | Z odrzucenia Pythii | Z odrzucenia DPS |
|---|---|---|---|---|---|---|
| DpDp | 39.3 | 40.5 | 0.492 ± 0.353 | 0.200 ± 0.103 | 0.75 | 0.11 |
- łącznie (kanały): Z_reject_pythia_mc_combined = 0.74σ
- łącznie (kanały): Z_reject_dps_flat_combined = 0.07σ
- Kontrola: μ wspólne we wszystkich (kanał × bin) = 1.50 [0.30, 2.56], Z = 1.29σ (oczekiwane ≈ nominał dla kanałów fitowanych).

Interpretacja: `pythia_mc` = zgodność kształtu z Pythią (SPS+MPI), `dps_flat` = zgodność z czystym DPS (Δφ płasko / like:opposite 1:1). Małe p odrzuca hipotezę. f_MC z N par MC minbias — jego błąd nie jest wliczony w μ_i (patrz kolumna f_MC ±).
