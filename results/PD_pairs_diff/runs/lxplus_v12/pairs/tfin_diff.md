# Pary D–D różniczkowo — znak (like/opposite) (ODSŁONIĘTE, etap PD)

- Zmienna: `sign`, biny: like, opposite; normalizacja μ_i: `mc` (SS_i = μ_i · N_par · f_i); kanały fitowane: D0D0, DpD0, DpDp.
- Pary MC minbias (Pythia) per bin: D0D0 [8, 14]; DpD0 [13, 28]; DpDp [3, 12]; nieokreślone w danych (poza binami): 0.
- Nominał (łącznie, bez binów): — (brak fit2d_unblind.json).

| kanał | bin | N par | S oczek. | f_fit (tryb) | f_MC ± | μ_i | 68 % | Z_i | dev/ndf | f_i dane ± |
|---|---|---|---|---|---|---|---|---|---|---|
| D0D0 | like | 4728 | 15.2 | 0.36 (mc) | 0.36 ± 0.10 | 0.96 ± 1.42 | -0.41–2.43 | 0.69 | 0.55 | 0.22 ± 0.27 |
| D0D0 | opposite | 4699 | 26.6 | 0.64 (mc) | 0.64 ± 0.10 | 1.95 ± 0.91 | 1.09–2.89 | 2.37 | 0.57 | 0.78 ± 0.27 |
| D0D0 | μ(like)/μ(opposite) | | | | | 0.49 ± 0.76 | | | | |
| D0D0 | test `dps_flat` h=[0.5, 0.5] | | | | | μ=1.59 | 0.84–2.36 | | χ²/ndf=1.35/1 | p=0.25 |
| D0D0 | test `pythia_mc` h=[0.36, 0.64] | | | | | μ=1.75 | 0.96–2.46 | | χ²/ndf=0.35/1 | p=0.55 |
| DpD0 | like | 14378 | 20.7 | 0.32 (mc) | 0.32 ± 0.07 | 1.78 ± 1.25 | 0.61–3.08 | 1.56 | 0.86 | 0.34 ± 0.18 |
| DpD0 | opposite | 14323 | 44.5 | 0.68 (mc) | 0.68 ± 0.07 | 1.57 ± 0.62 | 0.98–2.21 | 2.83 | 0.85 | 0.66 ± 0.18 |
| DpD0 | μ(like)/μ(opposite) | | | | | 1.13 ± 0.91 | | | | |
| DpD0 | test `dps_flat` h=[0.5, 0.5] | | | | | μ=1.71 | 1.08–2.23 | | χ²/ndf=0.79/1 | p=0.37 |
| DpD0 | test `pythia_mc` h=[0.32, 0.68] | | | | | μ=1.50 | 1.08–2.19 | | χ²/ndf=0.07/1 | p=0.80 |
| DpDp | like | 15160 | 3.7 | 0.20 (mc) | 0.20 ± 0.10 | 4.48 ± 4.47 | 0.20–5.00 | 1.05 | 1.03 | 0.49 ± 0.35 |
| DpDp | opposite | 16646 | 14.9 | 0.80 (mc) | 0.80 ± 0.10 | 1.15 ± 1.15 | 0.06–2.35 | 1.05 | 1.05 | 0.51 ± 0.35 |
| DpDp | μ(like)/μ(opposite) | | | | | 3.89 ± 5.48 | | | | |
| DpDp | test `dps_flat` h=[0.5, 0.5] | | | | | μ=1.80 | 0.58–3.15 | | χ²/ndf=0.01/1 | p=0.91 |
| DpDp | test `pythia_mc` h=[0.2, 0.8] | | | | | μ=1.50 | 0.30–2.56 | | χ²/ndf=0.56/1 | p=0.46 |

Łączenie po kanałach per bin (suma profili, jak nominalne `combine`):

| bin | kanały | μ | 68 % | Z |
|---|---|---|---|---|
| like | D0D0, DpD0, DpDp | 1.50 | 0.69–2.53 | 1.83 |
| opposite | D0D0, DpD0, DpDp | 1.50 | 1.16–2.09 | 3.79 |

- Test łączny `pythia_mc`: μ globalne = 1.50 [1.20, 2.04] (Z = 4.21σ); Δχ²(μ per bin wspólne w kanałach vs jedno μ z kształtem) = 0.00 / ndf 1 → p = 1.00.
- Test łączny `dps_flat`: μ globalne = 1.64 [1.22, 2.09] (Z = 4.07σ); Δχ²(μ per bin wspólne w kanałach vs jedno μ z kształtem) = 1.74 / ndf 1 → p = 0.19.

Like-sign (INT): f_LS = N_LS/(N_LS+N_OS) w sygnale; LHCb LS/OS 0.10–0.13 → f_LS 0.091–0.115; DPS czyste → 0,5; SPS czyste → 0.

| kanał | N_LS sig | N_OS sig | f_LS dane ± | f_LS Pythia ± | Z odrzucenia Pythii | Z odrzucenia DPS |
|---|---|---|---|---|---|---|
| D0D0 | 31.5 | 112.5 | 0.219 ± 0.266 | 0.364 ± 0.103 | 0.59 | 1.16 |
| DpD0 | 83.2 | 158.2 | 0.345 ± 0.181 | 0.317 ± 0.073 | 0.26 | 0.89 |
| DpDp | 39.4 | 40.6 | 0.493 ± 0.353 | 0.200 ± 0.103 | 0.75 | 0.11 |
- łącznie (kanały): Z_reject_pythia_mc_combined = 0.00σ
- łącznie (kanały): Z_reject_dps_flat_combined = 1.32σ
- Kontrola: μ wspólne we wszystkich (kanał × bin) = 1.50 [1.20, 2.04], Z = 4.21σ (oczekiwane ≈ nominał dla kanałów fitowanych).

Interpretacja: `pythia_mc` = zgodność kształtu z Pythią (SPS+MPI), `dps_flat` = zgodność z czystym DPS (Δφ płasko / like:opposite 1:1). Małe p odrzuca hipotezę. f_MC z N par MC minbias — jego błąd nie jest wliczony w μ_i (patrz kolumna f_MC ±).
