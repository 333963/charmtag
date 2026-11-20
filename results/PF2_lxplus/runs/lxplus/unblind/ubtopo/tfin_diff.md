# PD — pary różniczkowe `dy` (raport końcowy, 2026-11-19 04:52 UTC)

- Źródła: `/eos/user/k/klyczak/charm_q/ub_topo/dy` (pairs/fit2d --diff) i `/eos/user/k/klyczak/charm_q/pf_pairs/unblind` (nominał odsłonięty).
- Kontrola spójności z nominałem: Σ_i μ_i f_i (sygnał per kanał w jednostkach oczekiwania nominalnego) vs μ z fit2d_unblind.json (różnice ≈ pary nieokreślone/poza binami + niezależne nuisance per bin).

| kanał | Σ μ_i f_i | μ nominalne ± | pary nieokreślone (dane) |
|---|---|---|---|
| D0D0 | 1.57 | 1.16 ± 0.55 | 0 (łącznie) |
| DpD0 | 1.60 | 1.16 ± 0.42 | 0 (łącznie) |
| DpDp | 1.97 | 1.39 ± 0.93 | 0 (łącznie) |

- Nominał wspólne μ = 1.25 [0.8999999999999999, 1.5125000000000002], Z = 4.06σ.

---

# Pary D–D różniczkowo — Δy≈|Δη| (ODSŁONIĘTE, etap PD)

- Zmienna: `dy`, biny: [0,1), [1,5]; normalizacja μ_i: `mc` (SS_i = μ_i · N_par · f_i); kanały fitowane: D0D0, DpD0, DpDp.
- Pary MC minbias (Pythia) per bin: D0D0 [9, 13]; DpD0 [26, 15]; DpDp [12, 3]; nieokreślone w danych (poza binami): 0.
- Nominał (łącznie, bez binów): wspólne μ = 1.25 [0.90, 1.51], Z = 4.06σ.

| kanał | bin | N par | S oczek. | f_fit (tryb) | f_MC ± | μ_i | 68 % | Z_i | dev/ndf | f_i dane ± |
|---|---|---|---|---|---|---|---|---|---|---|
| D0D0 | [0,1) | 5691 | 17.1 | 0.41 (mc) | 0.41 ± 0.10 | 2.71 ± 1.62 | 1.17–4.39 | 1.80 | 0.62 | 0.71 ± 0.22 |
| D0D0 | [1,5] | 3736 | 24.7 | 0.59 (mc) | 0.59 ± 0.10 | 0.78 ± 0.67 | 0.18–1.49 | 1.28 | 0.51 | 0.29 ± 0.22 |
| D0D0 | μ([0,1))/μ([1,5]) | | | | | 3.46 ± 3.59 | | | | |
| D0D0 | test `dps_flat` h=[0.5, 0.5] | | | | | μ=1.23 | 0.64–2.04 | | χ²/ndf=0.76/1 | p=0.38 |
| D0D0 | test `pythia_mc` h=[0.41, 0.59] | | | | | μ=1.00 | 0.49–1.80 | | χ²/ndf=1.33/1 | p=0.25 |
| DpD0 | [0,1) | 18610 | 41.3 | 0.63 (mc) | 0.63 ± 0.08 | 2.36 ± 0.81 | 1.58–3.19 | 3.23 | 0.92 | 0.93 ± 0.15 |
| DpD0 | [1,5] | 10091 | 23.8 | 0.37 (mc) | 0.37 ± 0.08 | 0.29 ± 0.68 | -0.30–1.04 | 0.45 | 0.77 | 0.07 ± 0.15 |
| DpD0 | μ([0,1))/μ([1,5]) | | | | | 8.13 ± 19.36 | | | | |
| DpD0 | test `dps_flat` h=[0.5, 0.5] | | | | | μ=0.95 | 0.51–1.56 | | χ²/ndf=6.23/1 | p=0.01 |
| DpD0 | test `pythia_mc` h=[0.63, 0.37] | | | | | μ=1.50 | 0.83–1.98 | | χ²/ndf=3.60/1 | p=0.06 |
| DpDp | [0,1) | 24025 | 14.9 | 0.80 (mc) | 0.80 ± 0.10 | 1.42 ± 1.38 | 0.10–2.85 | 1.07 | 1.09 | 0.58 ± 0.30 |
| DpDp | [1,5] | 7781 | 3.7 | 0.20 (mc) | 0.20 ± 0.10 | 4.16 ± 3.26 | 1.08–5.00 | 1.39 | 1.01 | 0.42 ± 0.30 |
| DpDp | μ([0,1))/μ([1,5]) | | | | | 0.34 ± 0.43 | | | | |
| DpDp | test `dps_flat` h=[0.5, 0.5] | | | | | μ=1.80 | 0.73–3.00 | | χ²/ndf=0.07/1 | p=0.80 |
| DpDp | test `pythia_mc` h=[0.8, 0.2] | | | | | μ=2.00 | 0.67–3.23 | | χ²/ndf=0.64/1 | p=0.42 |

Łączenie po kanałach per bin (suma profili, jak nominalne `combine`):

| bin | kanały | μ | 68 % | Z |
|---|---|---|---|---|
| [0,1) | D0D0, DpD0, DpDp | 2.25 | 1.62–2.87 | 3.79 |
| [1,5] | D0D0, DpD0, DpDp | 0.75 | 0.21–1.17 | 1.45 |

- Test łączny `pythia_mc`: μ globalne = 1.25 [0.94, 1.76] (Z = 3.56σ); Δχ²(μ per bin wspólne w kanałach vs jedno μ z kształtem) = 3.75 / ndf 1 → p = 0.05.
- Test łączny `dps_flat`: μ globalne = 1.23 [0.83, 1.62] (Z = 3.33σ); Δχ²(μ per bin wspólne w kanałach vs jedno μ z kształtem) = 5.77 / ndf 1 → p = 0.02.
- Kontrola: μ wspólne we wszystkich (kanał × bin) = 1.25 [0.94, 1.76], Z = 3.56σ (oczekiwane ≈ nominał dla kanałów fitowanych).

Interpretacja: `pythia_mc` = zgodność kształtu z Pythią (SPS+MPI), `dps_flat` = zgodność z czystym DPS (Δφ płasko / like:opposite 1:1). Małe p odrzuca hipotezę. f_MC z N par MC minbias — jego błąd nie jest wliczony w μ_i (patrz kolumna f_MC ±).
