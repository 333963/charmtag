# N_D* pełna lumi — skan wariantów (NDS, naprawiony `dstar_refit.fit_single`)

| wariant | kat. | okno D0 | tło | N ± stat | skala σ | μ [MeV] | valid | χ²/ndf |
|---|---|---|---|---|---|---|---|---|
| nominal | nominal | w40 | threshold_pol2 | 4436 ± 671 | 1.02 | 145.40 | ✔ | 336/383 |
| bkg_expo_pol2_fullmodel | bkg | w40 | expo_pol2 | 4432 ± 683 | 1.02 | 145.40 | ✔ | 336/382 |
| bkg_sqrt_pol3 | bkg | w40 | sqrt_pol3 | 4970 ± 660 | 1.09 | 145.37 | ✔ | 338/383 |
| bkg_pure_expo | bkg_rigid | w40 | pure_expo | 5553 ± 623 | 1.17 | 145.40 | ✔ | 341/384 |
| range_hi158 | range | w40 | threshold_pol2 | 4226 ± 684 | 1.00 | 145.41 | ✔ | 308/343 |
| range_hi165 | range | w40 | threshold_pol2 | 4436 ± 671 | 1.02 | 145.40 | ✔ | 336/383 |
| range_lo141 | range | w40 | threshold_pol2 | 4317 ± 667 | 1.01 | 145.39 | ✔ | 326/373 |
| range_lo141_hi165 | range | w40 | threshold_pol2 | 4317 ± 667 | 1.01 | 145.39 | ✔ | 326/373 |
| bin_0.10 | bin | w40 | threshold_pol2 | 4443 ± 675 | 1.02 | 145.40 | ✔ | 179/188 |
| bin_0.20 | bin | w40 | threshold_pol2 | 4330 ± 664 | 1.01 | 145.39 | ✔ | 101/90 |
| sig_gauss_cb | sig | w40 | threshold_pol2 | 4278 ± 714 | 1.02 | 145.40 | ✔ | 336/383 |
| sig_scale_fixed | sig | w40 | threshold_pol2 | 4339 ± 406 | 1.00 | 145.40 | ✔ | 336/384 |
| window_d0_36 | window | w36 | threshold_pol2 | 4000 ± 613 | 1.00 | 145.39 | ✔ | 348/383 |
| window_d0_55 | window | w55 | threshold_pol2 | 5164 ± 899 | 1.16 | 145.46 | ✔ | 373/383 |

| bin pT | N ± stat | skala σ | valid | tło |
|---|---|---|---|---|
| 3-5 | 641 ± 282 | 0.80 | ✔ | threshold_pol2 |
| 5-8 | 1345 ± 303 | 1.08 | ✔ | threshold_pol2 |
| 8-30 | 2173 ± 192 | 0.95 | ✔ | threshold_pol2 |

**N_D* = 4436 ± 671 (stat) ± 605 (syst)**; syst per kat.: {"bkg": 534, "range": 211, "bin": 107, "sig": 159}; z „czystym wykładnikiem” (Δχ²=+5, σ×1.17): ±1270. Mediana w40 4339, rozrzut 190 (11 wariantów). Suma per bin 4158.
Stare: fit_dm (pf_pairs, użyte w unblind) 5588 ± 586; fit_dm (V) 6607; refit V (invalid) 4472. Korekta μ: ×1.260.
