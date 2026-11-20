# N_D* pełna lumi — skan wariantów (NDS, naprawiony `dstar_refit.fit_single`)

| wariant | kat. | okno D0 | tło | N ± stat | skala σ | μ [MeV] | valid | χ²/ndf |
|---|---|---|---|---|---|---|---|---|
| nominal | nominal | w40 | threshold_pol2 | 4436 ± 595 | 1.03 | 145.40 | ✔ | 336/383 |
| bkg_expo_pol2_fullmodel | bkg | w40 | expo_pol2 | 4433 ± 580 | 1.03 | 145.41 | ✔ | 336/382 |
| bkg_sqrt_pol3 | bkg | w40 | sqrt_pol3 | 4964 ± 642 | 1.09 | 145.37 | ✔ | 338/383 |
| bkg_pure_expo | bkg_rigid | w40 | pure_expo | 5553 ± 633 | 1.17 | 145.40 | ✔ | 341/384 |
| range_hi158 | range | w40 | threshold_pol2 | 4240 ± 583 | 1.00 | 145.40 | ✔ | 308/343 |
| range_hi165 | range | w40 | threshold_pol2 | 4436 ± 595 | 1.03 | 145.40 | ✔ | 336/383 |
| range_lo141 | range | w40 | threshold_pol2 | 4317 ± 584 | 1.01 | 145.39 | ✔ | 326/373 |
| range_lo141_hi165 | range | w40 | threshold_pol2 | 4317 ± 584 | 1.01 | 145.39 | ✔ | 326/373 |
| bin_0.10 | bin | w40 | threshold_pol2 | 4458 ± 680 | 1.02 | 145.40 | ✔ | 179/188 |
| bin_0.20 | bin | w40 | threshold_pol2 | 4342 ± 645 | 1.01 | 145.39 | ✔ | 101/90 |
| sig_gauss_cb | sig | w40 | threshold_pol2 | 4295 ± 601 | 1.02 | 145.39 | ✖ | 336/383 |
| sig_scale_fixed | sig | w40 | threshold_pol2 | 4347 ± 400 | 1.00 | 145.40 | ✔ | 336/384 |
| window_d0_36 | window | w36 | threshold_pol2 | 3999 ± 546 | 1.00 | 145.39 | ✔ | 348/383 |
| window_d0_55 | window | w55 | threshold_pol2 | 5188 ± 891 | 1.16 | 145.46 | ✔ | 373/383 |

| bin pT | N ± stat | skala σ | valid | tło |
|---|---|---|---|---|
| 3-5 | 646 ± 305 | 0.80 | ✔ | threshold_pol2 |
| 5-8 | 1353 ± 301 | 1.08 | ✔ | threshold_pol2 |
| 8-30 | 2174 ± 175 | 0.95 | ✔ | threshold_pol2 |

**N_D* = 4436 ± 595 (stat) ± 588 (syst)**; syst per kat.: {"bkg": 527, "range": 196, "bin": 95, "sig": 141}; z „czystym wykładnikiem” (Δχ²=+5, σ×1.17): ±1262. Mediana w40 4390, rozrzut 190 (10 wariantów). Suma per bin 4172.
Stare: fit_dm (pf_pairs, użyte w unblind) 5588 ± 586; fit_dm (V) 6607; refit V (invalid) 4472. Korekta μ: ×1.260.
