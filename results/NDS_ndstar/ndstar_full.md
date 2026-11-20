# N_D* pełna lumi — skan wariantów (NDS, naprawiony `dstar_refit.fit_single`)

| wariant | kat. | okno D0 | tło | N ± stat | skala σ | μ [MeV] | valid | χ²/ndf |
|---|---|---|---|---|---|---|---|---|
| nominal | nominal | w40 | threshold_pol2 | 4443 ± 612 | 1.03 | 145.40 | ✔ | 336/383 |
| bkg_expo_pol2_fullmodel | bkg | w40 | expo_pol2 | 4442 ± 567 | 1.03 | 145.41 | ✔ | 336/382 |
| bkg_sqrt_pol3 | bkg | w40 | sqrt_pol3 | 4963 ± 609 | 1.09 | 145.37 | ✔ | 338/383 |
| bkg_pure_expo | bkg_rigid | w40 | pure_expo | 5553 ± 612 | 1.17 | 145.40 | ✔ | 341/384 |
| range_hi158 | range | w40 | threshold_pol2 | 4240 ± 615 | 1.00 | 145.40 | ✔ | 308/343 |
| range_hi165 | range | w40 | threshold_pol2 | 4443 ± 612 | 1.03 | 145.40 | ✔ | 336/383 |
| range_lo141 | range | w40 | threshold_pol2 | 4317 ± 581 | 1.01 | 145.39 | ✔ | 326/373 |
| range_lo141_hi165 | range | w40 | threshold_pol2 | 4317 ± 581 | 1.01 | 145.39 | ✔ | 326/373 |
| bin_0.10 | bin | w40 | threshold_pol2 | 4471 ± 608 | 1.03 | 145.40 | ✖ | 179/188 |
| bin_0.20 | bin | w40 | threshold_pol2 | 4342 ± 603 | 1.01 | 145.39 | ✔ | 101/90 |
| sig_gauss_cb | sig | w40 | threshold_pol2 | 4279 ± 640 | 1.02 | 145.39 | ✔ | 336/383 |
| sig_scale_fixed | sig | w40 | threshold_pol2 | 4348 ± 406 | 1.00 | 145.40 | ✔ | 336/384 |
| window_d0_36 | window | w36 | threshold_pol2 | 3999 ± 554 | 1.00 | 145.39 | ✔ | 348/383 |
| window_d0_55 | window | w55 | threshold_pol2 | 5187 ± 878 | 1.16 | 145.46 | ✔ | 373/383 |

| bin pT | N ± stat | skala σ | valid | tło |
|---|---|---|---|---|
| 3-5 | 644 ± 305 | 0.80 | ✔ | threshold_pol2 |
| 5-8 | 1353 ± 302 | 1.08 | ✔ | threshold_pol2 |
| 8-30 | 2174 ± 175 | 0.95 | ✔ | threshold_pol2 |

**N_D* = 4443 ± 612 (stat) ± 591 (syst)**; syst per kat.: {"bkg": 521, "range": 202, "bin": 100, "sig": 163}; z „czystym wykładnikiem” (Δχ²=+5, σ×1.17): ±1258. Mediana w40 4345, rozrzut 195 (10 wariantów). Suma per bin 4171.
Stare: fit_dm (pf_pairs, użyte w unblind) 5588 ± 586; fit_dm (V) 6607; refit V (invalid) 4472. Korekta μ: ×1.258.

**Uzupełnienie 2026-11-09 (etap ndstar, `PD_pairs_diff/ndstar/ndstar.md`):** dodatkowa kategoria syst „źródło kształtu MC” — ten sam fit z kształtem
z MC potoku par (σ 0,707, k 2,11, f 0,699; użyty w biegu PD → 4022) daje na `norm_all` 4097 (χ² identyczne) → **±339**.
Syst łącznie 591 ⊕ 339 = **681**; δN/N (stat ⊕ syst) = 20,6 %.
