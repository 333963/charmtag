# SYST_TRIG — wiersz „trigger / selekcja zdarzeń” (reguły: TRIG_A/SYST_TRIG.md, prerejestracja 19.11)

Status: **efekt (zmiana nominału tylko osobną decyzją)**


| pomiar | μ (wejście) | d = (μ_nom − μ_czyste)/μ_nom | σ_d | Z |
|---|---|---|---|---|
| d₁ UB427 (427929 vs 428770; f_jk = 1.13, zmierzone (UB427JK, ±5 %)) | 0.53 vs 1.50 | +0.648 | 0.461 | +1.40 |
| d₂ UBTRIG (J12only vs reszta); r_J12 = 0.769 | 1.63 ± 0.46 vs -0.53 ± 0.51 (μ_nom 1.19) | +1.394 | 0.443 | +3.15 |
| d₂ kontrola (J12only vs MBTS_noj20 ⊕ other); r_J12 = 0.769 | 1.63 ± 0.46 vs -0.18 ± 0.85 (μ_nom 1.19) | +1.166 | 0.625 | +1.87 |
| **d (średnia)** | | **+1.036** | 0.319 | **+3.24** |

Emulacja J12 (TRIGCAL): κ = 1.24 → d_emu = 0.156 (zakres 0.101–0.230) przy r_J12 = 0.769; zgodność danych z emulacją: +2.76σ

Wiersz (μ = 1.5): `{"selekcja_zdarzen": -1.553835963860555, "jednostronny": true, "zastepuje_wag": true, "wariant_mu_corr": [-0.053835963860555, 0.4791607745991543]}`
