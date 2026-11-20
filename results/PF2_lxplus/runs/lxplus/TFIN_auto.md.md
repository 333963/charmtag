# T-FIN — pary D–D po poprawce weta wspólnych śladów (PF) — ŚLEPO

Wygenerowano 2026-11-04 20:33 UTC z `PF2_lxplus/runs/lxplus`. Odniesienia: ½ lumi (Mac) 1,65σ · wrzesień v1.0 2,63σ · pełna lumi tilt+rho (V, stare weto) 2,27σ.

## 1. Weto wspólnych śladów (pairs_counts.json)
track_id `d0z0`: par przed wetem 296782, usuniętych 226944 (76.5 %), w histogramach 56657 (region blind).

## 2. Tło i czułość — ślepo (compare_pf)
| kanał | SB×SB nominal | SB×SB PF | Δ | Z_A nominal | Z_A PF | σ_μ nominal | σ_μ PF | dev/ndf nom | dev/ndf PF | ρ PF | tilt PF |
|---|---|---|---|---|---|---|---|---|---|---|---|
| D0D0 | — | — | — | 1.38 | 1.90 | 0.79 | 0.57 | 0.79 | 0.65 | — | — |
| D⁺D0 | — | — | — | 1.60 | 2.71 | 0.69 | 0.41 | 1.79 | 0.94 | — | — |
| D⁺D⁺ | — | — | — | 0.47 | 1.10 | 2.72 | 0.90 | 3.05 | 1.07 | — | — |
Z_A łączne (kwadratura): nominal 2.17 → PF **3.49**.

## 3. Toys (PF, ślepo, μ=1)
| kanał | σ(pull) | bias/σ_μ | pokrycie 68 % | n_toys | G1 | G2 | G3 |
|---|---|---|---|---|---|---|---|
| D0D0 | — | -0.03 | 0.69 | 10000 | — | ✔ | ✔ |
| D⁺D0 | — | 0.03 | 0.68 | 10000 | — | ✔ | ✔ |
| D⁺D⁺ | — | -0.08 | 0.68 | 10000 | — | ✔ | ✔ |

## 4. Systematyki po PF (U + Y + S8 + N)
| kanał | σ_stat | kształt+tło+κ | SF | normalizacja | exp/pol2 (Y) | SALAD (S8) | suma syst | syst/stat | G7 |
|---|---|---|---|---|---|---|---|---|---|
| D0D0 | 0.57 | — | — | — | — | — | — | — | — |
| D⁺D0 | 0.41 | — | — | — | — | — | — | — | — |
| D⁺D⁺ | 0.90 | — | — | — | — | — | — | — | — |
Źródła: `pairs/syst/*.json` (7 plików), `syst_bkg_model_Dplus.json` (Y), `bkg_shape_salad_*.json` (S8). Brakujące → „—” i G7 = ?.

## 5. Mixing (S5 nom_mix)
D0D0: χ²/ndf same vs mix —, corr same — / mix — · D⁺D0: χ²/ndf same vs mix —, corr same — / mix — · D⁺D⁺: χ²/ndf same vs mix —, corr same — / mix —

## 6. Checklist GO/NO-GO (automatyczna część; reszta ręcznie w unblinding_checklist.md)
| # | kryterium | D0D0 | D⁺D0 | D⁺D⁺ |
|---|---|---|---|---|
| G1 | toys σ(pull) 0,9–1,1 | — | — | — |
| G2 | toys |bias| < 0,1σ | ✔ | ✔ | ✔ |
| G3 | pokrycie 68 % 0,62–0,74 | ✔ | ✔ | ✔ |
| G4 | dev/ndf ślepy < 1,5 | ✔ | ✔ | ✔ |
| G6 | udział wspólnych śladów w pasie (łączny) | ✖ | ✖ | ✔ |
| G7 | syst/stat < 1 | — | — | — |
| G11 | Z_A PF ≥ Z_A nom − 0,1 | ✔ | ✔ | ✔ |
Werdykt automatyczny: **NIEKOMPLETNE** (tylko G1–G4, G6, G11; G5/G7–G10/G12 ręcznie).

## 7. Wynik odsłonięty
— (nie odsłonięte; `submit_unblind.sh submit` po decyzji)

## 8. Porównanie z wcześniejszymi liczbami
| | Z | uwagi |
|---|---|---|
| ½ lumi (Mac, 23/47) | 1,65 | stare weto |
| wrzesień v1.0 (35 plików) | 2,63 (toys 2,8) | stare weto, inny klasyfikator |
| pełna lumi tilt+rho (V) ślepo | 2,27 | stare weto |
| **PF ślepo (Asimov)** | **3.49** | poprawne weto |
| **PF odsłonięte** | **—** | — |
