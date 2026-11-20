# T-FIN — pary D–D po poprawce weta wspólnych śladów (PF) — ODSŁONIĘTE

Wygenerowano 2026-11-05 19:03 UTC z `/eos/user/k/klyczak/charm_q/pf_pairs` i `/eos/user/k/klyczak/charm_q/pf_pairs/unblind`. Odniesienia: ½ lumi (Mac) 1,65σ · wrzesień v1.0 2,63σ · pełna lumi tilt+rho (V, stare weto) 2,27σ.

## 1. Weto wspólnych śladów (pairs_counts.json)
track_id `d0z0`: par przed wetem 296782, usuniętych 226944 (76.5 %), w histogramach 56657 (region blind).

## 2. Tło i czułość — ślepo (compare_pf)
| kanał | SB×SB nominal | SB×SB PF | Δ | Z_A nominal | Z_A PF | σ_μ nominal | σ_μ PF | dev/ndf nom | dev/ndf PF | ρ PF | tilt PF |
|---|---|---|---|---|---|---|---|---|---|---|---|
| D0D0 | 2382 | 1385 | -42 % | 1.38 | 1.90 | 0.79 | 0.57 | 0.79 | 0.65 | 0.073 | — |
| D⁺D0 | 17744 | 7263 | -59 % | 1.60 | 2.71 | 0.69 | 0.41 | 1.79 | 0.94 | -0.005 | — |
| D⁺D⁺ | 92077 | 14005 | -85 % | 0.47 | 1.10 | 2.72 | 0.90 | 3.05 | 1.07 | 0.003 | — |
Z_A łączne (kwadratura): nominal 2.17 → PF **3.49**.

## 3. Toys (PF, ślepo, μ=1)
| kanał | σ(pull) | bias/σ_μ | pokrycie 68 % | n_toys | G1 | G2 | G3 |
|---|---|---|---|---|---|---|---|
| D0D0 | 1.01 | -0.03 | 0.69 | 10000 | ✔ | ✔ | ✔ |
| D⁺D0 | 1.01 | 0.03 | 0.68 | 10000 | ✔ | ✔ | ✔ |
| D⁺D⁺ | 1.02 | -0.08 | 0.68 | 10000 | ✔ | ✔ | ✔ |

## 4. Systematyki po PF (U + Y + S8 + N)
| kanał | σ_stat | kształt+tło+κ | SF | normalizacja | exp/pol2 (Y) | SALAD (S8) | suma syst | syst/stat | G7 |
|---|---|---|---|---|---|---|---|---|---|
| D0D0 | 0.57 | 0.11 | 0.18 | 0.24 | — | — | 0.32 | 0.56 | ✔ |
| D⁺D0 | 0.41 | 0.15 | 0.13 | 0.19 | — | — | 0.27 | 0.67 | ✔ |
| D⁺D⁺ | 0.90 | 0.12 | 0.20 | 0.28 | — | — | 0.37 | 0.41 | ✔ |
Źródła: `pairs/syst/*.json` (7 plików), `syst_bkg_model_Dplus.json` (Y), `bkg_shape_salad_*.json` (S8). Brakujące → „—” i G7 = ?.

## 5. Mixing (S5 nom_mix)
D0D0: χ²/ndf SB×SB same vs mix 0.79, corr same 0.073 / mix -0.010, ramiona same/mix 0.988 · D⁺D0: χ²/ndf SB×SB same vs mix 0.94, corr same -0.005 / mix 0.002, ramiona same/mix 0.974 · D⁺D⁺: χ²/ndf SB×SB same vs mix 1.20, corr same 0.003 / mix 0.003, ramiona same/mix 0.976

## 6. Checklist GO/NO-GO (automatyczna część; reszta ręcznie w unblinding_checklist.md)
| # | kryterium | D0D0 | D⁺D0 | D⁺D⁺ |
|---|---|---|---|---|
| G1 | toys σ(pull) 0,9–1,1 | ✔ | ✔ | ✔ |
| G2 | toys |bias| < 0,1σ | ✔ | ✔ | ✔ |
| G3 | pokrycie 68 % 0,62–0,74 | ✔ | ✔ | ✔ |
| G4 | dev/ndf ślepy < 1,5 | ✔ | ✔ | ✔ |
| G6 | udział usuniętych par SB×SB w pasie (per kanał) | ✔ | ✖ | ✔ |
| G7 | syst/stat < 1 | ✔ | ✔ | ✔ |
| G11 | Z_A PF ≥ Z_A nom − 0,1 | ✔ | ✔ | ✔ |
G6 źródło: SB×SB nominal→PF (NIEizolowane: miesza weto ze zmianą kodu); udział usuniętych SB×SB: D0D0 42 % (pas 30–65 %), D⁺D0 59 % (pas 21–50 %), D⁺D⁺ 85 % (pas 72–93 %)

Werdykt automatyczny: **NO-GO** (G1–G4, G6, G7, G11 wszystkie kanały; G5, G8–G10, G12 ręcznie).

## 7. Wynik odsłonięty
| kanał | μ̂ | 68 % | Z lok. | dev/ndf | μ bez tagu | κ=channel μ (Z) |
|---|---|---|---|---|---|---|
| D0D0 | 1.16 ± 0.55 | [0.64, 1.72] | 2.25 (exp 1.90) | 0.69 | 1.32 ± 0.76 | 1.17 (2.28) |
| D⁺D0 | 1.16 ± 0.42 | [0.75, 1.58] | 3.01 (exp 2.71) | 0.96 | 0.94 ± 0.53 | 1.08 (2.79) |
| D⁺D⁺ | 1.39 ± 0.93 | [0.49, 2.33] | 1.57 (exp 1.10) | 1.07 | 1.39 ± 0.93 | 1.37 (1.55) |

**Wspólne μ = 1.25 [0.90, 1.51], Z = 4.06** (stat); z systematykami profilowanymi: μ = 1.25, Z = 4.11.
Wynik vs ślepy Asimov: Z_A PF = 3.49; dane/MC μ i przedział jak wyżej (normalizacja: direct).

## 8. Porównanie z wcześniejszymi liczbami
| | Z | uwagi |
|---|---|---|
| ½ lumi (Mac, 23/47) | 1,65 | stare weto |
| wrzesień v1.0 (35 plików) | 2,63 (toys 2,8) | stare weto, inny klasyfikator |
| pełna lumi tilt+rho (V) ślepo | 2,27 | stare weto |
| **PF ślepo (Asimov)** | **3.49** | poprawne weto |
| **PF odsłonięte** | **4.06** | stat; z syst: 4.11 |
