# T-FIN — pary D–D po poprawce weta wspólnych śladów (PF), pełna lumi — ŚLEPO (2026-11-04 21:40 CET)

Źródło: `charm_work/PF2_lxplus/runs/lxplus/` (pf.dag 13:08→20:18; `compare_pf.json` sha256 0df5cce3886dbe0b…, `pairs/fit2d_blind.json` 8cd559cec09ca20c…). Centrum SS danych nieoglądane. Normalizacja `direct`, N_D* = 5588±586 (**stary `fit_dm`** — pf.dag biegł ze źródłem sprzed 2.0.2; dotyczy skali μ, nie Z ślepego; węzeł `unblind` używa 2.0.2/`dstar_refit`). Model fit2d: kategorie score + κ (plan B, T/U) — bez członu tilt+rho (V).

## 1. Weto wspólnych śladów (`pairs_counts.json`)
`track_id=d0z0`: par przed wetem 296 782, usuniętych **226 944 (76,5 %)**, w histogramach 56 657 (region ślepy).

## 2. Tło i czułość — ślepo
| kanał | SB×SB nominal → PF | usunięte | Z_A nom → PF | σ_μ nom → PF | dev/ndf nom → PF | corr(m1,m2) SB×SB same (mix) |
|---|---|---|---|---|---|---|
| D0D0 | 2 382 → 1 385 | 42 % | 1,38 → **1,90** | 0,79 → 0,57 | 0,79 → 0,65 | 0,073 (−0,010) |
| D⁺D0 | 17 744 → 7 263 | 59 % | 1,60 → **2,71** | 0,69 → 0,41 | 1,79 → 0,94 | −0,005 (0,002) |
| D⁺D⁺ | 92 077 → 14 005 | 85 % | 0,47 → **1,10** | 2,72 → 0,90 | 3,05 → 1,07 | 0,003 (0,003) |
| **łącznie (kwadratura)** | | | **2,17 → 3,49** | | | |

Oczekiwane S (pełna płaszczyzna): D0D0 125,7 / D⁺D0 204,7 / D⁺D⁺ 61,2 par; w centrum 58,2 / 90,5 / 25,8.
Odniesienia (stare weto): ½ lumi 1,65σ · wrzesień 2,63σ (toys 2,8) · pełna lumi tilt+rho (V) 2,27σ.
Więz κ partnera: σ=0,3 → 3,29σ; σ=0,5 → 3,19σ. Warianty WP/arbitrażu/tagu: 3,35–3,50σ.

## 3. Toys (50 jobów × 200, μ_true = 1)
| kanał | σ(pull) | bias/σ_μ | pokrycie 68 % | Z med [16–84 %] | valid |
|---|---|---|---|---|---|
| D0D0 | 1,012 | −0,03 | 0,693 | 1,83 [0,80–2,84] | 0,996 |
| D⁺D0 | 1,006 | +0,03 | 0,684 | 2,67 [1,63–3,68] | 0,987 |
| D⁺D⁺ | 1,025 | −0,08 | 0,681 | 1,06 [0,04–2,06] | 1,000 |

## 4. Systematyki po PF (U; `pairs/syst/syst_blind.md`)
| kanał | σ_stat | kształty | SF | normalizacja | κ/kategorie | suma | syst/stat |
|---|---|---|---|---|---|---|---|
| D0D0 | 0,573 | 0,092 | 0,177 | 0,240 | 0,066 | 0,319 | 0,56 |
| D⁺D0 | 0,411 | 0,078 | 0,132 | 0,190 | 0,123 | 0,274 | 0,67 |
| D⁺D⁺ | 0,898 | 0,092 | 0,201 | 0,280 | 0,083 | 0,367 | 0,41 |
| wspólne | 0,313 | 0,071 | 0,132 | 0,216 | 0,094 | 0,279 | 0,89 |

Z_A z systematykami kształtu/tła/κ: 3,44 → 3,30σ. SF i normalizacja zmieniają tylko skalę μ. **Brak**: tune Pythii (MPI/CR/Monash — joby `gen_syst` czekają). Dodatkowo z projektu B: model tła exp/pol2 (Y) ±3–4 % @WP70, SALAD (S8) ±5/10 %, SF v3 (N2/B2) ±8–13 % — mieszczą się w kolumnach „kształty”/„SF”.
Dyspersja plików: σ_jk/σ D0D0 1,16–1,31, D⁺D0 1,75–1,85, **D⁺D⁺ 2,8–3,2** (klastrowanie zostało mimo weta) → `mu_err` = max(Hesse·infl., profil) jak w preregistracji.

## 5. Mixing (S5 `nom_mix`, po wecie)
χ²/ndf SB×SB same vs mix: D0D0 0,79 · D⁺D0 0,94 · D⁺D⁺ 1,20; ramiona same/mix 0,988±0,016 / 0,974±0,008 / 0,976±0,008. Korelacja m1–m2 w D⁺D⁺ **0,26 → 0,003** — potwierdzenie PF na pełnej lumi; **tilt+rho (V) nie jest już potrzebne** jako nominał (zostaje cross-check).

## 6. Bramki (preregistracja T-PREP)
| # | kryterium | D0D0 | D⁺D0 | D⁺D⁺ | uwagi |
|---|---|---|---|---|---|
| G1 | σ(pull) 0,90–1,10 | ✔ 1,01 | ✔ 1,01 | ✔ 1,03 | |
| G2 | bias < 0,1 σ | ✔ 0,03 | ✔ 0,03 | ✔ 0,08 | |
| G3 | pokrycie 0,62–0,74 | ✔ 0,69 | ✔ 0,68 | ✔ 0,68 | |
| G4 | dev/ndf < 1,5 | ✔ 0,65 | ✔ 0,94 | ✔ 1,07 | |
| G5 | ρ po wecie ≈ MC (< 0,10) | ✔ 0,07 (1393 par, ±0,03) | ✔ −0,005 | ✔ 0,003 | |
| G6 | udział usuniętych w pasie | ✔ 42 % (35–50) | **✖ 59 % (25–40)** | ✔ 85 % (70–90) | patrz niżej |
| G7 | syst/stat < 1 (wspólne < 1,2) | ✔ 0,56 | ✔ 0,67 | ✔ 0,41 | wspólne 0,89; tune czeka |
| G8 | SF w zakresie, błąd < 10 % | ✔ v2 0,88 / 0,95 / 0,70 | ✔ | ✔ | v3: unc 0,07–0,10 (N2), ±0,05 @8–30 (B2) |
| G9 | N_D* dstar_refit vs fit_dm | ◐ | ◐ | ◐ | fit_dm 5588; dstar_refit (B2/N2) ≈ 4,7–5,4 k → Δ ≈ −3…−16 %; wiersz syst 10,5 % jest; unblind liczy 2.0.2 |
| G10 | mixing χ²/ndf < 2 | ✔ 0,79 | ✔ 0,94 | ✔ 1,20 | |
| G11 | Z_A PF ≥ nom − 0,1 | ✔ | ✔ | ✔ | |
| G12 | hashe zamrożone | ◐ | ◐ | ◐ | hash wyników zapisany; **src na lxplus w czasie biegu = AA+PF sprzed 2.0.1** (lokalnie już 2.0.2) — hash biegu do odtworzenia z EOS `pf_pairs/src` |

**Werdykt automatyczny + ręczny: NIE-GO jeszcze (1 punkt do wyjaśnienia).** G1–G5, G7, G8, G10, G11 PASS we wszystkich kanałach. **G6 D⁺D0: usunięto 59 % par SB×SB wobec oczekiwanych 25–40 %** (PF na ½ lumi: 33 % danych, 47 % MC). Możliwe przyczyny: (a) „nominal” to bieg głównego DAG na starszym kodzie (L bez S5/AA) — różnice poza wetem (weto D* dla D⁺, dedup, progi), więc „nominal→PF” nie izoluje weta; (b) w SB×SB po WP udział wspólnych śladów jest wyższy niż na wszystkich parach. Test rozstrzygający (ślepo, ~1 h, fit-only): ten sam `pf_pairs/src` z `pairs.track_id: pt_eta` → różnica = czyste weto. Jeśli po tym D⁺D0 mieści się w 25–45 % lub nadwyżka ma wytłumaczenie w (a) → GO.

## 7. Rekomendacja
1. Przed odsłonięciem: (i) węzeł `pf_pteta` (fit-only, `track_id: pt_eta`) → G6 per kanał z tego samego kodu; (ii) zapisać hash `pf_pairs/src` z EOS w checkliście; (iii) `unblind` i tak przelicza N_D* `dstar_refit` (2.0.2) → G9 zamknięte w tym samym biegu.
2. Po (i): **GO** — odsłonięcie 3 kanałów naraz, raport μ±σ, Z lokalne + łączne (Asimov ślepy 3,49σ; z syst 3,30σ), przedział toys [1,6–3,7] dla D⁺D0 itd.
3. tilt+rho (V) → cross-check, nie nominał. Tune Pythii dopisać do tabeli syst po `gen_syst`.
