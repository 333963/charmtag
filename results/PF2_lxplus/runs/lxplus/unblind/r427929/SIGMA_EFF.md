# Dolna granica σ_eff z μ(427929) — 19.11.2026 (projekt A)

Wynik roboczy do noty (§8, wariant B). Zależy od etapu `ubndsw` (szerokość piku D* → N_D* 427929). Skrypt: obliczenie inline (numpy/scipy, toy 4·10⁵, ziarno 1).

## Model (istniejący, bez nowych założeń fizycznych)
μ = (1 − f_DPS) + f_DPS · σ_eff^Py / σ_eff,  σ_eff^Py ≈ 28 mb (Pythia A3; notatki projektu, „Interpretacja”).
- f_DPS = udział DPS w parach Pythii: **0,40** (≈ 2·f_LS generatora 0,20; `MULT_charm/PROCEDURA.md` §5); wariant **0,54** (2 × f_LS 0,27 po poprawce znaku D⁰, `unblind/fls/FLS.md`); 0,20 jako skrajny dół.
- SPS w Pythii przyjęty jako dokładny; DPS ma w danych tę samą kinematykę/wydajność co DPS Pythii; całe odchylenie μ od 1 przypisane DPS.
- Kontrola mapy: μ = 1,5 → σ_eff 12,4 mb (f 0,40) / 14,5 mb (f 0,54) — w przedziale „6–16 mb” z notatek projektu.
- σ_eff,min = 28·f / (μ_up − 1 + f); jeśli μ_up ≤ 1 − f → brak granicy.

## Wejście
μ(427929) = 0,50 +0,41/−0,37 (stat) ± 0,19 (syst) ± 0,21 (N_D*) (prawda mb_nd; `truthnorm/TRUTHNORM.md`, `BUDZET_427929.md`). Do górnej granicy wchodzi błąd górny stat 0,41.
Wariant „szerokość piku z MC” μ ≈ 0,72, te same błędy względne.

## Wyniki (górna granica μ jednostronna → dolna granica σ_eff [mb] dla f = 0,40 / 0,54 / 0,20)

| μ | metoda błędów | CL | μ_up | σ_eff > |
|---|---|---|---|---|
| 0,50 | **syst i N_D* multiplikatywne (log-normal, rel. 0,566), stat gauss** — nominał | 95 % | **2,12** | **7,4** / 9,1 / 4,2 |
| 0,50 | jw. | 68 % | 0,82 | 50 / 42 / 226 |
| 0,50 | wszystko addytywne gauss (σ = 0,50) | 95 % | 1,32 | 15,6 / 17,6 / 10,8 |
| 0,50 | jw. | 68 % | 0,73 | 84 / 55 / brak |
| 0,72 | log-normal | 95 % | 3,05 | 4,6 / 5,8 / 2,5 |
| 0,72 | log-normal | 68 % | 1,19 | 19 / 21 / 14 |
| 0,72 | addytywne gauss (σ = 0,72) | 95 % | 1,90 | 8,6 / 10,5 / 5,1 |

Nominał = log-normal: syst (38,7 %) i N_D* (42 %) są względne (μ ∝ 1/N_D*), a addytywny gauss zaniża górną granicę. Wariant „syst względne liczone w granicy” z rozkładem gaussa rozbiega się (błąd względny × 1,645 ≈ 0,94) — odrzucony.

## Odczyt
- **σ_eff > 7 mb (95 % CL, f_DPS = 0,40)** przy μ = 0,50; przy μ = 0,72 (pik z MC) > 5 mb. Rozpiętość wg f_DPS i metody: 4–18 mb.
- Granica jest słaba: wyklucza tylko bardzo silny DPS (σ_eff ≲ 5–7 mb); nie rozróżnia typowych wartości 10–28 mb. Ogranicza ją błąd N_D* (±42 %) i stat 427929.
- Etap `ubndsw`: jeśli szerokość prawdziwa (μ = 0,50, mniejszy syst N_D*) → granica wzrośnie; jeśli pik z MC (μ ≈ 0,72) → spadnie do ~5 mb.

## Aktualizacja 19.11 wieczór — syst końcowy (decyzja: tune ↑ = fragmentacja ALICE)
Wejście końcowe: μ = 0,50 +0,41/−0,37 (stat) +0,28/−0,19 (syst; względnie +0,551/−0,374) ± 0,21 (N_D*, 0,414) — `truthnorm/TRUTHNORM.md`, `BUDZET_427929.md`.
Skrypt `seff_final.py` (toy 10⁶, ziarno 1): μ_true = (μ + stat·g₁)·exp(ln(1+syst↑)·g₂)·exp(ln(1+N_D*)·g₃). Uwaga: toy z wersji wyżej (μ_up 2,12) nie dał się
odtworzyć (kod inline nie zapisany); ten jawny toy przy starych wejściach daje μ_up 1,63 (σ_eff > 10,9 mb), czyli jest mniej zachowawczy.

| wariant | μ_up 95 % (jawny toy) | σ_eff > [mb] f = 0,40 / 0,54 / 0,20 | μ_up 68 % |
|---|---|---|---|
| 0,50, syst stary sym. 0,387 (kontrola metody) | 1,63 | 10,9 / 13,0 / 6,8 | 0,70 |
| **0,50, syst końcowy +0,551** | **1,79** | **9,4** / 11,4 / 5,6 | 0,70 |
| 0,72 (skala σ z MC) | 2,58 | 5,7 / 7,1 / 3,1 | 1,01 |
| 0,70 (fragmentacja ALICE, μ/0,71) | 2,52 | 5,8 / 7,3 / 3,3 | 0,99 |

Zapis do noty (zachowawczo): poprzedni toy przeskalowany przyrostem μ_up z jawnego toya (×1,14 na μ_up − μ) → μ_up ≈ 2,35, **σ_eff > ~6 mb (95 % CL, f_DPS = 0,40)**;
rozpiętość wg metody i f_DPS 4–11 mb; warianty interpretacyjne (pik z MC, fragmentacja ALICE) > ~5 mb (jawny toy) / ~4,7 mb (FRAG.md, stary toy).
