# UBJK — jackknife μ po plikach danych (ODSŁONIĘTE, pełna lumi run 428770)

- grupy: events:200 (grupy zdarzeń: hash(run, event) mod 200); nogi: przebudowa cmd_pairs w scratch (src bieżący); syntetyczne: False
- grup (plików): 200; tryby: fixed; główny: **fixed**
- precyzja σ_jk ≈ ±5 % (1/√(2(n−1)), gauss.)

## tryb `fixed`

| | μ nom. | σ (Hesse/68 %) | σ_jk | σ_jk ważony | σ_jk/σ | σ_jk,w/σ | μ_J (bias-corr.) | Z nom. | Z bez pliku (min–max) | max |Δμ| (plik) |
|---|---|---|---|---|---|---|---|---|---|---|
| D0D0 | 0.546 | 0.705 | 0.901 | 0.901 | 1.28 | 1.28 | -0.073 | 0.84 | 0.00–1.41 | 0.551 (g192) |
| DpD0 | 0.767 | 0.589 | 0.680 | 0.680 | 1.16 | 1.16 | 0.812 | 1.51 | 0.61–1.67 | 0.424 (g146) |
| DpDp | -0.071 | 0.938 | 1.178 | 1.178 | 1.26 | 1.26 | 0.085 | 0.00 | 0.00–0.06 | 0.948 (g189) |
| combined | 0.567 | 0.401 | 0.454 | 0.454 | 1.13 | 1.13 | 0.378 | 1.57 | 0.84–1.87 | 0.243 (g146) |
| wmean | 0.535 | 0.407 | 0.488 | 0.488 | 1.20 | – | 0.499 | – | – | 0.223 (g146) |

Korelacje kanałów z jackknife (μ bez pliku k): D0D0–DpD0 -0.05, D0D0–DpDp +0.04, DpD0–DpDp -0.05

Wspólne μ (ciągłe minimum sumy profili, replika nominalna): 0.567, 68 % [0.18817204400010362, 0.9895507214901751], Z = 1.57; 05.11 (siatka 0,25): 0.50, Z = 1.56

σ_jk/σ > 1 ⇒ korelacje par w zdarzeniach (lub niejednorodność plików) poszerzają błąd ponad Poissona w komórkach;
jackknife ma rozdzielczość ±1/√(2(n−1)) ≈ ±5 % na σ (n = 200) — różnice poniżej ~2× tego nieistotne.
Tryb fixed = definicja `unblinded_jackknife` (dataleg/SB×SB z pełnej próbki); rederive = wszystkie wejścia z danych bez pliku k.
N_D* (normalizacja) i wejścia MC — skalowane ∝ zdarzeniom (fluktuacja N_D* w systematyce normalizacji, nie tutaj).
