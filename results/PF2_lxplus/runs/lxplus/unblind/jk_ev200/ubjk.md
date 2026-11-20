# UBJK — jackknife μ po plikach danych (ODSŁONIĘTE, pełna lumi run 428770)

- grupy: events:200 (grupy zdarzeń: hash(run, event) mod 200); nogi: legs_data.parquet z odsłonięcia (05.11); syntetyczne: False
- grup (plików): 200; tryby: fixed; główny: **fixed**
- precyzja σ_jk ≈ ±5 % (1/√(2(n−1)), gauss.)

## tryb `fixed`

| | μ nom. | σ (Hesse/68 %) | σ_jk | σ_jk ważony | σ_jk/σ | σ_jk,w/σ | μ_J (bias-corr.) | Z nom. | Z bez pliku (min–max) | max |Δμ| (plik) |
|---|---|---|---|---|---|---|---|---|---|---|
| D0D0 | 1.158 | 0.552 | 0.618 | 0.618 | 1.12 | 1.12 | 1.181 | 2.25 | 1.59–2.40 | 0.317 (g064) |
| DpD0 | 1.165 | 0.424 | 0.467 | 0.467 | 1.10 | 1.10 | 1.162 | 3.01 | 2.51–3.19 | 0.158 (g064) |
| DpDp | 1.386 | 0.925 | 1.122 | 1.122 | 1.21 | 1.21 | 1.478 | 1.57 | 0.95–1.73 | 0.559 (g079) |
| combined | 1.189 | 0.317 | 0.374 | 0.374 | 1.18 | 1.18 | 1.283 | 4.07 | 3.33–4.20 | 0.188 (g064) |
| wmean | 1.188 | 0.316 | 0.378 | 0.378 | 1.20 | – | 1.242 | – | – | 0.190 (g064) |

Korelacje kanałów z jackknife (μ bez pliku k): D0D0–DpD0 +0.16, D0D0–DpDp -0.00, DpD0–DpDp -0.01

Wspólne μ (ciągłe minimum sumy profili, replika nominalna): 1.189, 68 % [0.8759847333213475, 1.5091999367836362], Z = 4.07; 05.11 (siatka 0,25): 1.25, Z = 4.06

σ_jk/σ > 1 ⇒ korelacje par w zdarzeniach (lub niejednorodność plików) poszerzają błąd ponad Poissona w komórkach;
jackknife ma rozdzielczość ±1/√(2(n−1)) ≈ ±5 % na σ (n = 200) — różnice poniżej ~2× tego nieistotne.
Tryb fixed = definicja `unblinded_jackknife` (dataleg/SB×SB z pełnej próbki); rederive = wszystkie wejścia z danych bez pliku k.
N_D* (normalizacja) i wejścia MC — skalowane ∝ zdarzeniom (fluktuacja N_D* w systematyce normalizacji, nie tutaj).
