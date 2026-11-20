# UBJK — jackknife μ po plikach danych (ODSŁONIĘTE, pełna lumi run 428770)

- grupy: files (h5 (charmtag.syst.files.file_table)); nogi: legs_data.parquet z odsłonięcia (05.11); syntetyczne: False
- grup (plików): 39; tryby: fixed, rederive; główny: **rederive**
- precyzja σ_jk ≈ ±11 % (1/√(2(n−1)), gauss.)
- kontrola: `charmtag.syst.files.unblinded_jackknife` vs tryb fixed: max |Δσ_jk| = 3.42e-09

## tryb `fixed`

| | μ nom. | σ (Hesse/68 %) | σ_jk | σ_jk ważony | σ_jk/σ | σ_jk,w/σ | μ_J (bias-corr.) | Z nom. | Z bez pliku (min–max) | max |Δμ| (plik) |
|---|---|---|---|---|---|---|---|---|---|---|
| D0D0 | 1.158 | 0.552 | 0.530 | 0.550 | 0.96 | 1.00 | 1.168 | 2.25 | 1.73–2.52 | 0.253 (000027) |
| DpD0 | 1.165 | 0.424 | 0.459 | 0.451 | 1.08 | 1.06 | 1.131 | 3.01 | 2.35–3.26 | 0.233 (000036) |
| DpDp | 1.386 | 0.925 | 1.229 | 1.238 | 1.33 | 1.34 | 1.471 | 1.57 | 0.98–1.94 | 0.526 (000032) |
| combined | 1.189 | 0.317 | 0.318 | 0.316 | 1.00 | 1.00 | 1.190 | 4.07 | 3.54–4.25 | 0.132 (000020) |
| wmean | 1.188 | 0.316 | 0.321 | 0.320 | 1.02 | – | 1.229 | – | – | 0.136 (000020) |

Korelacje kanałów z jackknife (μ bez pliku k): D0D0–DpD0 +0.05, D0D0–DpDp -0.01, DpD0–DpDp -0.25

## tryb `rederive`

| | μ nom. | σ (Hesse/68 %) | σ_jk | σ_jk ważony | σ_jk/σ | σ_jk,w/σ | μ_J (bias-corr.) | Z nom. | Z bez pliku (min–max) | max |Δμ| (plik) |
|---|---|---|---|---|---|---|---|---|---|---|
| D0D0 | 1.158 | 0.552 | 0.524 | 0.545 | 0.95 | 0.99 | 1.166 | 2.25 | 1.74–2.54 | 0.247 (000027) |
| DpD0 | 1.165 | 0.424 | 0.463 | 0.454 | 1.09 | 1.07 | 1.130 | 3.01 | 2.33–3.26 | 0.236 (000036) |
| DpDp | 1.386 | 0.925 | 1.229 | 1.239 | 1.33 | 1.34 | 1.465 | 1.57 | 0.98–1.94 | 0.527 (000032) |
| combined | 1.189 | 0.317 | 0.317 | 0.315 | 1.00 | 0.99 | 1.189 | 4.07 | 3.53–4.25 | 0.137 (000020) |
| wmean | 1.188 | 0.316 | 0.321 | 0.319 | 1.01 | – | 1.226 | – | – | 0.141 (000020) |

Korelacje kanałów z jackknife (μ bez pliku k): D0D0–DpD0 +0.04, D0D0–DpDp -0.02, DpD0–DpDp -0.26

Wspólne μ (ciągłe minimum sumy profili, replika nominalna): 1.189, 68 % [0.8759847333213475, 1.5091999367836362], Z = 4.07; 05.11 (siatka 0,25): 1.25, Z = 4.06

σ_jk/σ > 1 ⇒ korelacje par w zdarzeniach (lub niejednorodność plików) poszerzają błąd ponad Poissona w komórkach;
jackknife ma rozdzielczość ±1/√(2(n−1)) ≈ ±11 % na σ (n = 39) — różnice poniżej ~2× tego nieistotne.
Tryb fixed = definicja `unblinded_jackknife` (dataleg/SB×SB z pełnej próbki); rederive = wszystkie wejścia z danych bez pliku k.
N_D* (normalizacja) i wejścia MC — skalowane ∝ zdarzeniom (fluktuacja N_D* w systematyce normalizacji, nie tutaj).
