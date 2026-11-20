# TRUTH_NPAIR — tryb truth n_pair (PROCEDURA §1–§4), 19.11.2026

Wejścia: `runs/lxplus_v20` (gen_norm ROOT v2.0, wariant 78 plików §3b, akceptacja PASS po przetagowaniu nóg 19.11), prawda A3 500 plików
(353,7 M zdarzeń ND), warianty tune po 100 plików. Kod: `truth_npair.py` → `truth_npair.json`, `truth_norm_v20.json` (format `charmtag.pairs.truthnorm`);
wariant informacyjny `truth_mbnd.py` → `truth_mbnd.json`. Bez danych par, bez μ (reguła 0.2). Test na 5 plikach: `truth_npair_test5.json`.

## Bramki (§3)
| | wynik | |
|---|---|---|
| B1 (§3a/§3b) | `acceptance.json` verdict PASS (12/12, w tym B1a_retag_bijective) | PASS |
| B2 | `ref_match` true (1344 / 9534 / 8665 / 22 / 41 / 15), `dstar_scan_equals_mc_rates` true | PASS |
| B3 | unmatched D0 0,62 %, D⁺ 0,11 %, D* 0 %; k_sp − 1 = 0,0063 / 0,0011 / 0 | PASS |
| B4 | split-half singli: R_D0 = 1,049 (próg max(0,05; 0,138)), R_D⁺ = 1,127 (próg 0,142) | PASS (D⁺ blisko progu) |
| B5 | split-half par: R = 1,04 / 1,31 / 1,19 (2σ = 0,69 / 0,56 / 0,74); R̄ = 1,218, 2σ̄ = 0,396 | PASS |
| B6 | c_w,truth / c_w,h5 (ε_trk 0,85): D0 0,939 / 0,904 = 1,040; D⁺ 0,781 / 0,769 = 1,015 | PASS |
| B7 (info) | prawdziwe D* (mode 1, pT > 3, \|η\| < 2,5) / zdarzenie: mb_nd 5,42·10⁻⁴ vs A3 ND 3,99·10⁻⁴ → ×1,36 | **FAIL (nie blokuje)** |

## Wynik A3 (§1)
| | D0D0 | D⁺D0 | D⁺D⁺ |
|---|---|---|---|
| [pary/D*]_truth (w2) | 0,02019 | 0,04715 | 0,02583 |
| direct mb_nd (h5, w2) | 0,02236 | 0,03657 | 0,01099 |
| **truth / direct** | 0,90 | **1,29** | **2,35** |
| stat. truth (jackknife 50 grup ⊕ bootstrap map) | 8,8 % | 6,7 % | 8,6 % |
| prawdziwe pary w akceptacji | 4846 | 12 906 | 8993 |

## Rozbieżność A3 ↔ mc23 mb_nd (poziom prawdy, bez map i wag; 60 plików A3 vs 78 mb_nd)
| na D* (mode 1, pT > 2,5, \|η\| < 2,7) | D0 | D⁺ | D0D0 | D⁺D0 | D⁺D⁺ |
|---|---|---|---|---|---|
| mb_nd (mc23) | 3,77 | 3,63 | 0,032 | 0,065 | 0,031 |
| A3 (pythia_dd) | 3,80 | **5,15** | 0,024 | 0,068 | **0,046** |
Udział mode = 1 wśród D⁺ jest ten sam (0,099), więc różnica jest w produkcji D⁺ względem D* (×1,42; por. VER (C): D⁺/D* = 1,33 w combined_*,
domknięcie fd→mb ×1,32). D⁺D0 ∝ D⁺/D*, D⁺D⁺ ∝ (D⁺/D*)² → truth/direct 1,29 / 2,35 to głównie ta różnica, nie statystyka.
**Wniosek:** tryb truth z A3 zmienia definicję odniesienia μ (dane/Pythia A3 standalone zamiast dane/mc23 mb_nd). Bramki B1–B6 przechodzą,
ale B7 i ten test pokazują, że A3 nie jest tą samą „Pythią” co nominał → nie rekomenduję A3 jako nominału bez decyzji o definicji odniesienia.

## Wariant informacyjny (ODSTĘPSTWO od §1): prawda z mb_nd × mapy direct
Te same zdarzenia mc23, reco n_trk i lead_pt (bez proxy), w2; 446 / 909 / 443 prawdziwych par (n_eff po w2·μμ 39 / 72 / 34 vs direct 8,7 / 13,3 / 7,5).
| | D0D0 | D⁺D0 | D⁺D⁺ | wspólne |
|---|---|---|---|---|
| truth_mb / direct | 1,03 | 1,27 | 1,06 | 1,17 |
| stat. (jackknife 26 grup plików ⊕ bootstrap map) | 16,7 % | 15,8 % | 16,5 % | 0,117 |
| model wag w2↔w1 | 9,9 % | 9,9 % | 39 % | 0,135 |
| próbka danych w2 | 0,5 % | 0,3 % | 0,2 % | 0,003 |

## Wiersze syst §4 (wspólne Σ w·δ, × 1,495 na μ)
| wiersz | A3 truth (prerejestracja) | mb_nd truth (odstępstwo) |
|---|---|---|
| stat. MC (prawda) | 0,064 → ±0,10 | 0,117 → ±0,18 |
| model wag w2↔w1 | 0,162 → ±0,24 | 0,135 → ±0,20 |
| próbka danych w2 | 0,002 → ±0,003 | 0,003 → ±0,005 |
| proxy n_ch→n_trk | 0,011 → ±0,02 | — (reco n_trk) |
| proxy lead_pt | 0,018 → ±0,03 | — (reco lead_pt) |
| transfer/binning map | 0,162 → ±0,24 | — (mapy z tej samej próbki; stat. map w bootstrapie) |
| faktoryzacja par | 0,039 → ±0,06 | 0,039 → ±0,06 |
| k_sp | 0,004 → ±0,006 | 0,004 → ±0,006 |
| tune (info; wiersz ±0,17 bez zmian) | max wariantów 18 / 12 / 15 % | — |
