# VER_budget — budżet μ po normalizacji z prawdy (gen_norm v2.0), 19.11.2026

Źródła: `GEN_mcstat/out/TRUTH_NPAIR.md` (bramki, liczby), `truth_npair.json`, `truth_mbnd.json`; konwencja jak `RAPORT.md` (16.11): wspólne Σ w·δ,
w = 0,298 / 0,580 / 0,122, × 1,495; „pozostałe” = √(0,284² − 0,188²) = 0,213 (SF, kształty, tune ±0,17, …) bez zmian; κ_cls 0,042 bez zmian.
Nominał NIE zmieniony (UNB_WYNIK.md bez zmian). Zmiana trybu normalizacji = osobna decyzja (PROCEDURA §5.4); bez fitu par danych (reguła 0.2).

## Tabela przed / po (syst. wspólne, bez N_D*)
| wiersz | 16.11 (direct) | A3 truth (prerejestracja §1) | mb_nd truth (odstępstwo, rekomendowane) |
|---|---|---|---|
| pozostałe (SF, kształty, tune, …) | ±0,32 | ±0,32 | ±0,32 |
| κ_cls | ±0,06 | ±0,06 | ±0,06 |
| stat. MC par | ±0,46 (n_eff 8,7/13,3/7,5) | ±0,10 | ±0,18 (n_eff 39/72/34) |
| model wag w2↔w1 | ±0,24 (proxy singli) | ±0,24 | ±0,20 |
| transfer/binning map | — | ±0,24 | — |
| faktoryzacja par, k_sp, proxy, próbka w2 | — | ±0,07 | ±0,06 |
| **syst razem** | **±0,61** | ±0,49 | **±0,42** |
| centralne n_pair vs direct (truth/direct, wspólne) | 1 | 1,30 (0,90 / 1,29 / 2,35) | 1,17 (1,03 / 1,27 / 1,06) |
| trigger / selekcja zdarzeń | osobno — **czeka na UBTRIGNULL / ub427jk** | | |

## Werdykt
1. **Bramki B1–B6 przechodzą** (B4 dla D⁺ blisko progu: 1,127 przy progu 1,142). Formalnie tryb truth A3 może zastąpić direct (§0.1).
2. **A3 nie nadaje się na nominał:** B7 FAIL (D*/zdarzenie mb_nd/A3 = 1,36) i na poziomie prawdy D⁺/D* = 5,15 (A3) vs 3,63 (mc23 mb_nd), ×1,42.
   Pary z D⁺ w A3 są przez to zawyżone (D⁺D⁰ ×1,29, D⁺D⁺ ×2,35 vs direct) — to zmiana definicji odniesienia μ, nie poprawa precyzji.
3. **Rekomendacja (odstępstwo od prerejestracji, do decyzji):** normalizacja z prawdy **mb_nd** × mapy direct — ta sama próbka mc23 co nominał
   (definicja μ bez zmian), bez proxy krotności, n_eff par ×5. Syst ±0,61 → **±0,42**; centralne n_pair rośnie (D⁺D⁰ +27 %), więc μ zmalałoby
   — przelicznik liniowy ~÷1,17 (to NIE jest fit; fit `--norm truth` z `truth_norm` z mb_nd dopiero po decyzji).
4. Bez zmiany trybu (direct) zapis zostaje: **μ = 1,50 +0,48/−0,46 (stat) ± 0,61 (syst) ± 0,33 (N_D*)**; prawda mb_nd potwierdza wiersz modelu wag
   (0,135 wspólne vs 0,163 proxy) i nie wnosi korekty poza statystyką.

## Co niepewne
- Jackknife prawdy mb_nd (26 grup po 3 pliki) zdominowany przez kilka ciężkich wag w2 (w do ~43); D⁺D⁺ model wag 39 % (w1 źle domyka ogon krotności).
- DpD0 truth_mb/direct = 1,27 (1σ statystyki direct) — który jest bliżej prawdy, rozstrzyga statystyka, nie test.
- Źródło różnicy A3 ↔ mc23 (D⁺/D*, D*/zdarzenie) niezbadane: wersja Pythii / dekaje EvtGen / filtr próbki mb_nd „Dpt3to150” (nie sprawdzone).
- Trigger (UBTRIG: d = 1,04 ± 0,32) może dominować nad wszystkim powyżej — wiersz osobny, czeka na test zerowy i jackknife 427929.
