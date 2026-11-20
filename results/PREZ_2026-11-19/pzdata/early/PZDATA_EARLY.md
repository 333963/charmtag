# PZDATA — wynik priorytetowy (20.11.2026 ~03:45): 427929, D⁰ i D* (noga D⁰), sieć MLP high-level na pasmach danych (sb_net_hi) vs nominalny XGB H

Źródło: etap ORCH `pzdata` (Q_lxplus/pz_data). Modele sb_net_hi wytrenowane na lxplus (tło = pasma boczne danych 428770 jak nominał, 5 foldów;
AUC na prawdzie MC w oknie 0,955 OOF). Zastosowanie (`pzd_apply.py`) uruchomione LOKALNIE na pełnych danych 427929 (458 plików h5, klasa MB, jak UB427)
i pełnym fd_nd v1.1 (24 pliki), bo węzły apply na lxplus czekały w kolejce (klaster zapchany); te same skrypty i modele. Raport `pzd_report.py`.
Progi: ten sam ε_S na sygnale MC fd per bin pT. Fit D⁰: model nominału. D*: Δm z nogą D⁰ dopasowaną do D* (jak charmtag dstar_control), `massfit.fit_dm`.

| cięcie na nodze D⁰ | wariant | N_D* ± σ | σ/N | N_D⁰ (pT 3–30) ± σ | σ/N D⁰ |
|---|---|---|---|---|---|
| brak | — | 2151 ± 332 | 0,154 | 14310 ± 1844 | 0,129 |
| WP90 | XGB H (nominał) | 1405 ± 116 | 0,083 | 12239 ± 868 | 0,071 |
| WP90 | sieć sb_net_hi | 1470 ± 112 | 0,076 | 8121 ± 721 | 0,089 |
| WP70 | XGB H (nominał) | 1138 ± 62 | 0,055 | 6014 ± 375 | 0,062 |
| WP70 | sieć sb_net_hi | 1090 ± 60 | 0,055 | 6380 ± 369 | 0,058 |

- Kontrola: N_D* bez cięcia 2151 ± 332 zgodne z NDS 427929 (2046 ± 485).
- D*: sieć = XGB w granicach ~8 % σ/N (WP90 0,92×, WP70 1,00×) — brak zysku dla σ(N_D*).
- D⁰ WP70: sieć 0,93× σ/N; WP90 fit przy S/B ≈ 0,03 niestabilny (N nominału 12239 vs sieć 8121 przy tym samym ε_S MC; ε_data/ε_MC 0,95 vs 0,63) — nie interpretować.
- Rzeźbienie tła danych (R pasmo wewn./zewn.): nom 1,043 / 1,058, sieć 1,035 / 1,054 (WP90 / WP70) — jednakowe.
- Uwaga: σ(N_D*) w budżecie wyniku głównego (±41 %) to głównie syst. kształtu/szerokości piku D* (NDS), nie stat. — cięcie na nodze D⁰ go nie usuwa; ten test mierzy tylko stat. precyzję.
Pozostałe warianty (mc_*, sb_net_lo), D⁺ i 428770 — pełny raport etapu (`../PZDATA.md` po zakończeniu DAG).
