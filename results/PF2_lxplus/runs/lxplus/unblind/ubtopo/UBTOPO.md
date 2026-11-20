# UB_TOPO — kształt `dy` sygnału par D–D vs Pythia/DPS (odsłonięte, 428770)

Biny ['[0,1)', '[1,5]']; hipotezy z generatora Pythii A3 (ważone ε nóg: `w`), rodzina h(λ) = λ·h_DPS + (1−λ)·h_Pythia.
λ = 0 — kształt Pythii (nominał MC), λ = 1 — czysty DPS (pary niezależne); „nadwyżka = DPS” ⇒ λ = (μ−1)/μ = 0.33 (μ = 1.5).

| kanał | h_Pythia | h_DPS | λ̂ | 68 % | Z(λ=0, Pythia) | Z(λ=1, DPS) | Z(λ = nadwyżka-DPS) |
|---|---|---|---|---|---|---|---|
| D0D0 | 0.601/0.399 | 0.406/0.594 | -0.60 | [-1.50; 0.72] | 0.47 | 1.16 | 0.74 |
| DpD0 | 0.601/0.399 | 0.422/0.578 | -1.50 (krawędź) | [-1.50; -0.95] | 2.01 | 2.71 | 2.28 |
| DpDp | 0.617/0.383 | 0.444/0.556 | 0.18 | [-1.40; 2.50] | 0.12 | 0.39 | 0.08 |
| łącznie |  |  | -1.11 | [-1.50; -0.51] | 1.71 | 2.80 | 2.13 |

Interpretacja (ustalona przed wynikiem, PROGNOZA.md): σ_λ ≈ 0,8 → test rozróżnia najwyżej czysty DPS od Pythii (~1,2σ);
pytania „nadwyżka DPS czy dżetowa” ten test nie rozstrzyga (Δλ ≈ 0,3 ≪ σ_λ; obciążenie dżetowe daje kształt SPS-podobny, λ ≤ 0).

Prognoza (PROGNOZA.md): σ_λ ≈ 0,85; czysty DPS vs Pythia ~1,2σ; nadwyżka-DPS (λ=0,33) vs Pythia ~0,5–0,75σ (ref/forecast.json)

Syst. wag ε (hipotezy bez wag ε): λ̂ = -1.15 [-1.5, -0.65].
