# Budżet systematyk μ par D–D z tune Pythii (MPI, CR z CRqcd, Monash)

Skrypt `TUNE_pythia/budget_with_tune.py` (z `charm_work/`). Wejścia: `PF2_lxplus/runs/lxplus/pairs/syst/syst_blind.json`, `TUNE_pythia/out/tune_syst.json`, `truth_norm.json` z `../runs/lxplus/truth` (ścieżka `norm.collect`).
Wartości względne (Δμ/μ, tj. Δμ przy μ=1). Nominalny fit, Z i μ centralne bez zmian.

## Reguły (z kodu `charmtag.syst`)
- Normalizacja: `norm.py` l. 97–110 — δ grupy tune = max po wariantach |[par/D*]_war/[par/D*]_A3 − 1|; l. 119 — **Δμ/μ = |δ|** (1. rząd). `dmu_rel(δ) = 1/(1+δ) − 1` (l. 26–28) służy tylko do pola asymetrii (l. 120), nie do tabeli.
- „wspólne”: `report.py` l. 40–41, 76–77 — Σ w·Δμ / Σ w, w = 1/σ²_stat (Asimov ślepy: D⁰D⁰ 0.573, D⁺D⁰ 0.411, D⁺D⁺ 0.898), brak w kanale = 0.
- Sumy: l. 86–93 — kwadratura wierszy (pomijane: pending oraz `kappa_partner*`); suma „wspólne” = kwadratura wierszy „wspólne”; l. 112–113 — „wszystko” = kwadratura grup.

## Kontrola: odtworzenie starych sum
Max |Δ| „wspólne” w wierszach: 0.0e+00; max |Δ| sum grup/„wszystko”: 4.8e-04 (próg 0.002).

## Wiersze tune (Normalizacja)

| grupa | warianty | Δμ D⁰D⁰ | Δμ D⁺D⁰ | Δμ D⁺D⁺ | Δμ wspólne | asym. wspólne (↑ / ↓) |
|---|---|---|---|---|---|---|
| tune_MPI | MPIup/MPIdown | 0.040 | 0.041 | 0.037 | 0.040 | +0.033 / -0.039 |
| tune_CR | CRup/CRdown/CRqcd | 0.112 | 0.109 | 0.096 | 0.108 | +0.121 / -0.040 |
| tune_Monash | Monash | 0.020 | 0.019 | 0.017 | 0.019 | +0.000 / -0.019 |
| **tune razem (kwadr.)** | | 0.121 | 0.118 | 0.104 | 0.117 | |

Asymetria: Δμ/μ = 1/(1+δ) − 1. CRqcd obniża oczekiwaną liczbę par (δ ≈ −10…−11 %), więc przesuwa μ **tylko w górę** (CR: +0.121 w górę vs -0.040 w dół, wspólne). Tabela (zgodnie z kodem) używa symetrycznego |δ|.

## Budżet przed / tune / po (względne)

| suma | D⁰D⁰ przed | tune | po | D⁺D⁰ przed | tune | po | D⁺D⁺ przed | tune | po | wspólne przed | tune | po |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Normalizacja | 0.240 | 0.121 | **0.269** | 0.190 | 0.118 | **0.224** | 0.280 | 0.104 | **0.299** | 0.216 | 0.117 | **0.245** |
| wszystko | 0.319 | 0.121 | **0.341** | 0.274 | 0.118 | **0.298** | 0.367 | 0.104 | **0.381** | 0.279 | 0.117 | **0.302** |
| Normalizacja bez N_D* | 0.216 | 0.121 | **0.248** | 0.159 | 0.118 | **0.198** | 0.260 | 0.104 | **0.280** | 0.188 | 0.117 | **0.222** |
| wszystko bez N_D* | 0.301 | 0.121 | **0.325** | 0.253 | 0.118 | **0.279** | 0.351 | 0.104 | **0.366** | 0.258 | 0.117 | **0.284** |

Wiersz N_D* z tabeli ślepej (`ndstar_fit`, fit Δm danych N_D*=5588±586): 0.105 we wszystkich kanałach — odjęty kwadratowo w wierszach „bez N_D*”, bo wynik końcowy podaje N_D* osobno (powtórka NDS 12.11: N_D* = 4436, 21,9 %).

## Przeskalowanie do μ_wsp. = 1,495 (≈ 1,50)

Baza N_D* = 4436 (21,9 %) z `NDS_ndstar/rerun_2026-11-12/rerun.md` (decyzja z 12.11.2026).

| wielkość | względna | × 1,495 |
|---|---|---|
| (a) syst. tune (MPI ⊕ CR ⊕ Monash) | 0.117 | ± 0,17 |
| (b) syst. łączna bez N_D*, z tune | 0.284 | ± 0,42 |
| (c) syst. łączna bez N_D*, bez tune | 0.258 | ± 0,39 |
| info: syst. łączna z N_D* (0,105), z tune | 0.302 | ± 0,45 |
| info: N_D* (NDS, rerun 12.11) | 0.219 | ± 0,33 |
| info: tune CR asym. (↑ / ↓) | +0.121 / -0.040 | +0,18 / -0,06 |

Per kanał tylko wartości względne (tabela wyżej).

Poprzednia baza (historycznie, zapis z 12.11 przed zmianą): μ = 1,49, N_D* = 4443 (20,6 %) → tune ± 0,17, syst bez N_D* ± 0,42 (asym. +0,43/−0,39), N_D* ± 0,31.

## Zapis wyniku (baza N_D* od 12.11)

> **μ_wsp. = 1,50 +0,48/−0,46 (stat; profil +0,41/−0,39 × 1,18 z jackknife po zdarzeniach, 13.11) ± 0,42 (syst, bez N_D*) ± 0,33 (N_D*)**

Składnik tune w ± 0,42 jest symetryzowany (|δ|); w rzeczywistości CRqcd przesuwa μ tylko w górę (CR: +0,18 / -0,06). Z = 4,06σ (4,11σ) bez zmian — normalizacja nie zmienia Z.
Wariant asymetryczny (informacyjnie, tune jako skrajne 1/(1+δ)−1 per grupa: +0.126 / −0.059): syst. bez N_D* = +0,43 / −0,40.

