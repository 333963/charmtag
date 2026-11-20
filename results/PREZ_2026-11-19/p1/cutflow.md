# Cutflow (S = truth-matched c/b, B = reszta poza refleksjami K↔π)

## Surowe ntuple (finder CharmAnalysis = preselekcja), lokalne pliki h5 v1.2 (mb_nd 24, fd_nd 8)

| próbka | gatunek | krok | S | B |
|---|---|---|---|---|
| mb_nd | D0 | candidates in ntuple (finder preselection) | 4,317 | 952,367 |
| mb_nd | D0 | mass range | 4,317 | 952,360 |
| mb_nd | D0 | 3 < pT < 30 GeV | 4,313 | 952,273 |
| mb_nd | Dplus | candidates in ntuple (finder preselection) | 2,631 | 3,330,776 |
| mb_nd | Dplus | mass range | 2,631 | 3,330,766 |
| mb_nd | Dplus | 3 < pT < 30 GeV | 2,629 | 3,330,531 |
| mb_nd | Dstar | candidates in ntuple (finder preselection) | 510 | 324,817 |
| mb_nd | Dstar | (no m cut; Δm range in fit) | 510 | 324,817 |
| mb_nd | Dstar | 3 < pT < 30 GeV | 510 | 324,800 |
| fd_nd | D0 | candidates in ntuple (finder preselection) | 124,215 | 501,503 |
| fd_nd | D0 | mass range | 124,214 | 501,502 |
| fd_nd | D0 | 3 < pT < 30 GeV | 124,049 | 501,434 |
| fd_nd | Dplus | candidates in ntuple (finder preselection) | 58,001 | 1,981,398 |
| fd_nd | Dplus | mass range | 58,001 | 1,981,394 |
| fd_nd | Dplus | 3 < pT < 30 GeV | 57,929 | 1,981,135 |
| fd_nd | Dstar | candidates in ntuple (finder preselection) | 25,427 | 206,448 |
| fd_nd | Dstar | (no m cut; Δm range in fit) | 25,427 | 206,448 |
| fd_nd | Dstar | 3 < pT < 30 GeV | 25,378 | 206,414 |

## Pełna statystyka (scores v11_full; wejście klasyfikatora = zakres masy i pT już nałożone)

### D0, MC minimum bias (mb_nd, Pythia A3)

| krok | S | B | S/B |
|---|---|---|---|
| mass range + 3<pT<30 GeV (classifier input) | 9,534 | 3,170,411 | 0.00301 |
| signal window |m − m0| < 55 MeV | 8,825 | 1,404,246 | 0.00628 |
| score ≥ WP80 | 7,098 | 122,334 | 0.058 |
| score ≥ WP70 | 6,187 | 55,933 | 0.111 |
| score ≥ WP50 | 4,432 | 11,360 | 0.39 |

### D0, MC forced decay (fd_nd)

| krok | S | B | S/B |
|---|---|---|---|
| mass range + 3<pT<30 GeV (classifier input) | 216,175 | 1,297,662 | 0.167 |
| signal window |m − m0| < 55 MeV | 200,101 | 574,737 | 0.348 |
| score ≥ WP80 | 161,210 | 66,015 | 2.44 |
| score ≥ WP70 | 141,293 | 35,849 | 3.94 |
| score ≥ WP50 | 101,321 | 11,570 | 8.76 |

### Dplus, MC minimum bias (mb_nd, Pythia A3)

| krok | S | B | S/B |
|---|---|---|---|
| mass range + 3<pT<30 GeV (classifier input) | 8,727 | 11,078,038 | 0.000788 |
| D* veto (min(m − m_Kπ) ≥ 160 MeV) | 8,665 | 10,846,291 | 0.000799 |
| signal window |m − m0| < 45 MeV | 7,912 | 1,910,643 | 0.00414 |
| score ≥ WP80 | 6,391 | 47,687 | 0.134 |
| score ≥ WP70 | 5,600 | 21,196 | 0.264 |
| score ≥ WP50 | 3,972 | 5,358 | 0.741 |

### Dplus, MC forced decay (fd_nd)

| krok | S | B | S/B |
|---|---|---|---|
| mass range + 3<pT<30 GeV (classifier input) | 149,705 | 5,135,226 | 0.0292 |
| D* veto (min(m − m_Kπ) ≥ 160 MeV) | 148,440 | 4,975,861 | 0.0298 |
| signal window |m − m0| < 45 MeV | 135,929 | 858,947 | 0.158 |
| score ≥ WP80 | 109,277 | 52,874 | 2.07 |
| score ≥ WP70 | 95,749 | 30,013 | 3.19 |
| score ≥ WP50 | 68,580 | 10,750 | 6.38 |

