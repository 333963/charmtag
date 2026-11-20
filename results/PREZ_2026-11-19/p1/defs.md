# (1) Definicje sygnału i tła (MC, truth-matching)

Źródło etykiet: `Dmeson_truth_category` w ntuplach CharmAnalysis (h5 v1.1/v1.2). Mapowanie liczba ↔ `Dmeson_truth_channel`
sprawdzone na `mb_nd_000001` (19.11); grupowanie jak `charmpy.core.cuts_utils.categorize_dmeson` + refleksje z charmtag.

| grupa (legenda wykresów) | kody | truth_channel | znaczenie |
|---|---|---|---|
| **Signal** | 0, 1 | MatchedC, MatchedB | wszystkie córki dopasowane do prawdziwego D tego gatunku i kanału; 0 = prompt (c), 1 = z rozpadu b |
| Reflection K↔π (tylko D⁰) | 100, 101 | (MatchedC/B + zamiana hipotez) | te same ślady co sygnał, zamienione przypisanie K↔π; v1.1: hipoteza dalsza od m₀ (`closest_mass`), truth + 100 |
| Other charm / other decay mode | 2, 3, 6, 7 | MatchedCharmC/B, DecayModeMismatchC/B | ślady z prawdziwego hadronu powabnego, ale innego gatunku albo innego kanału rozpadu (np. D⁺→KKπ, D_s) |
| Partially reco | 8–11 | MisMatchNoB/SomeB/AllB (+pdg) | część córek z hadronu c/b, reszta przypadkowa |
| Combinatorial | 4, 5, 12–15 | MatchedBottom, MatchedLight, Other* | brak związku z hadronem powabnym (w tym ślady z b/lekkich) |

„Tło” w treningu i na wykresach score = wszystko poza sygnałem i refleksjami (0, 1, 100, 101).
D*: sygnał = prawdziwe D* i prawdziwe D⁰ (truth_ds ∈ {0,1} i truth ∈ {0,1}).

Skład okna sygnału (|m − m₀| < 55 MeV D⁰, 45 MeV D⁺), MC mb_nd, bez cięcia MVA → po WP70:
- D⁰: sygnał 0,6 % → 9,7 %; refleksje 0,15 % → 2,9 %; inny charm 0,35 % → 2,3 %; partially reco 5,3 % → 10,3 %; kombinatoryka 93,6 % → 74,8 %.
- D⁺: sygnał 0,4 % → 20,9 %; inny charm 0,14 % → 2,3 %; partially reco 8,7 % → 34,8 %; kombinatoryka 90,7 % → 42,0 %.
(fd_nd, po WP70: D⁰ sygnał 64 %, refleksje 20 %; D⁺ sygnał 76 %, partially reco 21 %.)
