# gen_norm — wyniki (MC + rozkład zdarzeń danych; BEZ par danych, bez μ)

MC mb_nd: 78 plików, 19220986 zdarzeń; dane 428770: 97 plików (ref ['anchors:303497']: 441671 zdarzeń, wszystkie 7989131).

## (a) w2 — test akceptacji: **PASS** (tol |Δc_w| < 0.005)

| | odtworzone | ref | Δ | wszystkie dane |
|---|---|---|---|---|
| c_w D0 | 0.9179 | 0.9348 | -0.0168 | 0.9182 |
| c_w Dplus | 0.8572 | 0.7930 | +0.0643 | 0.8585 |
| c_w D0D0 | 1.4009 | 1.3743 | +0.0266 | 1.4016 |
| c_w DpD0 | 1.6772 | 1.2009 | +0.4763 | 1.6700 |
| c_w DpDp | 1.0585 | 0.9810 | +0.0776 | 1.0741 |

Liczniki: events 19220986 (≠ 19714175), dstar_n 1746 (≠ 1314), D0_n 12777 (≠ 9448), Dplus_n 8466 (≠ 8671), pairs_D0D0_n 40 (≠ 21), pairs_DpD0_n 55 (≠ 41), pairs_DpDp_n 15 (OK)

B1 v2.0 (blokujące): B1b_anchors_found OK (303497), B1b_ref_events OK (441690), B1a78_Dplus_n OK (8466), B1a78_pairs_DpDp_n OK (15), B1b78_Dplus_correct_per_event_unw OK (0.00044045607233676775), B1b78_Dplus_correct_per_event_w OK (0.0023675488773733377), B1b78_pairs_DpDp_per_event_unw OK (7.803970098100066e-07), B1b78_pairs_DpDp_per_event_w OK (5.179821073397761e-06), B1a_legs_without_event OK (0), B1a_dstar_without_event OK (0), B1a_event_key_duplicates OK (0), B1a_retag_bijective OK ({'000003': '000009', '000004': '000010', '000005': '000003', '000006': '000004', '000007': '000005', '000008': '000006', '000009': '000007', '000010': '000008'})

D0/D* i c_w (v2.0 ≠ v1.1 w kandydatach D0/D*) — tylko informacyjnie (B1c).

## (b) nogi direct ↔ prawda

Nogi (gn_legs): {'n_dstar': 1344, 'n_single': {'10': 9534, '4': 8665}, 'n_pairs': {'D0D0': 22, 'DpD0': 41, 'DpDp': 15}}; zgodność z nominałem: True.
Dopasowanie: D0 9304 nóg, bez prawdy 58, poza mode 0, poza akc. mapy 0; Dplus 8460 nóg, bez prawdy 9, poza mode 0, poza akc. mapy 0; Dstar 1312 nóg, bez prawdy 0, poza mode 0, poza akc. mapy 0

c_w selekcji h5 (w2 ref): D0 0.904, Dplus 0.769, D0D0 1.333, DpD0 1.199, DpDp 0.961

| domknięcie (split-half) | direct n (n_eff) | truth/direct (mult) | truth/direct (any) | pull |
|---|---|---|---|---|
| D0D0 | 22 (8.7) | 1.039 | 1.039 | 0.114 |
| DpD0 | 40 (13.3) | 1.315 | 1.314 | 1.148 |
| DpDp | 15 (7.5) | 1.193 | 1.193 | 0.530 |

## (c) n_ch → n_trk: POMINIĘTE — ntuple mb_nd v1.1 nie mają TruthParticles (tylko TruthDmesons) — n_ch→n_trk niedostępne; zostaje proxy n_trk = ε_trk·n_ch (0,80/0,85/0,90)
