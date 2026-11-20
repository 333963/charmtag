# charmtag — klasyfikator pojedynczego D + liczba D z fitu masy (etap L, 2026-11-02)

Lokalne narzędzie (nic na GitLabie). Nominalny model: **XGBoost (hist), wariant H** = 7 cech „frozen” + cosθ_P, |cosθ*|, |IP sig|, min|d0| córek, max|z0| córek, max|η| córek, dE/dx K i π z pędami. **Izolacja (ptcone) zakazana** — `features.forbidden`, sprawdzane w kodzie i testach.

- Sygnał: MC forcedDecay ND (`fd_nd`), truth ∈ `signal.truth_values` (0/1).
- Tło: sidebandy danych 428770 `none` (D0 |m−m0|>75 MeV, D+ >70 MeV; dla D+ weto D*: m(Kππ)−m(Kπ)<153 MeV), ważone do (pT,|η|) sygnału.
- k-fold = 5 po splitmix64(run, event_number); każdy kandydat (także w `apply`) dostaje wynik modelu foldu, na którym NIE był trenowany → wyniki out-of-fold. Biny pT 3–5, 5–8, 8–30; progi WP 80/70/50% per bin.

## Komendy
```
charmtag inspect  -c config/mac_v11.yaml                 # mapowanie kolumn → inspect.json
charmtag train    -c ... --species D0|Dplus [--variant H|frozen7]
charmtag apply    -c ... --species D0                     # scores/<sp>_<var>/<próbka>.h5 (row = indeks w `events`)
charmtag fit      -c ... --species D0                     # iminuit: Gauss/2-Gauss + exp/pol2 + refleksja K↔π (D0)
charmtag validate -c ... --species D0                     # ε_B@WP, R_win, kontrola D0 z D*, H vs frozen7
charmtag report   -c ... --species D0                     # report/<sp>.md, *_summary.json, PNG
charmtag all      -c ... --species D0                     # wszystko, oba warianty
charmtag synth    --h5-dir /tmp/x                         # syntetyczne h5 (tylko testy)
python -m pytest -q tests
```
Pełny bieg na Macu: dwuklik `~/Documents/charm/START_charmtag.command` → log `~/Documents/charm/START_charmtag.log`, wyniki `runs/v11_full/`.

## Konfiguracja (`config/default.yaml`, nadpisanie `-c plik.yaml`)
- `columns`: nazwa logiczna → lista kandydatów nazw kolumn (pierwsza istniejąca; `{i}` = córka; `re:` = regex).
- `species`: m0, zakres masy, okna, sidebandy, pasmo wewnętrzne do R_win, indeks kaonu/pionów, weto D*.
- `fit`: `signal_model` gauss|gauss2 (kształt z MC truth po cięciu), `background_model` exp|pol2, `sigma_mode` free|fixed|constrained (domyślnie więz 10% do σ z MC), `reflection_mode` template|gauss|off, `reflection_ratio` (null = z MC).

## Walidacja (validate)
- ε_B przy WP 80/70/50 (OOF, ważone i surowe), AUC per bin pT.
- R_win: (a) dane — [N(pasmo 50–75 MeV)/N(sidebandy)] po/przed cięciem per bin pT (pasmo nieużywane w treningu; ogon sygnału i refleksje mogą dawać R>1), (b) dane — lewy/prawy sideband (kształt tła), (c) MC mb_nd tło (truth) — okno sygnału vs sidebandy.
- D0 z D*: kandydat D0 dopasowany do D* w tym samym (run, event) (|m_D0 − (m_D*−Δm)|<1 MeV, pT, η), Δm 143.5–147.5 MeV minus sideband 150–165 MeV (skala szerokością) → ε(WP) dane vs MC.
- ε klasyfikatora z fitu: N(WP)/N(bez cięcia) w danych vs nominalne ε z MC.
- Porównanie H vs frozen7: zysk odrzucania 1/ε_B.

## Znane ograniczenia
- Lista 7 cech „frozen” i nazwy kolumn h5 v1.1 — ustawione na podstawie dokumentów; zweryfikować z `single_d_classifier_v20260925_01/code/common.py` i `charmtag inspect`.
- fd_nd scalony starą wersją merge: event_number nieunikalny (1…10000/job) → foldy MC nieidealnie rozłączne między jobami (sygnał MC i tak nie jest w danych); dopasowanie D* w MC wymaga zgodności m, pT, η.
- Skala sidebandu Δm w kontroli D* = stosunek szerokości (przybliżenie płaskiego tła).
- Dane = ½ lumi (joby 25–47).

## Pary D–D (etap O, 2026-11-02) — `charmtag pairs` / `charmtag fit2d`
Moduł `charmtag/pairs/` (addytywny; w `cli.py` tylko przekierowanie komend `pairs`, `fit2d`). Port M7 v3 + plan B z
`unblinding_freeze_v20260927_01` / `odsloniecie_wynik_v20260927_01` (+ systematyki kształtu z `wynik_stabilnosc`, więz κ partnera z T12).
```
charmtag pairs -c config/mac_v11.yaml            # nogi, kategorie, szablony MC, normalizacja, histogramy BEZ centrum 2D
charmtag fit2d -c config/mac_v11.yaml [--ntoys 30]   # ślepo: jakość fitu bez centrum, κ, Z Asimov, toys, przesunięcia systematyk
charmtag pairs|fit2d ... --unblind               # odsłonięcie (oba kroki wymagają flagi)
charmtag fit2d --closure [~/Documents/charm/unblinding_freeze_v20260927_01]   # domknięcie 09-27 (PASS: μ 1,37/1,62/1,57; Z 3,6σ)
```
Pełny bieg: dwuklik `~/Documents/charm/START_charmtag_O.command` → `runs/v11_full/pairs/` (fit2d_blind.md/json, pair_inputs.*).
Opis: notatka `charm-O-pary.md`.
