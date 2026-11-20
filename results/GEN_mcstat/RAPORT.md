# GEN_mcstat: czy prawda pythia_dd zastąpi n_pair z rekonstrukcji (16.11.2026)

Obliczenia lokalne, tylko odczyt wejść. Bez danych, bez μ. lxplus niedostępny (brak socketu ControlMastera), więc EOS sprawdzony tylko przez lokalne kopie.
Kod: `gen_inventory.py` → `out/gen_inventory_A3_trk{0.80,0.85,0.90}.json`.

## Werdykt: **częściowo**
- **Statystyka prawdy wystarcza z zapasem** (×100 par, n_eff po wagach ×60–100). Nowa produkcja 5·10⁸ nie jest potrzebna.
- **Zmienne wag są w prawdzie:** `E <ev> <proc> <nch> <leadpt>`, gdzie n_ch to naładowane stabilne cząstki z pT > 0,5 GeV i |η| < 2,5 (generatorowy odpowiednik `nTracksPt500Eta25`), a leadpt to pT wiodącej z nich.
  **Brakuje samej mapy w2** (`w2d.npy`, 15×13 komórek n_trk × lead_pt). Nie ma jej lokalnie, nie ma też per-zdarzeniowych (n_trk, lead_pt) danych ani mb_nd.
- **Brakuje map ε dla selekcji z `direct`.** Mapy fd_nd nie domykają się na mb_nd dla D⁺ względem D*, rozbieżność wynosi 32 % (§4).
- Do zrobienia nie jest więc nowy generator, tylko dwa małe zadania odczytu na lxplus oraz opcjonalny mały bieg walidacyjny (§5).

## 1. Inwentarz
| próbka (EOS `charm_q/truth/gen/`) | zdarzenia | lokalnie | tune |
|---|---|---|---|
| A3 | 5,0·10⁸ (500 plików × 10⁶) | wszystkie 500 (`AFP_ddpred/gen/A3`) | A3 |
| syst_MPIup / MPIdown / Monash | po 1,0·10⁸ | po 100 plików | A3 ± pT0Ref / Monash |
| syst_CRup / CRdown / CRqcd | 1,0 / 1,0 / 0,88·10⁸ | tylko `truth_norm.json` (`runs/lxplus/truth/combined_*`) | A3 + CR |
| `truth/effmap/fd_{nd,sd,dd,cd}` | ntuple 1,97 / 0,99 / 1,00 / 0,98 M | tak (`AFP_ddpred/effmap`) | pełna symulacja: prawdziwe D + `ntrk` reco + `matched` |

Skład A3 według procesu (`info.code()`): ND 353,7 M, SD 40,9 + 40,9 M, DD 56,3 M, CD 8,2 M.
W pliku są tylko zdarzenia z D o pT > 2 GeV i |η| < 3. Dla wszystkich zdarzeń zapisany jest tylko 1D `nch_hist` (wszystkie procesy razem), bez lead_pt i bez podziału na proces.
**Wag zdarzeń w plikach nie ma**, Pythia generuje je bez ważenia.

Pary w akceptacji (pT > 3, |η| < 2,2; obie nogi `mode==1`, czyli naturalne BR do Kπ / Kππ; D⁰ także z D*). A3, tylko ND, bo mb_nd zawiera tylko ND:

| | D⁰D⁰ | D⁺D⁰ | D⁺D⁺ | D* |
|---|---|---|---|---|
| prawdziwe pary (w tym LS) | 2310 (416) | 6091 (1166) | 4306 (834) | 125 705 |
| n_eff po ε×ρ | 1916 | 4685 | 3098 | — |
| n_eff po ε×ρ×w1 (proxy n_trk = 0,85·n_ch) | 555–495 | 1400–1340 | 920–830 | ~10⁴ |
| **reco mb_nd (direct)**: n / n_eff(w2) | 22 / 8,2 | 41 / 13,5 | 15 / 7,5 | 1344 |
| stat. względna: direct (n_eff) → prawda (n_eff w1) | 0,35 → **0,045** | 0,27 → **0,027** | 0,37 → **0,034** | |

Mediana n_ch par wynosi 44–46, a mediana n_trk par w reco 36–47: to te same zdarzenia o dużej krotności. Udział lead > 4 GeV to 0,69 / 0,63 / 0,53, podczas gdy w danych 0,645. Wagi dla par są więc umiarkowane.
Wszystkie procesy razem dają +2 % par. Warianty tune mają po 1/5 statystyki A3 (~460 / 1220 / 860 par), co wystarcza na warianty z dokładnością ~5 %.

## 2. Wagi 2D na poziomie generatora
- Definicja (`etap1_definicje_masy_v20260925_01/code/w2d.py`, `data_events.py`): w2(n_trk, lead_pt) = p_dane/p_MC(mb_nd, reco).
  Zmienne: n_trk = `EventInfo_nTracksPt500Eta25`, lead_pt = max pT śladu przy |η| < 2,5.
  Biny: NE = [−1, 2, 5, …, 85, 10⁴], PE = [0, 0,5, …, 15, 10⁴]. Komórki z < 20 zdarzeniami MC dostają wagę 0.
  `c_w = (X/D*)_w / (X/D*)_unw` (`build.py: weight_factors`).
- Na generatorze proxy jest bezpośrednie: n_trk ≈ ε_trk·n_ch (tak jak `truth_combine.py`, ε_trk = 0,85 ± 0,05), lead_pt ≈ leadpt.
  Do stosunków na D* nie potrzeba normalizacji Σw po wszystkich zdarzeniach, bo się skraca.
- **Test proxy na w1** (mapa 1D jest lokalnie, `weights_ntrk.csv`):
  | c_w1 | D⁰ | D⁺ | D⁰D⁰ | D⁺D⁰ | D⁺D⁺ |
  |---|---|---|---|---|---|
  | prawda × ε, proxy 0,85 (0,80–0,90) | 1,010 | 1,047 (1,043–1,049) | 1,58 (1,53–1,61) | 1,50 (1,49–1,52) | 1,58 (1,52–1,61) |
  | reco mb_nd (`t4_rates` w1) | 0,960 | 0,957 | 1,23 (n_eff 4) | 1,12 (n_eff 19) | 2,04 (n_eff 5) |

  Pojedyncze D zgadzają się na poziomie 5–9 %, czyli w granicach statystyki ważonych D* w reco (n_eff kilkaset). Pary zgadzają się w granicach szumu reco.
  Zależność od ε_trk wynosi ±3 % dla par.
- Proxy dla w2: ta sama procedura z mapą w2. Niedomknięcie typu n_trk ≈ ε·n_ch (zamiast macierzy odpowiedzi binomial + fałszywe ślady, przy w do 43 w ogonie) trzeba sprawdzić na mb_nd (§5 b).

## 3. Tryb `truth` w kodzie: już istnieje, ale nie rozwiązuje problemu
- `charmtag/pairs/truthnorm.py` + `build.py:556`, `charmtag pairs --norm truth`, test `tests/test_truthnorm.py`.
  Wzór: n_pair = N_D* · [pary/D*]_reco,Pythia(`fd_nd|trk0.85`) · c_w. Wejście `truth_norm.json` z `Q_lxplus/truth/truth_combine.py`: 5·10⁸ A3, |η| < 2,5, mapa ε(pT, |η|, n_trk) z fd_nd, ρ = 1,060 / 1,032 / 1,020.
- Różnice względem direct: (i) pary/D* z prawdy × ε zamiast z liczenia 22/41/15 par w mb_nd; (ii) **c_w pozostaje z `mb_rates_2d.json`**, czyli z tych samych 21/41/15 par.
  W efekcie zostaje szum wag (n_eff 8 / 14 / 8), a `rel_stat` (< 0,2 %) go pomija. To tylko pół rozwiązania.
- Wartości: pary/D* = 0,01257 / 0,0225 / 0,01049. Direct na mb_nd (h5) daje 22 / 41 / 15 na 1344 D*, czyli 0,0164 / 0,0305 / 0,0112.
  Przy N_D* = 5584 i tym samym c_w truth daje ≈ 96 / 151 / 58 par, a direct 126 / 205 / 61, czyli **truth/direct = 0,77 / 0,74 / 0,94**.
  Tryb truth nie był użyty na danych (brak `truth_norm_json` w configach PF2). Na prawdzie oparto tylko f_LS (`PD_fls_truth`).

## 4. Domknięcie na mb_nd (istniejące dane, bez nowego liczenia na reco)
- `truth_norm.json` → `closure.fd_nd_on_mb_nd`: mapa fd_nd na prawdzie mb_nd (19,7 M zdarzeń), zmierzone/przewidziane.
  Pojedyncze: D⁰ 0,725 ± 0,010, **D⁺ 0,949 ± 0,011**, D* 0,721 ± 0,025. Pary: 0,66 ± 0,14 / 0,80 ± 0,13 / 0,90 ± 0,23.
- W stosunku do D*, którego używa normalizacja:
  - D⁰/D* = 1,01 ± 0,03 (OK);
  - **D⁺/D* = 1,32 ± 0,04** (mapa fd zaniża D⁺ w mb_nd względem D*);
  - pary/D* = 0,91 ± 0,20 / 1,11 ± 0,18 / 1,25 ± 0,32 (ważona średnia 1,05 ± 0,12).
- Direct (reco, D* w oknie fitu) vs moja prawda × mapa przy |η| < 2,2 (0,01142 / 0,0232 / 0,01216 na D*): obserwowane 22 / 41 / 15 wobec przewidywanych 15,3 / 31,2 / 16,3, **razem 78 wobec 63 (1,24 ± 0,14)**.
  Po uwzględnieniu D⁺/D* = 1,32 przewidywane jest 85, czyli 0,92 ± 0,10.
  Dodatkowo D* w direct (1344, okno Δm/D⁰) to nie to samo co D* „matched” w mapie (1613). Tej selekcji brak w mapie.
- **Wniosek:** zgodność jest w granicach ~15–25 %, ale moc testu ogranicza ta sama statystyka 78 par.
  Mapy ε muszą odpowiadać dokładnie selekcji direct, czyli poprawne D z `legs_for_truth` i D* z oknem `dstar_norm`. Najlepiej budować je z mb_nd (9,5 k D⁰, 8,7 k D⁺, 1,3 k D*), a z fd_nd brać tylko ρ.
  Bez tego tryb truth wnosi nowy wiersz „transfer ε fd→mb” rzędu 30 % dla D⁺.

## 5. Czego brakuje (zamiast nowego generatora)
a) **Odtworzenie mapy w2 na lxplus** (odczyt ntupli, bez condora albo 1 mały job):
   - `data_events.py` na plikach danych 428770: zmienne (n_trk, lead_pt) i histogram 2D;
   - (n_trk, lead_pt) na 80 plikach mb_nd z gałęzi `Tracks`.
   Koszt ~2–5 CPU·h, wynik `w2d.npy` + krawędzie (kB). Sprawdzenie: odtworzone c_w(D⁰, D⁺, pary) = `mb_rates_2d.json` (0,935 / 0,793 / 1,374 / 1,201 / 0,981).

b) **Odpowiedź n_ch → n_trk:** jeśli ntuple mb_nd mają `TruthParticles` (do sprawdzenia), to 2D macierz (n_ch, lead_gen) → (n_trk, lead_reco) z 2–3 plików mb_nd.
   Jeśli nie, zostaje proxy ε_trk·n_ch z wariantami 0,80 / 0,90 (±3 % c_w).

c) **Mapy ε dla selekcji direct z mb_nd** (`effmap.py` jest gotowe; plik testowy `Q_lxplus/_test/effmap_v9_mbnd_000027.npz` istnieje), plus flaga „w oknie D*”. Koszt ~10 CPU·h (80 plików).

d) Opcjonalnie: **bieg walidacyjny** `pythia_dd` z dopisanym 2D histogramem (n_ch, lead) wszystkich zdarzeń per proces, 2·10⁷ ND, ~10 CPU·h.
   Pozwala porównać rozkład zdarzeń A3 8.317 z mb_nd mc23 (dziś możliwe tylko 1D: mediana n_ch·0,85 wszystkich procesów ≈ 6 wobec 10 w reco ND; to spójne, bo jest tam ~29 % dyfrakcji).
   **Nowej produkcji 5·10⁸ nie trzeba.**

## 6. Oczekiwane nowe wiersze (konwencja wspólna Σw·δ, w = 0,298 / 0,580 / 0,122, × 1,495)
- **Stat. MC par:** prawda (n_eff w1 jako proxy w2; w reco n_eff(w2)/n ≥ n_eff(w1)/n) ⊕ ρ fd (0,9 / 0,8 / 1,8 %) ⊕ ε_D* z mb_nd (~2,7 %).
  Daje to 5,2 / 3,8 / 4,6 %, wspólne 0,043, czyli **±0,06 na μ (było ±0,46)**.
- **Model wag (w2↔w1):** szum (dziś dominujący, n_eff 4–19) znika, bo w1↔w2 liczone jest na ~500–1400 n_eff par (szum ~5 %, ±0,02 na μ).
  Zostaje część fizyczna, której bez mapy w2 nie umiem policzyć. Górną granicą jest proxy pojedynczych D z VER (0,05 / 0,19 / 0,31, czyli ±0,24).
  Szacunek: **±0,1–0,24**, prawdopodobnie bliżej dolnej wartości dla D⁰D⁰ i D⁺D⁰.
- **Nowe wiersze trybu truth:**
  - transfer ε (po mapach z mb_nd: stat ε ~3 %, ±0,05);
  - proxy n_ch → n_trk (±3 %, ±0,05);
  - tune generatora (już jest, ±0,17; warianty mają tę samą strukturę).
- Razem stat MC + model + nowe: ok. ±0,15–0,27 zamiast ±0,52 (0,46 ⊕ 0,24).
