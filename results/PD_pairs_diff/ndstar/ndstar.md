# N_D* — spójna normalizacja par D–D i μ na wspólnej bazie (etap ndstar, 2026-11-09)

> **Nota 2026-11-13:** od 13.11 stat +0,48/−0,46 (× 1,18 z jackknife po zdarzeniach; Z = 4,06σ bez zmian, ≈ 3,4σ z korelacjami zdarzeń), zob. `PF2_lxplus/runs/lxplus/unblind/UNB_WYNIK.md` §1b.

Obliczenia lokalne, bez lxplus i bez nowych biegów par. Pliki w tym katalogu:
- `ndstar_rescale.py` → `ndstar_rescale_<N>.json`: przeliczenie μ z JSON-ów (bez fitów).
- `ndstar_decomp.py` → `ndstar_decomp.json`: lokalny refit samego Δm, histogram V `dstar_nom.npz`, kilka sekund.

Liczby ATLAS — tylko do użytku wewnętrznego.

## TL;DR
- **Rekomendacja: N_D* = 4436 +715/−645 (stat, MINOS; HESSE 671) ± 692 (syst) ≈ ±22 %.** Syst = NDS 605 ⊕ „źródło kształtu MC” ±336 (§1.3).
  Aktualizacja 12.11.2026 po poprawce NLL w `fit_single` (10.11): było 4443 ± 612 ± 681 (±20,6 %); zob. `NDS_ndstar/rerun_2026-11-12/rerun.md`.
  Tabele §3 i JSON-y `ndstar_rescale_*.json` zostają na bazie 4443 (różnica ×1,0014).
- μ ∝ 1/N_D* **dokładnie** (§2). Przeliczenie jest analityczne, nie trzeba nic liczyć ponownie.
- Bieg PD dostał 4022, bo kształt sygnału Δm (2 Gaussy) wziął z MC `signal_mc` potoku par: k = 2,11 zamiast 2,73 z V/NDS, węższe ogony.
  Kształt odpowiada za ≈ −8 %, binning 0,1 MeV za +0,5 %, reszta (−1,8 %) to histogram wejściowy i model tła (§1.3).
  Te same dane, ten sam typ fitu (`fit_single`).
- **Na wspólnej bazie 4443 nominał i PD się zgadzają.** Σμ_i f_i / μ_nom = 0,99 / 1,01 / 0,95 (było 1,37 / 1,41 / 1,31).
- **Nowa uwaga:** μ wspólne z JSON-ów (1,25 nominał, 1,50 PD) jest skwantowane do siatki profilu 0,25. `fit2d.combine` interpoluje liniowo,
  więc minimum ląduje zawsze na węźle. Minimum z gładkiej krzywej (spline / parabola przez 3 i 5 punktów) daje
  nominał 1,19–1,20 i PD 1,61 (w bazach użytych w biegach).
  **Na bazie 4436: nominał μ = 1,50 [1,10; 1,90], PD μ = 1,46 [1,08; 1,84].** Na bazie 4443 było 1,49 i 1,45. „1,57” jest zawyżone o ~5 % przez siatkę.
- f_LS, Z i testy kształtu (Δχ² Pythia/DPS) są niezmiennicze. Potwierdzone numerycznie: f_LS po przeliczeniu jest identyczne do 4 cyfr.

## 1. Skąd każde N_D*

| N_D* | metoda / plik | dane Δm | kształt sygnału | status |
|---|---|---|---|---|
| **5584 ± 585** (5588 ± 586 w NDS) | `massfit.fit_dm` — **fallback** w `pairs/build.py:n_dstar_data` (pf_pairs, unblind) | dane, okno D0 ±40 MeV, [139,6; 165], 656 845 wejść | Gauss σ≤2,5 MeV (`fit_dm`) | `pair_inputs.json` `meta.dstar_fit_data`: `method: massfit.fit_dm (fallback)`, `fallback_reason: refit niepoprawny (valid=False, N=4003.9)`. Znany artefakt V (sztywne tło + limit σ) — zawyżone o ~25 %, głównie bin pT 3–5 (stary 2782 vs nowy 580–790, B2b/NDS). **Odrzucone.** |
| **4436 ± 671 ± 605** (MINOS −645/+715; przed poprawką NLL 4443 ± 612 ± 591) | NDS: `dstar_refit.fit_single` (threshold_pol2), `NDS_ndstar/rerun_2026-11-12/ndstar_variants.py` → `rerun_2026-11-12/new/ndstar_full.json` | V `dstar_nom.npz` `norm_all` = ta sama selekcja (656 845 wejść w [139,6;165]), bin 0,05 MeV, [m_π+0,9; 160] → 482 979 wejść | 2G z MC V (`dstar_post.json` `norm.shape_mc`): σ=0,743, k=2,73, f=0,718 | valid, χ² 336/383; 14 wariantów, wszystkie valid; syst = bkg 534, range 211, bin 107, sig 159 (przed poprawką: 521, 202, 100, 163). Wartość 4443 użyta do korekty ×1,257 → „1,57”. |
| **4220 ± 424** | B2b: suma per bin pT pf95 (pass/fail przy WP95), okno D0 **±55**, próbka EV par D0←D* po dedup (`B2b_dstar/b2b_ndstar.json` `per_bin_nom55`) | inna próbka (EV) i **inne okno** | `fit_passfail` z kształtem MC | Inna obserwabla: `mc_rates` liczy D* w MC w oknie ±40, a ±55 nie pasuje do tego mianownika. Kontrola zgodna: 8–30 GeV 2110 vs V 2247; ±40 per bin daje sumę 3697. **Nie jako normalizacja**, tylko kontrola (−5 % od 4443, w granicach syst). |
| **4022 ± 562** | bieg PD `diff_v12`: `pairs/build.py:n_dstar_data` → `dstar_refit.fit_single` (nowy kod NDS z 5.11, bez fallbacku) | dane z potoku par (`D_d`, te same 656 845 wejść co nominał), bin 0,1 MeV, 482 264 wejść | 2G z MC potoku par (`legs_for("signal_mc")`, truth∈{0,1}, `shape_from_mc` w oknie 143–148,5): σ=0,707, **k=2,11**, f=0,699 | valid, χ² 188/188 (`runs/lxplus_v12/pairs/pair_inputs.json`). Ten sam MC kształt w nominalnym biegu pf_pairs dał 4004 (invalid → fallback 5584). |

### 1.1 Jak N_D* wchodzi w oczekiwane S (`pairs/build.py` l. 509–527, `model.py`)
- `norm_mode: direct`: n_pair[ch] = N_D* × n_pairs_MC[ch] / n_dstar_MC × c_w[ch], a n_real[dt] = N_D* × n_single_MC / n_dstar_MC × c_w.
  Alternatywy `frozen` i `truth` też są ∝ N_D*. Obie wielkości są **liniowe w N_D*, we wszystkich kanałach, tym samym czynnikiem**.
  Sprawdzone: n_pair i n_real nominał/PD = 1,3884 = 5583,6/4021,6.
- PD per bin: n_pair_i = n_pair × f_i, n_real_i = n_real × g_i (`diff.py:diff_inputs_json`). f_i i g_i nie zależą od N_D*.
- `model.expected`: SS = n_pair × Y_a Y_b (sygnał), SB/BS = n_real × Y × NF (prawdziwa × fałszywa noga), BB z danych SB×SB (nie zależy od N_D*).

### 1.2 Dlaczego bieg PD dostał 4022, a nie 4443
W kodzie NDS z 5.11 `fit_single` jest nominałem. Nominalny bieg pf_pairs (2.0.2 z poprawką ID śladu z 5.11) też wołał refit, ale ten był invalid → fallback `fit_dm`.
Bieg PD miał nowy, odporny `fit_single` (valid). Kształt sygnału brał jednak, jak każdy bieg par, z MC potoku par, a nie z kształtu V, którego użył NDS.

### 1.3 Rozkład 4443 → 4022 (`ndstar_decomp.py`, histogram V `norm_all`)

| bin | kształt NDS (V MC, k=2,73) | kształt PD (MC par, k=2,11) |
|---|---|---|
| 0,05 MeV | 4436 ± 671 (threshold_pol2, χ² 336/383) | 4100 ± 639 (threshold_pol2, χ² 336/383) |
| 0,10 MeV | 4443 ± 675 | 4107 ± 641 |

- Kształt: −7,7 %. RMS sygnału 1,25 → 1,01 MeV, a σ szerokiego Gaussa 2,03 → 1,49 MeV. Skala s ≈ 1,03 wolna w obu przypadkach.
- Binning: +0,2 % (po poprawce NLL; było +0,5 %).
- Reszta 4097 → 4022 (−1,8 %): inny histogram wejściowy (482 264 vs 482 979 wejść, inne krawędzie binów przy LO) i inny model tła z kolejki prób.
- **Dane nie rozróżniają kształtów** (χ² 336 vs 336). Różnica to czysta systematyka źródła kształtu MC: |4436 − 4100| = **±336 (7,6 %)** (przed poprawką NLL: 4097, ±339).
  NDS jej nie zawierał, bo jego kategoria „sig” to tylko GCB i stała skala.
- Uwaga z B2b: `shape_from_mc` w oknie 143–148,5 ma closure −0,08 (8–30, WP80) i obcina ogony, co zaniża k. Oba kształty są tym dotknięte.
  Szersze okno (141–152) prawdopodobnie podniesie k, a więc i N. Bardziej wiarygodny wydaje się więc kształt z większym k (V/NDS).

## 2. Rekomendacja i niezmienniczość

**N_D* = 4436 +715/−645 (stat) ± 692 (syst) [605 NDS ⊕ 336 kształt MC], δN/N = 21,9 %** (stat symetryzowane 680; z HESSE 671: 21,7 %; przed poprawką NLL 4443 ± 612 ± 681, 20,6 %; 19,1 % bez składnika kształtu). Uzasadnienie:
1. Valid fit odpornym modelem (`fit_single`). Pełny skan wariantów, mediana w40 = 4339 (rozrzut 190).
2. Ta sama selekcja co mianownik `mc_rates` (okno ±40, Δm < 165, bez cięcia score).
3. Kształt z większym k jest mniej obcięty (por. closure B2b).
4. Nominał już był przeliczony na tę bazę (×1,257), więc wyniki zostają spójne z dokumentacją.
5. 4022 i 4220 mieszczą się w ±1σ syst.

**Niezmienniczość (dokładna, nie przybliżona).** W `PairFit.lam` mamy λ = μ·SS·(…) + κ_a·d_sb·SB + κ_b·d_bs·BS + bb·BB, gdzie:
- κ jest wolne per noga (planB), z ograniczeniem tylko κ ≥ 0, bez priora;
- d_sb i d_bs mają prior 1 ± 0,10 (względny, niezależny od skali);
- `partner` karze tylko log(κ_i/κ_j), co też nie zależy od skali.

Przeskalowanie N_D* → N_D*/r daje SS → SS/r i SB/BS → SB/BS/r. Wtedy μ → r·μ i κ → r·κ odtwarzają identyczne λ, czyli identyczny NLL.
Stąd μ̂·N_D* = const. To samo dotyczy błędów Hesse, przedziałów profilu i krzywej profilu (rozciągnięta w μ o r).
Wartość ΔNLL w minimum i w μ=0 się nie zmienia, więc **Z, q_shape, Δχ² testów `pythia_mc`/`dps_flat` i Z odrzucenia Pythii/DPS są niezmiennicze**.
f_LS = μ_LS·n_LS / (μ_LS·n_LS + μ_OS·n_OS) to stosunek, w którym r się skraca. Numerycznie: 0,2189 / 0,3447 / 0,4927 przed i po.
Jedyny wyjątek to numeryka: punkt startowy minimalizatora (κ₀ = 1) i obcięcie siatek profilu (§3).

## 3. Tabela μ na wspólnej bazie N_D* = 4443
Czynniki: nominał k = 5583,6/4442,5 = **1,2569**, PD k = 4021,6/4442,5 = **0,9052**. μ_ref = k·μ_użyte. Z bez zmian.

### Per kanał / bin

| | μ (baza biegu) | **μ (4443) ± Hesse** | 68 % profil (4443) | Z | S oczek. (4443, μ=1) |
|---|---|---|---|---|---|
| nominał D0D0 | 1,158 | **1,46 ± 0,69** | [0,80; 2,16] | 2,25 | 46,3 |
| nominał D⁺D0 | 1,165 | **1,46 ± 0,53** | [0,95; 1,99] | 3,01 | 72,0 |
| nominał D⁺D⁺ | 1,386 | **1,74 ± 1,16** | [0,61; 2,93] | 1,57 | 20,5 |
| PD D0D0 like | 0,958 | **0,87 ± 1,29** | [−0,37; 2,20] | 0,69 | 16,8 |
| PD D0D0 opposite | 1,954 | **1,77 ± 0,82** | [0,99; 2,62] | 2,37 | 29,4 |
| PD D⁺D0 like | 1,781 | **1,61 ± 1,13** | [0,55; 2,79] | 1,56 | 22,8 |
| PD D⁺D0 opposite | 1,572 | **1,42 ± 0,56** | [0,89; 2,00] | 2,83 | 49,1 |
| PD D⁺D⁺ like | 4,478 | **4,05 ± 4,05** | [0,18; ≥4,53]* | 1,05 | 4,1 |
| PD D⁺D⁺ opposite | 1,153 | **1,04 ± 1,04** | [0,06; 2,12] | 1,05 | 16,4 |

\* Górna granica jest obcięta przez zakres `FINE` (μ ≤ 5,0 w bazie biegu). D⁺D⁺ like jest praktycznie niemierzalne.

Spójność Σ_i μ_i f_i / μ_nom (oba na 4443): **D0D0 0,99, D⁺D0 1,01, D⁺D⁺ 0,95**. Było 1,374 / 1,407 / 1,312 = k_nom/k_PD ≈ 1,388.

### μ globalne

„JSON” to wartość z raportu (minimum na węźle siatki 0,25). „gładkie” to minimum spline sumy profili (`ndstar_rescale.py`); parabola przez 3 i 5 węzłów daje to samo w granicach ±0,01.

| wielkość | baza biegu: JSON | baza biegu: gładkie | **4443: gładkie** | 4443: JSON×k | Z |
|---|---|---|---|---|---|
| nominał (3 kanały) | 1,25 [0,90; 1,51] | 1,19 [0,88; 1,51] | **1,49 [1,10; 1,90]** | 1,57 [1,13; 1,90] | 4,06 (4,07) |
| nominał + syst profil. | 1,25 [0,92; 1,54] | — (brak profili w JSON) | ≈ 1,49 (przesunięcie jak wyżej) | 1,57 [1,15; 1,93] | 4,11 |
| PD wspólne (kanały × biny, `pythia_mc`) | 1,50 [1,20; 2,04] | 1,61 [1,20; 2,03] | **1,45 [1,08; 1,84]** | 1,36 [1,09; 1,84] | 4,21 |
| PD `dps_flat` | 1,64 [1,22; 2,09] | — | — | 1,48 [1,10; 1,89] | 4,07 |
| PD like (biny łącznie) | 1,50 [0,69; 2,53] | 1,58 [0,69; 2,53] | **1,43 [0,62; 2,29]** | 1,36 | 1,83 |
| PD opposite | 1,50 [1,16; 2,09] | 1,61 [1,16; 2,09] | **1,46 [1,05; 1,89]** | 1,36 | 3,79 |

Wniosek: **nominał 1,49 i PD 1,45 na tej samej bazie są zgodne** (Δ ≈ 3 %). Po przeliczeniu „μ_PD 1,50 vs μ_nom 1,25” znika całkowicie.

**Proponowany zapis wyniku (do decyzji):** μ = **1,50 +0,41/−0,39 (stat) ± 0,33 (N_D*, 21,9 %)**, Z = 4,06σ (4,11 z syst kształtu).
Zastępuje „1,57 ± 0,38 ± 0,30”. Składniki:
- stat: przedział z gładkiego profilu × 1,257;
- 0,33 = 1,495 × 0,219 (przed poprawką NLL: 0,31 = 1,49 × 0,206; NDS 0,30 było 1,57 × 0,191).

Inne bazy, μ_glob (gładkie) = 6633/N_D* (nominał), 6459/N_D* (PD):

| N_D* | 5584 | 4443 | **4436** | 4220 | 4097 | 4022 |
|---|---|---|---|---|---|---|
| nominał | 1,19 | 1,49 | **1,50** | 1,57 | 1,62 | 1,65 |
| PD | 1,16 | 1,45 | **1,46** | 1,53 | 1,58 | 1,61 |

## 4. Co wymagałoby ponownego liczenia (nic z tego nie jest potrzebne do tabel wyżej)
1. **Poprawka kształtu MC** (okno `shape_from_mc` 141–152, wspólny kształt dla NDS i potoku par): tylko refit Δm, czyli `pairs/build.py:n_dstar_data` /
   `dstar_refit.fit_single`. Wejścia: Δm danych w oknie ±40 i Δm MC truth D*. Lokalnie wystarczy V `dstar_nom.npz` (dane), ale potrzebny jest
   **niebinowany** Δm MC lub histogram truth z szerszym oknem. μ przeliczyć potem tym samym czynnikiem, bez ponownego fitu par.
2. **N_D* per bin pT** (suma per bin 4171 NDS / B2b): wtedy n_pair = Σ_b N_D*(b)·r_b(MC) zależy od rozkładu pT i normalizacja **nie jest już jednym
   czynnikiem**. Każdy kanał dostaje inny czynnik (liniowy, ale per kanał). Trzeba przebudować `pair_inputs` (`mc_rates` per bin pT) i powtórzyć fit2d.
   Kroki: `pairs build`, potem `fit2d --unblind [--diff]`.
3. **Kwantyzacja siatki**: w `fit2d.combine`/`diff.combine_bins` profil liczony jest na siatce 0,25 z interpolacją liniową. Albo zagęścić `GRID` (0,05 wokół
   minimum, ponowny `profile` bez nowego builda), albo raportować minimum z paraboli/spline (jak tu). Dotyczy też przedziałów 68 % (±0,01–0,02).
   Profile `syst_profiled` nie są zapisane w JSON, więc ich gładkiego minimum nie da się policzyć lokalnie.
4. Błąd f_MC (±0,07–0,10, z 8–28 par MC) nadal nie jest wliczony w μ_i. To niezależne od N_D*.

## 5. Do przeniesienia (poza tym katalogiem — nie edytowane tutaj)
- `charm-PD-pairs-diff.md` §9 „Normalizacja”: zastąpić „nie porównywalne” zdaniem „na bazie N_D* = 4443: μ_PD = 1,45 [1,08; 1,84] vs μ_nom = 1,49 [1,10; 1,90];
  spójność 0,99/1,01/0,95” i dodać link do `ndstar/ndstar.md`.
- Nagłówek wyniku: rozważyć zmianę „μ = 1,57 ± 0,38 ± 0,30” na „1,50 +0,41/−0,39 ± 0,33” (siatka + syst kształtu MC). Interpretacja σ_eff nieznacznie się przesuwa (μ ×0,95).
- `NDS_ndstar` syst: dopisać kategorię „źródło kształtu MC” ±336 (rerun 2026-11-12).
