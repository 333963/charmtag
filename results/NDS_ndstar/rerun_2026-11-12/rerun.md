# N_D* — ponowny skan NDS po poprawce NLL w `dstar_refit.fit_single` (2026-11-12)

> **Nota 2026-11-13:** od 13.11 stat +0,48/−0,46 (× 1,18 z jackknife po zdarzeniach; Z = 4,06σ bez zmian, ≈ 3,4σ z korelacjami zdarzeń), zob. `PF2_lxplus/runs/lxplus/unblind/UNB_WYNIK.md` §1b.

Liczby ATLAS — do użytku wewnętrznego. Oryginały (`ndstar_full.json/.md`, `ndstar.md`, AA_merge, REPO) nietknięte.

## TL;DR
- **N_D* = 4436 ± 671 (stat, HESSE; MINOS −645/+715) ± 692 (syst) [NDS 605 ⊕ kształt MC 336]**.
- **δN/N = 21,7 % (HESSE) / 21,9 % (MINOS sym. 680)** — było 20,6 %.
- **μ = 1,495 (≈1,50) +0,41/−0,39 (stat) ± 0,33 (N_D*, 21,9 %)** — było 1,49 … ± 0,31. Przesunięcie μ centralnego +0,14 % (tylko zaokrąglenie przechodzi z 1,49 na 1,50).
- Zmiana pochodzi niemal wyłącznie z błędu stat (612 → 671, +10 %); centralne N i syst zmieniają się w granicach szumu numerycznego starego fitu.
- Wszystkie warianty valid (było `all_valid: false` przez `bin_0.10`). Żaden wariant nie przesunął się o > 0,1σ (67).

## Pliki (wszystko w tym katalogu)
| plik | co |
|---|---|
| `ndstar_variants.py` | kopia skryptu NDS + przełącznik `NLL_MODE=old/new` (opis niżej) i zapis pośredni `partial.jsonl` |
| `old_repro/` (+`old_repro.log`) | odtworzenie: AA_merge (stary `fit_single`) + `NLL_MODE=old` |
| `new/` (+`new.log`) | **nowy wynik**: charmtag z poprawką NLL z 10.11 + `NLL_MODE=new` |
| `nominal_minos.py` → `nominal_minos.json/.log` | nominał: HESSE, MINOS, skan profilu (33 pkt, ±1,6σ) |
| `jitter_test.py` → `jitter_old/new.json` | test szumu numerycznego (przesunięcie krawędzi binów o ±1e-9 MeV) |
| `ndstar_decomp_rerun.py` → `ndstar_decomp_new.json` / `_old.json` | kopia `ndstar_decomp.py` (dstar_refit z PYTHONPATH, wynik do pliku z argv) |

Wejścia: `V_daneMC/lxplus/dstar/dstar_nom.npz` (`norm_all`, `norm_<bin>`, MC sygnału, okna w36/w55), `V_daneMC/lxplus/summary/dstar_post.json`
(`norm.shape_mc` = identyczny z `shape_2g_mc` w `ndstar_full.json`). Środowisko: `.venv_h5`, numpy 2.5.4, scipy 1.18.1, iminuit 2.33.0.

## Zmiany w kopii `ndstar_variants.py`
Skrypt ma dwa własne fity z `Minuit` liczące NLL po staremu, Σ[λ − h·ln λ]:
- `fit_generic` (l. 42–43 oryg.) — warianty: tło `sqrt_pol3`, `pure_expo`, sygnał G+CB, skala stała;
- `gcb_shape_from_mc` (l. 64–65 oryg.) — fit kształtu G+CB do MC (wartości NLL mniejsze, ale ta sama postać).

W obu NLL zastąpione przez `_pnll(lam, hh)`: `NLL_MODE=new` → Σ[λ − h − h·ln(λ/h)] (h=0 → ln λ/1, jak w REPO `hs = where(h>0, h, 1)`),
`NLL_MODE=old` → stara postać. Stała addytywna — minimum się nie zmienia. Poza tym: zapis każdego wariantu do `partial.jsonl` w locie oraz pola
`nll_mode`, `dstar_refit_file` w JSON. Warianty z modelem tła A/B i sygnałem 2G idą przez `R.fit_single` (wersję wybiera PYTHONPATH).
Kształt GCB po poprawce: praktycznie identyczny (sc 0,63484→0,63487, k 1,71519→1,71523).

## 1. Odtworzenie starego wyniku
| | NDS `ndstar_full.json` (lxplus) | stary fit lokalnie (`old_repro`) | 
|---|---|---|
| N_nominal | 4442,5 ± 611,6 | 4436,3 ± 595,2 |
| N_syst | 590,7 | 587,6 |
| syst bkg / range / bin / sig | 520,9 / 202,3 / 100,1 / 163,4 | 527,3 / 196,0 / 94,8 / 140,8 |
| mediana / rozrzut w40 | 4345 / 195 | 4390 / 190 |
| all_valid | false (bin_0.10) | false (sig_gauss_cb) |

**4442,5 nie odtwarza się dokładnie — ale wejścia są te same, a różnica to szum numeryczny starego fitu.** Dowody:
1. Warianty, w których fit jest dobrze uwarunkowany, zgadzają się do 0,1–1: `bkg_sqrt_pol3` 4963,5/4963,6, `bkg_pure_expo` 5552,7/5552,7,
   `range_lo141` 4316,9/4316,8, `window_d0_36` 3998,9/3998,8, `window_d0_55` 5187,4/5187,5, bin 5–8 1352,6/1352,6. Identyczne n_entries (482 979) i χ² nominału (335,81).
2. `jitter_test`: przesunięcie krawędzi o ±1e-9 MeV (fizycznie zero) daje w starym ficie N = 4436,3 / 4445,4 / 4442,8 i HESSE = 595 / 626 / 587.
   Nowy fit: 4436,39 / 4436,39 / 4436,39, HESSE 671 / 673 / 675. Stary fit kończy z EDM 5e-4…1e-3 (tuż pod progiem), nowy 4e-6.
3. ΔNLL między 4442,5 a 4436,4 wynosi ~4e-5 (profil), poniżej rozdzielczości starego NLL (~−3·10⁶).
Zatem 612 (NDS) vs 588 (wcześniejsze sprawdzenie) vs 595 (tu) — **to ten sam szum HESSE starego fitu**; zależy od wersji numpy/kolejności
sumowania, nie od innego histogramu. Wersji numpy na lxplus z czasu NDS nie da się sprawdzić lokalnie (nie zapisana) — to jedyna niesprawdzona część.
(Przesunięcia ±1e-7/1e-6 MeV zmieniają liczbę binów na brzegu zakresu przez tolerancję 1e-9 w `_rebin` → N 4412 w obu wersjach; to nie jest test szumu.)

## 2. Nominał: HESSE vs MINOS (nowy fit)
| | wartość |
|---|---|
| N | 4436,39 (valid, EDM 4e-6, żaden parametr nie na limicie) |
| HESSE | ±670,9 |
| MINOS | −645,0 / +715,1 (oba valid) |
| skan profilu ΔNLL=0,5 (niezależny) | −644,1 / +714,4 |

MINOS potwierdzony skanem. Asymetria ±5 %. **Rekomendacja:** podawać N_D* stat jako **−645/+715 (MINOS)**; do symetrycznej
składowej normalizacyjnej μ brać średnią 680 (różnica względem HESSE 671 to 0,1 p.p. w δN/N). Uwaga: μ ∝ 1/N, więc asymetria odwraca się
(N w górę → μ w dół: −13,9 % / +17,0 % dla samego stat N_D*).

## 3. Stare vs nowe (NDS lxplus → nowy fit)
| wielkość | NDS (stary) | **nowy** | Δ |
|---|---|---|---|
| N_nominal | 4442,5 | **4436,4** | −6,2 (−0,01σ) |
| N_stat (HESSE) | 611,6 | **670,9** | +59 (+9,7 %) |
| N_syst (NDS) | 590,7 | **605,0** | +14 |
| syst bkg | 520,9 | 533,8 | +13 |
| syst range | 202,3 | 210,6 | +8 |
| syst bin | 100,1 | 106,9 | +7 |
| syst sig | 163,4 | 158,8 | −5 |
| syst z „czystym wykładnikiem” | 1257,6 | 1270,3 | +13 |
| rigid_expo N / Δχ² / skala | 5552,7 / 5,11 / 1,173 | 5553,3 / 5,12 / 1,173 | — |
| mediana w40 | 4345,2 | 4338,7 | −6,5 |
| rozrzut w40 | 195,2 | 190,1 | −5 |
| n wariantów w40 (valid) | 10 | 11 | `bin_0.10` teraz valid |
| all_valid | false | **true** | |
| suma per bin pT | 4170,7 | 4158,4 | −12 |
| kształt MC (decomp, |NDS − PD kształt|, bin 0,05) | 339 (4436 − 4097) | **336 (4436,4 − 4100,3)** | −3 |
| syst łączny (NDS ⊕ kształt) | 681 | **692** | +11 |
| δN/N | 20,6 % | **21,7 % (HESSE) / 21,9 % (MINOS 680)** | +1,1…1,3 p.p. |

### Warianty (N ± HESSE, valid, model tła)
| wariant | NDS (stary) | nowy | ΔN | uwagi |
|---|---|---|---|---|
| nominal | 4442,5 ± 612 ✔ thr_pol2 | 4436,4 ± 671 ✔ thr_pol2 | −6,2 | |
| bkg_expo_pol2_fullmodel | 4442,2 ± 567 ✔ (posdef ✖) | 4432,1 ± 683 ✔ (posdef ✔) | −10,0 | posdef ✖→✔ |
| bkg_sqrt_pol3 | 4963,5 ± 609 ✔ | 4970,2 ± 660 ✔ | +6,8 | |
| bkg_pure_expo | 5552,7 ± 612 ✔ | 5553,3 ± 623 ✔ | +0,6 | |
| range_hi158 | 4240,3 ± 615 ✔ | 4225,7 ± 684 ✔ | −14,5 | |
| range_hi165 | = nominal | = nominal | −6,2 | npz kończy się na 160 MeV → identyczny z nominałem (jak w NDS) |
| range_lo141 | 4316,9 ± 581 ✔ | 4316,5 ± 667 ✔ | −0,4 | |
| range_lo141_hi165 | = range_lo141 | = range_lo141 | −0,4 | j.w. |
| bin_0.10 | 4470,7 ± 608 **✖** | 4443,2 ± 675 **✔** | −27,5 | **valid ✖→✔** (największe ΔN, 0,04σ) |
| bin_0.20 | 4342,5 ± 603 ✔ | 4329,5 ± 664 ✔ | −12,9 | |
| sig_gauss_cb | 4279,2 ± 640 ✔ | 4277,6 ± 714 ✔ | −1,6 | |
| sig_scale_fixed | 4347,9 ± 406 ✔ | 4338,7 ± 406 ✔ | −9,2 | |
| window_d0_36 | 3998,9 ± 554 ✔ | 4000,4 ± 613 ✔ | +1,5 | |
| window_d0_55 | 5187,4 ± 878 ✔ | 5163,6 ± 899 ✔ | −23,9 | |

**Żaden wariant nie zmienił N o > 0,1σ (67).** Zmiana statusu valid: tylko `bin_0.10` (✖→✔). Model tła: bez zmian we wszystkich wariantach.
Posdef: `bkg_expo_pol2_fullmodel` ✖→✔. Wszystkie nowe fity: EDM ≤ 1e-4.

### Per bin pT
| bin | NDS | nowy | ΔN |
|---|---|---|---|
| 3–5 | 643,8 ± 305 ✔ | 640,5 ± 282 ✔ (s=0,80 — na dolnym limicie skali) | −3,3 |
| 5–8 | 1352,6 ± 302 ✔ | 1344,7 ± 303 ✔ | −7,9 |
| 8–30 | 2174,3 ± 175 ✔ | 2173,1 ± 192 ✔ | −1,2 |

### Dekompozycja kształtu MC (`ndstar_decomp_new.json`, stary `_old.json` = odtwarza `ndstar_decomp.json` co do 0,1)
| bin | kształt NDS (stary → nowy) | kształt PD (stary → nowy) |
|---|---|---|
| 0,05 | 4436,3 ± 595 → **4436,4 ± 671** | 4097,4 ± 674 (expo_pol2_s1, po 2 invalid) → **4100,3 ± 639 (threshold_pol2, valid od razu)** |
| 0,10 | 4457,6 ± 680 → 4443,2 ± 675 | 4081,3 ± 518 → 4107,1 ± 641 |

Składnik „źródło kształtu MC”: 339 → **336** (7,6 %, bez zmian w praktyce). Binning przy kształcie PD: −0,4 % → +0,2 %.

## 4. Przeliczenie μ
- δN/N = hypot(stat, hypot(605,0; 336,1)) / 4436,4 = hypot(stat, 692,0) / 4436,4:
  HESSE 671 → **21,73 %**; MINOS sym. 680 → **21,87 %**; (MINOS górny 715 → 22,4 %, dolny 645 → 21,3 %).
- μ centralne = 6633/4436,4 = **1,4951** (było 6633/4442,5 = 1,4931). Przedział stat ×(4442,5/4436,4 = 1,0014): [1,10; 1,90] → [1,10; 1,90] — **+0,41/−0,39 bez zmian**.
- Składnik normalizacyjny: 1,4951 × 0,2187 = **0,327 → ±0,33** (MINOS sym.); z HESSE 1,4951 × 0,2173 = 0,325 → ±0,32. Było 0,31.
- PD wspólne: 6459/4436,4 = 1,456 (było 1,45) — przesunięcie tym samym czynnikiem 1,0014.

**Proponowany zapis:** N_D* = 4436 +715/−645 (stat) ± 692 (syst) [NDS 605 ⊕ kształt MC 336], δN/N = 21,9 %;
**μ = 1,50 +0,41/−0,39 (stat) ± 0,33 (N_D*, 21,9 %)**, Z = 4,06σ (4,11 z syst kształtu) — Z bez zmian (niezmienniczość, `ndstar.md` §2).
Uwaga: 1,4951 leży na granicy zaokrąglenia; jeśli wolisz nie ruszać wartości centralnej, „1,49” z 3 cyframi to 1,495 — różnica nieistotna.

## 5. Proponowane zmiany (NIE wprowadzone)

### `charm_work/PD_pairs_diff/ndstar/ndstar.md`
| linia | stary tekst | nowy tekst |
|---|---|---|
| 10 | `**Rekomendacja: N_D* = 4443 ± 612 (stat) ± 681 (syst) ≈ ±21 %.** Syst = NDS 591 ⊕ nowy składnik „źródło kształtu MC” ±339 (§1.3).` | `**Rekomendacja: N_D* = 4436 +715/−645 (stat, MINOS) ± 692 (syst) ≈ ±22 %.** Syst = NDS 605 ⊕ „źródło kształtu MC” ±336 (§1.3). Po poprawce NLL w `fit_single` (10.11), zob. `NDS_ndstar/rerun_2026-11-12/rerun.md`.` |
| 19 | `**Na bazie 4443: nominał μ = 1,49 [1,10; 1,90], PD μ = 1,45 [1,08; 1,84].**` | `**Na bazie 4436: nominał μ = 1,50 [1,10; 1,90], PD μ = 1,46 [1,08; 1,84].**` (opcjonalne; zmiana ×1,0014) |
| 27 | `\| **4443 ± 612 ± 591** \| NDS: …` oraz w tym wierszu `valid, χ² 336/383; 14 wariantów; syst = bkg 521, range 202, bin 100, sig 163.` | `\| **4436 ± 671 ± 605** (MINOS −645/+715) \| NDS rerun 2026-11-12 (NLL wzgl. modelu nasyconego): …` oraz `valid, χ² 336/383; 14 wariantów, wszystkie valid; syst = bkg 534, range 211, bin 107, sig 159.` |
| 46 | `\| 0,05 MeV \| 4436 ± 595 (threshold_pol2, χ² 336/383) \| 4097 ± 674 (expo_pol2_s1, χ² 336/382) \|` | `\| 0,05 MeV \| 4436 ± 671 (threshold_pol2, χ² 336/383) \| 4100 ± 639 (threshold_pol2, χ² 336/383) \|` |
| 47 | `\| 0,10 MeV \| 4458 ± 680 \| 4081 ± 518 \|` | `\| 0,10 MeV \| 4443 ± 675 \| 4107 ± 641 \|` |
| 50 | `- Binning: +0,5 %.` | `- Binning: +0,2 % (po poprawce NLL; było +0,5 %).` |
| 52 | `\|4436 − 4097\| = **±339 (7,6 %)**.` | `\|4436 − 4100\| = **±336 (7,6 %)**.` |
| 59 | `**N_D* = 4443 ± 612 (stat) ± 681 (syst) [591 NDS ⊕ 339 kształt MC], δN/N = 20,6 %** (było 19,1 % bez składnika kształtu).` | `**N_D* = 4436 +715/−645 (stat) ± 692 (syst) [605 NDS ⊕ 336 kształt MC], δN/N = 21,9 %** (stat sym. 680; z HESSE 671: 21,7 %; przed poprawką NLL 20,6 %).` |
| 60 | `mediana w40 = 4345 (rozrzut 195).` | `mediana w40 = 4339 (rozrzut 190), 11 wariantów valid.` |
| 113 | `μ = **1,49 +0,41/−0,39 (stat) ± 0,31 (N_D*, 20,6 %)**` | `μ = **1,50 +0,41/−0,39 (stat) ± 0,33 (N_D*, 21,9 %)**` |
| 116 | `- 0,31 = 1,49 × 0,206 (NDS 0,30 było 1,57 × 0,191).` | `- 0,33 = 1,495 × 0,219 (przed poprawką NLL: 0,31 = 1,49 × 0,206; NDS 0,30 = 1,57 × 0,191).` |
| 120–122 | kolumna `4443` z wartościami `1,49` / `1,45` | kolumna `4436` z `1,50` / `1,46` (opcjonalnie) |
| 140 | `na „1,49 +0,41/−0,39 ± 0,31”` | `na „1,50 +0,41/−0,39 ± 0,33”` |
| 141 | `dopisać kategorię „źródło kształtu MC” ±339.` | `dopisać kategorię „źródło kształtu MC” ±336 (rerun 2026-11-12).` |

Tabele §3 (μ per kanał na bazie 4443) — opcjonalnie przeskalować ×1,0014; w zaokrągleniu do 0,01 zmieniają się co najwyżej o 0,01.

## 6. Czego nie sprawdzono / zastrzeżenia
- Dokładne 4442,5 ± 611,6 z lxplus nie jest odtwarzalne lokalnie (szum starego NLL; wersja numpy na lxplus nieznana). Wejścia uznane za identyczne na podstawie §1.
- AA_merge nadal ma starą wersję `fit_single` — każdy nowy bieg z AA_merge/Q_lxplus da zaszumiony HESSE (±~40 na nominale).
- `gcb_shape_from_mc` zmieniony tak samo; wpływ na kształt GCB pomijalny (5. cyfra).
- 6633 (stała μ·N z `ndstar.md`) przepisana z dokumentu, nieprzeliczana; przedział stat μ przeskalowany analitycznie, profile par nie były liczone.
