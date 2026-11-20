# `results/` — wyniki projektu A z liczbami

Stan na 20.11.2026. Raporty, tabele, małe pliki JSON z liczbami, nota i prezentacja. Skrypty, które je policzyły, są w
[`../analysis/`](../analysis/README.md); ścieżki w obu katalogach są takie same jak w `~/Documents/charm/charm_work/`
(`charm_work/<ścieżka>` → `results/<ścieżka>` albo `analysis/<ścieżka>`).

Pliki są kopiami stanu z 20.11.2026.
Raporty opisują stan z dnia powstania. Stan bieżący: `NOTA_WYNIK/nota.md` (19.11), `PF2_lxplus/runs/lxplus/unblind/UNB_WYNIK.md`
§0 i `truthnorm/TRUTHNORM.md`. Wcześniejsze raporty (ślepy bieg PF2, odsłonięcie 5.11, budżet 16.11) zostają jako historia.

Skrót w tabelach: `U/` = `PF2_lxplus/runs/lxplus/unblind/`.

## Wynik główny — run 427929 (czysty minimum bias, HLT_mb_mbts_L1MBTS_1)

> μ = 0,50 +0,41/−0,37 (stat, × 1,13 jackknife) +0,28/−0,19 (syst) ± 0,21 (N_D*);
> Z = 1,60 względem braku par (≈ 1,42 z korelacjami zdarzeń); względem Pythii 8 A3 (μ = 1): −0,93σ;
> σ_eff > ~6 mb (95 % CL, przy f_DPS = 0,40).

μ = [N_par / N_D*]_dane / [N_par / N_D*]_MC, normalizacja par z prawdy mb_nd. Zatwierdzony przez autora 19.11.2026.
Trzy znane niespójności tego wyniku (N_D* 1901 ↔ 2046, `tag_sf`/`sf_json`, skala szerokości piku D* 1,42) są opisane
w [`../PRZEGLAD.md`](../PRZEGLAD.md#znane-problemy); liczby w tym katalogu **nie** zostały z ich powodu zmienione.

### Indeks: od liczby do pliku i skryptu

Pliki względem `results/`, skrypty względem `analysis/`.

| liczba | wartość dokładna | plik → klucz | skrypt |
|---|---|---|---|
| μ wspólne | 0,49894733917589 | `U/truthnorm/res/tn_nom427.json` → `fit.combined` (μ, przedział 68 %, Z); to samo w `U/truthnorm/truthnorm.json` → `nom427` | `U/truthnorm/truthnorm_run.py --what nom427` |
| przedział 68 % z profilu | [0,17454; 0,86467] | jw. | jw. |
| Z = 1,60 | 1,60262 | jw.; `truthnorm.json` → `write_427929.Z` | jw. |
| μ per kanał D⁰D⁰ / D⁺D⁰ / D⁺D⁺ | 0,529 / 0,606 / −0,067 | `tn_nom427.json` → `fit.<kanał>.mu_hat` | jw. |
| czynniki f_ch (prawda mb_nd / direct) | 1,0317 / 1,2659 / 1,0597 | `GEN_mcstat/out/truth_mbnd.json` → `truth_self_over_direct`; kopia w `tn_nom427.json` → `f` | `GEN_mcstat/truth_mbnd.py` |
| n_pair przed i po × f_ch | 30,42 / 44,22 / 21,14 → 31,39 / 55,97 / 22,40 | `tn_nom427.json` → `n_pair_nominal`, `n_pair`; źródło: `U/r427929/pairs/pair_inputs.json` → `n_pair` | `charmtag pairs --unblind` (lxplus, snapshot charmtag 2.0.3 z 10.11), config `U/r427929/pairs/cfg427.yaml` |
| histogramy par i szablony (wejście fitu) | — | `U/r427929/pairs/pairs_full.npz.json`, `pair_inputs.npz.json`, `pair_inputs.json` | jw.; zadania `Q_lxplus/ub_427929/u4_job.sh` |
| czynnik jackknife 1,13 | 1,133895 | `U/r427929/jk_ev200/ubjk.json` → `per_mode.fixed.combined.ratio_68`; tabela `ubjk.md` | `Q_lxplus/ub_jk/ubjk_run.py` (zadania `Q_lxplus/ub427_jk/`); w zapisie użyto stałej `JK427 = 1.13` z `truthnorm_report.py` |
| stat +0,41 / −0,37 | +0,41327 / −0,36657 | `U/truthnorm/truthnorm.json` → `write_427929.stat_up`, `stat_dn` = (granica 68 % − μ) × 1,13 | `U/truthnorm/truthnorm_report.py` |
| syst +0,28 / −0,19 | +0,27506 / −0,18677 | `truthnorm.json` → `write_427929.syst_up`, `syst_dn` = μ × (0,5513 / 0,3743) | jw. |
| syst względny +0,551 / −0,374 | 0,551281 / 0,374330 | `U/r427929/budzet_427929.json` → `syst_rel_up`, `syst_rel_dn`; wiersze: `rows`, grupy: `groups`; tabela `BUDZET_427929.md` | `U/r427929/budzet_427929.py` |
| N_D* ± 0,21 | 0,20656 | `truthnorm.json` → `write_427929.nds` = μ × 0,4140 | `truthnorm_report.py` |
| N_D* względnie 41,4 % (stat 23,7 % ⊕ warianty 33,1 % ⊕ kształt MC 7,6 %) | 0,413997 | `budzet_427929.json` → `ndstar`; `U/r427929/nds/ndstar_variants_427929.json` → `summary.rel` | `Q_lxplus/ub427_nds/n8_ndstar.py` |
| N_D* = 2046 ± 485 (stat) ± 677 (14 wariantów) | 2045,52 ± 484,90 ± 676,76 | `U/r427929/ndstar_427929.json` → `fit`; `nds/ndstar_variants_427929.json` → `summary`; `nds/NDS_427929.md` | `Q_lxplus/ub_427929/u4_ndstar.py`, `Q_lxplus/ub427_nds/n8_ndstar.py` |
| N_D* użyte w ficie par | 1901,49 ± 444,53 | `U/r427929/pairs/pair_inputs.json` → `meta.dstar_fit_data` | `charmtag pairs` (`dstar_refit.fit_single` w wersji charmtag 2.0.3 z 10.11) |
| skala szerokości piku D* | 1,4204 (bieg par), 1,4043 (NDS) | `pair_inputs.json` → `meta.dstar_fit_data.scale`; `ndstar_427929.json` → `fit.scale`; `U/r427929/ndsw/NDSW.md` | jw.; `Q_lxplus/ub_ndsw/ndsw.py` |
| `tag_sf` = 0,7, `sf_json` = null | — | `pair_inputs.json` → `meta.config` | config `U/r427929/pairs/cfg427.yaml` |
| wagi zdarzeń c_w(427929) | D⁰D⁰ 0,977, D⁺D⁰ 0,762, D⁺D⁺ 0,996 | `U/r427929/cw/cw.json` → `c_w_427929` | `Q_lxplus/ub_427929/u4_cw.py` |
| Z ≈ 1,42 z korelacjami | 1,60262 / 1,13 | `U/truthnorm/TRUTHNORM.md` (sekcja „Zapis wyników”) | `truthnorm_report.py` |
| względem Pythii −0,93σ | −0,93185 | `truthnorm.json` → `write_427929.vs_pythia_sigma` | `truthnorm_report.py` |
| σ_eff > ~6 mb (95 % CL) | — | `U/r427929/SIGMA_EFF.md` (sekcja „Aktualizacja 19.11 wieczór”) | `U/r427929/seff_final.py` |
| górna strona tune (fragmentacja ALICE), +40,9 % | k = 0,65 / 0,73 / 0,81 | `U/frag/FRAG.md`; `budzet_427929.json` → `tune_up` | `U/frag/frag.py` |
| kontrola: normalizacja direct, μ = 0,567 (fit gładki) / 0,527 (zapis UB427 na N_D* 2046) | 0,567396 / 0,527445 | `truthnorm.json` → `nom427_f1`; `budzet_427929.json` → `mu`; `U/r427929/UB427.md`, `ub427.json` | `truthnorm_run.py --what nom427 --f-one`; `Q_lxplus/ub_427929/u4_report.py` |
| fit kanałów z lxplus (siatka profilu, μ wspólne 0,50 [0,20; 1,00]) | — | `U/r427929/pairs/fit2d_unblind.json` → `combined` | `charmtag fit2d --unblind` |
| sieć par jako kontrola schematu | μ = 0,10 [−0,11; 0,40] | `U/r427929/pnet/UBPNET427.md`, `ubpnet.json` | `Q_lxplus/ub_pnet/upn_fit.py` |
| propozycja: stat. MC par z prawdy mb_nd (0,162 zamiast 0,316) — nie przyjęta | — | `U/r427929/mcstat_truth_prop.json`; końcówka `BUDZET_427929.md` | `U/r427929/mcstat_truth_prop.py` |

Jak sprawdzić pierwszy wiersz na własnym komputerze: `python analysis/_odtworzenie/odtworz_wynik_glowny.py` (zob. `analysis/README.md`).

## Run 428770 — obserwacja obciążona triggerem L1_J12 (pomiar pomocniczy)

> μ = 1,27 +0,41/−0,40 (stat, × 1,18) ± 0,36 (syst, prawda mb_nd) ± 0,28 (N_D*); Z = 4,05 (≈ 3,43 z korelacjami zdarzeń).
> Zapis direct z 16.11: 1,50 +0,48/−0,46 ± 0,61 ± 0,33. Odsłonięcie 5.11: Z = 4,06σ.

| liczba | plik → klucz | skrypt |
|---|---|---|
| μ = 1,27, stat, syst, N_D*, Z = 4,05 | `U/truthnorm/truthnorm.json` → `write_428770` (fit: `res/tn_nom428.json`, przeskalowany z N_D* 5584 na 4436: `nom428_scaled`); `U/truthnorm/TRUTHNORM.md` | `U/truthnorm/truthnorm_run.py --what nom428`, `truthnorm_report.py` |
| odsłonięcie 5.11: Z = 4,06, μ per kanał | `U/pairs/fit2d_unblind.json` → `combined`, `channels`; `U/UNB_WYNIK.md` §1, `UNB_WYNIK.json`; rysunki `U/report_UNBLIND_*.png` | `charmtag pairs/fit2d --unblind`; zadania `Q_lxplus/pf_pairs/unblind_job.sh`, raport `Q_lxplus/pf_pairs/tools/fill_unblind.py` |
| wejścia fitu 428770 | `U/pairs/pair_inputs.json`, `pair_inputs.npz.json`, `pairs_full.npz.json` | jw. |
| zapis direct 1,50 +0,48/−0,46 ± 0,61 ± 0,33 | `U/UNB_WYNIK.md` §1b i §7; budżet: `VER_budget/RAPORT.md`, `VER_budget/RAPORT_gen_norm.md`, `TUNE_pythia/budget_with_tune.md` | `TUNE_pythia/budget_with_tune.py`, `VER_budget/kcls_check.py` |
| N_D* = 4436 ± 671 (NDS) | `NDS_ndstar/rerun_2026-11-12/` (`rerun.md`, `new/ndstar_full.json`) | `NDS_ndstar/rerun_2026-11-12/` w `analysis/` |
| czynnik jackknife 1,18 | `U/jk_ev200/ubjk.md`, `ubjk.json`; po plikach: `U/jk/` | `Q_lxplus/ub_jk/ubjk_run.py` |
| μ per klasa triggera (J12only 1,63; prawda mb_nd 1,37) | `U/ubtrig/UBTRIG.md`, `ubtrig.json`; `truthnorm.json` → `cls_*`, `J12_vs_rest` | `Q_lxplus/ub_trig/ubtrig_run.py`; `truthnorm_run.py --what cls:<i>` |
| test zerowy klas (24 ziarna) | `U/ubtrig/null/UBTRIG_NULL.md` | `Q_lxplus/ub_trig/ubtrig_null_agg.py` |
| miara triggera d = 1,04 ± 0,32 (prerejestracja) | `TRIG_A/SYST_TRIG.md`, `TRIG_A/syst_trig_out/` | `TRIG_A/syst_trig.py` |
| emulacja L1_J12 na generatorze, wydajność J12 z danych, SF w kategoriach par | `TRIG_A/trigmc/`, `TRIG_A/j12eff/`, `TRIG_A/pairsf/` | `TRIG_A/` w `analysis/` |
| topologia \|Δη\|: λ̂ = −1,11 | `U/ubtopo/UBTOPO.md` | `Q_lxplus/ub_topo/topo_report.py`, `topo_gen.py` |
| udział par tego samego znaku f_LS = 0,33 ± 0,14 | `U/fls/FLS.md` (zestawienie źródeł); `PD_pairs_diff/` | `charmtag pairs/fit2d --unblind --diff sign` (moduł `charmtag/pairs/diff.py`) |
| sieć par (428770) | `U/ubpnet/UBPNET.md`, `U/ubpnet/ubpnet.json` | `Q_lxplus/ub_pnet/upn_fit.py` |
| systematyka tune Pythii | `TUNE_pythia/tune_syst.md`, `TUNE_pythia/out/` | `TUNE_pythia/tune_syst.py` |

## Katalogi

| katalog | zawartość |
|---|---|
| `NOTA_WYNIK/` | `nota.md` — szkic noty wynikowej (19.11); `INNE_RUNY_MB.md`, `STRUMIEN_428770_v11.md` — notatki o próbkach |
| `PREZ_2026-11-19/` | tabele `p1/`, `p3/`, `p6/`; `pzdata/early/` |
| `PF2_lxplus/` | bieg par na lxplus: ślepy (`runs/lxplus/pairs/`, `pairs_syst/`, `mix/`, `TFIN_report.md`) i odsłonięty (`runs/lxplus/unblind/`); opis biegu `README_PF2.md` |
| `PF2_lxplus/runs/lxplus/unblind/` | `UNB_WYNIK.md` (zapis wyniku, historia od 5.11), `UNBLIND_OK_zapis_2026-11-05.txt` (treść flagi decyzji o odsłonięciu z sumami sha256; nazwa zmieniona, żeby kopia nie działała jako flaga), `truthnorm/` (normalizacja z prawdy), `r427929/` (run 427929), `pairs/` (428770), `jk/`, `jk_ev200/`, `ubtrig/`, `ubtopo/`, `ubpnet/`, `fls/`, `frag/`, `tfin.md` |
| `GEN_mcstat/` | normalizacja z prawdy: `RAPORT.md`, `out/truth_mbnd.json`, `out/truth_npair.json`, kontrola A3 z 20.11 (`A3_KOREKTA_2026-11-20.md`) |
| `TRIG_A/` | badanie triggera (prerejestracja `SYST_TRIG.md` i wyniki) |
| `NDS_ndstar/` | N_D* dla 428770 (warianty, ponowny fit 12.11) |
| `PD_pairs_diff/` | pary różniczkowe LS/OS, Δφ, Δy; przeskalowanie N_D* |
| `TUNE_pythia/`, `VER_budget/` | budżet systematyk 428770, weryfikacja budżetu (κ_cls, gen_norm) |
| `N2_sf/` | SF klasyfikatora dla D⁰ i D⁺ (`sf_for_pairs_v3.json` i propozycje) |
| `PF_pairs/out/` | diagnostyka weta wspólnych śladów (4.11) |

## Pliki `*.npz.json`

Repo nie przechowuje plików binarnych. Pięć małych tablic potrzebnych do powtórzenia fitu zapisano tekstowo
(`analysis/_odtworzenie/npzjson.py`, format `npz-as-json/1`: nazwa, dtype, kształt i dane każdej tablicy, suma md5 pliku źródłowego):
`U/r427929/pairs/pair_inputs.npz.json`, `U/r427929/pairs/pairs_full.npz.json`, `U/pairs/pair_inputs.npz.json`,
`U/pairs/pairs_full.npz.json`, `U/r427929/nds/dstar_hist_427929.npz.json`. Odtworzenie `.npz` jest dokładne
(sprawdzone `npzjson.same` dla wszystkich pięciu); robi je `analysis/_odtworzenie/zloz_charm_work.py`.

## Czego tu nie ma

- Danych i dużych plików: ntupli ROOT, h5, wyników klasyfikatora (`scores`), nóg par (parquet), map (`.npz`, `.npy`), modeli.
  Pliki `perfile.npz` z jackknife (histogramy per grupa) też nie; do fitu głównego wystarcza ich suma = `pairs_full`.
- Wyników kontroli krotności `MULT_charm` (niezamknięte) i wszystkiego z projektu B (repo `charmbkg`).
- Katalogu podejść odrzuconych `charm-a-nieudane` (cytowany w nocie §8; leży poza `charm_work`).
