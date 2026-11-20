# `analysis/` — skrypty, które policzyły wyniki projektu A

Stan na 20.11.2026. Tu leży kod spoza pakietu `charmtag`: skrypty wyniku głównego (run 427929), zadania lxplus/HTCondor,
skrypty raportów i budżetów. Wyniki, które te skrypty dały, są w [`../results/`](../results/README.md).

## Zasada układu

Ścieżki w tym katalogu i w `results/` są takie same jak w katalogu roboczym `~/Documents/charm/charm_work/`:

    charm_work/<ścieżka>   →   analysis/<ścieżka>   (kod: .py .sh .sub .cc .yaml .diff; całe paczki Q_lxplus)
                           →   results/<ścieżka>    (wyniki: .md .json .log .png .pdf)

Dzięki temu każda ścieżka cytowana w nocie, w `UNB_WYNIK.md` i w raportach (`PF2_lxplus/runs/lxplus/unblind/...`,
`GEN_mcstat/...`) wskazuje plik w repo: wystarczy dopisać z przodu `analysis/` albo `results/`.

Skrypty są kopiami stanu z dysku z 20.11.2026, więc ich ścieżki względne zakładają jedno wspólne drzewo
`charm_work`. Żeby je uruchomić, trzeba obie połowy złożyć:

    python analysis/_odtworzenie/zloz_charm_work.py /tmp/cw      # → /tmp/cw/charm_work/ (kopie; repo nietknięte)

W ścieżkach skryptów występuje login CERN autora (`/eos/user/k/...`, `/afs/cern.ch/user/k/...`).

## Odtworzenie wyniku głównego z samego repo

    pip install numpy scipy pandas iminuit
    python analysis/_odtworzenie/odtworz_wynik_glowny.py                   # 427929, kod zamrożony; kilka minut
    python analysis/_odtworzenie/odtworz_wynik_glowny.py --charmtag repo   # ten sam fit bieżącym pakietem charmtag
    python analysis/_odtworzenie/odtworz_wynik_glowny.py --what nom428     # 428770
    python analysis/_odtworzenie/odtworz_wynik_glowny.py --f-one           # kontrola: normalizacja direct (f ≡ 1)

Skrypt składa drzewo, uruchamia `truthnorm_run.py` i porównuje wynik z `results/.../truthnorm/res/tn_*.json`.
Sprawdzone 20.11.2026 (uruchomione tym narzędziem, różnica względem zapisu 0 na wszystkich cyfrach):

| bieg | kod | μ wspólne [68 %] | Z |
|---|---|---|---|
| `nom427` | zamrożony | 0,49894733917589 [0,17454499249837; 0,86467482555756] | 1,60261744409864 |
| `nom427` | pakiet z repo (2.1.1) | 0,49894733917589 [0,17454499249837; 0,86467482555756] | 1,60261744409864 |
| `nom428` | zamrożony | 1,01097126383925 [0,74366213694360; 1,28777065559511] | 4,04521627740340 |

`nom428` to fit na bazie N_D* = 5584; zapis 1,27 powstaje w `truthnorm_report.py` przez przeskalowanie do N_D* = 4436.
Wariantu `--f-one` tym narzędziem nie uruchamiano.

W złożonym drzewie działają też skrypty raportów (sprawdzone 20.11.2026, wynik bajt w bajt równy plikom w `results/`):

| skrypt (w `charm_work/PF2_lxplus/runs/lxplus/unblind/`) | odtwarza |
|---|---|
| `truthnorm/truthnorm_report.py` | `truthnorm/TRUTHNORM.md`, `truthnorm/truthnorm.json` |
| `r427929/budzet_427929.py` | `r427929/budzet_427929.json` i pierwsze 6856 bajtów `BUDZET_427929.md` (sekcja „PROPOZYCJA 19.11” na końcu pliku była dopisana ręcznie) |
| `r427929/mcstat_truth_prop.py` | `r427929/mcstat_truth_prop.json` |
| `r427929/seff_final.py`, `frag/frag.py` | liczby do `r427929/SIGMA_EFF.md` i `frag/FRAG.md` (wypisują na ekran; uruchamiają się, tekstu raportów nie porównywano) |

**Czego z repo odtworzyć się nie da:** wszystkiego przed ostatnim fitem. Ntuple ROOT, h5, wyniki klasyfikatora, nogi par,
mapy wag i próbki prawdy (parquet) są na EOS/dysku, nie w repo. Skrypty tych kroków są tu jako zapis tego, co wykonano.

## Która wersja `charmtag` policzyła co

| krok | kod | gdzie w repo |
|---|---|---|
| nogi, histogramy par, szablony, N_D* biegu par (1901), normalizacja direct — oba runy (lxplus) | snapshot `Q_lxplus/pf_pairs/src` = charmtag 2.0.3 z 10.11.2026 (deklaruje „2.0.2”) | kodu nie ma w repo |
| fit końcowy wyniku głównego i 428770 z prawdy mb_nd (lokalnie, 19.11) | `charm_work/PF_pairs/src` (deklaruje „0.1.0”, stan z 4.11) | `PF_pairs/src/` (kopia części A) |
| N_D* = 2046 i 14 wariantów (NDS) | `dstar_refit_nds.py` = `charmtag/dstar_refit.py` z charmtag 2.1.0 (stan z 16.11; md5 zgodne) | `Q_lxplus/ub_427929/`, `Q_lxplus/ub427_nds/` |
| odsłonięcie 428770 z 5.11 (Z = 4,06σ) | charmtag 2.0.2 z poprawką ID śladu z 5.11 | kodu nie ma w repo; sumy sha256 kodu par: `results/PF2_lxplus/runs/lxplus/unblind/UNBLIND_OK_zapis_2026-11-05.txt` |
| bieżący pakiet | `charmtag/` 2.1.1 (2.1.0 ze stanu z 16.11 + poprawki z 20.11) | katalog główny; fit końcowy daje te same liczby co kod zamrożony |

## Łańcuch wyniku głównego (run 427929) — skrypt po skrypcie

Ścieżki względem `analysis/`; wyniki względem `results/`. Skrót `U/` = `PF2_lxplus/runs/lxplus/unblind/`.

| # | co liczy | skrypt | wejścia | wynik w repo |
|---|---|---|---|---|
| 1 | filtr h5 v2.0 → selekcja nominału, klasa MB (GRL ∧ MBTS) | `Q_lxplus/ub_427929/u4_filter.py` (zadanie `u4_job.sh`) | h5 v2.0 run 427929 (EOS) | statystyki filtra w `U/r427929/ub427.json` → `filter` |
| 2 | wagi zdarzeń c_w dla 427929 (metoda stosunku R_X) | `Q_lxplus/ub_427929/u4_cw.py`, `gn_root.py`, `gn_combine.py`, `effmap_gn.py` | ROOT v2.0 mb_nd / 428770 / 427929 (EOS), `Q_lxplus/inputs/mb_rates_2d.json` | `U/r427929/cw/cw.json` |
| 3 | nogi, kategorie score, histogramy par, normalizacja direct, fit kanałów | `charmtag pairs --unblind`, `charmtag fit2d --unblind` (snapshot charmtag 2.0.3 z 10.11); config `U/r427929/pairs/cfg427.yaml` | h5 po filtrze, modele XGB H, score, `mb_rates_427929.json` z kroku 2 (poza repo) | `U/r427929/pairs/` (`pair_inputs.json`, `pair_inputs.npz.json`, `pairs_full.npz.json`, `fit2d_unblind.json`, `charmtag.log`) |
| 4 | N_D* metodą NDS (2046 ± 485) | `Q_lxplus/ub_427929/u4_ndstar.py` + `dstar_refit_nds.py` | h5 po filtrze | `U/r427929/ndstar_427929.json` |
| 5 | raport UB427 (μ direct, porównanie z 428770, przeskalowanie do N_D* = 2046) | `Q_lxplus/ub_427929/u4_report.py` | wyniki kroków 2–4, odsłonięcie 428770 | `U/r427929/UB427.md`, `ub427.json` |
| 6 | 14 wariantów N_D* → systematyka N_D* 41 % | `Q_lxplus/ub427_nds/n8_ndstar.py`, `nds_var.py`, `dstar_refit_nds.py` | h5 po filtrze; `ndstar_full_428770.json` | `U/r427929/nds/` (`NDS_427929.md`, `ndstar_variants_427929.json`, `dstar_hist_427929.npz.json`) |
| 7 | szerokość piku D* w klasach triggera (kryterium `ubndsw`) | `Q_lxplus/ub_ndsw/ndsw.py` | h5 obu runów, tabela triggerów | `U/r427929/ndsw/` (`NDSW.md`, `ndsw.json`) |
| 8 | jackknife po zdarzeniach (200 grup) → czynnik 1,13; gładkie minimum wspólne 0,567 | `Q_lxplus/ub_jk/ubjk_run.py` (zadania `Q_lxplus/ub427_jk/`) | krok 3 + przebudowa nóg | `U/r427929/jk_ev200/` (`ubjk.md`, `ubjk.json`, `perfile.json`) |
| 9 | czynniki f_ch = [pary/D*] z prawdy mb_nd / direct | `GEN_mcstat/truth_mbnd.py` (wejścia z `Q_lxplus/gen_norm/`: `gn_root.py`, `gn_legs.py`, `gn_combine.py`) | `GEN_mcstat/runs/lxplus_v20/*.parquet`, mapy `*.npz` (poza repo) | `GEN_mcstat/out/truth_mbnd.json` |
| 10 | **fit wyniku głównego**: n_pair × f_ch, pełny ponowny fit, profil wspólny | `U/truthnorm/truthnorm_run.py --what nom427` (minimalizator: `Q_lxplus/ub_jk/ubjk_run.py:replica_fit`; `charmtag` z `PF_pairs/src`) | kroki 3 i 9 | `U/truthnorm/res/tn_nom427.json` (kontrola f ≡ 1: `tn_nom427_f1.json`) |
| 11 | budżet systematyk 427929 | `U/r427929/budzet_427929.py` (ok. 30 liczb wpisanych w skrypcie; wiersze kształtów, SF, κ przeniesione z 428770) | `ub427.json`, `jk_ev200/ubjk.json`, `cw/cw.json`, `ndstar_427929.json`, `nds/ndstar_variants_427929.json`, `PF2_lxplus/runs/lxplus/pairs/syst/syst_blind.json` | `U/r427929/BUDZET_427929.md`, `budzet_427929.json` |
| 12 | zapis wyniku: stat × 1,13, syst i N_D* względnie z budżetu | `U/truthnorm/truthnorm_report.py` (stałe `JK427 = 1.13`, `JK428 = 1.18` wpisane w skrypcie) | `U/truthnorm/res/*.json`, `budzet_427929.json` | `U/truthnorm/TRUTHNORM.md`, `truthnorm.json` |
| 13 | granica σ_eff | `U/r427929/seff_final.py` | liczby z kroku 12 (wpisane) | `U/r427929/SIGMA_EFF.md` |
| 14 | fragmentacja (górna strona wiersza tune) | `U/frag/frag.py` | liczby wpisane (ALICE, Pythia) | `U/frag/FRAG.md` |

## Pozostałe katalogi

| katalog | co zawiera | wyniki |
|---|---|---|
| `_odtworzenie/` | `zloz_charm_work.py`, `odtworz_wynik_glowny.py`, `npzjson.py` (dopisane 20.11) | — |
| `PF_pairs/` | zamrożony `src/` (część A), diagnostyka weta wspólnych śladów `pf_chan.py`, `pf_diag.py`, `pf_impact.py`, łata `patch/pf_trackid.diff` | `PF_pairs/out/` |
| `PF2_lxplus/runs/lxplus/unblind/` | skrypty leżące przy wynikach: `truthnorm/*.py`, `r427929/*.py`, `frag/frag.py`, `r427929/pairs/cfg427.yaml` | ten sam katalog w `results/` |
| `Q_lxplus/pf_pairs/` | zadania par na lxplus (ślepy bieg PF2 i odsłonięcie 5.11): `pf_job.sh`, `unblind_job.sh`, `tools/` (raporty, bramki), `cfg/sf_v3.yaml` | `PF2_lxplus/` |
| `Q_lxplus/ub_427929/`, `ub427_nds/`, `ub427_jk/`, `ub_ndsw/` | run 427929 (kroki 1–8 wyżej) | `U/r427929/` |
| `Q_lxplus/ub_jk/` | jackknife μ po plikach i grupach zdarzeń, `replica_fit` (gładkie minimum wspólne) | `U/jk/`, `U/jk_ev200/` |
| `Q_lxplus/ub_trig/` | μ per klasa triggera w 428770 (`ubtrig_run.py`), test zerowy (`ubtrig_null_agg.py`), tabela triggerów z ROOT (`trg_root_extract.py`) | `U/ubtrig/` |
| `Q_lxplus/ub_comb/` | wspólny fit próbek bez triggera dżetowego — analiza wykonalności, bieg nie wysłany | opis w nocie §7 (C) |
| `Q_lxplus/ub_pnet/`, `pair_planes/`, `prez_net/` | sieć par jako kontrola (`upn_fit.py`), płaszczyzny mas par, porównanie drzewo/sieć do prezentacji | `U/ubpnet/`, `U/r427929/pnet/`, `PREZ_2026-11-19/p3/` |
| `Q_lxplus/ub_topo/` | topologia \|Δη\| par vs hipotezy generatora | `U/ubtopo/` |
| `Q_lxplus/gen_norm/` | normalizacja z prawdy: mapy wag w2, mapy wydajności, nogi MC (`gn_*.py`), testy lokalne | `GEN_mcstat/runs/lxplus_v20/` |
| `Q_lxplus/h5_v20/` | zadania konwersji ROOT v2.0 → h5 (bez konwertera, zob. niżej) | — |
| `Q_lxplus/ev_dstar/` | ekstrakcja par D⁰←D* z kluczem zdarzenia (`ev_extract.py`) i test klucza | — |
| `Q_lxplus/common.sh`, `config/`, `inputs/` | wspólne ścieżki i funkcje zadań, config lxplus, zamrożone `mb_rates_2d.json` i `sf_for_pairs_v2.json` | — |
| `GEN_mcstat/` | `truth_mbnd.py`, `truth_npair.py`, `a3_korekta.py`, `gen_inventory.py` | `GEN_mcstat/` |
| `TRIG_A/` | trigger: emulacja L1_J12 na generatorze (`trigmc/`, `pythia_dd_jet.cc`), wydajność J12 z danych (`j12eff/`), SF w kategoriach par (`pairsf/`), sieć zamiast drzewa (`nnh/`), prerejestrowana miara `syst_trig.py` | `TRIG_A/` |
| `NDS_ndstar/` | N_D* dla 428770: warianty, ponowny fit z 12.11, test `fit_single` | `NDS_ndstar/` |
| `PD_pairs_diff/` | pary różniczkowe LS/OS: przeskalowanie N_D* (`ndstar/`), łaty kodu `*.diff` (już w pakiecie) | `PD_pairs_diff/` |
| `TUNE_pythia/` | systematyka tune Pythii i budżet 428770 | `TUNE_pythia/` |
| `VER_budget/` | weryfikacja budżetu: κ_cls (`kcls_check.py`) | `VER_budget/` |
| `N2_sf/` | SF D⁰ z tagu D* i D⁺ (`n2_*.py`) | `N2_sf/` |
| `MULT_charm/` | kontrola krotności D⁰ (procedura, skrypty, zadania); wyniki niezamknięte | — (nie włączone) |
| `PREZ_2026-11-19/` | `deck/build_deck.py`, wykresy `p1/make_p1.py`, `p6/make_p6.py` | `PREZ_2026-11-19/` |
| `ORCH/stages/` | etapy pętli ORCH dla projektu A (wysyłka i pobieranie zadań lxplus) | — |

## Czego tu nie ma

- `Q_lxplus/h5_v20/charmpp_h5_v20.py` — konwerter ROOT → h5 napisany jako odtworzenie 1:1 funkcji z kodu analizy ATLAS na GitLabie CERN;
  pominięty z ostrożności (zasada: bez kopii kodu z GitLaba CERN). Jest lokalnie w `charm_work/Q_lxplus/h5_v20/`.
- `ORCH/pipeline.yaml` — wspólny dla projektów A i B.
- Snapshoty `src/` w paczkach `Q_lxplus/*` (kopie `charmtag` z modułami projektu B i modelami ONNX); `gen_norm/src/SNAPSHOT.txt` mówi, z którego stanu pakietu pochodzą.
- Pliki-flagi decyzji autora (`UNBLIND_OK`): repo ich nie zawiera pod nazwą, której szukają skrypty. Treść flagi odsłonięcia z 5.11
  (data decyzji i sumy sha256 configu oraz kodu par) jest zachowana jako zapis w `results/PF2_lxplus/runs/lxplus/unblind/UNBLIND_OK_zapis_2026-11-05.txt`.
- Kopie zapasowe `*.bak_*`, `*.przed_*`, pliki `diag.txt` (zrzuty `condor_q`), archiwa `.tgz`, `__pycache__`.
- Paczki `Q_lxplus` projektu B i paczki, których przynależności nie rozstrzygnięto (`tagtrig`, `trkjet*`, `pnet_v4`, `pz_data`, `pz_improve`).
