# charmtag — double charm: pary D–D w danych ATLAS minimum bias (pp, 13,6 TeV, 2022)

Projekt A. Pomiar par mezonów powabnych (D⁰D⁰, D⁺D⁰, D⁺D⁺) w jednym zderzeniu pp: klasyfikator pojedynczego D (XGBoost,
k-fold OOF, tło z sidebandów danych), fit masy, pary D–D z wetem wspólnych śladów i fitem 2D, normalizacja do N_D*, systematyki.
Mierzona wielkość: μ = [N_par / N_D*]_dane / [N_par / N_D*]_Pythia 8 A3.

Stan opisu: **20.11.2026**. Repo zawiera kod (`charmtag/`), skrypty analizy (`analysis/`), wyniki z liczbami (`results/`)
i dokumentację (`docs/`).

Projekt B (tło podobne do sygnału) to osobne repo [charmbkg](https://github.com/333963/charmbkg), pakiet `charmbkg`.

## Wynik główny — run 427929 (czysty minimum bias)

Zatwierdzony przez autora 19.11.2026. Próbka: run 427929, trigger HLT_mb_mbts_L1MBTS_1, normalizacja par z prawdy mb_nd.

> **μ = 0,50 +0,41/−0,37 (stat, × 1,13 jackknife) +0,28/−0,19 (syst) ± 0,21 (N_D*)**
>
> Z = 1,60 względem braku par (≈ 1,42 z korelacjami zdarzeń). Względem Pythii A3 (μ = 1): −0,93σ.
> Interpretacja: σ_eff > ~6 mb (95 % CL, przy f_DPS = 0,40) — granica słaba.

W czystym minimum bias nie ma nadwyżki par względem Pythii; sam sygnał par nie jest istotny (Z = 1,6).

- Nota: [`results/NOTA_WYNIK/nota.md`](results/NOTA_WYNIK/nota.md).
- Każda liczba → plik i skrypt: [`results/README.md`](results/README.md).
- Trzy znane niespójności tego wyniku: [Znane problemy](#znane-problemy).

## Run 428770 — obserwacja obciążona triggerem dżetowym (pomiar pomocniczy)

Strumień physics_Main; 77 % N_D* pochodzi ze zdarzeń HLT_noalg_L1J12. Pary odsłonięto 5.11.2026.

> μ = 1,27 +0,41/−0,40 (stat, × 1,18) ± 0,36 (syst, prawda mb_nd) ± 0,28 (N_D*); Z = 4,05 (≈ 3,43 z korelacjami zdarzeń).
> Zapis z normalizacją direct (16.11): 1,50 +0,48/−0,46 ± 0,61 ± 0,33. Odsłonięcie 5.11: Z = 4,06σ.

Nadwyżka siedzi w zdarzeniach wybranych wyłącznie triggerem L1_J12 (μ = 1,37 [1,05; 1,71] wobec −0,40 ± 0,37 w pozostałych),
pary są bardziej kolinearne niż w Pythii (λ̂ = −1,11). Wyniku nie interpretuje się jako DPS.
Źródła: `results/PF2_lxplus/runs/lxplus/unblind/` (`UNB_WYNIK.md`, `truthnorm/TRUTHNORM.md`, `ubtrig/UBTRIG.md`, `ubtopo/UBTOPO.md`), `results/TRIG_A/`.

## Znane problemy

Trzy niespójności wyniku głównego. **Nie zostały poprawione i żadna liczba wyniku nie została zmieniona** — to decyzja autora
(20.11.2026); tutaj są tylko opisane.

1. **Dwie bazy N_D\*.** Fit z normalizacją z prawdy dla 427929 wczytał wejścia par z N_D* = 1901 ± 445
   (`results/PF2_lxplus/runs/lxplus/unblind/r427929/pairs/pair_inputs.json` → `meta.dstar_fit_data`), a nota i budżet podają
   N_D* = 2046 ± 485 (`r427929/ndstar_427929.json`, metoda NDS). Dla 428770 raport przeskalował wynik do bazy NDS (5584 → 4436),
   dla 427929 nie (`analysis/PF2_lxplus/runs/lxplus/unblind/truthnorm/truthnorm_report.py`). μ ∝ 1/N_D*, więc na bazie 2046
   wynik to μ ≈ 0,46 zamiast 0,50. Te dwie liczby N_D* pochodzą z różnych wersji fitu Δm.
2. **Współczynniki dane/MC w nominale.** Wejścia par mają `sf_json: null` (SF klasyfikatora = 1) i `tag_sf: 0.7`
   (`pair_inputs.json` → `meta.config`), choć zmierzony SF dla D⁰ to 0,88–0,95 (`results/N2_sf/`). Według przeglądu kodu
   zastosowanie zmierzonych SF podniosłoby μ o ok. 5–13 % (oszacowanie z przeglądu, niepowtórzone tutaj).
3. **Szerokość piku D\*.** Fit N_D* ma wolną skalę szerokości i daje 1,42 × MC w biegu par (1,40 w NDS); w 428770, tym samym
   detektorem, 1,02. Ze skalą 1 N_D* = 1417 ± 200, co odpowiada μ ≈ 0,72 (`results/.../r427929/nds/NDS_427929.md`,
   `r427929/ndsw/NDSW.md`; kryterium rozstrzygające zapisane przed wynikiem dało „nierozstrzygnięte”).

Żadna z tych rzeczy nie zmienia wniosku fizycznego (błąd całkowity ok. ±0,5), ale przesuwają wartość centralną w zakresie ok. 0,46–0,72.

Inne ograniczenia, o których trzeba wiedzieć:

- Normalizacja MC opiera się na 22 / 41 / 15 parach w MC (`pair_inputs.json` → `meta.mc_rates_unweighted.n_pairs`).
- Z samego repo odtwarza się **ostatni krok** (fit i raporty), nie cały łańcuch: dane, h5, wyniki klasyfikatora i nogi par są poza repo.
- Fit końcowy wykonano 19.11 kodem spoza pakietu (`charm_work/PF_pairs/src`, „0.1.0”), a wejścia zbudował na lxplus snapshot
  pakietu w stanie charmtag 2.0.3 z 10.11.2026. Oba są udokumentowane, a fit końcowy bieżącym pakietem daje identyczne liczby — zob. niżej.
- Część wierszy budżetu systematyk 427929 jest przeniesiona z 428770, a ok. 30 liczb jest wpisanych w skrypcie budżetu.
- CI nie instaluje `charmbkg` i nie uruchamia testów projektu B (lokalnie przechodzą).

## Mapa repo

| ścieżka | co to jest |
|---|---|
| `charmtag/` | pakiet Pythona (2.1.1), CLI `charmtag …`: rdzeń (io, features, samples, model, pipeline, massfit, dstar_refit, report), `pairs/` (budowa par, model i fit 2D, pary różniczkowe), `syst/`, `hand/`; `_bkg_alias.py` — aliasy modułów projektu B → `charmbkg` |
| `tests/` | pytest (marker `slow`); `conftest.py` pomija pliki testów projektu B, gdy brak `charmbkg` |
| `analysis/` | skrypty, które policzyły wyniki: zadania lxplus, raporty, budżety; zamrożony kod fitu końcowego; narzędzia odtworzenia. Opis: [`analysis/README.md`](analysis/README.md) |
| `results/` | raporty, tabele, JSON z liczbami, nota, prezentacja. Indeks liczb: [`results/README.md`](results/README.md) |
| `docs/` | opis metody ([`double_charm.md`](docs/double_charm.md)), format danych, lxplus |
| `configs/examples/` | przykładowe configi |
| `scripts/studies/`, `inputs/` | dwa starsze skrypty pomocnicze i jeden plik wejściowy z listopada (sprzed `analysis/`) |

`analysis/` i `results/` mają ten sam układ ścieżek co katalog roboczy `charm_work`, więc ścieżki cytowane w nocie i raportach
(`PF2_lxplus/runs/lxplus/unblind/...`) wskazują pliki w repo po dopisaniu `results/` albo `analysis/`.

## Jak odtworzyć wynik główny

```bash
python -m venv .venv && source .venv/bin/activate
pip install numpy scipy pandas iminuit
python analysis/_odtworzenie/odtworz_wynik_glowny.py                   # kod zamrożony; kilka minut
python analysis/_odtworzenie/odtworz_wynik_glowny.py --charmtag repo   # ten sam fit bieżącym pakietem charmtag
```

Skrypt składa w katalogu tymczasowym drzewo `charm_work` z `analysis/` i `results/`, odtwarza pliki `.npz` z ich zapisu tekstowego,
uruchamia `truthnorm_run.py` i porównuje μ, przedział 68 % i Z z zapisem w
`results/PF2_lxplus/runs/lxplus/unblind/truthnorm/res/tn_nom427.json`. Sprawdzone 20.11.2026: różnica 0 na wszystkich cyfrach,
zarówno kodem zamrożonym, jak i pakietem z repo; tak samo fit 428770 (`--what nom428`). Raporty (`TRUTHNORM.md`, `budzet_427929.json`) odtwarzają się w złożonym drzewie
bajt w bajt — zob. [`analysis/README.md`](analysis/README.md).

Która wersja kodu policzyła co:

| krok | kod |
|---|---|
| nogi, histogramy par, N_D* biegu par, normalizacja direct (lxplus) | snapshot charmtag 2.0.3 z 10.11.2026 |
| fit końcowy z normalizacją z prawdy (19.11, lokalnie) | `analysis/PF_pairs/src/` — zamrożona kopia |
| odsłonięcie 428770 (5.11) | charmtag 2.0.2 z poprawką ID śladu z 5.11 |
| bieżący pakiet | `charmtag/` 2.1.1 — fit końcowy daje te same liczby co kod zamrożony |

## Instalacja i testy

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,onnx]"          # instaluje też charmbkg (zależność charmbkg@v1.3.0)
pytest -q -m "not slow"
```

Bez dostępu do `charmbkg` (tak działa CI):

```bash
pip install numpy pandas h5py scipy scikit-learn pyyaml pytest "iminuit>=2.20" matplotlib xgboost onnxruntime
pip install --no-deps -e .
pytest -q -m "not slow"               # pliki testów projektu B są pomijane, z komunikatem w nagłówku i podsumowaniu
```

Stan testów 20.11.2026 (`-m "not slow"`, Python 3.14): z `charmbkg` 1.3.0 — 179 passed, 11 skipped (6 z nich: brak `pyarrow` w środowisku); bez `charmbkg` — 96 passed, 7 skipped.

Na prawdziwych ntuplach: skopiuj `configs/examples/mac_v11.yaml`, ustaw `paths.h5_dir`, `paths.out_dir` i nazwy plików próbek, potem
`charmtag all -c twoj.yaml --species D0` (`--species Dplus`), `charmtag pairs -c twoj.yaml`, `charmtag fit2d -c twoj.yaml --ntoys 50`.
Uruchamianie na lxplus/HTCondor: [docs/running_on_lxplus.md](docs/running_on_lxplus.md).

## Dane wejściowe

h5 „v1.1”/„v1.2”/„v2.0” (compound dataset `events`, kolumny `Dmeson_*`; D⁰ = decayType 10, D⁺ = 4, D* = 1). Mapowanie kolumn
i jednostki (masa w GeV, pT córek w MeV) są w `charmtag/data/default.yaml`; opis formatu: [docs/data_format.md](docs/data_format.md).
Surowe ntuple, h5, modele, wyniki klasyfikatora i mapy wag nie są w repo. Małe zamrożone wejścia zadań lxplus:
`analysis/Q_lxplus/inputs/`.

## Reguły

- Bez danych i plików binarnych w repo (ntuple, `.root`, `.h5`, `.npz`, modele). Małe tablice wejściowe fitu są zapisane tekstowo (`*.npz.json`).
- `charmtag pairs` / `fit2d` bez `--unblind` nie zapisują i nie rysują centrum płaszczyzny mas par.
- Kod analizy ATLAS na GitLabie CERN (`charmpy`, `CharmAnalysis`) jest tylko czytany; to repo nie zawiera jego kopii.
- Projekt B osobno: nic z repo `charmbkg` nie trafia tutaj poza zależnością od pakietu `charmbkg`.

## Historia

Kolejne zapisy wyniku (5.11: μ 1,25 → 1,57; 9.11: 1,49; 13.11: 1,50; 19.11: wynik główny 427929) i ich uzasadnienia:
`results/PF2_lxplus/runs/lxplus/unblind/UNB_WYNIK.md`.
