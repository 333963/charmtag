# Format danych wejściowych — h5 v1.1 (i v1.2)

Pliki HDF5 produkowane z mini-ntupli ROOT (CharmAnalysis/charmpp, produkcja v1.1 „Aug7”) skryptem konwersji
(`charmpp_h5.py`, etap H; rozszerzenie v1.2: `charmpp_h5_v12.py`, etap H2). Jeden plik = jedna próbka
(dane `data_<run>`, MC inkluzywne `mb_<proces>`, MC forced-decay `fd_<proces>`), ewentualnie w częściach („parts”).
Lokalizacje plików i zasady dostępu: `DATA_ACCESS.md` (nie w repozytorium).

## Struktura

Jeden dataset `events` typu compound (**jeden wiersz = jeden kandydat D-mezonu**, nie zdarzenie), 60 kolumn w v1.1.
Wiersze tego samego zdarzenia są ciągłe (blok), kolejność jak w ROOT. Kandydaci różnych kanałów rozpadu są
w tym samym datasecie i rozróżniani polem `Dmeson_decayType`:

| `decayType` | kanał | córki | uwagi |
|---|---|---|---|
| 10 | D0 → K⁻π⁺ | 2 | **obie hipotezy K/π** tej samej pary torów są zapisane jako osobne wiersze (ok. 10 % wierszy w danych, 24 % w fd) |
| 1 | D*⁺ → D0 π⁺ₛ | 3 | `Dmeson_delta_mass` = m(D0πₛ) − m(D0); D0 z D* **nie** jest zapisywany dodatkowo jako wiersz `decayType 10` (v1.1; w starszej produkcji v9 był — 10,8 % duplikatów) |
| 4 | D⁺ → K⁻π⁺π⁺ | 3 | kaon = córka 0 (znak przeciwny do pionów) |
| 11 | Λc⁺ → pK⁻π⁺ | 3 | nieużywany w tej analizie |

Pakiet mapuje kolumny na nazwy logiczne przez `columns:` w `charmtag/data/default.yaml` (lista kandydatów na nazwę,
pierwsza istniejąca wygrywa), więc inne wersje ntupli wymagają tylko zmiany configu.

## Kolumny (v1.1) i jednostki

**Uwaga na jednostki:** większość w GeV, ale pT córek i m(Kπ) w **MeV**. Konwersję do GeV robi `features.compute`
(mnożniki w `units:` configu: `d_pt: 0.001`, `m_kpi: 0.001`); dalej w pakiecie wszystko jest w GeV
(`charmtag/constants.py`). Moduł `pairs` liczy histogramy w MeV (krawędzie binów `pairs.edges`), konwersja na brzegu.

| grupa | kolumny | jednostka / typ | znaczenie |
|---|---|---|---|
| zdarzenie | `run_number`, `event_number`, `event_weight`, `actual_mu`, `average_mu` | — | w **fd** `event_number` nie jest unikalny (numer w pliku; 10 k unikalnych na 187 k bloków) — k-fold grupuje nadmiarowo, OOF pozostaje poprawny |
| kandydat | `Dmeson_decayType`, `Dmeson_species`, `Dmeson_charge` | int | `charge` = ładunek kandydata (D0: 0) |
| masa | `Dmeson_mass` | GeV | masa pod hipotezą wiersza |
| | `Dmeson_delta_mass` | GeV | Δm dla D*, **−0,999** gdy nie-D* |
| | `Dmeson_mKpi1`, `Dmeson_mKpi2` | **MeV** | m(Kπ) dla D⁺ (weto D*: min < 160 MeV); dla D0 = −0,999; liczone z innych pędów niż `daughter*_pt` (różnice ~2 MeV) |
| kinematyka | `Dmeson_pt`, `Dmeson_eta`, `Dmeson_phi` | GeV, —, rad | |
| wierzchołek | `Dmeson_Lxy` | mm | preselekcja produkcji: Lxy ≥ 0,25 mm (wszystkie próbki v1.1) |
| | `Dmeson_chi2`, `Dmeson_cosThetaP`, `Dmeson_cosThetaStar`, `Dmeson_impact`, `Dmeson_impactSignificance`, `Dmeson_z0SinTheta`, `Dmeson_maxdR`, `Dmeson_min_daughter_d0` | — | cosθ* = kąt K w układzie spoczynkowym D (zgodny z przeliczeniem z córek, corr 0,99996) |
| córki `i = 0,1,2` | `Dmeson_daughter{i}_pt` | **MeV** | córka 0 = K (hipoteza), dalej piony; dla D0 trzecia córka = 0 |
| | `Dmeson_daughter{i}_eta`, `_phi` | —, rad | **po refit wierzchołka kandydata** — ten sam ślad w dwóch kandydatach ma inne wartości (~0,1 %) |
| | `Dmeson_daughter{i}_d0`, `_z0SinTheta`, `_z0SinThetaPV` | mm | parametry toru względem PV, **niezależne od refitu** → identyfikator śladu w parach (`pairs.track_id: d0z0`) |
| | `Dmeson_daughter{i}_dEdx`, `_npts_dEdx` | jedn. MIP (~0,87 dla π przy p ≫ 1 GeV), int | dE/dx pikseli; `npts = 0` ⇔ brak pomiaru |
| D* | `Dmeson_soft_pion_d0`, `Dmeson_soft_pion_dR` | mm, — | |
| izolacja | `Dmeson_ptcone40/60/80/100` | GeV | **zakazane jako cechy** (kontrola tylko) |
| prawda (MC) | `Dmeson_truth_category` | int | 0 = c, 1 = b (sygnał), ≥ 8 tło; brak = NaN/−1 w danych |
| | `reco_channel`, `truth_channel`, `presel_channel` | tekst | |
| generator | `pdfInfo_pdgid1/2` | int | tylko MC |

Czego **nie ma** w v1.1: ładunków i `pdgId` córek, `trackId`, `truthPdgId` córek, pędów `daughter*_p`
(liczymy p = pT·coshη), liczby śladów / pT wiodącego śladu w zdarzeniu (wagi krotności z MC przenosimy jako czynnik).

## Konsekwencje dla analizy

- **Hipotezy K/π (D0).** Dopasowanie do prawdy w v1.1 ignoruje hipotezę: dla par z oboma wierszami oba są
  `truth_category ∈ {0,1}` w 65 % (nigdy dokładnie jeden). Właściwą hipotezę wskazuje heurystyka `closest_mass`
  (`samples.mark_swaps`: dalsza od m0 → swap, `truth += 100`); jest definicyjna, nie fizyczna — w v1.2 zastępuje ją
  `Dmeson_hypo_swap` (`training.truth_swap_source: pdgid`).
- **Wspólne ślady w parach.** Identyfikator śladu z (pT, η) córki **nie działa** po refit; używamy bitów
  `(d0, z0SinThetaPV)` — zob. `docs/double_charm.md` (sekcja PF).
- **Duplikaty.** v1.1: 0 duplikatów na 5,4 M sprawdzonych wierszy. Strażnik `io.dedup: false|report|true`
  (`charmtag/dedup.py`, klucz: run, event, decayType, bity pT/η/φ córek, bity masy) — włączyć `report` przy każdej
  nowej produkcji (v9 miał 10,8 % duplikatów D0 z finderu D*, bez cięcia Lxy).
- **Dwie hipotezy ≠ duplikat** (inna masa) — zbija je arbitraż par (`pairs.build.arbitrate`), nie dedup.

## v1.2 (etap H2, konwersja `charmpp_h5_v12.py`)

Te same 60 kolumn bit w bit + 19 nowych (int32, sentinel −999):
`Dmeson_daughter{0,1,2}_{charge,pdgId,trackId,truthPdgId,truthParentPdgId}`, `Dmeson_fit_charge` (D⁺ ±1; D0/D* 0),
`Dmeson_truthPdgId`, `Dmeson_n_daughter_truth_ok`, **`Dmeson_hypo_swap`** (0 zgodna / 1 zamieniona / −1 nieokreślona).
Ładunki potwierdzają: kandydatów wrong-sign w ntuplach nie ma (odrzucane w finderach przed zapisem).
Loader toleruje brak nowych kolumn (v1.1 → stare zachowanie); przy v1.2 klucz dedup można przełączyć na `trackId`.

## Minimalny przykład odczytu

```python
import h5py, numpy as np
with h5py.File("mb_nd_part0.h5") as f:
    ev = f["events"]
    dt = ev["Dmeson_decayType"][:]
    m = ev["Dmeson_mass"][:]                       # GeV
    d_pt0 = ev["Dmeson_daughter0_pt"][:] * 1e-3    # MeV → GeV
    d0 = (dt == 10)
```

W pakiecie: `charmtag.io.open_events` (jeden plik lub lista części, `row` = wiersz globalny), `charmtag.samples.load_species`
(filtr `decayType`, jednostki, cechy, swapy, dedup).
