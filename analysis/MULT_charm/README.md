# MULT_charm — uzysk D i par D–D vs krotność zdarzenia, dane vs Pythia (16.11.2026)

Projekt A (pary D–D, charmtag) — **nie** projekt B. Nic nie uruchamiane na plikach v2.0 (jeszcze się pobierają); pary danych zablokowane (G5).

| plik | co robi |
|---|---|
| `V20_FORMAT.md` | format ntupli ROOT v2.0 (drzewa, gałęzie, zmiany v1.1→v2.0) + **lista zmian dla `Q_lxplus/gen_norm`** + uwaga o triggerze J12 |
| `PROCEDURA.md` | prerejestracja: obserwable, nominał/warianty, macierz przebiegów, bramki G0–G5, κ (wpływ na μ), reguły decyzji |
| `mult_common.py` | definicje MULT-v1: cięcia (PRE/NOM), selekcja zdarzeń, schematy triggera, estymatory (opp/blx/bl/trk500), biny, w2, fit α |
| `mult_extract.py` | ROOT (v2.0; zgodność v1.1) → npz: zdarzenia (wszystkie), kandydaci po PRE, prawdziwe D (MC) |
| `mult_fit.py` | binowane fity masy (D*: Δm, D⁺, D0 z refleksją z MC) — niezależne od charmtag |
| `mult_analyze.py` | Ŷ(x), R, ε_rel, α/Δα, pary MC, κ; `--pairs-data` → STOP bez `APPROVAL_PAIRS_DATA.txt` (tworzony ręcznie) |
| `mult_gates.py` | G0 (format/liczniki), G1a (PS MBTS w 427929), G1b (plateau sptrk), G1c (pile-up) |
| `mult_truth_gen.py` | przewidywanie Pythii z prawdy generatora `AFP_ddpred/gen` (A3, MPIup/down, Monash): Ŷ(n_ch/⟨n_ch⟩), pary/D |
| `lxplus/` | paczka HTCondor wzorem gen_norm: `mu_common.sh`, `mu_job.sh`, `make_dag_mu.sh`, `submit_mu.sh`, `mu.sub` |
| `tests/` | `make_synth_v20.py` (atrapa v2.0 z lokalnego v1.1 + emulacja triggerów/PS/GRL), `run_local.sh` (zamknięcie, wstrzyknięcie, G5) |
| `out/truth_gen/` | wynik `mult_truth_gen.py` (pełne A3 + warianty) |

## Użycie
Lokalnie: `.venv_h5/bin/python` (uproot 5.7.7, awkward 2.14, scipy). Na lxplus: LCG_110_swan + venv (`q_env`), sprawdzane w `mu_env`.
```
# lxplus (ręcznie/ORCH; nic nie wysłane): ~/private/charm_q/Q_lxplus/mult_charm = lxplus/* + mult_{common,extract,fit,analyze,gates}.py
./submit_mu.sh check                 # G0: pliki dpt3to150, format v2.0 (nSelBL) w 1. pliku każdej próbki, kolejka, AFS
DRYRUN=1 ./submit_mu.sh submit       # mu.dag + komendy
./submit_mu.sh submit                # extract × N plików (+fd_nd, mb_sd/dd/cd jeśli w SAMPLES) → analyze (macierz PROCEDURA §3 + gates)
# lokalnie
python3 mult_analyze.py --data DANE_NPZ --mc MC_NPZ[:w,…] [--eps-mc FD_NPZ] --scheme mbts --est opp --w2 none --out res
python3 mult_truth_gen.py --out out/truth_gen
tests/run_local.sh 16 8 4            # test na lokalnych v1.1 mb_nd (MC) + atrapach v2.0 („dane”)
```
EOS: wejście `/eos/user/k/klyczak/charm_raw_v20/{mb_nd,data_427929,data_428770}` (filtr `PTPROD=dpt3to150`), wyjście `charm_q/mult_charm/{parts,res.tgz,DONE}`.
Szacunek: extract ~1–3 min/plik MC (246k zd.), dane 428770 ~1–2 min/plik (82k zd.), 427929 ~583 pliki; analyze 64 przebiegi × ~2–5 min.

## Testy lokalne (16.11, `tests/run_local.sh 16 8 4`, wyniki `tests/local/res*`)
MC = 16 prawdziwych plików v1.1 mb_nd (3,94 M zd., tryb compat), „dane” = 8 atrap v2.0 z innych plików mb_nd (1,97 M zd.).
- Format: v1.1 (compat: nselbl/leadbl ← n500/lead z Tracks) i atrapa v2.0 czytane tym samym kodem; filtr `dpt3to150` w paczce lxplus działa (atrapa lxplus, `DRYRUN`, węzeł extract → „EOS”).
- **Zamknięcie (dane ≡ MC), D0, mbts:** Δα_fit = −0,13 ± 0,11 (`blx`), +0,05 ± 0,09 (`opp`) — zgodne z 0.
- **Wstrzyknięcie x^0,3:** Δα_fit = +0,17 ± 0,11 (`blx`) = −0,13 + 0,30 → czułość odtworzona.
- **G2 działa i jest potrzebna:** fit D0 przy cięciach NOM w czystym minbias MC ma α_fit − α_prawda = +1,0 (χ² 2162/10) — porównanie
  „dane-fit vs MC-prawda” dałoby fałszywe Z ≈ 8. Stąd wynik główny = Δα_fit (ta sama metoda) i reguła G2 (PROCEDURA §1, §4).
  D⁺ G2 PASS (blx), D* FAIL przy ~250 prawdziwych D* w 16 plikach (statystyka; pełne 80 plików do sprawdzenia).
- D⁺/D* w atrapie danych: brak binów z N/σ > 3 (G3) — atrapa ma sygnał minbias MC z 2 M zdarzeń; w danych 428770/427929 sygnał jest wielokrotnie większy.
- G5: `--pairs-data` → STOP bez `APPROVAL_PAIRS_DATA.txt` ✓.
- `tests/local/synth` (atrapy ROOT) i `tests/local/{mc,data}` (npz) to pliki testowe — można usunąć.
