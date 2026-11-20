# gen_norm — wagi w2 i mapy ε trybu `direct` z mb_nd (lxplus/HTCondor), 16.11.2026

Kontekst: `charm_work/GEN_mcstat/RAPORT.md` §5 (a–c), `VER_budget/RAPORT.md` (A). Procedura użycia wyników (prerejestracja):
`charm_work/GEN_mcstat/PROCEDURA.md`. **Pakiet liczy tylko MC i rozkład (n_trk, lead_pt) zdarzeń danych — żadnych par danych, żadnego μ.**

## v2.0 (17.11.2026) — domyślnie `GN_FMT=v20` (kopie v1.1: `*.bak_v11`)
- Wejście: `RAW_GN=/eos/user/k/klyczak/charm_raw_v20/{mb_nd,data_428770}`, tylko `*${PTPROD}*` (`dpt3to150`; `find`, nie `ls`), oczekiwane 80 / 97 plików
  (check = STOP przy innej liczbie). `../common.sh` (RAW v1.1) bez zmian — nadpisanie w `gn_common.sh`. Wyjście `GN_OUT=charm_q/gen_norm_v20`,
  ORCH `stages/gen_norm.sh` → `GEN_mcstat/runs/lxplus_v20`. `GN_FMT=v11` = dawne zachowanie (ścieżka v1.1 na EOS 17.11 niedostępna).
- `gn_root.py`: D0 w `sel` tylko 1750 < m < 2000 MeV (v2.0 ma 1,70–2,05; na v1.1 bez wpływu); `ev_*` + `nselbl, leadbl [GeV], npu, npv, grl, run, lb`
  i w danych `tp_/ps_<łańcuch>` (11 łańcuchów = `MULT_charm/mult_common.TRIG_STORE`) — gotowe pod w2(nSelBL, leadBL) i w2 per trigger.
- Nogi (`legs`) nadal z h5 v1.1 + scores nominalne (brak h5 v2.0). Join nóg ↔ ROOT v2.0 po (tag joba, eventNumber): podział mb_nd v2.0 = v1.1
  (plik `_000002`: te same eventNumber/n_trk/Tracks), sprawdzane w B1a (`legs_without_event = 0`).
- `gn_combine.py --fmt v20 --ref-anchors h5,h5`: referencja w2 = bloki zdarzeń v2.0 wokół kotwic (eventNumber z h5 v1.1 `_000021/_000030`),
  `ref_anchors.json`; akceptacja = B1a/B1b (blokujące) + B1c (D0/D*/c_w informacyjnie) — `GEN_mcstat/PROCEDURA.md` §3a.
- Weryfikacja 17.11 (odczyt lxplus, 1 plik/próbkę): wszystkie gałęzie używane przez `gn_root`/`effmap_gn` są (z `_NOSYS`), `TruthDmesons_*` w MC,
  11/11 `trigPassed`/`trigPS` w danych, `metadata` + `DSelection_Vertices` obecne; `passGRL` tylko w danych.
- Testy: `tests/run_local.sh` (tryb v11, bit w bit z etap1 — OK), `tests/test_v20.py` (kotwice, kolumny v2.0 na atrapie, okno D0 — OK).

## Co liczy
| węzeł | wejście | wyjście (EOS `charm_q/gen_norm/`) |
|---|---|---|
| `mc_<job>` × 80 | ROOT v1.1 mb_nd (`$RAW/<mb_nd>/…root`, lista z `config/samples.txt`) | `parts/mc/mb_nd_<job>.npz`: zdarzenia (evnum, n_trk, lead_pt, n_true_dstar, n_corr_d0/dp = etap1 `mb_truth.py`), poprawne kandydaty etap1 (trackId), prawdziwe D (`effmap_gn.py`) |
| `da_<job>` × 39 (`DATA_SCOPE=ref`: 2) | ROOT v1.1 428770 | `parts/data/data_428770_<job>.npz`: evnum, n_trk, lead_pt (= etap1 `data_events.py`) |
| `legs` | h5 v1.1 mb_nd + nominalne `scores/{D0_H,Dplus_H}/incl_mc.h5` (`BASE_CFG=config/lxplus_v11.yaml`) | `legs/legs_mb.parquet`, `legs/dstar_mb.parquet`, `legs_meta.json`, `mc_rates.json` — **kodem charmtag** `legs_for_truth` + `mc_rates` (snapshot `src/`), kontrola: 1344 D*, 9534/8665 D, 22/41/15 par |
| `combine` | wszystko powyżej | `w2d.npy`, `w2d_alldata.npy`, `w2d_edges.json`, `hist2d.npz`, `acceptance.json`, `mb_rates_2d_repro.json`, `mb_rates_2d_alldata.json`, `w2d_validation_repro.json`, `cw_h5.json`, `effmap_direct.npz`, `truthD_mb.parquet`, `legs_mb_w.parquet`, `dstar_mb_w.parquet`, `events_D_mb.parquet`, `closure_mb.json`, `response_nch.json`, `summary.md`, `run_info.json`, `DONE` |

(a) **w2**: `w2(n_trk, lead_pt) = p_dane/p_MC` dokładnie jak `etap1_definicje_masy_v20260925_01/code/w2d.py` (biny 15×13, komórki MC < 20 → 0).
Referencja = pliki danych `_000021`, `_000030` (441 691 zdarzeń; te same co etap1 i t4 — `kontrole_t1_t3_t4…/code/t4_rates.py`); wariant: wszystkie 39 plików.
**Test akceptacji** (`acceptance.json`): `analyze2.py 2d` odtworzony 1:1 → c_w (`weight_factors`) vs `ref/mb_rates_2d.json`:
0,935 / 0,793 / 1,374 / 1,201 / 0,981, tolerancja |Δ| < 0,005, plus liczniki 19 714 175 zdarzeń, 441 691 zdarzeń danych, 1314 D*, 9448/8671 D, 21/41/15 par — wszystkie równe → PASS.

(b) **mapy ε trybu direct**: nogi z h5 (selekcja `direct` bez zmian) dopasowane do prawdziwych D z `TruthDmesons` po (job, eventNumber, gatunek,
|Δη| < 0,05, |pT_r/pT_t − 1| < 0,15; D* także |Δφ| < 0,05). Na każdym prawdziwym D: `n_legs` (liczba nóg direct), `matched` (stara flaga Donly).
D*: noga = D* w **oknie fitu trybu direct** (`mc_rates.dsel`: |m−Δm−m_D0| < 40 MeV, Δm < 165 MeV, truth 0/1) — to jest „flaga okna D*”.
`effmap_direct.npz`: liczniki `<sp>__{den,num_any,num_mult,num_matched}` (nieważone i `w2__…`) w binach pT [2,5;3;4;5;6;8;12;30;∞) × |η| [0;0,5;1;1,5;2;2,2;2,5;2,7) × n_trk
[0;10;20;30;45;70;∞) (mianownik: mode == 1) + gotowe `eps_{any,mult,matched}__<sp>` z cofaniem binów jak `truth_combine.EffMap`.
ρ z fd: **nie liczone** (fd ma nieunikalne event_number → `check_event_key_unique` zatrzymałby nogi fd); używać istniejącego ρ fd_nd 1,060/1,032/1,020 jako wariantu,
a ρ selekcji direct sprawdza `closure_mb.json` (faktoryzacja na mb_nd, split-half).

(c) **n_ch → n_trk**: ntuple v1.1 mb_nd **nie mają TruthParticles** (podgląd `root_podglad_v20260929_01/mc_minbiasND_STRUKTURA.txt`: tylko TruthDmesons).
`gn_root.py` i tak sprawdza drzewa/gałęzie `TruthParticles*` w każdym pliku; brak → `response_nch.json = {"skipped": …}`; są → `response_nch.npz`.

## Diff względem istniejącego kodu
- `effmap_gn.py` = kopia `Q_lxplus/truth/effmap.py` + kolumna `evnum` (EventInfo_eventNumber) w wierszach prawdy (klucz joinu z h5); reszta bez zmian:
  ```
  +    evn = ak.to_numpy(a[k["evnum"]]).astype(np.uint64) if "evnum" in k else np.zeros(n_ev, np.uint64)
  +                evnum=evn[ev_of_truth][fsel],
  -             ntrk=resolve(keys, "EventInfo_nTracksPt500Eta25"))
  +             ntrk=resolve(keys, "EventInfo_nTracksPt500Eta25"), evnum=resolve(keys, "EventInfo_eventNumber"))
  ```
  Flagi okna D* nie dopisano do effmap (selekcja z ROOT nie byłaby „dokładnie direct”): okno D* bierze się z h5 kodem charmtag (`gn_legs.dstar_scan`
  = `scan_source` + `mc_rates.dsel`, z kontrolą liczby = `mc_rates.n_dstar`).
- `src/` = snapshot `charmtag` 2.1.0 (15.11) + lokalna zmiana `check_event_key_unique` w `pairs/build.py`, **bez zmian**
  (`src/SNAPSHOT.txt`). AA_merge nietknięte.
- `gn_root.py` (część MC) = logika `etap1 …/code/mb_truth.py` 1:1 (+ odczyt Tracks porcjami); `gn_combine.analyze2` = `analyze2.py 2d` 1:1
  (te same typy float32/int16) — test lokalny: wynik **bit w bit** równy oryginalnym skryptom etap1.

## Użycie (lxplus, `~/private/charm_q/Q_lxplus/gen_norm`; normalnie przez ORCH, etap `gen_norm`, domyślnie ZABLOKOWANY)
```
./submit_gn.sh check                         # ROOT 80+39 (w tym _000021/_000030), h5 v1.1, scores incl_mc, src, kolejka, DONE, AFS < 85 %
DRYRUN=1 ./submit_gn.sh submit               # gn.dag + komendy węzłów (nic nie wysyła)
./submit_gn.sh submit                        # condor_submit_dag -f -batch-name charm_gn -maxjobs 40 (czysty env: env -u LD_LIBRARY_PATH …)
./submit_gn.sh status | fetch
```
Zmienne: `GN_OUT` (domyślnie `/eos/user/k/klyczak/charm_q/gen_norm`), `DATA_SCOPE=all|ref`, `DATA_REF_TAGS=000021,000030`, `MAX_EVENTS=N` (próba),
`FORCE=1` (przelicz części; DONE → DONE.old_*), `NOMINAL`, `BASE_CFG`, `REF_RATES`, `MAXJOBS`. Części istniejące na EOS są pomijane (wznowienie po awarii).
Szacunek: MC 80 × ~3–5 min, dane 39 × ~5–10 min (Tracks 650 MB), legs ~0,5–1 h (80 h5 przez FUSE), combine ~10–20 min / 12 GB → ≈ 6–10 CPU·h, ściana ~1–2 h.
EOS: parts ≈ 0,5 GB, wyniki ≈ 30–80 MB (pobierane: bez `parts/`). AFS: tylko logi `logs/gen_norm/` (~120 × 2 pliki, < 5 MB).

## Test lokalny
`tests/run_local.sh [KATALOG] [N=3]` — atrapy ROOT z prawdziwych lokalnych h5 v1.2 mb_nd (FEAT_charmnet), oryginalne skrypty etap1 jako referencja,
łańcuch węzłów DAG (`make_dag_gn.sh` + `gn_job.sh`) z lokalnym `common.sh`, porównanie bit w bit, `tests/test_closure.py` (zabawka, truth/direct ≈ 1).
