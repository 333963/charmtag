# UBJK — jackknife μ po plikach danych na odsłoniętej próbce par D–D (12.11.2026)

**Po co.** Fit 2D par zakłada niezależny Poisson w komórkach, a pary z jednego zdarzenia są skorelowane (n partnerów na nogę → n wpisów;
`charm_work/BIN1885/bin1885.md` §4 i uwaga uboczna). Preregistracja (`T_PREP/unblinding_checklist.md` §2) przewidywała dla wyniku odsłoniętego
bootstrap/jackknife na poziomie zdarzeń. Jackknife „delete-one-file” (39 plików danych run 428770, zdarzenie nigdy nie jest dzielone między pliki)
mierzy σ(μ̂) z korelacjami w zdarzeniach i niejednorodnością plików; σ_jk/σ_Hesse > 1 ⇒ błąd z fitu zaniżony.

**Odsłonięcie.** Decyzja 05.11 (flaga `UNBLIND_OK`); węzeł sprawdza `$UB_OUT/UNBLIND_OK`, `$UB_OUT/DONE` i hash `pf_pairs/cfg/pf_unblind.yaml`
(= `cfg_sha256` z flagi), ustawia `CHARM_UNBLIND=1`. Nic nie jest zmieniane w odsłoniętym wyniku (wyjścia w osobnym `$UB_OUT/jk/`).

## Kod i wejścia
- `PYTHONPATH=$Q/pf_pairs/src` (przez `pf_env`; src NIE kopiowany), cfg `pf_pairs/cfg/pf_unblind.yaml`.
- `ubjk_run.py` importuje z charmtag (bez zmian snapshotu): `syst.files.{file_table, assign_files, per_file_hists, scaled_inputs, unblinded_jackknife}`,
  `pairs.fit2d.fit_channel`, `pairs.model`, `pairs.frozen.sbsb_counts`.
- Wejścia (EOS `charm_q/pf_pairs/unblind/pairs/`): `pair_inputs.{json,npz}`, `pairs_full.npz`, `fit2d_unblind.json`, `legs_data.parquet`
  (pełne nogi z `pairs --unblind`; plik → przez `cand` i liczby wierszy h5). Brak parquetu → prep przebudowuje pary w scratch
  (`cmd_pairs(unblind=True)`, 16 GB) i bierze tylko nogi; wejścia modelu zawsze z 05.11.

## DAG (`make_dag_ubjk.sh`)
`prep` → `fit_0 … fit_{N−1}` (+ `ref`) → `finish` (N = `UJ_NCHUNK`, domyślnie 8)
- **prep**: nogi → per plik: histogramy 2D (pełna płaszczyzna), `dataleg`, SB×SB. Kontrole (STOP przy niezgodności):
  Σ plików == `pairs_full.npz` bitowo; fit nominalny == `fit2d_unblind.json` (|Δμ| < 2·10⁻³). Σ dataleg/SB×SB == `pair_inputs` → `rederive_ok`.
- **fit_i**: repliki `[-1, 0..38][i::N]` (−1 = nominał z przedziałem 68 % wspólnego μ). Per replika i tryb: 3 × `fit_channel` (μ̂, Hesse, Z)
  + wspólne μ = ciągłe minimum Σ profili kanałów (Brent, ~15–25 × 3 fity z ustalonym μ) i Z = √(2[ΣNLL(0) − min]).
  Tryby: `fixed` = `files.scaled_inputs(I, 1−w_k)` (definicja `unblinded_jackknife`), `rederive` = też dataleg/SB×SB/nev_legs bez pliku k.
- **ref**: wprost `charmtag.syst.files.unblinded_jackknife` (odniesienie: tryb fixed musi dać to samo σ_jk).
- **finish**: `ubjk.json` + `ubjk.md`: μ̂ bez każdego pliku, σ_jk (standard i ważony Busing 1999 — pliki 109–248 k zdarzeń),
  σ_jk/σ_Hesse, σ_jk/σ_profil68, μ_J (bias-corr.), Z min–max, korelacje kanałów, wspólne μ/Z; `DONE`.

Koszt (pomiar lokalny fitów: D0D0 3,6 s, DpD0 5,1 s, DpDp 0,3 s + fit μ=0): kanały ~10–20 s/replikę/tryb, wspólne μ ~20–25 × 3 fity
z ustalonym μ ≈ 2–4 min → 79 replik·trybów ≈ 3–5 h CPU; 8 porcji (domyślnie) → ~0,5–1 h każda (workday 8 h). prep ~5–15 min
(z przebudową ~40 min), ref ~10–15 min, finish < 1 min.

Opcja `UJ_GROUPS=events:G` (np. 200; wyjście `unblind/jk_ev200/`): grupy zdarzeń z hasha (run, event) zamiast plików — zdarzenie nadal
niepodzielone, precyzja σ_jk ±1/√(2(G−1)) (G=200: ±5 %), ale bez niejednorodności między plikami; bez węzła ref. Koszt ∝ G
(G=200, tylko `UJ_MODES=fixed`: ~200 × 2–4 min ≈ 7–13 h CPU → `UJ_NCHUNK=25`).

## Użycie (lxplus, `~/private/charm_q/Q_lxplus/ub_jk`)
`./submit_ubjk.sh check|submit|status|fetch` — submit w czystym env (`env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:… condor_submit_dag -f
-batch-name charm_ubjk`). Zmienne: `UJ_GROUPS` (files), `UJ_NCHUNK` (8), `UJ_MODES` (fixed,rederive), `UJ_REF` (1), `UJ_REBUILD` (auto). `DRYRUN=1 ./submit_ubjk.sh submit` = tylko komendy.
Logi AFS `Q_lxplus/logs/ubjk/`, wyniki tylko EOS `charm_q/pf_pairs/unblind/jk/`. ORCH: `stages/ub_jk.sh` → `PF2_lxplus/runs/lxplus/unblind/jk/`.

Lokalny test: `CHARM_UNBLIND=1 PYTHONPATH=Q_lxplus/pf_pairs/src ../.venv_h5/bin/python Q_lxplus/ub_jk/ubjk_run.py selftest
--pairs-dir PF2_lxplus/runs/lxplus/unblind/pairs --nfiles 10 --out <scratch>` (z `charm_work/`; wielomianowy podział → oczekiwane σ_jk/σ ≈ 1).

## Ograniczenia
- σ_jk z 39 plików: precyzja ±1/√(2·38) ≈ ±11 % (gauss.); różnice < ~20 % nieistotne.
- N_D* (normalizacja) i wejścia MC skalowane ∝ zdarzeniom — ich fluktuacja należy do systematyki normalizacji.
- Wagi plików = `n_events_read` z h5 (jak `syst files`); diagnostyka: korelacja z liczbą zdarzeń z nogami (`perfile.json: w_corr`).
