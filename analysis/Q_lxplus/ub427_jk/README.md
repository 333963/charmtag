# UB427JK — jackknife μ po grupach zdarzeń dla runu 427929 (19.11.2026, projekt A)

**Po co.** UB427 (18.11) podał σ(μ₄₂₇₉₂₉) z profilu 68 % bez korelacji par w zdarzeniu; w 428770 jackknife po 200 grupach zdarzeń dał ×1,18.
Ten bieg mierzy ten czynnik wprost w 427929 → wejście `--ub427jk` do `TRIG_A/syst_trig.py` (prerejestracja `TRIG_A/SYST_TRIG.md`, d₁).

**Jak.** `../ub_jk/ubjk_run.py` bez zmian: `prep --groups events:200` (hash(run, event) mod 200), `fit` (tryb fixed, 25 porcji), `finish`.
Wejścia modelu = wynik UB427 na EOS `charm_q/ub_427929/pairs/` (pair_inputs.{json,npz}, pairs_full.npz, fit2d_unblind.json).
Nóg 427929 nie zapisano w UB427 → prep przebudowuje pary (`cmd_pairs --unblind`) z tym samym cfg co węzeł `pairs` UB427
(merge_cfg z `pf_pairs/cfg/pf_unblind.yaml`, h5 `ub_427929/h5/*_mb.h5`, scores MC nominalne + dane z `ub_427929/scores`, `cw/mb_rates_427929.json`)
i zapisuje `legs_data_rebuilt.parquet` na EOS (powtórka prep go używa). Kontrole ubjk (STOP): Σ grup == pairs_full.npz bitowo,
fit nominalny == fit2d_unblind.json 427929 (|Δμ| < 2·10⁻³).

**Użycie (lxplus `~/private/charm_q/Q_lxplus/ub427_jk`):** `./submit_u7.sh check|submit|status|dag` (czysty env, `-batch-name charm_ub427jk`).
Zmienne: `U7_G` (200), `U7_NCHUNK` (25), `U7_MODES` (fixed), `U7_REBUILD=1` (wymuś przebudowę), `FORCE=1`. Wyjście EOS `charm_q/ub_427929/jk_ev200/`,
logi AFS `Q_lxplus/logs/ub427jk/`. ORCH: `stages/ub427jk.sh` → `PF2_lxplus/runs/lxplus/unblind/r427929/jk_ev200/` (etap NIE dopisany do pipeline.yaml).
Koszt: prep ~20–40 min (16 GB), fit 201 replik × ~2–4 min / 25 porcji ≈ 0,5 h każda, finish < 1 min.

**Test lokalny (19.11):** `ubjk_run.py selftest` na wejściach 427929 (wielomianowy podział na 8 „plików”, oczekiwane σ_jk/σ ≈ 1) — wynik w raporcie studium;
`DRYRUN=1 ./u7_job.sh prep|fit 3|finish` i `submit_u7.sh dag` (27 węzłów) OK.
