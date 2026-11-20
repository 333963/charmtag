# prez_net — PREZ-2 (19.11.2026, projekt A): klasyfikator pojedynczych D, BDT vs sieć × high-level vs low-level

Prezentacja wg kartki, punkt 3 („też sieć, nie tylko drzewo”). TYLKO MC (fd_nd sygnał truth, mb_nd tło kombinatoryczne truth 8–99), h5 v1.1,
config i foldy nominalnego XGB H (`$RUNS/models/D0_H/config_used.yaml`, fold_of(run, event), 5). Opis modeli do slajdu: `model.md`.

DAG (`submit_pz.sh check|submit|status`, batch `charm_preznet`, wynik EOS `charm_q/prez_net`):
  ext × 4 (`pz_ext.py`: D0/D+ × fd/mb → npz: 18 cech H, ślady [n, nd, 11], L_xy, masa, wagi, fold, score nominału OOF)
  → train × 4 (`pz_train.py`: gatunek × {xgb, net}; każdy liczy hi i lo, 5 foldów OOF, ważność permutacyjna + SHAP)
  → report (`pz_report.py`: ROC, AUC, ε_B @ WP90, ε_S @ ε_B(nominał), rozkłady score, ważność, rzeźbienie masy; TABLE.md, PNG+PDF).
Etap ORCH: `ORCH/stages/preznet.sh` → `charm_work/PREZ_2026-11-19/p3/runs`.

Decyzje:
- trening i ewaluacja w oknie sygnału (|m − m0| < 55/45 MeV): sieć na śladach umie policzyć masę niezmienniczą, a sygnał leży tylko w piku —
  trening „pik vs pasma” nauczyłby masy (rzeźbienie). Tło poza oknem służy tylko do testu rzeźbienia (ε_B w oknie / poza oknem przy WP90).
- tło z MC (nie pasma danych jak nominał): ta sama definicja dla 4 wariantów, prawda dostępna; nominał (trenowany na pasmach danych) jako odniesienie.
- Deep Sets dla śladów (zbiór nieuporządkowany, 2 piony D+ wymienne, rola K/π jako cecha), MLP dla 18 cech.

Testy lokalne (19.11): `pz_selftest.py` (syntetyczne npz; xgb + net z torch CPU) OK; łańcuch na syntetycznych h5 v1.1 (`charmtag.synth`)
ext → train (xgb, net) → report OK (wszystkie wykresy). Pełny trening sieci (selftest, 20 k): MLP AUC 0,983 ≈ XGB; Deep Sets 0,915 vs XGB-lo 0,908.
Niesprawdzone lokalnie: czytanie prawdziwych h5 v1.1 z EOS i modele foldów nominału (sprawdza `submit_pz.sh check` + pierwszy węzeł ext).
