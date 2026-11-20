# PF2 — powtórka par D–D z poprawką weta wspólnych śladów (lxplus/HTCondor)

**Co:** `pairs → fit2d (ślepo) → 50×200 toys → toys_merge` [T] + `syst core/files/wp-build×6/wp/norm/report` [U] + `mix_pairs → mix_fit2d` [S5 nom_mix] + `compare`.
**Bez treningu**: modele/score/validation z nominalnego biegu (`$RUNS`), kod = `src/` (AA_merge + PF: `pairs.track_id: d0z0`), przez PYTHONPATH (editable `charmtag_src` głównego DAG nietknięty).
**Wyjścia:** `/eos/user/k/klyczak/charm_q/pf_pairs/{pairs,pairs/syst,pairs/toys,mix,compare_pf.md,DONE}`; AFS tylko logi `logs/pf/`; scratch kasowany (trap).
**Zawsze ślepo** (bez `--unblind`). NORM=direct domyślnie (porównywalne z ½ lumi/wrześniem); `NORM=truth` gdy `pairs.truth_norm_json` w `config/lxplus_v11.yaml`.

Mac: dwuklik `~/Documents/charm/START_lxplus_pf.command [check|submit|status|fetch]` (domyślnie: check, a gdy OK → submit; gdy na EOS jest DONE → fetch).
lxplus ręcznie: `cd ~/private/charm_q/Q_lxplus/pf_pairs && ./submit_pf.sh check|submit|status|fetch`.
Warunki `check`: główny DAG nie w kolejce (`IGNORE_MAIN=1` gdy zostały tylko toys), nominalne pliki na EOS, PF nie w kolejce, brak DONE, AFS < 85 %, src z poprawką.
Opcje `make_dag_pf.sh`: `NOTOYS=1 NOSYST=1 NOMIX=1 NTOYJOBS= NTOY=`.

`legs_regen` (nogi frozen z września nowym kluczem dla `--closure`): POMINIĘTE — na lxplus `inputs/` ma tylko `mb_rates_2d.json`, `sf_for_pairs_v2.json`; frozen `legs_*.parquet` są tylko na Macu (`unblinding_freeze_v20260927_01/frozen_inputs/`), a `--closure` i tak nie działa na lxplus (komentarz w lxplus_v11.yaml). Closure z nowym kluczem → osobny bieg na Macu (`charmtag pairs --closure` po regeneracji nóg z h5 — wymaga pełnych h5 lokalnie lub na lxplus z frozen skopiowanym do `inputs/`).

Szacunek czasu: pairs ~1–2 h (2 skany h5 po 8 GB/próbkę), fit2d 10 min, toys 50 jobów × ~4 h równolegle (~5 h), syst wp-build 6 × ~1–2 h, mix 1–2 h → **~6–8 h od submitu** (toys dominują; `NOTOYS=1` → ~3 h).
