# NDS — patch dla `Q_lxplus/pf_pairs/src` (NIE nakładać w trakcie biegu `unblind.dag`)

`Q_lxplus/pf_pairs/src/charmtag/` celowo zostało na 2.0.2 z poprawką ID śladu z 5.11 (hash cfg/src w `UNBLIND_OK` na EOS). Węzeł `unblind`
liczy N_D* starym fallbackiem `fit_dm` (5588) — korekta post hoc: `UNB_DRY/fill_unblind.py --ndstar NDS_ndstar/ndstar_full.json`.

Po zakończeniu unblind (DONE na EOS, fetch na Mac), żeby pf_pairs/src = kod NDS z 5.11:
```bash
cd ~/Documents/charm/charm_work
cp AA_merge/src/charmtag/dstar_refit.py Q_lxplus/pf_pairs/src/charmtag/dstar_refit.py
cp AA_merge/src/charmtag/pairs/build.py  Q_lxplus/pf_pairs/src/charmtag/pairs/build.py   # tylko docstring n_dstar_data
# albo: cd Q_lxplus/pf_pairs/src && patch -p2 < ../../../NDS_ndstar/patch_pf_pairs_src_nds.diff
```
Potem `submit_unblind.sh check` zgłosi zmianę hasha src — to oczekiwane; nowy bieg par (jeśli będzie) ma już N_D* z `fit_single`
(valid=True, bkg_model=threshold_pol2) i nie potrzebuje `--ndstar`.
