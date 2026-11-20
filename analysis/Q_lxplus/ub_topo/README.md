# ub_topo — kształt |Δη| sygnału par D–D vs Pythia/DPS (projekt A, 19.11.2026)

Decyzje i prognoza mocy: `PROGNOZA.md`. Lokalnie: `topo_gen.py` (hipotezy z generatora → `ref/topo_gen.json`), `topo_report.py --selftest`,
`topo_report.py --forecast-from PD_pairs_diff/runs/lxplus_v12/pairs` (→ `ref/forecast.json`). PYTHONPATH z `AA_merge/src` (diff.py identyczne md5 z `pf_pairs/src`).
lxplus (ORCH `stages/ubtopo.sh`): `submit_topo.sh check|submit|status` → DAG pairs → fit2d → report (batch `charm_ubtopo`),
mechanika `pd_diff` (pd_common.sh, tools/diff_fill.py; tylko odczyt), cfg `cfg/ub_topo.yaml` = pf_unblind.yaml + out_dir + pairs.diff.
Wyjście EOS `charm_q/ub_topo/dy/{UBTOPO.md, ubtopo.json, tfin_diff.md, DONE}` → `PF2_lxplus/runs/lxplus/unblind/ubtopo`.
