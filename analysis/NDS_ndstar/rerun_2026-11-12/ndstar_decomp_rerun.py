"""rerun 2026-11-12: kopia PD_pairs_diff/ndstar/ndstar_decomp.py z poprawionym fitem (dstar_refit z PYTHONPATH, np. charmtag z poprawką NLL z 10.11).
Zmiany: import dstar_refit z PYTHONPATH (nie AA_merge), wynik do pliku podanego w argv[1] (nie nadpisuje ndstar_decomp.json).
Kształt NDS brany z ORYGINALNEGO ndstar_full.json (identyczny z dstar_post.json norm.shape_mc)."""
import sys, os, json, numpy as np
W = "/Users/k/Documents/charm/charm_work"
from charmtag import dstar_refit as R
print("dstar_refit:", R.__file__)
z = np.load(os.path.join(W, "V_daneMC/lxplus/dstar/dstar_nom.npz")); e = z["edges"]; h = z["norm_all"].astype(float)
nds = json.load(open(os.path.join(W, "NDS_ndstar/ndstar_full.json")))["shape_2g_mc"]
pd_ = json.load(open(os.path.join(W, "PD_pairs_diff/runs/lxplus_v12/pairs/pair_inputs.json")))["meta"]["dstar_fit_data"]["shape_mc"]
def rebin2(h, e):
    n = (len(h) // 2) * 2; return h[:n].reshape(-1, 2).sum(1), e[:n + 1:2]
out = {"dstar_refit_file": R.__file__}
for bn, (hh, ee) in {"0.05": (h, e), "0.10": rebin2(h, e)}.items():
    for sn, sh in {"NDS_shape(V MC)": nds, "PD_shape(pairs MC)": pd_}.items():
        f = R.fit_single(hh, ee, sh)
        out[f"{bn}|{sn}"] = {k: f[k] for k in ("N", "N_err", "scale", "mu", "valid", "chi2", "ndf", "n_entries", "bkg_model", "edm", "posdef", "attempts")}
        print(f"bin {bn}  {sn:20s} N={f['N']:7.1f} ± {f['N_err']:5.1f}  s={f['scale']:.3f} mu={f['mu']:.2f} valid={f['valid']} "
              f"chi2/ndf={f['chi2']:.0f}/{f['ndf']} n={f['n_entries']} bkg={f['bkg_model']} att={[(a['model'], a['valid'], round(a['N'],1)) for a in f['attempts']]}", flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
