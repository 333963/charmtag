"""ndstar: rozkład różnicy N_D* NDS (4443) vs bieg PD (4022) na kształt MC i binning.
Lokalnie, tylko odczyt: V dstar_nom.npz (norm_all, 0,05 MeV) + dstar_refit.fit_single z AA_merge/src.
Użycie: .venv_h5/bin/python ndstar_decomp.py  (ścieżki względem charm_work)."""
import sys, os, json, numpy as np
W = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(W, "AA_merge", "src"))
from charmtag import dstar_refit as R
z = np.load(os.path.join(W, "V_daneMC/lxplus/dstar/dstar_nom.npz")); e = z["edges"]; h = z["norm_all"].astype(float)
nds = json.load(open(os.path.join(W, "NDS_ndstar/ndstar_full.json")))["shape_2g_mc"]
pd_ = json.load(open(os.path.join(W, "PD_pairs_diff/runs/lxplus_v12/pairs/pair_inputs.json")))["meta"]["dstar_fit_data"]["shape_mc"]
def rebin2(h, e):
    n = (len(h) // 2) * 2; return h[:n].reshape(-1, 2).sum(1), e[:n + 1:2]
out = {}
for bn, (hh, ee) in {"0.05": (h, e), "0.10": rebin2(h, e)}.items():
    for sn, sh in {"NDS_shape(V MC)": nds, "PD_shape(pairs MC)": pd_}.items():
        f = R.fit_single(hh, ee, sh); f.pop("attempts", None)
        out[f"{bn}|{sn}"] = {k: f[k] for k in ("N", "N_err", "scale", "mu", "valid", "chi2", "ndf", "n_entries", "bkg_model")}
        print(f"bin {bn}  {sn:20s} N={f['N']:7.1f} ± {f['N_err']:5.1f}  s={f['scale']:.3f} mu={f['mu']:.2f} valid={f['valid']} "
              f"chi2/ndf={f['chi2']:.0f}/{f['ndf']} n={f['n_entries']} bkg={f['bkg_model']}")
# efektywna szerokość sygnału (RMS 2G) — skala zmiany N przy tym samym tle
for sn, sh in {"NDS": nds, "PD": pd_}.items():
    rms = sh["sc"] * np.sqrt(sh["f"] + (1 - sh["f"]) * sh["k"] ** 2)
    print(f"{sn}: sc={sh['sc']:.3f} k={sh['k']:.2f} f={sh['f']:.3f} σ_wide={sh['sc']*sh['k']:.2f} RMS={rms:.3f} MeV")
json.dump(out, open(os.path.join(os.path.dirname(__file__), "ndstar_decomp.json"), "w"), indent=1)
