"""rerun 2026-11-12: czułość numeryczna fit_single (nominał) na fizycznie pomijalne zaburzenia: przesunięcie krawędzi o δ MeV.
Użycie: PYTHONPATH=<AA_merge/src | REPO/CHARM_git> python jitter_test.py OUT.json"""
import sys, json, numpy as np
from charmtag import dstar_refit as R
W = "/Users/k/Documents/charm/charm_work"
z = np.load(W + "/V_daneMC/lxplus/dstar/dstar_nom.npz"); e = z["edges"].astype(float); h = z["norm_all"].astype(float)
shape = json.load(open(W + "/V_daneMC/lxplus/summary/dstar_post.json"))["norm"]["shape_mc"]
out = []
for d in (0.0, 1e-9, -1e-9, 1e-7, -1e-7, 1e-6):
    f = R.fit_single(h, e + d, shape)
    out.append({"shift_mev": d, "N": f["N"], "N_err": f["N_err"], "valid": f["valid"], "edm": f["edm"], "bkg_model": f["bkg_model"]})
    print(R.__file__.split("/")[-4], out[-1], flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
