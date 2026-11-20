#!/usr/bin/env python3
"""UB427 N_D* (18.11.2026): N_D* w 427929 metodą NDS — ta sama, którą znormalizowano nominał 428770 (4436 ± 671, NDS_ndstar/rerun_2026-11-12):
D* z h5 po filtrze (= D_d potoku par: decayType 1, 0 < Δm < 0,165 GeV, |m(D*) − Δm − m_D0| < 40 MeV, bez cięcia score), histogram Δm [MeV]
na binach V (139,6–170, 0,05 MeV), `dstar_refit.fit_single` z poprawką NLL (charmtag z 10.11 = AA_merge; kopia dstar_refit_nds.py)
i kształtem 2G z MC V (shape_2g_mc z NDS). Bieg par w pf_pairs/src liczy N_D* starszym kodem (fallback fit_dm, +25 %) — μ przeskalowuje raport.
  u4_ndstar.py 'h5glob' --out ndstar_427929.json"""
import argparse
import glob
import json
import sys
import os

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dstar_refit_nds as R  # noqa: E402

SHAPE_NDS = {"mu": 145.54556094560178, "sc": 0.7428973985874197, "k": 2.728596842934596, "f": 0.7182675043379269}
N_NOM_428770 = {"N": 4436.390108391369, "N_err": 670.8561890721969, "src": "NDS_ndstar/rerun_2026-11-12/new/ndstar_full.json variants[0]"}
EDGES = np.round(np.arange(139.6, 170.0 + 1e-9, 0.05), 6)
M0_D0 = 1.86484
PT_BINS = {"3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("h5"); ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    fs = sorted(glob.glob(a.h5))
    dm, pt = [], []
    for f in fs:
        with h5py.File(f, "r") as g:
            e = g["events"]
            dt = e["Dmeson_decayType"][:]
            q = dt == 1
            if not q.any():
                continue
            x = e.fields(["Dmeson_delta_mass", "Dmeson_mass", "Dmeson_pt"])[:][q]
        s = (x["Dmeson_delta_mass"] > 0) & (x["Dmeson_delta_mass"] < 0.165) & (np.abs(x["Dmeson_mass"] - x["Dmeson_delta_mass"] - M0_D0) < 0.040)
        dm.append(x["Dmeson_delta_mass"][s] * 1e3); pt.append(x["Dmeson_pt"][s])
    dm = np.concatenate(dm) if dm else np.zeros(0); pt = np.concatenate(pt) if pt else np.zeros(0)
    h, _ = np.histogram(dm, EDGES)
    r = R.fit_single(h, EDGES, SHAPE_NDS); r.pop("attempts", None)
    out = {"files": len(fs), "n_cand_window": int(len(dm)), "n_entries_hist": int(h.sum()), "fit": r, "shape": SHAPE_NDS, "nominal_428770": N_NOM_428770,
           "method": "NDS: fit_single (NLL poprawione), kształt 2G MC V, biny 0,05 MeV, okno D0 ±40 MeV", "per_bin": {}}
    for b, (lo, hi) in PT_BINS.items():
        hb, _ = np.histogram(dm[(pt >= lo) & (pt < hi)], EDGES)
        try:
            rb = R.fit_single(hb, EDGES, SHAPE_NDS); rb.pop("attempts", None)
        except Exception as ex:  # noqa: BLE001
            rb = {"error": repr(ex)}
        out["per_bin"][b] = rb
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print(f"N_D* 427929 (NDS) = {r['N']:.1f} ± {r['N_err']:.1f}, valid {r['valid']}, model {r.get('bkg_model')}, wejść {int(h.sum())}; "
          + ", ".join(f"{b}: {v.get('N', float('nan')):.0f}" for b, v in out["per_bin"].items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
