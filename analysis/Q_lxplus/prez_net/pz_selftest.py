#!/usr/bin/env python3
"""PREZ-2 selftest: syntetyczne npz w formacie pz_ext (sygnał/tło z rozróżnialnymi rozkładami w kilku cechach śladów, cechy H liczone
ze śladów) → pz_train (xgb, net jeśli jest torch) → pz_report. Sprawdza: przebieg, AUC > 0,6 dla wszystkich wariantów, brak NaN w OOF w oknie,
spójność wierszy, wykresy. Użycie: pz_selftest.py --out DIR [--n 30000] [--no-net]"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def make(sp, sample, n, seed, out):
    import pz_ext as E
    rng = np.random.default_rng(seed)
    nd = 2 if sp == "D0" else 3; m0 = 1.86484 if sp == "D0" else 1.86966; sw = 0.055 if sp == "D0" else 0.045
    lo, hi = (1.75, 2.00) if sp == "D0" else (1.70, 2.20)
    sig = sample == "fd"
    m = rng.normal(m0, 0.015, n) if sig else rng.uniform(lo, hi, n)
    m = np.clip(m, lo, hi - 1e-6)
    pt = 3 + rng.exponential(2.0 if sig else 1.6, n); pt = np.clip(pt, 3, 29.9)
    eta = rng.uniform(-2.5, 2.5, n)
    T = np.zeros((n, nd, E.NTF), np.float32)
    for s in range(nd):
        tpt = pt / nd * rng.uniform(0.5, 1.5, n); teta = eta + rng.normal(0, 0.3, n)
        T[:, s, 0] = np.log(tpt); T[:, s, 1] = teta; T[:, s, 2] = rng.normal(0, 0.25, n); T[:, s, 3] = teta - eta
        T[:, s, 4] = np.log(tpt * np.cosh(teta)); ded = rng.normal(1.6 if (sig and s == 0) else 1.2, 0.3, n)
        miss = rng.random(n) < 0.05; T[:, s, 5] = np.where(miss, np.nan, ded); T[:, s, 6] = miss
        T[:, s, 7] = rng.integers(1, 5, n); T[:, s, 8] = rng.normal(0, 0.08 if sig else 0.04, n); T[:, s, 9] = rng.normal(0, 0.3, n)
        T[:, s, 10] = float(s == 0)
    Lxy = rng.exponential(1.2 if sig else 0.7, n) + 0.25
    ad0 = np.abs(T[:, :, 8])
    Xh = np.stack([np.log1p(pt), np.abs(eta), rng.chisquare(2 if sig else 4, n), Lxy, rng.normal(3 if sig else 1, 1.5, n), np.abs(rng.normal(0, 1, n)),
                   ad0.max(1), rng.exponential(3 if sig else 1.5, n), np.abs(rng.uniform(-1, 1, n)) * (0.8 if sig else 1), np.abs(rng.normal(3 if sig else 1, 1.5, n)),
                   ad0.min(1), np.abs(T[:, :, 9]).max(1), np.abs(T[:, :, 1]).max(1), T[:, 0, 5], T[:, 0, 4], T[:, 0, 7], np.nanmax(T[:, 1:, 5], 1),
                   T[:, 1:, 4].min(1)], 1).astype(np.float32)
    run = np.zeros(n, np.int64); event = rng.integers(0, 2**40, n)
    from charmtag import samples as S
    fold = S.fold_of(run, event, 5)
    feats = ["log_pt", "abs_eta", "chi2", "Lxy", "impact_sig", "abs_impact_z0", "d_max_abs_d0", "pointing", "abs_cos_theta_star", "abs_impact_sig",
             "d_min_abs_d0", "d_max_abs_z0", "d_max_abs_eta", "dedx_K", "log_p_K", "npts_K", "dedx_pi_max", "log_p_pi_min"]
    meta = {"species": sp, "sample": sample, "n": n, "n_inwin": int((np.abs(m - m0) < sw).sum()), "features_hi": feats, "track_features": E.TRK_NAMES,
            "m0": m0, "signal_window": sw, "mass_range": [lo, hi], "n_folds": 5,
            "reweight_bins": {"pt": [3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 12, 15, 20, 30], "abs_eta": [0, 0.5, 1.0, 1.5, 2.0, 2.7]}, "s_nom_available": True}
    s_nom = 1 / (1 + np.exp(-(0.8 * (Lxy - 1) + 0.5 * (Xh[:, 4] - 2) + rng.normal(0, 1, n)))).astype(np.float32)
    np.savez_compressed(out, meta=json.dumps(meta), Xh=Xh, T=T, G=Lxy[:, None].astype(np.float32), m=m, pt=pt, abs_eta=np.abs(eta), weight=np.ones(n),
                        run=run, event=event, truth=np.zeros(n) if sig else np.full(n, 20.0), fold=fold, inwin=np.abs(m - m0) < sw, s_nom=s_nom)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--n", type=int, default=30000)
    ap.add_argument("--no-net", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    try:
        import torch  # noqa: F401
        have_torch = not a.no_net
    except ImportError:
        have_torch = False
    dirs = []
    for i, sp in enumerate(("D0", "Dplus")):
        d = os.path.join(a.out, sp); os.makedirs(d, exist_ok=True); dirs.append(d)
        make(sp, "fd", a.n, 10 + i, os.path.join(d, "sig.npz")); make(sp, "mb", 2 * a.n, 20 + i, os.path.join(d, "bkg.npz"))
        for fam in (["xgb", "net"] if have_torch else ["xgb"]):
            r = subprocess.run([sys.executable, os.path.join(HERE, "pz_train.py"), "--sig", os.path.join(d, "sig.npz"), "--bkg", os.path.join(d, "bkg.npz"),
                                "--family", fam, "--out", d, "--threads", "4", "--quick"], capture_output=True, text=True)
            print(r.stdout[-1500:], r.stderr[-3000:])
            assert r.returncode == 0, f"pz_train {fam} {sp} padł"
    r = subprocess.run([sys.executable, os.path.join(HERE, "pz_report.py"), *sum([["--in", d] for d in dirs], []), "--out", os.path.join(a.out, "report")],
                       capture_output=True, text=True)
    print(r.stdout[-3000:], r.stderr[-3000:])
    assert r.returncode == 0, "pz_report padł"
    T = json.load(open(os.path.join(a.out, "report", "table.json")))
    ok = True
    for sp, t in T.items():
        for k, row in t["rows"].items():
            if not (row["auc"] > 0.6):
                print("ZA NISKIE AUC", sp, k, row["auc"]); ok = False
    for d in dirs:
        z = np.load(os.path.join(d, "scores_xgb.npz"))
        for k in ("s_hi", "s_lo"):
            if not np.isfinite(z[k][z["inwin"]]).all():
                print("NaN w OOF", d, k); ok = False
    figs = [f for f in os.listdir(os.path.join(a.out, "report")) if f.endswith(".png")]
    print("wykresy:", sorted(figs)); print("torch:", have_torch)
    print("SELFTEST:", "OK" if ok and len(figs) >= 8 else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
