"""VER_budget: niezależne odtworzenie κ_cls (fd_nd, v11 H) + test mechanizmu: κ przewidziane z ε(pT) faktoryzowanego vs obserwowane.
Także κ dla pełnego zestawu kategorii D⁰ do WP90 (bez tagu), i bootstrap po zdarzeniach."""
import os, json, numpy as np, pandas as pd, h5py
SC = "/home/k/Documents/charm/charm_work/L_charmtag/runs/v11_full/scores"
E = {"D0": (1.75, 2.0), "Dplus": (1.70, 2.20)}
CH = {"D0D0": ("D0", "D0"), "DpD0": ("Dplus", "D0"), "DpDp": ("Dplus", "Dplus")}
WPmax = {"D0": 0.65, "Dplus": 0.65}
legs = {}
for sp in ("D0", "Dplus"):
    with h5py.File(os.path.join(SC, f"{sp}_H", "signal_mc.h5")) as f:
        d = pd.DataFrame({k: f[k][:] for k in ("row", "run", "event", "m", "pt", "eta", "score", "truth", "veto_ok")})
    d = d[np.isin(d.truth, [0, 1]) & (d.m >= E[sp][0]) & (d.m < E[sp][1])]
    if sp == "Dplus":
        d = d[d.veto_ok.astype(bool)]
    else:
        d = d.assign(k=d.pt.round(4).astype(str) + d.eta.round(4).astype(str)).sort_values("score", ascending=False)
        d = d.drop_duplicates(["run", "event", "k"]).drop(columns="k")
    s = np.sort(d.score.to_numpy())
    for e in (0.65, 0.9):
        d[f"p{int(e*100)}"] = d.score >= np.quantile(s, 1 - e)
    pb = np.quantile(d.pt, np.linspace(0, 1, 21)); pb[0] = 0; pb[-1] = 1e9
    d["pbin"] = np.clip(np.digitize(d.pt, pb) - 1, 0, 19)
    legs[sp] = d.reset_index(drop=True)
res = {}
rng = np.random.default_rng(1)
for ch, (i, j) in CH.items():
    A, B = legs[i], legs[j]
    P = A.merge(B, on=["run", "event"], suffixes=("1", "2"))
    P = P[np.abs(P.row1 - P.row2) < 400]
    if i == j:
        P = P[P.row1 < P.row2]
    out = {"n_pairs": int(len(P))}
    for wp in ("p65", "p90"):
        eA, eB = A[wp].mean(), B[wp].mean()
        obs = (P[wp + "1"] & P[wp + "2"]).mean()
        epA = A.groupby("pbin")[wp].mean(); epB = B.groupby("pbin")[wp].mean()
        pred = (P.pbin1.map(epA) * P.pbin2.map(epB)).mean()
        # bootstrap po zdarzeniach
        ev = P.run.astype(str) + ":" + P.event.astype(str)
        g = pd.factorize(ev)[0]; ng = g.max() + 1
        num = np.bincount(g, (P[wp + "1"] & P[wp + "2"]).to_numpy(float), ng); cnt = np.bincount(g, None, ng)
        bs = []
        for _ in range(300):
            w = rng.poisson(1.0, ng)
            bs.append((w @ num) / (w @ cnt) / (eA * eB))
        out[wp] = {"kappa_obs": obs / (eA * eB), "kappa_obs_err_boot": float(np.std(bs)), "kappa_pT_factorized": pred / (eA * eB),
                   "kappa_residual": obs / pred}
    out["pt_mean_pair"] = [float(P.pt1.mean()), float(P.pt2.mean())]; out["pt_mean_single"] = [float(A.pt.mean()), float(B.pt.mean())]
    res[ch] = out
print(json.dumps(res, indent=1, default=float))
json.dump(res, open("/home/k/Documents/charm/charm_work/VER_budget/kcls_check.json", "w"), indent=1, default=float)
