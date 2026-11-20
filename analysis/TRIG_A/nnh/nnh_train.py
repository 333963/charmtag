#!/usr/bin/env python3
"""TRIG_A/nnh — sieć neuronowa zamiast drzew? Porównanie na TYCH SAMYCH 18 cechach H, tych samych kandydatach i foldach (projekt A).

Próbka: D0 z par D0←D* (wyniki b_trig: bt_fd_nd, bt_mb_nd — MC; bt_data_* — dane). Pojedyncze D, ŻADNYCH par D–D.
Modele (OOF, 5 foldów z kolumny `fold`; model k nie widział foldu k; dane oceniane modelem swojego foldu):
  nn64   — MLP 3×64 ReLU, Adam lr 0,01, 20 epok, batch 10 000 (przepis CharmNet bez członu dCor), standaryzacja
  nn256  — MLP 4×256 SiLU, dropout 0,1, AdamW lr 2e-3 (cosinus), wd 1e-4, 60 epok, batch 2048, wejścia po transformacji rang → Gauss
  xgb    — XGBoost na tej samej próbce (kontrola „ta sama statystyka treningu”)
  avg    — średnia nn256 i xgb
  H      — istniejący score nominalny (kolumna `score`; trenowany na pełnej próbce inkluzywnej — punkt odniesienia, nie ta sama statystyka)
Sygnał = truth 0/1, tło = truth 8–15; reszta wyłączona. Metryki MC w oknie |m − m0| < 55 MeV per bin pT: AUC, ε_B przy ε_S = 30/50/65/90 %,
Z/Z_H = √(ε_B^H/ε_B), bootstrap Poissona (sparowany), corr(score, m | tło). Dane: ε_B z pasma bocznego Δm, SF z fitu Δm pass/fail (WP50, WP90).
Wyjście (--out): nnh.json, NNH.md, nnh_scores_<próbka>.npz.
"""
import argparse, glob, json, os, sys, time
from multiprocessing import Pool
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import btrig_common as C  # noqa: E402

FEATS = ["log_pt", "abs_eta", "chi2", "Lxy", "impact_sig", "abs_impact_z0", "d_max_abs_d0", "pointing", "abs_cos_theta_star", "abs_impact_sig",
         "d_min_abs_d0", "d_max_abs_z0", "d_max_abs_eta", "dedx_K", "log_p_K", "npts_K", "dedx_pi_max", "log_p_pi_min"]
BINS = {"3-5": (3.0, 5.0), "5-8": (5.0, 8.0), "8-30": (8.0, 30.0)}
WPS = (0.3, 0.5, 0.65, 0.9)
SFWP = (0.5, 0.9)
MODELS = ["H", "xgb", "nn64", "nn256", "avg"]
KEEP = ("run", "event", "row", "file_idx", "m", "pt", "dm", "score", "fold", "truth", "truth_ds")
T0 = time.time()
log = lambda *a: print(f"[{time.time() - T0:6.0f}s]", *a, flush=True)  # noqa: E731


def load(paths, max_rows=0):
    cols = {}
    for p in paths:
        z = np.load(p, allow_pickle=False)
        for k in z.files:
            if k in KEEP or k in FEATS or k.startswith("e_tp_") or k == "e_grl":
                cols.setdefault(k, []).append(z[k][:max_rows] if max_rows else z[k])
    miss = [f for f in FEATS if f not in cols]
    if miss:
        raise SystemExit(f"STOP: brak cech {miss} w {paths[0]}")
    D = {k: np.concatenate(v) for k, v in cols.items()}
    D["X"] = np.stack([D.pop(f).astype(np.float32) for f in FEATS], axis=1)
    return D


def prep(X, med=None):
    X = np.array(X, np.float32, copy=True); bad = ~np.isfinite(X)
    if med is None:
        med = np.nanmedian(np.where(bad, np.nan, X), axis=0); med = np.where(np.isfinite(med), med, 0.0)
    if bad.any():
        X[bad] = np.broadcast_to(med.astype(np.float32), X.shape)[bad]
    return X, med


class RankGauss:
    def __init__(self, X, nq=1000):
        from scipy.special import erfinv
        q = np.linspace(0, 1, nq)
        self.q = [np.unique(np.quantile(X[:, j], q)) for j in range(X.shape[1])]
        self.erfinv = erfinv

    def __call__(self, X):
        out = np.empty_like(X, dtype=np.float32)
        for j, qs in enumerate(self.q):
            if len(qs) < 2:
                out[:, j] = 0.0; continue
            u = np.interp(X[:, j], qs, np.linspace(0, 1, len(qs)))
            out[:, j] = np.sqrt(2) * self.erfinv(np.clip(2 * u - 1, -0.999, 0.999))
        return out


def class_w(y):
    w = np.ones(len(y)); w[y] = 0.5 * len(y) / y.sum(); w[~y] = 0.5 * len(y) / (~y).sum(); return w.astype(np.float32)


def train_nn(kind, Xtr, ytr, val, seed, quick):
    """Zwraca funkcję predict(X_po_prep) → p. torch; bez torch tylko w trybie --allow-sklearn (test lokalny)."""
    if kind == "nn64":
        mean, std = Xtr.mean(0), Xtr.std(0); std = np.where(std > 0, std, 1.0)
        tf = lambda X: ((X - mean) / std).astype(np.float32)  # noqa: E731
        hp = dict(h=64, nl=3, act="relu", drop=0.0, lr=0.01, wd=0.0, ep=20, bs=10000, sched="plateau")
    else:
        tf = RankGauss(Xtr)
        hp = dict(h=256, nl=4, act="silu", drop=0.1, lr=2e-3, wd=1e-4, ep=60, bs=2048, sched="cos")
    if quick:
        hp["ep"] = 3
    Z = tf(Xtr); w = class_w(ytr)
    try:
        import torch
        from torch import nn
    except ImportError:
        if not ALLOW_SK:
            raise
        from sklearn.neural_network import MLPClassifier
        clf = MLPClassifier(hidden_layer_sizes=(hp["h"],) * hp["nl"], max_iter=hp["ep"], batch_size=min(hp["bs"], 2048), random_state=seed).fit(Z[~val], ytr[~val])
        return (lambda X: clf.predict_proba(tf(X))[:, 1]), {"backend": "sklearn(TEST)", "hp": hp}
    import copy
    torch.manual_seed(seed)
    act = nn.ReLU if hp["act"] == "relu" else nn.SiLU
    layers, d = [], Z.shape[1]
    for _ in range(hp["nl"]):
        layers += [nn.Linear(d, hp["h"]), act()] + ([nn.Dropout(hp["drop"])] if hp["drop"] else []); d = hp["h"]
    net = nn.Sequential(*layers, nn.Linear(d, 1))
    opt = torch.optim.AdamW(net.parameters(), lr=hp["lr"], weight_decay=hp["wd"])
    Zt, yt, wt = torch.tensor(Z), torch.tensor(ytr, dtype=torch.float32).unsqueeze(1), torch.tensor(w).unsqueeze(1)
    tr, va = torch.tensor(np.nonzero(~val)[0]), torch.tensor(np.nonzero(val)[0])
    steps = hp["ep"] * ((len(tr) + hp["bs"] - 1) // hp["bs"])
    sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, steps) if hp["sched"] == "cos" else torch.optim.lr_scheduler.ReduceLROnPlateau(opt, patience=2, factor=0.5)
    lossf = nn.BCEWithLogitsLoss(reduction="none")
    best, best_w, best_ep, hist = 1e9, None, -1, []
    for ep in range(hp["ep"]):
        net.train(); order = tr[torch.randperm(len(tr))]; tl = []
        for i in range(0, len(order), hp["bs"]):
            idx = order[i:i + hp["bs"]]
            ls = (lossf(net(Zt[idx]), yt[idx]) * wt[idx]).mean()
            opt.zero_grad(set_to_none=True); ls.backward(); opt.step(); tl.append(ls.item())
            if hp["sched"] == "cos":
                sch.step()
        net.eval()
        with torch.no_grad():
            vl = float((lossf(net(Zt[va]), yt[va]) * wt[va]).mean())
        if hp["sched"] != "cos":
            sch.step(vl)
        hist.append([float(np.mean(tl)), vl])
        if vl < best:
            best, best_w, best_ep = vl, copy.deepcopy(net.state_dict()), ep
    net.load_state_dict(best_w); net.eval()

    def predict(X, chunk=500000):
        out = []
        with torch.no_grad():
            for i in range(0, len(X), chunk):
                out.append(torch.sigmoid(net(torch.tensor(tf(X[i:i + chunk])))).numpy()[:, 0])
        return np.concatenate(out) if out else np.zeros(0, np.float32)
    return predict, {"backend": "torch", "hp": hp, "best_epoch": best_ep, "val": best, "hist": hist}


def train_xgb(Xtr, ytr, val, seed, quick, threads):
    import xgboost as xgb
    clf = xgb.XGBClassifier(n_estimators=60 if quick else 1500, max_depth=6, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, min_child_weight=5,
                            tree_method="hist", n_jobs=threads, random_state=seed, early_stopping_rounds=50, eval_metric="logloss")
    w = class_w(ytr)
    clf.fit(Xtr[~val], ytr[~val], sample_weight=w[~val], eval_set=[(Xtr[val], ytr[val])], sample_weight_eval_set=[w[val]], verbose=False)
    return (lambda X: clf.predict_proba(X)[:, 1]), {"best_iteration": int(getattr(clf, "best_iteration", -1))}


def auc(s, y, w=None):
    o = np.argsort(s, kind="mergesort"); s, y = s[o], y[o]; w = np.ones(len(s)) if w is None else w[o]
    r = np.cumsum(w) - 0.5 * w                      # rangi ważone (remisy pomijalne dla score ciągłych)
    ws, wb = w[y].sum(), w[~y].sum()
    return float(((r[y] * w[y]).sum() - 0.5 * ws * ws) / (ws * wb)) if ws > 0 and wb > 0 else float("nan")


def wq(x, w, q):
    o = np.argsort(x); c = np.cumsum(w[o]); return x[o][np.searchsorted(c, q * c[-1])]


def _fit(job):
    key, dp, df_, shape = job
    return key, C.fit_eps(dp, df_, shape, interval=True)


def main():
    global ALLOW_SK
    ap = argparse.ArgumentParser()
    ap.add_argument("--fd", required=True); ap.add_argument("--mb", required=True); ap.add_argument("--data", action="append", default=[], help="GLOB@RUN")
    ap.add_argument("--out", required=True); ap.add_argument("--threads", type=int, default=8); ap.add_argument("--nboot", type=int, default=200)
    ap.add_argument("--quick", action="store_true"); ap.add_argument("--max-rows", type=int, default=0); ap.add_argument("--allow-sklearn", action="store_true")
    ap.add_argument("--no-fit", action="store_true")
    a = ap.parse_args(); ALLOW_SK = a.allow_sklearn
    os.makedirs(a.out, exist_ok=True)
    try:
        import torch; torch.set_num_threads(a.threads)
    except ImportError:
        pass
    fd, mb = load(sorted(glob.glob(a.fd)), a.max_rows), load(sorted(glob.glob(a.mb)), a.max_rows)
    M = {k: np.concatenate([fd[k], mb[k]]) for k in fd if k in mb}
    M["is_fd"] = np.concatenate([np.ones(len(fd["pt"]), bool), np.zeros(len(mb["pt"]), bool)])
    y = np.isin(M["truth"], [0, 1]); bkg = (M["truth"] >= 8) & (M["truth"] <= 15); use = y | bkg
    Xall, med = prep(M["X"])
    log(f"MC: {len(y)} kandydatów, sygnał {y.sum()}, tło {bkg.sum()}, wyłączone {(~use).sum()}; foldy {np.bincount(M['fold'].astype(int))}")
    DATA = {}
    for spec in a.data:
        g, run = spec.rsplit("@", 1); fl = sorted(glob.glob(g))
        if not fl:
            raise SystemExit(f"STOP: brak plików {g}")
        D = load(fl, a.max_rows); D["X"], _ = prep(D["X"], med); DATA[int(run)] = D
        log(f"dane {run}: {len(D['pt'])} kandydatów z {len(fl)} plików")
    S = {m: np.full(len(y), np.nan, np.float32) for m in ("xgb", "nn64", "nn256")}
    SD = {r: {m: np.full(len(D["pt"]), np.nan, np.float32) for m in S} for r, D in DATA.items()}
    info = {}
    rng = np.random.default_rng(20261117)
    hv = rng.random(len(y)) < 0.05                              # 5 % walidacji (wczesne zatrzymanie), wewnątrz próbki treningowej foldu
    for k in sorted(set(M["fold"].astype(int))):
        tr = use & (M["fold"] != k); te = M["fold"] == k
        Xtr, ytr, val = Xall[tr], y[tr], hv[tr]
        for m in S:
            t = time.time()
            if m == "xgb":
                pred, inf = train_xgb(Xtr, ytr, val, 100 + k, a.quick, a.threads)
            else:
                pred, inf = train_nn(m, Xtr, ytr, val, 100 + k, a.quick)
            S[m][te] = pred(Xall[te])
            for r, D in DATA.items():
                dk = D["fold"] == k
                SD[r][m][dk] = pred(D["X"][dk])
            inf.pop("hist", None); info[f"{m}_fold{k}"] = inf
            log(f"fold {k} {m}: {time.time() - t:.0f}s {json.dumps({q: v for q, v in inf.items() if q != 'hp'}, default=float)}")
    S["H"] = M["score"].astype(np.float32); S["avg"] = 0.5 * (S["xgb"] + S["nn256"])
    for r, D in DATA.items():
        SD[r]["H"] = D["score"].astype(np.float32); SD[r]["avg"] = 0.5 * (SD[r]["xgb"] + SD[r]["nn256"])
        np.savez_compressed(os.path.join(a.out, f"nnh_scores_data_{r}.npz"), **{k: D[k] for k in ("run", "event", "row", "file_idx") if k in D}, **{"s_" + m: v for m, v in SD[r].items() if m != "H"})
    np.savez_compressed(os.path.join(a.out, "nnh_scores_mc.npz"), is_fd=M["is_fd"], **{k: M[k] for k in ("event", "row", "file_idx", "truth", "fold")}, **{"s_" + m: v for m, v in S.items() if m != "H"})
    # ------------------------------------------------------------------ metryki MC
    win = np.abs(M["m"] - C.M0_D0) < C.NOM_WIN
    R = {"mc": {}, "data": {}, "sf": [], "info": info, "n_mc": int(len(y)), "n_sig": int(y.sum()), "n_bkg": int(bkg.sum()), "feats": FEATS}
    THR = {}
    for b, (lo, hi) in BINS.items():
        inb = win & (M["pt"] >= lo) & (M["pt"] < hi)
        sg, bg = inb & y & M["is_fd"], inb & bkg & ~M["is_fd"]
        ev = np.nonzero(sg | bg)[0]; yy = y[ev]
        W = rng.poisson(1.0, size=(a.nboot, len(ev))).astype(np.float32)
        R["mc"][b] = {"n_sig": int(sg.sum()), "n_bkg": int(bg.sum())}
        base = {}
        for m in MODELS:
            s = S[m][ev]; rec = {"auc": auc(s, yy)}
            aucb = np.array([auc(s, yy, w) for w in W[:50]]); rec["auc_err"] = float(aucb.std())
            eb, ebb = {}, {}
            for wp in WPS:
                t = np.quantile(s[yy], 1 - wp); THR[(m, b, wp)] = float(t)
                eb[wp] = float((s[~yy] >= t).mean())
                tb = np.array([wq(s[yy], w[yy], 1 - wp) for w in W])                       # próg i ε_B w każdej replice
                ebb[wp] = np.array([(w[~yy] * (s[~yy] >= tt)).sum() / w[~yy].sum() for w, tt in zip(W, tb)])
            rec["eps_b"] = {str(wp): eb[wp] for wp in WPS}; base[m] = ebb
            mm = M["m"][inb & bkg]; ss = S[m][inb & bkg]
            rec["corr_m_bkg"] = float(np.corrcoef(ss, mm)[0, 1]) if len(mm) > 10 else float("nan")
            R["mc"][b][m] = rec
        for m in MODELS:
            zr = {}
            for wp in WPS:
                with np.errstate(all="ignore"):
                    rb = np.sqrt(base["H"][wp] / base[m][wp]); rb = rb[np.isfinite(rb)]
                    e0 = R["mc"][b][m]["eps_b"][str(wp)]
                    zr[str(wp)] = [float(np.sqrt(R["mc"][b]["H"]["eps_b"][str(wp)] / e0)) if e0 > 0 else float("nan"), float(rb.std()) if len(rb) else float("nan")]
            R["mc"][b][m]["Z_over_H"] = zr
    # ------------------------------------------------------------------ dane: ε_B z pasma bocznego Δm oraz SF z tagu D*
    jobs = []
    if DATA:
        from charmtag import dstar_refit as DR
        fsig = np.isin(fd["truth"], [0, 1]) & np.isin(fd["truth_ds"], [0, 1]) & (np.abs(fd["m"] - C.M0_D0) < C.NOM_WIN)
        nfd = len(fd["pt"])
        for r, D in DATA.items():
            bits = {k[2:]: D[k] for k in D if k.startswith("e_tp_") or k == "e_grl"}
            cls = C.classify(bits, r)
            names = ["all", "J12only"] if r != 427929 else ["MB427929"]
            nom = np.abs(D["m"] - C.M0_D0) < C.NOM_WIN
            sb = nom & (D["dm"] > 0.150) & (D["dm"] < 0.165)
            for c in names:
                for b, (lo, hi) in BINS.items():
                    inb = cls[c] & (D["pt"] >= lo) & (D["pt"] < hi)
                    key = f"{r}|{c}|{b}"; R["data"][key] = {"n_sb": int((inb & sb).sum())}
                    for m in MODELS:
                        s = SD[r][m]
                        R["data"][key][m] = {str(wp): float((s[inb & sb] >= THR[(m, b, wp)]).mean()) if (inb & sb).any() else float("nan") for wp in WPS}
                    if a.no_fit:
                        continue
                    shape = DR.shape_from_mc(fd["dm"][fsig & (fd["pt"] >= lo) & (fd["pt"] < hi)] * 1e3)
                    for m in ("H", "xgb", "nn256"):
                        for wp in SFWP:
                            # ε_MC: sygnał fd z prawdą D0 i D* (jak w B_trig), ten sam próg
                            sm = S[m][:nfd][fsig & (fd["pt"] >= lo) & (fd["pt"] < hi)]
                            emc = float((sm >= THR[(m, b, wp)]).mean())
                            pas = s_ok = SD[r][m] >= THR[(m, b, wp)]
                            jobs.append(((r, c, b, m, wp, emc, len(sm)), D["dm"][inb & nom & pas], D["dm"][inb & nom & ~s_ok], shape))
        if jobs:
            log(f"fity pass/fail: {len(jobs)}")
            with Pool(a.threads) as P:
                for (r, c, b, m, wp, emc, nmc), f in P.map(_fit, jobs, chunksize=1):
                    err = C.eps_sym_err(f); sf = f["eps"] / emc if emc > 0 and np.isfinite(f["eps"]) else float("nan")
                    R["sf"].append({"run": r, "cls": c, "bin": b, "model": m, "wp": wp, "eps_data": f["eps"], "eps_data_err": err, "eps_mc": emc, "n_mc": nmc,
                                    "determined": f["determined"], "valid": f["valid"], "N_sig": f["N"], "SF": sf,
                                    "SF_err": float(np.hypot(err / emc, sf * np.sqrt(emc * (1 - emc) / max(nmc, 1)) / emc)) if emc > 0 and np.isfinite(sf) else float("nan")})
    json.dump(R, open(os.path.join(a.out, "nnh.json"), "w"), indent=1, default=float)
    # ------------------------------------------------------------------ raport
    L = ["# TRIG_A/nnh — sieć neuronowa vs drzewa na 18 cechach H (D⁰ z par D⁰←D*, OOF 5 foldów)", "",
         f"MC: {len(y)} kandydatów (sygnał {y.sum()}, tło {bkg.sum()}). Te same kandydaty, cechy i foldy dla xgb / nn64 / nn256; `H` = score nominalny (trening na pełnej próbce inkluzywnej).",
         "Z/Z_H = √(ε_B^H / ε_B) przy tym samym ε_S; błąd z bootstrapu sparowanego. Okno |m − m₀| < 55 MeV; sygnał fd_nd, tło mb_nd.", ""]
    for b in BINS:
        L += [f"## MC, pT {b} GeV (N_S {R['mc'][b]['n_sig']}, N_B {R['mc'][b]['n_bkg']})", "",
              "| model | AUC | " + " | ".join(f"ε_B@{int(wp * 100)}" for wp in WPS) + " | " + " | ".join(f"Z/Z_H@{int(wp * 100)}" for wp in WPS) + " | corr(s, m \\| tło) |",
              "|---|---|" + "---|" * (2 * len(WPS) + 1)]
        for m in MODELS:
            r = R["mc"][b][m]
            L.append(f"| {m} | {r['auc']:.4f} ± {r['auc_err']:.4f} | " + " | ".join(f"{r['eps_b'][str(wp)]:.4f}" for wp in WPS) + " | "
                     + " | ".join(f"{r['Z_over_H'][str(wp)][0]:.3f} ± {r['Z_over_H'][str(wp)][1]:.3f}" for wp in WPS) + f" | {r['corr_m_bkg']:+.3f} |")
        L.append("")
    if R["data"]:
        L += ["## Dane: ε_B z pasma bocznego Δm (0,150–0,165 GeV, okno masy D⁰), progi z MC", "", "| run \\| klasa \\| bin | N_SB | model | " + " | ".join(f"ε_B@{int(wp * 100)}" for wp in WPS) + " | Z/Z_H@50 | Z/Z_H@90 |", "|---|---|---|" + "---|" * (len(WPS) + 2)]
        for key, d in R["data"].items():
            for m in MODELS:
                z = lambda wp: (np.sqrt(d["H"][str(wp)] / d[m][str(wp)]) if d[m][str(wp)] > 0 else float("nan"))  # noqa: E731
                L.append(f"| {key} | {d['n_sb']} | {m} | " + " | ".join(f"{d[m][str(wp)]:.4f}" for wp in WPS) + f" | {z(0.5):.3f} | {z(0.9):.3f} |")
        L.append("")
    if R["sf"]:
        L += ["## Dane: SF z tagu D* (fit Δm pass/fail)", "", "| run | klasa | bin | WP | model | ε_dane | ε_MC | SF | N_sig | określony |", "|---|---|---|---|---|---|---|---|---|---|"]
        for r in sorted(R["sf"], key=lambda r: (r["run"], r["cls"], r["bin"], r["wp"], r["model"])):
            L.append(f"| {r['run']} | {r['cls']} | {r['bin']} | {int(r['wp'] * 100)} | {r['model']} | {r['eps_data']:.3f} ± {r['eps_data_err']:.3f} | {r['eps_mc']:.3f} | {r['SF']:.3f} ± {r['SF_err']:.3f} | {r['N_sig']:.0f} | {r['determined']} |")
    open(os.path.join(a.out, "NNH.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L[:60]))


ALLOW_SK = False
if __name__ == "__main__":
    main()
