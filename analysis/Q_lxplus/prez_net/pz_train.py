#!/usr/bin/env python3
"""PREZ-2 / pz_train (19.11.2026, projekt A) — trening klasyfikatorów pojedynczych D: rodzina (xgb | net) × wejścia (hi | lo).

  hi = 18 cech H (jak nominalny XGB H: log pT, |η|, χ², L_xy, istotności parametru zderzenia, pointing, |cosθ*|, d0/z0/η córek, dE/dx i p K/π…)
  lo = ślady córek (log pT, η, Δφ i Δη względem D, log p, dE/dx + flaga braku, N pkt, d0, z0, rola K/π) + L_xy kandydata
  xgb: hiperparametry nominału (1500 drzew, głębokość 6, lr 0,05, subsample/colsample 0,8, min_child_weight 5, early stopping 50);
       lo → płaski wektor ślady (K, π malejąco po pT) + L_xy, NaN natywnie.
  net: hi → MLP 18→128→128→64→1 (SiLU, dropout 0,1); lo → Deep Sets (koder śladu φ 11→64→64→64, pooling suma ⊕ max, ⊕ L_xy → ρ 129→128→64→1).
       Wejścia: transformacja rang → Gauss dopasowana na foldach treningowych, braki → mediana (+ flaga). BCE z wagami, AdamW, OneCycle, ≤ 30 epok,
       early stopping (cierpliwość 5) na stracie walidacyjnej.
Próbka: okno sygnału |m − m0| < signal_window (jak region sygnału nominału) — sygnał fd truth, tło mb kombinatoryka; poza oknem: tło do testu rzeźbienia
(model foldu kandydata, kandydaci nigdy w treningu). Foldy: charmtag fold_of(run, event), 5 (= nominał), walidacja 10 % z foldów treningowych.
Wagi: tło przeważone do (pT, |η|) sygnału (charmtag.samples.reweight, biny nominału), klasy zrównoważone (jak model.train).
Ważność (hi): permutacyjna ΔAUC na foldach testowych (oba modele) + SHAP (xgb: pred_contribs).
Wyjście (--out): scores_<fam>.npz (s_hi, s_lo na wszystkich wierszach sig i bkg; OOF w oknie), imp_<fam>.json, info_<fam>.json.
Użycie: pz_train.py --sig sig.npz --bkg bkg.npz --family xgb|net --out DIR [--threads 8] [--epochs 30] [--quick]"""
import argparse
import json
import os
import sys
import time

import numpy as np

T0 = time.time()
log = lambda *a: print(f"[{time.time() - T0:6.0f}s]", *a, flush=True)  # noqa: E731
XGB = dict(n_estimators=1500, learning_rate=0.05, max_depth=6, min_child_weight=5, subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0)


# ----------------------------------------------------------------------------- dane
def load(sig, bkg):
    S, B = np.load(sig, allow_pickle=False), np.load(bkg, allow_pickle=False)
    ms, mb = json.loads(str(S["meta"])), json.loads(str(B["meta"]))
    D = {}
    for k in ("Xh", "T", "G", "m", "pt", "abs_eta", "weight", "fold", "inwin", "s_nom", "run", "event"):
        D[k] = np.concatenate([S[k], B[k]])
    D["y"] = np.concatenate([np.ones(len(S["m"]), np.int8), np.zeros(len(B["m"]), np.int8)])
    return D, ms, mb


def train_weights(D, bins):
    """w: tło przeważone do (pT, |η|) sygnału w oknie; klasy zrównoważone (suma S = suma B = n/2). Tylko wiersze w oknie (reszta 0)."""
    w = np.zeros(len(D["y"]))
    win = D["inwin"]
    s, b = win & (D["y"] == 1), win & (D["y"] == 0)
    bp, be = np.asarray(bins["pt"]), np.asarray(bins["abs_eta"])
    hs, _, _ = np.histogram2d(D["pt"][s], D["abs_eta"][s], [bp, be], weights=D["weight"][s])
    hb, _, _ = np.histogram2d(D["pt"][b], D["abs_eta"][b], [bp, be], weights=D["weight"][b])
    hs /= max(hs.sum(), 1e-12); hb /= max(hb.sum(), 1e-12)
    ratio = np.divide(hs, hb, out=np.zeros_like(hs), where=hb > 0)
    ip = np.clip(np.digitize(D["pt"], bp) - 1, 0, len(bp) - 2); ie = np.clip(np.digitize(D["abs_eta"], be) - 1, 0, len(be) - 2)
    w[s] = D["weight"][s]
    w[b] = D["weight"][b] * ratio[ip[b], ie[b]]
    n = win.sum()
    w[s] *= 0.5 * n / max(w[s].sum(), 1e-12); w[b] *= 0.5 * n / max(w[b].sum(), 1e-12)
    return w


def flat_lo(D):
    T = D["T"].reshape(len(D["T"]), -1)
    return np.concatenate([T, D["G"]], 1).astype(np.float32)


def wauc(s, y, w):
    o = np.argsort(-s, kind="stable"); s, y, w = s[o], y[o], w[o]
    tp = np.cumsum(w * (y == 1)); fp = np.cumsum(w * (y == 0))
    if tp[-1] <= 0 or fp[-1] <= 0:
        return np.nan
    tpr, fpr = np.r_[0, tp / tp[-1]], np.r_[0, fp / fp[-1]]
    return float(np.trapezoid(tpr, fpr) if hasattr(np, "trapezoid") else np.trapz(tpr, fpr))


# ----------------------------------------------------------------------------- preprocessing (rank → Gauss)
class RankGauss:
    def __init__(self, X, nq=400):
        from scipy.special import erfinv
        self.erfinv = erfinv
        X = X.reshape(-1, X.shape[-1])
        self.med = np.nan_to_num(np.nanmedian(X, 0)).astype(np.float32)
        Xf = np.where(np.isfinite(X), X, self.med)
        self.q = [np.unique(np.quantile(Xf[:, j], np.linspace(0, 1, nq))) for j in range(X.shape[1])]

    def __call__(self, X):
        sh = X.shape; X = X.reshape(-1, sh[-1])
        X = np.where(np.isfinite(X), X, self.med)
        Y = np.zeros(X.shape, np.float32)
        for j, qs in enumerate(self.q):
            if len(qs) < 2:
                continue
            u = np.interp(X[:, j], qs, np.linspace(0, 1, len(qs)))
            Y[:, j] = np.sqrt(2) * self.erfinv(np.clip(2 * u - 1, -0.998, 0.998))
        return Y.reshape(sh)


def nan_flags(X):
    """Kolumny z brakami → flagi 0/1 (dla sieci; xgb ma NaN natywnie)."""
    cols = np.flatnonzero(~np.isfinite(X).all(0))
    return cols


# ----------------------------------------------------------------------------- sieci
def make_mlp(nin):
    import torch.nn as nn
    return nn.Sequential(nn.Linear(nin, 128), nn.SiLU(), nn.Dropout(0.1), nn.Linear(128, 128), nn.SiLU(), nn.Dropout(0.1),
                         nn.Linear(128, 64), nn.SiLU(), nn.Linear(64, 1))


def make_deepsets(ntf, nglob):
    import torch
    import torch.nn as nn

    class DeepSets(nn.Module):
        """Koder śladu φ (wspólny dla wszystkich śladów) → suma ⊕ max po śladach (niezmiennicze względem permutacji) ⊕ cechy globalne → ρ."""
        def __init__(self):
            super().__init__()
            self.phi = nn.Sequential(nn.Linear(ntf, 64), nn.SiLU(), nn.Linear(64, 64), nn.SiLU(), nn.Linear(64, 64), nn.SiLU())
            self.rho = nn.Sequential(nn.Linear(128 + nglob, 128), nn.SiLU(), nn.Dropout(0.1), nn.Linear(128, 64), nn.SiLU(), nn.Linear(64, 1))

        def forward(self, t, g):
            h = self.phi(t)
            return self.rho(torch.cat([h.sum(1), h.max(1).values, g], 1))

    return DeepSets()


def nparams(net):
    return int(sum(p.numel() for p in net.parameters()))


def train_net(kind, Xtr, ytr, wtr, Xva, yva, wva, epochs, seed, batch=None):
    """kind: 'mlp' (X = [n, f]) | 'ds' (X = (T [n, nd, f], G [n, g])). Zwraca (predict(X) → p, info)."""
    import torch
    torch.manual_seed(seed); np.random.seed(seed)
    if kind == "mlp":
        net = make_mlp(Xtr.shape[1]); tens = lambda X, i: (torch.from_numpy(X[i]),)  # noqa: E731
    else:
        net = make_deepsets(Xtr[0].shape[2], Xtr[1].shape[1]); tens = lambda X, i: (torch.from_numpy(X[0][i]), torch.from_numpy(X[1][i]))  # noqa: E731
    fwd = lambda X, i: net(*tens(X, i)).squeeze(1)  # noqa: E731
    ntr = len(ytr)
    batch = batch or int(min(4096, max(256, ntr // 50)))                         # ≥ 50 kroków/epokę (małe próbki); 4096 przy pełnej statystyce
    yt, wt = torch.from_numpy(ytr.astype(np.float32)), torch.from_numpy(wtr.astype(np.float32))
    yv, wv = torch.from_numpy(yva.astype(np.float32)), torch.from_numpy(wva.astype(np.float32))
    opt = torch.optim.AdamW(net.parameters(), lr=2e-3, weight_decay=1e-4)
    spe = max(1, ntr // batch)
    sch = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=2e-3, total_steps=epochs * spe + 1)
    lossf = torch.nn.BCEWithLogitsLoss(reduction="none")
    best, best_state, bad, hist = np.inf, None, 0, []
    for ep in range(epochs):
        net.train()
        perm = np.random.permutation(ntr)
        for i in range(spe):
            idx = perm[i * batch:(i + 1) * batch] if ntr >= batch else perm
            opt.zero_grad()
            loss = (lossf(fwd(Xtr, idx), yt[idx]) * wt[idx]).mean()
            loss.backward(); opt.step(); sch.step()
        net.eval()
        with torch.no_grad():
            nv = len(yva); vl = 0.0
            for i in range(0, nv, 65536):
                idx = np.arange(i, min(nv, i + 65536))
                vl += float((lossf(fwd(Xva, idx), yv[idx]) * wv[idx]).sum())
            vl /= max(float(wv.sum()), 1e-12)
        hist.append(vl)
        if vl < best - 1e-5:
            best, best_state, bad = vl, {k: v.clone() for k, v in net.state_dict().items()}, 0
        else:
            bad += 1
            if bad >= 5:
                break
    net.load_state_dict(best_state); net.eval()

    def predict(X):
        n = len(X) if kind == "mlp" else len(X[0])
        out = np.empty(n, np.float32)
        with torch.no_grad():
            for i in range(0, n, 65536):
                idx = np.arange(i, min(n, i + 65536))
                out[idx] = torch.sigmoid(fwd(X, idx)).numpy()
        return out

    return predict, {"batch": batch, "epochs_run": len(hist), "best_val_loss": float(best), "val_hist": hist, "n_params": nparams(net)}


# ----------------------------------------------------------------------------- xgb
def train_xgb(Xtr, ytr, wtr, Xva, yva, wva, threads, quick):
    import xgboost as xgb
    kw = dict(XGB); kw["n_estimators"] = 60 if quick else kw["n_estimators"]
    m = xgb.XGBClassifier(objective="binary:logistic", eval_metric="logloss", tree_method="hist", early_stopping_rounds=50, n_jobs=threads, **kw)
    m.fit(Xtr, ytr, sample_weight=wtr, eval_set=[(Xva, yva)], sample_weight_eval_set=[wva], verbose=False)
    bi = int(getattr(m, "best_iteration", -1) or -1)
    return (lambda X: m.predict_proba(X)[:, 1].astype(np.float32)), {"best_iteration": bi, "n_trees": bi + 1}, m


# ----------------------------------------------------------------------------- główna pętla
def run(D, w, family, inputs, nfold, threads, epochs, quick, imp_names, perm_cap=150_000):
    """OOF w oknie (y), predykcja modelem foldu dla wierszy spoza okna. Zwraca (score[n], info, importance)."""
    rng = np.random.default_rng(17)
    n = len(D["y"]); score = np.full(n, np.nan, np.float32)
    win = D["inwin"]; fold = D["fold"].astype(int); y = D["y"]
    if inputs == "hi":
        Xraw = D["Xh"].astype(np.float32)
    elif family == "xgb":
        Xraw = flat_lo(D)
    else:
        Xraw = None
    info, imp_acc, shap_acc = [], [], []
    for k in range(nfold):
        te, trn = fold == k, (fold != k) & win
        va = trn & (rng.random(n) < 0.1); fit = trn & ~va
        tew = te & win
        if family == "xgb":
            pred, inf, model = train_xgb(Xraw[fit], y[fit], w[fit], Xraw[va], y[va], w[va], threads, quick)
            score[te] = pred(Xraw[te])
            if inputs == "hi":
                import xgboost as xgb
                sub = np.flatnonzero(tew); sub = sub if len(sub) <= 100_000 else rng.choice(sub, 100_000, replace=False)
                c = model.get_booster().predict(xgb.DMatrix(Xraw[sub]), pred_contribs=True)[:, :-1]
                shap_acc.append(np.abs(c).mean(0))
        else:
            if inputs == "hi":
                bad = nan_flags(Xraw[fit])
                rg = RankGauss(Xraw[fit])
                prep = lambda X: np.concatenate([rg(X), (~np.isfinite(X[:, bad])).astype(np.float32)], 1)  # noqa: E731
                pred_, inf = train_net("mlp", prep(Xraw[fit]), y[fit], w[fit], prep(Xraw[va]), y[va], w[va], epochs if not quick else 3, 100 + k)
                pred = lambda X: pred_(prep(X))  # noqa: E731
                score[te] = pred(Xraw[te])
            else:
                T, G = D["T"], D["G"]
                rgT = RankGauss(T[fit]); rgG = RankGauss(G[fit])
                prepT = lambda Tx: rgT(Tx)  # noqa: E731
                X_fit = (prepT(T[fit]), rgG(G[fit])); X_va = (prepT(T[va]), rgG(G[va]))
                pred_, inf = train_net("ds", X_fit, y[fit], w[fit], X_va, y[va], w[va], epochs if not quick else 3, 200 + k)
                te_i = np.flatnonzero(te)
                for i in range(0, len(te_i), 500_000):
                    j = te_i[i:i + 500_000]
                    score[j] = pred_((prepT(T[j]), rgG(G[j])))
                pred = None
        inf.update({"fold": k, "n_fit": int(fit.sum()), "n_val": int(va.sum()), "n_test_win": int(tew.sum()),
                    "auc_test": wauc(score[tew], y[tew], w[tew])})
        info.append(inf); log(f"{family}/{inputs} fold {k}: AUC {inf['auc_test']:.4f}", {kk: vv for kk, vv in inf.items() if kk in ('n_trees', 'epochs_run', 'n_params')})
        if inputs == "hi" and pred is not None:                                   # ważność permutacyjna (ΔAUC)
            sub = np.flatnonzero(tew); sub = sub if len(sub) <= perm_cap else rng.choice(sub, perm_cap, replace=False)
            Xs = Xraw[sub]; base = wauc(pred(Xs), y[sub], w[sub]); d = []
            for j in range(Xs.shape[1]):
                Xp = Xs.copy(); Xp[:, j] = Xp[rng.permutation(len(sub)), j]
                d.append(base - wauc(pred(Xp), y[sub], w[sub]))
            imp_acc.append(d)
    imp = {}
    if imp_acc:
        a = np.asarray(imp_acc); imp["perm_dauc"] = dict(zip(imp_names, a.mean(0).tolist())); imp["perm_dauc_err"] = dict(zip(imp_names, (a.std(0) / np.sqrt(len(a))).tolist()))
    if shap_acc:
        imp["shap_mean_abs"] = dict(zip(imp_names, np.asarray(shap_acc).mean(0).tolist()))
    return score, info, imp


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sig", required=True); ap.add_argument("--bkg", required=True); ap.add_argument("--family", required=True, choices=["xgb", "net"])
    ap.add_argument("--out", required=True); ap.add_argument("--threads", type=int, default=8); ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--quick", action="store_true"); ap.add_argument("--inputs", default="hi,lo")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    if a.family == "net":
        import torch
        torch.set_num_threads(a.threads)
    D, ms, mb = load(a.sig, a.bkg)
    w = train_weights(D, ms["reweight_bins"])
    nfold = int(ms["n_folds"])
    log(f"{ms['species']}: sygnał {int((D['y'] == 1).sum())} (okno {int(((D['y'] == 1) & D['inwin']).sum())}), tło {int((D['y'] == 0).sum())} "
        f"(okno {int(((D['y'] == 0) & D['inwin']).sum())})")
    res, infos, imps = {}, {}, {}
    for inp in a.inputs.split(","):
        s, inf, imp = run(D, w, a.family, inp, nfold, a.threads, a.epochs, a.quick, ms["features_hi"])
        res[f"s_{inp}"] = s; infos[inp] = inf; imps[inp] = imp
    np.savez_compressed(os.path.join(a.out, f"scores_{a.family}.npz"), y=D["y"], w=w, inwin=D["inwin"], m=D["m"], pt=D["pt"], fold=D["fold"],
                        s_nom=D["s_nom"], **res)
    json.dump({"species": ms["species"], "family": a.family, "folds": infos, "meta_sig": ms, "meta_bkg": mb, "xgb_params": XGB, "epochs_max": a.epochs},
              open(os.path.join(a.out, f"info_{a.family}.json"), "w"), indent=1, default=float)
    json.dump(imps, open(os.path.join(a.out, f"imp_{a.family}.json"), "w"), indent=1, default=float)
    log("OK", a.out)


if __name__ == "__main__":
    sys.exit(main())
