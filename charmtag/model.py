"""k-fold XGBoost (hist) po hash(run, event); wynik out-of-fold; progi WP per bin pT."""
from __future__ import annotations

import json
import logging
import os
import pickle
import time

import numpy as np
import pandas as pd

from .samples import fold_of

log = logging.getLogger("charmtag")


def _backend():
    try:
        import xgboost  # noqa: F401
        return "xgboost"
    except Exception as e:  # macOS bez libomp itp.
        log.warning("xgboost niedostępny (%s) → sklearn HistGradientBoosting", e)
        return "sklearn"


def _make(params: dict, backend: str):
    if backend == "xgboost":
        import xgboost as xgb
        p = dict(params)
        es = p.pop("early_stopping_rounds", None)
        nj = p.pop("n_jobs", 0) or None
        return xgb.XGBClassifier(objective="binary:logistic", eval_metric="logloss", early_stopping_rounds=es,
                                 n_jobs=nj, **p)
    from sklearn.ensemble import HistGradientBoostingClassifier
    return HistGradientBoostingClassifier(max_iter=params.get("n_estimators", 300),
                                          learning_rate=params.get("learning_rate", 0.05),
                                          max_depth=params.get("max_depth", 6), early_stopping=True,
                                          validation_fraction=0.1, n_iter_no_change=params.get("early_stopping_rounds", 30))


def _proba(model, X):
    return model.predict_proba(X)[:, 1].astype(np.float32)


def train(S: pd.DataFrame, B: pd.DataFrame, feats: list[str], cfg: dict, outdir: str) -> dict:
    tr = cfg["training"]
    n_folds = int(tr["n_folds"])
    backend = _backend()
    D = pd.concat([S, B], ignore_index=True)
    D["fold"] = fold_of(D["run"].fillna(0).astype(np.int64), D["event"].astype(np.int64), n_folds)
    X = D[feats].to_numpy(dtype=np.float32)
    y = D["y"].to_numpy()
    w = D["w_train"].to_numpy(dtype=float)
    # równoważenie klas (suma wag S = suma wag B)
    w = np.where(y == 1, w / w[y == 1].sum(), w / w[y == 0].sum()) * len(y) / 2
    oof = np.full(len(D), np.nan, np.float32)
    rng = np.random.default_rng(cfg["io"]["seed"])
    os.makedirs(outdir, exist_ok=True)
    info = []
    for k in range(n_folds):
        t0 = time.time()
        te = D["fold"].to_numpy() == k
        trn = ~te
        va = trn & (rng.random(len(D)) < 0.1)
        fit = trn & ~va
        m = _make(tr["xgb"], backend)
        if backend == "xgboost":
            m.fit(X[fit], y[fit], sample_weight=w[fit], eval_set=[(X[va], y[va])], sample_weight_eval_set=[w[va]],
                  verbose=False)
            m.save_model(os.path.join(outdir, f"fold{k}.json"))
            best = int(getattr(m, "best_iteration", -1) or -1)
        else:
            m.fit(X[trn], y[trn], sample_weight=w[trn])
            with open(os.path.join(outdir, f"fold{k}.pkl"), "wb") as f:
                pickle.dump(m, f)
            best = int(getattr(m, "n_iter_", -1))
        oof[te] = _proba(m, X[te])
        info.append({"fold": k, "n_train": int(fit.sum()), "n_test": int(te.sum()), "best_iteration": best,
                     "sec": round(time.time() - t0, 1)})
        log.info("fold %d: n_train=%d best_it=%d %.0fs", k, fit.sum(), best, time.time() - t0)
    D["score"] = oof
    imp = feature_importance(outdir, feats, n_folds, backend)
    meta = {"features": feats, "n_folds": n_folds, "backend": backend, "folds": info, "importance": imp}
    return D, meta


def feature_importance(outdir, feats, n_folds, backend):
    if backend != "xgboost":
        return {}
    import xgboost as xgb
    acc = {f: 0.0 for f in feats}
    for k in range(n_folds):
        b = xgb.Booster()
        b.load_model(os.path.join(outdir, f"fold{k}.json"))
        g = b.get_score(importance_type="gain")
        for i, f in enumerate(feats):
            acc[f] += g.get(f"f{i}", g.get(f, 0.0)) / n_folds
    tot = sum(acc.values()) or 1
    return {f: round(v / tot, 4) for f, v in sorted(acc.items(), key=lambda x: -x[1])}


class Ensemble:
    """Modele foldów; kandydat dostaje wynik modelu foldu, w którym NIE był trenowany."""

    def __init__(self, mdir: str):
        with open(os.path.join(mdir, "meta.json")) as f:
            self.meta = json.load(f)
        self.feats = self.meta["features"]
        self.n = self.meta["n_folds"]
        self.models = []
        for k in range(self.n):
            if self.meta["backend"] == "xgboost":
                import xgboost as xgb
                m = xgb.XGBClassifier()
                m.load_model(os.path.join(mdir, f"fold{k}.json"))
            else:
                with open(os.path.join(mdir, f"fold{k}.pkl"), "rb") as f:
                    m = pickle.load(f)
            self.models.append(m)

    def score(self, L: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        fold = fold_of(L["run"].fillna(0).astype(np.int64), L["event"].astype(np.int64), self.n)
        X = L[self.feats].to_numpy(dtype=np.float32)
        s = np.full(len(L), np.nan, np.float32)
        for k, m in enumerate(self.models):
            sel = fold == k
            if sel.any():
                s[sel] = _proba(m, X[sel])
        return s, fold


def wp_thresholds(score, w, pt, pt_bins, wps) -> dict:
    """Próg cięcia dający ε_S = wp w każdym binie pT (ważone)."""
    out = {}
    for lo, hi in zip(pt_bins[:-1], pt_bins[1:]):
        sel = (pt >= lo) & (pt < hi) & np.isfinite(score)
        s, ww = score[sel], w[sel]
        o = np.argsort(s)[::-1]
        c = np.cumsum(ww[o]) / max(ww.sum(), 1e-12)
        out[f"{lo}-{hi}"] = {f"{wp:.2f}": float(s[o][min(np.searchsorted(c, wp), len(o) - 1)]) if len(o) else 1.0
                             for wp in wps}
    return out


def pass_wp(score, pt, thr: dict, wp: float) -> np.ndarray:
    ok = np.zeros(len(score), bool)
    for key, d in thr.items():
        lo, hi = map(float, key.split("-"))
        sel = (pt >= lo) & (pt < hi)
        ok |= sel & (score >= d[f"{wp:.2f}"])
    return ok
