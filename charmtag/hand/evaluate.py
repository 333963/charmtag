"""Metryki: AUC, ε_B @ ε_S w binach pT, test rzeźbienia masy (R_win, kształt m po cięciu)."""
from __future__ import annotations
import numpy as np, pandas as pd

from .. import constants as _C
from sklearn.metrics import roc_auc_score

PT_BINS = [(3, 5), (5, 8), (8, 30)]
EFFS = (0.8, 0.7, 0.5)
M_NOM = {"D0": _C.M_D0, "Dp": _C.M_DPLUS}

def eps_b_at(y, s, w=None, effs=EFFS):
    y = np.asarray(y).astype(bool); s = np.asarray(s, float)
    w = np.ones(len(y)) if w is None else np.asarray(w, float)
    out = {}
    for e in effs:
        # próg = kwantyl (ważony) scoru sygnału przy zadanej ε_S
        ss, ws = s[y], w[y]
        o = np.argsort(ss); cw = np.cumsum(ws[o]) / ws.sum()
        thr = ss[o][np.searchsorted(cw, 1 - e)]
        out[e] = float(w[~y][s[~y] >= thr].sum() / max(w[~y].sum(), 1e-12))
        out[f"thr{e}"] = float(thr)
    return out

def per_ptbin(df, y, s, pt_col="pt", w=None, bins=PT_BINS):
    rows = []
    pt = df[pt_col].to_numpy(float)
    w = np.ones(len(df)) if w is None else np.asarray(w, float)
    for lo, hi in [("all", None)] + bins:
        sel = np.ones(len(df), bool) if lo == "all" else (pt >= lo) & (pt < hi)
        yy, ss, ww = np.asarray(y)[sel], np.asarray(s)[sel], w[sel]
        if yy.sum() < 20 or (~yy.astype(bool)).sum() < 20:
            continue
        eb = eps_b_at(yy, ss, ww)
        rows.append({"pt_bin": "all" if lo == "all" else f"{lo}-{hi}", "n_sig": int(yy.sum()), "n_bkg": int(len(yy) - yy.sum()),
                     "AUC": roc_auc_score(yy, ss, sample_weight=ww),
                     **{f"epsB@{int(e*100)}": eb[e] for e in EFFS}})
    return pd.DataFrame(rows)

def sculpting(m, s, thr, species="D0", win=0.025, sb=(0.05, 0.12), nbins=40, w=None):
    """R_win = (N_win/N_SB)_po / (N_win/N_SB)_przed  dla próbki TŁA (MC tła lub SB danych).
    1.0 = brak rzeźbienia; 1.24–1.3 = poziom XGB D+ z docu S. Zwraca też znormalizowane histogramy m."""
    m = np.asarray(m, float); s = np.asarray(s, float)
    w = np.ones(len(m)) if w is None else np.asarray(w, float)
    m0 = M_NOM[species]
    inwin = np.abs(m - m0) < win
    insb = (np.abs(m - m0) > sb[0]) & (np.abs(m - m0) < sb[1])
    cut = s >= thr
    def ratio(sel):
        a, b = w[sel & inwin].sum(), w[sel & insb].sum()
        return a / b if b > 0 else np.nan
    r_before, r_after = ratio(np.ones_like(cut)), ratio(cut)
    nwa, nsa = (cut & inwin).sum(), (cut & insb).sum()
    r_err = (r_after / r_before) * np.sqrt(1 / max(nwa, 1) + 1 / max(nsa, 1)) if (r_before and nwa and nsa) else np.nan
    edges = np.linspace(m0 - 0.2, m0 + 0.2, nbins + 1)
    h0, _ = np.histogram(m, edges, weights=w); h1, _ = np.histogram(m[cut], edges, weights=w[cut])
    h0n = h0 / max(h0.sum(), 1); h1n = h1 / max(h1.sum(), 1)
    # odległość kształtów: max |CDF różnica| (KS-like) + χ² kształtu
    ks = float(np.max(np.abs(np.cumsum(h0n) - np.cumsum(h1n))))
    return {"R_win": float(r_after / r_before) if r_before else np.nan, "R_win_err": float(r_err), "n_win_after": int(nwa), "n_sb_after": int(nsa), "eff_cut": float(w[cut].sum() / w.sum()),
            "KS_shape": ks, "edges": edges, "h_before": h0n, "h_after": h1n}

def sculpting_table(m, s, y, species, w=None, effs=EFFS):
    """R_win tła przy progach zdefiniowanych przez ε_S na sygnale (jeśli y podane) albo kwantylach scoru."""
    rows = []
    m = np.asarray(m); s = np.asarray(s)
    if y is not None and np.asarray(y).sum() > 20:
        eb = eps_b_at(y, s, w, effs); thrs = {e: eb[f"thr{e}"] for e in effs}
        bkg = ~np.asarray(y).astype(bool)
    else:
        thrs = {e: float(np.quantile(s, 1 - e)) for e in effs}; bkg = np.ones(len(s), bool)
    for e, t in thrs.items():
        r = sculpting(m[bkg], s[bkg], t, species, w=None if w is None else np.asarray(w)[bkg])
        rows.append({"epsS": e, "thr": t, "R_win": r["R_win"], "R_win_err": r["R_win_err"], "n_win": r["n_win_after"], "n_sb": r["n_sb_after"], "eff_bkg": r["eff_cut"], "KS_shape": r["KS_shape"]})
    return pd.DataFrame(rows)
