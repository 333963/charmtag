"""Fit masy (iminuit, binned extended NLL): Gauss/2-Gauss + exp/pol2 + refleksja K↔π (D0)."""
from __future__ import annotations

import math

import numpy as np
from iminuit import Minuit
from iminuit.cost import ExtendedBinnedNLL
from scipy.special import erf

from .constants import DM_DSTAR, M_PI


def _gcdf(x, mu, s):
    return 0.5 * (1 + erf((x - mu) / (math.sqrt(2) * s)))


def fit_mc_shape(m, w, lo, hi, m0, model="gauss2", nbins=100):
    """Kształt sygnału z MC truth: mu, s1, (f, r=s2/s1)."""
    m = np.asarray(m, dtype=np.float64)   # float32 → krawędzie histogramu float32 → szum NLL (NEP 50)
    sel = (m >= lo) & (m < hi)
    m, w = m[sel], np.asarray(w, dtype=np.float64)[sel]
    h, xe = np.histogram(m, nbins, (lo, hi), weights=w)
    sd = float(np.sqrt(np.cov(m, aweights=w))) if len(m) > 2 else 0.02
    sd = min(max(sd, 0.004), 0.05)
    if model == "gauss":
        def cdf(x, n, mu, s1):
            return n * _gcdf(x, mu, s1)
        mi = Minuit(ExtendedBinnedNLL(h, xe, cdf), n=h.sum(), mu=m0, s1=sd)
    else:
        def cdf(x, n, mu, s1, f, r):
            return n * (f * _gcdf(x, mu, s1) + (1 - f) * _gcdf(x, mu, s1 * r))
        mi = Minuit(ExtendedBinnedNLL(h, xe, cdf), n=h.sum(), mu=m0, s1=sd * 0.7, f=0.7, r=2.0)
        mi.limits["f"] = (0.05, 1.0)
        mi.limits["r"] = (1.05, 6.0)
    mi.limits["s1"] = (0.002, 0.08)
    mi.limits["mu"] = (m0 - 0.03, m0 + 0.03)
    mi.limits["n"] = (0, None)
    mi.migrad()
    out = {k: float(mi.values[k]) for k in mi.parameters}
    out.update({f"e_{k}": float(mi.errors[k]) for k in mi.parameters})
    out["valid"] = bool(mi.valid)
    out["model"] = model
    out["sigma_eff"] = _sigma_eff(out)
    return out


def _sigma_eff(sh):
    if sh.get("model") == "gauss2":
        f, s1, r = sh["f"], sh["s1"], sh["r"]
        return float(math.sqrt(f * s1 ** 2 + (1 - f) * (r * s1) ** 2))
    return float(sh["s1"])


def fit_data(m, lo, hi, m0, shape: dict, cfg_fit: dict, refl_template=None, refl_ratio=0.0, w=None):
    """Fit danych. shape = wynik fit_mc_shape. refl_template = histogram (nbins) z MC (kształt)."""
    nb = int(cfg_fit["nbins"])
    m = np.asarray(m, dtype=np.float64)   # float32 → szum NLL → zły MIGRAD/HESSE
    sel = (m >= lo) & (m < hi)
    h, xe = np.histogram(m[sel], nb, (lo, hi), weights=None if w is None else np.asarray(w, dtype=np.float64)[sel])
    span = hi - lo
    bkgm = cfg_fit["background_model"]
    gauss2 = shape.get("model") == "gauss2"
    f2, r2 = (shape["f"], shape["r"]) if gauss2 else (1.0, 1.0)
    if refl_template is not None and refl_ratio > 0:
        rt = np.asarray(refl_template, float)
        rc = np.concatenate([[0.0], np.cumsum(rt) / max(rt.sum(), 1e-12)])
    else:
        rc = None

    def sig_cdf(x, mu, s1):
        return f2 * _gcdf(x, mu, s1) + (1 - f2) * _gcdf(x, mu, s1 * r2)

    def bkg_cdf(x, c1, c2):
        # pol2: 1 + c1 t + c2 t² ≥ 0 na t∈[-1,1] wymuszone limitami c1∈[-0.95,0.95], c2∈[-0.04,0.5] (min ≥ 0,01, sprawdzone numerycznie)
        if bkgm == "exp":
            if abs(c1) < 1e-6:
                return (x - lo) / span
            return (1 - np.exp(-c1 * (x - lo))) / (1 - math.exp(-c1 * span))
        t = 2 * (x - lo) / span - 1  # pol2 na [-1,1]: 1 + c1 t + c2 t^2 (całka)
        F = lambda u: u + c1 * u ** 2 / 2 + c2 * u ** 3 / 3
        return (F(t) - F(-1.0)) / max(F(1.0) - F(-1.0), 1e-12)

    def cdf(x, ns, nb_, mu, s1, c1, c2):
        out = ns * (sig_cdf(x, mu, s1) - sig_cdf(lo, mu, s1)) + nb_ * bkg_cdf(x, c1, c2)
        if rc is not None:
            out = out + refl_ratio * ns * np.interp(x, xe, rc)
        return out

    cost = ExtendedBinnedNLL(h, xe, cdf)
    smc = shape["s1"]
    mode = cfg_fit["sigma_mode"]
    if mode == "constrained":
        from iminuit.cost import NormalConstraint
        cost = cost + NormalConstraint("s1", smc, cfg_fit["sigma_constraint_rel"] * smc)
    ntot = max(h.sum(), 1.0)
    mi = Minuit(cost, ns=0.05 * ntot, nb_=0.95 * ntot, mu=shape["mu"], s1=smc, c1=1.0 if bkgm == "exp" else 0.0, c2=0.0)
    mi.limits["ns"] = (-0.2 * ntot, 2 * ntot)
    mi.limits["nb_"] = (0, 2 * ntot)
    mi.limits["mu"] = (m0 - 0.02, m0 + 0.02)
    mi.limits["s1"] = (0.3 * smc, 3 * smc)
    if bkgm == "exp":
        mi.fixed["c2"] = True
        mi.limits["c1"] = (-50, 50)
    else:   # dodatniość tła pol2 na [-1,1]: prostokąt c1∈[-0.95,0.95], c2∈[-0.04,0.5] ⇒ min(1 + c1 t + c2 t²) ≥ 0,01 (P2-12, recenzja kodu 2.0.1)
        mi.limits["c1"] = (-0.95, 0.95)
        mi.limits["c2"] = (-0.04, 0.5)
    if mode == "fixed":
        mi.fixed["s1"] = True
    mi.migrad()
    mi.hesse()
    v = {k: float(mi.values[k]) for k in mi.parameters}
    e = {k: float(mi.errors[k]) for k in mi.parameters}
    mu_, s1_ = v["mu"], v["s1"]
    seff = _sigma_eff({**shape, "s1": s1_})
    a, b = mu_ - 2 * seff, mu_ + 2 * seff
    nbw = v["nb_"] * (bkg_cdf(b, v["c1"], v["c2"]) - bkg_cdf(a, v["c1"], v["c2"]))
    sw = v["ns"] * (sig_cdf(b, mu_, s1_) - sig_cdf(a, mu_, s1_))
    pred = np.diff(cdf(xe, *[v[k] for k in mi.parameters]))
    good = pred > 0
    chi2 = float(np.sum((h[good] - pred[good]) ** 2 / np.maximum(pred[good], 1)))
    ndf = int(good.sum() - mi.nfit)
    return {
        "N_sig": v["ns"], "N_sig_err": e["ns"], "N_bkg": v["nb_"], "mu": mu_, "mu_err": e["mu"],
        "sigma1": s1_, "sigma_eff": seff, "sigma_mc": smc, "S_2sig": sw, "B_2sig": nbw,
        "signif": sw / math.sqrt(max(sw + nbw, 1e-9)), "chi2": chi2, "ndf": ndf, "valid": bool(mi.valid),
        "refl_ratio": refl_ratio if rc is not None else 0.0, "n_entries": int(h.sum()),
        "_hist": h.tolist(), "_edges": xe.tolist(), "_pred": pred.tolist(),
        "_comp_bkg": (v["nb_"] * np.diff(bkg_cdf(xe, v["c1"], v["c2"]))).tolist(),
        "_comp_refl": (refl_ratio * v["ns"] * np.diff(np.interp(xe, xe, rc))).tolist() if rc is not None else None,
    }


def fit_dm(dm, lo=0.1396, hi=0.165, nbins=64, mu0=DM_DSTAR, s0=0.0007, fix_shape=False):
    """Fit Δm(D*): Gauss + tło progowe (x−x0)^a·exp(−b(x−x0)), x0 = m_π. Zwraca N_D* ± σ."""
    x0 = M_PI
    dm = np.asarray(dm, dtype=np.float64)   # float32 → szum NLL
    h, xe = np.histogram(dm[(dm >= lo) & (dm < hi)], nbins, (lo, hi))
    grid = np.linspace(lo, hi, 2001)

    def bcdf(x, a, b):
        u = np.clip(grid - x0, 1e-7, None)
        dens = u ** a * np.exp(-b * u)
        c = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(grid))])
        return np.interp(x, grid, c / max(c[-1], 1e-300))

    def cdf(x, ns, nb, mu, s, a, b):
        return ns * (_gcdf(x, mu, s) - _gcdf(lo, mu, s)) + nb * bcdf(x, a, b)

    n = max(h.sum(), 1)
    mi = Minuit(ExtendedBinnedNLL(h, xe, cdf), ns=0.1 * n, nb=0.9 * n, mu=mu0, s=s0, a=0.5, b=10.0)
    mi.limits["ns"] = (-0.1 * n, 1.5 * n)
    mi.limits["nb"] = (0, 2 * n)
    mi.limits["mu"] = (0.1446, 0.1462)
    mi.limits["s"] = (0.0002, 0.0025)
    mi.limits["a"] = (0.01, 3.0)
    mi.limits["b"] = (-50, 300)
    if fix_shape:
        mi.fixed["mu"] = True
        mi.fixed["s"] = True
    mi.migrad()
    mi.hesse()
    return {"N": float(mi.values["ns"]), "N_err": float(mi.errors["ns"]), "mu": float(mi.values["mu"]),
            "sigma": float(mi.values["s"]), "valid": bool(mi.valid), "n_entries": int(h.sum())}
