"""NDS: skan wariantów N_D* na pełnej lumi (histogramy Δm z `charmtag` validate / dstar_hist) → ndstar_full.json/.md.
Użycie: python scripts/studies/ndstar_variants.py dstar_nom.npz dstar_post.json OUT_DIR  (npz: edges, norm_all, norm_<bin>, hist. MC sygnału).
Warianty: model tła (nominał próg×pol2; pełny z e^{-bu}; czysty e^{-bu}; √u×pol3), zakres, binowanie, sygnał (2G / G+CB; skala σ stała),
okno masy D0 (±36/±55 — inna selekcja, raport), per bin pT. Systematyka = suma kwadratowa max |ΔN| per kategoria (okno wyłączone,
model „czysty wykładnik” liczony osobno — Δχ²=+5 i σ×1,17 = artefakt fit_dm)."""
import sys, json, os, numpy as np
from iminuit import Minuit
from scipy.stats import norm, crystalball
import dstar_refit_nds as R  # kopia NDS (nie charmtag z PF/src)
MPI = R.MPI
# rerun 2026-11-12: NLL_MODE=new (domyślnie) → NLL względem modelu nasyconego, Σ[λ − h − h·ln(λ/h)] (jak dstar_refit z 10.11);
# NLL_MODE=old → stare Σ[λ − h·ln λ] (do odtworzenia ndstar_full.json). Stała addytywna — minimum bez zmian.
NLL_MODE = os.environ.get("NLL_MODE", "new")
def _pnll(lam, hh):
    if NLL_MODE == "old":
        return float(np.sum(lam - hh * np.log(lam)))
    hs = np.where(hh > 0, hh, 1.0)
    return float(np.sum(lam - hh - hh * np.log(lam / hs)))

def sig_cdf_gcb(x, mu, sc, k, f, beta, mcb):
    return f * norm.cdf(x, mu, sc) + (1 - f) * (1 - crystalball.cdf(-(x - mu) / (sc * k), beta, mcb))

def bkg_cdf_gen(xe, kind, a, b, c1, c2):
    xe = np.asarray(xe, float); lo, hi = xe[:-1], xe[1:]
    x = 0.5 * (hi + lo)[:, None] + 0.5 * (hi - lo)[:, None] * R._GL_X[None, :]
    u = np.clip(x - MPI, 1e-9, None); t = u / 20.0
    if kind == "C":   fb = u ** a * np.exp(-b * u)
    elif kind == "E": fb = np.sqrt(u) * (1 + a * t + b * t * t + c1 * t ** 3)
    fb = np.clip(fb, 1e-12, None); ib = 0.5 * (hi - lo) * (fb * R._GL_W[None, :]).sum(1)
    cb = np.concatenate([[0.0], np.cumsum(ib)]); return cb / cb[-1]

def rebin(h, e, lo, hi, merge=1):
    h = np.asarray(h, float); e = np.asarray(e, float)
    if merge > 1:
        n = (len(h) // merge) * merge; h = h[:n].reshape(-1, merge).sum(1); e = e[:n + 1:merge]
    sel = (e[:-1] >= lo - 1e-9) & (e[1:] <= hi + 1e-9)
    return h[sel], np.r_[e[:-1][sel], e[1:][sel][-1]]

def fit_generic(h, e, shape, bkg, sig="2G", lo=R.LO, hi=160.0, merge=1, s_free=True):
    """Warianty poza fit_single: inne modele tła (C, E), sygnał G+CB, skala stała, zakres/binowanie (dla nominalnego modelu tła → fit_single)."""
    hh, ee = rebin(h, e, lo, hi, merge)
    if bkg in ("A", "B") and sig == "2G" and s_free:
        f = R.fit_single(hh, ee, shape, attempts=(("threshold_pol2", False, 2),) if bkg == "B" else (("expo_pol2", True, 2), ("expo_pol2_s1", True, 1)))
        return {**f, "bin_mev": float(np.diff(ee)[0])}
    def cdf_s(xe, mu, s):
        c = R._sig_cdf(xe, mu, shape["sc"] * s, shape["k"], shape["f"]) if sig == "2G" else sig_cdf_gcb(xe, mu, shape["sc"] * s, shape["k"], shape["f"], shape["beta"], shape["mcb"])
        return (c - c[0]) / (c[-1] - c[0])
    def cdf_b(xe, a, b, c1, c2):
        return R._bkg_cdf_t(xe, a, b, c1, c2, bkg == "A") if bkg in ("A", "B") else bkg_cdf_gen(xe, bkg, a, b, c1, c2)
    def nll(N, mu, s, nb, a, b, c1, c2):
        lam = np.clip(np.diff(N * cdf_s(ee, mu, s) + nb * cdf_b(ee, a, b, c1, c2)), 1e-9, None); return _pnll(lam, hh)
    best = None
    for a0, b0 in ((0.5, 0.03), (0.6, 0.1), (0.4, 0.0)):
        m = Minuit(nll, N=max(0.01 * hh.sum(), 50), mu=shape["mu"], s=1.0, nb=hh.sum(), a=a0, b=b0, c1=0.0, c2=0.0); m.errordef = 0.5; m.strategy = 2
        m.limits["N"] = (0, None); m.limits["nb"] = (0, None); m.limits["s"] = (0.8, 1.6); m.limits["mu"] = (shape["mu"] - 0.5, shape["mu"] + 0.5)
        m.limits["a"] = (0.01, 3); m.limits["b"] = (-0.5, 1.5); m.limits["c1"] = (-5, 20); m.limits["c2"] = (-10, 40)
        if bkg == "B": m.fixed["b"] = True; m.values["b"] = 0.0
        if bkg == "C": m.fixed["c1"] = True; m.fixed["c2"] = True
        if bkg == "E": m.limits["a"] = (-5, 20); m.limits["b"] = (-20, 40); m.limits["c1"] = (-40, 40); m.fixed["c2"] = True; m.values["a"] = 0.0; m.values["b"] = 0.0
        if not s_free: m.fixed["s"] = True
        m.simplex(); m.migrad(ncall=50000); m.hesse()
        if best is None or (m.valid and (not best.valid or m.fval < best.fval)): best = m
        if best.valid: break
    m = best; v = {k: float(m.values[k]) for k in m.parameters}
    lam = np.diff(v["N"] * cdf_s(ee, v["mu"], v["s"]) + v["nb"] * cdf_b(ee, v["a"], v["b"], v["c1"], v["c2"]))
    return {"N": v["N"], "N_err": float(m.errors["N"]), "scale": v["s"], "mu": v["mu"], "valid": bool(m.valid), "posdef": bool(m.fmin.has_posdef_covar),
            "edm": float(m.fmin.edm), "chi2": float(np.sum((hh - lam) ** 2 / np.maximum(lam, 1))), "ndf": int(len(hh) - m.nfit), "n_entries": int(hh.sum()),
            "bkg_model": {"A": "expo_pol2", "B": "threshold_pol2", "C": "pure_expo", "E": "sqrt_pol3"}[bkg], "bin_mev": float(np.diff(ee)[0])}

def gcb_shape_from_mc(hmc, e, s2g):
    hh, ee = rebin(hmc, e, 143.0, 148.5)
    def nll(n, mu, sc, k, f, beta, mcb):
        c = sig_cdf_gcb(ee, mu, sc, k, f, beta, mcb); lam = np.clip(np.diff(n * c), 1e-9, None); return _pnll(lam, hh)
    m = Minuit(nll, n=hh.sum(), mu=s2g["mu"], sc=s2g["sc"], k=2.0, f=0.6, beta=1.5, mcb=3.0); m.errordef = 0.5; m.strategy = 2
    m.limits["sc"] = (0.2, 2); m.limits["k"] = (1.0, 6); m.limits["f"] = (0.05, 0.98); m.limits["beta"] = (0.3, 5); m.limits["mcb"] = (1.1, 20); m.limits["mu"] = (144.8, 146.2)
    m.simplex(); m.migrad(ncall=50000); m.hesse()
    return {p: float(m.values[p]) for p in ("mu", "sc", "k", "f", "beta", "mcb")} | {"valid": bool(m.valid)}

