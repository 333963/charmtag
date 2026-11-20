"""Etap V: odporny fit Δm(D*) do kontroli D0 z tagu D* (zastępuje massfit.fit_dm w tej roli).
Model: sygnał = 2 Gaussy (kształt z MC truth D*: σ_core, k=σ_wide/σ_core, f_core; w danych wolna skala s≥0.8 wspólna dla pass/fail),
tło = (u)^a·exp(−b·u)·(1 + c1·u + c2·u²), u = Δm − m_π [MeV], osobne parametry tła w próbce pass i fail;
zakres [m_π+0.9 MeV, 160 MeV] (bez progu, gdzie model tła najsłabszy), biny 0.1 MeV.
ε_data = N_pass/(N_pass+N_fail) — JEDNOCZESNY fit pass/fail (ε i N_tot to parametry), więc brak osobnego fitu „bez cięcia”.
Wejście: npz z dstar_hist.py (h_<var>_<bin>_<wp>_pass/fail, edges) lub toy.
"""
import numpy as np
from iminuit import Minuit
from scipy.stats import norm

MPI = 139.57
LO, HI = MPI + 0.9, 160.0

def _sig_cdf(x, mu, sc, k, f):
    return f * norm.cdf(x, mu, sc) + (1 - f) * norm.cdf(x, mu, sc * k)

def _bkg_cdf(xe, a, b, c1, c2, n=4001):
    x = np.linspace(xe[0], xe[-1], n); u = np.clip(x - MPI, 1e-9, None)
    fb = np.clip(u ** a * np.exp(-b * u) * (1 + c1 * u + c2 * u * u), 1e-12, None)
    cb = np.concatenate([[0], np.cumsum(0.5 * (fb[1:] + fb[:-1]) * np.diff(x))])
    return np.interp(xe, x, cb / cb[-1])

def shape_from_mc(dm_mev_or_hist, edges=None):
    """Kształt sygnału z MC truth D* (Δm w MeV): 2 Gaussy, zwraca dict(mu, sc, k, f)."""
    from iminuit.cost import UnbinnedNLL, ExtendedBinnedNLL
    if edges is None:
        x = np.asarray(dm_mev_or_hist); x = x[(x > 143) & (x < 148.5)]
        pdf = lambda x, mu, sc, k, f: (f * norm.pdf(x, mu, sc) + (1 - f) * norm.pdf(x, mu, sc * k)) / (_sig_cdf(148.5, mu, sc, k, f) - _sig_cdf(143, mu, sc, k, f))
        m = Minuit(UnbinnedNLL(x, pdf), mu=145.43, sc=0.6, k=3, f=0.7)
    else:
        h = np.asarray(dm_mev_or_hist, float); e = np.asarray(edges, float)
        sel = (e[:-1] >= 143) & (e[1:] <= 148.5); hh = h[sel]; ee = np.r_[e[:-1][sel], e[1:][sel][-1]]
        cdf = lambda xe, n, mu, sc, k, f: n * (_sig_cdf(xe, mu, sc, k, f) - _sig_cdf(ee[0], mu, sc, k, f))
        m = Minuit(ExtendedBinnedNLL(hh, ee, cdf), n=hh.sum(), mu=145.43, sc=0.6, k=3, f=0.7)
    m.limits["sc"] = (0.15, 2); m.limits["k"] = (1.2, 8); m.limits["f"] = (0.05, 1); m.limits["mu"] = (144.8, 146.2)
    m.migrad(); m.hesse()
    return {p: float(m.values[p]) for p in ("mu", "sc", "k", "f")}

def _rebin(h, e):
    e = np.asarray(e, float); h = np.asarray(h, float)
    sel = (e[:-1] >= LO - 1e-9) & (e[1:] <= HI + 1e-9)
    return h[sel], np.r_[e[:-1][sel], e[1:][sel][-1]]

def fit_passfail(hp, hf, edges, shape, float_scale=True, mu_free=True):
    """hp/hf: histogramy Δm [MeV] pass/fail na wspólnych krawędziach. Zwraca ε, N_tot, błędy, χ²."""
    hp, e = _rebin(hp, edges); hf, _ = _rebin(hf, edges)
    def cdf_s(xe, mu, s):
        c = _sig_cdf(xe, mu, shape["sc"] * s, shape["k"], shape["f"]); return (c - c[0]) / (c[-1] - c[0])
    def nll(Nt, eps, mu, s, nbp, ap, bp, c1p, c2p, nbf, af, bf, c1f, c2f):
        cs = cdf_s(e, mu, s)
        lp = np.diff(Nt * eps * cs + nbp * _bkg_cdf(e, ap, bp, c1p, c2p))
        lf = np.diff(Nt * (1 - eps) * cs + nbf * _bkg_cdf(e, af, bf, c1f, c2f))
        lp = np.clip(lp, 1e-9, None); lf = np.clip(lf, 1e-9, None)
        return float(np.sum(lp - hp * np.log(lp)) + np.sum(lf - hf * np.log(lf)))
    best = None
    for b0 in (0.05, 0.15, 0.4):
        m = Minuit(nll, Nt=max(0.01 * (hp.sum() + hf.sum()), 50), eps=0.7, mu=shape["mu"], s=1.0,
                   nbp=hp.sum(), ap=0.5, bp=b0, c1p=0.0, c2p=0.0, nbf=hf.sum(), af=0.5, bf=b0, c1f=0.0, c2f=0.0)
        m.errordef = 0.5
        m.limits["eps"] = (-0.2, 1.5); m.limits["Nt"] = (0, None); m.limits["s"] = (0.8, 1.6); m.limits["mu"] = (shape["mu"] - 0.5, shape["mu"] + 0.5)
        for p in ("ap", "af"): m.limits[p] = (0.01, 3)
        for p in ("bp", "bf"): m.limits[p] = (-1, 3)
        for p in ("c1p", "c1f"): m.limits[p] = (-0.5, 2)
        for p in ("c2p", "c2f"): m.limits[p] = (-0.05, 0.2)
        for p in ("nbp", "nbf"): m.limits[p] = (0, None)
        if not float_scale: m.fixed["s"] = True
        if not mu_free: m.fixed["mu"] = True
        m.migrad(ncall=20000); m.hesse()
        if best is None or (m.valid and m.fval < best.fval) or (not best.valid and m.valid):
            best = m
        if best.valid:
            break   # kolejne starty tylko gdy fit nie zbiegł
    m = best
    v = {k: float(m.values[k]) for k in m.parameters}
    cs = cdf_s(e, v["mu"], v["s"])
    lp = np.diff(v["Nt"] * v["eps"] * cs + v["nbp"] * _bkg_cdf(e, v["ap"], v["bp"], v["c1p"], v["c2p"]))
    lf = np.diff(v["Nt"] * (1 - v["eps"]) * cs + v["nbf"] * _bkg_cdf(e, v["af"], v["bf"], v["c1f"], v["c2f"]))
    chi2 = float(np.sum((hp - lp) ** 2 / np.maximum(lp, 1)) + np.sum((hf - lf) ** 2 / np.maximum(lf, 1)))
    return {"eps": v["eps"], "eps_err": float(m.errors["eps"]), "N_tot": v["Nt"], "N_tot_err": float(m.errors["Nt"]),
            "scale": v["s"], "mu": v["mu"], "valid": bool(m.valid), "chi2": chi2, "ndf": int(2 * len(hp) - m.nfit),
            "n_pass": int(hp.sum()), "n_fail": int(hf.sum())}

def sb_count(hp, hf, edges, mu=145.43, half=1.0, sb=(150.0, 158.0), fit_lo=MPI + 0.9):
    """Kontrola bez kształtu sygnału: N(okno) − tło interpolowane z fitu samego tła do (zakres poza oknem ±2.5 MeV)."""
    out = {}
    for nm, h in (("pass", hp), ("fail", hf)):
        h, e = _rebin(h, edges); c = 0.5 * (e[1:] + e[:-1])
        side = (np.abs(c - mu) > 2.5)
        def nll(nb, a, b, c1, c2):
            lam = np.clip(nb * np.diff(_bkg_cdf(e, a, b, c1, c2)), 1e-9, None)[side]
            return float(np.sum(lam - h[side] * np.log(lam)))
        m = Minuit(nll, nb=h.sum(), a=0.5, b=0.1, c1=0.0, c2=0.0); m.errordef = 0.5
        m.limits["a"] = (0.01, 3); m.limits["b"] = (-1, 3); m.limits["c1"] = (-0.5, 2); m.limits["c2"] = (-0.05, 0.2); m.limits["nb"] = (0, None)
        m.migrad()
        lam = m.values["nb"] * np.diff(_bkg_cdf(e, *[m.values[k] for k in ("a", "b", "c1", "c2")]))
        w = np.abs(c - mu) < half
        out[nm] = float(h[w].sum() - lam[w].sum()); out[nm + "_err"] = float(np.sqrt(h[w].sum()))
    t = out["pass"] + out["fail"]
    out["eps"] = out["pass"] / t if t > 0 else float("nan")
    return out

def fit_single(h, edges, shape):
    """Jedna próbka (np. N_D* do normalizacji par): 2 Gaussy (kształt MC, skala wolna) + tło jak w fit_passfail."""
    h, e = _rebin(h, edges)
    def cdf_s(xe, mu, s):
        c = _sig_cdf(xe, mu, shape["sc"] * s, shape["k"], shape["f"]); return (c - c[0]) / (c[-1] - c[0])
    def nll(N, mu, s, nb, a, b, c1, c2):
        lam = np.clip(np.diff(N * cdf_s(e, mu, s) + nb * _bkg_cdf(e, a, b, c1, c2)), 1e-9, None)
        return float(np.sum(lam - h * np.log(lam)))
    best = None
    for b0 in (0.05, 0.15, 0.4):
        m = Minuit(nll, N=max(0.01 * h.sum(), 50), mu=shape["mu"], s=1.0, nb=h.sum(), a=0.5, b=b0, c1=0.0, c2=0.0)
        m.errordef = 0.5
        m.limits["N"] = (0, None); m.limits["s"] = (0.8, 1.6); m.limits["mu"] = (shape["mu"] - 0.5, shape["mu"] + 0.5)
        m.limits["a"] = (0.01, 3); m.limits["b"] = (-1, 3); m.limits["c1"] = (-0.5, 2); m.limits["c2"] = (-0.05, 0.2); m.limits["nb"] = (0, None)
        m.migrad(ncall=20000); m.hesse()
        if best is None or (m.valid and m.fval < best.fval) or (not best.valid and m.valid):
            best = m
        if best.valid:
            break
    m = best; v = {k: float(m.values[k]) for k in m.parameters}
    lam = np.diff(v["N"] * cdf_s(e, v["mu"], v["s"]) + v["nb"] * _bkg_cdf(e, v["a"], v["b"], v["c1"], v["c2"]))
    return {"N": v["N"], "N_err": float(m.errors["N"]), "scale": v["s"], "mu": v["mu"], "valid": bool(m.valid),
            "chi2": float(np.sum((h - lam) ** 2 / np.maximum(lam, 1))), "ndf": int(len(h) - m.nfit), "n_entries": int(h.sum())}

def old_fit(h, edges, **kw):
    """Odtworzenie starej metody (massfit.fit_dm, Gauss σ≤2,5 MeV + próg 2-parametrowy) z histogramu [MeV]."""
    from charmtag.massfit import fit_dm
    c = 0.5 * (np.asarray(edges[1:]) + np.asarray(edges[:-1]))
    return fit_dm(np.repeat(c, np.asarray(h, int)) / 1e3, **kw)


# ----------------------------------------------------------------------------- AA: hak do pipeline.dstar_control
DEFAULT_SHAPE = {"mu": 145.43, "sc": 0.6, "k": 3.0, "f": 0.7}   # start; właściwy kształt: shape_from_mc (V dstar_post)


def control_row(j, thr, wps, cfg, species, M, step_mev=0.5, hi_mev=160.0):
    """Wiersz kontroli D* dla pipeline.dstar_control (opcja validate.dstar_refit: true).
    j: kandydaci D0↔D* w binie pT (kolumny dm [GeV], score, pt). Jednoczesny fit pass/fail (etap V) per WP."""
    vc = cfg.get("validate", {})
    shape = dict(DEFAULT_SHAPE, **(vc.get("dstar_shape") or {}))
    edges = np.arange(MPI + 0.9, hi_mev + 1e-9, step_mev)
    dm = j["dm"].to_numpy(float) * 1e3
    row = {"method": "V: jednoczesny fit pass/fail (2 Gaussy MC × skala + u^a e^{-bu}(1+c1u+c2u²))", "shape": shape}
    h_all, _ = np.histogram(dm, edges)
    for wp in wps:
        p = M.pass_wp(j["score"].to_numpy(), j["pt"].to_numpy(), thr, wp)
        hp, _ = np.histogram(dm[p], edges)
        hf = h_all - hp
        try:
            r = fit_passfail(hp, hf, edges, shape, float_scale=vc.get("dstar_float_scale", True), mu_free=True)
        except Exception as ex:  # fit nie zbiegł → NaN, nie wywracać validate
            r = {"eps": float("nan"), "eps_err": float("nan"), "valid": False, "error": str(ex)}
        row[f"refit_{int(wp*100)}"] = r
        row[f"eff_{int(wp*100)}"] = float(r["eps"])
        row[f"err_{int(wp*100)}"] = float(r["eps_err"])
    return row
