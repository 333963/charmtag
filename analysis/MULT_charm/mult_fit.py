"""mult_fit.py — MULT_charm: binowane fity masy (rozszerzona NLL Poissona, scipy), niezależne od charmtag.

Modele (x w GeV):
  dstar : x = Δm ∈ [0,1395; 0,165]; S = Gauss(μ, σ₁) ⊕ Gauss(μ, σ₂=k·σ₁) z ułamkiem f; B = (x−m_π)^p·exp(−c(x−m_π))
  dplus : x = m(Kππ) ∈ [1,73; 2,00]; S = Gauss(μ, σ); B = exp(c·(x−1,87))·(1 + a·(x−1,87))
  d0    : x = m(Kπ) ∈ [1,75; 2,00]; S = Gauss(μ, σ); R (refleksja K↔π) = r·N_S·Gauss(μ_r, σ_r), r/μ_r/σ_r z MC (stałe); B jak dplus
Procedura MULT-v1: fit zintegrowany (wszystkie biny krotności) → μ, σ (i f, k dla D*) ustalone; w binach wolne N_S, N_B i parametry tła.
Błąd N_S z odwrotnego hesjanu (numerycznego) NLL po wszystkich wolnych parametrach.
Poprawka 19.11 (bieg lxplus 18.11: kształt bez więzów uciekał — D⁺ μ = 135 GeV, D0 σ = 41–55 MeV, D* μ na krawędzi okna — w MC i w danych):
  kształt referencyjny = fit samego sygnału z kandydatów MC z prawdą (`truth_shape`); fit zintegrowany ma μ i σ₁ w przedziałach
  wokół referencji (transformacja sigmoidalna: μ_ref ± DMU, σ₁ ∈ [SLO; SHI]·σ_ref; D*: f, k z referencji), flaga `at_bound`.
  Bez referencji — przedziały wokół masy PDG (DEF_BOUNDS).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr

from mult_common import DM_DSTAR, M_D0, M_DP, M_PI

WINDOWS = {"dstar": (0.1395, 0.165, 51), "dplus": (1.73, 2.00, 54), "d0": (1.75, 2.00, 50)}
SUB = 6
DMU = {"dstar": 0.0004, "dplus": 0.006, "d0": 0.006}            # ± przesunięcie μ dane/MC względem referencji [GeV]
SLO, SHI = 0.75, 1.5                                             # σ₁ / σ₁,ref
DEF_BOUNDS = {"dstar": (0.0002, 0.0012), "dplus": (0.006, 0.030), "d0": (0.006, 0.030)}   # σ₁ bez referencji


def _sig(t, lo, hi):
    return lo + (hi - lo) / (1.0 + np.exp(-np.clip(t, -50, 50)))


def _isig(v, lo, hi):
    u = np.clip((v - lo) / (hi - lo), 1e-6, 1 - 1e-6)
    return float(np.log(u / (1 - u)))


def _gbin(e, mu, s):
    return ndtr((e[1:] - mu) / s) - ndtr((e[:-1] - mu) / s)


def _bkg_shape(model, e, p):
    xs = (e[:-1, None] + (np.arange(SUB) + 0.5)[None, :] / SUB * np.diff(e)[:, None])
    if model == "dstar":
        q = np.clip(xs - M_PI, 1e-6, None)
        f = q ** p[0] * np.exp(-p[1] * q * 100.0)
    else:
        u = xs - 1.87
        f = np.exp(p[0] * u) * np.clip(1.0 + p[1] * u, 1e-9, None)
    b = f.mean(axis=1) * np.diff(e)
    return b / b.sum()


class Fitter:
    def __init__(self, model, refl=None):
        self.model = model
        self.lo, self.hi, self.nb = WINDOWS[model]
        self.e = np.linspace(self.lo, self.hi, self.nb + 1)
        self.refl = refl or {"r": 1.0, "mu": M_D0, "sigma": 0.09}

    def hist(self, x, w=None):
        h, _ = np.histogram(x, self.e, weights=w)
        return h

    def sig_shape(self, sh):
        s = self.sig_shape_core(sh)
        if self.model == "d0":
            s = s + self.refl["r"] * _gbin(self.e, self.refl["mu"], self.refl["sigma"])
        return s

    def bounds(self, ref=None):
        """Przedziały (μ, σ₁) dla fitu zintegrowanego."""
        m0 = DM_DSTAR if self.model == "dstar" else (M_DP if self.model == "dplus" else M_D0)
        if ref is None:
            return (m0 - 3 * DMU[self.model], m0 + 3 * DMU[self.model]), DEF_BOUNDS[self.model]
        return (ref["mu"] - DMU[self.model], ref["mu"] + DMU[self.model]), (SLO * ref["s1"], SHI * ref["s1"])

    def _unpack(self, p, free_shape, sh0, bnd=None):
        sh = dict(sh0)
        if free_shape:
            (ml, mh), (sl, shh) = bnd
            sh["mu"], sh["s1"] = _sig(p[4], ml, mh), _sig(p[5], sl, shh)
            if self.model == "dstar" and len(p) > 6:
                sh["f"], sh["k"] = _sig(p[6], 0.05, 0.95), _sig(p[7], 1.2, 5.0)
        return sh

    def truth_shape(self, x_true):
        """Kształt sygnału z kandydatów MC z prawdą (sam sygnał, bez tła; μ, σ₁ [, f, k] wolne w szerokich przedziałach)."""
        h = np.asarray(self.hist(x_true), float)
        if h.sum() < 50:
            return None
        m0 = DM_DSTAR if self.model == "dstar" else (M_DP if self.model == "dplus" else M_D0)
        bnd = ((m0 - 0.02 * (0.1 if self.model == "dstar" else 1)), (m0 + 0.02 * (0.1 if self.model == "dstar" else 1))), \
            ((0.0001, 0.003) if self.model == "dstar" else (0.003, 0.06))
        p0 = [_isig(m0, *bnd[0]), _isig(0.0006 if self.model == "dstar" else 0.015, *bnd[1])] + ([0.0, 0.0] if self.model == "dstar" else [])

        def nll(q):
            sh = self._unpack(np.r_[0, 0, 0, 0, q], True, {"mu": 0, "s1": 0}, bnd)
            s = np.clip(self.sig_shape_core(sh), 1e-12, None); lam = h.sum() * s / s.sum()
            return float((lam - h * np.log(lam)).sum())
        r = minimize(nll, np.array(p0, float), method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-8, "fatol": 1e-8})
        sh = self._unpack(np.r_[0, 0, 0, 0, r.x], True, {"mu": 0, "s1": 0}, bnd)
        return {k: float(v) for k, v in sh.items()}

    def sig_shape_core(self, sh):
        """Kształt sygnału bez refleksji D0."""
        if self.model == "dstar":
            return sh["f"] * _gbin(self.e, sh["mu"], sh["s1"]) + (1 - sh["f"]) * _gbin(self.e, sh["mu"], sh["s1"] * sh["k"])
        return _gbin(self.e, sh["mu"], sh["s1"])

    def fit(self, h, shape=None, ref=None):
        """h: liczniki (nieważone). shape=None → fit kształtu (zintegrowany, μ/σ₁ w przedziałach wokół ref). Zwraca dict (N_S, err, shape, nll, …)."""
        h = np.asarray(h, float); ntot = max(h.sum(), 1.0)
        free = shape is None
        bnd = self.bounds(ref) if free else None
        if free:
            sh0 = dict(ref) if ref else ({"mu": DM_DSTAR, "s1": 0.0006, "f": 0.6, "k": 2.2} if self.model == "dstar"
                                         else {"mu": M_DP if self.model == "dplus" else M_D0, "s1": 0.012})
        else:
            sh0 = shape
        b0 = [0.5, 8.0] if self.model == "dstar" else [-2.0, 0.0]
        p0 = [0.05 * ntot, 0.95 * ntot] + b0
        if free:
            p0 += [_isig(sh0["mu"], *bnd[0]), _isig(sh0["s1"], *bnd[1])]
            if self.model == "dstar" and ref is None:
                p0 += [_isig(sh0["f"], 0.05, 0.95), _isig(sh0["k"], 1.2, 5.0)]

        def nll(p):
            sh = self._unpack(p, free, sh0, bnd)
            lam = p[0] * self.sig_shape(sh) + abs(p[1]) * _bkg_shape(self.model, self.e, p[2:4])
            lam = np.clip(lam, 1e-12, None)
            return float((lam - h * np.log(lam)).sum())

        best = None
        for s0 in (0.02, 0.05, 0.15):
            p0[0], p0[1] = s0 * ntot, (1 - s0) * ntot
            r = minimize(nll, np.array(p0, float), method="Nelder-Mead", options={"maxiter": 20000, "maxfev": 20000, "xatol": 1e-7, "fatol": 1e-7})
            r = minimize(nll, r.x, method="Powell", options={"maxiter": 20000, "xtol": 1e-8, "ftol": 1e-10})
            if best is None or r.fun < best.fun:
                best = r
        p = best.x
        H = _hess(nll, p)
        err = float("nan")
        for k in (len(p), 4):                                      # pełny hesjan; gdy osobliwy (parametr kształtu przy granicy) — przy ustalonym kształcie
            try:
                v = np.linalg.inv(H[:k, :k])[0, 0]
            except np.linalg.LinAlgError:
                continue
            if np.isfinite(v) and v > 0:
                err = float(np.sqrt(v)); break
        err_fallback = not np.isfinite(err)
        if err_fallback:
            err = float(np.sqrt(max(p[0], 1.0)))
        out = {"N_S": float(p[0]), "err": err, "N_B": float(abs(p[1])), "shape": {k: float(v) for k, v in self._unpack(p, free, sh0, bnd).items()},
               "nll": float(best.fun), "ok": bool(np.isfinite(best.fun)) and not err_fallback, "n_entries": float(h.sum()), "err_fallback": err_fallback}
        if free:
            out["at_bound"] = bool(np.any(np.abs(p[4:6]) > 6.0))     # parametr przy granicy przedziału (|t| > 6 ⇔ < 0,25 % od krawędzi)
            out["bounds"] = {"mu": list(bnd[0]), "s1": list(bnd[1])}
            out["ok"] = out["ok"] and not out["at_bound"]
        return out


def _hess(f, p, rel=1e-4):
    n = len(p); H = np.zeros((n, n)); f0 = f(p)
    hstep = np.maximum(np.abs(p) * rel, 1e-7)
    for i in range(n):
        for j in range(i, n):
            pp = p.copy(); pp[i] += hstep[i]; pp[j] += hstep[j]; fpp = f(pp)
            pm = p.copy(); pm[i] += hstep[i]; pm[j] -= hstep[j]; fpm = f(pm)
            mp = p.copy(); mp[i] -= hstep[i]; mp[j] += hstep[j]; fmp = f(mp)
            mm = p.copy(); mm[i] -= hstep[i]; mm[j] -= hstep[j]; fmm = f(mm)
            H[i, j] = H[j, i] = (fpp - fpm - fmp + fmm) / (4 * hstep[i] * hstep[j])
    return H


def refl_from_mc(m_refl, m_true):
    """Parametry refleksji D0 z MC: Gauss dopasowany momentami do m(Kπ) wierszy z zamienioną hipotezą; r = N_refl/N_true."""
    m_refl = np.asarray(m_refl); m_true = np.asarray(m_true)
    lo, hi, _ = WINDOWS["d0"]
    s = (m_refl > lo) & (m_refl < hi)
    if s.sum() < 20 or len(m_true) < 20:
        return None
    return {"r": float(s.sum() / max(((m_true > lo) & (m_true < hi)).sum(), 1)), "mu": float(np.mean(m_refl[s])),
            "sigma": float(np.std(m_refl[s])), "n_refl": int(s.sum())}
