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

SHAPE_K_MIN = 1.2


def shape_from_mc(dm_mev_or_hist, edges=None, k_max: float = 3.0):
    """Kształt sygnału z MC truth D* (Δm w MeV): 2 Gaussy, zwraca dict(mu, sc, k, f, k_at_limit).
    Fit tylko w oknie 143–148,5 MeV, k = σ_wide/σ_core ∈ (1,2; k_max). Limit k_max = 3 (było 8, 2026-11-13): w oknie 5,5 MeV szeroki
    Gauss z dużym k jest nieokreślony (modeluje płaski cokół), a w ficie pass/fail na [m_π+0,9; 160] MeV pochłania tło (closure MC 8–30 GeV:
    ε_fit/ε_truth 0,957–0,989 przy k = 7,9; 0,978–0,995 przy k ≤ 4; 0,991–0,996 przy k ≤ 3). Kształty z k ≈ 2,2–2,74 bez zmian
    (NUM_audit_2026-11-12/kmax). k_at_limit: k w odległości < 1 % zakresu od k_max (kształt do sprawdzenia)."""
    from iminuit.cost import UnbinnedNLL, ExtendedBinnedNLL
    if not k_max > SHAPE_K_MIN:
        raise ValueError(f"shape_from_mc: k_max = {k_max} musi być > {SHAPE_K_MIN}")
    k0 = min(3.0, k_max - 0.1 * (k_max - SHAPE_K_MIN))   # start k = 3 jak dotąd (dla k_max ≥ 3,2); niżej 10 % zakresu pod limitem
    if edges is None:
        x = np.asarray(dm_mev_or_hist); x = x[(x > 143) & (x < 148.5)]
        pdf = lambda x, mu, sc, k, f: (f * norm.pdf(x, mu, sc) + (1 - f) * norm.pdf(x, mu, sc * k)) / (_sig_cdf(148.5, mu, sc, k, f) - _sig_cdf(143, mu, sc, k, f))
        m = Minuit(UnbinnedNLL(x, pdf), mu=145.43, sc=0.6, k=k0, f=0.7)
    else:
        h = np.asarray(dm_mev_or_hist, float); e = np.asarray(edges, float)
        sel = (e[:-1] >= 143) & (e[1:] <= 148.5); hh = h[sel]; ee = np.r_[e[:-1][sel], e[1:][sel][-1]]
        cdf = lambda xe, n, mu, sc, k, f: n * (_sig_cdf(xe, mu, sc, k, f) - _sig_cdf(ee[0], mu, sc, k, f))
        m = Minuit(ExtendedBinnedNLL(hh, ee, cdf), n=hh.sum(), mu=145.43, sc=0.6, k=k0, f=0.7)
    m.limits["sc"] = (0.15, 2); m.limits["k"] = (SHAPE_K_MIN, k_max); m.limits["f"] = (0.05, 1); m.limits["mu"] = (144.8, 146.2)
    m.migrad(); m.hesse()
    out = {p: float(m.values[p]) for p in ("mu", "sc", "k", "f")}
    out["k_at_limit"] = bool(out["k"] > k_max - 0.01 * (k_max - SHAPE_K_MIN))
    return out

def _rebin(h, e):
    e = np.asarray(e, float); h = np.asarray(h, float)
    sel = (e[:-1] >= LO - 1e-9) & (e[1:] <= HI + 1e-9)
    return h[sel], np.r_[e[:-1][sel], e[1:][sel][-1]]

def _nll_sat(lam, n, pos):
    """Poissonowskie NLL względem modelu nasyconego w postaci dokładnej przy dużych zliczeniach:
    Σ n·(d − log1p d), d = λ/n − 1 (dla n = 0: λ). Równe Σ[λ − n − n ln(λ/n)], ale bez kasowania λ − n przy ≳10⁴ wpisów/bin
    (postać λ − n − n ln(λ/n) daje tam szum HESSE — audyt NLL 2026-11-12 §3). pos = n > 0 (liczone raz)."""
    d = lam[pos] / n[pos] - 1.0
    return float(np.sum(n[pos] * (d - np.log1p(d))) + np.sum(lam[~pos]))


def _rebin_pf(h, e, tol=1e-6):
    """Jak `_rebin`, ale z tolerancją 1e-6 MeV (≪ szerokość binu): przesunięcie krawędzi o ±1e-9 nie zmienia liczby binów
    (N2: krawędź leży 6e-14 MeV od LO, więc tolerancja 1e-9 z `_rebin` gubiła bin przy przesunięciu −1e-9)."""
    e = np.asarray(e, float); h = np.asarray(h, float)
    sel = (e[:-1] >= LO - tol) & (e[1:] <= HI + tol)
    return h[sel], np.r_[e[:-1][sel], e[1:][sel][-1]]


# fit_passfail (2026-11-12): kolejka modeli tła jak w fit_single — nominał „próg × pol2” (b ≡ 0, zagnieżdżony), wykładnik tylko
# gdy żaden start nominału nie zbiegł. Pełny model ma płaski kierunek b ↔ (c1, c2): z nim ~50 % fitów pass/fail było niepoprawnych.
FIT_PASSFAIL_ATTEMPTS = (("threshold_pol2", False, 2), ("expo_pol2", True, 2), ("expo_pol2_s1", True, 1))
PF_EPS_LIMITS = (-0.2, 1.5)    # celowo poza [0, 1]: przy małym S nie obcinać ε (obciążenie); ε przy granicy → flaga eps_at_limit
PF_TOL = 0.01                  # Minuit tol → EDM < 2e-3·tol·errordef = 1e-5 (domyślnie 1e-4 ⇒ ε ruchome o ~1 % σ)
PF_SIDE = 2.5                  # [MeV] okno sygnału wyłączone w fitach samego tła (starty; jak 2,5 w sb_count)
PF_MIN_SIGNIF = 3.0            # N_tot/σ(N_tot) poniżej tego ⇒ sygnał nieistotny ⇒ ε nieokreślone (toy N2 3-5: przy Z ≈ 1,5–2 ε z „zamkniętym”
                               # przedziałem to głównie fluktuacje w dół — pull −2,6; granica ε = 1,5 obcina fluktuacje w górę)
PF_EXPO_DNLL = 4.6             # test modelu (tylko diagnostyka, wynik bez zmian): e^{-bu} (b_pass, b_fail) poprawia NLL o > 4,6
                               # (2ΔNLL > 9,2 ⇒ p < 1 % dla 2 st. sw.) → bkg_tension. Przełączanie modelu badane i odrzucone: na MC fd
                               # (mało tła) minimum modelu z wykładnikiem jest wielokrotne i wybór modelu zależał od perturbacji 1e-6.
PF_DNLL_AMBIG = 0.5            # drugie minimum bliżej niż ΔNLL 0,5 (1σ) i dalej niż 1σ w ε ⇒ ε nieokreślone


def _pf_param_limits(m, shape, expo):
    m.limits["eps"] = PF_EPS_LIMITS; m.limits["Nt"] = (0, None); m.limits["s"] = (0.8, 1.6)
    m.limits["mu"] = (shape["mu"] - 0.5, shape["mu"] + 0.5)
    for q in ("p", "f"):
        m.limits["a" + q] = (0.01, 3); m.limits["b" + q] = (-0.5, 1.5); m.limits["c1" + q] = (-5, 20)
        m.limits["c2" + q] = (-10, 40); m.limits["nb" + q] = (0, None)
        if not expo:
            m.fixed["b" + q] = True


def _pf_bkg_start(h, e, mu, expo):
    """Fit samego tła (próg × pol2[·e^{-bu}]) do binów poza |Δm − μ| < PF_SIDE: start tła i nadwyżka w oknie (start sygnału)."""
    c = 0.5 * (e[1:] + e[:-1]); side = np.abs(c - mu) > PF_SIDE; hs = h[side]; pos = hs > 0
    def nll(nb, a, b, c1, c2):
        lam = np.clip(nb * np.diff(_bkg_cdf_t(e, a, b, c1, c2, expo))[side], 1e-9, None)
        return _nll_sat(lam, hs, pos)
    best = None
    for a0 in (0.5, 1.5, 2.8):   # tło przy małej statystyce ma kilka minimów (parametry pol2 na granicach) — kilka startów
        m = Minuit(nll, nb=max(h.sum(), 1.0), a=a0, b=0.03 if expo else 0.0, c1=0.0, c2=0.0)
        m.errordef = 0.5; m.errors["nb"] = max(np.sqrt(h.sum()), 1.0)
        m.limits["nb"] = (0, None); m.limits["a"] = (0.01, 3); m.limits["b"] = (-0.5, 1.5); m.limits["c1"] = (-5, 20); m.limits["c2"] = (-10, 40)
        if not expo:
            m.fixed["b"] = True
        m.simplex(); m.migrad(ncall=20000)
        if best is None or m.fval < best.fval:
            best = m
    m = best
    v = {k: float(m.values[k]) for k in m.parameters}
    lam = v["nb"] * np.diff(_bkg_cdf_t(e, v["a"], v["b"], v["c1"], v["c2"], expo))
    return v, float(h[~side].sum() - lam[~side].sum())


def _pf_starts(hp, hf, e, shape, expo):
    """Starty fitu pass/fail (wykonywane WSZYSTKIE): tło z fitów sidebandów + ε, N_tot z nadwyżki w oknie; ten sam start tła
    z ε0 ∈ {0,2; 0,5; 0,9}; start „ogólny” (płaskie tło, ε0 = 0,7) jak w starej wersji."""
    bp, sp = _pf_bkg_start(hp, e, shape["mu"], expo)
    bf, sf = _pf_bkg_start(hf, e, shape["mu"], expo)
    c = _sig_cdf(np.array([shape["mu"] - PF_SIDE, shape["mu"] + PF_SIDE, e[0], e[-1]]), shape["mu"], shape["sc"], shape["k"], shape["f"])
    fw = (c[1] - c[0]) / (c[3] - c[2])                       # ułamek sygnału w oknie ±PF_SIDE
    ntot = hp.sum() + hf.sum()
    nt0 = max((sp + sf) / fw, 50.0)
    e0 = float(np.clip(sp / (sp + sf), 0.05, 0.95)) if sp + sf > 0 else 0.5
    bkg = {"nbp": bp["nb"], "ap": bp["a"], "bp": bp["b"], "c1p": bp["c1"], "c2p": bp["c2"],
           "nbf": bf["nb"], "af": bf["a"], "bf": bf["b"], "c1f": bf["c1"], "c2f": bf["c2"]}
    sig = {"mu": shape["mu"], "s": 1.0}
    st = [("sideband", {"Nt": nt0, "eps": e0, **sig, **bkg})]
    st += [(f"sideband_eps{x:g}", {"Nt": nt0, "eps": x, **sig, **bkg}) for x in (0.2, 0.5, 0.9)]
    b0 = 0.03 if expo else 0.0
    st.append(("generic", {"Nt": max(0.01 * ntot, 50.0), "eps": 0.7, **sig, "nbp": max(hp.sum(), 1.0), "ap": 0.5, "bp": b0, "c1p": 0.0,
                           "c2p": 0.0, "nbf": max(hf.sum(), 1.0), "af": 0.5, "bf": b0, "c1f": 0.0, "c2f": 0.0}))
    return st


PF_BKG_GRID = tuple((a0, c1, c2) for a0 in (0.3, 1.0, 2.0, 3.0) for c1, c2 in ((0.0, 0.0), (5.0, -3.0), (-2.0, 1.0)))


def _pf_polish_start(m, hp, hf, pp, pf_, e, cdf_s, expo):
    """Przy sygnale ustalonym z fitu m: tło pass i fail dopasowane osobno z siatki PF_BKG_GRID (i z wartości m). Zwraca start pełnego
    fitu (sygnał z m, najlepsze tła), gdy NLL którejś części spadło o > 1e-3; inaczej None."""
    v = {q: float(m.values[q]) for q in m.parameters}
    cs = cdf_s(v["mu"], v["s"]); st, gain = dict(v), 0.0
    for q, h, pos, sig in (("p", hp, pp, v["Nt"] * v["eps"] * cs), ("f", hf, pf_, v["Nt"] * (1 - v["eps"]) * cs)):
        def nll(nb, a, b, c1, c2):
            return _nll_sat(np.clip(sig + nb * np.diff(_bkg_cdf_t(e, a, b, c1, c2, expo)), 1e-9, None), h, pos)
        own = (v["nb" + q], v["a" + q], v["b" + q], v["c1" + q], v["c2" + q]); f0 = nll(*own); best = (f0, own)
        for a0, c10, c20 in PF_BKG_GRID:
            mm = Minuit(nll, nb=max(v["nb" + q], 1.0), a=a0, b=0.03 if expo else 0.0, c1=c10, c2=c20)
            mm.errordef = 0.5; mm.errors["nb"] = max(np.sqrt(max(v["nb" + q], 1.0)), 1.0)
            mm.limits["nb"] = (0, None); mm.limits["a"] = (0.01, 3); mm.limits["b"] = (-0.5, 1.5); mm.limits["c1"] = (-5, 20); mm.limits["c2"] = (-10, 40)
            if not expo:
                mm.fixed["b"] = True
            mm.migrad(ncall=10000)
            if mm.fval < best[0]:
                best = (float(mm.fval), tuple(float(x) for x in mm.values))
        if best[0] < f0 - 1e-3:
            gain += f0 - best[0]
            st.update(dict(zip(("nb" + q, "a" + q, "b" + q, "c1" + q, "c2" + q), best[1])))
    return st if gain > 0 else None


def _pf_profile_side(prof, eps0, err, side, f0):
    """Jedna strona przedziału z profilu NLL(ε) (ΔNLL = 0,5): kroki k·σ_HESSE na zewnątrz aż do przekroczenia 0,5 albo granicy ε,
    potem brentq. prof(x) → (NLL_min przy ε = x, wektor parametrów). Zwraca (ε_kraniec | NaN, min ΔNLL widziane, x tego minimum)."""
    from scipy.optimize import brentq
    lim = PF_EPS_LIMITS[1] if side > 0 else PF_EPS_LIMITS[0]
    x_in, d_in, dmin, xmin = eps0, 0.0, 0.0, eps0
    for kk in (0.5, 1.0, 1.5, 2.0, 3.0, 4.5, 7.0, 10.0):
        x = eps0 + side * kk * err
        x = min(x, lim) if side > 0 else max(x, lim)
        d = prof(x)[0] - f0
        if d < dmin:
            dmin, xmin = d, x
        if d >= 0.5:
            r = brentq(lambda t: prof(t)[0] - f0 - 0.5, x_in, x, xtol=1e-3 * err) if d_in < 0.5 else x_in
            return float(r), dmin, xmin
        x_in, d_in = x, d
        if x == lim:
            break
    return float("nan"), dmin, xmin


def _fit_passfail_full(hp, hf, edges, shape, float_scale=True, mu_free=True, attempts=FIT_PASSFAIL_ATTEMPTS, interval=True,
                      model_check=True):
    """Jak fit_passfail, ale zwraca (wynik, Minuit wybranego minimum) — do skanów profilu w badaniach i testach."""
    return _pf_core(hp, hf, edges, shape, float_scale, mu_free, attempts, interval, model_check)


def fit_passfail(hp, hf, edges, shape, float_scale=True, mu_free=True, attempts=FIT_PASSFAIL_ATTEMPTS, interval=True,
                 model_check=True):
    """Jednoczesny fit Δm pass/fail (sygnał: 2 Gaussy z MC, wspólne μ i skala σ; tło próg × pol2 osobno pass i fail) → ε, N_tot.
    hp/hf: histogramy Δm [MeV] na wspólnych krawędziach.

    Odporność (2026-11-12): tło u^a(1 + c1·t + c2·t²), t = u/20 MeV, całki GL (`_bkg_cdf_t`), nominał bez wykładnika (b ≡ 0; kolejka
    `attempts` jak w fit_single); NLL nasycone w postaci log1p; wszystkie starty z `_pf_starts` (simplex → migrad strategy 2 → hesse)
    i wybór najniższego NLL spośród poprawnych; przedział ε ze skanu profilu NLL(ε) (ΔNLL = 0,5, brentq) — nie z MINOS, bo MINOS
    potrafi zwrócić „poprawny” kraniec w złym miejscu.

    Klucze jak dotąd: eps, eps_err (HESSE), N_tot, N_tot_err, scale, mu, valid, chi2, ndf, n_pass, n_fail. Nowe:
    eps_lo/eps_hi (przedział z profilu; NaN gdy nie zamknięty przed granicą ε), eps_err_lo/eps_err_hi, posdef, edm, nll, bkg_model, eps_at_limit,
    params_at_limit (inne parametry na granicy — wtedy HESSE bywa niestabilne, np. skala σ = 0,8; miarodajny przedział z profilu),
    n_minima (różne minima ze startów), dnll_second / eps_second (drugie minimum), determined (ε wyznaczone przez dane),
    undetermined_reason, attempts (log modeli i startów).
    determined = valid ∧ posdef ∧ ε nie przy granicy ∧ przedział z profilu zamknięty w granicach ∧ brak drugiego minimum o ΔNLL < 0,5
    i |Δε| > σ (ze startów albo ze skanu profilu) ∧ N_tot > PF_MIN_SIGNIF·σ(N_tot). Bin nieokreślony zwraca najlepszy fit, ale NIE należy go używać jako pomiaru ε.
    interval=False: bez skanu profilu (szybciej; eps_lo/eps_hi = NaN, determined bez warunku przedziału).
    model_check=True: po nominale (próg × pol2) test wykładnika: expo_pol2 ze startu w minimum nominału (b0 = 0; 0,05; 0,15) + polerowanie;
    zwraca dnll_expo (NLL_nominał − NLL_expo), eps_expo (ε z modelu z wykładnikiem — miara syst. modelu tła) i bkg_tension
    (dnll_expo > PF_EXPO_DNLL: nominał niewystarczający, ε może być obciążone). Wynik ε zawsze z nominału."""
    return _pf_core(hp, hf, edges, shape, float_scale, mu_free, attempts, interval, model_check)[0]


def _pf_core(hp, hf, edges, shape, float_scale, mu_free, attempts, interval, model_check):
    hp, e = _rebin_pf(hp, edges); hf, _ = _rebin_pf(hf, edges)
    pp, pf_ = hp > 0, hf > 0
    sc, k, f = shape["sc"], shape["k"], shape["f"]
    def cdf_s(mu, s):
        c = _sig_cdf(e, mu, sc * s, k, f); return np.diff(c) / (c[-1] - c[0])
    def lams(v, expo):
        cs = cdf_s(v["mu"], v["s"])
        lp = v["Nt"] * v["eps"] * cs + v["nbp"] * np.diff(_bkg_cdf_t(e, v["ap"], v["bp"], v["c1p"], v["c2p"], expo))
        lf = v["Nt"] * (1 - v["eps"]) * cs + v["nbf"] * np.diff(_bkg_cdf_t(e, v["af"], v["bf"], v["c1f"], v["c2f"], expo))
        return np.clip(lp, 1e-9, None), np.clip(lf, 1e-9, None)
    names = ("Nt", "eps", "mu", "s", "nbp", "ap", "bp", "c1p", "c2p", "nbf", "af", "bf", "c1f", "c2f")

    def make(st, expo, strategy):
        def nll(Nt, eps, mu, s, nbp, ap, bp, c1p, c2p, nbf, af, bf, c1f, c2f):
            lp, lf = lams(dict(zip(names, (Nt, eps, mu, s, nbp, ap, bp, c1p, c2p, nbf, af, bf, c1f, c2f))), expo)
            return _nll_sat(lp, hp, pp) + _nll_sat(lf, hf, pf_)
        m = Minuit(nll, **{p: float(st[p]) for p in names})
        m.errordef = 0.5; m.strategy = strategy; m.tol = PF_TOL
        m.errors["Nt"] = max(np.sqrt(max(st["Nt"], 1.0)) * 2, 10.0); m.errors["eps"] = 0.1; m.errors["mu"] = 0.05; m.errors["s"] = 0.05
        for q in ("p", "f"):
            m.errors["nb" + q] = max(np.sqrt(max(st["nb" + q], 1.0)), 1.0); m.errors["a" + q] = 0.05
            m.errors["c1" + q] = 0.1; m.errors["c2" + q] = 0.1; m.errors["b" + q] = 0.01
        _pf_param_limits(m, shape, expo)
        if not float_scale:
            m.fixed["s"] = True
        if not mu_free:
            m.fixed["mu"] = True
        return m

    def polish(mb, expo, strategy, res):
        """Polerowanie najlepszego fitu (zob. _pf_polish_start), do 3 obiegów; dopisuje próby do res, zwraca najlepszy Minuit."""
        for it in range(3):
            st = _pf_polish_start(mb, hp, hf, pp, pf_, e, cdf_s, expo)
            if st is None:
                break
            m = make(st, expo, strategy)
            try:
                m.migrad(ncall=50000); m.hesse()
            except (RuntimeError, ValueError) as ex:
                res.append((f"polish{it}", None, str(ex))); break
            res.append((f"polish{it}", m, ""))
            if not (m.fval < mb.fval - 1e-3 and (m.valid or not mb.valid)):
                break
            mb = m
        return mb

    def run(expo, strategy):
        res = []
        for label, st in _pf_starts(hp, hf, e, shape, expo):
            m = make(st, expo, strategy)
            try:
                m.simplex(); m.migrad(ncall=50000); m.hesse()
            except (RuntimeError, ValueError) as ex:   # np. NaN w NLL z dziwnego startu — pomijamy start
                res.append((label, None, str(ex))); continue
            res.append((label, m, ""))
        # polerowanie: tła pass i fail są niezależne przy ustalonym sygnale, a przy małej statystyce tła mają kilka lokalnych minimów
        # (parametry pol2 / a na granicach), do których starty trafiają chaotycznie (simplex). Dla najlepszego fitu: tło pass i fail
        # dopasowane osobno (sygnał stały) z siatki startów PF_BKG_GRID → najlepsze → pełny migrad. Powtarzane, póki NLL spada.
        ok = [m for _, m, _ in res if m is not None]
        if ok:
            polish(min([m for m in ok if m.valid] or ok, key=lambda m: m.fval), expo, strategy, res)
        return res

    log_att, chosen, all_runs = [], None, []
    for name, expo, strategy in attempts:
        runs = run(expo, strategy)
        ok = [(lb, m) for lb, m, _ in runs if m is not None]
        val = [(lb, m) for lb, m in ok if m.valid]
        good = [(lb, m) for lb, m in val if m.fmin.has_posdef_covar]   # valid i HESSE dodatnio określona — preferowane
        pick = min(good or val or ok, key=lambda x: x[1].fval) if ok else None
        log_att.append({"model": name, "n_starts": len(runs), "n_valid": len(val),
                        "starts": [{"start": lb, "valid": bool(m.valid) if m else False, "fval": float(m.fval) if m else None,
                                    "eps": float(m.values["eps"]) if m else None, "edm": float(m.fmin.edm) if m else None,
                                    "error": err or None} for lb, m, err in runs]})
        if pick is not None and (chosen is None or (pick[1].valid and not chosen[1].valid)):
            chosen, all_runs = (name, pick[1], expo), (val or ok)
        if chosen is not None and chosen[1].valid:
            break
    if chosen is None:
        raise RuntimeError("fit_passfail: żaden start nie dał wyniku")
    name, m, expo = chosen
    strategy = dict((nm, st) for nm, _, st in attempts)[name]
    dnll_expo = eps_expo = float("nan")
    if model_check and not expo and m.valid:   # test modelu tła: czy wykładnik (zagnieżdżony, b=0 ↔ nominał) istotnie poprawia fit
        v0 = {q: float(m.values[q]) for q in m.parameters}; best_x = None
        for b0 in (0.0, 0.05, 0.15):
            mx = make({**v0, "bp": b0, "bf": b0}, True, strategy)
            try:
                mx.migrad(ncall=50000); mx.hesse()
            except (RuntimeError, ValueError):
                continue
            if mx.valid and mx.fmin.has_posdef_covar and (best_x is None or mx.fval < best_x.fval):
                best_x = mx
        if best_x is not None:   # tło z wykładnikiem też ma lokalne minima (MC fd: mało tła) → to samo polerowanie
            px = []
            best_x = polish(best_x, True, strategy, px)
            best_x = min([best_x] + [mm for _, mm, _ in px if mm is not None and mm.valid and mm.fmin.has_posdef_covar], key=lambda mm: mm.fval)
            dnll_expo, eps_expo = float(m.fval - best_x.fval), float(best_x.values["eps"])
            log_att.append({"model": "expo_pol2_test", "dnll": dnll_expo, "eps": eps_expo, "valid": True})

    def distinct_minima(runs, m):
        """Różne minima wśród startów wybranego modelu (ten sam punkt = |Δε| < 0,1σ i ΔNLL < 0,05); [0] = wybrany fit."""
        err0 = float(m.errors["eps"]) if np.isfinite(m.errors["eps"]) and m.errors["eps"] > 0 else 1.0
        mins = []   # wybrany fit pierwszy (może nie być najniższy, gdy niższe minimum ma HESSE niedodatnią → dnll_second < 0)
        for mm in [m] + sorted([r for _, r in runs if r is not m], key=lambda r: r.fval):
            x, fv = float(mm.values["eps"]), float(mm.fval)
            if not any(abs(x - x0) < 0.1 * err0 and abs(fv - f0) < 0.05 for x0, f0 in mins):
                mins.append((x, fv))
        return mins

    prof_info, ends = {}, (float("nan"), float("nan"))
    for it in range(2):   # 2. obieg tylko gdy skan profilu znalazł niższe minimum
        if not (interval and m.valid):
            break
        f0, e0 = float(m.fval), float(m.values["eps"])
        err = float(m.errors["eps"]) if np.isfinite(m.errors["eps"]) and m.errors["eps"] > 0 else 0.1
        best_v = {p: float(m.values[p]) for p in m.parameters}; cache = {e0: (f0, best_v)}

        def prof(x):
            x = float(x)
            if x in cache:
                return cache[x]
            near = min(cache, key=lambda t: abs(t - x))
            out = None
            for st in (cache[near][1], best_v) if near != e0 else (best_v,):
                mm = make({**st, "eps": x}, expo, 1); mm.fixed["eps"] = True
                mm.migrad(ncall=20000)
                if out is None or mm.fval < out[0]:
                    out = (float(mm.fval), {p: float(mm.values[p]) for p in mm.parameters})
            cache[x] = out
            return out
        lo, dlo, xlo = _pf_profile_side(prof, e0, err, -1, f0)
        hi, dhi, xhi = _pf_profile_side(prof, e0, err, +1, f0)
        prof_info = {"n_eval": len(cache) - 1, "dnll_min_scan": float(min(dlo, dhi))}
        if min(dlo, dhi) < -0.01 and it == 0:   # profil zszedł poniżej minimum: start z tego punktu, eps wolne
            xn = xlo if dlo < dhi else xhi
            mm = make(cache[xn][1], expo, strategy); mm.migrad(ncall=50000); mm.hesse()
            prof_info["new_min"] = {"eps": float(mm.values["eps"]), "dnll": float(mm.fval - f0), "valid": bool(mm.valid)}
            if mm.valid and mm.fval < f0 - 1e-3:
                all_runs.append(("profile_new_min", mm)); m = mm
                continue
        ends = (lo, hi)
        break
    mins = distinct_minima(all_runs, m)
    fbest = float(m.fval)
    second = min(mins[1:], key=lambda t: t[1]) if len(mins) > 1 else None   # najniższe z pozostałych minimów
    v = {p: float(m.values[p]) for p in m.parameters}
    lp, lf = lams(v, expo)
    chi2 = float(np.sum((hp - lp) ** 2 / np.maximum(lp, 1)) + np.sum((hf - lf) ** 2 / np.maximum(lf, 1)))
    eps, err = v["eps"], float(m.errors["eps"])
    lo_lim, hi_lim = PF_EPS_LIMITS
    at_lim = bool(eps - lo_lim < 1e-3 or hi_lim - eps < 1e-3)
    closed = bool(np.isfinite(ends[0]) and np.isfinite(ends[1]))
    posdef = bool(m.fmin.has_posdef_covar)
    at_lims = [p for p in m.parameters if not m.fixed[p] and p != "eps" and any(
        np.isfinite(b) and abs(v[p] - b) < 1e-3 * max(1.0, abs(b)) for b in m.limits[p])]   # HESSE wtedy mniej wiarygodne niż profil
    ambig = bool(second is not None and second[1] - fbest < PF_DNLL_AMBIG and abs(second[0] - eps) > err)
    why = [w for w, bad in (("invalid", not m.valid), ("not_posdef", not posdef), ("eps_at_limit", at_lim),
                            ("interval_open", interval and not closed), ("ambiguous_minima", ambig),
                            ("low_signal", not v["Nt"] > PF_MIN_SIGNIF * float(m.errors["Nt"]))) if bad]
    return {"eps": eps, "eps_err": err, "N_tot": v["Nt"], "N_tot_err": float(m.errors["Nt"]),
            "scale": v["s"], "mu": v["mu"], "valid": bool(m.valid), "chi2": chi2, "ndf": int(2 * len(hp) - m.nfit),
            "n_pass": int(hp.sum()), "n_fail": int(hf.sum()),
            "eps_lo": ends[0], "eps_hi": ends[1], "eps_err_lo": eps - ends[0], "eps_err_hi": ends[1] - eps, "profile": prof_info,
            "posdef": posdef, "edm": float(m.fmin.edm), "nll": fbest, "bkg_model": name, "dnll_expo": dnll_expo, "eps_expo": eps_expo,
            "bkg_tension": bool(dnll_expo > PF_EXPO_DNLL), "eps_at_limit": at_lim, "params_at_limit": at_lims,
            "n_minima": len(mins), "dnll_second": float(second[1] - fbest) if second else float("nan"),
            "eps_second": float(second[0]) if second else float("nan"),
            "determined": not why, "undetermined_reason": ",".join(why), "attempts": log_att}, m

def sb_count(hp, hf, edges, mu=145.43, half=1.0, sb=(150.0, 158.0), fit_lo=MPI + 0.9):
    """Kontrola bez kształtu sygnału: N(okno) − tło interpolowane z fitu samego tła do (zakres poza oknem ±2.5 MeV)."""
    out = {}
    for nm, h in (("pass", hp), ("fail", hf)):
        h, e = _rebin(h, edges); c = 0.5 * (e[1:] + e[:-1])
        side = (np.abs(c - mu) > 2.5)
        hs = h[side]; pos = hs > 0
        def nll(nb, a, b, c1, c2):
            lam = np.clip(nb * np.diff(_bkg_cdf(e, a, b, c1, c2)), 1e-9, None)[side]
            return _nll_sat(lam, hs, pos)   # postać nasycona log1p (wartości jak dotąd, NLL O(n_bin))
        m = Minuit(nll, nb=h.sum(), a=0.5, b=0.1, c1=0.0, c2=0.0); m.errordef = 0.5
        m.limits["a"] = (0.01, 3); m.limits["b"] = (-1, 3); m.limits["c1"] = (-0.5, 2); m.limits["c2"] = (-0.05, 0.2); m.limits["nb"] = (0, None)
        m.migrad()
        lam = m.values["nb"] * np.diff(_bkg_cdf(e, *[m.values[k] for k in ("a", "b", "c1", "c2")]))
        w = np.abs(c - mu) < half
        out[nm] = float(h[w].sum() - lam[w].sum()); out[nm + "_err"] = float(np.sqrt(h[w].sum()))
    t = out["pass"] + out["fail"]
    out["eps"] = out["pass"] / t if t > 0 else float("nan")
    return out

_GL_X, _GL_W = np.polynomial.legendre.leggauss(6)


def _bkg_cdf_t(xe, a, b, c1, c2, expo=True):
    """Tło dla fit_single (NDS): u^a·e^{−bu}·(1 + c1·t + c2·t²), t = u/20 MeV → współczynniki O(1) (lepsze uwarunkowanie Hesse).
    Całki po binach kwadraturą Gaussa–Legendre'a (6 węzłów/bin; na [m_π+0,9; 160] całka jest gładka) — NLL gładkie
    do precyzji maszynowej, więc Migrad osiąga EDM < 1e-4 (siatka + np.interp w `_bkg_cdf` dawała EDM ~1e-3 → valid=False).
    expo=False: b≡0 (model „próg × pol2”, zagnieżdżony w pełnym). Zwraca CDF na krawędziach (0 … 1)."""
    xe = np.asarray(xe, float); lo, hi = xe[:-1], xe[1:]
    x = 0.5 * (hi + lo)[:, None] + 0.5 * (hi - lo)[:, None] * _GL_X[None, :]
    u = np.clip(x - MPI, 1e-9, None); t = u / 20.0
    fb = u ** a * (1 + c1 * t + c2 * t * t)
    if expo:
        fb = fb * np.exp(-b * u)
    fb = np.clip(fb, 1e-12, None)
    ib = 0.5 * (hi - lo) * (fb * _GL_W[None, :]).sum(1)
    cb = np.concatenate([[0.0], np.cumsum(ib)])
    return cb / cb[-1]


# Kolejność prób fit_single (NDS, G9): nominał „próg × pol2” (u^a(1+c1t+c2t²)) → pełny model z wykładnikiem → pełny, strategia 1.
# Pełny model (u^a e^{-bu} × pol2, nominał V) ma na 20 MeV zakresu płaski kierunek b ↔ (c1, c2): Migrad znajduje minimum
# (N stabilne), ale EDM nie spada poniżej progu / Hesse niedodatnio określona → valid=False na pełnej lumi (S/B ≈ 1 %).
# Model bez wykładnika jest w nim zagnieżdżony (b=0), daje ten sam χ² i N (NDS: 4443 vs 4437, Δχ² < 0,1) i zbiega → nominał.
# Czysty wykładnik (bez pol2) jest za sztywny: zawyża σ (×1,17) i N (+25 %) — ten sam artefakt co w massfit.fit_dm → nie w kolejce.
FIT_SINGLE_ATTEMPTS = (("threshold_pol2", False, 2), ("expo_pol2", True, 2), ("expo_pol2_s1", True, 1))


def fit_single(h, edges, shape, attempts=FIT_SINGLE_ATTEMPTS):
    """Jedna próbka (np. N_D* do normalizacji par): 2 Gaussy (kształt MC, skala wolna) + tło u^a(1+c1t+c2t²)[·e^{-bu}], t=u/20 MeV.
    Odporność (NDS): simplex→migrad (strategy 2)→hesse; gdy fit niepoprawny — kolejny model z `attempts` (zob. FIT_SINGLE_ATTEMPTS).
    Zwraca N, N_err, scale, mu, valid, chi2, ndf, n_entries + bkg_model (użyty model), edm, posdef, attempts (lista prób)."""
    h, e = _rebin(h, edges)
    def cdf_s(xe, mu, s):
        c = _sig_cdf(xe, mu, shape["sc"] * s, shape["k"], shape["f"]); return (c - c[0]) / (c[-1] - c[0])
    hs = np.where(h > 0, h, 1.0)   # NLL względem modelu nasyconego (stała addytywna): wartości O(n_bin) zamiast O(10⁶) —
    def run(expo, strategy):       # bez tego różnice NLL ~1e-4 giną w zaokrągleniach i EDM nie spada poniżej progu (valid=False)
        def nll(N, mu, s, nb, a, b, c1, c2):
            lam = np.clip(np.diff(N * cdf_s(e, mu, s) + nb * _bkg_cdf_t(e, a, b, c1, c2, expo)), 1e-9, None)
            return float(np.sum(lam - h - h * np.log(lam / hs)))
        best = None
        for a0, b0 in ((0.5, 0.03), (0.6, 0.1), (0.4, 0.0)):
            m = Minuit(nll, N=max(0.01 * h.sum(), 50), mu=shape["mu"], s=1.0, nb=h.sum(), a=a0, b=b0 if expo else 0.0, c1=0.0, c2=0.0)
            m.errordef = 0.5; m.strategy = strategy
            m.limits["N"] = (0, None); m.limits["s"] = (0.8, 1.6); m.limits["mu"] = (shape["mu"] - 0.5, shape["mu"] + 0.5)
            m.limits["a"] = (0.01, 3); m.limits["b"] = (-0.5, 1.5); m.limits["c1"] = (-5, 20); m.limits["c2"] = (-10, 40); m.limits["nb"] = (0, None)
            if not expo:
                m.fixed["b"] = True
            m.simplex(); m.migrad(ncall=50000); m.hesse()
            if best is None or (m.valid and (not best.valid or m.fval < best.fval)):
                best = m
            if best.valid:
                break
        return best
    log_att, m, used = [], None, None
    for name, expo, strategy in attempts:
        mm = run(expo, strategy)
        log_att.append({"model": name, "valid": bool(mm.valid), "N": float(mm.values["N"]), "edm": float(mm.fmin.edm)})
        if m is None or (mm.valid and not m.valid):
            m, used = mm, (name, expo)
        if m.valid:
            break
    v = {k: float(m.values[k]) for k in m.parameters}
    lam = np.diff(v["N"] * cdf_s(e, v["mu"], v["s"]) + v["nb"] * _bkg_cdf_t(e, v["a"], v["b"], v["c1"], v["c2"], used[1]))
    return {"N": v["N"], "N_err": float(m.errors["N"]), "scale": v["s"], "mu": v["mu"], "valid": bool(m.valid),
            "chi2": float(np.sum((h - lam) ** 2 / np.maximum(lam, 1))), "ndf": int(len(h) - m.nfit), "n_entries": int(h.sum()),
            "bkg_model": used[0], "edm": float(m.fmin.edm), "posdef": bool(m.fmin.has_posdef_covar), "attempts": log_att}

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
