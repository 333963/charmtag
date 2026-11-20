"""`dstar_refit.fit_passfail` (2026-11-12): jednoczesny fit Δm pass/fail → ε, N_tot.
Toy pass/fail (tło z wykładnikiem, b ≠ 0, czyli spoza nominalnego modelu „próg × pol2”): zbieżność, ε w 3σ, przedział z profilu,
stabilność na perturbacje nieistotne fizycznie, flaga nieokreśloności bez sygnału. Szybkie (~30 s łącznie)."""
import numpy as np
import pytest
from charmtag import dstar_refit as R

SHAPE = {"mu": 145.5456, "sc": 0.7429, "k": 2.7286, "f": 0.7183}   # kształt MC truth D* (V dstar_post.json)
EDGES = np.arange(R.LO, R.HI + 1e-9, 0.1)
OLD_KEYS = {"eps", "eps_err", "N_tot", "N_tot_err", "scale", "mu", "valid", "chi2", "ndf", "n_pass", "n_fail"}


def _lam(nt, eps, nbp, nbf, b=0.03, ap=0.55, af=0.45, c1=0.2, edges=EDGES):
    c = R._sig_cdf(edges, SHAPE["mu"], SHAPE["sc"], SHAPE["k"], SHAPE["f"]); cs = np.diff(c) / (c[-1] - c[0])
    bp = nbp * np.diff(R._bkg_cdf_t(edges, ap, b, c1, 0.0, True))
    bf = nbf * np.diff(R._bkg_cdf_t(edges, af, b, -0.1, 0.05, True))
    return nt * eps * cs + bp, nt * (1 - eps) * cs + bf


def _toy(seed, **kw):
    rng = np.random.default_rng(seed)
    lp, lf = _lam(**kw)
    return rng.poisson(lp).astype(float), rng.poisson(lf).astype(float)


@pytest.fixture(scope="module")
def toy_fit():
    hp, hf = _toy(11, nt=3000, eps=0.7, nbp=20000, nbf=100000)
    return hp, hf, R.fit_passfail(hp, hf, EDGES, SHAPE)


def test_passfail_toy_converges_and_covers(toy_fit):
    """Toy jak V 5-8 GeV (S/B pod pikiem ~0,3/0,05): valid, posdef, określony, ε w 3σ, przedział z profilu zawiera ~±1σ HESSE."""
    _, _, f = toy_fit
    assert OLD_KEYS <= set(f)                          # API wstecznie zgodne (pipeline.dstar_control, t4_fit, n2_fit)
    assert f["valid"] and f["posdef"] and f["determined"], (f["undetermined_reason"], f["attempts"])
    assert f["bkg_model"] == "threshold_pol2"
    assert abs(f["eps"] - 0.7) < 3 * f["eps_err"], (f["eps"], f["eps_err"])
    assert f["eps_lo"] < f["eps"] < f["eps_hi"]
    assert 0.6 < (f["eps_hi"] - f["eps_lo"]) / (2 * f["eps_err"]) < 1.5
    assert f["n_minima"] >= 1 and not f["eps_at_limit"]
    assert abs(f["N_tot"] - 3000) < 3 * f["N_tot_err"]


def test_passfail_stable_under_irrelevant_perturbations(toy_fit, monkeypatch):
    """Krawędzie ±1e-9 MeV i starty ×(1 + 1e-6·N(0,1)): |Δε| < 2 % σ, |Δσ|/σ < 5 %, ten sam stan flag."""
    hp, hf, f0 = toy_fit
    runs = [R.fit_passfail(hp, hf, EDGES + d, SHAPE, interval=False) for d in (1e-9, -1e-9)]
    orig = R._pf_starts
    rng = np.random.default_rng(5)
    monkeypatch.setattr(R, "_pf_starts", lambda *a, **k: [(lb, {p: v * (1 + 1e-6 * rng.normal()) for p, v in st.items()})
                                                           for lb, st in orig(*a, **k)])
    runs.append(R.fit_passfail(hp, hf, EDGES, SHAPE, interval=False))
    for r in runs:
        assert r["valid"] and r["posdef"]
        assert abs(r["eps"] - f0["eps"]) < 0.02 * f0["eps_err"], (r["eps"], f0["eps"], f0["eps_err"])
        assert abs(r["eps_err"] / f0["eps_err"] - 1) < 0.05


def test_passfail_no_signal_flagged_undetermined():
    """Bez sygnału (N_tot = 0): ε nie jest wyznaczone przez dane → determined=False z podanym powodem (nie „naprawiane” na siłę)."""
    hp, hf = _toy(3, nt=0, eps=0.7, nbp=20000, nbf=100000)
    f = R.fit_passfail(hp, hf, EDGES, SHAPE)
    assert OLD_KEYS <= set(f)
    assert not f["determined"] and f["undetermined_reason"], f
    assert f["N_tot"] < 3 * max(f["N_tot_err"], 1.0)


def test_nll_sat_matches_poisson_and_is_exact_at_large_counts():
    """Postać log1p = Σ[λ − n − n ln(λ/n)] (n = 0: λ) i nie traci precyzji przy 10⁶ wpisów/bin."""
    n = np.array([0.0, 1.0, 5.0, 1e6, 1e6]); lam = np.array([0.3, 1.2, 4.0, 1e6 + 1.0, 1e6 * (1 + 1e-9)])
    pos = n > 0
    ref = [0.3, 1.2 - 1 - np.log(1.2), 4.0 - 5 - 5 * np.log(0.8), None, None]
    for i in range(3):
        assert np.isclose(R._nll_sat(lam[i:i + 1], n[i:i + 1], pos[i:i + 1]), ref[i], rtol=1e-12)
    assert np.isclose(R._nll_sat(lam[3:4], n[3:4], pos[3:4]), 0.5 / 1e6, rtol=1e-5)          # (λ−n)²/(2n)
    assert np.isclose(R._nll_sat(lam[4:5], n[4:5], pos[4:5]), 0.5 * 1e6 * 1e-18, rtol=1e-4)


def test_passfail_model_check_flags_strong_exponent():
    """Tło z silnym wykładnikiem (b = 0,15/MeV, wysoka statystyka): test modelu zgłasza bkg_tension (ΔNLL > PF_EXPO_DNLL);
    ε z modelu z wykładnikiem (eps_expo) w 3σ prawdy. Wynik nominalny bez zmian (próg × pol2)."""
    hp, hf = _toy(21, nt=3000, eps=0.7, nbp=50000, nbf=300000, b=0.15)
    f = R.fit_passfail(hp, hf, EDGES, SHAPE, interval=False)
    assert f["valid"] and f["bkg_tension"] and f["dnll_expo"] > R.PF_EXPO_DNLL, (f["dnll_expo"], f["attempts"])
    assert f["bkg_model"] == "threshold_pol2"
    assert abs(f["eps_expo"] - 0.7) < 3 * f["eps_err"], (f["eps_expo"], f["eps_err"])
