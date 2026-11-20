"""`dstar_refit.shape_from_mc` — limit k = σ_wide/σ_core (2026-11-13, k_max domyślnie 3 zamiast 8).
Syntetyk (binowany, ustalony seed): wąski Gauss + płaski cokół w oknie 143–148,5 → bez limitu (k_max = 8) szeroki Gauss modeluje cokół
i k ucieka do granicy; z domyślnym k_max k ≤ k_max i k_at_limit = True. Czysty 2-Gauss z k = 2,5: wynik niezależny od k_max. Szybkie (< 2 s)."""
import numpy as np
import pytest
from charmtag import dstar_refit as R

E = np.arange(R.LO, R.HI + 1e-9, 0.05)
C = 0.5 * (E[1:] + E[:-1])
PARS = ("mu", "sc", "k", "f")


def _hist(seed, n, sc, k, f, pedestal=0.0, mu=145.5):
    """Poisson z: n·(2 Gaussy) + pedestal·n zdarzeń płasko w oknie 143–148,5 MeV."""
    c = R._sig_cdf(E, mu, sc, k, f)
    lam = n * np.diff(c) / (c[-1] - c[0]) + pedestal * n * ((C > 143) & (C < 148.5)) * np.diff(E) / 5.5
    return np.random.default_rng(seed).poisson(lam).astype(float)


def test_shape_k_limit_respected_and_flagged():
    h = _hist(1, 20000, 0.8, 1.5, 1.0, pedestal=0.1)       # sygnał 1 Gauss σ = 0,8 MeV + 10 % płaskiego tła w oknie
    free = R.shape_from_mc(h, E, k_max=8.0)
    assert free["k"] > 7.0 and free["k_at_limit"], free   # stary limit: k ucieka do granicy 8
    sh = R.shape_from_mc(h, E)
    assert sh["k"] <= 3.0 + 1e-9 and sh["k_at_limit"], sh
    assert set(sh) == {*PARS, "k_at_limit"} and isinstance(sh["k_at_limit"], bool)
    assert abs(sh["mu"] - 145.5) < 0.05


@pytest.mark.parametrize("k_max", [None, 4.0])
def test_shape_clean_2g_independent_of_kmax(k_max):
    h = _hist(1, 50000, 0.75, 2.5, 0.7)
    ref = R.shape_from_mc(h, E, k_max=8.0)
    sh = R.shape_from_mc(h, E) if k_max is None else R.shape_from_mc(h, E, k_max=k_max)
    assert not ref["k_at_limit"] and not sh["k_at_limit"]
    assert 2.2 < sh["k"] < 2.8
    for p in PARS:                                        # różnice tylko na poziomie zbieżności MIGRAD (~1e-4 wzgl.)
        assert abs(sh[p] / ref[p] - 1) < 2e-3, (p, sh[p], ref[p])


def test_shape_unbinned_branch_limit():
    rng = np.random.default_rng(3)
    x = np.r_[rng.normal(145.5, 0.8, 20000), rng.uniform(143, 148.5, 2000)]
    sh = R.shape_from_mc(x)
    assert sh["k"] <= 3.0 + 1e-9 and isinstance(sh["k_at_limit"], bool)


def test_shape_kmax_invalid():
    with pytest.raises(ValueError):
        R.shape_from_mc(np.ones(len(E)), E, k_max=1.0)
