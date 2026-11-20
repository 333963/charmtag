"""NDS (G9): `dstar_refit.fit_single` musi zbiegać na pełnej lumi (valid=True, Hesse OK) i odtwarzać N na toyach.
Test na prawdziwych histogramach uruchamia się tylko gdy CHARM_DSTAR_NPZ wskazuje npz z V (dane nie są w repo)."""
import os
import numpy as np
import pytest
from charmtag import dstar_refit as R

SHAPE = {"mu": 145.5456, "sc": 0.7429, "k": 2.7286, "f": 0.7183}   # kształt MC truth D* (V dstar_post.json)


def _toy(n_sig, n_bkg, seed, a=0.55, b=0.03, c1=0.0, c2=0.0, edges=None):
    rng = np.random.default_rng(seed)
    edges = np.arange(R.LO, R.HI + 1e-9, 0.1) if edges is None else edges
    sig = np.where(rng.random(n_sig) < SHAPE["f"], rng.normal(SHAPE["mu"], SHAPE["sc"], n_sig), rng.normal(SHAPE["mu"], SHAPE["sc"] * SHAPE["k"], n_sig))
    x = np.linspace(edges[0], edges[-1], 20001); u = x - R.MPI; t = u / 20
    pdf = u ** a * np.exp(-b * u) * (1 + c1 * t + c2 * t * t); pdf /= pdf.sum()
    bkg = rng.choice(x, size=n_bkg, p=pdf) + rng.uniform(-2.5e-4, 2.5e-4, n_bkg)
    h, _ = np.histogram(np.r_[sig, bkg], edges)
    return h, edges


def test_fit_single_toy_full_lumi_like_converges_and_is_unbiased():
    """S/B ≈ 1 % jak na pełnej lumi (4,4k sygnału na 480k tła): valid, Hesse dodatnio określona, pull < 3."""
    h, e = _toy(4400, 480000, seed=1)
    f = R.fit_single(h, e, SHAPE)
    assert f["valid"] and f["posdef"], f["attempts"]
    assert abs(f["N"] - 4400) < 3 * f["N_err"], (f["N"], f["N_err"])
    assert 0.9 < f["scale"] < 1.15
    assert f["bkg_model"] in ("expo_pol2", "threshold_pol2")
    assert set(f) >= {"N", "N_err", "scale", "mu", "valid", "chi2", "ndf", "n_entries"}   # API bez zmian (pairs/build.py)


def test_fit_single_toy_ensemble_pull():
    """10 toyów: średni pull N w (−1, 1), wszystkie fity poprawne."""
    pulls, ok = [], []
    for s in range(10):
        h, e = _toy(4400, 480000, seed=100 + s, c1=0.3, c2=-0.05)
        f = R.fit_single(h, e, SHAPE)
        ok.append(f["valid"]); pulls.append((f["N"] - 4400) / f["N_err"])
    assert all(ok)
    assert abs(np.mean(pulls)) < 1.0 and np.std(pulls) < 1.6, pulls


def test_fit_single_small_sample_still_works():
    """Mała próbka (jak testy pairs na MC): zbiega, N w 3σ."""
    h, e = _toy(300, 5000, seed=7)
    f = R.fit_single(h, e, SHAPE)
    assert f["valid"] and abs(f["N"] - 300) < 3 * f["N_err"]


@pytest.mark.skipif(not os.environ.get("CHARM_DSTAR_NPZ"), reason="pełna lumi: ustaw CHARM_DSTAR_NPZ=.../dstar_nom.npz")
def test_fit_single_full_lumi_real_histograms():
    """Prawdziwe histogramy V (norm_all = selekcja normalizacji par): valid=True, N zgodne z NDS (4,4k ± 0,6k), per bin też zbiega."""
    z = np.load(os.environ["CHARM_DSTAR_NPZ"]); e = z["edges"]
    f = R.fit_single(z["norm_all"], e, SHAPE)
    assert f["valid"] and f["posdef"], f["attempts"]
    assert 3500 < f["N"] < 5300 and f["N_err"] < 800, (f["N"], f["N_err"])
    assert f["scale"] < 1.1                      # bez artefaktu σ (fit_dm: σ 1,12 MeV ≈ ×1,5)
    for b in ("3-5", "5-8", "8-30"):
        fb = R.fit_single(z[f"norm_{b}"], e, SHAPE)
        assert fb["valid"], (b, fb["attempts"])
