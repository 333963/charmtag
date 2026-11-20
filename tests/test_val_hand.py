"""VAL: charmtag.hand (kod W) — score na kanonicznej ramce i adapter z surowych kolumn h5 v1.1."""
import numpy as np
from charmtag import hand
from charmtag.hand import adapt


def _raw(n=500, nd=3, seed=0):
    rng = np.random.default_rng(seed)
    g = {"Dmeson_mass": rng.uniform(1.7, 2.2, n), "Dmeson_pt": rng.uniform(3, 30, n), "Dmeson_eta": rng.uniform(-2, 2, n),
         "Dmeson_Lxy": rng.exponential(0.5, n), "Dmeson_cosThetaP": 1 - rng.exponential(0.01, n), "Dmeson_chi2": rng.exponential(2, n),
         "Dmeson_cosThetaStar": rng.uniform(-1, 1, n), "Dmeson_impactSignificance": rng.exponential(1, n)}
    for i in range(nd):
        g[f"Dmeson_daughter{i}_pt"] = rng.uniform(500, 5000, n); g[f"Dmeson_daughter{i}_eta"] = rng.uniform(-2, 2, n)
        g[f"Dmeson_daughter{i}_d0"] = rng.normal(0, 0.05, n); g[f"Dmeson_daughter{i}_dEdx"] = rng.uniform(0.8, 1.5, n)
        g[f"Dmeson_daughter{i}_npts_dEdx"] = rng.integers(1, 5, n).astype(float)
    return g


def test_adapter_and_score_dplus():
    df = adapt.to_canonical(_raw(), "Dplus")
    s = hand.score(df, species="Dplus", variant="B", rev=1, dedx_scale=adapt.DEDX_SCALE)
    col = [c for c in s.columns if c.startswith("score_")][0]
    assert s[col].between(0, 1).all() and np.isfinite(s[col]).all()


def test_adapter_d0():
    df = adapt.to_canonical(_raw(nd=2), "D0")
    assert "p_pi2" not in df and np.isfinite(df["lxy_sig"]).all()
    s = hand.score(df, species="D0", variant="A", rev=1)
    assert len(s) == len(df)


def test_variant_B_is_mass_safe():
    assert "cos_star" not in hand.MASS_SAFE and "abs_eta" not in hand.MASS_SAFE and "logp_K" not in hand.MASS_SAFE
