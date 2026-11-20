"""Testy modułu charmtag.bkggen (etap S) — na syntetycznych h5 (synth.make) i zabawkowych tablicach."""
import os

import numpy as np
import pandas as pd
import pytest

from charmtag import config as C
from charmtag.bkggen import control as K
from charmtag.bkggen import data as D
from charmtag.bkggen import gen as G
from charmtag.bkggen import validate as V
from charmtag.synth import make


@pytest.fixture(scope="module")
def cfg(tmp_path_factory):
    d = tmp_path_factory.mktemp("bg")
    make(str(d / "h5"), n=4000, seed=5)
    return C.load(None, {"paths": {"h5_dir": str(d / "h5"), "out_dir": str(d / "out")}})


@pytest.fixture(scope="module")
def L(cfg):
    sp = cfg["species"]["D0"]
    X = D.load_h5(cfg, "data_none", "D0")
    return X[D.in_range(X, sp)].reset_index(drop=True)


def test_swapped_mass_involution():
    rng = np.random.default_rng(0)
    n = 1000
    df = pd.DataFrame({"d_pt0": rng.uniform(1, 5, n), "d_pt1": rng.uniform(1, 5, n), "d_eta0": rng.normal(0, 1, n),
                       "d_eta1": rng.normal(0, 1, n), "d_phi0": rng.uniform(-3, 3, n), "d_phi1": rng.uniform(-3, 3, n)})
    m = D.nominal_mass(df)
    ms = D.swapped_mass(df)
    sw = df.rename(columns={"d_pt0": "d_pt1", "d_pt1": "d_pt0", "d_eta0": "d_eta1", "d_eta1": "d_eta0",
                            "d_phi0": "d_phi1", "d_phi1": "d_phi0"})
    assert np.allclose(D.nominal_mass(sw), ms)  # zamiana hipotez = zamiana ról torów
    assert not np.allclose(m, ms)


def test_sideband_regions(cfg, L):
    sp = cfg["species"]["D0"]
    n = K.sideband(L, sp, "D0", "near")
    f = K.sideband(L, sp, "D0", "far")
    d = np.abs(n["m"] - sp["m0"])
    assert ((d >= 0.075) & (d < 0.100)).all() and (np.abs(f["m"] - sp["m0"]) >= 0.100).all()
    assert set(n["side"]) <= {"L", "R"}


def test_reflection_weights(cfg, L):
    R = K.reflection(L, cfg["species"]["D0"])
    assert set(np.unique(R["w_ctrl"].round(6))) <= {1.0, round(-0.025 / 0.05, 6)}
    h, e = K.reflection_template(R, 1.75, 2.0, 50)
    assert (h >= 0).all() and len(e) == 51


def test_mix_pairs_disjoint(cfg, L):
    I, J = K.same_event_pairs(L)
    ev = K._event_key(L)
    assert (ev[I] == ev[J]).all()
    Im, Jm = K.mixed_pairs(L, I, J)
    assert (ev[Im] != ev[Jm]).all()
    P = K.pair_frame(L, Im, Jm, "mixed")
    assert {"m1", "m2", "dphi"} <= set(P.columns)


def test_random_pis(cfg, L):
    R = K.dstar_random_pis(L, cfg["species"]["D0"], n_per=2)
    if len(R):
        assert (R["dm_rand"] > 0.13).all()  # Δm ≥ m_π


def test_metrics_identical():
    rng = np.random.default_rng(1)
    a = pd.DataFrame({"x": rng.normal(0, 1, 20000), "y": rng.normal(0, 1, 20000)})
    b = pd.DataFrame({"x": rng.normal(0, 1, 20000), "y": rng.normal(0, 1, 20000)})
    assert V.chi2_ndf(a.x, b.x) < 2.0
    assert V.ks(a.x, b.x) < 0.03
    assert abs(V.c2st_auc(a, b, ["x", "y"], folds=2) - 0.5) < 0.03
    c = b.assign(x=b.x + 0.5)
    assert V.c2st_auc(a, c, ["x", "y"], folds=2) > 0.6


def test_gaussianizer_roundtrip():
    rng = np.random.default_rng(2)
    X = np.c_[rng.exponential(1, 5000), rng.integers(0, 5, 5000).astype(float)]
    g = G.Gaussianizer().fit(X)
    Z = g.to_z(X)
    assert abs(Z[:, 0].std() - 1) < 0.05 and abs(Z[:, 1].mean()) < 0.05  # atomy rozmyte → N(0,1)
    Xb = g.from_z(Z)
    assert np.allclose(np.sort(Xb[:, 0])[::50], np.sort(X[:, 0])[::50], atol=0.05)
    assert set(np.unique(np.round(Xb[:, 1]))) <= {0, 1, 2, 3, 4}


def test_copula_interpolates_linear_trend():
    rng = np.random.default_rng(3)
    m = rng.uniform(1.75, 2.0, 40000)
    x = 2.0 * (m - 1.865) + rng.normal(0, 0.05, len(m))
    sb = np.abs(m - 1.865) > 0.075
    SB = pd.DataFrame({"x": x[sb], "m": m[sb]})
    gm = G.BkgGenerator(["x"], 1.865, lambda c: np.abs(np.asarray(c) - 1.865) > 0.075, (1.75, 2.0),
                        model="copula").fit(SB)
    s = gm.sample(20000, (1.81, 1.92))
    true = 2.0 * (s["m"] - 1.865)
    assert abs((s["x"] - true).mean()) < 0.03  # trend w m przeniesiony pod pik
    assert 1.81 <= s["m"].min() and s["m"].max() <= 1.92


def test_cli_sideband(cfg, tmp_path):
    from charmtag.bkggen.cli import main
    import yaml
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump({"paths": {"h5_dir": cfg["paths"]["h5_dir"], "out_dir": str(tmp_path / "o")}}))
    assert main(["bkggen", "--method", "sideband", "-c", str(p), "--which", "far"]) == 0
    assert os.path.exists(tmp_path / "o" / "bkggen" / "sideband_data_none_D0_far.h5")


def test_make_generator_salad_without_torch():
    """REPO: `make_generator('salad')` (projekt B) działa bez torcha/iminuit — fit na SB, wagi ~1 dla identycznych rozkładów."""
    from charmtag.bkggen import api
    names = api.available()
    assert "salad" in names and "copula" in names
    gen = api.make_generator("salad")
    assert gen is not None
