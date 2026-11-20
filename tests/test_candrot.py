"""S6: testy candrot (próbki kontrolne na poziomie kandydata) i kin3 (Dalitz D⁺)."""
import numpy as np
import pandas as pd
import pytest

from charmtag import config as C
from charmtag import features as F
from charmtag.bkggen import candrot as R
from charmtag.bkggen import kin3


@pytest.fixture(scope="module")
def cfg():
    return C.load(None)


def toy(nd, n=4000, seed=0, m_target=None):
    """Ramka logiczna charmtag (GeV) z losowymi córkami; zmienne wierzchołka losowe."""
    rng = np.random.default_rng(seed)
    d = {k: rng.normal(size=n) for k in F.SCALARS}
    d["run"] = np.full(n, 1)
    d["event"] = rng.integers(0, n // 3, n)
    d["weight"] = np.ones(n)
    d["truth"] = rng.choice([0, 8, 12], n)
    d["Lxy"] = rng.exponential(1, n) + 0.25
    d["cos_theta_p"] = rng.uniform(0.9, 1, n)
    for i in range(nd):
        d[f"d_pt{i}"] = rng.uniform(0.5, 6, n)
        d[f"d_eta{i}"] = rng.uniform(-2, 2, n)
        d[f"d_phi{i}"] = rng.uniform(-np.pi, np.pi, n)
        d[f"d_p{i}"] = np.nan
        d[f"d_d0{i}"] = rng.normal(size=n)
        d[f"d_z0{i}"] = rng.normal(size=n)
        d[f"d_dedx{i}"] = rng.uniform(0.5, 2, n)
        d[f"d_npts{i}"] = rng.integers(0, 5, n).astype(float)
    return pd.DataFrame(d)


def sp_wide(cfg, name):
    sp = dict(cfg["species"][name])
    sp["mass_range"] = [0.0, 100.0]
    return sp


def test_kin3_roundtrip():
    rng = np.random.default_rng(1)
    n = 5000
    pts = [rng.uniform(0.3, 8, n) for _ in range(3)]
    etas = [rng.uniform(-2.5, 2.5, n) for _ in range(3)]
    phis = [rng.uniform(-np.pi, np.pi, n) for _ in range(3)]
    a = kin3.angles3(pts, etas, phis)
    b = kin3.daughters3(a["m"], a["pt"], a["eta"], a["phi"], a["s12"], a["s13"], a["cts"], a["phis"], a["psi"])
    for i in range(3):
        assert np.allclose(b[0][i], pts[i], atol=1e-8)
        assert np.allclose(b[1][i], etas[i], atol=1e-8)
        assert np.allclose(np.angle(np.exp(1j * (b[2][i] - phis[i]))), 0, atol=1e-8)
    u, v = kin3.dalitz_norm(a["m"], a["s12"], a["s13"])
    assert (u >= 0).all() and (u <= 1).all() and (v >= 0).all() and (v <= 1).all()


def test_recompute_cts_matches_kin3(cfg):
    L = toy(3)
    sp = sp_wide(cfg, "Dplus")
    R.recompute(L, sp, "Dplus")
    a = kin3.angles3([L[f"d_pt{i}"].to_numpy() for i in range(3)], [L[f"d_eta{i}"].to_numpy() for i in range(3)],
                     [L[f"d_phi{i}"].to_numpy() for i in range(3)])
    assert np.allclose(L["m"], a["m"], atol=1e-9)
    assert np.allclose(L["cos_theta_star"], a["cts"], atol=1e-9)
    assert np.allclose(L["s12"], a["s12"], atol=1e-9)


@pytest.mark.parametrize("name,nd", [("D0", 2), ("Dplus", 3)])
def test_rotation_identities(cfg, name, nd):
    L = toy(nd)
    sp = sp_wide(cfg, name)
    R.recompute(L, sp, name)
    X = R.make(L, sp, name, "rot_pi", daughter=0, presel=False)
    assert np.allclose(X["d_pt0"], L["d_pt0"]) and np.allclose(X["d_eta0"], L["d_eta0"])
    assert np.allclose(np.cos(X["d_phi0"] - L["d_phi0"]), -1)
    assert not np.allclose(X["m"], L["m"])
    assert np.allclose(X["Lxy"], L["Lxy"])  # wierzchołek: keep
    Y = R.make(X, sp, name, "rot_pi", daughter=0, presel=False)
    assert np.allclose(Y["m"], L["m"], atol=1e-9)
    Z = R.make(R.make(L, sp, name, "flip_eta", daughter=1, presel=False), sp, name, "flip_eta", daughter=1,
               presel=False)
    assert np.allclose(Z["m"], L["m"], atol=1e-9)


def test_rot_rand_angle_range(cfg):
    L = toy(2)
    sp = sp_wide(cfg, "D0")
    X = R.make(L, sp, "D0", "rot_rand", daughter=1, presel=False, n_rep=2)
    dphi = np.abs(np.angle(np.exp(1j * (X["d_phi1"].to_numpy() - np.r_[L["d_phi1"], L["d_phi1"]]))))
    assert (dphi >= 0.5 * np.pi - 1e-9).all()
    assert np.isclose(X["w_ctrl"].sum(), len(L))


def test_swap_d0_equals_reflected_mass(cfg):
    L = toy(2)
    sp = sp_wide(cfg, "D0")
    R.recompute(L, sp, "D0")
    X = R.make(L, sp, "D0", "swap", presel=False)
    assert np.allclose(X["m"], L["m_swap"], atol=1e-9)
    assert np.allclose(X["d_dedx0"], L["d_dedx1"])  # zmienne torowe wędrują z torem


def test_swap_dplus_pion_choice(cfg):
    L = toy(3)
    sp = sp_wide(cfg, "Dplus")
    X = R.make(L, sp, "Dplus", "swap", presel=False, swap_pion=2)
    assert np.allclose(X["d_pt0"], L["d_pt2"]) and np.allclose(X["d_pt2"], L["d_pt0"])
    assert np.allclose(X["d_pt1"], L["d_pt1"])


def test_mix_other_event(cfg):
    L = toy(2, n=6000)
    sp = sp_wide(cfg, "D0")
    pool = R.track_pool(L, 2)
    X = R.make(L, sp, "D0", "mix", daughter=1, presel=False, pool=pool, seed=3)
    assert len(X) > 0.95 * len(L)
    # tor podmieniony pochodzi z innego zdarzenia, kolumny toru spójne (pT ↔ dE/dx tego samego toru)
    key = dict(zip(np.round(pool["d_pt"], 10), np.round(pool["d_dedx"], 6)))
    assert all(key[a] == b for a, b in zip(np.round(X["d_pt1"], 10), np.round(X["d_dedx1"], 6)))
    ev_of = dict(zip(np.round(pool["d_pt"], 10), pool["ev"]))
    evX = R._event_key(X)
    assert np.mean([ev_of[a] != e for a, e in zip(np.round(X["d_pt1"], 10), evX)]) == 1.0
    # pT podmienionego toru blisko oryginału (ranking ±50 w puli)
    src = L.iloc[X["src_row"].to_numpy()]
    assert np.median(np.abs(np.log(X["d_pt1"].to_numpy() / src["d_pt1"].to_numpy()))) < 0.1


def test_presel_and_vertex_resample(cfg):
    L = toy(2, n=8000)
    sp = dict(cfg["species"]["D0"])
    R.recompute(L, sp, "D0")
    ref = toy(2, n=3000, seed=5)
    R.recompute(ref, sp_wide(cfg, "D0"), "D0")
    ref["Lxy"] = 1000 + np.arange(len(ref))
    X = R.make(L, sp_wide(cfg, "D0"), "D0", "rot_rand", presel=False, vertex_mode="resample", vertex_ref=ref)
    assert (X["Lxy"] >= 1000).all()
    Y = R.make(L, sp, "D0", "rot_rand", n_rep=3)
    lo, hi = sp["mass_range"]
    assert ((Y["m"] >= lo) & (Y["m"] < hi) & (Y["pt"] >= 3)).all()
    with pytest.raises(ValueError):
        R.make(L, sp, "D0", "nope")
