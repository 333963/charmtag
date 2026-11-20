import json
import os

import numpy as np
import pandas as pd
import pytest

from charmtag import config as C
from charmtag import features as F
from charmtag import massfit as MF
from charmtag import pipeline as P
from charmtag import samples as SA
from charmtag import store
from charmtag.io import resolve_columns
from charmtag.synth import make


@pytest.fixture(scope="session")
def cfg(tmp_path_factory):
    d = tmp_path_factory.mktemp("ct")
    make(str(d / "h5"), n=4000, seed=3)
    return C.load(None, {"paths": {"h5_dir": str(d / "h5"), "out_dir": str(d / "out")},
                         "training": {"xgb": {"n_estimators": 60}}})


def test_resolve_columns():
    cm = resolve_columns(["dmeson_M", "Dmeson_daughter0_pt", "Dmeson_daughter1_pt"],
                         {"m": ["Dmeson_m"], "d_pt": ["Dmeson_daughter{i}_pt"], "x": ["nope"]}, 3)
    assert cm["m"] == "dmeson_M" and cm["d_pt1"] == "Dmeson_daughter1_pt"
    assert "d_pt2" in cm.missing and "x" in cm.missing


def test_fold_deterministic_balanced():
    ev = np.arange(100000)
    f1 = SA.fold_of(np.zeros_like(ev), ev, 5)
    f2 = SA.fold_of(np.zeros_like(ev), ev, 5)
    assert (f1 == f2).all()
    c = np.bincount(f1, minlength=5) / len(ev)
    assert np.all(np.abs(c - 0.2) < 0.01)
    assert (SA.fold_of(np.full(10, 428770), np.arange(10), 5) != SA.fold_of(np.zeros(10), np.arange(10), 5)).any()


def test_no_isolation_features(cfg):
    for v in ("H", "frozen7"):
        assert not any("ptcone" in f for f in F.feature_list(cfg, v))
    bad = json.loads(json.dumps(cfg))
    bad["features"]["extra_H"].append("ptcone")
    with pytest.raises(ValueError):
        F.feature_list(bad, "H")
    assert len(F.feature_list(cfg, "frozen7")) == 7
    assert "log_pt" in F.feature_list(cfg, "frozen7")


def test_reweight():
    rng = np.random.default_rng(0)
    S = pd.DataFrame({"pt": 3 + rng.exponential(3, 20000), "abs_eta": rng.uniform(0, 2, 20000), "weight": 1.0})
    B = pd.DataFrame({"pt": 3 + rng.exponential(2, 40000), "abs_eta": rng.uniform(0, 2, 40000), "weight": 1.0})
    w = SA.reweight(B, S, {"pt": [3, 4, 5, 6, 8, 10, 15, 30], "abs_eta": [0, 1, 2]})
    e = [3, 4, 5, 6, 8, 10, 15, 30]
    hs = np.histogram(S["pt"], e)[0] / np.histogram(S["pt"], e)[0].sum()
    hb = np.histogram(B["pt"], e, weights=w)[0] / w[B["pt"] < 30].sum()
    assert np.allclose(hs, hb, atol=0.01)
    assert abs(w.sum() - len(S)) < 1e-6 * len(S)


def _toy(rng, ns, nb, lo=1.75, hi=2.0, m0=1.865, s=0.02, refl=0):
    b = lo + rng.exponential(0.3, 3 * nb)
    b = b[b < hi][:nb]
    m = np.concatenate([rng.normal(m0, s, ns), b, rng.normal(m0, 0.07, refl)])
    return m[(m >= lo) & (m < hi)]


def test_massfit_yield():
    rng = np.random.default_rng(1)
    sh = MF.fit_mc_shape(rng.normal(1.865, 0.02, 20000), np.ones(20000), 1.75, 2.0, 1.865, "gauss")
    pulls = []
    fc = {"nbins": 100, "background_model": "exp", "sigma_mode": "constrained", "sigma_constraint_rel": 0.1}
    for i in range(20):
        r = MF.fit_data(_toy(rng, 500, 5000), 1.75, 2.0, 1.865, sh, fc)
        pulls.append((r["N_sig"] - 500) / r["N_sig_err"])
        assert r["valid"]
    assert abs(np.mean(pulls)) < 0.8 and 0.5 < np.std(pulls) < 1.6


def test_massfit_reflection():
    rng = np.random.default_rng(2)
    sh = MF.fit_mc_shape(rng.normal(1.865, 0.02, 20000), np.ones(20000), 1.75, 2.0, 1.865, "gauss2")
    tmpl, _ = np.histogram(rng.normal(1.865, 0.07, 50000), 100, (1.75, 2.0))
    fc = {"nbins": 100, "background_model": "pol2", "sigma_mode": "constrained", "sigma_constraint_rel": 0.1}
    nr = []
    nw = []
    for i in range(10):
        m = _toy(rng, 2000, 6000, refl=1000)
        nr.append(MF.fit_data(m, 1.75, 2.0, 1.865, sh, fc, tmpl, 0.5)["N_sig"])
        nw.append(MF.fit_data(m, 1.75, 2.0, 1.865, sh, fc)["N_sig"])
    assert abs(np.mean(nr) / 2000 - 1) < abs(np.mean(nw) / 2000 - 1) + 0.02
    assert abs(np.mean(nr) / 2000 - 1) < 0.06


@pytest.fixture(scope="session")
def trained(cfg):
    out = {}
    for sp in ("D0", "Dplus"):
        for v in ("H", "frozen7"):
            out[(sp, v)] = P.cmd_train(cfg, sp, v)
            P.cmd_apply(cfg, sp, v)
            P.cmd_fit(cfg, sp, v)
        P.cmd_validate(cfg, sp)
    return out


def test_train_metrics(trained):
    for (sp, v), (meta, met) in trained.items():
        assert meta["n_folds"] == 5 and len(meta["folds"]) == 5
        assert not any("ptcone" in f for f in meta["features"])
        for key, m in met.items():
            assert m["auc"] > 0.85, (sp, v, key, m["auc"])
            assert 0.75 < m["effS_80"] < 0.85
            assert m["effB_50"] <= m["effB_70"] <= m["effB_80"]


def test_apply_equals_oof(cfg, trained):
    """Score z `apply` dla próbki treningowej = wynik out-of-fold z treningu (ten sam fold-hash)."""
    oof = store.load(os.path.join(P.mdir(cfg, "D0", "H"), "oof.h5"))
    sc = P.load_scores(cfg, "D0", "H", "signal_mc")
    s_oof = oof[oof["y"] == 1].sort_values(["event", "m"])["score"].to_numpy()
    tm = sc[sc["truth"].isin([0, 1])]
    s_app = tm.sort_values(["event", "m"])["score"].to_numpy()
    assert len(s_oof) == len(s_app)
    assert np.allclose(np.sort(s_oof), np.sort(s_app), atol=1e-5)


def test_apply_row_index(cfg, trained):
    import h5py
    sc = P.load_scores(cfg, "D0", "H", "data_none")
    with h5py.File(C.sample_path(cfg, "data_none")) as f:
        ev = f["events"].fields(["Dmeson_mass"])[:]["Dmeson_mass"]
    assert np.allclose(ev[sc["row"].to_numpy()], sc["m"].to_numpy())


def test_fit_output(cfg, trained):
    fr = json.load(open(os.path.join(cfg["paths"]["out_dir"], "fits", "D0_H.json")))
    e = fr["bins"]["3-5"]["data_none_wp70"]
    assert e["N_sig"] > 0 and e["N_sig_err"] > 0 and e["refl_ratio"] > 0


def test_validation_and_report(cfg, trained):
    from charmtag import report as R
    v = json.load(open(os.path.join(cfg["paths"]["out_dir"], "validation", "D0.json")))
    d = v["variants"]["H"]["dstar_data"]["3-5"]
    assert d["fit_all"]["N"] > 50 and np.isfinite(d["eff_70"])
    assert abs(d["eff_70"] - 0.7) < 0.15
    assert "H_vs_frozen7" in v
    md, pngs = R.make(cfg, "D0")
    assert "Fit masy" in md and all(os.path.exists(p) for p in pngs)
    md, _ = R.make(cfg, "Dplus")
    assert "Dplus" in md


def test_dplus_veto(cfg):
    L = SA.load_species(cfg, "signal_mc", "Dplus", cache=False)
    keep = SA.dstar_veto_mask(L, cfg["species"]["Dplus"])
    assert keep.mean() > 0.9


def test_cli_inspect(cfg, capsys):
    from charmtag.cli import main
    assert main(["inspect", "--h5-dir", cfg["paths"]["h5_dir"], "--out-dir", cfg["paths"]["out_dir"]]) == 0
    rep = json.load(open(os.path.join(cfg["paths"]["out_dir"], "inspect.json")))
    assert rep["data_none"]["missing"] == [] or "truth" not in rep["signal_mc"]["missing"]


def test_mark_swaps():
    sp = {"m0": 1.86484}
    cfg_ = {"signal": {"truth_values": [0, 1], "swap_resolution": "closest_mass"}}
    L = pd.DataFrame({"run": 0, "event": [1, 1, 2, 3, 3], "d_pt0": [1.0, 2.0, 1.0, 5.0, 6.0],
                      "d_pt1": [2.0, 1.0, 3.0, 6.0, 5.0], "m": [1.80, 1.87, 1.9, 1.86, 1.95],
                      "truth": [0.0, 0.0, 1.0, 12.0, 12.0]})
    SA.mark_swaps(L, sp, cfg_)
    assert L["truth"].tolist() == [100.0, 0.0, 1.0, 12.0, 12.0]


def test_fit_dm():
    rng = np.random.default_rng(5)
    x0 = 0.13957
    b = x0 + rng.gamma(1.5, 0.008, 60000)
    s = rng.normal(0.14543, 0.0007, 3000)
    r = MF.fit_dm(np.concatenate([s, b]))
    assert abs(r["N"] - 3000) < 3 * r["N_err"] + 100
