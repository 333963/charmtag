"""Testy modułu par (etap O): model 2D, ślepota, arbitraż, normalizacja, domknięcie z 09-27 (jeśli są wejścia)."""
import json
import os

import numpy as np
import pandas as pd
import pytest

from charmtag import config as C
from charmtag import pipeline as P
from charmtag.pairs import CHANNELS
from charmtag.pairs.model import Inputs, Model, M0, region_mask
from charmtag.pairs.hists import build_pair_hists
from charmtag.pairs import build as B
from charmtag.pairs import fit2d as F

FREEZE = os.path.expanduser(os.environ.get("CHARM_FREEZE", ""))    # zamrożone wejścia (niepubliczne) — testy closure skipują bez nich
STAGEB = os.path.expanduser(os.environ.get("CHARM_STAGEB", ""))
MB = {"dstar_true_per_event": {"unw": 6.67e-05, "w": 4.51e-04}, "D0_correct_per_event": {"unw": 4.79e-4, "w": 3.03e-3},
      "Dplus_correct_per_event": {"unw": 4.40e-4, "w": 2.36e-3},
      "pairs_D0D0": {"per_event_unw": 1.07e-6, "per_event_w": 9.9e-6}, "pairs_DpD0": {"per_event_unw": 2.08e-6, "per_event_w": 1.69e-5},
      "pairs_DpDp": {"per_event_unw": 7.6e-7, "per_event_w": 5.05e-6}}


# ----------------------------------------------------------------------------- toy wejścia
def toy_inputs(seed=1, nleg=40000, npair=400.0):
    rng = np.random.default_rng(seed)
    E = {10: np.linspace(1750, 2000, 26), 4: np.linspace(1700, 2200, 26)}
    C_ = {k: 0.5 * (v[1:] + v[:-1]) for k, v in E.items()}
    ncat = {10: 5, 4: 3}
    sig, dl, Y = {}, {}, {}
    for dt in (10, 4):
        for c in range(ncat[dt]):
            s = np.exp(-0.5 * ((C_[dt] - M0[dt]) / (12 + 3 * c)) ** 2) + (0.15 * np.exp(-0.5 * ((C_[dt] - M0[dt]) / 60) ** 2) if dt == 10 else 0)
            sig[(dt, c)] = s / s.sum()
            u = (C_[dt] - C_[dt][0]) / (C_[dt][-1] - C_[dt][0])
            b = np.exp(-1.2 * u) * nleg * (c + 1) / 6
            dl[(dt, c)] = rng.poisson(b).astype(float)
            Y[(dt, c)] = 0.2
    ncat_pairs = {}
    I0 = Inputs(edges=E, ncat=ncat, sig=sig, sig_comp={}, Y=Y, dataleg=dl, nev_legs=1e6, ndata=1e6, n_real={10: 3000.0, 4: 2000.0},
                n_pair={"D0D0": npair, "DpD0": npair, "DpDp": npair}, sbsb={}, meta={})
    # sbsb z niezależnych nóg × czynnik korelacji 1
    from charmtag.pairs.model import Shapes, cat_pairs
    S = Shapes(I0)
    for ch, (a, b) in CHANNELS.items():
        for i, j in cat_pairs(ncat, a, b):
            BB2 = np.outer(S.fshape(a, S.TH[(a, i)]), S.fshape(b, S.TH[(b, j)]))
            fr = BB2[np.ix_(I0.SBM[a], I0.SBM[b])].sum()
            ncat_pairs[(ch, i, j)] = int(rng.poisson(1e6 * S.NF[(a, i)] * S.NF[(b, j)] * (0.5 if (a == b and i == j) else 1) * fr * 3))
    I0.sbsb = ncat_pairs
    return I0


def test_inputs_roundtrip(tmp_path):
    I = toy_inputs()
    I.to_files(str(tmp_path / "pi"))
    J = Inputs.from_files(str(tmp_path / "pi"))
    assert J.ncat == I.ncat and np.allclose(J.sig[(10, 3)], I.sig[(10, 3)]) and J.sbsb == I.sbsb


@pytest.mark.slow
def test_asimov_unbiased():
    M = Model(toy_inputs())
    a = F.asimov(M)
    for ch in CHANNELS:
        assert abs(a[ch]["mu"] - 1.0) < 0.02, a[ch]
        assert a[ch]["Z"] > 0.5


@pytest.mark.slow
def test_toys_pull():
    M = Model(toy_inputs(npair=800.0))
    cf = M.channel("DpDp")
    t = F.toys(M, "DpDp", cf.lam(cf.x0()), 12, seed=5)
    assert abs(t["pull_mean"]) < 0.8 and 0.4 < t["pull_std"] < 1.8


@pytest.mark.slow
def test_kappa_variants_and_syst_run():
    M = Model(toy_inputs())
    cf = M.channel("DpD0")
    n = cf.lam(cf.x0())
    for kw in ({"kappa": "channel"}, {"kappa": "partner", "partner_sigma": 0.5}, {"syst": True}):
        c2 = M.channel("DpD0", **kw)
        m = c2.fit(n)
        assert abs(m.values["mu"] - 1.0) < 0.05, (kw, m.values["mu"])


def test_region_masks():
    I = toy_inputs()
    for ch in CHANNELS:
        b, f, s, a = (region_mask(I, ch, r) for r in ("blind", "full", "sbsb", "arms"))
        assert f.all() and not (b & region_mask(I, ch, "centre")).any() and (s | a | region_mask(I, ch, "centre")).all()


# ----------------------------------------------------------------------------- nogi i pary
def toy_legs(seed=2, nev=3000):
    rng = np.random.default_rng(seed)
    rows = []
    cand = 0
    for e in range(nev):
        for _ in range(rng.integers(1, 4)):
            dt = int(rng.choice([10, 4]))
            m = rng.uniform(1750, 2000) if dt == 10 else rng.uniform(1700, 2200)
            t = rng.integers(0, 50, 3)
            rows.append(("428770", e, cand, dt, m, int(rng.integers(0, 5 if dt == 10 else 3)), t[0], t[1], t[2] if dt == 4 else -1))
            cand += 1
    return pd.DataFrame(rows, columns=["source_file", "entry", "cand", "decay_type", "m", "cat", "trk0", "trk1", "trk2"])


def test_blind_centre_never_filled():
    I = toy_inputs()
    L = toy_legs()
    Hb, cb = build_pair_hists(L, I, "blind")
    Hf, cf = build_pair_hists(L, I, "full")
    for ch in CHANNELS:
        cen = region_mask(I, ch, "centre")
        assert Hb[ch][:, cen].sum() == 0
        assert np.array_equal(Hb[ch][:, ~cen], Hf[ch][:, ~cen])
    assert cf["pairs_kept"] > cb["pairs_kept"]


def test_shared_track_removed():
    I = toy_inputs()
    L = pd.DataFrame([("1", 7, 0, 4, 1869.0, 0, 11, 12, 13), ("1", 7, 1, 10, 1865.0, 0, 13, 20, -1),
                      ("1", 7, 2, 10, 1800.0, 0, 30, 31, -1)],
                     columns=["source_file", "entry", "cand", "decay_type", "m", "cat", "trk0", "trk1", "trk2"])
    H, c = build_pair_hists(L, I, "full")
    assert c["pairs_kept"] == 2          # (0,1) dzieli ślad 13


def test_arbitration_mass_blind():
    df = pd.DataFrame({"source_file": ["1"] * 4, "entry": [5, 5, 6, 6], "cand": [0, 1, 2, 3], "trk0": [1, 2, 3, 4], "trk1": [2, 1, 4, 3],
                       "m": [1800.0, 1900.0, 1850.0, 1870.0]})
    a = B.arbitrate(df)
    b = B.arbitrate(df.assign(m=df.m[::-1].to_numpy()))
    assert len(a) == 2 and list(a.cand) == list(b.cand)


def test_categories_monotonic():
    sig = pd.DataFrame({"score": np.random.default_rng(1).random(100000)})
    thr = B.thresholds(sig, [0.3, 0.45, 0.65, 0.9])
    assert thr["0.3"] > thr["0.45"] > thr["0.65"] > thr["0.9"]
    c = B.cat_of(sig.score.to_numpy(), thr, B.DEFAULTS["cats"])
    frac = [np.mean(c == k) for k in range(3)]
    assert np.allclose(frac, [0.30, 0.15, 0.20], atol=0.01)


def test_track_id():
    pt = np.array([1000.5, 1000.5, -999.0], np.float32)
    eta = np.array([0.5, 0.5, -999.0], np.float32)
    t = B._trk_id(pt, eta)
    assert t[0] == t[1] and t[2] == -1


def test_weight_factors():
    cs, cp, rw = B.weight_factors(MB)
    assert 0.8 < cs[10] < 1.1 and cp["D0D0"] > 1.2


# ----------------------------------------------------------------------------- end-to-end na syntetycznych h5
@pytest.fixture(scope="module")
def synth_run(tmp_path_factory):
    from charmtag.synth import make
    d = tmp_path_factory.mktemp("pairs")
    make(str(d / "h5"), n=6000, seed=4)
    json.dump(MB, open(d / "mb.json", "w"))
    cfg = C.load(None, {"paths": {"h5_dir": str(d / "h5"), "out_dir": str(d / "out")}, "training": {"xgb": {"n_estimators": 40}},
                        "pairs": {"mb_rates_json": str(d / "mb.json")}})
    for sp in ("D0", "Dplus"):
        P.cmd_train(cfg, sp, "H")
        P.cmd_apply(cfg, sp, "H")
    return cfg


def test_pairs_end_to_end_blind(synth_run):
    I, H, counts = B.cmd_pairs(synth_run, unblind=False)
    out = B.pdir(synth_run)
    assert os.path.exists(os.path.join(out, "pair_inputs.json")) and os.path.exists(os.path.join(out, "pairs_blind.npz"))
    assert not os.path.exists(os.path.join(out, "pairs_full.npz"))
    for ch in CHANNELS:
        assert H[ch][:, region_mask(I, ch, "centre")].sum() == 0
    assert all(np.isfinite(v) for v in I.n_pair.values())
    assert set(I.Y) == {(10, c) for c in range(5)} | {(4, c) for c in range(3)}


def test_fit2d_cli_refuses_unblind_without_hists(synth_run):
    from charmtag.pairs.cli import main
    B.cmd_pairs(synth_run, unblind=False)
    with pytest.raises(SystemExit):
        main(["fit2d", "--unblind", "--out-dir", synth_run["paths"]["out_dir"]])


# ----------------------------------------------------------------------------- domknięcie z 09-27 (tylko gdy są zamrożone wejścia)
@pytest.mark.skipif(not os.path.isdir(FREEZE), reason="brak unblinding_freeze_v20260927_01")
def test_closure_0927_fast():
    from charmtag.pairs.frozen import load_frozen, frozen_legs_all
    I = load_frozen(FREEZE)
    H, _ = build_pair_hists(frozen_legs_all(FREEZE), I, "full")
    if os.path.exists(STAGEB):
        old = np.load(STAGEB)
        assert all(np.array_equal(old[ch], H[ch]) for ch in CHANNELS)
    M = Model(I)
    cf = M.channel("DpDp")
    m1 = cf.fit(H["DpDp"])
    assert abs(m1.values["mu"] - 1.5748) < 2e-3
    cf = M.channel("D0D0")
    m1 = cf.fit(H["D0D0"])
    m0 = cf.fit(H["D0D0"], fix_mu=0.0, start=np.array(m1.values))
    assert abs(m1.values["mu"] - 1.3664) < 2e-3 and abs(F.z_from(m1, m0) - 1.884) < 0.01


def test_check_event_key_unique():
    """Audyt 16.11 (errata v2.0): powtórzony event_number w różnych blokach wierszy → SystemExit; bloki rozłączne → OK."""
    import pandas as pd
    import pytest
    from charmtag.pairs.build import check_event_key_unique
    ok = pd.DataFrame({"source_file": ["0"] * 5, "entry": [1, 1, 2, 3, 3], "cand": [0, 1, 2, 3, 5]})
    check_event_key_unique(ok)
    bad = pd.DataFrame({"source_file": ["0"] * 5, "entry": [1, 1, 2, 1, 3], "cand": [0, 1, 2, 3, 5]})
    with pytest.raises(SystemExit):
        check_event_key_unique(bad)


def test_mc_pair_rel_stat_neff():
    """VER_budget 16.11: mc_stat z n_eff (config > mb_rates_json > fallback 1/√n)."""
    import numpy as np
    from charmtag.pairs.build import mc_pair_rel_stat
    n = {"D0D0": 21, "DpD0": 41, "DpDp": 15}
    r = mc_pair_rel_stat(n, {}, {"mc_pairs_neff": {"D0D0": 8.2, "DpD0": 13.5, "DpDp": 7.5}})
    assert np.isclose(r["D0D0"], 1 / np.sqrt(8.2)) and np.isclose(r["DpDp"], 1 / np.sqrt(7.5))
    r = mc_pair_rel_stat(n, {"pairs_D0D0": {"n_eff": 4.0}}, {})
    assert np.isclose(r["D0D0"], 0.5) and np.isclose(r["DpD0"], 1 / np.sqrt(41))
