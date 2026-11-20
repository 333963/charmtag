"""Testy etapu S8: SALAD (reweighting MC tła → dane w SB) — SaladReweighter8 (sklearn), bkg_source synthetic/salad,
`charmtag s8salad` (systematyka kształtu po cięciu), dispatch CLI. Syntetyczny h5 (charmtag.synth) — mechanika, nie fizyka."""
import json

import numpy as np
import pytest

from charmtag import bkgsrc as BS
from charmtag import config as C
from charmtag import pipeline as P
from charmtag import samples as SA
from charmtag.bkggen.salad8 import SaladReweighter8
from charmtag.synth import make

MB = {"dstar_true_per_event": {"unw": 6.67e-05, "w": 4.51e-04}, "D0_correct_per_event": {"unw": 4.79e-4, "w": 3.03e-3},
      "Dplus_correct_per_event": {"unw": 4.40e-4, "w": 2.36e-3},
      "pairs_D0D0": {"per_event_unw": 1.07e-6, "per_event_w": 9.9e-6}, "pairs_DpD0": {"per_event_unw": 2.08e-6, "per_event_w": 1.69e-5},
      "pairs_DpDp": {"per_event_unw": 7.6e-7, "per_event_w": 5.05e-6}}


@pytest.fixture(scope="module")
def h5dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("s8")
    make(str(d / "h5"), n=5000, seed=8)
    json.dump(MB, open(d / "mb.json", "w"))
    return d


def _cfg(d, name, **ov):
    base = {"paths": {"h5_dir": str(d / "h5"), "out_dir": str(d / name)}, "training": {"xgb": {"n_estimators": 40}},
            "pairs": {"mb_rates_json": str(d / "mb.json")}}
    return C.load(None, C._merge(base, ov))


def _toy(n, rng, shift=0.0, m0=1.87):
    m = rng.uniform(1.70, 2.20, n)
    t = (m - m0) / 0.1
    x1 = rng.normal(shift * t, 1.0, n)           # zależne od m tylko w „danych”
    x2 = rng.exponential(1.0, n) + shift
    return np.column_stack([x1, x2]), m


@pytest.mark.parametrize("kind", ["sk_mlp", "hgb", "hgb_nom"])
def test_reweighter_closes_sb_shift_and_interpolates(kind):
    rng = np.random.default_rng(1)
    sp = {"m0": 1.87}
    Xm, mm = _toy(40000, rng, 0.0)
    Xd, md = _toy(40000, rng, 0.6)
    sb = lambda m: np.abs(m - sp["m0"]) > 0.07  # noqa: E731
    rw = SaladReweighter8(["a", "b"], sp, kind=kind, seed=0, max_iter=30).fit(Xm[sb(mm)], mm[sb(mm)], Xd[sb(md)], md[sb(md)], max_n=30000)
    assert 0.5 < rw.val_auc < 1.0
    w = rw.weights(Xm, mm)
    assert w.shape == (len(Xm),) and abs(w.mean() - 1) < 1e-6 and w.min() > 0
    # ważona średnia x2 w SR ≈ średnia danych w SR (interpolacja); surowe MC jest przesunięte o −0,6
    sr = ~sb(mm); srd = ~sb(md)
    mean_w = np.average(Xm[sr, 1], weights=w[sr]); mean_raw = Xm[sr, 1].mean(); mean_d = Xd[srd, 1].mean()
    assert abs(mean_w - mean_d) < 0.5 * abs(mean_raw - mean_d), (mean_w, mean_raw, mean_d)


def test_bad_kind_rejected():
    with pytest.raises(ValueError):
        SaladReweighter8(["a"], {"m0": 1.87}, kind="nope")


def test_api_salad_uses_sklearn_default():
    from charmtag.bkggen import api
    g = api.make_generator("salad", species="D0")
    assert g.name == "salad" and "salad" in api.available()


def test_training_background_salad_source(h5dir):
    cfg = _cfg(h5dir, "sal", training={"bkg_source": "synthetic", "bkg_options": {"synthetic_model": "salad", "salad_kind": "hgb",
                                                                                   "salad_max_n": 20000}})
    assert BS.source(cfg) == "synthetic"
    B = BS.training_background(cfg, "D0")
    sp = cfg["species"]["D0"]
    assert len(B) > 0 and (np.abs(B["m"] - sp["m0"]) < sp["sideband_min"]).all()
    assert "w_salad" in B and B.attrs["synthetic"]["model"] == "salad" and 0 < B.attrs["synthetic"]["neff_over_n"] <= 1.0
    assert (B["weight"] > 0).all()
    S, B2 = SA.training_set(cfg, "D0")
    assert (B2["y"] == 0).all() and np.isfinite(B2["w_train"]).all() and B2["w_train"].sum() > 0


def test_default_synthetic_is_still_copula(h5dir):
    cfg = _cfg(h5dir, "cop", training={"bkg_source": "synthetic"})
    assert BS.options(cfg)["synthetic_model"] == "copula"


def test_s8salad_shape_syst_and_cli(h5dir, tmp_path):
    from charmtag import s8salad as S8
    cfg = _cfg(h5dir, "s8sys", training={"bkg_options": {"salad_kind": "hgb", "salad_max_n": 20000}})
    P.cmd_train(cfg, "D0", "H")
    S8.shape_syst(cfg, "D0", "H", out=str(tmp_path))
    p = tmp_path / "bkg_shape_salad_D0.json"
    assert p.exists()
    j = json.load(open(p))
    for k in ("effB_salad_SR", "effB_mcraw_SR", "effB_data_SB", "rel_salad_vs_SB", "definition", "c2st"):
        assert k in j
    assert set(j["effB_salad_SR"]) >= {"3-5", "5-8", "8-30", "all"}
    assert all(0 <= v <= 1 for v in j["effB_salad_SR"]["all"].values())
    # dispatch CLI
    from charmtag import cli
    assert "s8salad" in cli._SUBCOMMANDS   # REPO: dispatch tabelą, nie łańcuchem if


def test_s8salad_max_rows_subsamples(h5dir, tmp_path):
    """S5MEM: limit salad_max_rows → keep_frac < 1, mniej wierszy w pamięci, ε_B zgodne w granicach statystyki."""
    from charmtag import s8salad as S8
    full = _cfg(h5dir, "s8full", training={"bkg_options": {"salad_kind": "hgb", "salad_max_n": 20000, "salad_max_rows": 0}})
    P.cmd_train(full, "D0", "H")
    n_d, _ = SA.count_species(full, "data_none", "D0")
    assert S8.keep_frac_for(full, "data_none", "D0") == 1.0                    # 0 → bez limitu
    lim = max(200, n_d // 3)
    sub = _cfg(h5dir, "s8full", training={"bkg_options": {"salad_kind": "hgb", "salad_max_n": 20000, "salad_max_rows": lim}})
    f = S8.keep_frac_for(sub, "data_none", "D0")
    assert 0 < f < 1 and abs(f - lim / n_d) < 1e-9
    a = S8.shape_syst(full, "D0", "H", out=str(tmp_path / "a"))
    b = S8.shape_syst(sub, "D0", "H", out=str(tmp_path / "b"))
    assert b["n_data"] < a["n_data"] and b["keep_frac"]["data_none"] < 1
    assert abs(b["n_data"] - lim) < 5 * np.sqrt(lim) + 0.05 * lim
    for wp, ea in a["effB_data_SB"]["all"].items():                            # stosunki niezmienione w granicach stat.
        eb = b["effB_data_SB"]["all"][wp]
        n_sb = max(b["n_data"] * 0.3, 1)
        assert abs(eb - ea) < 5 * np.sqrt(max(ea * (1 - ea), 1e-3) / n_sb) + 0.02
