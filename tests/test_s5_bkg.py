"""Testy etapu S5: training.bkg_source, szablon refleksji z danych w fit masy, extra_wps/output_tag, scores_dir/models_dir,
mixing nóg w pairs, s5compare. Na syntetycznym h5 (charmtag.synth) — tylko mechanika, nie fizyka."""
import json
import os

import numpy as np
import pandas as pd
import pytest

from charmtag import bkgsrc as BS
from charmtag import config as C
from charmtag import pipeline as P
from charmtag import refldata as RD
from charmtag import samples as SA
from charmtag.synth import make

MB = {"dstar_true_per_event": {"unw": 6.67e-05, "w": 4.51e-04}, "D0_correct_per_event": {"unw": 4.79e-4, "w": 3.03e-3},
      "Dplus_correct_per_event": {"unw": 4.40e-4, "w": 2.36e-3},
      "pairs_D0D0": {"per_event_unw": 1.07e-6, "per_event_w": 9.9e-6}, "pairs_DpD0": {"per_event_unw": 2.08e-6, "per_event_w": 1.69e-5},
      "pairs_DpDp": {"per_event_unw": 7.6e-7, "per_event_w": 5.05e-6}}


@pytest.fixture(scope="module")
def h5dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("s5")
    make(str(d / "h5"), n=5000, seed=5)
    json.dump(MB, open(d / "mb.json", "w"))
    return d


def _cfg(d, name, **ov):
    base = {"paths": {"h5_dir": str(d / "h5"), "out_dir": str(d / name)}, "training": {"xgb": {"n_estimators": 40}},
            "pairs": {"mb_rates_json": str(d / "mb.json")}}
    return C.load(None, C._merge(base, ov))


def test_default_is_sideband_and_unknown_rejected(h5dir):
    cfg = _cfg(h5dir, "x")
    assert BS.source(cfg) == "sideband"
    assert "bkg_source" not in cfg["training"]  # default.yaml nietknięty
    with pytest.raises(ValueError):
        BS.source({"training": {"bkg_source": "nope"}})


def test_swapped_mass_is_involution_like():
    rng = np.random.default_rng(0)
    n = 1000
    L = pd.DataFrame({"d_pt0": rng.uniform(1, 5, n), "d_pt1": rng.uniform(1, 5, n), "d_eta0": rng.normal(0, 1, n),
                      "d_eta1": rng.normal(0, 1, n), "d_phi0": rng.uniform(-3, 3, n), "d_phi1": rng.uniform(-3, 3, n)})
    sp = {"kaon_index": 0, "pion_indices": [1]}
    a = BS.swapped_mass(L, sp)
    Lb = L.rename(columns={"d_pt0": "d_pt1", "d_pt1": "d_pt0", "d_eta0": "d_eta1", "d_eta1": "d_eta0", "d_phi0": "d_phi1", "d_phi1": "d_phi0"})
    b = BS.swapped_mass(Lb, {"kaon_index": 1, "pion_indices": [0]})  # ta sama fizyka, inne indeksy
    assert np.allclose(a, b)
    s = BS._fv(L.d_pt0.values, L.d_eta0.values, L.d_phi0.values, BS.M_K) + BS._fv(L.d_pt1.values, L.d_eta1.values, L.d_phi1.values, BS.M_PI)
    nom = np.sqrt(s[:, 0] ** 2 - (s[:, 1:] ** 2).sum(1))
    assert not np.allclose(a, nom)


@pytest.mark.parametrize("species", ["D0", "Dplus"])
def test_far_veto_source(h5dir, species):
    cfg = _cfg(h5dir, f"fv_{species}", training={"bkg_source": "sideband_far_veto"})
    sp = cfg["species"][species]
    S, B = SA.training_set(cfg, species)
    far = BS.options(cfg)["far_min"][species]
    assert len(B) > 50 and (np.abs(B["m"] - sp["m0"]) > far).all()
    if species == "D0":
        assert (np.abs(BS.swapped_mass(B, sp) - sp["m0"]) >= 0.025).all()
    assert np.isclose(B["w_train"].sum(), S["weight"].sum())


def test_mc_source(h5dir):
    cfg = _cfg(h5dir, "mc", training={"bkg_source": "mc"})
    S, B = SA.training_set(cfg, "D0")
    t = B["truth"].to_numpy()
    assert len(B) > 50 and ((t >= 8) & (t < 100)).all()
    assert (np.abs(B["m"] - 1.86484) < 0.075).any()  # region: all → także pod pikiem


def test_synthetic_source_and_train(h5dir):
    cfg = _cfg(h5dir, "syn", training={"bkg_source": "synthetic"})
    S, B = SA.training_set(cfg, "D0")
    sp = cfg["species"]["D0"]
    assert len(B) > 50 and (np.abs(B["m"] - sp["m0"]) <= sp["sideband_min"] + 1e-9).all()
    from charmtag import features as F
    for f in F.feature_list(cfg, "H"):
        if f in B:
            x = B[f].to_numpy(float)
            assert np.isfinite(x).mean() > 0.3, f
    meta, metrics = P.cmd_train(cfg, "D0", "H")
    assert all(np.isfinite(v["auc"]) for v in metrics.values())


def test_refl_template_recovers_shape():
    rng = np.random.default_rng(1)
    m0, lo, hi = 1.86484, 1.75, 2.0
    n_ref, n_comb = 20000, 200000
    m_ref = rng.normal(m0, 0.06, n_ref)                   # refleksje: m szerokie, m_swap w piku
    ms_ref = rng.normal(m0, 0.010, n_ref)
    m_c = rng.uniform(lo, hi, n_comb)                      # kombinatoryka: m_swap płaskie
    ms_c = rng.uniform(1.6, 2.2, n_comb)
    o = dict(RD.DEFAULTS, smooth_bins=0)
    h, npk, tot = RD.template(np.r_[m_ref, m_c], np.r_[ms_ref, ms_c], m0, lo, hi, 50, o)
    exp_ref = ((np.abs(ms_ref - m0) < 0.025) & (m_ref >= lo) & (m_ref < hi)).sum()
    assert abs(tot - exp_ref) < 5 * np.sqrt(npk)
    shape = {"model": "gauss", "mu": m0, "s1": 0.0125}
    assert abs(RD.f_window(shape, m0, 0.025) - 0.9545) < 1e-3


def test_fit_refl_data_extra_wps_and_overrides(h5dir):
    cfg = _cfg(h5dir, "nom")
    P.cmd_train(cfg, "D0", "H")
    P.cmd_apply(cfg, "D0", "H")
    r0 = P.cmd_fit(cfg, "D0", "H")
    nom = cfg["paths"]["out_dir"]
    for norm in ("mc_ratio", "data"):
        c2 = _cfg(h5dir, f"refl_{norm}", paths={"models_dir": os.path.join(nom, "models"), "scores_dir": os.path.join(nom, "scores")},
                  fit={"reflection_source": "data", "refl_data": {"norm": norm, "min_peak": 20}, "extra_wps": [0.8, 0.5],
                       "output_tag": "_refldata"})
        r = P.cmd_fit(c2, "D0", "H")
        p = os.path.join(c2["paths"]["out_dir"], "fits", "D0_H_refldata.json")
        assert os.path.exists(p)
        assert not os.path.exists(os.path.join(c2["paths"]["out_dir"], "scores"))  # wyniki z nominalnego
        e = next(iter(r["bins"].values()))
        assert {"data_none_wp80", "data_none_wp70", "data_none_wp50"} <= set(e)
        assert "refl_data_wp70" in e
        f = e["data_none_wp70"]
        if "N_sig" in f:
            assert f.get("refl_source") == "data" and f["refl_norm"] == norm
    # nominalny fit bez zmian (brak kluczy S5)
    e0 = next(iter(r0["bins"].values()))
    assert "refl_data_wp70" not in e0 and "data_none_wp80" not in e0


def test_mixing_independent_legs_closure():
    from charmtag.pairs import CHANNELS
    from charmtag.pairs.mixctl import mixed_pairs, hist_pairs, _masks
    from charmtag.pairs.hists import leg_bins, build_pair_hists
    from charmtag.pairs.model import Inputs
    rng = np.random.default_rng(3)
    E = {10: np.linspace(1750, 2000, 26), 4: np.linspace(1700, 2200, 26)}
    ncat = {10: 5, 4: 3}
    I = Inputs(edges=E, ncat=ncat, sig={}, sig_comp={}, Y={}, dataleg={}, nev_legs=1, ndata=1, n_real={}, n_pair={}, sbsb={})
    nev = 40000
    k = rng.integers(2, 5, nev)
    ev = np.repeat(np.arange(nev), k)
    dt = np.where(rng.random(len(ev)) < 0.5, 10, 4)
    m = np.where(dt == 10, 1750 + rng.exponential(120, len(ev)), 1700 + rng.exponential(250, len(ev)))
    L = pd.DataFrame({"source_file": "428770", "entry": ev, "cand": np.arange(len(ev)), "decay_type": dt, "m": m,
                      "cat": np.where(dt == 10, rng.integers(0, 5, len(ev)), rng.integers(0, 3, len(ev))),
                      "trk0": np.arange(len(ev)) * 3, "trk1": np.arange(len(ev)) * 3 + 1, "trk2": -1})
    Lb = leg_bins(L, I)
    P_ = mixed_pairs(Lb, I, factor=2, seed=1)
    assert len(P_) > 0 and not ((P_["pool1"] == P_["pool2"]) & (P_["pos1"] == P_["pos2"])).any()
    mix = hist_pairs(P_, I)
    same, _ = build_pair_hists(L, I, "sbsb")
    for ch in CHANNELS:
        sb, arms, centre = _masks(I, ch)
        assert (same[ch][:, ~sb] == 0).all()
        k_ = same[ch][:, sb].sum() / mix[ch][:, sb].sum()
        pred = mix[ch] * k_
        g = sb[None] & (pred > 5)
        chi2 = ((same[ch][g] - pred[g]) ** 2 / pred[g]).sum() / g.sum()
        assert chi2 < 1.6, (ch, chi2)


def test_pairs_mix_control_end_to_end(h5dir):
    from charmtag.pairs import build as B
    cfg = _cfg(h5dir, "pmix", pairs={"mix_control": True, "mix_factor": 1})
    for sp in ("D0", "Dplus"):
        P.cmd_train(cfg, sp, "H")
        P.cmd_apply(cfg, sp, "H")
    I, H, counts = B.cmd_pairs(cfg, unblind=False)
    out = B.pdir(cfg)
    js = json.load(open(os.path.join(out, "mix_control.json")))
    assert set(js["channels"]) == {"D0D0", "DpD0", "DpDp"}
    assert os.path.exists(os.path.join(out, "mix_hists.npz"))
    assert not os.path.exists(os.path.join(out, "pairs_full.npz"))
    # s5compare na tym biegu (nominal = ten sam katalog → zysk 1)
    from charmtag.s5cmp import compare
    res, md = compare(cfg, {"nom": cfg["paths"]["out_dir"], "kopia": cfg["paths"]["out_dir"]}, "nom", species=("D0",), with_eff=True)
    assert "| nom |" in md and "Pary D–D" in md
    r = res["variants"]["kopia"]["species"]["D0"]["bins"]
    assert all(v.get("effB_dataSB_70") is not None for v in r.values())


def test_cli_dispatch_s5compare(h5dir, tmp_path):
    from charmtag.cli import main
    cfg = _cfg(h5dir, "nom")
    out = tmp_path / "cmp"
    rc = main(["s5compare", "--run", f"a={cfg['paths']['out_dir']}", "--out", str(out), "--species", "D0", "--no-eff",
               "-c", os.path.join(os.path.dirname(C.DEFAULT), "default.yaml")])
    assert rc == 0 and (tmp_path / "cmp.md").exists()


def test_refl_hybrid_centre_from_mc():
    e = np.linspace(1.75, 2.0, 26)
    m0 = 1.86484
    hm = np.exp(-0.5 * ((0.5 * (e[1:] + e[:-1]) - m0) / 0.06) ** 2) * 1000
    hd = hm * 0.5
    c = np.abs(0.5 * (e[1:] + e[:-1]) - m0) < 0.025
    hd[c] = 0  # dane: centrum zniszczone (mieszanie z sygnałem)
    h, tails, fc = RD.hybrid(hd, hm, e, m0, 0.025)
    assert np.allclose(h, hm * 0.5) and np.isclose(tails, 0.5 * hm[~c].sum()) and 0 < fc < 1
