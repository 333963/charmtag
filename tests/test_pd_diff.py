"""PD: pary różniczkowe (LS/OS, Δφ) — zmienne par, biny, bitowa identyczność starych wyjść, odzyskanie μ_i, end-to-end na syntetyku."""
import hashlib
import json
import os
import shutil

import h5py
import numpy as np
import pandas as pd
import pytest

from charmtag import config as C
from charmtag import pipeline as P
from charmtag.pairs import CHANNELS
from charmtag.pairs import build as B
from charmtag.pairs import diff as DF
from charmtag.pairs.frozen import sbsb_counts
from charmtag.pairs.hists import build_pair_hists
from charmtag.pairs.model import Model, region_mask
from test_pairs import MB, toy_inputs, toy_legs

# hash histogramów (toy_legs seed 2 × toy_inputs seed 1) z kodu PRZED etapem PD (snapshot AA_merge, 2.0.3) — bitowa identyczność
BASE_HASH = {"blind": "b6d566d96fe5be89d1643549484618f6d17d9c1f6a2e69f41b7c1dc2006979dd",
             "full": "f5234dcda00567b5d1a87ccc71f1a729671e049dfd912e39ee8146f8b53369f2"}


def _hash(H):
    return hashlib.sha256(b"".join(H[ch].astype("<f8").tobytes() for ch in CHANNELS)).hexdigest()


def _kin_legs(seed=2, v11=False):
    """Nogi toy + kinematyka; v11=True → bez ładunku córki 0 (jak h5 v1.1: D0 nieokreślone w `sign`)."""
    L = toy_legs(seed)
    rng = np.random.default_rng(seed + 100)
    n = len(L)
    return L.assign(eta=rng.uniform(-2.5, 2.5, n), phi=rng.uniform(-np.pi, np.pi, n), charge=np.where(L.decay_type == 4, rng.choice([-1, 1], n), 0),
                    dq0=np.where((L.decay_type == 10) & (not v11), rng.choice([-1, 1], n), np.nan))


# ----------------------------------------------------------------------------- zmienne par
def test_wrap_dphi_and_sign():
    assert np.allclose(DF.wrap_dphi([0.1, 3.0, -3.0, 0.0], [0.1, -3.0, 3.0, np.pi]), [0.0, 2 * np.pi - 6.0, 2 * np.pi - 6.0, np.pi])
    s = DF.charm_sign([4, 4, 10, 10, 10, 10], [1, -1, 0, 0, 0, 0], [np.nan, np.nan, -1, 1, np.nan, -999])
    assert list(s) == [1, -1, 1, -1, 0, 0]      # D⁺ → c, D⁻ → c̄, D0 (K⁻) → c, anty-D0 (K⁺) → c̄, bez ładunku K → nieokreślone


def test_diff_bin_variants():
    P = pd.DataFrame({"decay_type1": [4, 4, 4, 10, 10], "decay_type2": [4, 4, 10, 10, 10],
                      "phi1": [0.0, 0.0, 0.0, 1.0, 0.0], "phi2": [0.5, 3.0, np.pi, -1.0, 0.0], "eta1": [0.0, 0.0, 1.0, 0.0, 0.0], "eta2": [0.5, 2.0, -1.0, 0.0, 0.0],
                      "charge1": [1, 1, 1, 0, 0], "charge2": [1, -1, 0, 0, 0], "dq01": [np.nan] * 3 + [-1, -1], "dq02": [np.nan, np.nan, -1, -1, np.nan]})
    sd = DF.spec_from_cfg({"diff": {"var": "dphi"}})
    assert list(DF.diff_bin(P, sd)) == [0, 1, 1, 1, 0]                   # π dokładnie → ostatni bin domknięty
    sy = DF.spec_from_cfg({"diff": {"var": "dy", "edges": [0, 1, 5]}})
    assert list(DF.diff_bin(P, sy)) == [0, 1, 1, 0, 0]
    ss = DF.spec_from_cfg({"diff": {"var": "sign"}})
    assert list(DF.diff_bin(P, ss)) == [0, 1, 0, 0, -1]                  # D⁺D⁺ like, D⁺D⁻ opp, D⁺D0(K⁻) = c c like, D0D0(K⁻K⁻) like, bez K → nieokr.
    sx = DF.spec_from_cfg({"diff": {"var": "sign_dphi"}})
    assert DF.labels(sx) == ["like", "opp:[0,1.57)", "opp:[1.57,3.14]"]
    assert list(DF.diff_bin(P, sx)) == [0, 2, 0, 0, -1]
    assert np.allclose(DF.flat_fractions(sx), [0.5, 0.25, 0.25]) and np.allclose(DF.flat_fractions(sd), [0.5, 0.5])


def test_spec_defaults_and_errors():
    assert DF.spec_from_cfg({}) is None and B.pcfg({})["diff"] is None
    s = DF.spec_from_cfg({"diff": True})
    assert s["var"] == "sign" and s["channels"] == list(CHANNELS) and s["expect"] == "mc"
    s = DF.spec_from_cfg({}, "dphi")
    assert s["edges"] == [0.0, np.pi / 2, np.pi] and s["channels"] == ["D0D0", "DpD0"]
    s = DF.spec_from_cfg({"diff": {"var": "dphi", "edges": [0, 1, 2, 3.1416]}}, "cfg")
    assert len(s["edges"]) == 4
    for bad in ({"var": "x"}, {"var": "dphi", "edges": [1, 0]}, {"expect": "y"}, {"channels": ["Q"]}):
        with pytest.raises(SystemExit):
            DF.spec_from_cfg({"diff": bad})


# ----------------------------------------------------------------------------- histogramy: stare wyjścia bez zmian, biny sumują się do całości
def test_hists_bitwise_identical_without_diff():
    I = toy_inputs()
    for reg in ("blind", "full"):
        H, c = build_pair_hists(toy_legs(), I, reg)
        assert _hash(H) == BASE_HASH[reg], reg
        H2, c2 = build_pair_hists(_kin_legs(), I, reg)          # dodatkowe kolumny nóg nie zmieniają histogramów
        assert _hash(H2) == BASE_HASH[reg] and c2["pairs_kept"] == c["pairs_kept"]


def test_hists_diff_bins_partition_total():
    I = toy_inputs()
    L = _kin_legs()
    H, c = build_pair_hists(L, I, "full")
    for var in ("dphi", "sign", "sign_dphi"):
        spec = DF.spec_from_cfg({"diff": {"var": var}})
        Hd_tot, cd, Hd = build_pair_hists(L, I, "full", diff=spec)
        assert _hash(Hd_tot) == BASE_HASH["full"] and cd["pairs_kept"] == c["pairs_kept"]
        nb = DF.nbins(spec)
        assert set(Hd) == {DF.hkey(ch, k) for ch in CHANNELS for k in range(nb)}
        tot = sum(cd["diff_kept"].values()) + cd["diff_undefined"]
        assert tot == c["pairs_kept"]
        for ch in CHANNELS:
            S = sum(Hd[DF.hkey(ch, k)] for k in range(nb))
            assert (S <= H[ch]).all()
            if var == "dphi":
                assert np.array_equal(S, H[ch])                 # Δφ zawsze określone
        if var == "sign":
            assert cd["diff_undefined"] == 0                    # „v1.2": wszystkie nogi z ładunkiem → wszystko określone
            _, cv, Hv = build_pair_hists(_kin_legs(v11=True), I, "full", diff=spec)
            assert cv["diff_undefined"] > 0                     # v1.1: D0 nogi bez dq0 → nieokreślone, D⁺D⁺ kompletne
            assert np.array_equal(sum(Hv[DF.hkey("DpDp", k)] for k in range(2)), H["DpDp"])
            assert all(Hv[DF.hkey(ch, k)].sum() == 0 for ch in ("D0D0", "DpD0") for k in range(2))
    # ślepo: centrum nigdy w binach
    spec = DF.spec_from_cfg({"diff": {"var": "dphi"}})
    _, _, Hb = build_pair_hists(L, I, "blind", diff=spec)
    for ch in CHANNELS:
        for k in range(2):
            assert Hb[DF.hkey(ch, k)][:, region_mask(I, ch, "centre")].sum() == 0


def test_sbsb_and_mc_rates_per_bin():
    L = _kin_legs()
    ncat = {10: 5, 4: 3}
    spec = DF.spec_from_cfg({"diff": {"var": "dphi", "edges": [0, 1, 2, 3.1416]}})
    tot = sbsb_counts(L, ncat)
    per = sbsb_counts(L, ncat, diff=spec)
    assert set(per) == {0, 1, 2}
    for key, v in tot.items():
        assert sum(per[k][key] for k in per) == v
    legs_c = {dt: L[L.decay_type == dt].assign(truth=0) for dt in (10, 4)}
    D = pd.DataFrame({"m": [], "dm": [], "truth": []})
    p = B.pcfg({})
    r = B.mc_rates(legs_c, D, p, diff=spec)
    for ch in CHANNELS:
        assert sum(r["n_pairs_diff"][ch]) + r["n_pairs_diff_undefined"][ch] == r["n_pairs"][ch]
    f, e = DF.fractions_mc(r["n_pairs_diff"], spec, "DpDp")
    assert abs(f.sum() - 1) < 1e-12 and (e > 0).all()
    ff, mode = DF.fit_fractions(spec, {"DpDp": [1, 0, 0]}, "DpDp")
    assert mode.startswith("flat") and np.allclose(ff, DF.flat_fractions(spec))   # pusty bin MC → płasko


# ----------------------------------------------------------------------------- fit per bin: odzyskanie μ_i i testy kształtu
def test_diff_fit_recovers_mu_per_bin(tmp_path):
    I = toy_inputs(npair=300.0)
    spec = DF.spec_from_cfg({"diff": {"var": "sign", "channels": ["DpDp"]}})
    sb = sbsb_counts(_kin_legs(), {10: 5, 4: 3}, diff=spec)
    n_mc = {ch: [30, 70] for ch in CHANNELS}
    D = DF.diff_inputs_json(I, spec, sb, n_mc, {"diff_kept": {"0": 0, "1": 0}})
    assert abs(D["bins"][0]["n_pair"]["DpDp"] - 0.3 * 300) < 1e-9 and D["bins"][0]["f_mode"]["DpDp"] == "mc"
    mu_true = {0: 2.5, 1: 0.6}
    rng = np.random.default_rng(7)
    Hd = {}
    for k in range(2):
        cf = Model(DF.inputs_for_bin(I, D, k)).channel("DpDp")
        x = cf.x0(); x[0] = mu_true[k]
        Hd[DF.hkey("DpDp", k)] = rng.poisson(cf.lam(x)).astype(float)
    rep = DF.run_diff(I, D, Hd, str(tmp_path))
    b = rep["channels"]["DpDp"]["bins"]
    for k, lab in enumerate(D["labels"]):
        r = b[lab]
        assert abs(r["mu_hat"] - mu_true[k]) < 2.5 * r["mu_err_hesse"], (lab, r["mu_hat"], r["mu_err_hesse"])
        assert r["mu_profile_68"][0] < mu_true[k] + 2 * r["mu_err_hesse"]
    t = rep["channels"]["DpDp"]["tests"]
    assert set(t) == {"dps_flat", "pythia_mc"} and 0 <= t["pythia_mc"]["p"] <= 1
    ls = rep["like_sign"]["channels"]["DpDp"]
    assert 0 < ls["f_LS"][0] < 1 and "Z_reject_pythia" in ls
    # f_LS^dane = μ_LS f_LS^MC / Σ: tu ≈ 2.5·0.3/(2.5·0.3+0.6·0.7) = 0.64
    assert abs(ls["f_LS"][0] - 0.64) < 3 * ls["f_LS"][1] + 0.05
    assert "like" in rep["bins"] and rep["tests"]["pythia_mc"]["ndf"] == 1
    assert os.path.exists(tmp_path / "fit2d_diff.json")
    path = DF.write_report(str(tmp_path), I, D, rep)
    txt = open(path).read()
    assert "Like-sign" in txt and "DpDp | like" in txt and "0.09" in txt


def test_combine_bins_shape_scaling():
    scan = {f"{m:g}": float((m - 1.0) ** 2) for m in DF.GRID}          # parabola z minimum w μ=1
    r = DF.combine_bins({"a": scan})
    assert abs(r["mu"] - 1.0) < 0.01 and abs(r["q_shape"]) < 1e-9 and abs(r["mu_68"][1] - r["mu_68"][0] - 2 * np.sqrt(0.5)) < 0.02
    r2 = DF.combine_bins({"a": scan}, {"a": 2.0})                     # ΔNLL(2μ) → minimum w μ=0.5
    assert abs(r2["mu"] - 0.5) < 0.01
    r3 = DF.combine_bins({"a": scan, "b": {f"{m:g}": float((m - 3.0) ** 2) for m in DF.GRID}})
    assert abs(r3["mu"] - 2.0) < 0.01 and abs(r3["q_shape"] - 4.0) < 0.05    # 2·(1+1)


# ----------------------------------------------------------------------------- end-to-end na syntetyku (v1.1 i „v1.2" z ładunkiem córki)
@pytest.fixture(scope="module")
def synth_diff(tmp_path_factory):
    from charmtag.synth import make
    d = tmp_path_factory.mktemp("pd")
    make(str(d / "h5"), n=6000, seed=4)
    json.dump(MB, open(d / "mb.json", "w"))
    cfg = C.load(None, {"paths": {"h5_dir": str(d / "h5"), "out_dir": str(d / "out")}, "training": {"xgb": {"n_estimators": 40}},
                        "pairs": {"mb_rates_json": str(d / "mb.json")}})
    for sp in ("D0", "Dplus"):
        P.cmd_train(cfg, sp, "H")
        P.cmd_apply(cfg, sp, "H")
    return cfg, d


def _make_v12(src_dir, dst_dir, shuffle=False, seed=3):
    """„v1.2": te same wiersze + Dmeson_daughter0_charge (losowo ±1 dla D0, znak jak ładunek dla D⁺); shuffle → inna kolejność (test kontroli)."""
    os.makedirs(dst_dir, exist_ok=True)
    rng = np.random.default_rng(seed)
    for fn in os.listdir(src_dir):
        with h5py.File(os.path.join(src_dir, fn)) as f:
            a = f["events"][:]
        n = len(a)
        dq = np.where(a["Dmeson_decayType"] == 4, -a["Dmeson_charge"], rng.choice([-1, 1], n)).astype(np.int32)
        dq[rng.random(n) < 0.02] = -999
        dt = np.dtype([(k, a.dtype.fields[k][0]) for k in a.dtype.names] + [("Dmeson_daughter0_charge", "<i4")])
        b = np.empty(n, dt)
        for k in a.dtype.names:
            b[k] = a[k]
        b["Dmeson_daughter0_charge"] = dq
        if shuffle:
            b = b[rng.permutation(n)]
        with h5py.File(os.path.join(dst_dir, fn), "w") as f:
            f.create_dataset("events", data=b)


def test_end_to_end_sign_v11_then_v12(synth_diff):
    cfg, d = synth_diff
    out = B.pdir(cfg)
    with pytest.raises(SystemExit):
        B.cmd_pairs(cfg, unblind=False, diff="sign")
    # v1.1: D⁺D⁺ określone (Dmeson_charge), D0 nogi nieokreślone
    I, H, counts = B.cmd_pairs(cfg, unblind=True, diff="sign")
    assert os.path.exists(os.path.join(out, "pairs_full_diff.npz")) and os.path.exists(os.path.join(out, "pair_inputs_diff.json"))
    D = json.load(open(os.path.join(out, "pair_inputs_diff.json")))
    Hd = dict(np.load(os.path.join(out, "pairs_full_diff.npz")))
    assert D["labels"] == ["like", "opposite"] and I.meta["diff"]["var"] == "sign"
    assert sum(Hd[DF.hkey("DpDp", k)].sum() for k in range(2)) == H["DpDp"].sum() > 0
    assert all(Hd[DF.hkey(ch, k)].sum() == 0 for ch in ("D0D0", "DpD0") for k in range(2))
    assert counts["diff_undefined"] > 0 and D["n_mc_pairs"]["DpDp"][0] + D["n_mc_pairs"]["DpDp"][1] == I.meta["mc_rates_unweighted"]["n_pairs"]["DpDp"]
    H0 = dict(np.load(os.path.join(out, "pairs_full.npz")))
    legs = pd.read_parquet(os.path.join(out, "legs_data.parquet")) if B._has_pa() else None
    if legs is not None:
        assert {"eta", "phi", "charge", "dq0"} <= set(legs.columns) and legs.dq0.isna().all()
    # bez diff: nominalne histogramy identyczne
    I2, H2, _ = B.cmd_pairs(cfg, unblind=True)
    assert all(np.array_equal(H0[ch], H2[ch]) for ch in CHANNELS) and "diff" not in I2.meta
    # v1.2 jako charge_source → D0 nogi określone; kontrola wierszy
    v12 = str(d / "h5_v12")
    _make_v12(str(d / "h5"), v12)
    cfg12 = C.load(None, {"paths": {"h5_dir": cfg["paths"]["h5_dir"], "out_dir": cfg["paths"]["out_dir"]},
                          "pairs": {"mb_rates_json": cfg["pairs"]["mb_rates_json"], "diff": {"var": "sign", "charge_source": v12}}})
    I3, H3, c3 = B.cmd_pairs(cfg12, unblind=True, diff="cfg")
    Hd3 = dict(np.load(os.path.join(out, "pairs_full_diff.npz")))
    assert all(np.array_equal(H0[ch], H3[ch]) for ch in CHANNELS)
    assert sum(Hd3[DF.hkey("D0D0", k)].sum() for k in range(2)) > 0 and c3["diff_undefined"] < counts["diff_undefined"]
    D3 = json.load(open(os.path.join(out, "pair_inputs_diff.json")))
    assert D3["n_mc_pairs"]["D0D0"][0] + D3["n_mc_pairs"]["D0D0"][1] > 0
    bad = str(d / "h5_v12_bad")
    _make_v12(str(d / "h5"), bad, shuffle=True)
    cfgbad = C.load(None, {"paths": {"h5_dir": cfg["paths"]["h5_dir"], "out_dir": cfg["paths"]["out_dir"]},
                           "pairs": {"mb_rates_json": cfg["pairs"]["mb_rates_json"], "diff": {"var": "sign", "charge_source": {"data_none": os.path.join(bad, "data_428770_d150_Donly_none.h5")}}}})
    with pytest.raises(SystemExit):
        B.cmd_pairs(cfgbad, unblind=True, diff="cfg")
    shutil.rmtree(bad, ignore_errors=True)


def test_fit2d_cli_diff_requires_unblind_and_inputs(synth_diff):
    from charmtag.pairs.cli import main
    cfg, d = synth_diff
    with pytest.raises(SystemExit):
        main(["fit2d", "--diff", "--out-dir", cfg["paths"]["out_dir"]])
    empty = str(d / "empty_out")
    os.makedirs(os.path.join(empty, "pairs"), exist_ok=True)
    with pytest.raises(SystemExit):
        main(["fit2d", "--unblind", "--diff", "--out-dir", empty])


def test_d0_flavour_closest_mass():
    """Znak D0 z hipotezy bliższej m0 (obie hipotezy tej samej pary torów: trk0/trk1 zamienione, dq0 przeciwne); brak drugiego rzędu → własny dq0."""
    L = pd.DataFrame({"source_file": ["1"] * 5, "entry": [1, 1, 2, 2, 3], "trk0": [10, 11, 20, 21, 30], "trk1": [11, 10, 21, 20, 31],
                      "m": [1865.0, 1900.0, 1950.0, 1868.0, 1800.0], "dq0": [-1, 1, 1, -1, np.nan]}, index=[7, 3, 9, 2, 5])
    out = B.d0_flavour_closest_mass(L)
    assert out.dq0.tolist()[:4] == [-1.0, -1.0, -1.0, -1.0] and np.isnan(out.dq0.iloc[4])
    assert list(out.index) == list(L.index) and out.m.tolist() == L.m.tolist()           # nic poza dq0 się nie zmienia
    assert DF.spec_from_cfg({"diff": True})["d0_flavour"] == "closest_mass"
    with pytest.raises(SystemExit):
        DF.spec_from_cfg({"diff": {"d0_flavour": "x"}})
    # make_legs: dq0 nogi po arbitrażu = dq0 hipotezy bliższej m0 (niezależnie od tego, który rząd wygrał arbitraż)
    sc = pd.DataFrame({"row": [0, 1], "run": [1, 1], "event": [5, 5], "m": [1.900, 1.865], "pt": [5.0, 5.0], "score": [0.9, 0.9],
                       "veto_ok": [True, True], "truth": [-1, -1]})
    T = pd.DataFrame({"row": [0, 1], "trk0": [10, 11], "trk1": [11, 10], "trk2": [-1, -1], "eta": [0, 0], "phi": [0, 0], "charge": [0, 0], "dq0": [1.0, -1.0]})
    thr = {"0.9": 0.5, "0.65": 0.6, "0.45": 0.7, "0.3": 0.8}
    E = {10: np.linspace(1750, 2000, 26), 4: np.linspace(1700, 2200, 26)}
    p = {**B.pcfg({}), "diff": DF.spec_from_cfg({"diff": True})}
    legs = B.make_legs(sc, T, 10, thr, p, None, E)
    assert len(legs) == 1 and legs.dq0.iloc[0] == -1.0                                   # hipoteza z m=1865 (bliższa m0) ma dq0=−1
    p2 = {**p, "diff": {**p["diff"], "d0_flavour": "hypothesis"}}
    legs2 = B.make_legs(sc, T, 10, thr, p2, None, E)
    assert legs2.cand.iloc[0] == legs.cand.iloc[0] and legs2.dq0.iloc[0] == T.dq0[legs2.cand.iloc[0]]   # ten sam arbitraż, dq0 własne


def test_charge_files_keeps_sample_subdir(tmp_path):
    """PD (06.11): katalog charge_source ma te same podkatalogi próbek co h5_dir (v1.2 = parts/<próbka>/...);
    wcześniej brano samą nazwę pliku → brak pliku → None po cichu → nogi D0 „nieokreślone"."""
    from charmtag.pairs import build
    h5, v12 = tmp_path / "v11", tmp_path / "v12"
    for base in (h5, v12):
        (base / "data_x").mkdir(parents=True)
        (base / "data_x" / "data_x_1_none.h5").write_bytes(b"")
    cfg = {"paths": {"h5_dir": str(h5), "samples": {"data_none": "data_x/data_x_*_none.h5"}}, "io": {"dataset": "events"}, "columns": {}}
    assert build.charge_files(cfg, {"charge_source": str(v12)}, "data_none") == [str(v12 / "data_x" / "data_x_1_none.h5")]
    (v12 / "data_x" / "data_x_1_none.h5").unlink()
    assert build.charge_files(cfg, {"charge_source": str(v12)}, "data_none") is None
    assert build.charge_files(cfg, {"charge_source": None}, "data_none") is None
