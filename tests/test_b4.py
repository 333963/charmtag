"""B4 — testy: ε_B pod pikiem (pol2 / SB sym / asym), systematyka bkg_shape_salad, SF_rel z B2.
Szybkie (<5 s); jeden test `slow` z większą statystyką."""

import json
import math
import os

import numpy as np
import pytest

from charmtag import sf_rel, val_epsb
from charmtag.syst import bkg_shape_salad as bss

M0 = 1.8696  # D⁺ [GeV]
R = val_epsb.default_regions("Dplus")


def _toy(n_bkg=400_000, n_sig=3_000, slope=-2.0, seed=0, eps0=0.03, sig_eps=0.8):
    """Tło exp(−2 m) w 1,70–2,20, ε_B(m) = eps0·exp(slope·(m−m0)) (slope −2/GeV → ×2,7 lewy/prawy jak X/S8);
    sygnał Gauss σ 8 MeV z ε_S = sig_eps. Zwraca m, pass, prawdziwe ε_B w SR."""
    rng = np.random.default_rng(seed)
    lo, hi = R["lo"], R["hi"]
    m = lo + (hi - lo) * rng.random(n_bkg * 2)
    m = m[rng.random(m.size) < np.exp(-2.0 * (m - lo))][:n_bkg]
    eps_m = np.clip(eps0 * np.exp(slope * (m - M0)), 1e-4, 0.5)
    p_b = rng.random(m.size) < eps_m
    ms = rng.normal(M0, 0.008, n_sig)
    p_s = rng.random(n_sig) < sig_eps
    sr = np.abs(m - M0) < R["sr_half"]
    eps_true = p_b[sr].mean()
    return np.concatenate([m, ms]), np.concatenate([p_b, p_s]), eps_true


def test_pol2_recovers_eps_b_under_peak_and_proxies_order():
    m, p, eps_true = _toy()
    res = val_epsb.eps_b_sr_pol2(m, p, M0, R["sr_half"], R["lo"], R["hi"])
    assert res["ok"]
    assert abs(res["eps"] / eps_true - 1) < 0.06, (res["eps"], eps_true)
    sym = val_epsb.eps_b_sb_sym(m, p, M0, R["sb_min"], R["sb_max"])["eps"]
    asym = val_epsb.eps_b_sb_asym(m, p, R["lo"], R["left_max"], R["right_min"], R["hi"])["eps"]
    # nachylenie ujemne: lewy SB (1,70–1,80) ma ε wyższe, prawy (1,94–2,20) niższe; asym (wiecej prawego) zaniża,
    # sym ~ zgodne/zawyża lekko (S8: asym 0,78, sym 1,06–1,14 przy slope −3)
    assert asym < eps_true < 1.2 * eps_true
    assert abs(sym / eps_true - 1) < abs(asym / eps_true - 1)


def test_nominal_block_fields_and_legacy_note():
    m, p, _ = _toy(n_bkg=150_000, n_sig=1_000, seed=1)
    blk = val_epsb.nominal_block(m, {80: p, 70: p & (np.arange(p.size) % 2 == 0)}, M0, R, method="pol2")
    assert blk["method"] == "pol2" and "legacy" in blk["legacy_note"].lower()
    for wp in ("80", "70"):
        row = blk["wps"][wp]
        assert row["method_used"] == "pol2"
        assert row["nominal"] > 0 and np.isfinite(row["rel_sym_vs_nom"]) and np.isfinite(row["rel_asym_vs_nom"])
        assert "sb_asym_legacy" in row and "sb_sym" in row
    assert blk["wps"]["70"]["nominal"] < blk["wps"]["80"]["nominal"]


def test_salad_method_with_weights_and_fallback():
    m, p, eps_true = _toy(n_bkg=100_000, n_sig=500, seed=2)
    rng = np.random.default_rng(3)
    mc_m = R["lo"] + (R["hi"] - R["lo"]) * rng.random(50_000)
    mc_sr = np.abs(mc_m - M0) < R["sr_half"]
    mc_pass = rng.random(mc_m.size) < eps_true
    w = np.clip(rng.lognormal(0, 0.3, mc_m.size), 0.05, 20)
    blk = val_epsb.nominal_block(
        m, {80: p, 50: p}, M0, R, method="salad", salad={80: {"mc_sr": mc_sr, "mc_passed": mc_pass, "w": w}}
    )
    assert blk["wps"]["80"]["method_used"] == "salad"
    assert abs(blk["wps"]["80"]["nominal"] / eps_true - 1) < 0.1
    assert blk["wps"]["50"]["method_used"] == "pol2" and "skipped" in blk["wps"]["50"]["salad"]


def test_bkg_shape_salad_defaults_and_pairs(tmp_path):
    dp = bss.build("Dplus")
    d0 = bss.build("D0")
    d0t8 = bss.build("D0", d0_truth_min2=False)
    assert dp["per"]["80|8-30"]["rel_unc"] == 0.05 and dp["per"]["50|3-5"]["rel_unc"] == 0.25
    assert d0["per"]["50|8-30"]["rel_unc"] == 0.15 and d0t8["per"]["50|8-30"]["rel_unc"] == 0.20
    assert dp["per"]["80|8-30"]["rel_salad_vs_SB"] is None
    assert abs(bss.default_unc("Dplus", 60) - 0.175) < 1e-9  # interpolacja 70→50
    pb = bss.pairs_block(d0, dp)
    assert math.isclose(pb["channels"]["DplusDplus"]["80|8-30"]["rel_unc_B"], math.sqrt(2) * 0.05)
    assert math.isclose(pb["channels"]["DplusD0"]["50|8-30"]["rel_unc_B"], math.hypot(0.25, 0.15))
    assert pb["channels"]["D0D0"]["80|8-30"]["rel_corr_B"] == 1.0
    path = bss.write(str(tmp_path), "Dplus")
    assert os.path.exists(path)
    with open(path) as f:
        assert json.load(f)["source"] == "bkg_shape_salad"


def test_bkg_shape_salad_reads_s8_json(tmp_path):
    # trzy układy JSON-u S8
    flat = {"rel_salad_vs_SB": {"80|8-30": 1.28, "70|8-30": 1.32}}
    nested = {"rel_salad_vs_SB": {"80": {"8-30": 1.28, "5-8": 1.20}}}
    recs = {"rel_salad_vs_SB": [{"wp": 80, "bin": "8-30", "rel_salad_vs_SB": 1.28}]}
    # prawdziwy układ S8 (pełna lumi 05.11): bin zewnętrzny, forma "−1", plus effB_* (stosunek liczony wprost)
    real = {"effB_salad_SR": {"8-30": {"80": 0.0267, "70": 0.013}, "all": {"80": 0.02}},
            "effB_data_SB": {"8-30": {"80": 0.0267 / 1.28, "70": 0.0111}, "all": {"80": 0.015}},
            "effB_salad_SB": {"8-30": {"80": 0.0267 / 1.28 * 1.16, "70": 0.0111}, "all": {"80": 0.015}},
            "rel_salad_vs_SB": {"8-30": {"80": 0.28, "70": 0.168}, "all": {"80": 0.32}}}
    real_rel_only = {"rel_salad_vs_SB": {"8-30": {"80": 0.28, "5-8": 0.58}, "all": {"80": 0.32}}}
    for i, js in enumerate((flat, nested, recs, real, real_rel_only)):
        p = tmp_path / f"s8_{i}.json"
        p.write_text(json.dumps(js))
        b = bss.build("Dplus", str(p))
        assert b["per"]["80|8-30"]["rel_salad_vs_SB"] == pytest.approx(1.28)
        assert b["per"]["80|8-30"]["method"] == "S8 json"
        assert b["per"]["50|3-5"]["rel_salad_vs_SB"] is None
        assert not any(k.startswith("all") or k.endswith("all") for k in bss.load_s8(str(p)))
        if js is real:  # zamknięcie SB 1,16 > B1 default 0,05 → niepewność z closure; 70: closure 1,0 → B1 default 0,10
            assert b["per"]["80|8-30"]["rel_unc"] == pytest.approx(0.16) and b["per"]["70|8-30"]["rel_unc"] == 0.10
        else:
            assert b["per"]["80|8-30"]["closure_SB"] is None and b["per"]["80|8-30"]["rel_unc"] == 0.05
    assert bss.load_s8(None) == {} and bss.load_s8(str(tmp_path / "nie_ma.json")) == {}


def test_sf_rel_lookup_transfer_and_annotate():
    assert sf_rel.sf_rel("D0", "8-30", 70) == (0.97, 0.05)
    assert sf_rel.sf_rel("D0", (8, 30), "WP50") == (0.94, 0.05)
    assert sf_rel.sf_rel("D0", "3_5", 0.8) == (1.0, 0.15)
    v, u = sf_rel.sf_rel("Dplus", "8-30", 70)
    assert v == 0.97 and u == pytest.approx(math.hypot(0.05, 0.05))
    assert sf_rel.sf_rel("D0", "8-30", 95) == (1.0, 0.0)
    with pytest.raises(KeyError):
        sf_rel.sf_rel("D0", "30-60", 70)
    rows = [{"bin": "8-30", "wp": 70, "eps_S": 0.70}, {"bin": "5-8", "wp": 50, "eps_S": 0.50}, {"bin": "x"}]
    out = sf_rel.annotate_rows(rows, "D0")
    assert out[0]["sf_rel_b2"]["syst_rel"] == 0.05 and out[0]["sf_rel_b2"]["eps_S_syst_abs"] == pytest.approx(0.035)
    assert out[0]["sf_rel_b2"]["central_shift_applied"] is False and out[0]["sf_rel_b2"]["eps_S_nominal"] == 0.70
    assert out[1]["sf_rel_b2"]["syst_rel"] == 0.15
    assert "skipped" in out[2]["sf_rel_b2"]


@pytest.mark.slow
def test_pol2_bias_vs_slope_slow():
    """Nachylenia ε(m) −1…−3 /GeV (×1,6…×4,5 lewy/prawy) i WP50-podobne ε0 — pol2 w ±10 %, asym zawsze zaniża."""
    for slope in (-1.0, -2.0, -3.0):
        m, p, eps_true = _toy(n_bkg=1_200_000, n_sig=5_000, slope=slope, eps0=0.01, seed=int(-slope))
        res = val_epsb.eps_b_sr_pol2(m, p, M0, R["sr_half"], R["lo"], R["hi"])
        asym = val_epsb.eps_b_sb_asym(m, p, R["lo"], R["left_max"], R["right_min"], R["hi"])["eps"]
        assert abs(res["eps"] / eps_true - 1) < 0.10, (slope, res["eps"], eps_true)
        assert asym < eps_true


def test_pipeline_eps_b_nominal_hook_and_report():
    """Hak B4 w pipeline.cmd_validate (`validate.eps_b_nominal`) + wiersz w report.sculpt_section."""
    import pandas as pd

    from charmtag import pipeline as P
    from charmtag import report as RP

    m, p, eps_true = _toy(n_bkg=200_000, n_sig=1_000, seed=5)
    d = pd.DataFrame({"m": m, "pt": np.full(m.size, 6.0), "score": p.astype(float)})
    sp = {"m0": M0, "mass_range": [1.70, 2.20], "signal_window": 0.045, "sideband_min": 0.070}
    orig = P._passed_by_wp
    P._passed_by_wp = lambda dd, thr, wps: {wp: dd["score"].to_numpy() > 0.5 for wp in wps}
    try:
        blk = P._eps_b_nominal({"validate": {"eps_b_nominal": "pol2"}}, "Dplus", d, sp, None, [0.8, 0.7])
    finally:
        P._passed_by_wp = orig
    assert set(blk["wps"]) == {"80", "70"} and blk["method"] == "pol2"
    assert abs(blk["wps"]["80"]["nominal"] / eps_true - 1) < 0.08
    L = RP.sculpt_section({"sculpt_data": None, "sculpt_mcbkg": None, "eps_B_SR": blk}, [0.8, 0.7])
    assert L == []   # sekcja renderuje się tylko razem z blokiem VAL
    L = RP.sculpt_section({"sculpt_data": {"all": {}}, "eps_B_SR": blk}, [0.8, 0.7])
    assert any("nominał B4" in x for x in L) and any(x.startswith("| 80 |") for x in L)
