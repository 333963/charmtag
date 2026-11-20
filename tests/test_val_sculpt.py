"""VAL: metryki rzeźbienia rozróżniają gładkie nachylenie ε_B(m) (R_win_window>1, R_win_sym≈1, pol2_closure≈1)
od prawdziwego piku w ε po cięciu (oba >1)."""
import numpy as np
import pandas as pd
from charmtag import sculpt as SC
from charmtag import pipeline as P

SP = {"m0": 1.86966, "mass_range": [1.70, 2.20], "signal_window": 0.045, "sideband_min": 0.070, "inner_band": [0.045, 0.070]}
THR = {"3-30": {"0.80": 0.2, "0.70": 0.3, "0.50": 0.5}}
WPS = [0.80, 0.70, 0.50]


def _bkg(n=400_000, seed=1):
    rng = np.random.default_rng(seed)
    m = rng.uniform(1.70, 2.20, n)                       # tło płaskie w m
    pt = rng.uniform(3, 30, n)
    return m, pt, rng


def test_slope_only_is_not_sculpting():
    m, pt, rng = _bkg()
    # ε(m) maleje liniowo ×2.5 na 0.5 GeV: score = u^(a(m)) → P(score>t) zależy od m gładko
    eps80 = 0.12 * (1 - 1.2 * (m - 1.70))                # 0.12 → 0.048
    u = rng.random(len(m))
    score = np.where(u < eps80, 0.9, 0.0)                 # przechodzi WP80 z prawd. eps80(m)
    d = pd.DataFrame({"m": m, "pt": pt, "score": score})
    old = P.rwin(d, SP, THR, [0.80], "window")["80"]["R"]
    blk = SC.validate_block(m, {0.80: score >= 0.2}, SP, has_signal=False)["80"]
    assert old > 1.10, old                                 # stara metryka „widzi rzeźbienie”
    assert abs(blk["R_win_sym"] - 1) < 3 * blk["R_win_sym_err"] + 0.01
    assert abs(blk["pol2_closure"] - 1) < 3 * blk["pol2_closure_err"] + 0.02
    assert blk["slope_perGeV"] < -0.8
    assert blk["proxy_ratio_SB_asym"] < 0.9 < blk["proxy_ratio_SB_sym"] < 1.1  # asymetryczne SB zaniżają ε_B pod pikiem


def test_real_peak_in_efficiency_is_flagged():
    m, pt, rng = _bkg()
    eps = 0.05 * (1 + 1.5 * np.exp(-0.5 * ((m - SP["m0"]) / 0.02) ** 2))   # pik ε pod m0
    score = np.where(rng.random(len(m)) < eps, 0.9, 0.0)
    d = pd.DataFrame({"m": m, "pt": pt, "score": score})
    old = P.rwin(d, SP, THR, [0.80], "window")["80"]["R"]
    blk = SC.validate_block(m, {0.80: score >= 0.2}, SP, has_signal=False)["80"]
    assert old > 1.3 and blk["R_win_sym"] > 1.3 and blk["pol2_closure"] < 0.8


def test_data_block_inner_band_and_regions():
    m, pt, rng = _bkg(100_000)
    score = np.where(rng.random(len(m)) < 0.1, 0.9, 0.0)
    blk = SC.validate_block(m, {0.80: score >= 0.2}, SP, has_signal=True)["80"]
    assert abs(blk["R_inner_sym"] - 1) < 4 * blk["R_inner_sym_err"]
    assert set(blk["eps"]) >= {"SR", "SB_asym", "SB_sym", "SB_left", "SB_right", "inner_sym"}


def test_validate_hooks_present():
    assert hasattr(P, "_sculpt") and hasattr(P, "_hand_control") and hasattr(P, "_salad_eps")
    assert "R_win_sym" in SC.LEGEND
