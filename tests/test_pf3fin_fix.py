"""PF3-FIN: ID śladu ≥ 0 także dla d0 < 0 (bug `<< 32`), oraz ślepe liczniki bez centrum."""
import numpy as np
import pandas as pd
from charmtag.pairs import build
from charmtag.pairs.frozen import shared_track


def test_negative_d0_track_is_matched():
    pt = np.array([1.0, 1.001], np.float32)
    d0 = np.array([-0.0123, -0.0123], np.float32)      # ujemne d0 → bit znaku float32
    z0 = np.array([-0.31, -0.31], np.float32)
    t = build._trk_id_d0z0(pt, d0, z0)
    assert (t >= 0).all() and t[0] == t[1]
    P = pd.DataFrame({"trk01": [t[0]], "trk11": [5], "trk21": [-1], "trk02": [t[1]], "trk12": [7], "trk22": [-1]})
    assert shared_track(P).all()


def test_sign_of_d0_still_distinguishes_from_other_tracks():
    pt = np.ones(4, np.float32)
    d0 = np.array([0.01, -0.01, 0.02, -0.02], np.float32)
    z0 = np.array([0.5, 0.6, 0.5, 0.7], np.float32)
    t = build._trk_id_d0z0(pt, d0, z0)
    assert len(set(t.tolist())) == 4 and (t >= 0).all()


def test_missing_daughter_is_minus_one():
    t = build._trk_id_d0z0(np.array([0.0], np.float32), np.array([-0.1], np.float32), np.array([0.2], np.float32))
    assert t[0] == -1
