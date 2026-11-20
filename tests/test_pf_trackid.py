"""PF: ID śladu z (d0, z0) musi łączyć ten sam ślad po refit (inne pT/η), a stary (pT, η) — nie."""
import numpy as np
import pandas as pd
from charmtag.pairs import build
from charmtag.pairs.frozen import shared_track


def _g(pt, eta, d0, z0):
    g = {}
    for i in range(3):
        g[f"d_pt{i}"] = np.asarray(pt[i], np.float32); g[f"d_eta{i}"] = np.asarray(eta[i], np.float32)
        g[f"d_d0{i}"] = np.asarray(d0[i], np.float32); g[f"d_z0{i}"] = np.asarray(z0[i], np.float32)
    return g


def test_d0z0_matches_refit_track_and_pt_eta_does_not():
    # dwa kandydaty D+ (3 córki) dzielące ślad #0: po refit pT różni się o 0,1 %, η o 1e-4; d0/z0 identyczne
    pt = [[1.000, 1.001], [2.0, 3.0], [0.7, 0.9]]
    eta = [[0.50, 0.5001], [1.0, -1.0], [0.1, 0.2]]
    d0 = [[0.0123, 0.0123], [0.05, -0.02], [0.001, 0.002]]
    z0 = [[-0.31, -0.31], [0.2, 0.4], [0.9, 0.8]]
    g = _g(pt, eta, d0, z0)
    new = build.track_ids({"pairs": {"track_id": "d0z0"}}, g)
    old = build.track_ids({"pairs": {"track_id": "pt_eta"}}, g)
    assert new[0][0] == new[0][1]          # ten sam ślad → ten sam ID
    assert old[0][0] != old[0][1]          # stary ID rozjeżdża się po refit
    assert len({new[i][0] for i in range(3)}) == 3  # różne ślady → różne ID
    P = pd.DataFrame({"trk01": [new[0][0]], "trk11": [new[1][0]], "trk21": [new[2][0]],
                      "trk02": [new[0][1]], "trk12": [new[1][1]], "trk22": [new[2][1]]})
    assert shared_track(P).all()
    P = pd.DataFrame({"trk01": [old[0][0]], "trk11": [old[1][0]], "trk21": [old[2][0]],
                      "trk02": [old[0][1]], "trk12": [old[1][1]], "trk22": [old[2][1]]})
    assert not shared_track(P).any()


def test_missing_d0_columns_is_error_not_silent_fallback():
    g = _g([[1.0], [2.0], [0.0]], [[0.1], [0.2], [0.0]], [[0.1], [0.2], [0.0]], [[0.1], [0.2], [0.0]])
    for k in list(g):
        if k.startswith("d_d0"):
            del g[k]
    import pytest
    with pytest.raises(SystemExit):        # recenzja kodu 2.0.1 P1-3: brak d0/z0 → twardy błąd, nie cichy powrót do (pT, η)
        build.track_ids({"pairs": {"track_id": "d0z0"}}, g)


def test_default_is_d0z0():
    assert build.pcfg({})["track_id"] == "d0z0"
