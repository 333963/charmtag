"""PF3: rozbicie par ze wspólnym śladem (liczba śladów, K/π) liczone przed wetem, per kanał i region."""
import pandas as pd
from charmtag.pairs.hists import shared_track_detail, _shared_detail


def _P():
    # 4 pary: D0D0 wspólny K (jeden ślad), D+D0 wspólny π z K (jeden), D+D+ dwa wspólne (K i π), D0D0 bez wspólnych
    return pd.DataFrame({
        "decay_type1": [10, 4, 4, 10], "decay_type2": [10, 10, 4, 10],
        "trk01": [1, 10, 20, 30], "trk11": [2, 11, 21, 31], "trk21": [-1, 12, 22, -1],
        "trk02": [1, 13, 20, 40], "trk12": [3, 10, 21, 41], "trk22": [-1, -1, 23, -1],
        "win1": [False, True, False, False], "win2": [False, False, False, False]})


def test_detail_counts():
    P = _P(); n, kk, kpi, pipi = shared_track_detail(P)
    assert list(n) == [1, 1, 2, 0]
    assert list(kk) == [True, False, True, False]
    assert list(kpi) == [False, True, False, False]
    assert list(pipi) == [False, False, True, False]


def test_accumulate_by_channel_region():
    P = _P(); sh = shared_track_detail(P)[0] > 0; counts = {}
    _shared_detail(P, sh, counts)
    d = counts["shared_detail"]
    assert d["D0D0"]["sbsb"] == {"n_pairs": 2, "n_shared": 1, "shared_1track": 1, "shared_2tracks": 0, "shared_K_K": 1, "shared_K_pi": 0, "shared_pi_pi": 0}
    assert d["DpD0"]["arms"]["n_shared"] == 1 and d["DpD0"]["arms"]["shared_K_pi"] == 1
    assert d["DpDp"]["sbsb"]["shared_2tracks"] == 1 and d["DpDp"]["sbsb"]["shared_pi_pi"] == 1
    _shared_detail(P, sh, counts)  # idempotentna akumulacja
    assert counts["shared_detail"]["D0D0"]["sbsb"]["n_pairs"] == 4
