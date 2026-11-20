"""19.11: klucz zdarzenia w ev_extract — (run, event) niejednoznaczne w fd_nd → indeks bloku evblk + porcje wyrównane do zdarzeń."""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import ev_extract as E  # noqa: E402


def test_blocks_split_repeated_event_numbers():
    ev = np.array([7, 7, 7, 9, 9, 7, 7, 5])          # zdarzenie 7 dwa razy (różne joby), run = 0
    w = np.array([1., 1., 1., 2., 2., 3., 3., 1.])
    b = E.event_blocks({"run_number": np.zeros(8, int), "event_number": ev, "event_weight": w})
    assert b.tolist() == [0, 0, 0, 1, 1, 2, 2, 3]


def test_blocks_same_number_adjacent_different_weight():
    b = E.event_blocks({"event_number": np.array([4, 4, 4, 4]), "event_weight": np.array([1., 1., 2., 2.])})
    assert b.tolist() == [0, 0, 1, 1]


def test_merge_no_cross_event_pairs():
    # D0 w bloku 0, D* w bloku 2 o tym samym (run, event): stary klucz paruje, nowy nie
    L = pd.DataFrame({"run": [0, 0], "event": [7, 7], "evblk": [0, 2], "m": [1.865, 1.865]})
    D = pd.DataFrame({"run": [0], "event": [7], "evblk": [2], "m0_from": [1.865]})
    assert len(L.merge(D, on=["run", "event"])) == 2
    assert len(L.merge(D, on=["run", "event", "evblk"])) == 1


def test_chunks_do_not_cut_events():
    blk = np.repeat(np.arange(10), 3)                 # 10 zdarzeń × 3 wiersze
    bnd = E.chunk_bounds(blk, 4)
    assert bnd[0][0] == 0 and bnd[-1][1] == len(blk)
    assert all(a2 == b1 for (_, b1), (a2, _) in zip(bnd, bnd[1:]))
    for a, b in bnd:
        assert b == len(blk) or blk[b] != blk[b - 1]


def test_real_fd_file_if_present():
    p = os.path.expanduser("~/Documents/charm/charm_work/FEAT_charmnet/data/fd_nd/fd_nd_000001_none.h5")
    if not os.path.exists(p):
        return
    import h5py
    with h5py.File(p, "r") as f:
        e = f["events"].fields([c for c in E.EV_COLS if c in f["events"].dtype.names])[:]
    blk = E.event_blocks({c: e[c] for c in e.dtype.names})
    nb = int(blk[-1]) + 1
    assert nb > len(np.unique(e["event_number"]))      # więcej zdarzeń niż numerów zdarzeń (powtórzenia)
