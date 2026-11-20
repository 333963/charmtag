#!/usr/bin/env python3
"""Test §3b (18.11.2026): wykluczenie plików mb_nd (gn_combine.b1_exclude) na danych syntetycznych.
Sprawdza: nogi/D* wykluczonych tagów usunięte; po złączeniu z zdarzeniami 78 plików legs_without_event = 0; referencje B1a = etap1 − wkład
wykluczonych tagów; testy count_events / count_Dplus_n / count_pairs_DpDp_n / B1b_* przeniesione do informacyjnych; tolerancja B1b78.
Uruchom: .venv_h5/bin/python tests/test_exclude.py"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gn_combine as C  # noqa: E402


def make(rng, tags, ev_per_tag=200):
    rows, evs = [], []
    for k, t in enumerate(tags):
        for e in range(ev_per_tag):
            evs.append({"file": k, "tag": t, "evnum": np.uint64(10**6 * int(t) + e)})
            nl = rng.integers(0, 3)
            for i in range(nl):
                rows.append({"tag": t, "evnum": np.uint64(10**6 * int(t) + e), "row": len(rows), "decay_type": int(rng.choice([10, 4])),
                             "trk0": 3 * len(rows), "trk1": 3 * len(rows) + 1, "trk2": 3 * len(rows) + 2 if rng.random() < .5 else -1})
    return pd.DataFrame(rows), pd.DataFrame(evs)


def main():
    rng = np.random.default_rng(3)
    all_tags = [f"{i:06d}" for i in range(1, 9)]
    excl = ["000003", "000007"]
    legs, ev_all = make(rng, all_tags)
    Ds = legs[legs.decay_type == 10].iloc[::4][["tag", "evnum", "row"]].copy()
    ev = ev_all[~ev_all.tag.isin(excl)].reset_index(drop=True)
    ev["file"] = pd.factorize(ev.tag)[0]
    keep = ~legs.tag.isin(excl)
    # „etap1” na 80 plikach ≙ wszystkie 8 tagów; rep = te same definicje na tagach bez wykluczonych
    n_dp_all = int((legs.decay_type == 4).sum()); n_dp_keep = int((legs[keep].decay_type == 4).sum())
    P_all = C.pairs_of_legs(legs); P_keep = C.pairs_of_legs(legs[keep])
    pp_all = int(((P_all.decay_type1 == 4) & (P_all.decay_type2 == 4)).sum()); pp_keep = int(((P_keep.decay_type1 == 4) & (P_keep.decay_type2 == 4)).sum())
    old_rc, old_n = dict(C.REF_COUNTS), C.N_MB_FILES_REF
    C.REF_COUNTS.update({"events": len(ev_all), "Dplus_n": n_dp_all, "pairs_DpDp_n": pp_all}); C.N_MB_FILES_REF = len(all_tags)
    try:
        f = len(excl) / len(all_tags)
        ref = {"events": len(ev_all), "sum_w": float(len(ev_all)),
               "Dplus_correct_per_event": {"unw": n_dp_all / len(ev_all), "w": n_dp_all / len(ev_all), "n": n_dp_all},
               "pairs_DpDp": {"per_event_unw": pp_all / len(ev_all), "per_event_w": pp_all / len(ev_all), "n": pp_all}}
        rep = {"events": len(ev), "sum_w": float(len(ev)),
               "Dplus_correct_per_event": {"unw": n_dp_keep / len(ev), "w": n_dp_keep / len(ev), "n": n_dp_keep},
               "pairs_DpDp": {"per_event_unw": pp_keep / len(ev), "per_event_w": pp_keep / len(ev), "n": pp_keep, "n_eff": float(pp_keep)}}
        checks = {"count_events": {"got": len(ev), "ref": len(ev_all), "pass": False},
                  "count_Dplus_n": {"got": n_dp_keep, "ref": n_dp_all, "pass": False},
                  "count_pairs_DpDp_n": {"got": pp_keep, "ref": pp_all, "pass": False},
                  "B1b_Dplus_correct_per_event_unw": {"pass": False}, "B1b_Dplus_correct_per_event_w": {"pass": False},
                  "B1b_pairs_DpDp_per_event_unw": {"pass": False}, "B1b_pairs_DpDp_per_event_w": {"pass": False},
                  "B1b_anchors_found": {"pass": True}}
        acc = {"checks": checks, "info_nonblocking": {}}
        L2, D2 = C.b1_exclude(acc, rep, ref, legs, Ds, ev, excl, neff={"Dplus_correct_per_event": float(n_dp_keep), "pairs_DpDp": float(pp_keep)})
    finally:
        C.REF_COUNTS.clear(); C.REF_COUNTS.update(old_rc); C.N_MB_FILES_REF = old_n
    assert not L2.tag.isin(excl).any() and not D2.tag.isin(excl).any()
    assert len(L2) == int(keep.sum())
    m = L2.merge(ev[["tag", "evnum"]].assign(hit=1), on=["tag", "evnum"], how="left")
    assert int(m.hit.isna().sum()) == 0, "legs_without_event po wykluczeniu ≠ 0"
    m0 = legs.merge(ev[["tag", "evnum"]].assign(hit=1), on=["tag", "evnum"], how="left")
    assert int(m0.hit.isna().sum()) == int((~keep).sum()), "bez wykluczenia nogi wykluczonych tagów = legs_without_event"
    ch = acc["checks"]
    for k in ("count_events", "count_Dplus_n", "count_pairs_DpDp_n", "B1b_Dplus_correct_per_event_unw", "B1b_pairs_DpDp_per_event_w"):
        assert k not in ch, k
    assert ch["B1a78_Dplus_n"]["ref"] == n_dp_keep and ch["B1a78_Dplus_n"]["pass"], ch["B1a78_Dplus_n"]
    assert ch["B1a78_pairs_DpDp_n"]["ref"] == pp_keep and ch["B1a78_pairs_DpDp_n"]["pass"], ch["B1a78_pairs_DpDp_n"]
    for k in ("B1b78_Dplus_correct_per_event_unw", "B1b78_Dplus_correct_per_event_w", "B1b78_pairs_DpDp_per_event_unw", "B1b78_pairs_DpDp_per_event_w"):
        assert ch[k]["pass"] and ch[k]["rel_diff"] < 1e-9, (k, ch[k])   # wagi = 1 → referencja z odjęciem dokładna
    tol = ch["B1b78_Dplus_correct_per_event_w"]["tol"]
    ex = acc["info_nonblocking"]["exclusion_78"]
    assert ex["events_implied_excluded"] == len(ev_all) - len(ev) and ex["events_implied_in_file_range"], ex
    print(f"test_exclude OK: nogi {len(legs)} → {len(L2)} (wykluczone {ex['legs_excluded']}, {100 * ex['legs_excluded_frac']:.1f} %), "
          f"D+ ref {ch['B1a78_Dplus_n']['ref']}, D+D+ ref {ch['B1a78_pairs_DpDp_n']['ref']}, B1b78 D+ rel {ch['B1b78_Dplus_correct_per_event_unw']['rel_diff']:.3f} < {tol:.3f}: "
          f"{ch['B1b78_Dplus_correct_per_event_unw']['pass']}")


if __name__ == "__main__":
    main()
