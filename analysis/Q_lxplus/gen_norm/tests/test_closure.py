#!/usr/bin/env python3
"""Test jednostkowy gn_combine: dopasowanie nóg ↔ prawda i domknięcie truth×mapy vs direct na zabawce z ε faktoryzowanym
(oczekiwane truth/direct = 1 w granicach statystyki). Uruchom: python3 tests/test_closure.py"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gn_combine as G  # noqa: E402


def toy(seed=3, n_ev=60000):
    rng = np.random.default_rng(seed)
    rows = []
    ntrk = rng.integers(5, 80, n_ev); w2 = rng.gamma(2.0, 0.5, n_ev)
    for e in range(n_ev):
        for _ in range(rng.poisson(1.6)):
            pdg = rng.choice([421, 411, 413], p=[0.45, 0.4, 0.15]) * rng.choice([-1, 1])
            rows.append((f"{e % 10:06d}", np.uint64(e), pdg, rng.uniform(3, 20), rng.uniform(-2.4, 2.4), rng.uniform(-3, 3), ntrk[e], w2[e]))
    T = pd.DataFrame(rows, columns=["tag", "evnum", "pdg", "pt", "eta", "phi", "ntrk", "w2"])
    T["mode"] = 1; T["matched"] = 0
    eps = 0.15 + 0.03 * np.minimum(T.pt.to_numpy(), 10) + 0.002 * T.ntrk.to_numpy()                # ε(pT, n_trk), faktoryzowane
    rec = rng.random(len(T)) < eps
    L = T[rec & (np.abs(T.pdg) != 413)].copy()
    L["decay_type"] = np.where(np.abs(L.pdg) == 421, 10, 4); L["row"] = np.arange(len(L))
    for i in range(3):
        L[f"trk{i}"] = np.arange(len(L)) * 3 + i                                                     # rozłączne ślady
    D = T[rec & (np.abs(T.pdg) == 413)].copy(); D["row"] = np.arange(len(D))
    return T, L.reset_index(drop=True), D.reset_index(drop=True)


def main():
    T, legs, Ds = toy()
    legs["tidx"] = -1
    for dt, pdg in ((10, 421), (4, 411)):
        m = (legs.decay_type == dt).to_numpy()
        legs.loc[m, "tidx"] = G.match_legs(legs[m], T, dt, pdg)
    Ds["tidx"] = G.match_legs(Ds, T, 1, 413, use_phi=True)
    assert (legs.tidx >= 0).all() and (Ds.tidx >= 0).all(), "dopasowanie nie znalazło prawdy"
    T["n_legs"] = 0
    for X in (legs, Ds):
        cnt = X.tidx.value_counts(); T.loc[cnt.index, "n_legs"] += cnt.to_numpy()
    assert int(T.n_legs.sum()) == len(legs) + len(Ds)
    P = G.pairs_of_legs(legs)
    P["ch"] = np.where((P.decay_type1 == 10) & (P.decay_type2 == 10), "D0D0", np.where((P.decay_type1 == 4) & (P.decay_type2 == 4), "DpDp", "DpD0"))
    M = G.build_maps(T)
    full = {sp: {"num_mult": G.eps_from(M, sp, "num_mult"), "num_any": G.eps_from(M, sp, "num_any")} for sp in ("D0", "Dplus", "Dstar")}
    c = G.closure(T, legs, Ds, P, lambda tg: full, "toy")
    ok = True
    for ch, v in c["channels"].items():
        print(ch, f"direct {v['direct_n']} (n_eff {v['direct_n_eff']:.0f})  truth/direct {v['ratio_truth_over_direct_mult']:.3f}  pull {v['pull_mult']:+.2f}")
        ok &= abs(v["pull_mult"]) < 3.5
    for sp, v in c["singles"].items():
        print(sp, "singles truth/direct", round(v["ratio_truth_over_direct_mult"], 4))
        ok &= abs(v["ratio_truth_over_direct_mult"] - 1) < 0.03
    print("TEST", "OK" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
