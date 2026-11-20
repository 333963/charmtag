"""HYP (2026-11-08): arbiter hipotezy D0→Kπ — toy MC, rzeźbienie, addytywność, integracja (samples/pipeline/pairs), ekstrakcja."""
from __future__ import annotations

import json
import os

import h5py
import numpy as np
import pandas as pd
import pytest

from charmtag import config as C
from charmtag import hyp as H
from charmtag import hyp_toy, store
from charmtag.constants import M_D0


def _prep(L):
    gid, n, rank = H.pair_index(L)
    return L.assign(gid=gid, n=n, rank=rank, twin=H.twin_index(gid, n))


def _train(L, feats="pid", n_folds=3):
    tw = L["twin"].to_numpy()
    hs = L["hypo_swap"].to_numpy()
    ok = (tw >= 0) & np.isin(hs, (0, 1))
    ok &= np.where(tw >= 0, hs[np.maximum(tw, 0)] + hs == 1, False)
    a = H.Arbiter(feats, n_folds)
    a.fit(L[ok], (hs[ok] == 0).astype(int), n_estimators=200)
    return a, ok


def _eps(L, keep, sel):
    c = sel & (L["hypo_swap"].to_numpy() == 0)
    return keep[c].mean(), c.sum()


@pytest.fixture(scope="module")
def toy():
    L = _prep(hyp_toy.make("fd_nd", 12000, 0, 11))
    a, ok = _train(L)
    p = a.predict(L)
    return L, a, ok, p


def test_pairs_and_antisymmetry(toy):
    L, a, ok, p = toy
    assert set(np.unique(L["n"])) <= {1, 2}
    tw = L["twin"].to_numpy()
    i = np.flatnonzero(tw >= 0)
    np.testing.assert_allclose(p[i] + p[tw[i]], 1.0, atol=1e-5)      # P(A) + P(B) = 1 (antysymetria)
    k = H.choose(L, "arbiter", p, gid=L.gid.values, n=L.n.values, rank=L["rank"].values)
    assert k[L["n"].to_numpy() == 1].all()                            # pojedyncze bez zmian
    assert (pd.Series(k[i]).groupby(L["gid"].to_numpy()[i]).sum() == 1).all()   # dokładnie 1 wiersz na parę


def test_arbiter_beats_hash(toy):
    L, a, ok, p = toy
    g = dict(gid=L.gid.values, n=L.n.values, rank=L["rank"].values)
    e_arb, n = _eps(L, H.choose(L, "arbiter", p, **g), ok)
    e_hash, _ = _eps(L, H.choose(L, "hash", **g), ok)
    assert abs(e_hash - 0.5) < 4 * 0.5 / np.sqrt(n)
    assert e_arb > e_hash + 4 * 0.5 / np.sqrt(n)


def test_arbiter_beats_closest_mass_at_low_p():
    """Reżim, w którym dE/dx rozdziela K/π (p_min < 0,7 GeV, centralnie), a rozdzielczość masy jest słaba: arbiter > closest_mass."""
    L = _prep(hyp_toy.make("fd_nd", 15000, 0, 3, res=0.05, pt_lo=0.3, d_pt_min=0.1, eta_max=0.6))
    a, ok = _train(L)
    p = a.predict(L)
    g = dict(gid=L.gid.values, n=L.n.values, rank=L["rank"].values)
    pmin = np.minimum(H.track_p(L, 0), H.track_p(L, 1))
    sel = ok & (pmin < 0.7)
    e_arb, n = _eps(L, H.choose(L, "arbiter", p, **g), sel)
    e_cm, _ = _eps(L, H.choose(L, "closest_mass", **g), sel)
    assert n > 400
    assert e_arb > e_cm + 3 * np.sqrt(0.25 / n)


def test_no_sculpting_pid_but_kin_and_closest_mass_sculpt(toy):
    L, a, ok, p = toy
    S = _prep(hyp_toy.make("fd_nd", 12000, 150000, 12))
    a_kin, _ = _train(S, "pid_kin")
    B = S[S["truth"].to_numpy() >= 8].reset_index(drop=True)
    B = _prep(B)
    B = B[B["n"] == 2].reset_index(drop=True)
    B = _prep(B)
    g = dict(gid=B.gid.values, n=B.n.values, rank=B["rank"].values)
    m = B["m"].to_numpy()
    z = {}
    for name, k in (("arbiter", H.choose(B, "arbiter", a.predict(B), **g)), ("hash", H.choose(B, "hash", **g)),
                    ("closest_mass", H.choose(B, "closest_mass", **g)),
                    ("kin", H.choose(B, "arbiter", a_kin.predict(B, twin=B["twin"].to_numpy()), **g))):
        z[name] = H.sculpt_test(m, k, B.gid.values, B.n.values, M_D0)["peak"]["z"]
    assert abs(z["arbiter"]) < 3 and abs(z["hash"]) < 3
    assert z["closest_mass"] > 20
    assert z["kin"] > 3, z          # (p0, p1, α) wyznaczają masę → pid_kin rzeźbi (dlatego wyłączony w charmtag)


def test_hash_choice_blind_and_deterministic(toy):
    L, *_ = toy
    k1 = H.choose(L, "hash")
    k2 = H.choose(L.iloc[::-1].reset_index(drop=True), "hash")[::-1]
    assert (k1 == k2).all()


def test_additive_defaults():
    cfg = C.load()
    assert cfg["signal"]["swap_resolution"] == "closest_mass"
    from charmtag.pairs.build import DEFAULTS, d0_arbitrate
    assert DEFAULTS["d0_arbitration"] == "hash"
    df = pd.DataFrame({"trk0": [1, 2, 5], "trk1": [2, 1, 6], "source_file": "1", "entry": 7, "cand": [10, 11, 12], "p_hyp": [0.2, 0.8, 0.5]})
    from charmtag.pairs.build import arbitrate
    pd.testing.assert_frame_equal(d0_arbitrate(df, DEFAULTS), arbitrate(df))
    out = d0_arbitrate(df, {**DEFAULTS, "d0_arbitration": "arbiter"})
    assert sorted(out.cand) == [11, 12]
    with pytest.raises(SystemExit):
        d0_arbitrate(df.drop(columns="p_hyp"), {**DEFAULTS, "d0_arbitration": "arbiter"})


def _toy_h5(tmp, sample="fd_nd", ns=1500, nb=1500, seed=5):
    L = hyp_toy.make(sample, ns, nb, seed)
    os.makedirs(os.path.join(tmp, "h5"), exist_ok=True)
    p = os.path.join(tmp, "h5", f"{sample}_000001_none.h5")
    o = hyp_toy.write_h5(L, p, 80, seed)
    return p, o


def _cfg(tmp, p, **sig):
    return C.load(overrides={"paths": {"h5_dir": os.path.dirname(p), "out_dir": os.path.join(tmp, "out"),
                                       "samples": {"signal_mc": os.path.basename(p)}}, "signal": sig})


def test_v12_columns_and_truth_modes(tmp_path):
    from charmtag import samples as SA
    p, o = _toy_h5(str(tmp_path))
    old = pd.concat([x for x, _ in SA.iter_species(_cfg(str(tmp_path), p), "signal_mc", "D0")], ignore_index=True)
    assert {"hypo_swap", "d_trk0", "d_trk1"} <= set(old.columns)              # v1.2 → kolumny opcjonalne obecne
    assert "p_hyp" not in old
    arb = pd.concat([x for x, _ in SA.iter_species(_cfg(str(tmp_path), p, swap_resolution="arbiter"), "signal_mc", "D0")],
                    ignore_index=True)
    assert "p_hyp" not in arb                                                  # bez model_dir → tylko prawda z hypo_swap
    t, hs = arb["truth"].to_numpy(), arb["hypo_swap"].to_numpy()
    assert (t[hs == 1] >= 100).all() and (t[(hs == 0) & (t < 8)] < 100).all()
    # v1.1 (bez kolumn v1.2) → ramka bez nowych kolumn (stara ścieżka bit w bit)
    with h5py.File(p) as f:
        a = f["events"][:]
    keep = [n for n in a.dtype.names if n not in ("Dmeson_hypo_swap",) and "trackId" not in n and "truthPdgId" not in n]
    p11 = os.path.join(str(tmp_path), "h5", "v11.h5")
    with h5py.File(p11, "w") as f:
        f.create_dataset("events", data=a[keep])
    v11 = pd.concat([x for x, _ in SA.iter_species(_cfg(str(tmp_path), p11), "signal_mc", "D0")], ignore_index=True)
    assert not ({"hypo_swap", "d_trk0"} & set(v11.columns))
    assert set(v11.columns) == set(old.columns) - {"hypo_swap", "d_trk0", "d_trk1"}


def test_iter_species_with_arbiter_model(tmp_path):
    from charmtag import samples as SA
    p, o = _toy_h5(str(tmp_path), ns=3000, nb=500)
    L = _prep(o.reset_index(drop=True))
    a, ok = _train(L)
    md = os.path.join(str(tmp_path), "arb")
    a.save(md)
    cfg = _cfg(str(tmp_path), p, swap_resolution="arbiter", arbiter={"model_dir": md})
    X = pd.concat([x for x, _ in SA.iter_species(cfg, "signal_mc", "D0")], ignore_index=True)
    assert {"p_hyp", "hyp_keep"} <= set(X.columns)
    gid, n, _ = H.pair_index(X, None)
    s = pd.Series(X["hyp_keep"].to_numpy()).groupby(gid).sum()
    assert (s == 1).all()
    a2 = H.Arbiter.load(md)
    np.testing.assert_allclose(a2.predict(X), X["p_hyp"].to_numpy(), atol=1e-6)


def test_fit_arbiter_requires_hyp_keep(tmp_path):
    from charmtag import pipeline as P
    cfg = C.load(overrides={"paths": {"h5_dir": str(tmp_path), "out_dir": str(tmp_path)},
                            "signal": {"swap_resolution": "arbiter", "arbiter": {"model_dir": "/x"}}})
    md = P.mdir(cfg, "D0", "H")
    os.makedirs(md, exist_ok=True)
    json.dump({"thresholds": {}}, open(os.path.join(md, "meta.json"), "w"))
    sd = P.sdir(cfg, "D0", "H")
    os.makedirs(sd, exist_ok=True)
    store.save(os.path.join(sd, "signal_mc.h5"), pd.DataFrame({"m": [1.86], "pt": [4.0], "score": [0.9], "veto_ok": [True], "truth": [0.0]}))
    with pytest.raises(SystemExit, match="hyp_keep"):
        P.cmd_fit(cfg, "D0", "H")


def test_extract_join_and_guards(tmp_path):
    pytest.importorskip("pyarrow")
    from charmtag import hyp_extract as X
    p, o = _toy_h5(str(tmp_path), ns=800, nb=800)
    with h5py.File(p) as f:
        n = f["events"].shape[0]
    S = pd.DataFrame({"row": o["h5row"].to_numpy(), "m": o["m"].to_numpy(), "score": o["score"].to_numpy(np.float32),
                      "fold": np.zeros(len(o), np.int8), "pt": o["pt"].to_numpy()}).sort_values("row")
    sc = os.path.join(str(tmp_path), "sc.h5")
    store.save(sc, S, {"source": "/eos/a/" + os.path.basename(p), "file_offsets": json.dumps([0, n])})
    out = os.path.join(str(tmp_path), "o", "fd_nd.parquet")
    assert X.main(["--sample", "fd_nd", "--h5", p, "--out", out, "--scores", sc, "--chunk", "700", "--score-mode", "join"]) == 0
    P = pd.read_parquet(out)
    meta = json.load(open(out.replace(".parquet", ".meta.json")))
    assert meta["score_join"]["mass_mismatch"] == 0 and meta["score_join"]["row_missing"] == 0
    ref = o[(o.pt >= 3) & (o.pt < 30)].sort_values("h5row")
    np.testing.assert_array_equal(P["row"].to_numpy(), ref["h5row"].to_numpy())
    np.testing.assert_allclose(P["score"].to_numpy(), ref["score"].to_numpy(), rtol=1e-6)
    assert {"hypo_swap", "d_trk0", "d_dedx1", "dedx_K", "d_tpdg0"} <= set(P.columns)
    # zła masa → join odmawia
    S2 = S.assign(m=S.m + 1e-3)
    store.save(sc, S2, {"source": "/eos/a/" + os.path.basename(p), "file_offsets": json.dumps([0, n])})
    with pytest.raises(SystemExit, match="niezgodnych mas"):
        X.main(["--sample", "fd_nd", "--h5", p, "--out", out, "--scores", sc, "--score-mode", "join"])
    # inne offsety → auto: fallback (bez modeli → score NaN)
    store.save(sc, S, {"source": "/eos/a/" + os.path.basename(p), "file_offsets": json.dumps([0, n + 1])})
    assert X.main(["--sample", "fd_nd", "--h5", p, "--out", out, "--scores", sc]) == 0
    assert json.load(open(out.replace(".parquet", ".meta.json")))["score_mode"] == "none"
    # podpróbka danych: całe zdarzenia
    out2 = os.path.join(str(tmp_path), "o", "data_428770.parquet")
    pd_, od = _toy_h5(str(tmp_path), "data_428770", 400, 3000, 9)
    assert X.main(["--sample", "data_428770", "--h5", pd_, "--out", out2, "--frac", "0.3", "--score-mode", "none"]) == 0
    D = pd.read_parquet(out2)
    allev = od[(od.pt >= 3) & (od.pt < 30)]
    assert 0.2 < D["event"].nunique() / allev["event"].nunique() < 0.4
    kept = set(D["event"])
    assert len(allev[allev.event.isin(kept)]) == len(D)
    # fd: te same event_number/trackId w dwóch plikach (jak v1.1/v1.2 fd) → pary rozdzielone blokami zdarzeń (evblk)
    import shutil
    p2 = p.replace("_000001_", "_000002_")
    shutil.copy(p, p2)
    out3 = os.path.join(str(tmp_path), "o", "fd2.parquet")
    assert X.main(["--sample", "fd_nd", "--h5", os.path.join(os.path.dirname(p), "fd_nd_*_none.h5"), "--out", out3,
                   "--score-mode", "none", "--chunk", "500"]) == 0
    P2 = pd.read_parquet(out3)
    assert len(P2) == 2 * len(P) and set(P2["file"]) == {0, 1}
    gid, n, _ = H.pair_index(P2)
    assert n.max() == 2 and (n == 2).sum() == 2 * (H.pair_index(P)[1] == 2).sum()


@pytest.mark.slow
def test_hyp_analysis_toy_end_to_end(tmp_path):
    import importlib.util
    here = os.path.dirname(__file__)
    f = os.path.join(here, "..", "..", "out", "hyp_analysis.py")
    if not os.path.exists(f):
        pytest.skip("hyp_analysis.py poza repo")
    spec = importlib.util.spec_from_file_location("hyp_analysis", f)
    A = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(A)
    R = A.run(None, str(tmp_path), toy=True, toy_scale=0.2)
    assert R["sculpt"]["arbiter"]["test"]["peak"]["z"] < 3
    assert R["eps"]["arbiter"]["all"]["eps"] > R["eps"]["hash"]["all"]["eps"]


def test_pertrack_v2_antisym_beats_hash_blind_and_roundtrip(tmp_path):
    """Arbiter v2 (T1, mode="pertrack"): PID per tor z centrowaniem na pionach — antysymetria, ε > hash, ślepota na masę
    z konstrukcji (para ππ: f(niższa m) ≈ 0,5), zapis/odczyt."""
    pytest.importorskip("xgboost")
    L = _prep(hyp_toy.make("fd_nd", 15000, 60000, 21, res=0.05, pt_lo=0.3, d_pt_min=0.1, eta_max=0.6))   # reżim, gdzie dE/dx rozdziela
    a = H.Arbiter("pid", 3, mode="pertrack")
    p = a.fit(L, None)
    tw = L["twin"].to_numpy()
    i = np.flatnonzero(tw >= 0)
    np.testing.assert_allclose(p[i] + p[tw[i]], 1.0, atol=1e-6)
    g = dict(gid=L.gid.values, n=L.n.values, rank=L["rank"].values)
    hs = L["hypo_swap"].to_numpy()
    ok = (tw >= 0) & np.isin(hs, (0, 1)) & (L["truth"].to_numpy() < 8)
    e_arb, n = _eps(L, H.choose(L, "arbiter", p, **g), ok)
    assert e_arb > 0.5 + 4 * 0.5 / np.sqrt(n), e_arb
    # tło ππ: wybór niezależny od masy
    pp = (L["d_tpdg0"].to_numpy() == 211) & (L["d_tpdg1"].to_numpy() == 211) & (tw >= 0)
    j = np.flatnonzero(pp & (np.arange(len(L)) < np.maximum(tw, 0)))
    m = L["m"].to_numpy()
    lower = np.where(p[j] >= 0.5, m[j] <= m[tw[j]], m[tw[j]] < m[j])
    assert len(j) > 2000 and abs(lower.mean() - 0.5) < 4 * 0.5 / np.sqrt(len(j)), (lower.mean(), len(j))
    a.save(str(tmp_path / "v2"))
    b = H.Arbiter.load(str(tmp_path / "v2"))
    assert b.mode == "pertrack"
    np.testing.assert_allclose(b.predict(L), p, atol=1e-6)
    old = H.Arbiter("pid", 3)                                     # zgodność wstecz: domyślnie mode="pairs"
    assert old.mode == "pairs"
