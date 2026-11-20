"""T4: tag D*⁺→D0π⁺_soft w ekstrakcji HYP (charmtag.dstar_tag + hyp_extract) — dane syntetyczne v1.2 (NIE fizyka)."""
import json

import h5py
import numpy as np
import pytest

from charmtag import dstar_tag as DT
from charmtag.constants import DM_DSTAR, M_D0

S = -999  # sentinel int v1.2


# ----------------------------------------------------------------------------- jednostkowe
def _T(rows):
    """rows: (evblk, t0, t1, t2, dm, qs, q0, q1, row)"""
    a = np.array(rows, dtype=float).T
    return DT.dstar_table(a[0], a[8], a[1], a[2], a[3], a[4], a[5], a[6], a[7], pis_pt=np.full(len(rows), 0.4))


def test_tag_basic_both_hypotheses():
    # D*+ (π_s +) z torami 5 (K⁻) i 7 (π⁺) → D0 → K⁻π⁺: wiersz z córką 0 = tor 5 (q=−1) zgodny, wiersz 7-jako-K refleksja
    T = _T([(3, 5, 7, 9, 0.1454, +1, -1, +1, 100)])
    ev = np.array([3, 3, 3, 4])
    t0 = np.array([5, 7, 5, 5])
    t1 = np.array([7, 5, 8, 7])
    q0 = np.array([-1, +1, -1, -1])
    r = DT.tag_d0(ev, t0, t1, T, q0)
    assert list(r["dstar_tag"]) == [1, 1, 0, 0]          # inna para torów / inne zdarzenie → brak tagu
    assert list(r["hypo_from_tag"]) == [0, 1, -1, -1]
    assert list(r["soft_pi_charge"]) == [1, 1, 0, 0]
    assert np.allclose(r["dm"][:2], 0.1454) and np.isnan(r["dm"][2:]).all()
    assert list(r["dstar_row"]) == [100, 100, -1, -1]
    assert list(r["dstar_n"]) == [1, 1, 0, 0]


def test_tag_dstar_minus_and_fallback_without_d0_charge():
    # D*− (π_s −): K⁺ = tor o ładunku +1. Wiersz D0 bez ładunków (sentinel) → tor K z ładunków córek D*
    T = _T([(1, 11, 12, 13, 0.1460, -1, -1, +1, 7)])   # D*: córka0 = 11 (q −1), córka1 = 12 (q +1) → K⁺ = 12
    r = DT.tag_d0(np.array([1, 1]), np.array([12, 11]), np.array([11, 12]), T, np.array([S, S]))
    assert list(r["hypo_from_tag"]) == [0, 1]
    r2 = DT.tag_d0(np.array([1, 1]), np.array([12, 11]), np.array([11, 12]), T, None)
    assert list(r2["hypo_from_tag"]) == [0, 1]
    # bez ładunków nigdzie → konwencja charmpp: córka 0 D* = K
    T0 = _T([(1, 11, 12, 13, 0.1460, -1, np.nan, np.nan, 7)])
    r3 = DT.tag_d0(np.array([1, 1]), np.array([12, 11]), np.array([11, 12]), T0, None)
    assert list(r3["hypo_from_tag"]) == [1, 0]


def test_tag_multiple_dstar_conflict_and_best_dm():
    T = _T([(2, 5, 7, 9, 0.1600, -1, +1, -1, 50),     # daleko od Δm0, π_s −
            (2, 5, 7, 10, 0.1455, +1, -1, +1, 51)])   # najbliżej 145,43 MeV, π_s +
    r = DT.tag_d0(np.array([2]), np.array([5]), np.array([7]), T, np.array([-1]))
    assert r["dstar_row"][0] == 51 and r["soft_pi_charge"][0] == 1 and r["hypo_from_tag"][0] == 0
    assert r["dstar_n"][0] == 2 and r["dstar_conflict"][0] == 1


def test_table_cuts():
    T = _T([(1, 5, 7, 9, 0.17, 1, -1, 1, 1),        # Δm ≥ dm_max
            (1, 5, -1, 9, 0.145, 1, -1, 1, 2),       # brak trackId
            (1, 5, 7, 9, -0.999, 1, -1, 1, 3),       # sentinel Δm
            (1, 5, 7, 9, 0.145, S, -1, 1, 4),        # ładunek π_s nieznany
            (1, 5, 7, 9, 0.145, 1, S, 1, 5)])        # q0 sentinel → NaN, wiersz zostaje
    assert list(T["ds_row"]) == [5] and np.isnan(T["ds_q0"].iloc[0])


def test_carry_and_confusion():
    T = _T([(1, 5, 7, 9, 0.145, 1, -1, 1, 1), (2, 5, 7, 9, 0.145, 1, -1, 1, 2)])
    assert list(DT.carry_last_event(T, 2)["ds_row"]) == [2]
    c = DT.confusion(np.array([0, 1, 1, -1, 0]), np.array([0, 1, 0, 0, -1]))
    assert c == {"tag0_true0": 1, "tag0_true1": 0, "tag1_true0": 1, "tag1_true1": 1}


# ----------------------------------------------------------------------------- end-to-end hyp_extract
COLS_F = ["Dmeson_mass", "Dmeson_pt", "Dmeson_eta", "Dmeson_phi", "Dmeson_delta_mass", "Dmeson_Lxy", "Dmeson_chi2", "Dmeson_impact",
          "Dmeson_impactSignificance", "Dmeson_z0SinTheta", "Dmeson_cosThetaStar", "Dmeson_cosThetaP", "Dmeson_ptcone40",
          "Dmeson_mKpi1", "Dmeson_mKpi2", "event_weight"]
PER_F = ["pt", "eta", "phi", "d0", "z0SinThetaPV", "dEdx"]
PER_I = ["npts_dEdx", "charge", "pdgId", "trackId", "truthPdgId", "truthParentPdgId"]


def _dtype():
    dt = [("run_number", "i4"), ("event_number", "i8"), ("Dmeson_decayType", "i4"), ("Dmeson_charge", "i4"),
          ("Dmeson_truth_category", "i4"), ("Dmeson_hypo_swap", "i4"), ("Dmeson_truthPdgId", "i4"), ("Dmeson_n_daughter_truth_ok", "i4")]
    dt += [(c, "f4") for c in COLS_F]
    for i in range(3):
        dt += [(f"Dmeson_daughter{i}_{q}", "f4") for q in PER_F] + [(f"Dmeson_daughter{i}_{q}", "i4") for q in PER_I]
    return np.dtype(dt)


def _row(rng, ev, dt, trk, q, m, dm=-0.999, hs=-1, truth=0, flav=1):
    """trk/q: trackId i ładunek córek (lista 2 lub 3). flav: +1 D0 / −1 D0bar (prawda)."""
    r = {c: 0.0 for c in COLS_F}
    pt = rng.uniform(4, 10)
    r.update({"run_number": 428770, "event_number": ev, "Dmeson_decayType": dt, "Dmeson_truth_category": truth, "Dmeson_hypo_swap": hs,
              "Dmeson_mass": m, "Dmeson_pt": pt, "Dmeson_eta": rng.uniform(-2, 2), "Dmeson_phi": rng.uniform(-3, 3),
              "Dmeson_delta_mass": dm, "Dmeson_Lxy": 1.0, "Dmeson_chi2": 1.0, "Dmeson_impact": 0.01, "Dmeson_impactSignificance": 0.5,
              "Dmeson_cosThetaStar": 0.1, "Dmeson_cosThetaP": 0.999, "Dmeson_ptcone40": 1.0, "Dmeson_mKpi1": 1000, "Dmeson_mKpi2": 1000,
              "event_weight": 1.0, "Dmeson_charge": int(q[-1]) if dt == 1 else 1, "Dmeson_truthPdgId": 421 * flav,
              "Dmeson_n_daughter_truth_ok": 2})
    for i in range(3):
        on = i < len(trk)
        r.update({f"Dmeson_daughter{i}_pt": pt * 400 if on else -999, f"Dmeson_daughter{i}_eta": r["Dmeson_eta"] if on else -999,
                  f"Dmeson_daughter{i}_phi": r["Dmeson_phi"] if on else -999, f"Dmeson_daughter{i}_d0": 0.02,
                  f"Dmeson_daughter{i}_z0SinThetaPV": 0.01, f"Dmeson_daughter{i}_dEdx": 1.5 if on else -999,
                  f"Dmeson_daughter{i}_npts_dEdx": 3 if on else -1, f"Dmeson_daughter{i}_charge": int(q[i]) if on else S,
                  f"Dmeson_daughter{i}_pdgId": (321 if i == 0 else 211) * int(q[i]) if on else S,
                  f"Dmeson_daughter{i}_trackId": int(trk[i]) if on else S,
                  f"Dmeson_daughter{i}_truthPdgId": S, f"Dmeson_daughter{i}_truthParentPdgId": S})
    return r


def _make_events(n_ev=60, seed=3):
    """Zdarzenia: (a) D* + D0 obu hipotez (poprawne zapachy, prawda hypo_swap), (b) D0 bez D*, co 7. — D* z błędnym Δm (poza oknem)."""
    rng = np.random.default_rng(seed)
    rows, truth = [], []
    for ev in range(n_ev):
        flav = 1 if rng.random() < 0.5 else -1          # D0 → K⁻π⁺ (π_s +) / D0bar → K⁺π⁻ (π_s −)
        tk, tp, ts = 10 + ev, 20 + ev, 30 + ev           # tor K, tor π, π_soft
        qk, qp, qs = -flav, flav, flav
        has_ds = ev % 3 != 2
        dm = 0.1454 if ev % 7 else 0.2
        if has_ds:
            rows.append(_row(rng, ev, 1, [tk, tp, ts], [qk, qp, qs], 2.01, dm=dm, hs=0, flav=flav))
        good = rng.normal(M_D0, 0.02)
        refl = good + rng.normal(0, 0.06)
        # wiersz poprawnej hipotezy (córka 0 = tor K) i refleksji (córka 0 = tor π); kolejność losowa
        pair = [_row(rng, ev, 10, [tk, tp], [qk, qp], good, hs=0, flav=flav), _row(rng, ev, 10, [tp, tk], [qp, qk], refl, hs=1, flav=flav)]
        if rng.random() < 0.5:
            pair = pair[::-1]
        rows += pair
        rows.append(_row(rng, ev, 4, [tk, tp, ts], [qk, qp, qs], 1.87))  # D+ (inny gatunek) — ignorowany
        truth.append((ev, has_ds and ev % 7 != 0))
    a = np.zeros(len(rows), _dtype())
    for j, r in enumerate(rows):
        for k, v in r.items():
            a[k][j] = v
    return a, dict(truth)


@pytest.fixture(scope="module")
def extracted(tmp_path_factory):
    pq = pytest.importorskip("pyarrow.parquet")
    from charmtag import hyp_extract as HX
    d = tmp_path_factory.mktemp("t4")
    a, has_tag = _make_events()
    ev = a["event_number"]                  # granica plików na granicy zdarzenia (części h5 = całe zdarzenia; plik = osobny blok)
    k = int(np.flatnonzero(ev != ev[len(a) // 2])[np.flatnonzero(ev != ev[len(a) // 2]) > len(a) // 2][0])
    paths = []
    for j, (lo, hi) in enumerate([(0, k), (k, len(a))]):
        p = d / f"fd_nd_{j:06d}_none.h5"
        with h5py.File(p, "w") as f:
            f.create_dataset("events", data=a[lo:hi])
        paths.append(str(p))
    out = str(d / "fd_nd.parquet")
    # chunk 5 < liczba wierszy zdarzenia (4) + przesunięcia → zdarzenia przecinane granicą porcji (test przeniesienia D*)
    rc = HX.main(["--sample", "fd_nd", "--h5", str(d / "fd_nd_*_none.h5"), "--out", out, "--score-mode", "none", "--chunk", "5"])
    assert rc == 0
    P = pq.read_table(out).to_pandas()
    meta = json.load(open(str(d / "fd_nd.meta.json")))
    return P, meta, has_tag


def test_e2e_columns(extracted):
    P, meta, _ = extracted
    for c in DT.TAG_COLS:
        assert c in P.columns, c
    assert P["dstar_row"].dtype == np.int64 and P["hypo_from_tag"].dtype == np.int8
    assert meta["dstar"]["enabled"] is True


def test_e2e_tag_matches_truth(extracted):
    P, meta, has_tag = extracted
    t = P[P["dstar_tag"] == 1]
    # każdy otagowany wiersz: hypo_from_tag == prawda hypo_swap (obie hipotezy), Δm w oknie
    assert len(t) > 0 and (t["hypo_from_tag"].to_numpy() == t["hypo_swap"].to_numpy()).all()
    assert ((t["dm"] - DM_DSTAR).abs() < 0.002).all()
    # zbiór otagowanych zdarzeń = zdarzenia z D* w oknie Δm (w zakresie masy/pT); zero tagów poza nimi
    ev_tag = set(t["event"].unique())
    ev_in = set(P["event"].unique())
    assert ev_tag == {e for e in ev_in if has_tag[e]}
    n = P[P["dstar_tag"] == 0]
    assert (n["hypo_from_tag"] == -1).all() and n["dm"].isna().all() and (n["soft_pi_charge"] == 0).all()
    # ładunek π_s = −ładunek toru K (wiersz zgodny: córka 0 = K)
    ok = t[t["hypo_from_tag"] == 0]
    assert (ok["soft_pi_charge"].to_numpy() == -ok["d_charge0"].to_numpy()).all()
    c = meta["dstar"]["confusion_vs_hypo_swap"]
    assert c["tag0_true1"] == 0 and c["tag1_true0"] == 0 and c["tag0_true0"] + c["tag1_true1"] == len(t)
    assert meta["dstar"]["counts"]["n_rows_tagged"] == len(t)


def test_e2e_chunk_invariance(extracted, tmp_path):
    """Duża porcja (bez cięcia zdarzeń) daje te same tagi co porcja 5."""
    pq = pytest.importorskip("pyarrow.parquet")
    from charmtag import hyp_extract as HX
    P, _, _ = extracted
    a, _ = _make_events()
    p = tmp_path / "fd_nd_000000_none.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset("events", data=a)
    out = str(tmp_path / "x.parquet")
    for ch in (5, 1000):
        HX.main(["--sample", "fd_nd", "--h5", str(p), "--out", out, "--score-mode", "none", "--chunk", str(ch)])
        Q = pq.read_table(out).to_pandas().sort_values("row").reset_index(drop=True)
        if ch == 5:
            Q5 = Q
    for c in ("dstar_tag", "hypo_from_tag", "soft_pi_charge", "dstar_row"):
        assert np.array_equal(Q5[c].to_numpy(), Q[c].to_numpy()), c


def test_e2e_off_v11(tmp_path):
    """--dstar-tag off (i v1.1 bez trackId) → kolumny są, wszystkie puste."""
    pq = pytest.importorskip("pyarrow.parquet")
    from charmtag import hyp_extract as HX
    a, _ = _make_events(20)
    p = tmp_path / "fd_nd_000000_none.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset("events", data=a)
    out = str(tmp_path / "y.parquet")
    HX.main(["--sample", "fd_nd", "--h5", str(p), "--out", out, "--score-mode", "none", "--dstar-tag", "off"])
    Q = pq.read_table(out).to_pandas()
    assert (Q["dstar_tag"] == 0).all() and (Q["hypo_from_tag"] == -1).all()


def test_e2e_dstar_only(tmp_path):
    pq = pytest.importorskip("pyarrow.parquet")
    from charmtag import hyp_extract as HX
    a, _ = _make_events(30)
    p = tmp_path / "data_428770_000000_none.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset("events", data=a)
    full, only = str(tmp_path / "f.parquet"), str(tmp_path / "o.parquet")
    HX.main(["--sample", "data_428770", "--h5", str(p), "--out", full, "--score-mode", "none", "--chunk", "7"])
    HX.main(["--sample", "data_428770", "--h5", str(p), "--out", only, "--score-mode", "none", "--chunk", "7", "--dstar-only",
             "--max-gb", "0.000001"])
    F = pq.read_table(full).to_pandas()
    O = pq.read_table(only).to_pandas()
    assert len(O) == int((F["dstar_tag"] == 1).sum()) > 0 and (O["dstar_tag"] == 1).all()
    assert np.array_equal(np.sort(O["row"].to_numpy()), np.sort(F.loc[F["dstar_tag"] == 1, "row"].to_numpy()))
    meta = json.load(open(str(tmp_path / "o.meta.json")))
    assert meta["dstar"]["only_tagged"] is True and meta["frac"] == 1.0
