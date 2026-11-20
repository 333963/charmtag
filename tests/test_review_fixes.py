"""Poprawki po recenzji kodu 2.0.1: ślepy eksport nóg, N_D* z dstar_refit, track_id bez cichego fallbacku, ID 64-bitowe."""
import numpy as np
import pandas as pd
import pytest

from charmtag.pairs import build
from charmtag.pairs.model import M0, SBW


def _legs():
    rows = []
    for dt, ms in ((10, [M0[10], M0[10] + 20, M0[10] - 100, M0[10] + 90]), (4, [M0[4], M0[4] - 80, M0[4] + 10])):
        for k, m in enumerate(ms):
            rows.append({"source_file": "1", "entry": k, "cand": 10 * dt + k, "decay_type": dt, "m": float(m), "cat": 0,
                         "trk0": k, "trk1": k + 100, "trk2": -1 if dt == 10 else k + 200})
    return pd.DataFrame(rows)


def test_blind_export_has_no_window_legs():
    L = _legs()
    sb = build.legs_for_export(L, unblind=False)
    assert len(sb) == 3   # D0: ±100/±90 poza 75; D⁺: −80 poza 70
    assert (np.abs(sb.m - sb.decay_type.map(M0)) > sb.decay_type.map(SBW)).all()
    assert len(build.legs_for_export(L, unblind=True)) == len(L)


def test_track_ids_missing_d0z0_is_hard_error():
    g = {f"d_pt{i}": np.ones(2, np.float32) for i in range(3)}
    g.update({f"d_eta{i}": np.zeros(2, np.float32) for i in range(3)})
    with pytest.raises(SystemExit):
        build.track_ids({"pairs": {"track_id": "d0z0"}}, g)
    t = build.track_ids({"pairs": {"track_id": "pt_eta"}}, g)   # jawnie stare → działa z ostrzeżeniem
    assert len(t) == 3


def test_track_id_is_injective_and_nonnegative():
    pt = np.array([1.0, 1.0], np.float32)
    d0 = np.array([0.5, -0.5], np.float32)                   # różnią się tylko bitem znaku
    z0 = np.array([0.25, 0.25], np.float32)
    a = build._trk_id_d0z0(pt, d0, z0)
    assert a[0] != a[1] and (a >= 0).all()   # PF3-FIN: ID ≥ 0 także dla d0 < 0
    # (d0, z0) i (z0, d0) nie mogą się zlać: górne 32 bity = d0, dolne = z0
    b = build._trk_id_d0z0(pt, z0, d0)
    assert a[0] != b[0]


def test_n_dstar_data_uses_refit_and_keeps_old(monkeypatch):
    rng = np.random.default_rng(0)
    n_sig, n_bkg = 3000, 20000
    dm_sig = rng.normal(0.14543, 0.0006, n_sig)
    u = rng.gamma(1.5, 0.004, n_bkg)
    dm_bkg = 0.13957 + u
    dm = np.r_[dm_sig, dm_bkg]
    m = M0[10] / 1000 + dm                      # D.m − D.dm = m(D0) → w oknie
    D = pd.DataFrame({"m": m, "dm": dm, "truth": np.r_[np.zeros(n_sig), np.full(n_bkg, 9.0)]})
    p = dict(build.DEFAULTS)
    f = build.n_dstar_data(D, p, D_mc=D)
    assert f["method"].startswith("dstar_refit")
    assert "old" in f and f["old"]["n_entries"] > 0
    assert abs(f["N"] - n_sig) < 5 * max(f["N_err"], 30)


def test_synth_has_same_event_and_shared_track_pairs(tmp_path):
    import h5py
    from charmtag import synth
    synth.make(str(tmp_path), n=1500, seed=3)
    with h5py.File(tmp_path / "data_428770_d150_Donly_none.h5") as f:
        ev = f["events"]
        key = ev["run_number"][:].astype(np.int64) * 10**9 + ev["event_number"][:].astype(np.int64)
        d0 = ev["Dmeson_daughter0_d0"][:]
        z0 = ev["Dmeson_daughter0_z0SinThetaPV"][:]
    df = pd.DataFrame({"k": key, "d0": d0, "z0": z0})
    multi = df.groupby("k").size()
    assert (multi >= 2).sum() > 50                       # zdarzenia z ≥2 kandydatami → pary same-event
    shared = df.groupby(["k", "d0", "z0"]).size()
    assert (shared >= 2).sum() > 10                      # ten sam ślad (bity d0,z0) w dwóch kandydatach tego samego zdarzenia
