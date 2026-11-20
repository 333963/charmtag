"""Etap DD: testy charmtag/dedup.py."""
import numpy as np
import pandas as pd

from charmtag import dedup


def _frame(n=50, seed=1):
    rng = np.random.default_rng(seed)
    L = pd.DataFrame({"run": np.zeros(n, np.int64), "event": rng.integers(0, 10, n), "decay_type": 10,
                      "m": rng.uniform(1.75, 2.0, n).astype(np.float32), "pt": rng.uniform(3, 30, n).astype(np.float32)})
    for i in range(3):
        for v in ("pt", "eta", "phi", "d0", "z0"):
            L[f"d_{v}{i}"] = (rng.normal(size=n).astype(np.float32) if i < 2 else np.float32(np.nan))
    for i in range(3):
        L[f"d_pt{i}"] = np.abs(L[f"d_pt{i}"]) + np.float32(0.5) if i < 2 else np.float32(0.0)
    return L


def test_no_duplicates_untouched():
    L = _frame()
    r = dedup.report(L)
    assert r["dup_rows"] == 0 and r["n_unique"] == len(L)
    assert dedup.dedup_mask(L).all()
    assert len(dedup.apply(L, True)) == len(L)


def test_exact_duplicates_removed_first_kept():
    L = _frame()
    D = pd.concat([L, L.iloc[[3, 7, 7]]], ignore_index=True)  # wiersz 3 ×2, wiersz 7 ×3
    r = dedup.report(D)
    assert r["dup_rows"] == 2 + 3 and r["max_mult"] == 3 and r["n_unique"] == len(L)
    keep = dedup.dedup_mask(D)
    assert keep[: len(L)].all() and not keep[len(L):].any()
    out = dedup.apply(D, True)
    assert len(out) == len(L) and out.attrs["dedup_report"]["dup_rows"] == 5


def test_other_hypothesis_is_not_duplicate():
    """Ta sama para torów, inna masa (K↔π) = dwa różne kandydaci — nie usuwać."""
    L = _frame(n=5)
    S = L.copy()
    S["m"] = (S["m"] + np.float32(0.02)).astype(np.float32)
    for v in ("pt", "eta", "phi", "d0", "z0"):  # zamiana kolejności córek = zamiana hipotezy K↔π
        S[f"d_{v}0"], S[f"d_{v}1"] = L[f"d_{v}1"].to_numpy(), L[f"d_{v}0"].to_numpy()
    D = pd.concat([L, S], ignore_index=True)
    assert dedup.report(D)["dup_rows"] == 0


def test_refit_copy_is_duplicate_with_d0z0_key():
    """DD2: ta sama para śladów (d0, z0 identyczne), pędy po refit różne o 0,1 % → duplikat (stary klucz by go nie widział)."""
    L = _frame(n=8)
    R = L.copy()
    for i in range(2):
        R[f"d_pt{i}"] = (R[f"d_pt{i}"] * np.float32(1.001)).astype(np.float32)
        R[f"d_eta{i}"] = (R[f"d_eta{i}"] + np.float32(1e-4)).astype(np.float32)
    R["m"] = (R["m"] + np.float32(1e-4)).astype(np.float32)
    D = pd.concat([L, R], ignore_index=True)
    assert dedup.report(D)["dup_rows"] == 2 * len(L)
    assert len(dedup.apply(D, True)) == len(L)


def test_fallback_key_without_d0z0(caplog):
    L = _frame(n=6).drop(columns=[c for c in dedup._TRACK_COLS])
    D = pd.concat([L, L.iloc[[1]]], ignore_index=True)
    assert dedup.report(D)["dup_rows"] == 2
    assert any("klucz zapasowy" in r.message for r in caplog.records)


def test_report_mode_does_not_drop():
    L = _frame()
    D = pd.concat([L, L.iloc[[0]]], ignore_index=True)
    out = dedup.apply(D, "report")
    assert len(out) == len(D) and out.attrs["dedup_report"]["dup_rows"] == 2
    assert dedup.apply(D, False) is D


def test_bits_nan_stable():
    a = dedup._bits(np.array([np.nan, 1.0], dtype=np.float32))
    b = dedup._bits(np.array([np.nan, 1.0], dtype=np.float32))
    assert (a == b).all()
