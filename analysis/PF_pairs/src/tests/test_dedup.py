"""Etap DD: testy charmtag/dedup.py."""
import numpy as np
import pandas as pd

from charmtag import dedup


def _frame(n=50, seed=1):
    rng = np.random.default_rng(seed)
    L = pd.DataFrame({"run": np.zeros(n, np.int64), "event": rng.integers(0, 10, n), "decay_type": 10,
                      "m": rng.uniform(1.75, 2.0, n).astype(np.float32), "pt": rng.uniform(3, 30, n).astype(np.float32)})
    for i in range(3):
        for v in ("pt", "eta", "phi"):
            L[f"d_{v}{i}"] = (rng.normal(size=n).astype(np.float32) if i < 2 else np.float32(np.nan))
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
    D = pd.concat([L, S], ignore_index=True)
    assert dedup.report(D)["dup_rows"] == 0


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
