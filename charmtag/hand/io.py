"""Odczyt próbek: h5 (v1.1 – założenie: grupa/kolumna lub pandas HDFStore) oraz npz (S6_kand/sub, S_bkggen/sub)."""
from __future__ import annotations
import numpy as np, pandas as pd, h5py

def read_h5(path: str, group: str | None = None) -> pd.DataFrame:
    try:
        with pd.HDFStore(path, "r") as st:
            keys = st.keys()
            if keys:
                k = group if group and f"/{group}" in keys else keys[0]
                return st[k]
    except Exception:
        pass
    with h5py.File(path, "r") as f:
        g = f[group] if group and group in f else f
        # jeśli brak datasetów na tym poziomie, zejdź do pierwszej grupy
        if not any(isinstance(g[k], h5py.Dataset) for k in g):
            g = g[next(iter(g))]
        cols = {}
        for k in g:
            d = g[k]
            if isinstance(d, h5py.Dataset):
                a = d[()]
                cols[k] = list(a) if a.ndim > 1 else a
        return pd.DataFrame(cols)

def read_npz(path: str) -> pd.DataFrame:
    """Obsługuje: (a) tablice nazwane po kolumnach, (b) X + feature_names (+ y/label/weight)."""
    z = np.load(path, allow_pickle=True)
    keys = list(z.keys())
    if "X" in keys:
        names = [str(s) for s in z["feature_names"]] if "feature_names" in keys else [f"f{i}" for i in range(z["X"].shape[1])]
        df = pd.DataFrame(np.asarray(z["X"]), columns=names)
        for k in keys:
            if k not in ("X", "feature_names") and np.ndim(z[k]) == 1 and len(z[k]) == len(df):
                df[k] = z[k]
        return df
    n = None
    cols = {}
    for k in keys:
        a = z[k]
        if a.ndim == 0:
            continue
        if n is None:
            n = len(a)
        if len(a) == n:
            cols[k] = list(a) if a.ndim > 1 else a
    return pd.DataFrame(cols)

def label_column(df: pd.DataFrame) -> str | None:
    for c in ("label", "y", "is_signal", "isSignal", "truth", "sig", "target"):
        if c in df.columns:
            return c
    return None
