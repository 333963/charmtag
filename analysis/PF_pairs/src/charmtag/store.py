"""Proste h5 (kolumna-na-dataset) z dopisywaniem porcjami — bez pytables."""
from __future__ import annotations

import h5py
import numpy as np
import pandas as pd


class Appender:
    def __init__(self, path: str, attrs: dict | None = None):
        self.f = h5py.File(path, "w")
        for k, v in (attrs or {}).items():
            self.f.attrs[k] = v

    def add(self, df: pd.DataFrame):
        if len(df) == 0:
            return
        for c in df.columns:
            a = df[c].to_numpy()
            if a.dtype == object:
                a = a.astype("S")
            if c not in self.f:
                self.f.create_dataset(c, data=a, maxshape=(None,), chunks=True, compression="gzip", compression_opts=1)
            else:
                d = self.f[c]
                n = d.shape[0]
                d.resize(n + len(a), axis=0)
                d[n:] = a

    def close(self):
        self.f.close()


def save(path: str, df: pd.DataFrame, attrs: dict | None = None):
    a = Appender(path, attrs)
    a.add(df)
    a.close()


def load(path: str, columns: list[str] | None = None) -> pd.DataFrame:
    with h5py.File(path, "r") as f:
        cols = columns or list(f.keys())
        return pd.DataFrame({c: f[c][:] for c in cols if c in f})
