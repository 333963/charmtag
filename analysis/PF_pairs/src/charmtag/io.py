"""Odczyt h5 v1.1 (dataset `events`, compound) lub kolumna-na-dataset; rozwiązywanie nazw kolumn."""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

import h5py
import numpy as np
import pandas as pd


def list_columns(path, dataset: str = "events") -> list[str]:
    if isinstance(path, (list, tuple)):  # etap P: wiele plików → kolumny z pierwszego
        path = path[0]
    with h5py.File(path, "r") as f:
        if dataset in f and isinstance(f[dataset], h5py.Dataset) and f[dataset].dtype.names:
            return list(f[dataset].dtype.names)
        cols = []

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and obj.ndim == 1 and obj.dtype.names is None:
                cols.append(name.split("/")[-1])

        f.visititems(visit)
        return cols


def n_rows(path, dataset: str = "events") -> int:
    if isinstance(path, (list, tuple)):  # etap P
        return int(sum(n_rows(q, dataset) for q in path))
    with h5py.File(path, "r") as f:
        if dataset in f and isinstance(f[dataset], h5py.Dataset):
            return int(f[dataset].shape[0])
        n = 0

        def visit(name, obj):
            nonlocal n
            if isinstance(obj, h5py.Dataset) and obj.ndim == 1:
                n = max(n, obj.shape[0])

        f.visititems(visit)
        return n


@dataclass
class ColumnMap:
    """Nazwa logiczna → nazwa kolumny w pliku (lub None)."""
    mapping: dict = field(default_factory=dict)
    missing: list = field(default_factory=list)

    def __getitem__(self, k):
        return self.mapping[k]

    def get(self, k, default=None):
        return self.mapping.get(k, default)

    def has(self, k):
        return self.mapping.get(k) is not None


def resolve_columns(available: list[str], spec: dict, n_daughters: int = 3) -> ColumnMap:
    avail = set(available)
    lower = {c.lower(): c for c in available}
    cm = ColumnMap()
    for logical, cands in spec.items():
        if isinstance(cands, str):
            cands = [cands]
        if any("{i}" in c for c in cands):
            for i in range(n_daughters):
                key = f"{logical}{i}"
                cm.mapping[key] = _first(avail, lower, [c.format(i=i) for c in cands])
                if cm.mapping[key] is None:
                    cm.missing.append(key)
        else:
            cm.mapping[logical] = _first(avail, lower, cands)
            if cm.mapping[logical] is None:
                cm.missing.append(logical)
    return cm


def _first(avail, lower, cands):
    for c in cands:
        if c in avail:
            return c
    for c in cands:
        if c.lower() in lower:
            return lower[c.lower()]
    for c in cands:  # regex (np. "re:Dmeson_.*Chi2")
        if c.startswith("re:"):
            pat = re.compile(c[3:])
            for a in sorted(avail):
                if pat.fullmatch(a):
                    return a
    return None


def read_rows(path: str, columns: list[str], dataset: str = "events", chunk_rows: int = 2_000_000,
              row_filter=None, max_rows: int | None = None) -> pd.DataFrame:
    """Czyta wybrane kolumny porcjami; row_filter(df_chunk) -> maska bool (np. decayType).
    path może być listą plików (etap P) — czytane kolejno, max_rows liczone łącznie."""
    if isinstance(path, (list, tuple)):
        parts, left = [], (int(max_rows) if max_rows else None)
        for q in path:
            if left is not None and left <= 0:
                break
            n_q = n_rows(q, dataset)
            parts.append(read_rows(q, columns, dataset, chunk_rows, row_filter, left))
            if left is not None:
                left -= n_q
        return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=list(dict.fromkeys(c for c in columns if c)))
    columns = list(dict.fromkeys(c for c in columns if c))
    out = []
    with h5py.File(path, "r") as f:
        compound = dataset in f and isinstance(f[dataset], h5py.Dataset) and f[dataset].dtype.names
        if compound:
            ds = f[dataset]
            n = ds.shape[0]
        else:
            dsets = {}

            def visit(name, obj):
                if isinstance(obj, h5py.Dataset) and obj.ndim == 1:
                    dsets.setdefault(name.split("/")[-1], obj)

            f.visititems(visit)
            n = max(d.shape[0] for d in dsets.values())
        if max_rows:
            n = min(n, int(max_rows))
        for start in range(0, n, chunk_rows):
            stop = min(n, start + chunk_rows)
            if compound:
                arr = ds.fields(columns)[start:stop] if len(columns) > 1 else ds[start:stop][columns]
                df = pd.DataFrame({c: np.asarray(arr[c]) for c in columns})
            else:
                df = pd.DataFrame({c: np.asarray(dsets[c][start:stop]) for c in columns})
            if row_filter is not None:
                df = df[row_filter(df)]
            out.append(df.reset_index(drop=True))
    if not out:
        return pd.DataFrame(columns=columns)
    return pd.concat(out, ignore_index=True)


def expand(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path))


# ----------------------------------------------------------------------------- etap P: wiele plików
class _FieldsView:
    def __init__(self, me, cols):
        self.me, self.cols = me, cols

    def __getitem__(self, sl):
        return self.me._read(sl, self.cols)


class MultiEvents:
    """Wirtualny compound `events` z listy plików (sklejonych w kolejności listy), indeks globalny.
    Wiersz globalny g ↔ plik i = searchsorted(offsets, g, 'right')-1, wiersz lokalny g - offsets[i].
    Obsługuje ds.shape, ds.fields(cols)[a:b], ds[a:b] (krok 1) — tyle, ile używa charmtag."""

    def __init__(self, paths, dataset: str = "events"):
        self.paths = list(paths)
        self.files = [h5py.File(p, "r") for p in self.paths]
        self.dsets = []
        for p, f in zip(self.paths, self.files):
            d = f[dataset] if dataset in f else None
            if not (isinstance(d, h5py.Dataset) and d.dtype.names):
                self.close()
                raise ValueError(f"{p}: odczyt wielu plików wymaga compound `{dataset}`")
            self.dsets.append(d)
        self.offsets = np.concatenate([[0], np.cumsum([d.shape[0] for d in self.dsets])]).astype(np.int64)
        self.shape = (int(self.offsets[-1]),)
        self.dtype = self.dsets[0].dtype if self.dsets else None

    def __len__(self):
        return self.shape[0]

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def close(self):
        for f in self.files:
            try:
                f.close()
            except Exception:
                pass

    def fields(self, cols):
        return _FieldsView(self, [cols] if isinstance(cols, str) else list(cols))

    def __getitem__(self, sl):
        return self._read(sl, None)

    def _read(self, sl, cols):
        if not isinstance(sl, slice):
            raise TypeError("MultiEvents: tylko wycinki a:b")
        start, stop, step = sl.indices(self.shape[0])
        if step != 1:
            raise ValueError("MultiEvents: krok != 1")
        out = []
        for i, d in enumerate(self.dsets):
            a, b = max(start, self.offsets[i]), min(stop, self.offsets[i + 1])
            if a >= b:
                continue
            a, b = int(a - self.offsets[i]), int(b - self.offsets[i])
            x = d.fields(cols)[a:b] if cols else d[a:b]
            if out and x.dtype != out[0].dtype:
                x = x.astype(out[0].dtype)
            out.append(x)
        if not out:
            d0 = self.dsets[0]
            return (d0.fields(cols)[0:0] if cols else d0[0:0])
        return out[0] if len(out) == 1 else np.concatenate(out)

    def locate(self, rows):
        """Globalne wiersze → (indeks pliku, wiersz lokalny)."""
        rows = np.asarray(rows, np.int64)
        i = np.searchsorted(self.offsets, rows, side="right") - 1
        return i, rows - self.offsets[i]


def open_events(paths, dataset: str = "events") -> MultiEvents:
    """Jeden plik lub lista → MultiEvents (context manager)."""
    if isinstance(paths, (str, os.PathLike)):
        paths = [paths]
    return MultiEvents(paths, dataset)
